"""ParamCheck：物理干预参数报告预检.

设计继承 vec-trust-kit 的 format_check.py（本地自查、分级 PASS/FAIL/SKIP/WARN、
输出"线索"而非官方结论），但检查对象从 h5ad 提交文件改为**物理干预参数报告**：
在没有湿实验条件时，先在纸面判断一份参数的报告能不能支撑它宣称的效应。

为什么需要它：
  A2 证据库反复出现"来源未报告剂量"——不是数据缺憾，而是复现性硬伤。
  一份只写"红光照射""低频刺激"的报告，无法判断落在剂量窗的哪一段：
  PBM 3 J/cm2 有效而 30 J/cm2 抑制，tFUS MI>1.9 从神经调控变成显著风险。
  ParamCheck 把这些硬阈值显式化，让"参数写全了没有"可以被机器判定。

判定契约（依被测对象而变，禁止恒真）：
  COMPLETE     关键参数齐全、剂量带数值单位、落在安全窗 → 可作为可复现报告
  PARTIAL      有剂量但缺部分参数 → 结论受限，需补报告再判定
  INSUFFICIENT 关键参数缺失或只有相对描述 → 不可复现
  OUT_OF_WINDOW 参数在数值上越过硬安全阈值 → 不是"效果弱"而是风险等级变化
  NO_EVIDENCE  没有任何可核对参数 → 等同宣称，不可进入证据库

免责：这是按公开物理参数规格实现的**本地自查线索**，不是监管结论，也不能
替代 IRB / 设备安全评估。参数在窗内 ≠ 有效；参数越界 ≠ 一定有害，只表示
"已越过公开的非显著风险阈值，需要额外证据"。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
import re

DISCLAIMER = (
    "ParamCheck 为本地自查线索（按公开物理参数规格实现），非监管/临床结论；"
    "参数落在窗口内不代表有效，越界也不代表必然有害，只表示越过公开的非显著风险阈值。"
)

# 关键参数规格：(展示名, 接收键的别名元组, 单位, 安全区间或 None, 是否必需)
# 区间取公开共识的"非显著风险/常用剂量"范围，不是疗效保证。
PARAM_SPECS: dict[str, tuple[str, ...]] = {
    "photobiomodulation": (
        "wavelength", "power_density", "fluence", "exposure_time", "spot_area",
    ),
    "vagus nerve stimulation": (
        "frequency", "pulse_width", "current", "session_duration", "site",
    ),
    "focused ultrasound": (
        "frequency", "mechanical_index", "intensity", "pulse_duration",
        "duty_cycle", "temperature_rise",
    ),
    "mechanical loading": (
        "strain_or_shear", "frequency", "session_duration", "substrate_stiffness",
    ),
}

# 键别名：报告里可能出现的字段名都归一化到规范键
ALIASES: dict[str, tuple[str, ...]] = {
    "wavelength": ("wavelength", "wavelength_nm", "波长"),
    "power_density": ("power_density", "irradiance", "辐照度"),
    "fluence": ("fluence", "dose", "energy_density", "剂量", "能量密度"),
    "exposure_time": ("exposure_time", "duration", "照射时间"),
    "spot_area": ("spot_area", "光斑面积"),
    "frequency": ("frequency", "freq", "频率"),
    "pulse_width": ("pulse_width", "脉宽"),
    "current": ("current", "amplitude", "电流"),
    "session_duration": ("session_duration", "duration", "时长", "min_per_day"),
    "site": ("site", "target", "刺激部位"),
    "mechanical_index": ("mechanical_index", "mi", "mitc", "机械指数"),
    "intensity": ("intensity", "isppa", "ispta", "声强"),
    "pulse_duration": ("pulse_duration", "脉冲宽度"),
    "duty_cycle": ("duty_cycle", "占空比"),
    "temperature_rise": ("temperature_rise", "delta_t", "温升"),
    "strain_or_shear": ("strain_or_shear", "strain", "shear_stress", "应变", "剪切力"),
    "substrate_stiffness": ("substrate_stiffness", "stiffness", "基质刚度"),
}

# 硬安全阈值：越界即 OUT_OF_WINDOW（不是扣分，是判定改级）
# (规范键, 比较方式, 阈值, 说明)
HARD_LIMITS: dict[str, tuple[str, float, str]] = {
    "photobiomodulation": ("fluence", ">", 200, "能量密度 >200 J/cm2 超出常用 PBM 剂量窗，可能进入抑制/损伤段"),
    "focused ultrasound": ("mechanical_index", ">", 1.9, "MI>1.9 越过 FDA/ITRUSST 非显著风险阈值"),
}
# 温升单列（键可能在报告里出现）
TEMP_LIMIT = ("temperature_rise", ">", 2.0, "温升>2°C 越过 ITRUSST 非显著热风险阈值")

# 归一化：单位写法的白名单（用于判断剂量是否"带单位"）
UNIT_HINTS = (
    "nm", "mw", "j/cm", "w/cm", "hz", "khz", "mhz", "ma", "us", "ms", "min",
    "s", "%", "dyn", "kpa", "mpa", "w", "cm2", "mm", "°c", "deg",
)
RELATIVE_WORDS = ("低剂量", "高剂量", "适量", "low", "high", "moderate", "relative", "some", "若干", "一定")


@dataclass
class Check:
    name: str
    status: str  # PASS / FAIL / SKIP / WARN
    detail: str

    def to_dict(self) -> dict:
        return {"name": self.name, "status": self.status, "detail": self.detail}


@dataclass
class ParamReport:
    modality: str
    source_id: str
    checks: list[Check] = field(default_factory=list)
    verdict: str = "NO_EVIDENCE"
    reason: str = ""
    falsification: str = ""

    @property
    def passed(self) -> bool:
        return self.verdict in ("COMPLETE",) and not any(c.status == "FAIL" for c in self.checks)

    def to_dict(self) -> dict:
        return {
            "modality": self.modality,
            "source_id": self.source_id,
            "verdict": self.verdict,
            "reason": self.reason,
            "falsification": self.falsification,
            "passed": self.passed,
            "checks": [c.to_dict() for c in self.checks],
            "disclaimer": DISCLAIMER,
        }


def normalize(report: dict) -> dict:
    """把报告字段按别名归一化；大小写不敏感，值转字符串。"""
    out: dict[str, str] = {}
    lower = {str(k).lower(): str(v).strip() for k, v in report.items() if v is not None}
    for canon, names in ALIASES.items():
        for n in names:
            if n.lower() in lower and lower[n.lower()]:
                out[canon] = lower[n.lower()]
                break
    return out


def _numeric(value: str):
    """从字符串里抽第一个数字；失败返回 None。"""
    m = re.search(r"-?\d+(?:\.\d+)?", value.replace(",", ""))
    return float(m.group()) if m else None


def _has_unit(value: str) -> bool:
    v = value.lower()
    return any(u in v for u in UNIT_HINTS)


def _is_relative(value: str) -> bool:
    v = value.lower()
    return any(w in v for w in RELATIVE_WORDS)


def check(report: dict, *, source_id: str = "") -> ParamReport:
    """对一份参数报告做预检。

    report 至少要有 modality；其余参数键见 ALIASES。
    """
    modality = str(report.get("modality", "")).strip().lower()
    rep = ParamReport(modality=modality, source_id=source_id or str(report.get("id", "")))

    # 0) 模态必须已知
    if modality not in PARAM_SPECS:
        rep.checks.append(Check("modality", "FAIL", f"未知模态 {modality!r}；支持：{', '.join(PARAM_SPECS)}"))
        rep.verdict = "NO_EVIDENCE"
        rep.reason = "模态未知，无法套用参数规格"
        rep.falsification = "补齐为一个已知模态后重新预检，判定必须随之改变"
        return rep
    rep.checks.append(Check("modality", "PASS", f"已知模态：{modality}"))

    norm = normalize(report)
    specs = PARAM_SPECS[modality]
    present = [k for k in specs if k in norm]
    missing = [k for k in specs if k not in norm]
    rep.checks.append(Check(
        "key_parameters",
        "PASS" if not missing else "FAIL",
        f"关键参数 {len(present)}/{len(specs)} 齐备" + (f"；缺：{', '.join(missing)}" if missing else ""),
    ))

    # 1) 剂量字段：有数值 + 单位 vs 只有相对描述 vs 完全缺失
    dose_key = {"photobiomodulation": "fluence",
                "vagus nerve stimulation": "current",
                "focused ultrasound": "intensity",
                "mechanical loading": "strain_or_shear"}[modality]
    dose = norm.get(dose_key, "")
    if not dose:
        rep.checks.append(Check("dose", "FAIL", f"缺核心剂量参数 {dose_key}（不可复现）"))
    elif _is_relative(dose) and _numeric(dose) is None:
        rep.checks.append(Check("dose", "FAIL", f"剂量 {dose!r} 仅为相对描述，无数值，不可复现"))
    elif not _has_unit(dose):
        rep.checks.append(Check("dose", "WARN", f"剂量 {dose!r} 有数字但未标单位，跨研究不可比"))
    else:
        rep.checks.append(Check("dose", "PASS", f"剂量带数值和单位：{dose}"))

    # 2) 硬阈值：越界直接改判 OUT_OF_WINDOW
    limit = HARD_LIMITS.get(modality)
    limit_hits: list[str] = []
    for key, op, threshold, message in ([limit] if limit else []) + [TEMP_LIMIT]:
        if not key or key not in norm:
            continue
        num = _numeric(norm[key])
        if num is None:
            continue
        if (op == ">" and num > threshold) or (op == "<" and num < threshold):
            limit_hits.append(f"{key}={num} → {message}")
    if limit_hits:
        rep.checks.append(Check("hard_limits", "FAIL", "；".join(limit_hits)))
    else:
        checked = [k for k, _, _, _ in ([limit] if limit else []) + [TEMP_LIMIT] if k and k in norm]
        rep.checks.append(Check("hard_limits", "PASS" if checked else "SKIP",
                                f"硬阈值核对 {len(checked)} 项，未越界" if checked else "报告未给出可核对的硬阈值参数"))

    # 3) 报告是否给出细胞/靶点与终点（来自原证据字段，若有）
    for key, label in (("cell_type", "靶细胞"), ("endpoint", "终点")):
        val = str(report.get(key, "")).strip()
        rep.checks.append(Check(key, "PASS" if val else "WARN",
                                f"{label}已标注：{val[:60]}" if val else f"未标注{label}，无法判断是否作用于目标"))

    # 4) 来源可核对
    url = str(report.get("source_url", "")).strip()
    rep.checks.append(Check("source_url", "PASS" if url.startswith("http") else "WARN",
                            "有可核对来源" if url.startswith("http") else "无可核对一手来源"))

    fails = [c for c in rep.checks if c.status == "FAIL"]
    warns = [c for c in rep.checks if c.status == "WARN"]

    # 判定：依报告内容而变，带明确推翻条件
    if limit_hits:
        rep.verdict = "OUT_OF_WINDOW"
        rep.reason = "参数越过公开硬安全阈值；这不是效应强弱问题而是风险等级变化"
    elif not present:
        rep.verdict = "NO_EVIDENCE"
        rep.reason = "没有任何可核对的物理参数，等同宣称"
    elif fails:
        rep.verdict = "INSUFFICIENT"
        rep.reason = f"关键参数/剂量不完整（{len(fails)} 项失败），当前不可复现"
    elif warns:
        rep.verdict = "PARTIAL"
        rep.reason = f"核心参数齐但 {len(warns)} 项待补，结论受限"
    else:
        rep.verdict = "COMPLETE"
        rep.reason = "关键参数齐备、剂量带数值单位、未越过硬阈值"

    rep.falsification = (
        f"若补充/更正以下任一项，判定必须改变：{', '.join(missing) if missing else '（关键参数已齐）'}；"
        "硬阈值越界时任何'有效'叙事都不能抵消风险改判。"
    )
    return rep


# 从非结构化 signal 文本抽取参数的正则（证据库把参数写在自由文本里）
# 每条：(规范键, 正则)。抽到就作为该键的候选值。
SIGNAL_PATTERNS: dict[str, tuple[str, ...]] = {
    "wavelength": (r"(\d+\s*nm)",),
    "power_density": (r"(\d+(?:\.\d+)?\s*mW\s*/\s*cm\s*2?)",),
    "fluence": (r"(\d+(?:\.\d+)?\s*J\s*/\s*cm\s*2?)",),
    "frequency": (r"(\d+(?:\.\d+)?\s*(?:kHz|MHz|Hz))", r"(\d+\s*Hz)"),
    "pulse_width": (r"(\d+(?:\.\d+)?\s*(?:us|µs|ms))",),
    "current": (r"(\d+(?:\.\d+)?\s*mA)",),
    "mechanical_index": (r"(?:MI|MItc)\s*[<=]?\s*(\d+(?:\.\d+)?)", r"MI\s*>?\s*(\d+(?:\.\d+)?)"),
    "intensity": (r"(\d+(?:\.\d+)?\s*W\s*/\s*cm\s*2?)",),
    "duty_cycle": (r"(\d+\s*%)",),
    "temperature_rise": (r"(\d+(?:\.\d+)?\s*(?:°C|C\b|度))",),
    "session_duration": (r"(\d+\s*min)", r"(\d+\s*sessions?)", r"(\d+\s*h/day)"),
    "exposure_time": (r"(\d+\s*(?:min|s)\b)",),
    "strain_or_shear": (r"(\d+(?:\.\d+)?\s*%)", r"(\d+\s*dyn\s*/\s*cm)", r"(\d+\s*kPa)"),
    "substrate_stiffness": (r"(\d+\s*kPa)",),
    "spot_area": (r"(\d+\s*cm\s*2)",),
}


def extract_from_signal(signal: str, modality: str) -> dict[str, str]:
    """从自由文本 signal 抽取可核对参数；只补空缺，不覆盖显式字段。"""
    found: dict[str, str] = {}
    for canon, patterns in SIGNAL_PATTERNS.items():
        # 只处理该模态关心的参数，避免跨模态误抽（如 % 被当成应变）
        if canon not in PARAM_SPECS.get(modality, ()):
            continue
        for pat in patterns:
            m = re.search(pat, signal, re.IGNORECASE)
            if m:
                val = m.group(1) if m.lastindex else m.group(0)
                found.setdefault(canon, val.strip())
                break
    return found


def check_evidence_row(row: dict) -> ParamReport:
    """从证据库一行生成参数报告，便于批量预检现有 103 条。

    证据库把参数写在 signal 自由文本 + dose/dose_unit 分列里；
    这里先按正则抽取，再让显式字段覆盖，兼容两种写法。
    """
    modality = str(row.get("modality", "")).strip().lower()
    payload = dict(row)
    extracted = extract_from_signal(str(row.get("signal", "")), modality)
    for k, v in extracted.items():
        payload.setdefault(k, v)
    # dose + dose_unit 合并给 ParamCheck 判单位
    if row.get("dose") and row.get("dose_unit"):
        payload["dose"] = f'{row["dose"]} {row["dose_unit"]}'
        # 数值剂量同时回填到该模态的核心剂量键（若 signal 未抽到）
        dose_key = {"photobiomodulation": "fluence",
                    "vagus nerve stimulation": "current",
                    "focused ultrasound": "intensity",
                    "mechanical loading": "strain_or_shear"}.get(modality)
        # 只有当 dose_unit 是"真剂量单位"而非波长/相对时回填
        unit = row.get("dose_unit", "").lower()
        if dose_key and any(u in unit for u in ("j/cm", "mw", "w/cm", "ma", "hz", "%", "dyn", "kpa")):
            payload.setdefault(dose_key, f'{row["dose"]} {row["dose_unit"]}')
    return check(payload, source_id=row.get("id", ""))


def render(rep: ParamReport) -> str:
    lines = [f"参数报告预检：{rep.source_id or '(未命名)'}  模态={rep.modality}",
             f"判定：{rep.verdict}  理由：{rep.reason}", ""]
    for c in rep.checks:
        lines.append(f"  [{c.status:<4}] {c.name}: {c.detail}")
    lines += ["", f"推翻条件：{rep.falsification}", "", DISCLAIMER]
    return "\n".join(lines)


def to_json(rep: ParamReport, *, indent: int = 2) -> str:
    return json.dumps(rep.to_dict(), ensure_ascii=False, indent=indent)
