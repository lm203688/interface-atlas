"""DoseAudit：物理干预剂量报告完整性引擎（L1 生态入口）.

为什么存在：
    证据库与论文普遍只报"设备参数"，且连设备参数都报不全（实测 103 条：
    PBM 37 条仅 19 条报波长、6 条报能量密度、2 条报光斑面积、1 条报照射时间；
    VNS 18 条 signal 几乎只有 "taVNS" 这类标签）。这使"不同研究能不能比较"
    在源头就不可判定——PEMF/EMTT 系统综述（Biomedicines 2026）因此无法做
    meta 分析并喊出 urgent call to action。

    根因不是数据不够，是没有剂量报告标准。放疗靠 ICRU 报告标准 + DICOM-RT
    + TPS 解决同一问题。本模块是物理干预版的第一步：定义四层剂量本体，

    并用它把"报全了没有 / 能不能复现"变成机器可判定的输出。

四层剂量本体（TIER）—— 这是核心贡献，也是 L2 计算引擎的接口契约：
    T1 设备参数 (device)      设备读数直接给出：波长、频率、功率、脉宽……
    T2 界面条件 (interface)   施加器到靶区的传递条件：距离、角度、线圈几何、
                              颅骨衰减、电极位置……（文献最常漏报的隐藏变量）
    T3 靶区剂量 (target)      真正到达靶区的剂量场：靶区 fluence / 声压 / 电场、
                              靶区深度、组织光学/声学/电学参数
    T4 个体化解剖 (subject)   是否用个体影像/体模，而非群体平均值

    T1 + T2 = 剂量场方程的完备边界条件；T3 是该方程的解。
    → 有 T1+T2 无 T3，可交给 L2 计算引擎求解（reproducible_with_model）；
    → 缺 T2，方程欠定（underdetermined），再多 T3 也无法外推到其他个体；
    → 缺 T1，无输入（not_reproducible）。

判定契约（依被测报告而变，禁止恒真）：
    reproducible_from_report  T1/T2/T3 全齐，报告自带靶区剂量，可直接复现
    reproducible_with_model   T1/T2 齐但缺 T3，边界条件完备，可由剂量场引擎求解
    underdetermined           T1 齐但 T2 缺，界面条件缺失，剂量场欠定
    not_reproducible          T1 缺，无可用输入，等同宣称

与 paramcheck 的关系：
    paramcheck 保留五级"判决"（COMPLETE/PARTIAL/...）不放宽；本模块是正交的
    "完整性诊断"——输出 completeness_score、missing_params、per-tier 分数。
    两者并存，前者回答"能不能过"，后者回答"缺什么、缺到什么程度"。

免责：本规范是基于公开文献的社区草案，非 ICRU/FDA/ISO 官方标准；分数是
报告完整性的度量，不是疗效或安全性的判定。
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import re
from pathlib import Path

DISCLAIMER = (
    "DoseAudit 基于公开文献的社区草案规范 v1，非 ICRU/FDA/ISO 官方标准；"
    "completeness_score 度量报告完整性，不代表疗效或安全性。"
)

TIER_NAMES = {
    1: "device（设备参数）",
    2: "interface（界面条件）",
    3: "target（靶区剂量）",
    4: "subject（个体化解剖）",
}
TIER_WEIGHTS = {1: 0.40, 2: 0.25, 3: 0.25, 4: 0.10}


@dataclass(frozen=True)
class DoseField:
    """一个剂量学必报字段。"""
    key: str
    label: str
    tier: int
    unit: str
    aliases: tuple[str, ...] = ()
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "key": self.key, "label": self.label, "tier": self.tier,
            "unit": self.unit, "aliases": list(self.aliases), "note": self.note,
        }


def _f(key, label, tier, unit, *aliases, note="") -> DoseField:
    return DoseField(key, label, tier, unit, aliases, note)


# ---------------------------------------------------------------------------
# 四层剂量本体：每模态的必报字段
# 单位与阈值取公开共识；区间不是疗效保证，只是"可复现所需的最小报告粒度"。
# ---------------------------------------------------------------------------
DOSE_SPEC: dict[str, tuple[DoseField, ...]] = {
    "photobiomodulation": (
        _f("wavelength", "波长", 1, "nm", "wavelength_nm", "波长"),
        _f("power_density", "辐照度/功率密度", 1, "mW/cm2", "irradiance", "功率密度", "辐照度"),
        _f("fluence", "能量密度", 1, "J/cm2", "energy_density", "剂量", "能量密度"),
        _f("exposure_time", "照射时长", 1, "s", "irradiation_time", "照射时间"),
        _f("spot_area", "光斑面积", 1, "cm2", "beam_area", "光斑面积"),
        _f("applicator_distance", "施加器-皮肤距离", 2, "cm", "source_distance", "照射距离"),
        _f("applicator_angle", "入射角", 2, "deg", "incidence_angle", "入射角度"),
        _f("beam_profile", "光束分布", 2, "enum: uniform|gaussian", "beam_uniformity"),
        _f("target_fluence", "靶区实际光通量", 3, "J/cm2", "target_dose", "靶区剂量"),
        _f("target_depth", "靶区深度", 3, "cm", "depth", "深度"),
        _f("tissue_optics", "组织光学参数", 3, "mu_a/mu_s (1/cm)", "optical_properties", "组织光学"),
        _f("anatomy_source", "解剖来源", 4, "enum: phantom|subject_imaging|population", "体模", "个体影像"),
    ),
    "vagus nerve stimulation": (
        _f("frequency", "频率", 1, "Hz", "freq", "频率"),
        _f("pulse_width", "脉宽", 1, "us", "pulse_duration", "脉宽"),
        _f("current", "电流强度", 1, "mA", "amplitude", "电流", "强度"),
        _f("session_duration", "单次时长", 1, "min", "duration_per_session", "时长"),
        _f("site", "刺激部位", 1, "enum: cervical|auricular|abdominal|other", "target", "部位", "刺激部位"),
        _f("electrode_geometry", "电极几何/触点配置", 2, "text", "lead_config", "电极配置"),
        _f("fascicle_target", "靶纤维束", 2, "text", "fiber_target", "纤维束"),
        _f("impedance", "接口阻抗", 2, "kOhm", "电极阻抗"),
        _f("neural_dose", "神经剂量(ECAP归一化)", 3, "mV or normalized ECAP", "ecap", "神经剂量"),
        _f("target_organ_engagement", "靶器官接 engagement", 3, "text", "靶器官激活", "脾神经激活"),
        _f("anatomy_source", "解剖来源", 4, "enum: population|subject_imaging", "个体化解剖"),
    ),
    "focused ultrasound": (
        _f("frequency", "频率", 1, "MHz", "freq", "频率"),
        _f("mechanical_index", "机械指数", 1, "MI", "mi", "mitc", "机械指数"),
        _f("intensity", "声强", 1, "W/cm2 or mW/cm2", "isppa", "ispta", "声强"),
        _f("pulse_duration", "脉冲宽度", 1, "ms or us", "脉宽"),
        _f("duty_cycle", "占空比", 1, "%", "占空比"),
        _f("temperature_rise", "温升", 1, "degC", "delta_t", "温升"),
        _f("transducer_geometry", "换能器几何", 2, "text", "探头几何", "换能器"),
        _f("skull_attenuation", "颅骨衰减", 2, "dB or text", "颅骨", "经颅衰减"),
        _f("target_pressure", "靶区声压", 3, "MPa", "target_pressure_field", "靶区声压"),
        _f("target_temperature", "靶区温度", 3, "degC", "靶区温升"),
        _f("tissue_acoustics", "组织声学参数", 3, "attenuation/speed (dB/cm/MHz)", "声学参数", "组织声学"),
        _f("anatomy_source", "解剖来源", 4, "enum: phantom|subject_imaging|population", "个体影像"),
    ),
    "mechanical loading": (
        _f("strain_or_shear", "应变/剪切应力", 1, "% or Pa or dyne/cm2", "strain", "shear_stress", "应变", "剪切力"),
        _f("frequency", "加载频率", 1, "Hz", "freq", "频率"),
        _f("session_duration", "单次时长", 1, "min", "duration", "时长"),
        _f("substrate_stiffness", "基质刚度", 1, "kPa", "stiffness", "刚度"),
        _f("loading_waveform", "加载波形", 2, "enum: static|cyclic|impact", "波形"),
        _f("contact_geometry", "接触几何", 2, "text", "加载面积", "接触面积"),
        _f("target_strain", "靶区实际应变", 3, "%", "target_deformation", "靶区应变"),
        _f("tissue_viscoelasticity", "组织粘弹性参数", 3, "text", "粘弹性", "组织力学参数"),
        _f("anatomy_source", "解剖来源", 4, "enum: phantom|subject_imaging|population", "个体影像"),
    ),
    # 代谢/微生物：无传统剂量单位，作为附录，用暴露量替代（v1 不纳入主评分）
    "metabolic / microbiota": (
        _f("exposure_amount", "摄入/暴露量", 1, "g/day or CFU or mmol/L", "dose", "摄入量"),
        _f("exposure_duration", "暴露时长", 1, "day or week", "duration", "疗程"),
        _f("compound_identity", "化合物/菌株身份", 2, "text", "菌株", "代谢物"),
        _f("target_concentration", "靶部位浓度", 3, "mmol/L or uM", "SCFA浓度", "靶区浓度"),
    ),
}

# 附录模态归一到主类
MODALITY_ALIASES = {
    "metabolic / bile acid": "metabolic / microbiota",
    "metabolic / incretin": "metabolic / microbiota",
}
APPENDIX_MODALITIES = {"metabolic / microbiota"}

# 各模态的 signal 正则（复用 paramcheck，避免重复维护）
SIGNAL_PATTERNS = {
    "photobiomodulation": {
        "wavelength": r"\d+\s*nm", "power_density": r"\d+(?:\.\d+)?\s*mW\s*/\s*cm",
        "fluence": r"\d+(?:\.\d+)?\s*J\s*/\s*cm", "exposure_time": r"\d+\s*(?:min|s)\b",
        "spot_area": r"\d+\s*cm\s*2",
    },
    "vagus nerve stimulation": {
        "frequency": r"\d+(?:\.\d+)?\s*(?:kHz|Hz)", "pulse_width": r"\d+(?:\.\d+)?\s*(?:us|µs|ms)",
        "current": r"\d+(?:\.\d+)?\s*mA", "session_duration": r"\d+\s*min",
    },
    "focused ultrasound": {
        "frequency": r"\d+(?:\.\d+)?\s*MHz", "mechanical_index": r"MI\s*[<= >]?\s*\d+(?:\.\d+)?",
        "intensity": r"\d+(?:\.\d+)?\s*(?:mW|W)\s*/\s*cm", "temperature_rise": r"\d+(?:\.\d+)?\s*(?:°C|C\b|度)",
        "duty_cycle": r"\d+\s*%",
    },
    "mechanical loading": {
        "strain_or_shear": r"\d+(?:\.\d+)?\s*(?:%|Pa|dyne|kPa)", "frequency": r"\d+\s*Hz",
        "substrate_stiffness": r"\d+\s*kPa",
    },
}

RELATIVE_WORDS = ("low", "high", "relative", "some", "若干", "一定", "lower", "higher", "适量", "低剂量", "高剂量")


def canonical_modality(modality: str) -> str:
    m = str(modality).strip().lower()
    return MODALITY_ALIASES.get(m, m)


@dataclass
class FieldCheck:
    key: str
    label: str
    tier: int
    unit: str
    status: str          # PRESENT / NO_UNIT / RELATIVE / MISSING
    value: str
    note: str = ""

    @property
    def score(self) -> float:
        return {"PRESENT": 1.0, "NO_UNIT": 0.5, "RELATIVE": 0.25, "MISSING": 0.0}[self.status]

    def as_dict(self) -> dict:
        return {
            "key": self.key, "label": self.label, "tier": self.tier, "unit": self.unit,
            "status": self.status, "value": self.value, "score": self.score, "note": self.note,
        }


@dataclass
class DoseAuditReport:
    source_id: str
    modality: str
    is_appendix: bool
    fields: list[FieldCheck]
    completeness_score: float
    tier_scores: dict[str, float]
    missing_params: list[str]
    reproducibility: str
    reason: str
    falsification: str

    def as_dict(self) -> dict:
        return {
            "schema": "dose_report_spec/v1",
            "source_id": self.source_id,
            "modality": self.modality,
            "is_appendix": self.is_appendix,
            "completeness_score": round(self.completeness_score, 3),
            "tier_scores": {k: round(v, 3) for k, v in self.tier_scores.items()},
            "missing_params": self.missing_params,
            "reproducibility": self.reproducibility,
            "reason": self.reason,
            "falsification": self.falsification,
            "fields": [f.as_dict() for f in self.fields],
            "disclaimer": DISCLAIMER,
        }


def _has_unit(value: str, field: DoseField) -> bool:
    """字段值是否带可核对单位（enum/text 字段只需非空）。

    通用单位表兜底：字段声明的单位写法可能与实际报告不同（如字段写
    "ms or us"，报告写 "10 ms"），因此用全局单位词表做子串命中。
    """
    if field.unit.startswith("enum:") or field.unit == "text":
        return True
    v = value.lower()
    units = {
        "nm", "mw", "j/cm", "w/cm", "mhz", "khz", "hz", "ma", "us", "µs", "ms",
        "min", "cm2", "cm", "mm", "deg", "°c", "%", "kpa", "mpa", "pa", "dyn",
        "db", "mmol", "um", "µm", "cfu", "g/day", "mv", "j", "w", "s", "v",
    }
    if any(u in v for u in units):
        return True
    tokens = re.findall(r"[a-zA-Zµ°]+", field.unit.lower())
    return any(t in v for t in tokens if len(t) >= 2 and t not in ("or", "text", "enum"))


def _collect_values(row: dict, modality: str) -> dict[str, str]:
    """从显式字段 + 别名 + signal 正则 + dose 列，收集原始值（不覆盖显式值）。"""
    out: dict[str, str] = {}
    lower = {str(k).strip().lower(): str(v).strip() for k, v in row.items() if v is not None}
    spec = DOSE_SPEC.get(modality, ())

    # 1) 显式字段与别名
    for f in spec:
        for name in (f.key, *f.aliases):
            v = lower.get(name.lower(), "")
            if v:
                out[f.key] = v
                break

    # 2) 正则抽取：signal 是主来源；dose / dose_unit / safety_boundary 也常藏参数
    #    （如 FUS 的 "MI<=1.9; ISPTA<=720 mW/cm2" 写在 dose 列，不是 signal）
    import re
    haystack = " ".join(str(row.get(k, "") or "") for k in
                        ("signal", "dose", "dose_unit", "safety_boundary"))
    for f in spec:
        if f.key in out:
            continue
        pat = SIGNAL_PATTERNS.get(modality, {}).get(f.key)
        if not pat:
            continue
        m = re.search(pat, haystack, re.IGNORECASE)
        if m:
            out[f.key] = m.group(0).strip()

    # 3) dose + dose_unit 回填到数值型剂量字段（优先级最低）
    #    不 break：一条 dose 文本可能含多个单位的剂量字段
    dose, unit = str(row.get("dose", "")).strip(), str(row.get("dose_unit", "")).strip()
    if dose and unit:
        combined = f"{dose} {unit}"
        for f in spec:
            if f.key in out or f.tier != 1:
                continue
            if f.key in ("fluence", "current", "strain_or_shear", "intensity",
                         "frequency", "exposure_amount", "session_duration", "exposure_time"):
                out[f.key] = combined
    return out


def audit(report: dict, *, source_id: str = "") -> DoseAuditReport:
    """对一份剂量报告做完整性审计。

    report 至少含 modality；其余键见 DOSE_SPEC（含别名）与 signal 文本。
    """
    modality = canonical_modality(report.get("modality", ""))
    spec = DOSE_SPEC.get(modality)
    rep_id = source_id or str(report.get("id", ""))
    if spec is None:
        return DoseAuditReport(
            source_id=rep_id, modality=modality, is_appendix=False, fields=[],
            completeness_score=0.0, tier_scores={}, missing_params=["<modality>"],
            reproducibility="not_reproducible",
            reason=f"未知模态 {modality!r}，无法套用剂量本体；支持：{', '.join(DOSE_SPEC)}",
            falsification="补齐为一个已知模态后重新审计，判定必须随之改变",
        )

    values = _collect_values(report, modality)
    checks: list[FieldCheck] = []
    for f in spec:
        v = values.get(f.key, "")
        if not v:
            status = "MISSING"
        elif any(w in v.lower() for w in RELATIVE_WORDS) and not any(c.isdigit() for c in v):
            status = "RELATIVE"
        elif any(c.isdigit() for c in v) and not _has_unit(v, f):
            status = "NO_UNIT"
        else:
            status = "PRESENT"
        checks.append(FieldCheck(f.key, f.label, f.tier, f.unit, status, v))

    # 分层分数：只计主表模态（附录模态不进主评分体系）
    tier_scores: dict[str, float] = {}
    for tier in sorted(TIER_NAMES):
        group = [c for c in checks if c.tier == tier]
        if group:
            tier_scores[str(tier)] = sum(c.score for c in group) / len(group)
        else:
            tier_scores[str(tier)] = 0.0
    score = sum(TIER_WEIGHTS[t] * tier_scores[str(t)] for t in TIER_NAMES)

    missing = [c.key for c in checks if c.status == "MISSING"]
    no_unit = [c.key for c in checks if c.status == "NO_UNIT"]

    # 可复现性判定用"层内是否仍有 MISSING 字段"，与分数解耦：
    # 分数度量完整性（连续），判定回答能否复现（离散），两者职责不同。
    def tier_complete(tier: int) -> bool:
        group = [c for c in checks if c.tier == tier]
        return bool(group) and all(c.status != "MISSING" for c in group)

    c1, c2, c3 = tier_complete(1), tier_complete(2), tier_complete(3)
    if not c1:
        verdict = "not_reproducible"
        reason = "T1 设备参数有字段完全缺失，无可用输入；等同宣称"
    elif not c2:
        verdict = "underdetermined"
        reason = "T1 齐但 T2 界面条件缺失（施加器几何/距离/组织参数），剂量场方程欠定"
    elif not c3:
        verdict = "reproducible_with_model"
        reason = "T1+T2 完备（边界条件齐），缺 T3 靶区剂量；可交 L2 剂量场引擎求解"
    else:
        verdict = "reproducible_from_report"
        reason = "T1/T2/T3 全齐，报告自带靶区剂量，可直接复现"

    falsification = (
        f"若补齐以下字段，completeness_score 必须上升且判定必须可改变："
        f"{', '.join(missing) if missing else '（必报字段已齐）'}。"
        f"{'；有数值无单位字段：' + ', '.join(no_unit) if no_unit else ''}"
        "；分数只度量报告完整性，不度量疗效或安全性。"
    )
    return DoseAuditReport(
        source_id=rep_id, modality=modality, is_appendix=modality in APPENDIX_MODALITIES,
        fields=checks, completeness_score=score, tier_scores=tier_scores,
        missing_params=missing, reproducibility=verdict, reason=reason,
        falsification=falsification,
    )


def audit_library(csv_path: str | Path) -> list[DoseAuditReport]:
    """批量审计证据库，用于生成字段缺失图谱。"""
    with open(csv_path, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [audit(row, source_id=row.get("id", "")) for row in rows]


def gap_map(reports: list[DoseAuditReport]) -> dict:
    """把批量审计汇总成字段缺失图谱——规范的实证基础。"""
    from collections import Counter, defaultdict
    by_modality: dict[str, list[DoseAuditReport]] = defaultdict(list)
    for r in reports:
        by_modality[r.modality].append(r)

    out: dict = {
        "schema": "dose_gap_map/v1",
        "total_rows": len(reports),
        "modality_count": {m: len(v) for m, v in sorted(by_modality.items())},
        "by_modality": {},
        "field_missing_rate": {},
        "overall": {},
    }
    all_field_miss: Counter = Counter()
    all_field_total: Counter = Counter()
    for modality, reps in sorted(by_modality.items()):
        spec = DOSE_SPEC.get(modality, ())
        rates = []
        for f in spec:
            n_miss = sum(1 for r in reps if any(x.key == f.key and x.status == "MISSING" for x in r.fields))
            rates.append({"key": f.key, "label": f.label, "tier": f.tier,
                          "missing": n_miss, "total": len(reps),
                          "missing_rate": round(n_miss / len(reps), 3)})
            if modality not in APPENDIX_MODALITIES:
                all_field_miss[f.key] += n_miss
                all_field_total[f.key] += len(reps)
        verdicts = Counter(r.reproducibility for r in reps)
        scores = [r.completeness_score for r in reps]
        out["by_modality"][modality] = {
            "n": len(reps),
            "field_missing": rates,
            "verdict_distribution": dict(verdicts),
            "mean_completeness": round(sum(scores) / len(scores), 3) if scores else 0.0,
        }

    out["field_missing_rate"] = [
        {"key": k, "missing": all_field_miss[k], "total": all_field_total[k],
         "missing_rate": round(all_field_miss[k] / all_field_total[k], 3)}
        for k in sorted(all_field_total, key=lambda x: -all_field_miss[x] / all_field_total[x])
    ]
    out["overall"] = {
        "verdict_distribution": dict(Counter(r.reproducibility for r in reports)),
        "mean_completeness": round(sum(r.completeness_score for r in reports) / len(reports), 3),
        "max_completeness": round(max(r.completeness_score for r in reports), 3),
    }
    return out


def json_schema() -> dict:
    """导出机器可读 JSON Schema（draft 2020-12）。"""
    props: dict = {
        "id": {"type": "string", "description": "报告唯一标识"},
        "modality": {"type": "string", "enum": list(DOSE_SPEC)},
        "signal": {"type": "string", "description": "自由文本描述；可被正则抽取但不替代必报字段"},
    }
    for spec in DOSE_SPEC.values():
        for f in spec:
            if f.key not in props:
                entry: dict = {"type": "string", "description": f"{f.label}（T{f.tier}，单位 {f.unit}）"}
                if f.unit.startswith("enum:"):
                    entry["enum"] = f.unit.split("enum:", 1)[1].split("|")
                props[f.key] = entry
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://github.com/lm203688/interface-atlas/blob/main/schema/dose_report_v1.schema.json",
        "title": "Physical Intervention Dose Report Spec v1",
        "description": "物理干预剂量报告规范 v1（社区草案）：四层剂量本体 T1 设备/T2 界面/T3 靶区/T4 个体化",
        "type": "object",
        "required": ["id", "modality"],
        "properties": props,
        "x-tier-weights": TIER_WEIGHTS,
        "x-disclaimer": DISCLAIMER,
    }


def render(rep: DoseAuditReport) -> str:
    lines = [
        f"剂量完整性审计：{rep.source_id or '(未命名)'}  模态={rep.modality}"
        + ("（附录）" if rep.is_appendix else ""),
        f"完整性分数：{rep.completeness_score:.3f}  "
        + "  ".join(f"T{k}={v:.2f}" for k, v in rep.tier_scores.items()),
        f"可复现性：{rep.reproducibility}  理由：{rep.reason}",
        f"缺失字段（{len(rep.missing_params)}）：{', '.join(rep.missing_params) if rep.missing_params else '无'}",
        "",
    ]
    for tier in sorted(TIER_NAMES):
        group = [f for f in rep.fields if f.tier == tier]
        if not group:
            continue
        lines.append(f"  [{TIER_NAMES[tier]}]")
        for f in group:
            lines.append(f"    {f.status:<8} {f.label}({f.key}, {f.unit}) = {f.value or '—'}")
    lines += ["", f"推翻条件：{rep.falsification}", "", DISCLAIMER]
    return "\n".join(lines)


def to_json(rep: DoseAuditReport, *, indent: int = 2) -> str:
    return json.dumps(rep.as_dict(), ensure_ascii=False, indent=indent)
