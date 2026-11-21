"""L1 规范可引用条文注册表（外部引用的稳定锚点）.

为什么存在：
    DOSE_REPORTING_SPEC.md 是给人读的规范，但没有稳定条文编号，外部论文只能
    写"按 interface-atlas 规范报告波长"——无法精确到条、无法版本冻结、无法
    机器核对。本模块把 dosimetry.DOSE_SPEC（唯一权威字段源）编译成带稳定 ID
    的条文注册表 schema/l1_clauses_v1.json，供外部逐条引用。

设计约束（诚实边界）：
    1. 唯一权威源是 DOSE_SPEC——注册表由它生成，绝不手工维护第二份字段清单。
    2. ID 用语义 key（如 IA-L1-PBM-T1-wavelength）而非序号——插入新字段不会
       让已发布 ID 漂移，这是"可引用"的前提。
    3. evidence_refs 是该模态的实证背景文献，只说明"字段为何必要"；不代表
       该字段被这些文献验证，更不构成疗效或安全性证据。此限定写进注册表。
    4. 条文文本是规范性报告要求（"须报什么/怎么报"），不含任何疗效声明。
    5. 社区草案，非 ICRU / FDA / ISO 官方标准。
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from .dosimetry import (
    APPENDIX_MODALITIES,
    DOSE_SPEC,
    TIER_NAMES,
    TIER_WEIGHTS,
    DoseField,
    canonical_modality,
)

SCHEMA_NAME = "l1_clauses/v1"
SPEC_VERSION = "1.0.0"
PINNED_REF = "v1.3.0"

_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = _ROOT / "schema" / "l1_clauses_v1.json"
EVIDENCE_PATH = _ROOT / "data" / "evidence_matrix.csv"

# 模态短码：与 DOSE_SPEC 插入顺序一致；短码一经发布不得更改。
MODALITY_CODES = {
    "photobiomodulation": "PBM",
    "vagus nerve stimulation": "VNS",
    "focused ultrasound": "FUS",
    "mechanical loading": "MECH",
    "metabolic / microbiota": "MET",
}

_TIER_EN = {1: "device", 2: "interface", 3: "target", 4: "subject"}

_RULES = (
    ("IA-L1-RULE-001", "not_reproducible",
     "T1 有必报字段完全缺失：无可用剂量学输入，等同宣称，不可复现。"),
    ("IA-L1-RULE-002", "underdetermined",
     "T1 齐、T2 缺：剂量场方程边界条件不完整，不可外推到其他个体。"),
    ("IA-L1-RULE-003", "reproducible_with_model",
     "T1+T2 齐、T3 缺：边界条件完备，可由剂量场引擎求解（L2 的定位）。"),
    ("IA-L1-RULE-004", "reproducible_from_report",
     "T1/T2/T3 全齐：报告自带靶区剂量，可直接复现。"),
    ("IA-L1-RULE-005", "completeness_score",
     "完整性分数 = Σ(层均分 × 层权重)；仅度量报告完整性，不度量疗效或安全性。"),
)

_SOURCE_KIND_RANK = {"pmc": 0, "doi": 1, "institutional": 2}
_REF_LIMIT = 5  # 每模态最多挂 5 条真实出处，足够引用且不把注册表撑大


def _source_kind(url: str) -> str:
    u = url.lower()
    if "pmc.ncbi.nlm.nih.gov" in u or "pubmed.ncbi.nlm.nih.gov" in u:
        return "pmc"
    if "doi.org" in u or u.startswith("10."):
        return "doi"
    return "institutional"


def evidence_refs(modality: str, *, limit: int = _REF_LIMIT,
                  evidence_path: Path = EVIDENCE_PATH) -> list[dict]:
    """从真实证据库取该模态的实证背景出处；不编造 DOI/URL。

    空文件/缺列时返回空列表（注册表仍能生成，但出处为空——不伪造）。
    排序按 (出处强度, URL) 确定性进行，保证注册表逐字节可复现。
    """
    if not evidence_path.exists():
        return []
    refs: list[dict] = []
    with evidence_path.open(encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            if canonical_modality(row.get("modality", "")) != modality:
                continue
            url = (row.get("source_url") or "").strip()
            if not url:
                continue
            refs.append({
                "evidence_id": row.get("id", ""),
                "url": url,
                "kind": _source_kind(url),
            })
    seen: set[str] = set()
    unique: list[dict] = []
    for ref in refs:
        if ref["url"] in seen:
            continue
        seen.add(ref["url"])
        unique.append(ref)
    unique.sort(key=lambda r: (_SOURCE_KIND_RANK[r["kind"]], r["url"]))
    return unique[:limit]


def clause_definition(field: DoseField, modality: str) -> str:
    """把权威字段定义编译成规范性条文文本（无疗效声明）。"""
    appendix = modality in APPENDIX_MODALITIES
    noun = {
        1: "暴露量" if appendix else "设备参数",
        2: "暴露条件" if appendix else "界面条件",
        3: "靶部位暴露读数" if appendix else "靶区剂量",
        4: "个体化解剖",
    }[field.tier]
    rule = {
        1: ("须带可核对单位与暴露物身份。" if appendix
            else "须为设备面板或说明书可直接读出的量，并带可核对单位。"),
        2: ("须描述暴露施加条件与材料/菌株信息。" if appendix
            else "须描述从设备到靶区的传递条件；该项是跨研究比较的关键，也是文献最常漏报项。"),
        3: ("须为实测读数，并标注检测方法。" if appendix
            else "须为测量值，或明确标注计算来源的估计值；模拟量不得冒充体内实测。"),
        4: "须区分群体平均、体模与个体影像，不得省略。",
    }[field.tier]
    text = (f"报告须给出{field.label}（{field.unit}），"
            f"作为{noun}（T{field.tier}）。{rule}")
    if field.note:
        text += f" 说明：{field.note}"
    return text


def build_registry(*, evidence_path: Path = EVIDENCE_PATH) -> dict:
    """从 DOSE_SPEC 编译条文注册表。字段增删自动同步，不产生第二真源。"""
    clauses: list[dict] = []
    by_modality: dict[str, int] = {}
    for modality, code in MODALITY_CODES.items():
        refs = evidence_refs(modality, evidence_path=evidence_path)
        for field in DOSE_SPEC[modality]:
            clauses.append({
                "id": f"IA-L1-{code}-T{field.tier}-{field.key}",
                "modality": modality,
                "modality_code": code,
                "tier": field.tier,
                "tier_en": _TIER_EN[field.tier],
                "tier_zh": TIER_NAMES[field.tier],
                "key": field.key,
                "label": field.label,
                "unit": field.unit,
                "aliases": list(field.aliases),
                "definition": clause_definition(field, modality),
                "note": field.note or None,
                "appendix_only": modality in APPENDIX_MODALITIES,
                "evidence_refs": refs,
            })
        by_modality[modality] = len(DOSE_SPEC[modality])

    rules = []
    for rid, outcome, text in _RULES:
        item = {"id": rid, "outcome": outcome, "definition": text}
        if outcome == "completeness_score":
            item["tier_weights"] = {f"T{k}": v for k, v in TIER_WEIGHTS.items()}
        rules.append(item)

    return {
        "schema": SCHEMA_NAME,
        "spec_version": SPEC_VERSION,
        "title": "物理干预剂量报告规范 v1 — 可引用条文注册表",
        "status": "community draft v1; not an ICRU / FDA / ISO official standard",
        "source_of_truth": "interface_atlas.dosimetry.DOSE_SPEC",
        "evidence_ref_scope": (
            "出处是该模态的实证背景文献，用于说明字段为何必要；"
            "不代表该字段已被这些文献验证，不构成疗效或安全性证据。"
        ),
        "citation": {
            "version_pin": PINNED_REF,
            "rule": "引用单条须写条文 ID 与规范版本；引用整份须写版本与仓库 URL。",
            "preferred": (
                "interface-atlas contributors. Physical Intervention Dose Reporting "
                "Spec v1 (interface-atlas L1), clause <ID>, "
                f"{PINNED_REF}. https://github.com/lm203688/interface-atlas"
            ),
        },
        "tiers": [
            {"id": f"T{k}", "en": _TIER_EN[k], "zh": TIER_NAMES[k],
             "weight": TIER_WEIGHTS[k]}
            for k in (1, 2, 3, 4)
        ],
        "rules": rules,
        "clauses": clauses,
        "counts": {
            "clauses": len(clauses),
            "rules": len(rules),
            "by_modality": by_modality,
        },
    }


def render(registry: dict) -> str:
    """确定性渲染：键排序 + 固定缩进，保证 --check 逐字节可复现。"""
    return json.dumps(registry, ensure_ascii=False, indent=2, sort_keys=True) + "\n"


def save_registry(path: Path = REGISTRY_PATH, *,
                  evidence_path: Path = EVIDENCE_PATH) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render(build_registry(evidence_path=evidence_path)), encoding="utf-8")
    return path


def check_registry(path: Path = REGISTRY_PATH) -> tuple[bool, str]:
    """磁盘注册表是否与权威源编译结果一致。"""
    if not path.exists():
        return False, f"注册表缺失：{path}"
    expected = render(build_registry())
    actual = path.read_text(encoding="utf-8")
    if actual != expected:
        return False, "注册表与 DOSE_SPEC 编译结果不一致，请运行 l1clauses --save 重建"
    return True, f"注册表一致：{len(build_registry()['clauses'])} 条条文 + {len(_RULES)} 条判定规则"


def fingerprint(path: Path = REGISTRY_PATH) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_clause(clause_id: str, *, path: Path = REGISTRY_PATH) -> dict | None:
    data = json.loads(path.read_text(encoding="utf-8"))
    target = clause_id.strip()
    for item in data["clauses"] + data["rules"]:
        if item["id"].lower() == target.lower():
            return item
    return None


def format_clause(item: dict) -> str:
    """单条条文的人类可读引用卡。"""
    lines = [f"[{item['id']}] {item['definition']}"]
    if "modality" in item:
        lines.append(
            f"    模态={item['modality']}  层=T{item['tier']} {item['tier_en']}  "
            f"字段={item['key']}  单位={item['unit']}"
        )
        refs = item.get("evidence_refs") or []
        lines.append(f"    实证背景出处={len(refs)} 条"
                     + ("：" + ", ".join(r["url"] for r in refs[:3]) if refs else "（证据库暂无此模态出处）"))
    if "tier_weights" in item:
        lines.append("    层权重=" + ", ".join(f"{k}={v}" for k, v in item["tier_weights"].items()))
    lines.append("    引用：" + "interface-atlas L1 " + item["id"] + f" ({PINNED_REF})")
    return "\n".join(lines)
