"""A2 物理干预证据与剂量图谱.

Schema v1.0.0：每条证据写成可机读记录，强制区分「机制」「人体终点」「剂量窗」「副作用」。
设计原则（继承 vec-trust-kit 的教训）：
  - 不接受无剂量、无来源、无证据等级的"有效"宣称；
  - 缺少 DOI/URL 的记录不删除，而是显式标记为 unresolved 并在审计中报出；
  - 证据等级 A/B/C/D 的含义见 README，不允许自定义等级。

用法：
  from interface_atlas.evidence import load, validate, query, audit
  rows = load()
  report = audit()
"""
from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
EVIDENCE_CSV = DATA_DIR / "evidence_matrix.csv"
SCHEMA_VERSION = "1.1.0"

FIELDS = (
    "id", "modality", "signal", "receiver_gene", "receiver", "cell_type",
    "endpoint", "effect_direction", "effect_size", "dose", "dose_unit",
    "model_system", "human_evidence", "evidence_grade", "adverse_event",
    "safety_boundary", "source_url", "notes",
)
REQUIRED = (
    "id", "modality", "receiver_gene", "cell_type", "endpoint",
    "effect_direction", "evidence_grade", "human_evidence",
)
GRADES = {"A", "B", "C", "D"}
DIRECTIONS = {"increase", "decrease", "biphasic", "mixed", "none", "physical_constraint"}
HUMAN = {"yes", "no", "mixed"}

# 种子证据：全部来自《非药物接口层深挖报告_2026-11-16.md》及其核查来源。
# source_url 为空者表示前轮调研记录但未保留可点击链接，audit() 会明确列出，不做静默填补。
SEED: list[dict[str, str]] = [
    dict(id="PBM-001", modality="photobiomodulation", signal="660 nm red", receiver_gene="COX4I1", receiver="CCO",
         cell_type="fibroblasts", endpoint="collagen synthesis", effect_direction="increase", effect_size="modest",
         dose="1-10", dose_unit="J/cm2", model_system="human skin", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="biphasic: high fluence inhibitory", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/",
         notes="umbrella review 纳入 meta 分析；效应量多为小到中等"),
    dict(id="PBM-002", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="keratinocytes", endpoint="proliferation", effect_direction="increase", effect_size="",
         dose="lower than fibroblasts", dose_unit="relative", model_system="in vitro", human_evidence="no", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC11630336/",
         notes="成纤维细胞所需剂量显著高于角质形成细胞，是细胞阈值差异的例子"),
    dict(id="PBM-003", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="macrophages", endpoint="M2 reparative polarization", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="preclinical + clinical", human_evidence="mixed", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://www.jkslms.or.kr/view.html?uid=395&vmd=Full",
         notes="PBM 免疫调节综述：巨噬/中性粒/DC/Treg 均被报道"),
    dict(id="PBM-004", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="Treg", endpoint="Treg expansion", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="preclinical", human_evidence="no", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://www.jkslms.or.kr/view.html?uid=395&vmd=Full", notes=""),
    dict(id="PBM-005", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="knee osteoarthritis", endpoint="disability (WOMAC)", effect_direction="increase", effect_size="moderate certainty",
         dose="", dose_unit="", model_system="human RCT meta", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/", notes="GRADE 中等确定性五项之一"),
    dict(id="PBM-006", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="fibromyalgia patients", endpoint="fatigue severity", effect_direction="decrease", effect_size="large SMD",
         dose="", dose_unit="", model_system="human RCT meta", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/", notes="GRADE 中等确定性五项之一"),
    dict(id="PBM-007", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="hair follicle", endpoint="hair density in androgenetic alopecia", effect_direction="increase", effect_size="moderate certainty",
         dose="", dose_unit="", model_system="human RCT meta", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/", notes=""),
    dict(id="PBM-008", modality="photobiomodulation", signal="transcranial NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="aging brain", endpoint="cognitive function", effect_direction="increase", effect_size="moderate certainty",
         dose="", dose_unit="", model_system="human RCT meta", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/",
         notes="唯一擦边抗衰的人体终点，但仍非死亡率/健康寿命"),
    dict(id="PBM-009", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="35 clinical endpoints", endpoint="23/35 endpoints no significant effect", effect_direction="none", effect_size="0 high-certainty",
         dose="", dose_unit="", model_system="umbrella review", human_evidence="yes", evidence_grade="A",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/",
         notes="最关键的反证：15 meta、>9000 人，无任何高确定性证据"),
    dict(id="PBM-010", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="cortical neurons", endpoint="ATP production", effect_direction="biphasic", effect_size="peak at 3 J/cm2, suppressed at 30 J/cm2",
         dose="25 mW/cm2 / 3 vs 30", dose_unit="J/cm2", model_system="in vitro", human_evidence="no", evidence_grade="C",
         adverse_event="mitochondrial impairment at high fluence", safety_boundary="Arndt-Schulz biphasic curve",
         source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC11671386/", notes="剂量窗是具体数字，不是'越多越好'"),
    dict(id="PBM-011", modality="photobiomodulation", signal="800 nm transcranial", receiver_gene="COX4I1", receiver="CCO",
         cell_type="human skull", endpoint="light transmittance", effect_direction="physical_constraint", effect_size="4.18-4.24% (mouse 40.1%)",
         dose="800", dose_unit="nm", model_system="ex vivo comparison", human_evidence="no", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC11671386/",
         notes="解释为何动物剂量不能直接外推到人"),
    dict(id="PBM-012", modality="photobiomodulation", signal="red/NIR", receiver_gene="COX4I1", receiver="CCO",
         cell_type="skin / retina", endpoint="phototoxicity", effect_direction="increase", effect_size="",
         dose="0.01-10", dose_unit="J/cm2 ocular", model_system="clinical guidance", human_evidence="yes", evidence_grade="B",
         adverse_event="retinal stress, transient photophobia", safety_boundary="eye protection; photosensitizing drugs; active tumor caution",
         source_url="https://eyewiki.org/w/index.php?oldid=124006", notes=""),
    dict(id="VNS-001", modality="vagus nerve stimulation", signal="transcutaneous auricular", receiver_gene="CHRNA7", receiver="alpha7 nAChR",
         cell_type="macrophages", endpoint="TNF-alpha / IL-6", effect_direction="decrease", effect_size="",
         dose="", dose_unit="", model_system="animal + non-CVD human", human_evidence="mixed", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC13598969/",
         notes="人体直接证据与动物推测必须分开；综述已明确标注"),
    dict(id="VNS-002", modality="vagus nerve stimulation", signal="taVNS", receiver_gene="CHRNA7", receiver="alpha7 nAChR",
         cell_type="human serum", endpoint="CRP", effect_direction="decrease", effect_size="",
         dose="", dose_unit="", model_system="human", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC13598969/", notes=""),
    dict(id="VNS-003", modality="vagus nerve stimulation", signal="implanted VNS", receiver_gene="CHRNA7", receiver="alpha7 nAChR",
         cell_type="rheumatoid arthritis patients", endpoint="approved indication", effect_direction="increase", effect_size="RESET-RA 242 例假刺激对照 RCT，3 月 ACR20 达主要终点",
         dose="1 分钟/日（植入式）", dose_unit="implanted, once daily", model_system="regulatory", human_evidence="yes", evidence_grade="A",
         adverse_event="hoarseness, cough, infection; 器械相关严重不良事件 1.7%", safety_boundary="implant-only; 禁用于迷走神经切除/脾切除/起搏器患者", source_url="https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpma/pma.cfm?id=P240039",
         notes="SetPoint System，FDA PMA P240039，2025-07-30 批准；临床登记 NCT04539964。一手来源 2026-11-17 核实"),
    dict(id="VNS-004", modality="vagus nerve stimulation", signal="VNS meta", receiver_gene="CHRNA7", receiver="alpha7 nAChR",
         cell_type="human inflammatory markers", endpoint="cytokine reduction", effect_direction="none", effect_size="数值下降但未达显著：TNF-alpha g=-0.21 p=0.359；IL-6 g=-0.94 p=0.112",
         dose="", dose_unit="", model_system="36 项系统综述 / 26 项纳入 meta", human_evidence="yes", evidence_grade="A",
         adverse_event="", safety_boundary="仅 4 项急性炎症长期研究亚组 CRP 显著", source_url="https://doi.org/10.1016/j.bbi.2023.12.008",
         notes="Schiweck et al. 2024, Brain Behav Immun 116:237-258；PMID 38070618。预注册 PRISMA meta，1135 人。抗炎抗衰宣称的关键反证。一手 DOI 2026-11-17 核实"),
    dict(id="VNS-005", modality="vagus nerve stimulation", signal="consumer auricular stimulation", receiver_gene="CHRNA7", receiver="alpha7 nAChR",
         cell_type="human autonomic", endpoint="HRV change", effect_direction="none", effect_size="g=0.014, BF01=24.678（强支持零假设）；正式发表版 16 项单盲研究",
         dose="", dose_unit="", model_system="16 项单盲假刺激对照 meta", human_evidence="yes", evidence_grade="A",
         adverse_event="", safety_boundary="", source_url="https://doi.org/10.1111/psyp.13933",
         notes="Wolf et al. 2021, Psychophysiology 58(11):e13933。更正：前轮记为『17 项』实为 bioRxiv 预印本数字，正式发表版为 16 项；结论不变（taVNS 不改变 HRV）。消费级耳夹可能根本未激活迷走神经。一手 DOI 2026-11-17 核实"),
    dict(id="FUS-001", modality="focused ultrasound", signal="low-intensity tFUS", receiver_gene="PIEZO1", receiver="mechanosensitive channels",
         cell_type="neurons", endpoint="excitation or inhibition", effect_direction="biphasic", effect_size="parameter-dependent",
         dose="ISPTA<=720 mW/cm2; ISPPA<=190 W/cm2; MI<=1.9", dose_unit="FDA diagnostic limits",
         model_system="in vitro + animal + human", human_evidence="mixed", evidence_grade="B",
         adverse_event="headache, fatigue (~3.4% mild transient in 677 participants)", safety_boundary="MI<=1.9; thermal rise<=2 C",
         source_url="https://pubmed.ncbi.nlm.nih.gov/41099683/", notes="机制未定：离子通道/膜电容/神经递质"),
    dict(id="FUS-002", modality="focused ultrasound", signal="tFUS safety", receiver_gene="PIEZO1", receiver="mechanosensitive channels",
         cell_type="brain", endpoint="microhemorrhage case", effect_direction="increase", effect_size="1 case, MI>1.9",
         dose="MI>1.9", dose_unit="mechanical index", model_system="clinical case", human_evidence="yes", evidence_grade="B",
         adverse_event="nucleus accumbens microhemorrhage and edema", safety_boundary="超出 FDA 限值即属 significant risk",
         source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC13000543/", notes="安全边界不是建议而是硬阈值"),
    dict(id="FUS-003", modality="focused ultrasound", signal="0.25-0.65 MHz transcranial", receiver_gene="PIEZO1", receiver="mechanosensitive channels",
         cell_type="human skull", endpoint="acoustic attenuation", effect_direction="physical_constraint", effect_size="65-90%",
         dose="0.25-0.65", dose_unit="MHz", model_system="human skull measurement", human_evidence="no", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://rd.springer.com/article/10.1186/s12984-025-01753-2",
         notes="动物剂量外推到人的第二道物理屏障"),
    dict(id="MECH-001", modality="mechanical loading", signal="fluid shear stress", receiver_gene="PIEZO1", receiver="PIEZO1",
         cell_type="osteocytes / BMSC", endpoint="osteoblast differentiation", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="animal + human genetics", human_evidence="mixed", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://doi.org/10.3389/fendo.2025.1658967",
         notes="PIEZO1 随年龄下降，导致骨力学响应减退"),
    dict(id="MECH-002", modality="mechanical loading", signal="osmotic / shear", receiver_gene="TRPV4", receiver="TRPV4",
         cell_type="chondrocytes", endpoint="anabolic signaling", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="in vitro / animal", human_evidence="no", evidence_grade="C",
         adverse_event="", safety_boundary="", source_url="https://www.sciencedirect.com/science/article/abs/pii/S1063458425012786", notes=""),
    dict(id="MECH-003", modality="mechanical loading", signal="overload", receiver_gene="PIEZO1", receiver="PIEZO1",
         cell_type="chondrocytes", endpoint="cartilage degeneration / apoptosis", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="gain-of-function model", human_evidence="no", evidence_grade="C",
         adverse_event="spontaneous cartilage degeneration", safety_boundary="mechanical overload harmful",
         source_url="https://www.sciencedirect.com/science/article/abs/pii/S1063458425012786", notes="机械信号同样有剂量窗"),
    dict(id="MET-001", modality="metabolic / microbiota", signal="dietary fiber -> SCFA", receiver_gene="FFAR2", receiver="GPR43",
         cell_type="intestinal L cells", endpoint="GLP-1 secretion", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="human metabolism", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12337796/",
         notes="最接近药物式分子特异的非药物路线，但效应远弱于 GLP-1RA"),
    dict(id="MET-002", modality="metabolic / microbiota", signal="SCFA", receiver_gene="FFAR3", receiver="GPR41",
         cell_type="intestinal L cells", endpoint="GLP-1 / PYY secretion", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="human metabolism", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12337796/", notes=""),
    dict(id="MET-003", modality="metabolic / bile acid", signal="secondary bile acids", receiver_gene="GPBAR1", receiver="TGR5",
         cell_type="intestinal L cells", endpoint="GLP-1 secretion", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="human metabolism", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12337796/", notes=""),
    dict(id="MET-004", modality="metabolic / bile acid", signal="bile acids", receiver_gene="NR1H4", receiver="FXR",
         cell_type="liver / intestine", endpoint="FGF15/19 production", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="human metabolism", human_evidence="yes", evidence_grade="B",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12337796/", notes=""),
    dict(id="MET-005", modality="metabolic / incretin", signal="endogenous GLP-1", receiver_gene="GLP1R", receiver="GLP-1R",
         cell_type="pancreatic beta cells / brain", endpoint="insulin secretion / satiety", effect_direction="increase", effect_size="",
         dose="", dose_unit="", model_system="human physiology", human_evidence="yes", evidence_grade="A",
         adverse_event="", safety_boundary="", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12337796/",
         notes="药物 GLP-1RA 的受体靶点；非药物路线只能提高内源分泌，强度不可等同"),
]


def _seed_rows() -> list[dict[str, str]]:
    return [{k: str(r.get(k, "")) for k in FIELDS} for r in SEED]


def save_seed(path: Path | str = EVIDENCE_CSV) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(_seed_rows())
    return p


def load(path: Path | str = EVIDENCE_CSV) -> list[dict[str, str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"未找到 {p}。请先运行 python -m interface_atlas.cli evidence --init")
    with open(p, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def validate_row(row: dict[str, str]) -> list[str]:
    errs = []
    for key in REQUIRED:
        if not str(row.get(key, "")).strip():
            errs.append(f"缺少必填字段 {key} (id={row.get('id', '?')})")
    if row.get("evidence_grade") and row["evidence_grade"] not in GRADES:
        errs.append(f"非法证据等级 {row['evidence_grade']} (id={row.get('id')})")
    if row.get("effect_direction") and row["effect_direction"] not in DIRECTIONS:
        errs.append(f"非法效应方向 {row['effect_direction']} (id={row.get('id')})")
    if row.get("human_evidence") and row["human_evidence"] not in HUMAN:
        errs.append(f"非法 human_evidence {row['human_evidence']} (id={row.get('id')})")
    return errs


def query(rows: list[dict[str, str]] | None = None, *, modality: str | None = None,
          gene: str | None = None, cell_type: str | None = None,
          grade_min: str | None = None, human_only: bool = False,
          keyword: str | None = None) -> list[dict[str, str]]:
    rows = rows if rows is not None else load()
    order = {"A": 0, "B": 1, "C": 2, "D": 3}
    out = []
    for r in rows:
        if modality and modality.lower() not in r["modality"].lower():
            continue
        if gene and gene.upper() != r["receiver_gene"].upper():
            continue
        if cell_type and cell_type.lower() not in r["cell_type"].lower():
            continue
        if grade_min and order.get(r["evidence_grade"], 9) > order[grade_min]:
            continue
        if human_only and r["human_evidence"] != "yes":
            continue
        if keyword and keyword.lower() not in " ".join(r.values()).lower():
            continue
        out.append(r)
    return out


def fingerprint(rows: list[dict[str, str]] | None = None) -> str:
    """证据库内容指纹（SHA256 前 12 位），用于 link 输出与跨版本可复现。"""
    rows = rows if rows is not None else load()
    canonical = "\n".join(
        "|".join(r.get(k, "") for k in FIELDS)
        for r in sorted(rows, key=lambda x: x.get("id", ""))
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12]


def _dose_band(row: dict[str, str]) -> str:
    """剂量分层：数值剂量 / 相对剂量 / 物理参数 / 缺失。

    区分「真的没报剂量」与「报的是相对关系或物理屏障参数」，避免把
    "lower than fibroblasts" 或 "800 nm 穿透率" 伪装成可比较的剂量。
    """
    dose = row.get("dose", "").strip()
    unit = row.get("dose_unit", "").strip()
    if not dose:
        return "missing"
    if unit.lower() in {"relative", ""} and any(
        t in dose.lower() for t in ("lower than", "higher than", "relative", "vs", "/")
    ):
        return "relative"
    if unit.lower() in {"nm", "mhz", "mechanical index", "fda diagnostic limits"}:
        return "physical_parameter"
    if unit:
        return "numeric"
    return "relative"


def audit(rows: list[dict[str, str]] | None = None) -> dict:
    rows = rows if rows is not None else load()
    errors: list[str] = []
    warnings: list[str] = []
    for r in rows:
        errors.extend(validate_row(r))
    ids = [r["id"] for r in rows]
    dupes = [k for k, v in Counter(ids).items() if v > 1]
    if dupes:
        errors.append(f"重复 id：{dupes}")

    no_source = [r["id"] for r in rows if not r["source_url"].strip()]
    # v1.1 退出门：A 级记录缺一手来源 = 硬失败。B/C/D 缺源仅警告。
    grade_a_no_source = [r["id"] for r in rows
                         if r["evidence_grade"] == "A" and not r["source_url"].strip()]
    for rid in grade_a_no_source:
        errors.append(f"A 级证据缺一手来源：{rid}（不得用于对外结论）")

    bands: dict[str, list[str]] = {"numeric": [], "relative": [],
                                   "physical_parameter": [], "missing": []}
    for r in rows:
        bands[_dose_band(r)].append(r["id"])
    no_dose = bands["missing"]

    unresolved = [r["id"] for r in rows if "待补" in r["notes"] or "unresolved" in r["notes"].lower()]
    # 人体阳性宣称（increase/decrease/biphasic）缺数值剂量：警告，不阻断。
    # 反证记录（none/physical_constraint）与监管记录本就不适用剂量窗，不计入。
    dose_gap = [r["id"] for r in rows
                if r["human_evidence"] == "yes" and r["evidence_grade"] in {"A", "B"}
                and r["effect_direction"] in {"increase", "decrease", "biphasic"}
                and _dose_band(r) in {"missing", "relative"}]
    if dose_gap:
        warnings.append(
            f"{len(dose_gap)} 条 A/B 级人体阳性宣称缺数值剂量（来源本身未报告，禁止编造）：{dose_gap}"
        )

    return {
        "status": "PASS" if not errors else "FAIL",
        "schema_version": SCHEMA_VERSION,
        "fingerprint": fingerprint(rows),
        "total_rows": len(rows),
        "grade_distribution": dict(Counter(r["evidence_grade"] for r in rows)),
        "human_evidence": dict(Counter(r["human_evidence"] for r in rows)),
        "missing_source_url": no_source,
        "grade_a_missing_source": grade_a_no_source,
        "dose_bands": {k: len(v) for k, v in bands.items()},
        "dose_band_detail": bands,
        "missing_dose": no_dose,
        "dose_warnings": warnings,
        "unresolved_notes": unresolved,
        "errors": errors,
        "caveat": "A 级缺一手来源为硬失败；人体宣称缺剂量为警告（来源未报告，禁止编造）。"
                  "nCPM 与剂量均不可跨研究直接比较；本库是假设生成器，不是疗效证据。",
    }
