"""L3 剂量-效应公开基准（种子版）：三层断裂 + 多指标 + 可证伪。

设计哲学（为什么不是 leaderboard）
----------------------------------
这个基准**不排名模型**、不给"哪个引擎最准"的单一分数。VCC2025 的教训
（单一 MAE 诱导过拟合）在剂量-效应领域会更致命，因为这里存在一条
**三层因果断裂**，每一层都有独立的失败模式：

    [设备参数] --T1/T2--> [模拟靶区剂量] --①--> [体内实测剂量] --②--> [生物效应]

* 断裂①（模拟 ≠ 体内）：Kop et al. 2026 (bioRxiv, DOI 10.64898/
  2026.09.27.754841) 在 24 名成人用 MR-ARFI 实测组织位移，发现**位移
  预测了 LGN/V1 的 BOLD 效应，而模拟的 in-situ 声强与两者都不相关**。
  即：当前声传播模拟在个体幅值上不可靠，颅骨与声耦合的差异没被建模。
* 断裂②（体内剂量 ≠ 效应）：Spivak et al. 2026 (Neuromodulation, DOI
  10.1016/j.neurom.2026.03.015) 把 tFUS 升到 10.08 W/cm²（14× FDA 诊断
  限）仍未见不良事件，但剂量-效应"mixed and inconclusive"。剂量上去了，
  效应没跟上是真实存在的，不是噪声。
* 暴露度量本身的多维性：Ghaempanah et al. 2025 证明换 CSF 分割管线后
  电场峰值几乎不变（0.236→0.224 V/m）但聚焦体积掉 39%（10.9→6.43 cm³）。
  单一峰值做"剂量"会漏掉最重要的差异。

因此本基准的产物是**一组带推翻条件的公开种子案例 + 多指标评测协议**，
供外部团队用同一口径报数、互相证伪，而不是一个可被刷分的排行榜。

种子案例（全部锚定已发表/预印本一手来源）
------------------------------------------
B1 THIELSCHER-2026-DOSE-HARM：focal tDCS 0.2 V/m 统一场强（43 人 / 8 靶区），
   跨样本泛化；DOI 10.1016/j.neuroimage.2026.121882
B2 KOP-2026-MRARFI：MR-ARFI 位移预测效应、模拟强度不预测；24 人双盲；
   DOI 10.64898/2026.09.27.754841
B3 SPIVAK-2026-ESCALATION：tFUS 剂量递增 0.72–10.08 W/cm² 安全但效应
   inconclusive；DOI 10.1016/j.neurom.2026.03.015
B4 GHAEMPANAH-2025-SEGMENTATION：CSF 分割 → 峰值不变/聚焦性 -39%；
   BeNe Brain Stimulation 2025 poster（DOI 见案例 provenance）

每个案例带 ``falsifier``：什么新证据会推翻这条基准结论。留不出推翻条件的
案例不收录——这是进入本基准的硬门。

MIT License. 案例文字为事实转述，版权归原作者；仅索引 DOI，不复制全文。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Iterable, Sequence

# ---------------------------------------------------------------------------

DISCLAIMER = (
    "L3 种子基准不是临床指南或疗效结论；它索引已发表研究的剂量-效应发现，"
    "用于统一评测口径与可证伪比较。任何条目不得被解读为治疗建议。"
)

THREE_GAPS: tuple[dict[str, str], ...] = (
    {
        "gap": "① 模拟剂量 vs 体内剂量",
        "anchor": "Kop et al. 2026, 24 人 MR-ARFI：组织位移预测 BOLD 效应，"
                  "模拟 in-situ 声强与位移/效应均不相关",
        "implication": "声压/场强点估计必须配体内剂量锚点（MR-ARFI / ECAP / 实测），"
                       "否则剂量-效应结论不可信",
        "doi": "10.64898/2026.09.27.754841",
    },
    {
        "gap": "② 体内剂量 vs 生物效应",
        "anchor": "Spivak et al. 2026：tFUS 升至 14× FDA 诊断限仍安全，"
                  "但剂量-效应 inconclusive（10 人，样本不足）",
        "implication": "剂量上升≠效应上升；非线性/倒 U 型可能存在，"
                       "基准须能容纳非单调剂量-效应，禁止预设线性",
        "doi": "10.1016/j.neurom.2026.03.015",
    },
    {
        "gap": "③ 暴露度量的多维性",
        "anchor": "Ghaempanah et al. 2025：换 CSF 分割管线峰值 0.236→0.224 V/m"
                  "（~5%）而聚焦体积 10.9→6.43 cm³（~39%）",
        "implication": "单一标量剂量不充分；基准强制多指标（峰值+聚焦+覆盖+分布）",
        "doi": "N/A (conference poster, BeNe Brain Stimulation 2025)",
    },
)

METRICS: tuple[str, ...] = (
    "peak_accuracy",        # 峰值与锚点一致度（各模态自己的锚点）
    "focality_agreement",   # 聚焦性排序/分位一致
    "coverage_calibration", # 覆盖区间校准（conformal 命中率）
    "dose_response_shape",  # 剂量-效应形状（单调/非单调）一致性
    "anatomy_sensitivity",  # 个体化解剖敏感性是否落在实测变异范围
)


@dataclass(frozen=True)
class BenchmarkCase:
    """一条种子基准案例。"""
    case_id: str
    modality: str
    title: str
    dose_anchor: str            # 基准剂量值/条件（带单位）
    effect_observation: str     # 实测效应发现（事实转述）
    measurement: str            # 体内剂量/效应的测量手段
    source_doi: str
    source_year: int
    evidence_tier: str          # peer-reviewed | preprint | conference
    gaps: tuple[str, ...]       # 涉及的三层断裂编号
    required_metrics: tuple[str, ...]  # 本案例应报的指标
    falsifier: str              # 什么新证据推翻它（硬门）
    caveats: tuple[str, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "case_id": self.case_id, "modality": self.modality, "title": self.title,
            "dose_anchor": self.dose_anchor, "effect_observation": self.effect_observation,
            "measurement": self.measurement, "source_doi": self.source_doi,
            "source_year": self.source_year, "evidence_tier": self.evidence_tier,
            "gaps": list(self.gaps), "required_metrics": list(self.required_metrics),
            "falsifier": self.falsifier, "caveats": list(self.caveats),
            "provenance": self.provenance, "disclaimer": DISCLAIMER,
        }


# ---------------------------------------------------------------------------
# 种子案例
# ---------------------------------------------------------------------------

CASES: tuple[BenchmarkCase, ...] = (
    BenchmarkCase(
        case_id="B1_THIELSCHER2026_DOSE_HARM",
        modality="transcranial electric",
        title="Focal tDCS 跨靶区统一场强 0.2 V/m",
        dose_anchor="平均电场 0.2 V/m（43 名参与者 × 8 靶区，center-surround montage 半径按靶区优化）",
        effect_observation="用建模框架把既往产生行为/生理效应的场强统一到 0.2 V/m，"
                           "独立样本跨样本泛化验证在群体平均上达到目标场强，同时保留个体差异；"
                           "已用于约 1000 数据集的多中心研究",
        measurement="SimNIBS 电场建模 + 跨样本泛化检验（组平均场强）",
        source_doi="10.1016/j.neuroimage.2026.121882",
        source_year=2026, evidence_tier="peer-reviewed",
        gaps=("③",),
        required_metrics=("peak_accuracy", "anatomy_sensitivity", "coverage_calibration"),
        falsifier="若独立大样本显示 0.2 V/m 在不同靶区并不能统一行为/生理效应量"
                  "（即场强对齐后效应仍系统性随靶区变化），则本基准的'统一场强'前提被推翻",
        caveats=("组平均对齐不等于个体对齐；个体差异被刻意保留",),
        provenance={"journal": "NeuroImage 331:121882",
                    "authors": "Thielscher, Hayek, Puonti, et al.",
                    "open_source": True},
    ),
    BenchmarkCase(
        case_id="B2_KOP2026_MRARFI",
        modality="focused ultrasound",
        title="MR-ARFI 实测位移预测效应、模拟强度不预测",
        dose_anchor="64 元相控阵 500 kHz、50 W/cm² 自由水 ISPPA、145 ms 连续脉冲，靶向 LGN；"
                    "以 MR-ARFI 组织位移（µm）作体内剂量",
        effect_observation="24 名成人双盲交叉：位移预测 LGN 与同侧 V1 的视觉诱发 BOLD 增幅"
                           "（个体间与体素内均成立）；模拟 in-situ 声强与位移及效应均无显著关联",
        measurement="MR-ARFI（运动编码梯度相位）位移 + 同步 TUS-fMRI；d′=−0.05 盲法验证",
        source_doi="10.64898/2026.09.27.754841",
        source_year=2026, evidence_tier="preprint",
        gaps=("①", "②"),
        required_metrics=("dose_response_shape", "anatomy_sensitivity", "focality_agreement"),
        falsifier="若后续预注册多中心重复实验显示模拟 in-situ 声强能稳定预测位移与 BOLD"
                  "（而非仅 MR-ARFI 位移预测），则'模拟剂量幅值不可靠'这条结论被推翻；"
                  "本案例为预印本，尚未同行评议",
        caveats=("单靶点（LGN）；预印本；位移是吸收+组织力学加权量，不是纯声能",),
        provenance={"journal": "bioRxiv 2026.09.27.754841",
                    "authors": "Kop, Mohammadjavadi, Feredoes, et al.",
                    "n_subjects": 24},
    ),
    BenchmarkCase(
        case_id="B3_SPIVAK2026_ESCALATION",
        modality="focused ultrasound",
        title="tFUS 剂量递增至 14× FDA 诊断限：安全但效应未定",
        dose_anchor="ISPTA.3 = 0.72–10.08 W/cm²（4 档），靶向右杏仁核，BX Pulsar 1002",
        effect_observation="10 名健康成人：全剂量档无不良事件、无 MRI 可检结构改变；"
                           "但情绪效价/唤醒的剂量-效应总体不显著、兴奋/抑制范式无差异，"
                           "趋势 mixed and inconclusive",
        measurement="BabelBrain 颅内能量沉积建模 + IAPS 情绪任务 + 结构 MRI 安全监测",
        source_doi="10.1016/j.neurom.2026.03.015",
        source_year=2026, evidence_tier="peer-reviewed",
        gaps=("②",),
        required_metrics=("dose_response_shape", "peak_accuracy", "coverage_calibration"),
        falsifier="若更大样本（预设效能）在这些剂量档显示清晰的单调或倒 U 型剂量-效应，"
                  "则'当前证据 inconclusive'被更新；本案例不支持'越高越好'，也不证明无效",
        caveats=("n=10；单次刺激后效应衰减；时间维度未充分刻画",),
        provenance={"journal": "Neuromodulation (2026), DOI 10.1016/j.neurom.2026.03.015",
                    "preprint": "10.1101/2025.10.21.25338018", "n_subjects": 10},
    ),
    BenchmarkCase(
        case_id="B4_GHAEMPANAH2025_SEGMENTATION",
        modality="transcranial electric",
        title="CSF 分割管线改变聚焦性 39% 而峰值近乎不变",
        dose_anchor="1 mA 额部 montage（anode F3 / cathode F4，5×5 cm 电极，4 mm 厚），47 人 T1",
        effect_observation="CHARM 默认 CSF 掩膜 vs 商用分割：E99 峰值 0.236→0.224 V/m"
                           "（中位 Δ=−0.004 V/m，p=0.071，不显著）；"
                           "但 ≥75% 峰值体积 10.9→6.43 cm³（中位 Δ=−4.16 cm³，p=3.6e−14，r=0.87）",
        measurement="SimNIBS 4.1；Wilcoxon 符号秩检验；E99 与 ≥75% 峰值皮层体积",
        source_doi="N/A (BeNe Brain Stimulation 2025 poster; VUB research portal)",
        source_year=2025, evidence_tier="conference",
        gaps=("③",),
        required_metrics=("focality_agreement", "peak_accuracy", "anatomy_sensitivity"),
        falsifier="若独立队列显示 CSF 分割差异对聚焦体积影响 <10% 且不显著，"
                  "则'分割管线主导聚焦性'被推翻；本条为会议摘要，未同行评议全文",
        caveats=("会议 poster；单一 montage；聚焦性指标为皮层体积阈值法",),
        provenance={"conference": "BeNe Brain Stimulation, Maastricht 2025-11-13",
                    "authors": "Ghaempanah, Van Schependom, Nagels",
                    "effect_size_r": 0.87},
    ),
)


def validate_case(case: BenchmarkCase) -> tuple[bool, list[str]]:
    """基准准入硬门：无可证伪推翻条件 / 无来源 / 指标不合法 → 拒绝收录。"""
    errs: list[str] = []
    if not case.falsifier or len(case.falsifier.strip()) < 20:
        errs.append("falsifier 缺失或过短（不可证伪的案例不收录）")
    if (not case.source_doi or
            (case.source_doi.startswith("N/A") and case.evidence_tier == "peer-reviewed")):
        errs.append("同行评议案例必须有 DOI；会议/预印本可显式标 N/A 并注明场合")
    if not case.measurement.strip():
        errs.append("measurement 缺失")
    bad = [m for m in case.required_metrics if m not in METRICS]
    if bad:
        errs.append(f"非法指标名：{bad}")
    if not case.gaps:
        errs.append("必须标注涉及的三层断裂编号")
    return (not errs), errs


# ---------------------------------------------------------------------------
# 多指标评测：拒绝单一 MAE 聚合
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MetricScore:
    metric: str
    value: float
    unit: str
    method: str
    reliability: str   # analytic | empirical | judgement

    def as_dict(self) -> dict:
        return {"metric": self.metric, "value": round(self.value, 6),
                "unit": self.unit, "method": self.method,
                "reliability": self.reliability}


@dataclass(frozen=True)
class CaseScore:
    case_id: str
    scores: tuple[MetricScore, ...]
    passed: bool | None    # None = 证据不足，不下结论（诚实留标）
    note: str

    def as_dict(self) -> dict:
        return {"case_id": self.case_id,
                "scores": [s.as_dict() for s in self.scores],
                "passed": self.passed, "note": self.note}


def evaluate_case(case: BenchmarkCase,
                  observed: dict[str, float]) -> CaseScore:
    """按案例要求的多指标评测一个提交结果。

    ``observed`` 是提交方按各指标口径报来的数值（键为 METRICS 名）。
    缺任一必需指标 → ``passed=None``（证据不足，不聚合、不猜），
    这正是与"单一 MAE 强制出分"的 leaderboard 的根本区别。
    """
    scores: list[MetricScore] = []
    missing: list[str] = []
    for m in case.required_metrics:
        if m not in observed:
            missing.append(m)
            continue
        val = float(observed[m])
        unit = {
            "peak_accuracy": "relative_error",
            "focality_agreement": "rank_correlation",
            "coverage_calibration": "empirical_coverage",
            "dose_response_shape": "shape_match_score",
            "anatomy_sensitivity": "max_min_ratio",
        }[m]
        scores.append(MetricScore(m, val, unit, "submitter-reported", "empirical"))
    if missing:
        return CaseScore(case.case_id, tuple(scores), None,
                         f"缺必需指标 {missing}：拒绝出分（证据不足，不聚合为单一分数）")
    # 多指标门槛：不做加权求和，逐项过阈值（避免 VCC 单一分数覆辙）
    thresholds = {
        "peak_accuracy": (0.0, 0.30),          # 相对误差 ≤30%
        "focality_agreement": (0.5, 1.0),      # 秩相关 ≥0.5
        "coverage_calibration": (0.85, 1.0),   # 经验覆盖 ≥85%
        "dose_response_shape": (0.6, 1.0),     # 形状一致 ≥0.6
        "anatomy_sensitivity": (0.0, 4.0),     # 个体变异比 ≤4（Berger ~2.1 的宽松上界）
    }
    failed = []
    for s in scores:
        lo, hi = thresholds[s.metric]
        if not (lo <= s.value <= hi):
            failed.append(f"{s.metric}={s.value:.3f}∉[{lo},{hi}]")
    passed = not failed
    note = "全部必需指标逐项过阈（无加权总分）" if passed else "未过阈：" + ", ".join(failed)
    return CaseScore(case.case_id, tuple(scores), passed, note)


def run_benchmark(submissions: dict[str, dict[str, float]]) -> dict[str, Any]:
    """批量评测。返回逐案例多指标分数 + 显式拒绝聚合总分。"""
    per_case: list[dict[str, Any]] = []
    for case in CASES:
        sub = submissions.get(case.case_id, {})
        cs = evaluate_case(case, sub)
        per_case.append({"case": case.as_dict(), "score": cs.as_dict()})
    n_pass = sum(1 for p in per_case if p["score"]["passed"] is True)
    n_fail = sum(1 for p in per_case if p["score"]["passed"] is False)
    n_na = sum(1 for p in per_case if p["score"]["passed"] is None)
    return {
        "three_gaps": list(THREE_GAPS),
        "metrics": list(METRICS),
        "cases": per_case,
        "summary": {"passed": n_pass, "failed": n_fail, "insufficient_evidence": n_na},
        "aggregate_score": None,
        "aggregate_refused_reason": "拒绝单一聚合分：三层断裂 + VCC2025 单一 MAE 教训；"
                                    "必须逐案例逐指标解读（见 L3_BENCHMARK.md）",
        "disclaimer": DISCLAIMER,
    }


def build_report() -> dict[str, Any]:
    """生成基准自检报告：案例准入、三层断裂、评测协议完整性。"""
    admissions = []
    for case in CASES:
        ok, errs = validate_case(case)
        admissions.append({"case_id": case.case_id, "admitted": ok, "errors": errs})
    # 协议自检：缺一指标就应 passed=None
    probe = evaluate_case(CASES[0], {"peak_accuracy": 0.1})  # 只给 1/3 指标
    # 全给但一项越界
    probe2 = evaluate_case(CASES[0], {
        "peak_accuracy": 0.1, "anatomy_sensitivity": 2.0, "coverage_calibration": 0.5,
    })
    gates = {
        "all_cases_admitted": all(a["admitted"] for a in admissions),
        "insufficient_refused": probe.passed is None,
        "out_of_range_failed": probe2.passed is False,
        "aggregate_refused": run_benchmark({})["aggregate_score"] is None,
        "three_gaps_complete": len(THREE_GAPS) == 3 and all(
            g.get("anchor") and g.get("implication") and g.get("doi") for g in THREE_GAPS
        ),
        "every_case_has_falsifier": all(bool(c.falsifier.strip()) for c in CASES),
        "every_case_has_doi_or_explicit_na": all(bool(c.source_doi) for c in CASES),
    }
    return {
        "admissions": admissions,
        "gates": gates,
        "overall_passed": all(gates.values()),
        "n_cases": len(CASES),
        "metrics": list(METRICS),
        "three_gaps": list(THREE_GAPS),
        "disclaimer": DISCLAIMER,
    }


# ---------------------------------------------------------------------------

def render_report(rep: dict[str, Any]) -> str:
    lines = [
        "L3 剂量-效应种子基准 — 自检报告",
        "=" * 46,
        f"案例数 {rep['n_cases']}   整体 {'PASS' if rep['overall_passed'] else 'FAIL'}",
        "三层断裂：",
    ]
    for g in rep["three_gaps"]:
        lines.append(f"  {g['gap']}")
        lines.append(f"    锚点：{g['anchor']}")
    lines.append("评测协议门禁：")
    for k, v in rep["gates"].items():
        lines.append(f"  [{'OK' if v else 'BAD'}] {k}")
    lines.append("案例准入：")
    for a in rep["admissions"]:
        lines.append(f"  [{'OK' if a['admitted'] else 'BAD'}] {a['case_id']}")
    lines.append(f"  {DISCLAIMER}")
    return "\n".join(lines)


def render_case(case: BenchmarkCase) -> str:
    lines = [
        f"[{case.case_id}] {case.title}",
        f"  模态={case.modality}  来源={case.source_doi} ({case.source_year}, {case.evidence_tier})",
        f"  剂量锚点：{case.dose_anchor}",
        f"  效应发现：{case.effect_observation}",
        f"  测量手段：{case.measurement}",
        f"  涉及断裂：{', '.join(case.gaps)}",
        f"  应报指标：{', '.join(case.required_metrics)}",
        f"  ✗ 推翻条件：{case.falsifier}",
    ]
    for c in case.caveats:
        lines.append(f"  ⚠ {c}")
    lines.append(f"  {DISCLAIMER}")
    return "\n".join(lines)
