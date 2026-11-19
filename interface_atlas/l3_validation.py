"""L2 个体化层 + L3 基准的独立验证。

可证伪精神：每项检验都设计成"如果实现错了就会 FAIL"，而不是复述函数返回值。
* 几何守恒用独立求和（不复用 layer_thicknesses 的内部逻辑）。
* 许可门禁用构造的 CC BY-NC / DUA 数据，验证 fail-closed。
* 群体抽样验证确定性（同 seed）与非退化（不同个体确实不同、变异在合理范围）。
* 特征抽取验证"缺项不聚合""逐项可信度标注""跨模态不归一化"。
* 基准协议验证拒分机制、三层断裂锚点带 DOI、每案有推翻条件。

MIT License.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import anatomy as an
from . import dose_features as df
from . import electric_field as ef
from . import l3bench as bench
from . import ultrasound_field as uf
from . import light_field as lf


# ---------------------------------------------------------------------------
# 个体化解剖层
# ---------------------------------------------------------------------------

def geometry_conservation() -> dict[str, Any]:
    """几何恒等式：各层厚之和必须等于最外层半径；层半径严格递增。"""
    rows = []
    ok = True
    for name, model in an.BUILTIN_MODELS.items():
        radii = [lay.outer_radius_cm for lay in model.layers]
        # 独立重算厚度：内层段 = 相邻半径差，最内层 = 首半径
        thickness_independent = [radii[0]] + [b - a for a, b in zip(radii, radii[1:])]
        s = sum(thickness_independent)
        strict = all(b > a for a, b in zip(radii, radii[1:]))
        match = abs(s - model.outer_radius_cm) < 1e-9
        row = {"model": name, "sum_thickness": round(s, 6),
               "outer_radius": round(model.outer_radius_cm, 6),
               "strictly_increasing": strict, "conserved": match}
        rows.append(row)
        ok = ok and strict and match
    return {"passed": bool(ok), "rows": rows,
            "note": "厚度和==外半径是层模型的几何恒等式；失败说明层半径定义错乱"}


def license_gate() -> dict[str, Any]:
    """许可门禁必须 fail-closed：NC/DUA 数据商用被拒，MIT 通过。"""
    mit_ok, _ = an.check_license(an.get_model("spherical"), commercial=True)
    nc_ok, nc_reason = an.check_license(an.get_model("ernie_extended_metadata"), commercial=True)
    dua_ok, dua_reason = an.check_license(an.get_model("hcp_population"), commercial=True)
    # 非商用一律放行（研究用途）
    nc_research_ok, _ = an.check_license(an.get_model("ernie_extended_metadata"), commercial=False)
    passed = bool(mit_ok and not nc_ok and not dua_ok and nc_research_ok)
    return {"passed": passed,
            "mit_commercial_ok": bool(mit_ok),
            "nc_commercial_rejected": bool(not nc_ok), "nc_reason": nc_reason,
            "dua_commercial_rejected": bool(not dua_ok), "dua_reason": dua_reason,
            "nc_research_allowed": bool(nc_research_ok),
            "note": "CC BY-NC 与 HCP DUA 不得进商用管线；MIT 代理体模可商用",
            "anchors": {k: an.SOURCES[k]["license"] for k in
                        ("analytic", "siminibs_ernie_extended", "hcp_population")}}


def population_sampling() -> dict[str, Any]:
    """群体抽样：同 seed 可复现、不同个体非退化、颅骨厚度在生理合理区间。"""
    a = an.population_sample(n=30, seed=7)
    b = an.population_sample(n=30, seed=7)
    c = an.population_sample(n=30, seed=8)
    ra = [m.skull_thickness_mm() for m in a]
    rb = [m.skull_thickness_mm() for m in b]
    rc = [m.skull_thickness_mm() for m in c]
    reproducible = ra == rb
    seed_sensitive = ra != rc
    nondegenerate = max(ra) - min(ra) > 1.0  # 至少 1mm 跨度
    physiological = all(2.0 < x < 15.0 for x in ra)
    passed = bool(reproducible and seed_sensitive and nondegenerate and physiological)
    return {"passed": passed, "n": len(a),
            "reproducible_same_seed": reproducible,
            "different_seed_differs": seed_sensitive,
            "skull_mm_min": round(min(ra), 3), "skull_mm_max": round(max(ra), 3),
            "spread_mm": round(max(ra) - min(ra), 3),
            "all_physiological_2_15mm": physiological,
            "illustrative": all(m.illustrative for m in a),
            "note": "分布为 illustrative 示例值；真实队列用 load_population_stats 替换"}


# ---------------------------------------------------------------------------
# 特征层
# ---------------------------------------------------------------------------

def feature_multimetric() -> dict[str, Any]:
    """特征向量必须多指标齐出，且每个指标带可信度来源；缺一即拒。"""
    lr = lf.dose_field(lf.TissueOptics(0.3, 70.0), max_depth_cm=2.0, target_depth_cm=1.0)
    vl = df.from_light_field(lr)
    tr = uf.TRANSDUCER_PRESETS["tfus_transcranial"]
    vu = df.from_ultrasound(uf.dose_field(tr, (uf.TISSUE_PRESETS["skull_thin"],)), transducer=tr)
    ve = df.from_electric(ef.dose_field(1.0))
    d = vu.as_dict()
    # 多指标：关键字段都存在
    needed = ("peak", "penetration_depth_cm", "decay_per_cm", "fwhm_cm",
              "focality_proxy", "coverage_depths_cm", "reliability")
    complete = all(k in d for k in needed)
    # 每个指标都有 reliability 标注
    rel_keys = set(vu.reliability)
    rel_complete = {"peak", "focality_proxy", "decay_per_cm"}.issubset(rel_keys)
    # 光模态聚焦性必须诚实标 unavailable（P1 无横向场）
    light_focality_honest = vl.reliability["focality_proxy"] == "unavailable"
    # 跨模态峰值单位不同，不归一化：三者 peak_unit 必须各不相同
    units = {vl.peak_unit, vu.peak_unit, ve.peak_unit}
    no_fake_normalization = len(units) == 3
    passed = bool(complete and rel_complete and light_focality_honest and no_fake_normalization)
    return {"passed": passed, "ultrasound_multimetric_complete": complete,
            "reliability_labeled": rel_complete,
            "light_focality_unavailable": light_focality_honest,
            "distinct_peak_units": sorted(units), "no_cross_modality_normalization": no_fake_normalization,
            "note": "单一标量会漏掉 Ghaempanah 2025 的 39% 聚焦性差异"}


def numeric_field_focality() -> dict[str, Any]:
    """数值场聚焦性：已知几何（中心 1 个体素=峰值）应给出可解析核对的体积占比。"""
    import numpy as np
    grid = np.zeros((11, 11, 11))
    grid[5, 5, 5] = 1.0              # 单个体素峰值
    vox = 0.5                         # 0.5cm 体素 → 总体积 11^3*0.125=166.375 cm3
    v = df.from_numeric_field(modality="focused ultrasound", subject_id="synthetic",
                              peak=1.0, peak_unit="norm", field_grid=grid, voxel_size_cm=vox)
    expected_cov = 1 * vox ** 3        # 仅峰值体素 ≥0.5·peak
    total = 11 ** 3 * vox ** 3
    expected_foc = expected_cov / total
    foc_ok = abs(v.focality_proxy - expected_foc) < 1e-9
    cov_ok = abs(v.coverage_depths_cm[0.5] - expected_cov) < 1e-9
    return {"passed": bool(foc_ok and cov_ok),
            "expected_focality": round(expected_foc, 9),
            "measured_focality": round(v.focality_proxy, 9),
            "expected_coverage_cm3": round(expected_cov, 6),
            "measured_coverage_cm3": round(v.coverage_depths_cm[0.5], 6),
            "note": "用单峰值体素的解析体积占比交叉核对数值体素求和（不同实现路径）"}


# ---------------------------------------------------------------------------
# 基准协议
# ---------------------------------------------------------------------------

def benchmark_protocol() -> dict[str, Any]:
    """基准协议：证据不足拒分、越界拒过、无聚合总分、每案有推翻条件与来源。"""
    # 空提交 → 全部 insufficient，无聚合分
    empty = bench.run_benchmark({})
    all_na = all(c["score"]["passed"] is None for c in empty["cases"])
    no_aggregate = empty["aggregate_score"] is None
    # 逐项过阈 / 越界
    good = bench.evaluate_case(bench.CASES[0], {
        "peak_accuracy": 0.1, "anatomy_sensitivity": 2.0, "coverage_calibration": 0.9})
    bad = bench.evaluate_case(bench.CASES[0], {
        "peak_accuracy": 0.1, "anatomy_sensitivity": 2.0, "coverage_calibration": 0.5})
    gate_ok = good.passed is True and bad.passed is False
    # 每案准入 + 推翻条件
    admitted = [bench.validate_case(c)[0] for c in bench.CASES]
    all_admitted = all(admitted) and len(bench.CASES) >= 4
    passed = bool(all_na and no_aggregate and gate_ok and all_admitted)
    return {"passed": passed, "empty_all_insufficient": all_na,
            "aggregate_refused": no_aggregate, "threshold_gate_ok": gate_ok,
            "all_cases_admitted": all_admitted, "n_cases": len(bench.CASES),
            "required_metrics": list(bench.METRICS),
            "note": "拒绝单一聚合分；缺指标=证据不足不猜；逐指标过阈不做加权求和"}


def three_gaps_anchored() -> dict[str, Any]:
    """三层断裂每条都必须有实证锚点与 DOI/出处，杜绝空泛断言。"""
    gaps = bench.THREE_GAPS
    has_anchor = all(g.get("anchor", "").strip() for g in gaps)
    has_implication = all(g.get("implication", "").strip() for g in gaps)
    has_doi = all(g.get("doi", "").strip() for g in gaps)
    key_dois = {"10.64898/2026.09.27.754841", "10.1016/j.neurom.2026.03.015"}
    dois = {g["doi"] for g in gaps}
    has_kop_spivak = key_dois.issubset(dois)
    passed = bool(len(gaps) == 3 and has_anchor and has_implication and has_doi and has_kop_spivak)
    return {"passed": passed, "n_gaps": len(gaps), "all_anchored": has_anchor,
            "all_with_implication": has_implication, "all_with_doi": has_doi,
            "kop_and_spivak_present": has_kop_spivak,
            "note": "锚点必须是一手研究（Kop 2026 MR-ARFI / Spivak 2026 剂量递增），不编造"}


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------

def build_report() -> dict[str, Any]:
    checks = {
        "geometry_conservation": geometry_conservation(),
        "license_gate": license_gate(),
        "population_sampling": population_sampling(),
        "feature_multimetric": feature_multimetric(),
        "numeric_field_focality": numeric_field_focality(),
        "benchmark_protocol": benchmark_protocol(),
        "three_gaps_anchored": three_gaps_anchored(),
    }
    gates = {
        "l3_benchmark_selfcheck": bench.build_report()["overall_passed"],
    }
    overall = all(c["passed"] for c in checks.values()) and all(gates.values())
    return {"checks": checks, "gates": gates, "overall_passed": bool(overall),
            "disclaimer": bench.DISCLAIMER}


def save_report(path: str | Path) -> str:
    rep = build_report()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rep, ensure_ascii=False, indent=2), encoding="utf-8")
    return str(p)


def render_report(rep: dict[str, Any]) -> str:
    lines = ["L2 个体化 + L3 基准 — 独立验证报告", "=" * 46,
             f"整体：{'PASS' if rep['overall_passed'] else 'FAIL'}", ""]
    for name, c in rep["checks"].items():
        lines.append(f"[{'PASS' if c['passed'] else 'FAIL'}] {name}")
        if not c["passed"]:
            lines.append(f"       {c}")
    for name, ok in rep["gates"].items():
        lines.append(f"[{'PASS' if ok else 'FAIL'}] {name}")
    lines.append(f"\n{rep['disclaimer']}")
    return "\n".join(lines)
