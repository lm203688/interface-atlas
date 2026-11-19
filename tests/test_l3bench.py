"""L3 基准的可证伪测试：准入硬门、拒分机制、三层断裂锚点、多指标协议。"""
import json

import pytest

from interface_atlas import l3bench as bench
from interface_atlas import l3_validation as lv


def test_all_cases_admitted():
    assert len(bench.CASES) >= 4
    for c in bench.CASES:
        ok, errs = bench.validate_case(c)
        assert ok, (c.case_id, errs)


def test_every_case_has_falsifier():
    """无可证伪推翻条件的案例不得进入基准——这是硬门。"""
    for c in bench.CASES:
        assert len(c.falsifier.strip()) >= 20, c.case_id
        assert "推翻" in c.falsifier or "更新" in c.falsifier, c.case_id


def test_case_without_falsifier_rejected():
    bad = bench.BenchmarkCase(
        case_id="BAD", modality="x", title="bad", dose_anchor="x",
        effect_observation="x", measurement="", source_doi="10.1/x",
        source_year=2020, evidence_tier="preprint", gaps=("①",),
        required_metrics=("peak_accuracy",), falsifier="短的",  # 太短 + measurement 空
    )
    ok, errs = bench.validate_case(bad)
    assert ok is False
    assert any("falsifier" in e for e in errs)


def test_insufficient_evidence_refused():
    """缺必需指标必须 passed=None（不猜、不聚合）。"""
    cs = bench.evaluate_case(bench.CASES[0], {"peak_accuracy": 0.1})
    assert cs.passed is None
    assert "缺必需指标" in cs.note


def test_threshold_gate():
    sub = {"peak_accuracy": 0.1, "anatomy_sensitivity": 2.0, "coverage_calibration": 0.9}
    assert bench.evaluate_case(bench.CASES[0], sub).passed is True
    sub["coverage_calibration"] = 0.5  # 越界
    cs = bench.evaluate_case(bench.CASES[0], sub)
    assert cs.passed is False
    assert "coverage_calibration" in cs.note


def test_aggregate_score_always_none():
    """拒绝单一聚合分——VCC2025 单一 MAE 教训。"""
    subs = {c.case_id: {"peak_accuracy": 0.1, "anatomy_sensitivity": 2.0,
                        "coverage_calibration": 0.9, "focality_agreement": 0.8,
                        "dose_response_shape": 0.9}
            for c in bench.CASES}
    r = bench.run_benchmark(subs)
    assert r["aggregate_score"] is None
    assert r["aggregate_refused_reason"]
    assert "summary" in r


def test_three_gaps_anchored():
    gaps = {g["gap"][:2]: g for g in bench.THREE_GAPS}
    assert len(gaps) == 3
    assert all(g["anchor"] and g["implication"] and g["doi"] for g in gaps.values())


def test_kop_and_spivak_dois_present():
    """三层断裂的一手锚点必须是这两篇，不编造。"""
    all_text = json.dumps(bench.THREE_GAPS) + json.dumps([c.as_dict() for c in bench.CASES])
    assert "10.64898/2026.09.27.754841" in all_text  # Kop 2026 MR-ARFI
    assert "10.1016/j.neurom.2026.03.015" in all_text  # Spivak 2026
    assert "10.1016/j.neuroimage.2026.121882" in all_text  # Thielscher 2026


def test_validation_report_all_pass():
    rep = lv.build_report()
    assert rep["overall_passed"] is True
    for name, c in rep["checks"].items():
        assert c["passed"] is True, (name, c)
