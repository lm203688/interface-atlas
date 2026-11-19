"""service.py 共享内核的契约测试：三入口（CLI/API/MCP）同源性的地基."""
from __future__ import annotations

import pytest

from interface_atlas import claimcheck, service


def test_meta_carries_version_fingerprint_boundary():
    m = service.meta()
    assert m["library"] == "interface-atlas"
    assert m["version"]
    assert len(m["evidence_fingerprint"]) == 12  # sha256 前 12 位
    assert m["evidence_count"] == 103
    assert "非临床" in m["boundary"] or "不构成医疗建议" in m["boundary"]


def test_claim_check_matches_cli_kernel():
    """服务层判定与内核完全一致（同源性是硬要求）."""
    text = "植入式迷走刺激能治疗类风湿关节炎"
    via_service = service.claim_check(text)["result"]
    via_kernel = claimcheck.check(text)
    for key in ("verdict", "reason", "matched_ids", "counter_ids", "falsification"):
        assert via_service[key] == via_kernel[key]
    assert via_service["verdict"] == "SUPPORTED"


def test_claim_check_verdict_depends_on_input():
    """判定依被测对象而变（防恒真，继承 v0.1 教训）."""
    supported = service.claim_check("植入式迷走刺激能治疗类风湿关节炎")["result"]["verdict"]
    refuted = service.claim_check("红光面板能抗衰老")["result"]["verdict"]
    none_ = service.claim_check("特定频率电磁波修复细胞")["result"]["verdict"]
    assert supported == "SUPPORTED"
    assert refuted == "REFUTED"
    assert none_ == "NO_EVIDENCE"
    assert len({supported, refuted, none_}) == 3


def test_claim_check_rejects_empty():
    with pytest.raises(ValueError):
        service.claim_check("   ")


def test_param_check_verdict_varies_with_input():
    complete = service.param_check({"modality": "photobiomodulation",
                                    "wavelength": "660nm", "power_density": "25mW/cm2",
                                    "fluence": "6J/cm2", "exposure_time": "240s",
                                    "spot_area": "1cm2"})["result"]
    assert complete["verdict"] in ("COMPLETE", "PARTIAL")  # 核心参数齐 + 数值单位
    oow = service.param_check({"modality": "focused ultrasound",
                               "mechanical_index": "2.4", "intensity": "15W/cm2"})["result"]
    assert oow["verdict"] == "OUT_OF_WINDOW"
    none_ = service.param_check({"modality": "electromagnetic_resonance"})["result"]
    assert none_["verdict"] == "NO_EVIDENCE"
    assert len({complete["verdict"], "OUT_OF_WINDOW", "NO_EVIDENCE"}) == 3


def test_dose_audit_four_reproducibility_tiers():
    """同一审计器随字段完备度给出不同档位（可证伪）."""
    only_t1 = service.dose_audit({"modality": "photobiomodulation",
                                  "wavelength": "660nm"})["result"]
    assert only_t1["reproducibility"] in ("not_reproducible", "underdetermined")
    t1_t2 = service.dose_audit({"modality": "photobiomodulation",
                                "wavelength": "660nm", "power_density": "25mW/cm2",
                                "fluence": "6J/cm2", "exposure_time": "240s",
                                "spot_area": "1cm2",
                                "applicator_distance": "2cm", "applicator_angle": "90deg",
                                "beam_profile": "gaussian"})["result"]
    assert t1_t2["reproducibility"] in ("reproducible_with_model", "reproducible_from_report")
    assert t1_t2["completeness_score"] > only_t1["completeness_score"]


def test_dose_field_three_modalities():
    light = service.dose_field("light", {"preset": "muscle_nir", "target_depth_cm": 1.0})["result"]
    assert light  # 非空结果
    us = service.dose_field("ultrasound", {"preset": "tfus_transcranial"})["result"]
    assert us
    es = service.dose_field("electric", {"current_ma": 1.0})["result"]
    assert es
    with pytest.raises(ValueError):
        service.dose_field("magnetism", {})


def test_evidence_query_filters():
    all_rows = service.evidence_query()["result"]
    assert all_rows["count"] == 103
    pbm = service.evidence_query(modality="photobiomodulation")["result"]
    assert 0 < pbm["count"] < 103
    grade_a = service.evidence_query(grade_min="A")["result"]
    assert all(r["evidence_grade"] == "A" for r in grade_a["rows"])


def test_headmodel_license_gate_fail_closed():
    models = service.headmodel_list()["result"]["models"]
    assert len(models) == 4
    nc = next(m for m in models if "NC" in m["license"] or "DUA" in m["license"])
    assert nc["commercial_ok"] is False
    denied = service.headmodel_license(nc["name"], commercial=True)["result"]
    assert denied["allowed"] is False
    with pytest.raises(Exception):
        service.headmodel_license("nonexistent_model")


def test_headmodel_sample_reproducible():
    a = service.headmodel_sample(n=50, seed=42)["result"]
    b = service.headmodel_sample(n=50, seed=42)["result"]
    assert a == b
    assert 2.0 <= a["skull_thickness_mm"]["min"] <= a["skull_thickness_mm"]["max"] <= 15.0


def test_l3_cases_have_falsifiers():
    cases = service.l3_cases()["result"]["cases"]
    assert len(cases) == 4
    assert all(len(c["falsifier"]) >= 20 for c in cases)  # 硬门：推翻条件 ≥20 字
    assert {c["evidence_tier"] for c in cases} >= {"peer-reviewed"}


def test_report_journey_combines_claim_and_evidence():
    rep = service.report("红光面板能抗衰老")["result"]
    assert rep["verdict"]["verdict"] == "REFUTED"
    assert rep["related_evidence"]  # 命中证据随报告返回
    assert rep["claim"] == "红光面板能抗衰老"
