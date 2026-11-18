"""剂量报告规范 v1 测试：重点是判定必须依被测报告而变（禁止恒真）。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from interface_atlas import dosimetry, evidence


FUS_FULL = {
    "id": "FUS-TEST", "modality": "focused ultrasound", "signal": "tFUS",
    "frequency": "0.5 MHz", "mechanical_index": "MI 1.2", "intensity": "100 mW/cm2",
    "pulse_duration": "10 ms", "duty_cycle": "50%", "temperature_rise": "0.5 degC",
    "transducer_geometry": "circular 10mm", "skull_attenuation": "3 dB",
    "target_pressure": "0.4 MPa", "target_temperature": "37.2 degC",
    "tissue_acoustics": "0.5 dB/cm/MHz", "anatomy_source": "subject_imaging",
}


def test_verdict_changes_with_report_stages():
    """核心可证伪性：逐层补全，判定必须阶梯式改变，分数单调递增。"""
    stages = [
        ({}, "not_reproducible"),
        ({k: FUS_FULL[k] for k in ("frequency", "mechanical_index", "intensity",
                                    "pulse_duration", "duty_cycle", "temperature_rise")},
         "underdetermined"),
        ({"transducer_geometry": "circular 10mm", "skull_attenuation": "3 dB"},
         "reproducible_with_model"),
        ({"target_pressure": "0.4 MPa", "target_temperature": "37.2 degC",
          "tissue_acoustics": "0.5 dB/cm/MHz"}, "reproducible_from_report"),
        ({"anatomy_source": "subject_imaging"}, "reproducible_from_report"),
    ]
    report: dict = {"id": "X", "modality": "focused ultrasound"}
    last_score = -1.0
    for patch, expected in stages:
        report.update(patch)
        rep = dosimetry.audit(report, source_id="X")
        assert rep.reproducibility == expected, (rep.reproducibility, expected, rep.reason)
        if expected == "reproducible_from_report" and last_score == pytest.approx(0.9):
            assert rep.completeness_score >= last_score, "T4 补全后分数不应下降"
        else:
            assert rep.completeness_score > last_score, "分数必须随补全单调上升"
        last_score = rep.completeness_score
    assert last_score == pytest.approx(1.0)


def test_relative_dose_is_not_present():
    rep = dosimetry.audit({"id": "X", "modality": "photobiomodulation", "fluence": "low dose"})
    fluence = next(f for f in rep.fields if f.key == "fluence")
    assert fluence.status == "RELATIVE"


def test_unitless_number_is_flagged():
    rep = dosimetry.audit({"id": "X", "modality": "photobiomodulation", "wavelength": "660"})
    w = next(f for f in rep.fields if f.key == "wavelength")
    assert w.status in ("NO_UNIT", "MISSING"), w.status


def test_signal_regex_extracts_device_params():
    rep = dosimetry.audit({
        "id": "X", "modality": "photobiomodulation",
        "signal": "660 nm red, 25 mW/cm2, 6 J/cm2, 240 s",
    })
    keys = {f.key for f in rep.fields if f.status != "MISSING"}
    assert {"wavelength", "power_density", "fluence", "exposure_time"} <= keys


def test_unknown_modality_fails_closed():
    rep = dosimetry.audit({"id": "X", "modality": "quantum healing"})
    assert rep.reproducibility == "not_reproducible"
    assert rep.completeness_score == 0.0
    assert "<modality>" in rep.missing_params


def test_appendix_modality_does_not_claim_reproducible():
    rep = dosimetry.audit({"id": "X", "modality": "metabolic / bile acid"})
    assert rep.is_appendix is True
    assert rep.reproducibility in ("not_reproducible", "underdetermined",
                                   "reproducible_with_model", "reproducible_from_report")
    assert rep.modality == "metabolic / microbiota"  # 归一化到附录主类


def test_library_gap_map_is_reproducible():
    reports = dosimetry.audit_library(evidence.EVIDENCE_CSV)
    gap = dosimetry.gap_map(reports)
    assert gap["total_rows"] == len(evidence.load())
    assert gap["schema"] == "dose_gap_map/v1"
    assert 0.0 <= gap["overall"]["mean_completeness"] <= 1.0
    # 两次生成必须逐字段一致
    assert gap == dosimetry.gap_map(dosimetry.audit_library(evidence.EVIDENCE_CSV))


def test_shipped_gap_map_matches_current_data():
    p = Path(__file__).resolve().parents[1] / "data" / "dose_gap_map.json"
    if not p.exists():
        pytest.skip("gap map 未生成")
    shipped = json.loads(p.read_text(encoding="utf-8"))
    current = dosimetry.gap_map(dosimetry.audit_library(evidence.EVIDENCE_CSV))
    assert shipped == current, "随库分发的缺失图谱必须与当前数据一致"


def test_json_schema_valid_shape():
    s = dosimetry.json_schema()
    assert s["$schema"].endswith("2020-12/schema")
    assert "photobiomodulation" in s["properties"]["modality"]["enum"]
    assert s["x-tier-weights"] == dosimetry.TIER_WEIGHTS
    assert len(s["properties"]) >= 40


def test_disclaimer_present_in_all_outputs():
    rep = dosimetry.audit({"id": "X", "modality": "photobiomodulation"})
    assert "社区草案" in dosimetry.DISCLAIMER
    assert dosimetry.DISCLAIMER in json.dumps(rep.as_dict(), ensure_ascii=False)
    assert dosimetry.DISCLAIMER in dosimetry.render(rep)
