"""L2 电场剂量场的可证伪测试：球谐求解、物理方向、线性、逆设计、验证门。"""
import math

import pytest

from interface_atlas import electric_field as ef, es_validation as ev


# ---------------- 球谐求解器 / 单位 ----------------

def test_single_sphere_matches_closed_form():
    """单层均匀球：球谐数值解必须逐点吻合独立解析闭式（扩散区）。"""
    r = ev.single_sphere_vs_closed_form()
    assert r["passed"] is True
    assert r["worst_relative_error"] < 1e-3


def test_center_field_magnitude():
    """球心场解析值 ~0.18 V/m（1mA, σ0.33, R9cm）——文献 0.2V/m 量级。"""
    e = ef.uniform_sphere_closed_form(0.01, 0.33, 9.0, 1.0)
    assert 0.15 < e < 0.25


def test_vm_per_mm_unit_identity():
    """1 V/m = 1 mV/mm（神经元兴奋性常用单位），必须单位恒等。"""
    r = ef.dose_field(1.0, target_depth_cm=2.2)
    assert r.field_mv_mm == pytest.approx(r.field_v_m, rel=1e-12)


# ---------------- 物理方向 / 线性 ----------------

def test_current_linearity():
    """E∝I（可证伪）：电流翻倍场翻倍。"""
    r = ev.current_linearity()
    assert r["passed"] is True
    assert r["field_ratio_2x"] == pytest.approx(2.0, rel=1e-9)


def test_conductivity_scales_inversely():
    """所有 σ×k → E×1/k（固定电流源+准静态；可证伪）。"""
    r = ev.conductivity_scale_invariance(k=2.0)
    assert r["passed"] is True
    assert r["measured_ratio"] == pytest.approx(0.5, rel=1e-9)


def test_skull_conductivity_monotone():
    """颅骨电导率↑ → 靶区场↑（颅骨阻抗降低），方向必须正确。"""
    r = ev.skull_monotonicity()
    assert r["passed"] is True
    assert r["monotone_increasing"] is True
    assert r["field_max_min_ratio"] > 1.5


def test_four_sphere_dose_magnitude():
    """四球 1mA 靶区场落在文献量级 0.1-0.5 V/m。"""
    r = ef.dose_field(1.0, target_depth_cm=2.2)
    assert 0.1 <= r.field_v_m <= 0.5
    assert r.target_layer == "brain"


def test_target_outside_brain_warns():
    """靶区落在脑层外（过浅）→ 必须 WARN 点电极近场不可信。"""
    r = ef.dose_field(1.0, target_depth_cm=0.5)  # r=9.5cm，在颅骨/头皮层
    assert r.target_layer != "brain"
    assert any("点电极" in w or "脑内" in w for w in r.warnings)


# ---------------- 逆设计 ----------------

def test_inverse_closed_form_roundtrip():
    """逆设计闭式解正向复核残差 ~0。"""
    r = ev.inverse_roundtrip()
    assert r["passed"] is True
    assert r["worst_residual_rel"] < 1e-9


def test_inverse_scales_linearly_with_field():
    """目标场翻倍 → 所需电流翻倍（线性，可证伪）。"""
    a = ef.inverse_current(0.2)
    b = ef.inverse_current(0.4)
    assert b["required_current_ma"] == pytest.approx(2 * a["required_current_ma"], rel=1e-9)


def test_inverse_flags_over_practice_window():
    """目标场需 >2mA → 标 within=False（advisory，非硬阈值）。"""
    inv = ef.inverse_current(1.0)  # 需 ~5mA
    assert inv["within_practice_window"] is False
    assert inv["required_current_ma"] > ef.PRACTICE_CURRENT_MAX_MA


# ---------------- 实践窗口 / advisory ----------------

def test_advisory_labeling():
    """2mA 标为 advisory，3mA 显式告警，不冒充 FDA 硬阈值。"""
    r = ev.advisory_labeling()
    assert r["passed"] is True
    assert r["advisory_labeled"] and r["over_window_flagged"]


def test_high_current_triggers_warning():
    """3mA → 越实践窗口告警（可证伪）。"""
    r = ef.dose_field(3.0, target_depth_cm=2.2)
    assert any("实践窗口" in w or "研究性" in w for w in r.warnings)


def test_no_brain_current_fraction_fabrication():
    """不输出入脑电流比例（同心球无法定义，25% 属真实头模 CSF 分流）。"""
    r = ef.dose_field(1.0, target_depth_cm=2.2)
    assert "brain_current_fraction" not in r.as_dict()


# ---------------- 验证全套 ----------------

def test_validation_report_all_gates_pass():
    """整套验证判据门全绿。"""
    rep = ev.build_report()
    assert rep["overall_passed"] is True
    for k in rep["gated_checks"]:
        assert rep["checks"][k].get("passed") is True, k


def test_near_field_truncation_shrinks_with_harmonics():
    """近场截断误差随阶数增加而减小（如实报告的 report 项）。"""
    r = ev.near_field_truncation()
    assert r["error_shrinks_with_more_harmonics"] is True
    errs = [abs(x["relative_error"]) for x in r["rows"]]
    assert errs[0] > errs[-1]
