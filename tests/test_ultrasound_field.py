"""L2 超声剂量场的可证伪测试：物理方向、单位、阈值、逆设计、独立验证。"""
import math

import pytest

from interface_atlas import paramcheck, ultrasound_field as uf, us_validation as uv


TR = uf.TRANSDUCER_PRESETS["tfus_transcranial"]


# ---------------- 基础声学量 / 单位 ----------------

def test_wavelength_units():
    """λ=c/f 单位核对：1MHz 水中应 0.148cm，0.5MHz 脑中应 0.312cm。"""
    assert uf.wavelength_cm(1.0, 1480.0) == pytest.approx(0.148, rel=1e-9)
    assert uf.wavelength_cm(0.5, 1560.0) == pytest.approx(0.312, rel=1e-9)


def test_focus_gain_monotone_and_f_number():
    """聚焦增益随孔径增大、随焦距减小而增大；F# 定义正确。"""
    lam = uf.wavelength_cm(0.5, 1560.0)
    small = uf.Transducer(1.0, 8.0, 0.5, 1.0)
    large = uf.Transducer(3.0, 8.0, 0.5, 1.0)
    long_f = uf.Transducer(3.0, 12.0, 0.5, 1.0)
    assert uf.focus_gain(large, lam) > uf.focus_gain(small, lam)
    assert uf.focus_gain(large, lam) > uf.focus_gain(long_f, lam)
    assert TR.f_number == pytest.approx(8.0 / (2 * 3.0))


def test_intensity_pressure_roundtrip():
    """声强-声压正逆变换闭合，且 1MPa 量级正确（水 ~34 W/cm²？核对）。"""
    rho, c = 1000.0, 1480.0
    p = uf.intensity_to_pressure_pa(10.0, rho, c)
    assert uf.pressure_to_intensity_w_cm2(p, rho, c) == pytest.approx(10.0, rel=1e-12)
    # 10 W/cm² → 峰值声压 ~17.2 MPa（I=p²/2ρc，p=√(2·1000·1480·1e5)）
    assert p / 1e6 == pytest.approx(math.sqrt(2 * rho * c * 10 * 1e4) / 1e6, rel=1e-9)


def test_mi_definition_and_frequency_scaling():
    """MI=p[MPa]/√f[MHz]：同声压下频率×4 → MI 减半。"""
    assert uf.mechanical_index(1.0, 1.0) == pytest.approx(1.0)
    assert uf.mechanical_index(1.0, 4.0) == pytest.approx(0.5)


def test_layer_attenuation_frequency_scaling():
    """α(f)=α(f0)·(f/f0)^ν，且分层衰减可加。"""
    lay = uf.TISSUE_PRESETS["skull_thin"]
    assert lay.attenuation_at(2.0) == pytest.approx(2 * lay.attenuation_at(1.0))
    trans, db = uf.layer_transmission((lay, lay), 0.5)
    # db = 层数 × α(f) × 层厚；skull_thin 厚 0.4cm
    assert db == pytest.approx(2 * lay.attenuation_at(0.5) * lay.thickness_cm)
    assert 0 < trans < 1


# ---------------- 正向剂量场 ----------------

def test_dose_field_dose_scales_linearly_with_surface_intensity():
    """剂量场对表面声强线性（可证伪）：翻倍 → 焦点声强翻倍、声压 √2 倍。"""
    base = uf.dose_field(TR, (uf.TISSUE_PRESETS["skull_thin"],))
    tr2 = uf.Transducer(TR.radius_cm, TR.focal_length_cm, TR.frequency_mhz,
                        TR.intensity_w_cm2 * 2, TR.efficiency)
    dbl = uf.dose_field(tr2, (uf.TISSUE_PRESETS["skull_thin"],))
    assert dbl.isppa_w_cm2 == pytest.approx(2 * base.isppa_w_cm2, rel=1e-9)
    assert dbl.peak_pressure_pa == pytest.approx(math.sqrt(2) * base.peak_pressure_pa, rel=1e-9)


def test_thicker_skull_lowers_focus_dose():
    """颅骨越厚 → 焦点剂量越低（可证伪，方向必须正确）。"""
    thin = uf.dose_field(TR, (uf.TISSUE_PRESETS["skull_thin"],))
    thick = uf.dose_field(TR, (uf.TISSUE_PRESETS["skull_thick"],))
    assert thick.peak_pressure_mpa < thin.peak_pressure_mpa
    assert thick.attenuation_db_total > thin.attenuation_db_total
    assert thick.transmission < thin.transmission


def test_dose_field_reports_thresholds_and_warns_on_high_intensity():
    """高声强越 MI 硬阈值 → 必须报 exceeded + WARN，不静默给'有效'剂量。"""
    hi = uf.Transducer(3.0, 8.0, 0.5, 5.0, 0.7)
    r = uf.dose_field(hi, (uf.TISSUE_PRESETS["skull_thin"],))
    assert r.mi_exceeded is True
    assert r.mi > uf.MI_LIMIT
    assert any("非线性" in w or "MI" in w for w in r.warnings)
    normal = uf.dose_field(TR, (uf.TISSUE_PRESETS["skull_thin"],))
    assert normal.mi_exceeded is False


def test_dose_field_result_is_finite_and_shaped():
    r = uf.dose_field(TR, (uf.TISSUE_PRESETS["skull_thin"],))
    d = r.as_dict()
    assert all(math.isfinite(v) for k, v in d.items()
               if isinstance(v, (int, float)) and not isinstance(v, bool))
    assert d["aperture_half_angle_deg"] > 0
    assert d["lateral_focus_mm"] > 0 and d["axial_focus_mm"] > 0
    assert 0 < d["transmission"] <= 1


def test_thresholds_single_source_with_paramcheck():
    """引擎阈值必须与 paramcheck 单一真源一致（改一处另一处跟着变）。"""
    assert uf.MI_LIMIT == paramcheck.HARD_LIMITS["focused ultrasound"][2] == 1.9
    assert uf.TEMP_LIMIT_C == paramcheck.TEMP_LIMIT[2] == 2.0


# ---------------- 逆设计 ----------------

def test_inverse_intensity_closed_form_and_roundtrip():
    """逆设计闭式解正向复核残差 ~0。"""
    for target in (0.2, 0.5, 0.9):
        inv = uf.inverse_intensity(target, TR, (uf.TISSUE_PRESETS["skull_thin"],))
        assert inv["residual_rel"] < 1e-9
        assert inv["required_surface_intensity_w_cm2"] > 0
        assert inv["recomputed_pressure_mpa"] == pytest.approx(target, rel=1e-9)


def test_inverse_flags_mi_infeasible_target():
    """目标声压越 MI 阈值 → 标不可行，不给'有效但危险'处方。"""
    inv = uf.inverse_intensity(3.0, TR, (uf.TISSUE_PRESETS["skull_thin"],))
    assert inv["feasible"] is False
    assert inv["mi_exceeded"] is True
    assert inv["infeasibility_note"]
    feasible = uf.inverse_intensity(0.5, TR, (uf.TISSUE_PRESETS["skull_thin"],))
    assert feasible["feasible"] is True


def test_inverse_scales_quadratically_with_pressure():
    """声压∝√I ⇒ 目标声压翻倍 → 所需表面声强 4 倍（可证伪）。"""
    a = uf.inverse_intensity(0.3, TR, (uf.TISSUE_PRESETS["skull_thin"],))
    b = uf.inverse_intensity(0.6, TR, (uf.TISSUE_PRESETS["skull_thin"],))
    assert b["required_surface_intensity_w_cm2"] == pytest.approx(
        4 * a["required_surface_intensity_w_cm2"], rel=1e-9
    )


# ---------------- 解剖敏感性 ----------------

def test_skull_sensitivity_direction_and_ratio():
    """同设备不同颅骨厚度 → 焦点声压差异（PNAS 论点超声版）。"""
    s = uf.skull_sensitivity(TR, (0.2, 0.4, 0.6, 0.8))
    vals = [c["peak_pressure_mpa"] for c in s["cases"]]
    assert all(a > b for a, b in zip(vals, vals[1:]))  # 越厚越低，单调
    assert s["peak_pressure_max_min_ratio"] > 1.5
    assert "设备参数" in s["interpretation"]


# ---------------- 独立验证基准 ----------------

def test_axis_rs_integral_matches_closed_form():
    """数值 RS 求积 vs 解析闭式：两种实现路径必须吻合（RMS<1%）。"""
    r = uv.axis_comparison()
    assert r["passed"] is True
    assert r["rms_diff"] < 0.01
    assert r["max_abs_diff"] < 0.02


def test_validation_report_all_gates_pass():
    """整套验证判据门全绿；surrogate 单列不混入门。"""
    rep = uv.build_report()
    assert rep["overall_passed"] is True
    for k in rep["gated_checks"]:
        assert rep["checks"][k].get("passed") is True, k
    assert "focus_gain_surrogate" not in rep["gated_checks"]
    assert rep["checks"]["focus_gain_surrogate"]["role"].startswith("measure")


def test_temperature_scaling_is_falsifiable():
    """温升粗估对声强/占空比/时长的线性必须成立；不满足即失败。"""
    r = uv.temperature_scaling()
    assert r["passed"] is True
    assert all(r["checks"].values())
