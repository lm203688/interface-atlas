"""L2 光剂量场与独立 MC 的可证伪测试。"""
import math

import pytest

from interface_atlas import light_field as lf
from interface_atlas import mc_light as mc


# ---------------- light_field: 正向 ----------------

def test_dose_field_monotonic_in_absorption():
    """吸收越强 → 深部剂量越低、漫反射越低（必须单调，不能恒真）。"""
    fs, rds = [], []
    for mu_a in (0.5, 1.0, 2.0, 4.0):
        t = lf.TissueOptics(mu_a, 100.0)
        r = lf.dose_field(t, target_depth_cm=0.5, depths_cm=(1e-6, 0.5))
        fs.append(r.target_fluence_rel)
        rds.append(r.rd_total)
    assert all(a > b for a, b in zip(fs, fs[1:])), fs
    assert all(a > b for a, b in zip(rds, rds[1:])), rds


def test_dose_field_surface_normalized_and_decaying():
    t = lf.TISSUE_PRESETS["skin_nir"]
    r = lf.dose_field(t, target_depth_cm=0.5, max_depth_cm=2.0, n_points=41)
    assert r.points[0].fluence_rel == pytest.approx(1.0, abs=1e-3)  # z→0+ 归一化
    assert all(p.fluence_rel >= q.fluence_rel for p, q in zip(r.points, r.points[1:]))
    assert r.points[-1].fluence_rel < r.points[0].fluence_rel  # 确实衰减


def test_mu_eff_analytic_matches_definition():
    t = lf.TissueOptics(2.0, 100.0)
    c = lf.optical_coefficients(t)
    assert c["mu_eff"] == pytest.approx(math.sqrt(3 * 2.0 * 102.0))
    assert c["a_prime"] == pytest.approx(100.0 / 102.0)


def test_diffusion_applicability_warns_on_high_absorption():
    strong = lf.TissueOptics(mu_a=30.0, mu_s_prime=100.0)  # mu_a/mu_s'=0.3
    r = lf.dose_field(strong, target_depth_cm=0.5, depths_cm=(1e-6, 0.5))
    assert any("扩散近似" in w for w in r.warnings)
    weak = lf.TissueOptics(mu_a=1.0, mu_s_prime=100.0)
    r2 = lf.dose_field(weak, target_depth_cm=0.5, depths_cm=(1e-6, 0.5))
    assert not any("扩散近似" in w for w in r2.warnings)


def test_reflectance_bounds():
    for t in lf.TISSUE_PRESETS.values():
        rd = lf.diffuse_reflectance_total(t)
        assert 0.0 <= rd <= 1.0


# ---------------- inverse ----------------

def test_inverse_design_closed_form_and_roundtrip():
    t = lf.TISSUE_PRESETS["muscle_nir"]
    for target in (0.005, 0.01, 0.05):
        inv = lf.inverse_irradiance(target, tissue=t, target_depth_cm=1.0)
        assert inv["required_surface_multiplier"] > 0
        assert inv["residual_rel"] < 1e-9
        assert inv["target_fluence_recomputed"] == pytest.approx(target, rel=1e-9)


def test_inverse_flags_physically_infeasible_targets():
    """深靶区 + 强吸收光学参数 → 表面倍数爆炸，必须标不可行而非给天文数字。"""
    t = lf.TissueOptics(4.0, 100.0)  # 强吸收皮肤
    inv = lf.inverse_irradiance(0.01, tissue=t, target_depth_cm=1.5)
    assert inv["feasible"] is False
    assert inv["infeasibility_note"]
    shallow = lf.inverse_irradiance(0.5, tissue=t, target_depth_cm=0.05)
    assert shallow["feasible"] is True


def test_inverse_scales_linearly_with_target():
    """目标剂量翻倍 → 所需表面倍数翻倍（线性，可证伪）。"""
    t = lf.TISSUE_PRESETS["skin_nir"]
    a = lf.inverse_irradiance(0.01, tissue=t, target_depth_cm=0.5)
    b = lf.inverse_irradiance(0.02, tissue=t, target_depth_cm=0.5)
    assert b["required_surface_multiplier"] == pytest.approx(
        2 * a["required_surface_multiplier"], rel=1e-9
    )


def test_gradient_is_nonzero_and_direction_correct():
    """数值梯度：μa↑深部剂量↓；μs'↑深部剂量也↓（散射增强使衰减更快）。

    后者是非平凡的物理结论——μeff=√(3μa(μa+μs')) 随 μs' 增大，穿透变浅。
    两个方向都必须为负且非零，否则说明可微接口失效。
    """
    def f(mu_a, mu_s_prime):
        return lf.dose_field(
            lf.TissueOptics(mu_a, mu_s_prime), target_depth_cm=0.5,
            depths_cm=(1e-6, 0.5),
        ).target_fluence_rel
    g = lf.gradient(f, {"mu_a": 2.0, "mu_s_prime": 100.0})
    assert g["mu_a"] < 0
    assert g["mu_s_prime"] < 0
    assert abs(g["mu_a"]) > 1e-9 and abs(g["mu_s_prime"]) > 1e-9


# ---------------- independent MC ----------------

def test_mc_pure_absorption_matches_beer_lambert():
    from interface_atlas.mc_light import Layer, run_mc
    r = run_mc((Layer(10.0, 0.0, 0.0, 0.1),), n_photons=30000,
               seed=1, n_bins=10, max_depth_cm=0.1)
    assert r.transmittance == pytest.approx(math.exp(-1.0), abs=0.01)


def test_mc_prahl_benchmark_passes():
    """独立 MC 必须复现权威透过率 Tt（剂量核心量），能量守恒。"""
    r = mc.validate_prahl_benchmark(n_photons=60000)
    assert r["passed"] is True
    assert abs(r["relative_error"]["Tt"]) <= 0.02
    assert r["energy_balance_error"] < 1e-6
    # Rd 偏差如实记录（不隐藏）
    assert "rd_known_limitation" in r


def test_mc_surrogate_error_is_reported_not_hidden():
    """P1 代理误差是测量输出，不是 pass/fail；必须带解释字段。"""
    r = mc.validate_patterson_mu_eff(n_photons=40000)
    assert r["role"].startswith("measure")
    assert "surrogate_relative_error" in r
    assert "interpretation" in r


def test_anatomy_sensitivity_direction():
    """解剖差异：脂肪更厚 → 靶区剂量更高（光通道），方向必须正确。"""
    r = mc.anatomy_sensitivity()
    vals = [c["target_absorption_rel"] for c in r["cases"]]
    assert vals[0] < vals[-1]  # 0mm 脂肪 < 5mm 脂肪
    assert r["max_min_ratio"] > 1
