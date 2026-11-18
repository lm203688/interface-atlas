"""L2 中阶校准区间测试：覆盖保证 + 反证 + 可证伪性。"""
import math

import pytest

from interface_atlas import calibration as cal
from interface_atlas.light_field import TissueOptics, TISSUE_PRESETS


@pytest.fixture(scope="module")
def report():
    return cal.load_calibration()


# --- 报告结构可复现 ---
def test_report_structure(report):
    for k in ("version", "bias_b", "q_by_alpha", "coverage_by_alpha",
              "falsification", "r_fit", "r_calib", "r_test", "domain", "seeds"):
        assert k in report
    assert report["bias_b"] > 1.0  # P1 系统性高估深部剂量
    assert len(report["r_fit"]) >= 4
    assert len(report["r_calib"]) >= 4
    assert len(report["r_test"]) >= 4


def test_mc_seed_reproducible():
    """同一 seed 的 μeff 比完全一致（无浮点非确定）。"""
    r1 = cal.mu_eff_ratio(0.1, 10.0, 0.9, n_photons=20000, seed=123)
    r2 = cal.mu_eff_ratio(0.1, 10.0, 0.9, n_photons=20000, seed=123)
    assert r1 is not None and r2 is not None
    assert r1 == r2


# --- conformal 覆盖率保证（核心可证伪断言） ---
def test_coverage_meets_target(report):
    """边缘覆盖保证的统计检验：n_test 个独立点的观测覆盖率不得显著低于 1-α。

    用二项分布的正态近似下界：若真实覆盖率=1-α，观测覆盖率低于
    (1-α) - 2·sqrt((1-α)·α/n) 的概率仅约 2.3%。这比硬阈值更诚实——小样本下
    观测值必然波动，conformal 保证的是边缘概率而非每个有限样本精确达标。
    """
    n_test = report["sample_sizes"]["test"]
    for a in cal.ALPHAS:
        cov = report["coverage_by_alpha"][str(a)]
        p = 1.0 - a
        lower = p - 2.0 * math.sqrt(p * a / n_test)
        assert cov >= lower, (
            f"alpha={a} 覆盖率 {cov} 显著低于目标 {p}（统计下界 {lower:.3f}）"
        )


def test_interval_widens_as_alpha_drops(report):
    """置信越高（alpha 越小）区间越宽——conformal 分位数单调性。"""
    q20 = report["q_by_alpha"]["0.2"]
    q10 = report["q_by_alpha"]["0.1"]
    q05 = report["q_by_alpha"]["0.05"]
    assert math.isfinite(q20) and math.isfinite(q10) and math.isfinite(q05)
    assert q05 >= q10 >= q20


# --- 反证：区间砍半覆盖率必须下降（防恒真工具） ---
def test_falsification_interval_halving(report):
    f = report["falsification"]
    assert f["passed"], f"反证失败：{f}"
    assert f["coverage_half_interval"] < f["coverage_full_interval"]


def test_independent_test_not_in_calibration(report):
    """test 与 calib 的 r 样本不完全重合（三段分离，非同集重放）。"""
    assert set(report["r_test"]) != set(report["r_calib"])


# --- 推理区间 ---
def test_calibrated_interval_ordering(report):
    t = TISSUE_PRESETS["muscle_nir"]
    r = cal.calibrated_interval(t, target_depth_cm=0.5, alpha=0.10, report=report)
    lo, hi = r["target_fluence_interval"]
    # 去偏点估计必须落在覆盖区间内
    assert 0 <= lo <= r["target_fluence_point"] <= hi <= 1.0
    mlo, mhi = r["mu_eff_interval_per_cm"]
    # 区间以去偏值 b·μeff_P1 为中心；原始 P1 低于区间（因 b>1 低估 μeff）
    mu_debiased = r["bias_b"] * r["mu_eff_p1_per_cm"]
    assert mlo <= mu_debiased <= mhi
    assert r["mu_eff_p1_per_cm"] < mlo
    # b>1 ⇒ 原始 P1 系统性高估深部剂量；去偏点 <= 原始 P1
    assert r["target_fluence_point"] <= r["target_fluence_p1_raw"]
    assert hi <= r["target_fluence_p1_raw"] + 1e-9


def test_deeper_target_wider_relative_interval(report):
    """靶区越深，μeff 误差放大越强 → 剂量区间相对宽度越大。"""
    t = TISSUE_PRESETS["muscle_nir"]
    shallow = cal.calibrated_interval(t, target_depth_cm=0.2, report=report)
    deep = cal.calibrated_interval(t, target_depth_cm=1.5, report=report)

    def rel_width(r):
        lo, hi = r["target_fluence_interval"]
        p = r["target_fluence_point"]
        return (hi - lo) / p if p > 1e-15 else float("inf")

    assert rel_width(deep) > rel_width(shallow)


def test_alpha_monotonic_interval_width(report):
    t = TISSUE_PRESETS["skin_nir"]
    kw = dict(target_depth_cm=0.4, report=report)
    w20 = cal.calibrated_interval(t, alpha=0.20, **kw)
    w05 = cal.calibrated_interval(t, alpha=0.05, **kw)
    span = lambda r: r["target_fluence_interval"][1] - r["target_fluence_interval"][0]
    assert span(w05) >= span(w20)


# --- 适用域守卫 ---
def test_in_domain_guard():
    # 域内：muscle_nir 量级 mu_a=1.5/mu_s'=75/g=0.9（弱吸收且在范围内）
    inside = cal.in_domain(TissueOptics(1.5, 75.0, 0.9, 1.4))
    assert inside["ok"] is True, inside["checks"]
    strong = cal.in_domain(TissueOptics(50.0, 15.0, 0.9, 1.4))
    assert strong["ok"] is False
    assert strong["checks"]["weak_absorption"] is False


def test_out_of_domain_flagged(report):
    t = TissueOptics(50.0, 15.0, 0.9, 1.4)  # 强吸收，扩散近似失效
    r = cal.calibrated_interval(t, target_depth_cm=0.5, report=report)
    assert r["in_domain"]["ok"] is False
    assert "不适用覆盖保证" in r["disclaimer"]


def test_render(report):
    t = TISSUE_PRESETS["muscle_nir"]
    r = cal.calibrated_interval(t, report=report)
    text = cal.render(r)
    assert "split conformal" in text
    assert "覆盖" in text
