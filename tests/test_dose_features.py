"""剂量场特征层的可证伪测试：多指标、可信度标注、数值场体积分交叉核对。"""
import numpy as np
import pytest

from interface_atlas import dose_features as df
from interface_atlas import electric_field as ef
from interface_atlas import light_field as lf
from interface_atlas import ultrasound_field as uf


def test_light_features_multimetric():
    r = lf.dose_field(lf.TissueOptics(0.3, 70.0), max_depth_cm=2.0, target_depth_cm=1.0)
    v = df.from_light_field(r)
    assert v.modality == "photobiomodulation"
    assert v.peak > 0 and v.decay_per_cm > 0
    assert 0.9 in v.coverage_depths_cm and 0.5 in v.coverage_depths_cm
    assert 1.0 / np.e == pytest.approx(next(k for k in v.coverage_depths_cm if 0.3 < k < 0.4))
    # 1/e 穿透深度必须来自插值而非 fallback（target_depth*2=2.0）
    assert v.penetration_depth_cm < 1.0


def test_light_focality_honestly_unavailable():
    r = lf.dose_field(lf.TissueOptics(0.3, 70.0), target_depth_cm=1.0)
    assert df.from_light_field(r).reliability["focality_proxy"] == "unavailable"


def test_ultrasound_focality_needs_transducer():
    tr = uf.TRANSDUCER_PRESETS["tfus_transcranial"]
    r = uf.dose_field(tr, (uf.TISSUE_PRESETS["skull_thin"],))
    without = df.from_ultrasound(r)
    with_ = df.from_ultrasound(r, transducer=tr)
    assert without.reliability["focality_proxy"] == "unavailable"
    assert math_isnan(without.focality_proxy)
    assert with_.reliability["focality_proxy"] == "proxy"
    assert 0 < with_.focality_proxy < 1
    assert with_.penetration_depth_cm == pytest.approx(tr.focal_length_cm)


def math_isnan(x):
    return isinstance(x, float) and np.isnan(x)


def test_electric_focality_ratio():
    r = ef.dose_field(1.0)
    v = df.from_electric(r)
    assert v.peak == pytest.approx(r.field_v_m)
    assert v.focality_proxy == pytest.approx(r.center_field_v_m / r.field_v_m, rel=1e-6)
    assert v.reliability["fwhm_cm"] == "unavailable"


def test_numeric_field_volume_integration():
    """单峰值体素：聚焦性=1/N，覆盖体积=voxel^3（解析可核对）。"""
    g = np.zeros((9, 9, 9))
    g[4, 4, 4] = 1.0
    v = df.from_numeric_field(modality="focused ultrasound", subject_id="s",
                              peak=1.0, peak_unit="norm", field_grid=g, voxel_size_cm=0.5)
    assert v.focality_proxy == pytest.approx(1 / (9 ** 3))
    assert v.coverage_depths_cm[0.5] == pytest.approx(0.125)
    assert v.reliability["focality_proxy"] == "analytic"


def test_cohort_no_single_score():
    vs = []
    for i in range(5):
        layers = tuple(ef.TissueLayer(l.outer_radius_cm,
                       l.conductivity_s_m * (0.6 + 0.4 * i), l.name) for l in ef.HEAD_4SPHERE)
        vs.append(df.from_electric(ef.dose_field(1.0, layers=layers)))
    s = df.cohort_summary(vs)
    assert s["n_subjects"] == 5
    assert s["single_metric_leaderboard_forbidden"] is True
    assert s["peak"]["max_min_ratio"] > 1.0  # 个体差异确实存在
    assert "sd" in s["peak"] and "q05" in s["peak"]  # 分布而非点估计


def test_as_dict_json_serializable():
    import json
    r = lf.dose_field(lf.TissueOptics(0.3, 70.0), target_depth_cm=1.0)
    json.dumps(df.from_light_field(r).as_dict(), allow_nan=False)
