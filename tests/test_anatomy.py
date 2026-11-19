"""L2 个体化解剖层的可证伪测试：几何守恒、许可门禁、群体抽样、导入契约。"""
import math

import pytest

from interface_atlas import anatomy as an


def test_geometry_conservation():
    for name, m in an.BUILTIN_MODELS.items():
        radii = [lay.outer_radius_cm for lay in m.layers]
        assert all(b > a for a, b in zip(radii, radii[1:])), name
        assert sum(m.layer_thicknesses().values()) == pytest.approx(m.outer_radius_cm), name


def test_skull_thickness_positive():
    for name, m in an.BUILTIN_MODELS.items():
        t = m.skull_thickness_mm()
        assert 2.0 < t < 20.0, (name, t)


def test_license_gate_fail_closed():
    """CC BY-NC / HCP DUA 商用必须被拒，MIT 通过；研究用途放行。"""
    assert an.check_license(an.get_model("spherical"), commercial=True)[0] is True
    assert an.check_license(an.get_model("ernie_extended_metadata"), commercial=True)[0] is False
    assert an.check_license(an.get_model("hcp_population"), commercial=True)[0] is False
    assert an.check_license(an.get_model("ernie_extended_metadata"), commercial=False)[0] is True


def test_require_license_raises():
    with pytest.raises(ValueError, match="许可门禁"):
        an.require_license(an.get_model("ernie_extended_metadata"), commercial=True)


def test_population_reproducible():
    a = an.population_sample(n=10, seed=123)
    b = an.population_sample(n=10, seed=123)
    assert [m.skull_thickness_mm() for m in a] == [m.skull_thickness_mm() for m in b]


def test_population_nondegenerate():
    s = an.population_sample(n=40, seed=9)
    skulls = [m.skull_thickness_mm() for m in s]
    assert max(skulls) - min(skulls) > 1.0
    assert all(m.illustrative for m in s)


def test_population_uses_custom_stats():
    """自定义统计（颅骨均值更大）应抬高抽样颅骨厚度，证明 load 钩子真的接进来了。"""
    st = {k: dict(v) for k, v in an.EXAMPLE_POPULATION_STATS.items()}
    st["skull_thickness_mm"]["mean"] = 12.0
    st["skull_thickness_mm"]["sd"] = 0.1
    base = an.population_sample(n=20, seed=3)
    custom = an.population_sample(n=20, seed=3, stats=st)
    assert sum(m.skull_thickness_mm() for m in custom) > sum(m.skull_thickness_mm() for m in base)
    assert all(not m.illustrative for m in custom)


def test_layer_geometry_validation():
    """半径非递增（内层 > 外层）必须 fail-fast。"""
    with pytest.raises(ValueError, match="严格递增"):
        an.HeadModel("x", "analytic", "MIT", (
            an.TissueLayerGeometry("brain", 9.0, 0.3, 0.5),   # 内层半径 9
            an.TissueLayerGeometry("skull", 8.0, 0.01, 20.0), # 外层半径 8 < 9 → 非递增
        ))


def test_anatomy_source_aligned_with_dosimetry():
    """tier-4 字段取值必须与 dosimetry 本体一致，不能另起炉灶。"""
    from interface_atlas import dosimetry as dos
    for m in an.BUILTIN_MODELS.values():
        assert m.anatomy_source in an.ANATOMY_SOURCE_CHOICES
    assert set(an.ANATOMY_SOURCE_CHOICES) == {"phantom", "subject_imaging", "population"}


def test_render_contains_license_and_disclaimer():
    txt = an.render(an.get_model("hcp_population"))
    assert "CC BY 4.0" in txt and "illustrative" in txt
