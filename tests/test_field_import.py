"""P2 真解导入 + 全指标的可证伪测试。

不依赖真实外部数据（受许可约束、本仓不打包）；用合成场做解析交叉核对，
依赖（nibabel/h5py/scipy）缺失时用 pytest.importorskip 跳过对应用例，
核心契约与 NPZ 通路零新增依赖、必测。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from interface_atlas import cli, dose_features as df, field_import as fi


# --- NumericField 契约 ------------------------------------------------------

def test_unit_is_mandatory():
    g = np.zeros((3, 3, 3))
    g[1, 1, 1] = 1.0
    with pytest.raises(ValueError, match="显式声明"):
        fi.from_array(g, field_unit="", voxel_size_cm=0.1)
    with pytest.raises(ValueError):
        fi.from_array(g, field_unit="V/m", voxel_size_cm=0.0)


def test_anisotropic_voxel_volume():
    g = np.zeros((3, 4, 5))
    g[1, 2, 2] = 1.0
    f = fi.from_array(g, field_unit="V/m", voxel_size_cm=(0.1, 0.2, 0.5))
    assert f.voxel_size_cm == (0.1, 0.2, 0.5)
    assert f.voxel_volume_cm3 == pytest.approx(0.01)
    assert f.shape == (3, 4, 5)


def test_scale_recorded_and_applied():
    g = np.ones((2, 2, 2))
    f = fi.from_array(g, field_unit="MPa", voxel_size_cm=0.1,
                      scale=1e-6, source_ref="pa_field")
    assert f.scale_applied == 1e-6
    assert float(f.grid.max()) == pytest.approx(1e-6)
    assert f.provenance["source"] == "in-memory array"


def test_validate_grid_rejects_bad_shapes():
    with pytest.raises(ValueError, match="3D"):
        fi.from_array(np.zeros((3, 3)), field_unit="V/m", voxel_size_cm=0.1)
    # 单个 NaN 容忍（抽取器只取 finite）；全 NaN 才拒绝
    one_nan = np.zeros((2, 2, 2)); one_nan[0, 0, 0] = np.nan
    assert fi.from_array(one_nan, field_unit="V/m", voxel_size_cm=0.1).shape == (2, 2, 2)
    all_nan = np.full((2, 2, 2), np.nan)
    with pytest.raises(ValueError, match="非有限"):
        fi.from_array(all_nan, field_unit="V/m", voxel_size_cm=0.1)


# --- NPZ 开放契约（零依赖，必测）-------------------------------------------

def test_npz_roundtrip(tmp_path):
    g = np.zeros((5, 5, 5)); g[2, 2, 2] = 2.0
    f = fi.from_array(g, field_unit="V/m", voxel_size_cm=(0.2, 0.2, 0.4),
                      source_tool="kwave", source_ref="synthetic")
    p = fi.save_npz(f, tmp_path / "f.npz", provenance_json=json.dumps({"license": "synthetic"}))
    loaded = fi.from_npz(p, field_unit="V/m")
    assert loaded.unit == "V/m"
    assert loaded.voxel_size_cm == (0.2, 0.2, 0.4)
    assert float(loaded.grid.max()) == 2.0
    assert loaded.provenance["license"] == "synthetic"


def test_npz_unit_mismatch_warns(tmp_path):
    g = np.ones((2, 2, 2))
    f = fi.from_array(g, field_unit="Pa", voxel_size_cm=0.1)
    p = fi.save_npz(f, tmp_path / "u.npz")
    loaded = fi.from_npz(p, field_unit="MPa")
    assert any("不一致" in w for w in loaded.warnings)


def test_npz_missing_key_raises(tmp_path):
    p = tmp_path / "bad.npz"
    np.savez(p, wrong=np.ones((2, 2, 2)), voxel_size_cm=0.1)
    with pytest.raises(KeyError, match="field"):
        fi.from_npz(p, field_unit="V/m")


def test_npz_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        fi.from_npz(tmp_path / "nope.npz", field_unit="V/m")


# --- 全指标抽取的解析交叉核对 ----------------------------------------------

def _gaussian():
    z, y, x = np.mgrid[0:21, 0:21, 0:21]
    g = np.exp(-(((x - 10) / 3.0) ** 2 + ((y - 10) / 4.0) ** 2 + ((z - 10) / 5.0) ** 2))
    return fi.from_array(g, field_unit="V/m", voxel_size_cm=(0.1, 0.2, 0.5),
                         source_tool="array", source_ref="synthetic-gauss")


def test_full_metrics_analytic_against_gaussian():
    """高斯场 exp(-(d/a)^2)：FWHM=2·sqrt(ln2)·a·cadence，1/e=a·cadence。"""
    f = _gaussian()
    v = df.from_numeric_field_full(f, modality="transcranial electric",
                                   subject_id="s", axis=2)  # x 轴: a=3 vox, cadence=0.5
    assert v.peak == pytest.approx(1.0)
    expected_fwhm = 2 * math.sqrt(math.log(2)) * 3.0 * 0.5
    expected_pen = 3.0 * 0.5
    assert v.fwhm_cm == pytest.approx(expected_fwhm, rel=0.02)
    assert v.penetration_depth_cm == pytest.approx(expected_pen, rel=0.02)
    assert v.decay_per_cm > 0
    assert v.reliability["fwhm_cm"] == "analytic"
    assert v.reliability["focality_proxy"] == "analytic"


def test_auto_axis_is_proxy():
    f = _gaussian()
    v = df.from_numeric_field_full(f, modality="transcranial electric", subject_id="s")
    assert v.provenance["profile_axis_auto"] is True
    assert v.reliability["fwhm_cm"] == "proxy"
    assert any("自动选择" in w for w in v.warnings)


def test_volume_integration_single_voxel():
    g = np.zeros((9, 9, 9)); g[4, 4, 4] = 1.0
    f = fi.from_array(g, field_unit="norm", voxel_size_cm=0.5)
    v = df.from_numeric_field_full(f, modality="focused ultrasound", subject_id="s", axis=0)
    assert v.focality_proxy == pytest.approx(1 / (9 ** 3))
    assert v.coverage_volumes_cm3[0.5] == pytest.approx(0.125)
    assert v.provenance["total_volume_cm3"] == pytest.approx(9 ** 3 * 0.125)


def test_peak_at_boundary_has_no_crossing():
    """峰值在轴边界：该侧无交点，指标标 unavailable，不编造。"""
    g = np.zeros((5, 5, 5)); g[0, 2, 2] = 1.0
    g[1, 2, 2] = 0.2
    f = fi.from_array(g, field_unit="norm", voxel_size_cm=0.1)
    v = df.from_numeric_field_full(f, modality="focused ultrasound", subject_id="s", axis=0)
    assert v.provenance["peak_index"][0] == 0
    assert v.reliability["fwhm_cm"] == "unavailable"


def test_full_metrics_json_serializable():
    v = df.from_numeric_field_full(_gaussian(), modality="transcranial electric",
                                   subject_id="s", axis=1)
    d = v.as_dict()
    assert "coverage_volumes_cm3" in d
    json.dumps(d, allow_nan=False)


def test_no_cross_modality_normalization():
    """真解保留来源单位，不做跨模态归一化。"""
    v = df.from_numeric_field_full(_gaussian(), modality="focused ultrasound",
                                   subject_id="s", axis=2)
    assert v.peak_unit == "V/m"  # 即便模态标超声，也不擅自改单位


# --- 可选依赖：MAT v5（scipy 有）-------------------------------------------

def test_kwave_mat_v5_via_scipy(tmp_path):
    pytest.importorskip("scipy.io")
    from scipy.io import savemat
    g = np.zeros((6, 6, 6)); g[3, 3, 3] = 5.0  # Pa 量级
    p = tmp_path / "p.mat"
    savemat(str(p), {"p": g})
    f = fi.from_kwave_mat(p, field_unit="MPa", variable="p",
                          voxel_size_cm=(0.1, 0.1, 0.2), scale=1e-6)
    assert float(f.grid.max()) == pytest.approx(5e-6)
    assert f.voxel_volume_cm3 == pytest.approx(0.002)
    assert "scipy.io" in f.provenance["read_backend"]


def test_kwave_rejects_bad_scale(tmp_path):
    p = tmp_path / "p.mat"
    p.write_bytes(b"dummy")
    with pytest.raises(ValueError, match="scale"):
        fi.from_kwave_mat(p, field_unit="MPa", voxel_size_cm=0.1, scale=0)


# --- CLI --------------------------------------------------------------------

def test_cli_fieldimport_npz(tmp_path):
    g = np.zeros((5, 5, 5)); g[2, 2, 2] = 1.0
    f = fi.from_array(g, field_unit="V/m", voxel_size_cm=0.1)
    p = fi.save_npz(f, tmp_path / "f.npz")
    assert cli.main(["fieldimport", "npz", str(p), "--unit", "V/m", "--describe"]) == 0
    assert cli.main(["fieldimport", "npz", str(p), "--unit", "V/m",
                     "--axis", "2", "--as-json"]) == 0
    # unit 必填：argparse 层拦截
    with pytest.raises(SystemExit):
        cli.main(["fieldimport", "npz", str(p)])


def test_cli_kwave_requires_voxel():
    with pytest.raises(ValueError, match="--voxel"):
        cli._parse_voxel(None)
    assert cli._parse_voxel("0.1") == 0.1
    assert cli._parse_voxel("0.1,0.2,0.3") == (0.1, 0.2, 0.3)
