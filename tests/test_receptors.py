import pandas as pd
import pytest

from interface_atlas import receptors


def test_cached_matrix_shape_and_coverage():
    df = receptors.load_matrix()
    assert set(df["gene"].unique()) == set(receptors.RECEPTOR_PANEL)
    assert df["cell_type"].nunique() == 155
    assert len(df) == 13 * 155
    assert df.isna().sum().sum() == 0


def test_audit_passes_on_cached_data():
    report = receptors.audit(receptors.load_matrix())
    assert report["status"] == "PASS"
    assert report["genes"] == 13
    assert report["cell_types"] == 155
    # 广布接收器必须被显式预警，不能假装细胞特异
    assert "COX4I1" in report["low_selectivity_genes_at_1ncpm"]


def test_threshold_sensitivity_is_monotonic():
    mat = pd.DataFrame({"G": [0.0, 0.4, 1.2, 6.0, 15.0]},
                       index=["c0", "c1", "c2", "c3", "c4"])
    counts = receptors.threshold_sensitivity(mat)["G"]
    values = [counts[t] for t in receptors.THRESHOLDS]
    assert values == sorted(values, reverse=True)
    assert values == [3, 3, 2, 1, 0]


def test_build_matrix_requires_real_source():
    if not receptors.HPA_ZIP.exists():
        pytest.skip("HPA zip 未下载")
    mat = receptors.build_matrix()
    assert mat.shape == (155, 13)
    assert mat.isna().sum().sum() == 0
    assert mat["COX4I1"].max() > mat["COX4I1"].min()


def test_load_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        receptors.load_matrix(tmp_path / "nope.csv")
