from interface_atlas import evidence, link


def test_crosslink_shape_and_status():
    rows = link.crosslink()
    # 交叉核验必须覆盖证据库全量，随数据规模自动伸缩
    assert len(rows) == len(evidence.load())
    assert len(rows) >= 100
    assert all(tuple(r.keys()) == link.FIELDS for r in rows)
    summary = link.summarize(rows)
    assert summary["total"] == len(rows)
    assert 0 < summary["matched_rate"] < 1  # 既不能全通过，也不能全失败


def test_cholinergic_receptor_gap_is_visible():
    rows = link.crosslink()
    hit = next(r for r in rows if r["id"] == "VNS-001")
    assert hit["status"] == "MATCHED"
    assert hit["matched_hpa_cell_type"] == "macrophages"
    assert float(hit["n_cpm"]) < 5.0  # 文献宣称与 RNA 水平之间的真实张力


def test_physical_barrier_not_falsely_verified():
    rows = link.crosslink()
    skull = [r for r in rows if r["evidence_cell_type"] == "human skull"]
    assert skull
    assert all(r["status"].startswith("UNMATCHED") for r in skull)


def test_alias_mapping_for_l_cells():
    rows = link.crosslink()
    hit = next(r for r in rows if r["id"] == "MET-001")
    assert hit["matched_hpa_cell_type"] == "enteroendocrine cells"
    assert hit["status"] == "MATCHED"
