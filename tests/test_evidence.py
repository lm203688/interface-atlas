from interface_atlas import evidence


def test_seed_schema_and_validation():
    rows = evidence.load()
    assert len(rows) == 28
    assert all(tuple(r.keys()) == evidence.FIELDS for r in rows)
    assert evidence.validate_row({"id": "X"})  # 必填缺失必须报错
    bad = dict(rows[0])
    bad["evidence_grade"] = "A+"
    assert any("证据等级" in e for e in evidence.validate_row(bad))


def test_query_filters():
    assert len(evidence.query(gene="CHRNA7")) == 5
    assert len(evidence.query(modality="photobiomodulation")) == 12
    assert evidence.query(grade_min="A")
    assert all(r["evidence_grade"] == "A" for r in evidence.query(grade_min="A"))
    assert all(r["human_evidence"] == "yes" for r in evidence.query(human_only=True))
    assert evidence.query(keyword="retina")[0]["id"] == "PBM-012"


def test_audit_v11_source_closed_and_dose_banded():
    report = evidence.audit()
    assert report["status"] == "PASS"
    assert report["schema_version"] == "1.1.0"
    # 三条 A 级来源已补一手 DOI/URL，缺源归零
    assert report["missing_source_url"] == []
    assert report["grade_a_missing_source"] == []
    assert report["unresolved_notes"] == []
    # 剂量分层：数值/相对/物理参数/缺失四类齐全
    assert set(report["dose_bands"]) == {"numeric", "relative", "physical_parameter", "missing"}
    assert report["dose_bands"]["numeric"] >= 4
    # 缺剂量不阻断，但必须有警告；且不得伪造数值
    assert report["dose_warnings"]
    assert len(report["fingerprint"]) == 12


def test_fingerprint_is_stable_and_content_sensitive():
    rows = evidence.load()
    assert evidence.fingerprint(rows) == evidence.fingerprint(rows)
    changed = [dict(r) for r in rows]
    changed[0]["effect_size"] = changed[0]["effect_size"] + " modified"
    assert evidence.fingerprint(changed) != evidence.fingerprint(rows)


def test_no_duplicate_ids():
    ids = [r["id"] for r in evidence.load()]
    assert len(ids) == len(set(ids))
