from interface_atlas import evidence


def test_seed_schema_and_validation():
    rows = evidence.load()
    # v1.2：证据库扩至 103 条，CSV 与内嵌 SEED 必须逐条一致（可复现重建）
    assert len(rows) >= 100
    assert len(rows) == len(evidence.SEED)
    assert all(tuple(r.keys()) == evidence.FIELDS for r in rows)
    assert evidence.validate_row({"id": "X"})  # 必填缺失必须报错
    bad = dict(rows[0])
    bad["evidence_grade"] = "A+"
    assert any("证据等级" in e for e in evidence.validate_row(bad))


def test_query_filters():
    # 断言不依赖固定条数，但必须随过滤条件变化，防止恒真/硬编码
    cherna = evidence.query(gene="CHRNA7")
    pbm = evidence.query(modality="photobiomodulation")
    assert len(cherna) > 5 and len(pbm) > len(cherna)
    assert all(r["receiver_gene"].upper() == "CHRNA7" for r in cherna)
    assert all("photobiomodulation" in r["modality"].lower() for r in pbm)
    assert evidence.query(grade_min="A")
    assert all(r["evidence_grade"] == "A" for r in evidence.query(grade_min="A"))
    assert all(r["human_evidence"] == "yes" for r in evidence.query(human_only=True))
    assert any(r["id"] == "PBM-012" for r in evidence.query(keyword="retina"))


def test_audit_v11_source_closed_and_dose_banded():
    report = evidence.audit()
    assert report["status"] == "PASS"
    assert report["schema_version"] == "1.2.0"
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
