from interface_atlas import claimcheck, evidence


def test_verdicts_are_falsifiable_not_hardcoded():
    """核心回归：四个判定都必须依被测宣称而变，不能恒返回同一值。"""
    claims = {
        "植入式迷走刺激能治疗类风湿关节炎": "SUPPORTED",
        "近红外光改善膝关节骨关节炎的功能": "SUPPORTED",
        "红光面板能抗衰老延寿": "REFUTED",
        "迷走神经刺激可以降低炎症因子": "REFUTED",
        "耳夹刺激能改善心率变异": "REFUTED",
        "特定频率电磁波能修复细胞": "NO_EVIDENCE",
    }
    results = {c: claimcheck.check(c)["verdict"] for c in claims}
    for claim, expected in claims.items():
        assert results[claim] == expected, f"{claim}: {results[claim]} != {expected}"
    assert len(set(results.values())) >= 3  # 至少三种不同判定


def test_supported_carries_traceable_evidence_and_falsification():
    r = claimcheck.check("近红外光改善膝关节骨关节炎的功能")
    assert r["verdict"] == "SUPPORTED"
    assert "PBM-005" in r["matched_ids"]
    assert r["falsification"]  # 必须写明推翻条件
    assert r["evidence_fingerprint"] == evidence.fingerprint()
    assert r["boundaries"]  # 必须带边界声明


def test_refuted_surfaces_counter_evidence_not_vague_uncertainty():
    r = claimcheck.check("迷走神经刺激可以降低炎症因子")
    assert r["verdict"] == "REFUTED"
    assert "VNS-004" in r["counter_ids"]
    assert not r["matched_ids"]


def test_anti_aging_cannot_be_upgraded_from_surrogate_endpoints():
    r = claimcheck.check("红光照射促进胶原蛋白合成从而抗衰老")
    assert r["verdict"] == "EXCEEDED"
    assert "抗衰" in r["reason"] or "延寿" in r["reason"]


def test_unparseable_claim_returns_no_evidence():
    r = claimcheck.check("量子共振")
    assert r["verdict"] == "NO_EVIDENCE"
    assert r["matched_ids"] == []


def test_markdown_output_contains_all_trace_fields():
    md = claimcheck.to_markdown(claimcheck.check("红光面板能抗衰老延寿"))
    for field in ("判定", "依据", "命中证据", "反证证据", "推翻条件", "边界"):
        assert field in md
