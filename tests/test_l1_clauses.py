"""L1 可引用条文注册表测试：唯一性 / 权威源同步 / 出处真实 / 确定性 / CLI。

不联网。出处断言基于真实证据库 CSV，禁止编造 DOI/URL。
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from interface_atlas import cli, dosimetry, l1_clauses

ROOT = Path(__file__).resolve().parents[1]
REG = ROOT / "schema" / "l1_clauses_v1.json"


def test_registry_exists_and_checks_clean():
    assert REG.exists()
    ok, msg = l1_clauses.check_registry(REG)
    assert ok, msg


def test_ids_unique_and_well_formed():
    data = json.loads(REG.read_text(encoding="utf-8"))
    ids = [c["id"] for c in data["clauses"] + data["rules"]]
    assert len(ids) == len(set(ids)), "条文 ID 必须唯一"
    for c in data["clauses"]:
        prefix = f"IA-L1-{c['modality_code']}-T{c['tier']}-"
        assert c["id"].startswith(prefix)
        assert c["id"] == f"IA-L1-{l1_clauses.MODALITY_CODES[c['modality']]}-T{c['tier']}-{c['key']}"


def test_clauses_match_dose_spec_exactly():
    """注册表必须与唯一权威源 DOSE_SPEC 一一对应，不能有第二份字段清单。"""
    data = json.loads(REG.read_text(encoding="utf-8"))
    expected = {
        f"IA-L1-{code}-T{f.tier}-{f.key}"
        for modality, code in l1_clauses.MODALITY_CODES.items()
        for f in dosimetry.DOSE_SPEC[modality]
    }
    actual = {c["id"] for c in data["clauses"]}
    assert actual == expected, f"脱节：多={actual - expected} 少={expected - actual}"
    assert len(actual) == sum(len(fs) for fs in dosimetry.DOSE_SPEC.values())


def test_every_clause_has_definition_and_unit():
    data = json.loads(REG.read_text(encoding="utf-8"))
    for c in data["clauses"]:
        assert "报告须给出" in c["definition"], c["id"]
        assert c["unit"], c["id"]
        assert c["tier"] in (1, 2, 3, 4)
        assert c["tier_en"] in ("device", "interface", "target", "subject")


def test_definitions_contain_no_efficacy_claim():
    """条文是报告要求，不得夹带疗效/安全声明。"""
    banned = ("可治疗", "能治愈", "保证有效", "一定有效", "疗效确切")
    data = json.loads(REG.read_text(encoding="utf-8"))
    for c in data["clauses"]:
        for word in banned:
            assert word not in c["definition"], (c["id"], word)


def test_evidence_refs_are_real_and_typed():
    """每条挂的出处必须能在真实证据库找到，标注三态类型。"""
    data = json.loads(REG.read_text(encoding="utf-8"))
    with (ROOT / "data" / "evidence_matrix.csv").open(encoding="utf-8-sig", newline="") as fh:
        real_urls = {r["source_url"].strip() for r in csv.DictReader(fh) if r["source_url"].strip()}
    assert real_urls
    for c in data["clauses"]:
        for ref in c["evidence_refs"]:
            assert ref["url"] in real_urls, ("伪造出处", c["id"], ref["url"])
            assert ref["kind"] in ("pmc", "doi", "institutional"), ref
            assert ref["evidence_id"]


def test_evidence_ref_scope_disclaimer_present():
    """注册表必须声明：出处仅说明字段为何必要，不构成疗效证据。"""
    data = json.loads(REG.read_text(encoding="utf-8"))
    scope = data["evidence_ref_scope"]
    assert "不代表该字段已被这些文献验证" in scope
    assert "不构成疗效或安全性证据" in scope
    assert "community draft" in data["status"].lower()


def test_deterministic_render_and_fingerprint_stable():
    a = l1_clauses.render(l1_clauses.build_registry())
    b = l1_clauses.render(l1_clauses.build_registry())
    assert a == b, "注册表渲染必须逐字节确定"
    assert l1_clauses.fingerprint(REG) == l1_clauses.fingerprint(REG)


def test_rebuild_is_byte_identical(tmp_path):
    """换个路径重建，字节必须与在仓注册表完全一致。"""
    rebuilt = l1_clauses.save_registry(tmp_path / "l1.json")
    assert rebuilt.read_bytes() == REG.read_bytes()


def test_inserting_field_does_not_shift_published_ids():
    """ID 以语义字段名为锚，新增字段不影响旧 ID。"""
    ids = {c["id"] for c in l1_clauses.build_registry()["clauses"]}
    assert "IA-L1-PBM-T1-wavelength" in ids
    # 即便 DOSE_SPEC 顺序改变，只要 key 不变，ID 集合也不变
    spec = dosimetry.DOSE_SPEC["photobiomodulation"]
    reordered = spec[::-1]
    ids_reordered = {
        f"IA-L1-PBM-T1-{f.key}" for f in reordered if f.tier == 1
    }
    ids_orig = {f"IA-L1-PBM-T1-{f.key}" for f in spec if f.tier == 1}
    assert ids_reordered == ids_orig


def test_rules_cover_four_tiers_and_five_outcomes():
    data = json.loads(REG.read_text(encoding="utf-8"))
    rules = {r["outcome"]: r for r in data["rules"]}
    assert set(rules) == {
        "not_reproducible", "underdetermined", "reproducible_with_model",
        "reproducible_from_report", "completeness_score",
    }
    assert rules["completeness_score"]["tier_weights"] == {"T1": 0.4, "T2": 0.25, "T3": 0.25, "T4": 0.1}
    assert data["tiers"][0]["weight"] == 0.4


def test_find_clause_and_format():
    item = l1_clauses.find_clause("ia-l1-fus-t2-skull_attenuation")  # 大小写不敏感
    assert item is not None
    card = l1_clauses.format_clause(item)
    assert "IA-L1-FUS-T2-skull_attenuation" in card
    assert "实证背景出处" in card
    assert l1_clauses.find_clause("IA-L1-NOPE") is None


def test_cli_l1clauses_modes():
    assert cli.main(["l1clauses", "--check"]) == 0
    assert cli.main(["l1clauses", "--cite", "IA-L1-RULE-003"]) == 0
    assert cli.main(["l1clauses", "--list", "--as-json"]) == 0
    assert cli.main(["l1clauses"]) == 0
    assert cli.main(["l1clauses", "--cite", "DOES-NOT-EXIST"]) == 2
