"""P1 证据产线测试：三态核验 / 去重 / 指纹 / 不自动入库 / CLI。

网络 I/O 全部用 monkeypatch 替换，测试不依赖外网、不发真实请求。
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from interface_atlas import cli, evidence, harvest, harvest_validation

ROOT = Path(__file__).resolve().parents[1]


def _rows() -> list[dict]:
    return evidence.load()


def _base_candidate(rows: list[dict]) -> dict:
    return dict(rows[0])


# --- 来源解析 ---------------------------------------------------------------

def test_doi_extraction_formats():
    assert harvest.extract_doi("https://doi.org/10.1016/J.BBI.2023.12.008") == "10.1016/j.bbi.2023.12.008"
    assert harvest.extract_doi("doi:10.1111/psyp.13933") == "10.1111/psyp.13933"
    assert harvest.extract_doi("10.1126/sciadv.abh2169") == "10.1126/sciadv.abh2169"
    assert harvest.extract_doi("https://example.org/no-doi") is None


def test_source_kind_all_types():
    assert harvest.source_kind("https://doi.org/10.1234/x") == "doi"
    assert harvest.source_kind("https://pmc.ncbi.nlm.nih.gov/articles/PMC123/") == "pmc"
    assert harvest.source_kind("https://pubmed.ncbi.nlm.nih.gov/41099683/") == "pubmed"
    assert harvest.source_kind("https://example.org/guide") == "unsupported"


# --- schema 与去重 ----------------------------------------------------------

def test_schema_accepts_real_row_and_rejects_bad():
    rows = _rows()
    assert harvest.check_schema(rows[0]) == []
    bad = {"id": "X", "effect_direction": "maybe", "evidence_grade": "Z", "source_url": "ftp://x"}
    errs = harvest.check_schema(bad)
    assert any("effect_direction" in e for e in errs)
    assert any("evidence_grade" in e for e in errs)
    assert any("source_url" in e for e in errs)


def test_duplicate_detection_three_ways():
    rows = _rows()
    assert any("ID 重复" in e for e in harvest.check_duplicate(rows[0], rows))
    by_url = {**rows[20], "id": "NEW-001"}
    assert any("来源 URL 重复" in e for e in harvest.check_duplicate(by_url, rows))
    by_triple = {**rows[20], "id": "NEW-002", "source_url": "https://example.org/v2"}
    assert any("已存在" in e for e in harvest.check_duplicate(by_triple, rows))


# --- 指纹 -------------------------------------------------------------------

def test_fingerprint_stable_and_tamper_evident():
    a = {"id": "X", "endpoint": "alpha", "modality": "m"}
    b = {"modality": "m", "endpoint": "alpha", "id": "X"}
    assert harvest.fingerprint(a) == harvest.fingerprint(b)
    assert harvest.fingerprint(a) != harvest.fingerprint({**a, "endpoint": "alpha!"})
    assert len(harvest.fingerprint(a)) == 64


# --- 三态核验（无真实网络） -------------------------------------------------

def test_tri_state_verification(monkeypatch):
    original = harvest._fetch_json
    cases = {
        "ok": ("verified", ("ok", {"result": {"1": {}}})),
        "not_found": ("refuted", ("not_found", None)),
        "rate_limited": ("unverified", ("rate_limited", None)),
        "error": ("unverified", ("error", "timeout")),
    }
    try:
        for label, (want, response) in cases.items():
            monkeypatch.setattr(harvest, "_fetch_json", lambda *a, r=response, **k: r)
            got = harvest.verify_source({"source_url": "https://doi.org/10.1234/probe"})["state"]
            assert got == want, label
    finally:
        monkeypatch.setattr(harvest, "_fetch_json", original)


def test_offline_never_touches_network(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("离线模式不得联网")
    monkeypatch.setattr(harvest, "_fetch_json", boom)
    v = harvest.verify_source({"source_url": "https://doi.org/10.1234/probe"}, offline=True)
    assert v["state"] == "unverified" and v["reason"] == "offline"


# --- 分流与不自动入库 --------------------------------------------------------

def test_process_candidates_routing(monkeypatch):
    rows = _rows()
    novel = {**rows[0], "id": "NEW-100", "endpoint": "novel endpoint X",
             "source_url": "https://doi.org/10.1234/novel"}
    dup = dict(rows[0])
    candidates = [novel, dup]

    rep = harvest.process_candidates(candidates, existing=rows, offline=True)
    assert rep["summary"]["needs_manual_verification"] == 1
    assert rep["summary"]["rejected"] == 1
    assert rep["summary"]["auto_committed"] == 0

    # 在线模式注入 404：应进入 rejected 而非 manual
    monkeypatch.setattr(harvest, "_fetch_json", lambda *a, **k: ("not_found", None))
    rep2 = harvest.process_candidates([novel], existing=rows)
    assert rep2["summary"]["rejected"] == 1
    assert rep2["rejected"][0]["reason"] == "source_refuted"


def test_no_auto_commit_preserves_evidence_file(tmp_path):
    csv_path = ROOT / "data" / "evidence_matrix.csv"
    before = csv_path.read_bytes()
    rows = _rows()
    rep = harvest.process_candidates([{**rows[0], "id": "NEW-999",
                                       "source_url": "https://example.org/x"}],
                                     existing=rows, offline=True)
    out = harvest.save_review_queue(rep, tmp_path / "queue.json")
    assert out.exists()
    assert json.loads(out.read_text(encoding="utf-8"))["summary"]["auto_committed"] == 0
    assert csv_path.read_bytes() == before
    assert len(_rows()) == len(rows)


# --- 独立验证与 CLI ---------------------------------------------------------

def test_independent_validation_report():
    rep = harvest_validation.build_report()
    assert rep["overall_passed"] is True
    assert all(v["passed"] for v in rep["checks"].values())


def test_cli_validate_and_harvest(tmp_path):
    assert cli.main(["harvest", "--validate"]) == 0
    payload = [{**_rows()[0], "id": "NEW-CLI", "endpoint": "cli novel",
                "source_url": "https://example.org/cli"}]
    p = tmp_path / "candidates.json"
    p.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    queue = tmp_path / "queue.json"
    assert cli.main(["harvest", "--input", str(p), "--offline", "--save", str(queue)]) == 0
    assert queue.exists()
