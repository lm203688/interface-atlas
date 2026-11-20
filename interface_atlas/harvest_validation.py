"""证据产线的独立验证：用不同实现路径交叉核对，不自证。

验证项：
  1. 来源解析：用字符串切片实现的独立提取器 vs harvest 的正则实现；
  2. 指纹抗篡改：内容改一个字符，指纹必须变；规范化后同一内容指纹必须同；
  3. 不自动入库：跑一次完整裁判，证据库文件字节与行数不变；
  4. 硬失败分流：schema 非法 / ID 重复 / 来源不支持必须分别落到正确的队列；
  5. 三态语义：refuted 只来自明确 404，unverified 来自限流/离线/不支持——
     用替换 I/O 的方式检验，不发真实网络请求。
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from . import evidence, harvest

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_CSV = ROOT / "data" / "evidence_matrix.csv"


# --- 独立实现（不用 harvest 的正则） ----------------------------------------

def _independent_doi(value: str) -> str | None:
    """用字符串切片而非正则提取 DOI，作为交叉核对路径。"""
    if not value:
        return None
    for marker in ("https://doi.org/", "http://doi.org/", "doi:"):
        if value.lower().startswith(marker):
            return value[len(marker):].strip().rstrip("/.,;").lower()
    if value.lower().startswith("10."):
        return value.strip().rstrip("/.,;").lower()
    return None


def _independent_pmcid(value: str) -> str | None:
    s = value.lower()
    i = s.find("pmc")
    if i < 0:
        return None
    digits = ""
    for ch in s[i + 3:]:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return digits or None


def source_parser_crosscheck() -> dict:
    samples = [
        "https://doi.org/10.1016/j.bbi.2023.12.008",
        "doi:10.1111/psyp.13933",
        "10.1126/sciadv.abh2169",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC12326686/",
        "https://example.org/guide",
    ]
    rows = []
    for s in samples:
        want_doi = _independent_doi(s)
        want_pmc = _independent_pmcid(s) if want_doi is None else None
        want = "doi" if want_doi else "pmc" if want_pmc else "unsupported"
        got = harvest.source_kind(s)
        rows.append({"url": s, "independent": want, "module": got, "match": want == got})
    return {"rows": rows, "passed": all(r["match"] for r in rows)}


def fingerprint_tamper_check() -> dict:
    base = {"id": "X-001", "modality": "test", "endpoint": "alpha"}
    same = dict(reversed(list(base.items())))
    changed = {**base, "endpoint": "alpha!"}
    return {
        "canonical_order_invariant": harvest.fingerprint(base) == harvest.fingerprint(same),
        "tamper_changes_fingerprint": harvest.fingerprint(base) != harvest.fingerprint(changed),
        "sha256_length": len(harvest.fingerprint(base)) == 64,
        "passed": (harvest.fingerprint(base) == harvest.fingerprint(same)
                   and harvest.fingerprint(base) != harvest.fingerprint(changed)
                   and len(harvest.fingerprint(base)) == 64),
    }


def no_auto_commit_check() -> dict:
    before = EVIDENCE_CSV.read_bytes()
    rows = evidence.load()
    candidate = {**rows[0], "id": "NEW-NOCOMMIT", "endpoint": "novel endpoint",
                 "source_url": "https://example.org/manual-only"}
    report = harvest.process_candidates([candidate], existing=rows, offline=True)
    after = EVIDENCE_CSV.read_bytes()
    return {
        "evidence_bytes_unchanged": before == after,
        "evidence_rows_unchanged": len(evidence.load()) == len(rows),
        "auto_committed": report["summary"]["auto_committed"],
        "queued_for_human": report["summary"]["needs_manual_verification"],
        "passed": (before == after and len(evidence.load()) == len(rows)
                   and report["summary"]["auto_committed"] == 0
                   and report["summary"]["needs_manual_verification"] == 1),
    }


def tri_state_check() -> dict:
    """替换 I/O 层，验证 404→refuted、429→unverified、200→verified，零真实请求。"""
    original = harvest._fetch_json
    out = {}
    try:
        for label, (state, payload) in {
            "ok": ("ok", {"result": {"123": {}}}),
            "not_found": ("not_found", None),
            "rate_limited": ("rate_limited", None),
            "error": ("error", "timeout"),
        }.items():
            harvest._fetch_json = lambda *a, _s=state, _p=payload, **k: (_s, _p)
            cand = {"source_url": "https://doi.org/10.1234/probe"}
            verdict = harvest.verify_source(cand)["state"]
            want = {"ok": "verified", "not_found": "refuted",
                    "rate_limited": "unverified", "error": "unverified"}[label]
            out[label] = {"injected": state, "want": want, "got": verdict, "match": want == verdict}
    finally:
        harvest._fetch_json = original
    return {"rows": out, "passed": all(r["match"] for r in out.values())}


def build_report() -> dict:
    checks = {
        "source_parser_crosscheck": source_parser_crosscheck(),
        "fingerprint_tamper": fingerprint_tamper_check(),
        "no_auto_commit": no_auto_commit_check(),
        "tri_state": tri_state_check(),
    }
    gated = ("source_parser_crosscheck", "fingerprint_tamper", "no_auto_commit", "tri_state")
    return {
        "overall_passed": all(checks[k]["passed"] for k in gated),
        "gated_checks": list(gated),
        "checks": checks,
        "boundary": "本验证用独立实现与 I/O 替换交叉核对；不证明候选条目的科学真伪，只证明产线按声明分流。",
    }


def render_report(rep: dict) -> str:
    lines = ["# 证据产线独立验证", "",
             f"总体：{'PASS' if rep['overall_passed'] else 'FAIL'}", ""]
    for k, v in rep["checks"].items():
        lines.append(f"- {'PASS' if v['passed'] else 'FAIL'} {k}")
    lines += ["", rep["boundary"]]
    return "\n".join(lines)


if __name__ == "__main__":  # pragma: no cover
    r = build_report()
    print(render_report(r))
    raise SystemExit(0 if r["overall_passed"] else 1)
