"""证据产线：LLM 抽取候选 → 确定性机器裁判 → 人审队列 → 指纹签发.

设计原则（继承项目质量门禁）：
  1. **不编造**：本模块不生成证据内容，只核验外部给的候选条目；
  2. **不自动入库**：核验通过也只是进入人审队列，绝不直接写 evidence_matrix.csv；
  3. **区分"证伪"与"无法验证"**：404 是 refuted，429 限流/离线/不支持的来源
     是 unverified，后者不判死刑——实测 CrossRef 会 429，把它当伪造会打瘫产线；
  4. **三入口无关**：纯内核模块，零第三方依赖（标准库 urllib/json/hashlib）。

核验三态：
  verified   来源可用且条目存在
  refuted    来源明确返回"不存在"（可推翻的硬失败）
  unverified 离线 / 限流 / 网络错误 / 不支持的来源类型（需人工核验）
"""
from __future__ import annotations

import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from . import evidence

USER_AGENT = "interface-atlas-evidence-harvest/1.2 (research use; DOI/source verification)"

REQUIRED_FIELDS = ("id", "modality", "signal", "endpoint", "effect_direction",
                   "evidence_grade", "source_url")
ALLOWED_DIRECTION = {"increase", "decrease", "physical_constraint", "mixed", "none", "biphasic"}
ALLOWED_GRADE = {"A", "B", "C", "D"}
ALLOWED_HUMAN = {"yes", "no", "mixed", ""}

CROSSREF = "https://api.crossref.org/works/"
NCBI_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


# ---------------------------------------------------------------------------
# 来源解析
# ---------------------------------------------------------------------------

def extract_doi(value: str) -> str | None:
    """从 DOI 链接 / doi: 前缀 / 裸 DOI 中提取规范 DOI（小写）。"""
    if not value:
        return None
    s = value.strip()
    m = re.search(r"(?:doi\.org/|doi:\s*)(10\.\d{4,9}/\S+)", s, re.I)
    if m:
        return m.group(1).rstrip("/.,;").lower()
    if re.match(r"^10\.\d{4,9}/\S+$", s, re.I):
        return s.rstrip("/.,;").lower()
    return None


def extract_pmcid(value: str) -> str | None:
    m = re.search(r"(?:PMC|articles/PMC)(\d+)", value, re.I)
    return m.group(1) if m else None


def extract_pmid(value: str) -> str | None:
    m = re.search(r"(?:pubmed\.ncbi\.nlm\.nih\.gov/|[\?&]id=)(\d+)", value, re.I)
    return m.group(1) if m else None


def source_kind(url: str) -> str:
    if extract_doi(url):
        return "doi"
    if extract_pmcid(url):
        return "pmc"
    if extract_pmid(url):
        return "pubmed"
    return "unsupported"


# ---------------------------------------------------------------------------
# 网络核验：可 monkeypatch 的单一 I/O 点
# ---------------------------------------------------------------------------

def _fetch_json(url: str, timeout: float = 10.0) -> tuple[str, Any]:
    """返回 (state, payload)，state ∈ ok | not_found | rate_limited | error。

    所有网络异常在此收敛，调用方只处理三态，不写散落的 try/except。
    """
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return "ok", json.load(r)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return "not_found", None
        if e.code == 429:
            return "rate_limited", None
        return "error", f"HTTP {e.code}"
    except Exception as e:  # 超时 / DNS / TLS
        return "error", f"{type(e).__name__}: {e}"


def verify_source(candidate: dict, *, offline: bool = False,
                  timeout: float = 10.0) -> dict:
    """三态核验来源。offline=True 时全部 unverified，不发任何网络请求。"""
    url = str(candidate.get("source_url", ""))
    kind = source_kind(url)
    if offline:
        unsupported = kind == "unsupported"
        return {"state": "unverified", "kind": kind,
                "reason": "unsupported_source" if unsupported else "offline",
                "detail": "来源类型无公开核验接口，需人工确认" if unsupported
                          else "离线模式：未联网核验，需人工确认"}
    if kind == "unsupported":
        return {"state": "unverified", "kind": kind, "reason": "unsupported_source",
                "detail": "来源类型无公开核验接口（机构页/指南页等），需人工确认"}

    if kind == "doi":
        doi = extract_doi(url)
        state, payload = _fetch_json(CROSSREF + urllib.parse.quote(doi, safe="/"), timeout)
        if state == "ok":
            return {"state": "verified", "kind": kind, "doi": doi,
                    "detail": "CrossRef 200"}
        if state == "not_found":
            return {"state": "refuted", "kind": kind, "doi": doi,
                    "detail": "CrossRef 404：DOI 不存在"}
        return {"state": "unverified", "kind": kind, "doi": doi,
                "reason": state, "detail": f"CrossRef 无法核验：{payload}"}

    # PMC / PubMed 共用 NCBI eutils
    db, uid = ("pmc", extract_pmcid(url)) if kind == "pmc" else ("pubmed", extract_pmid(url))
    url_q = f"{NCBI_SUMMARY}?db={db}&id={uid}&retmode=json"
    state, payload = _fetch_json(url_q, timeout)
    if state == "ok":
        result = (payload or {}).get("result", {})
        if result.get(str(uid), {}).get("error"):
            return {"state": "refuted", "kind": kind, "uid": uid,
                    "detail": f"NCBI 明确返回不存在：{result[str(uid)]['error']}"}
        return {"state": "verified", "kind": kind, "uid": uid, "detail": "NCBI eutils 命中"}
    if state == "not_found":
        return {"state": "refuted", "kind": kind, "uid": uid, "detail": "NCBI 404"}
    return {"state": "unverified", "kind": kind, "uid": uid,
            "reason": state, "detail": f"NCBI 无法核验：{payload}"}


# ---------------------------------------------------------------------------
# 确定性裁判
# ---------------------------------------------------------------------------

def check_schema(candidate: dict) -> list[str]:
    """字段与枚举校验（不判断内容的科学真伪）。"""
    errs = []
    if not isinstance(candidate, dict):
        return ["候选条目不是对象"]
    for key in REQUIRED_FIELDS:
        if not str(candidate.get(key, "")).strip():
            errs.append(f"缺必填字段 {key}")
    if candidate.get("effect_direction") not in ALLOWED_DIRECTION:
        errs.append(f"effect_direction 非法：{candidate.get('effect_direction')!r}")
    if candidate.get("evidence_grade") not in ALLOWED_GRADE:
        errs.append(f"evidence_grade 非法：{candidate.get('evidence_grade')!r}")
    if str(candidate.get("human_evidence", "")) not in ALLOWED_HUMAN:
        errs.append(f"human_evidence 非法：{candidate.get('human_evidence')!r}")
    url = str(candidate.get("source_url", ""))
    if url and not url.startswith(("https://", "http://")):
        errs.append("source_url 不是 http(s) 链接")
    return errs


def _norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(text).lower()).strip()


def check_duplicate(candidate: dict, existing: list[dict]) -> list[str]:
    """三重去重：ID / 来源 URL / 归一化(modality,endpoint,cell_type)。"""
    errs = []
    cid = str(candidate.get("id", "")).strip().lower()
    url = str(candidate.get("source_url", "")).strip()
    key = (_norm(candidate.get("modality", "")), _norm(candidate.get("endpoint", "")),
           _norm(candidate.get("cell_type", "")))
    for row in existing:
        if cid and str(row.get("id", "")).strip().lower() == cid:
            errs.append(f"ID 重复：{cid}")
        if url and str(row.get("source_url", "")).strip() == url:
            errs.append(f"来源 URL 重复：{url[:60]}")
        other = (_norm(row.get("modality", "")), _norm(row.get("endpoint", "")),
                 _norm(row.get("cell_type", "")))
        if all(k and k == o for k, o in zip(key, other)):
            errs.append(f"同模态+终点+细胞类型条目已存在：{row.get('id')}")
    return errs


def fingerprint(candidate: dict) -> str:
    """候选内容指纹（规范化 JSON 的 SHA-256），用于人审留痕与回归。"""
    canonical = json.dumps(candidate, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def process_candidates(candidates: list[dict], *, existing: list[dict] | None = None,
                       offline: bool = False) -> dict:
    """机器裁判主入口：分流 ready_for_review / needs_manual_verification / rejected。

    任何 schema 失败、重复或 refuted → rejected；unverified → 人工队列；
    verified 且无硬失败 → ready_for_review（仍需人签字才入库）。
    """
    if existing is None:
        existing = evidence.load()
    ready, manual, rejected = [], [], []
    for i, raw in enumerate(candidates):
        cand = dict(raw or {})
        fp = fingerprint(cand)
        errs = check_schema(cand)
        dups = check_duplicate(cand, existing) if not errs else []
        verification = {"state": "skipped", "reason": "schema_or_duplicate_failed"}
        if not errs and not dups:
            verification = verify_source(cand, offline=offline)
        item = {"index": i, "id": cand.get("id", f"candidate_{i}"),
                "fingerprint": fp, "candidate": cand,
                "schema_errors": errs, "duplicate_errors": dups,
                "verification": verification}
        if errs or dups:
            item["reason"] = "schema_or_duplicate_failed"
            rejected.append(item)
        elif verification["state"] == "refuted":
            item["reason"] = "source_refuted"
            rejected.append(item)
        elif verification["state"] == "unverified":
            item["reason"] = "source_unverified"
            manual.append(item)
        else:
            item["reason"] = "verified_pending_human_signoff"
            ready.append(item)

    return {
        "summary": {
            "input": len(candidates), "ready_for_review": len(ready),
            "needs_manual_verification": len(manual), "rejected": len(rejected),
            "offline": bool(offline), "auto_committed": 0,
        },
        "boundary": (
            "候选由 LLM 抽取，机器只做 schema/重复/来源核验；verified 不等于科学结论正确。"
            "本产线永不自动入库，全部条目须人工签字后由维护者合并。"
        ),
        "ready_for_review": ready,
        "needs_manual_verification": manual,
        "rejected": rejected,
    }


def save_review_queue(report: dict, path: str | Path = "data/harvest_review_queue.json") -> Path:
    """把报告写入人审队列（仅队列文件，不动证据库）。"""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return p


def render(report: dict) -> str:
    s = report["summary"]
    lines = [
        "# 证据产线裁判报告", "",
        f"- 输入候选：{s['input']}",
        f"- 待人工签字：{s['ready_for_review']}",
        f"- 需人工核验来源：{s['needs_manual_verification']}",
        f"- 拒绝：{s['rejected']}",
        f"- 自动入库：**{s['auto_committed']}**（本产线永不自动入库）", "",
        report["boundary"], "",
    ]
    for section in ("ready_for_review", "needs_manual_verification", "rejected"):
        if report[section]:
            lines.append(f"## {section}")
            for it in report[section]:
                lines.append(f"- {it['id']} | {it['reason']} | fp={it['fingerprint'][:12]}")
            lines.append("")
    return "\n".join(lines)
