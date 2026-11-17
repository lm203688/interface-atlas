"""A1 ↔ A2 交叉核验：证据库声明的细胞，在 HPA 单细胞数据里是否真有该接收器.

原则：匹配不上就写 UNMATCHED，绝不把"文献说的细胞"偷换成"看起来差不多的细胞"。
别名来自 data/aliases.csv（可审计、可版本控制），每条带 basis 依据。
空 hpa_cell_type 表示该证据不是单一细胞类型（患者群体/组织/生物标志物），显式跳过。

v0.3 起输出携带证据库版本指纹，保证 A1/A2/link 三者可复现对齐。
"""
from __future__ import annotations

import csv
from pathlib import Path

from . import evidence, receptors

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
ALIASES_CSV = DATA_DIR / "aliases.csv"

FIELDS = ("id", "gene", "evidence_cell_type", "matched_hpa_cell_type",
          "n_cpm", "pct_rank", "status", "note")


def load_aliases(path: Path | str = ALIASES_CSV) -> dict[str, list[tuple[str, str]]]:
    """读取别名表：cell_type -> [(hpa_cell_type, basis), ...]."""
    p = Path(path)
    out: dict[str, list[tuple[str, str]]] = {}
    with open(p, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            key = row["evidence_cell_type"].strip()
            target = row["hpa_cell_type"].strip()
            out.setdefault(key, []).append((target, row.get("basis", "").strip()))
    return out


def _match_cell_type(raw: str, available: set[str],
                     aliases: dict[str, list[tuple[str, str]]]) -> tuple[str, str]:
    key = raw.strip()
    if key in aliases:
        targets = [(t, b) for t, b in aliases[key] if t]
        if not targets:
            basis = aliases[key][0][1] if aliases[key] else "非单一细胞类型"
            return "", f"UNMATCHED: {basis}"
        missing = [t for t, _ in targets if t not in available]
        if missing:
            return "", f"UNMATCHED: 别名 {missing} 不在 HPA 细胞表中"
        if len(targets) > 1:
            return targets[0][0], f"AMBIGUOUS_CANDIDATES: {[t for t, _ in targets]}；取首位展示，需人工选择"
        return targets[0][0], f"ALIAS: {targets[0][1]}"
    hit = sorted(c for c in available if key.lower() in c.lower())
    if len(hit) == 1:
        return hit[0], "EXACT_SUBSTRING"
    if len(hit) > 1:
        return "", f"AMBIGUOUS: {len(hit)} 个候选，需人工选择"
    return "", "UNMATCHED: HPA 中无对应细胞类型"


def crosslink(evidence_rows: list[dict] | None = None,
              receptor_rows=None,
              aliases: dict[str, list[tuple[str, str]]] | None = None) -> list[dict[str, str]]:
    erows = evidence_rows if evidence_rows is not None else evidence.load()
    rrows = receptor_rows if receptor_rows is not None else receptors.load_matrix()
    aliases = aliases if aliases is not None else load_aliases()
    available = set(rrows["cell_type"].unique())
    lookup = {(r.gene, r.cell_type): (r.n_cpm, r.pct_rank) for r in rrows.itertuples()}

    out: list[dict[str, str]] = []
    for e in erows:
        gene = e["receiver_gene"]
        matched, note = _match_cell_type(e["cell_type"], available, aliases)
        if matched:
            n_cpm, pct = lookup.get((gene, matched), (None, None))
            status = "MATCHED" if n_cpm is not None else "UNMATCHED: 该细胞无此基因记录"
            out.append({
                "id": e["id"], "gene": gene, "evidence_cell_type": e["cell_type"],
                "matched_hpa_cell_type": matched,
                "n_cpm": "" if n_cpm is None else str(n_cpm),
                "pct_rank": "" if pct is None else str(pct),
                "status": status, "note": note,
            })
        else:
            out.append({
                "id": e["id"], "gene": gene, "evidence_cell_type": e["cell_type"],
                "matched_hpa_cell_type": "", "n_cpm": "", "pct_rank": "",
                "status": "UNMATCHED", "note": note,
            })
    return out


def summarize(rows: list[dict[str, str]],
              evidence_rows: list[dict] | None = None) -> dict:
    from collections import Counter
    counts = dict(Counter(r["status"].split(":")[0] for r in rows))
    total = len(rows)
    matched = counts.get("MATCHED", 0)
    ambiguous = counts.get("AMBIGUOUS_CANDIDATES", 0)
    non_single = sum(1 for r in rows if r["status"] == "UNMATCHED"
                     and "非单一细胞类型" in r["note"])
    return {
        "total": total,
        "status_counts": counts,
        "matched_rate": round(matched / max(total, 1), 3),
        "matched_plus_candidates_rate": round((matched + ambiguous) / max(total, 1), 3),
        "non_single_cell_records": non_single,
        "evidence_fingerprint": evidence.fingerprint(evidence_rows),
        "caveat": "nCPM 是 RNA 水平，不等于膜蛋白丰度；免疫受体表达高度依赖激活状态。"
                  "UNMATCHED 不等于证据错误，只表示本图谱无法核验。"
                  "患者群体/组织/生物标志物不是细胞类型，不计入可核验率分母的业务解读。",
    }
