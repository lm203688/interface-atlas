"""interface-atlas 统一命令行入口."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import __version__, claimcheck, dosimetry, evidence, link, paramcheck, receptors


def _print_rows(rows: list[dict], columns: list[str]) -> None:
    if not rows:
        print("（无匹配记录）")
        return
    widths = {c: max(len(c), max(len(str(r.get(c, ""))) for r in rows)) for c in columns}
    print(" | ".join(c.ljust(widths[c]) for c in columns))
    print("-+-".join("-" * widths[c] for c in columns))
    for r in rows:
        print(" | ".join(str(r.get(c, "")).ljust(widths[c]) for c in columns))


def _receptors(args: argparse.Namespace) -> int:
    if args.rebuild:
        mat = receptors.build_matrix(args.zip)
        path = receptors.save_matrix(mat)
        print(f"已生成 {path}（{len(mat)} 细胞类型 × {len(mat.columns)} 接收器）")
    df = receptors.load_matrix()
    if args.audit:
        mat = receptors.build_matrix(args.zip) if args.rebuild else None
        report = receptors.audit(df, mat)
        for k, v in report.items():
            if k != "threshold_sensitivity":
                print(f"{k}: {v}")
        if args.audit_out:
            out = receptors.write_audit(report, args.audit_out)
            print(f"审计报告已写入 {out}")
        return 0 if report["status"] == "PASS" else 1
    if args.heatmap:
        out = receptors.write_heatmap(df, args.heatmap, args.top_n)
        print(f"热图已写入 {out}")
        return 0
    rows = receptors.top_cells(df, gene=args.gene, modality=args.modality, n=args.top).to_dict("records")
    _print_rows(rows, ["gene", "receptor", "cell_type", "n_cpm", "pct_rank", "specificity_tau"])
    return 0


def _evidence(args: argparse.Namespace) -> int:
    if args.init:
        path = evidence.save_seed(args.csv)
        print(f"已生成种子证据库 {path}（{len(evidence.SEED)} 条）")
    rows = evidence.load(args.csv)
    if args.audit:
        report = evidence.audit(rows)
        for k, v in report.items():
            print(f"{k}: {v}")
        return 0 if report["status"] == "PASS" else 1
    rows = evidence.query(rows, modality=args.modality, gene=args.gene,
                          cell_type=args.cell_type, grade_min=args.grade,
                          human_only=args.human, keyword=args.keyword)
    _print_rows(rows, ["id", "modality", "cell_type", "endpoint", "effect_direction",
                       "evidence_grade", "source_url"])
    return 0


def _link(args: argparse.Namespace) -> int:
    rows = link.crosslink()
    if args.unmatched_only:
        rows = [r for r in rows if not r["status"].startswith("MATCHED")]
    _print_rows(rows, ["id", "gene", "evidence_cell_type", "matched_hpa_cell_type",
                       "n_cpm", "status", "note"])
    print(f"\n汇总：{link.summarize(link.crosslink())}")
    return 0


def _claim(args: argparse.Namespace) -> int:
    result = claimcheck.check(args.text)
    if args.markdown:
        print(claimcheck.to_markdown(result))
    else:
        print(f"verdict: {result['verdict']}")
        print(f"reason:  {result['reason']}")
        print(f"matched: {', '.join(result['matched_ids']) or '无'}")
        print(f"counter: {', '.join(result['counter_ids']) or '无'}")
        print(f"falsification: {result['falsification']}")
        print(f"fingerprint: {result['evidence_fingerprint']}")
    return 0


def _paramcheck(args: argparse.Namespace) -> int:
    if args.batch:
        rows = evidence.load()
        if args.modality:
            rows = evidence.query(rows, modality=args.modality)
        reports = [paramcheck.check_evidence_row(r) for r in rows]
        from collections import Counter
        dist = Counter(r.verdict for r in reports)
        print(f"批量预检 {len(reports)} 条：{dict(dist)}")
        for rep in reports:
            if rep.verdict in ("OUT_OF_WINDOW", "NO_EVIDENCE"):
                print(f"  [{rep.verdict}] {rep.source_id}: {rep.reason}")
        return 0
    if args.json:
        data = json.loads(Path(args.json).read_text(encoding="utf-8"))
        report = paramcheck.check(data, source_id=str(data.get("id", "")))
    else:
        payload = {"modality": args.modality}
        for kv in args.set or []:
            if "=" not in kv:
                print(f"参数格式错误（应为 key=value）：{kv}")
                return 2
            k, v = kv.split("=", 1)
            payload[k.strip()] = v.strip()
        for k in ("cell_type", "endpoint", "source_url"):
            val = getattr(args, k, None)
            if val:
                payload[k] = val
        report = paramcheck.check(payload)
    if args.as_json:
        print(paramcheck.to_json(report))
    else:
        print(paramcheck.render(report))
    return 0 if report.verdict != "OUT_OF_WINDOW" else 1


def _doseaudit(args: argparse.Namespace) -> int:
    if args.batch:
        rows = evidence.load(Path(args.csv) if args.csv else evidence.EVIDENCE_CSV)
        if args.modality:
            rows = evidence.query(rows, modality=args.modality)
        reports = [dosimetry.audit(r, source_id=r.get("id", "")) for r in rows]
        from collections import Counter
        dist = Counter(r.reproducibility for r in reports)
        scores = [r.completeness_score for r in reports]
        print(f"剂量完整性批量审计 {len(reports)} 条：{dict(dist)}")
        print(f"平均完整度 {sum(scores)/len(scores):.3f}，最高 {max(scores):.3f}")
        print("\n可复现性最佳的前 5 条：")
        for r in sorted(reports, key=lambda x: -x.completeness_score)[:5]:
            print(f"  {r.completeness_score:.3f}  {r.reproducibility:24} {r.source_id}")
        return 0
    if args.json:
        data = json.loads(Path(args.json).read_text(encoding="utf-8"))
        rep = dosimetry.audit(data, source_id=str(data.get("id", "")))
    else:
        if not args.modality:
            print("错误：需指定 --modality 或 --json 或 --batch")
            return 2
        payload: dict = {"modality": args.modality}
        for kv in args.set or []:
            if "=" not in kv:
                print(f"参数格式错误（应为 key=value）：{kv}")
                return 2
            k, v = kv.split("=", 1)
            payload[k.strip()] = v.strip()
        if args.signal:
            payload["signal"] = args.signal
        rep = dosimetry.audit(payload)
    print(dosimetry.to_json(rep) if args.as_json else dosimetry.render(rep))
    return 0


def _dosegap(args: argparse.Namespace) -> int:
    csv_path = Path(args.csv) if args.csv else evidence.EVIDENCE_CSV
    reports = dosimetry.audit_library(csv_path)
    gap = dosimetry.gap_map(reports)
    out = Path(args.out) if args.out else Path("data/dose_gap_map.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(gap, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"字段缺失图谱已写入 {out}")
    print(f"  总记录 {gap['total_rows']}，平均完整度 {gap['overall']['mean_completeness']}")
    print(f"  可复现性分布 {gap['overall']['verdict_distribution']}")
    print("\n缺失率最高的 5 个字段：")
    for x in gap["field_missing_rate"][:5]:
        print(f"  {x['missing_rate']:5.0%}  {x['key']}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="interface-atlas", description="非药物接口层：接收器图谱 + 证据剂量图谱")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("receptors", help="A1 接收器能力图谱")
    r.add_argument("--rebuild", action="store_true", help="从 HPA zip 重建表达矩阵")
    r.add_argument("--zip", default=str(receptors.HPA_ZIP), help="HPA tsv.zip 路径")
    r.add_argument("--top", type=int, default=12, help="显示前 N 行")
    r.add_argument("--gene", help="按基因过滤")
    r.add_argument("--modality", help="按模态关键词过滤")
    r.add_argument("--audit", action="store_true", help="运行完整性审计")
    r.add_argument("--audit-out", help="审计 Markdown 输出路径")
    r.add_argument("--heatmap", help="热图 PNG 输出路径")
    r.add_argument("--top-n", type=int, default=40, help="热图展示的细胞类型数")
    r.set_defaults(func=_receptors)

    e = sub.add_parser("evidence", help="A2 证据剂量图谱")
    e.add_argument("--init", action="store_true", help="从内置种子重建证据 CSV")
    e.add_argument("--csv", default=str(evidence.EVIDENCE_CSV), help="证据 CSV 路径")
    e.add_argument("--audit", action="store_true", help="运行 schema 审计")
    e.add_argument("--modality", help="按模态过滤")
    e.add_argument("--gene", help="按接收器基因过滤")
    e.add_argument("--cell-type", help="按细胞类型过滤")
    e.add_argument("--grade", help="最低证据等级 A/B/C/D")
    e.add_argument("--human", action="store_true", help="仅人体证据")
    e.add_argument("--keyword", help="全文关键词")
    e.set_defaults(func=_evidence)

    l = sub.add_parser("link", help="交叉核验：证据库细胞 ↔ HPA 实际表达")
    l.add_argument("--unmatched-only", action="store_true", help="只显示无法核验的记录")
    l.set_defaults(func=_link)

    c = sub.add_parser("claim", help="宣称-证据对齐判定（可证伪）")
    c.add_argument("text", help="一句宣称，如『红光面板能抗衰老』")
    c.add_argument("--markdown", action="store_true", help="Markdown 格式输出")
    c.set_defaults(func=_claim)

    pc = sub.add_parser("paramcheck", help="参数报告预检（本地自查，非监管结论）")
    pc.add_argument("--modality", help="模态，如 photobiomodulation / focused ultrasound")
    pc.add_argument("--set", action="append", metavar="KEY=VALUE", help="参数键值对，可重复")
    pc.add_argument("--cell-type", help="靶细胞")
    pc.add_argument("--endpoint", help="终点指标")
    pc.add_argument("--source-url", help="一手来源 URL")
    pc.add_argument("--json", help="从 JSON 文件读入完整报告")
    pc.add_argument("--as-json", action="store_true", help="输出 JSON")
    pc.add_argument("--batch", action="store_true", help="批量预检证据库全部记录")
    pc.set_defaults(func=_paramcheck)

    da = sub.add_parser("doseaudit", help="剂量报告完整性审计（规范 v1，社区草案）")
    da.add_argument("--modality", help="模态，如 photobiomodulation")
    da.add_argument("--set", action="append", metavar="KEY=VALUE", help="字段键值对，可重复")
    da.add_argument("--signal", help="自由文本描述（可正则抽取，不替代必报字段）")
    da.add_argument("--json", help="从 JSON 文件读入完整报告")
    da.add_argument("--as-json", action="store_true", help="输出 JSON")
    da.add_argument("--batch", action="store_true", help="批量审计证据库")
    da.add_argument("--csv", help="证据库 CSV 路径（批量时用）")
    da.set_defaults(func=_doseaudit)

    dg = sub.add_parser("dosegap", help="生成字段缺失图谱（规范的实证基础）")
    dg.add_argument("--csv", help="证据库 CSV 路径")
    dg.add_argument("--out", help="输出 JSON 路径，默认 data/dose_gap_map.json")
    dg.set_defaults(func=_dosegap)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
