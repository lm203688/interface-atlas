"""interface-atlas 统一命令行入口."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import __version__, calibration, claimcheck, dosimetry, evidence, link, mc_light, paramcheck, receptors
from . import light_field, ultrasound_field as usfield, us_validation
from . import electric_field as esfield, es_validation
from . import anatomy, dose_features, l3bench, l3_validation
from . import harvest, harvest_validation, l1_clauses


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


def _dosefield(args: argparse.Namespace) -> int:
    """L2 光剂量场：正向计算 + 逆设计（承接 L1 的 reproducible_with_model）。"""
    if args.preset:
        if args.preset not in light_field.TISSUE_PRESETS:
            print(f"未知预设 {args.preset}；可用：{', '.join(light_field.TISSUE_PRESETS)}")
            return 2
        tissue = light_field.TISSUE_PRESETS[args.preset]
    else:
        if args.mu_a is None or args.mu_sp is None:
            print("请提供 --preset 或同时提供 --mu-a 与 --mu-sp")
            return 2
        tissue = light_field.TissueOptics(
            float(args.mu_a), float(args.mu_sp), float(args.g), float(args.n),
            args.preset or "custom",
        )
    depths = None
    if args.depths:
        depths = [float(x) for x in args.depths.split(",")]

    if args.inverse is not None:
        inv = light_field.inverse_irradiance(
            float(args.inverse), tissue=tissue,
            target_depth_cm=float(args.target_depth), geometry=args.geometry,
        )
        print("## L2 逆设计（闭式解 + 正向复核）")
        print(f"- 光学参数: μa={tissue.mu_a}/cm, μs'={tissue.mu_s_prime}/cm, n={tissue.n}")
        print(f"- 靶区深度 {args.target_depth} cm, 目标相对剂量 {args.inverse}")
        print(f"- 需表面剂量倍数: {inv['required_surface_multiplier']:.4f}")
        print(f"- 正向复核残差: {inv['residual_rel']:.2e}")
        print(f"- 可行性: {'可行' if inv['feasible'] else '不可行'}")
        if inv.get("infeasibility_note"):
            print(f"  ⚠ {inv['infeasibility_note']}")
        print(f"> {inv['caveat']}")
        if args.as_json:
            print(json.dumps(inv, ensure_ascii=False, indent=2))
    else:
        res = light_field.dose_field(
            tissue, depths_cm=depths, max_depth_cm=float(args.max_depth),
            n_points=int(args.points), target_depth_cm=float(args.target_depth),
            geometry=args.geometry,
        )
        print(light_field.render(res))
        if args.as_json:
            print(json.dumps(res.as_dict(), ensure_ascii=False, indent=2))
    return 0


def _mcbench(args: argparse.Namespace) -> int:
    """L2 验证：独立 MC 对权威基准 + P1 代理误差 + 解剖敏感性。"""
    which = args.which
    n = int(args.photons)
    if which in ("all", "prahl"):
        r = mc_light.validate_prahl_benchmark(n_photons=n)
        print("=== Prahl/Wang-Jacques MCML 基准 ===")
        print(json.dumps(r, ensure_ascii=False, indent=2))
        print(f"[{'PASS' if r['passed'] else 'FAIL'}] 透过率 Tt 误差 {r['relative_error']['Tt']*100:+.3f}%")
    if which in ("all", "patterson"):
        r2 = mc_light.validate_patterson_mu_eff(n_photons=n)
        print("\n=== P1 扩散代理 vs 独立 MC 的误差测量 ===")
        print(json.dumps(r2, ensure_ascii=False, indent=2))
    if which in ("all", "anatomy"):
        r3 = mc_light.anatomy_sensitivity()
        print("\n=== 解剖敏感性（PNAS 论点）===")
        print(json.dumps(r3, ensure_ascii=False, indent=2))
    return 0


def _calibrate(args: argparse.Namespace) -> int:
    """L2 中阶：split conformal 剂量区间（给定组织 → 点估计 + 覆盖区间）。"""
    if args.rebuild:
        print("重建校准报告（跑独立 MC，约需数分钟）…")
        rep = calibration.build_calibration(
            n_fit=int(args.fit), n_calib=int(args.calib), n_test=int(args.test),
            n_photons=int(args.photons),
        )
        p = calibration.save_calibration(rep)
        print(f"校准报告已写入 {p}")
        print(f"  偏差因子 b={rep.bias:.4f}，q={ {k: round(v, 4) for k, v in rep.q_by_alpha.items()} }")
        print(f"  test 覆盖率={ {k: round(v, 3) for k, v in rep.coverage_by_alpha.items()} }")
        print(f"  反证区间砍半: {'PASS' if rep.falsification['passed'] else 'FAIL'}")
        if args.as_json:
            print(json.dumps(rep.as_dict(), ensure_ascii=False, indent=2))
        return 0 if rep.falsification["passed"] else 1

    if args.preset:
        if args.preset not in light_field.TISSUE_PRESETS:
            print(f"未知预设 {args.preset}；可用：{', '.join(light_field.TISSUE_PRESETS)}")
            return 2
        tissue = light_field.TISSUE_PRESETS[args.preset]
    elif args.mu_a is not None and args.mu_sp is not None:
        tissue = light_field.TissueOptics(float(args.mu_a), float(args.mu_sp),
                                          float(args.g), float(args.n), "custom")
    else:
        print("请提供 --preset 或同时 --mu-a/--mu-sp，或用 --rebuild 重建报告")
        return 2

    try:
        r = calibration.calibrated_interval(
            tissue, target_depth_cm=float(args.target_depth), alpha=float(args.alpha),
            geometry=args.geometry,
        )
    except FileNotFoundError as exc:
        print(f"错误：{exc}")
        return 2
    if args.as_json:
        print(json.dumps(r, ensure_ascii=False, indent=2))
    else:
        print(calibration.render(r))
    return 0


def _usfield(args: argparse.Namespace) -> int:
    """L2 超声剂量场：正向计算 + 逆设计（承接 L1 的 tFUS reproducible_with_model）。"""
    if args.preset:
        if args.preset not in usfield.TRANSDUCER_PRESETS:
            print(f"未知换能器预设 {args.preset}；可用：{', '.join(usfield.TRANSDUCER_PRESETS)}")
            return 2
        tr = usfield.TRANSDUCER_PRESETS[args.preset]
    else:
        if args.intensity is None:
            print("请提供 --preset 或 --intensity（表面声强 W/cm²）")
            return 2
        tr = usfield.Transducer(
            float(args.radius), float(args.focal), float(args.frequency),
            float(args.intensity), float(args.efficiency),
        )

    # 颅骨层：按名称预设或自定义厚度
    if args.skull_mm is not None:
        skull = usfield.TissueLayer(float(args.skull_mm) / 10.0, 2800.0, 1900.0, 25.0,
                                    1.0, 1.0, f"skull custom {args.skull_mm}mm")
    elif args.skull:
        if args.skull not in usfield.TISSUE_PRESETS:
            print(f"未知组织预设 {args.skull}；可用：{', '.join(usfield.TISSUE_PRESETS)}")
            return 2
        skull = usfield.TISSUE_PRESETS[args.skull]
    else:
        skull = usfield.TISSUE_PRESETS["skull_thin"]

    if args.inverse is not None:
        inv = usfield.inverse_intensity(
            float(args.inverse), tr, (skull,),
        )
        print("## L2 超声逆设计（闭式解 + 正向复核）")
        print(f"- 换能器: f={tr.frequency_mhz}MHz, a={tr.radius_cm}cm, F={tr.focal_length_cm}cm, F#={tr.f_number:.2f}")
        print(f"- 颅骨: {skull.name}, 衰减 {skull.attenuation_at(tr.frequency_mhz)*skull.thickness_cm:.2f}dB")
        print(f"- 靶区目标声压 {args.inverse} MPa → 所需表面声强 {inv['required_surface_intensity_w_cm2']:.4f} W/cm²")
        print(f"- 正向复核残差: {inv['residual_rel']:.2e}")
        print(f"- 目标处 MI: {inv['mi_at_target']:.3f}（阈值 {usfield.MI_LIMIT}）")
        print(f"- 可行性: {'可行' if inv['feasible'] else '不可行'}")
        if inv.get("infeasibility_note"):
            print(f"  ⚠ {inv['infeasibility_note']}")
        print(f"> {inv['caveat']}")
        if args.as_json:
            print(json.dumps(inv, ensure_ascii=False, indent=2))
        return 0

    result = usfield.dose_field(
        tr, (skull,),
        duty_cycle=float(args.duty), pulse_duration_s=float(args.duration),
    )
    print(usfield.render(result))
    if args.as_json:
        print(json.dumps(result.as_dict(), ensure_ascii=False, indent=2))
    return 0 if not result.mi_exceeded else 1


def _usbench(args: argparse.Namespace) -> int:
    """L2 超声验证：数值 RS 积分 / 独立验证全套 / 颅骨敏感性。"""
    which = args.which
    if which in ("all", "rs"):
        r = us_validation.axis_comparison()
        print("=== 平面活塞轴上场：数值 RS 求积 vs 解析闭式 ===")
        print(f"[{'PASS' if r['passed'] else 'FAIL'}] RMS={r['rms_diff']}, max={r['max_abs_diff']}")
        if which == "rs":
            print(json.dumps(r, ensure_ascii=False, indent=2))
    if which in ("all", "skull"):
        s = usfield.skull_sensitivity(
            usfield.TRANSDUCER_PRESETS["tfus_transcranial"], (0.2, 0.4, 0.6, 0.8)
        )
        print("\n=== 颅骨厚度敏感性（PNAS 论点超声版）===")
        print(json.dumps(s, ensure_ascii=False, indent=2))
    if which == "all":
        rep = us_validation.build_report()
        print("\n=== 独立验证全套 ===")
        print(us_validation.render_report(rep))
        return 0 if rep["overall_passed"] else 1
    return 0


def _esfield(args: argparse.Namespace) -> int:
    """L2 电场剂量场：电流 + 电极 → 靶区电场（含逆设计，四球解析基准）。"""
    if args.preset == "uniform":
        layers = esfield.UNIFORM_SPHERE
    elif args.preset == "4sphere":
        layers = esfield.HEAD_4SPHERE
    else:
        print(f"未知头模型预设 {args.preset}；可用：4sphere / uniform")
        return 2

    if args.inverse is not None:
        inv = esfield.inverse_current(
            float(args.inverse), layers=layers,
            target_depth_cm=float(args.target_depth),
            electrode_area_cm2=float(args.electrode_area),
            n_harmonics=int(args.harmonics),
        )
        print("## L2 电场逆设计（闭式解 + 正向复核）")
        print(f"- 头模型 {args.preset}，靶区深度 {args.target_depth} cm，目标场强 {args.inverse} V/m")
        print(f"- 所需注入电流 {inv['required_current_ma']:.4f} mA（场/电流 k={inv['field_per_ma_v_m']:.4f} V/m/mA）")
        print(f"- 正向复核残差: {inv['residual_rel']:.2e}")
        print(f"- 实践窗口(≤{inv['practice_current_max_ma']}mA, advisory): "
              f"{'内' if inv['within_practice_window'] else '超出'}")
        print(f"  ⚠ {inv['note']}")
        print(f"> {inv['caveat']}")
        if args.as_json:
            print(json.dumps(inv, ensure_ascii=False, indent=2, default=str))
        return 0

    result = esfield.dose_field(
        float(args.current), layers=layers,
        target_depth_cm=float(args.target_depth),
        electrode_area_cm2=float(args.electrode_area),
        n_harmonics=int(args.harmonics),
    )
    print(esfield.render(result))
    if args.as_json:
        print(json.dumps(result.as_dict(), ensure_ascii=False, indent=2))
    return 0


def _esbench(args: argparse.Namespace) -> int:
    """L2 电场验证：单层闭式交叉核对 / 独立验证全套 / 颅骨敏感性。"""
    which = args.which
    if which in ("all", "closed"):
        r = es_validation.single_sphere_vs_closed_form()
        print("=== 单层球谐 vs 解析闭式 ===")
        print(f"[{'PASS' if r['passed'] else 'FAIL'}] 最坏相对误差 {r['worst_relative_error']}")
        if which == "closed":
            print(json.dumps(r, ensure_ascii=False, indent=2, default=str))
    if which in ("all", "skull"):
        s = esfield.skull_conductivity_sensitivity()
        print("\n=== 颅骨电导率敏感性（个体差异论点）===")
        print(json.dumps(s, ensure_ascii=False, indent=2, default=str))
    if which == "all":
        rep = es_validation.build_report()
        print("\n=== 独立验证全套 ===")
        print(es_validation.render_report(rep))
        return 0 if rep["overall_passed"] else 1
    return 0


def _headmodel(args: argparse.Namespace) -> int:
    """L2 个体化解剖：内置代理体模 / 群体抽样 / 许可门禁。"""
    if args.action == "list":
        print("内置代理体模（非真实影像数据）：")
        for name, m in anatomy.BUILTIN_MODELS.items():
            ok, _ = anatomy.check_license(m, commercial=True)
            print(f"  {name:24s} 颅骨={m.skull_thickness_mm():.2f}mm "
                  f"许可={m.license} 商用={'允许' if ok else '禁止'} "
                  f"{'[illustrative]' if m.illustrative else ''}")
        return 0
    if args.action == "show":
        m = anatomy.get_model(args.model)
        print(anatomy.render(m))
        if args.as_json:
            print(json.dumps(m.as_dict(), ensure_ascii=False, indent=2))
        return 0
    if args.action == "sample":
        models = anatomy.population_sample(n=args.n, seed=args.seed)
        print(f"群体解剖抽样 n={args.n} seed={args.seed}（illustrative 示例分布）：")
        skulls = [m.skull_thickness_mm() for m in models]
        print(f"  颅骨厚度 mm：min={min(skulls):.2f} max={max(skulls):.2f} "
              f"均值={sum(skulls)/len(skulls):.2f}")
        if args.as_json:
            print(json.dumps([m.as_dict() for m in models], ensure_ascii=False, indent=2))
        return 0
    if args.action == "license":
        m = anatomy.get_model(args.model)
        ok, reason = anatomy.check_license(m, commercial=args.commercial)
        print(f"[{'PASS' if ok else 'BLOCK'}] {m.subject_id}: {reason}")
        return 0 if ok else 1
    return 1


def _dosefeatures(args: argparse.Namespace) -> int:
    """L3 供给：从三模态剂量场抽取多指标特征。"""
    if args.modality == "light":
        r = light_field.dose_field(light_field.TissueOptics(args.mu_a, args.mu_sp),
                                   max_depth_cm=2.0, target_depth_cm=args.depth)
        v = dose_features.from_light_field(r, subject_id=args.subject)
    elif args.modality == "ultrasound":
        tr = usfield.TRANSDUCER_PRESETS[args.preset]
        layers = (usfield.TISSUE_PRESETS["skull_thin"],) if args.skull else ()
        r = usfield.dose_field(tr, layers)
        v = dose_features.from_ultrasound(r, subject_id=args.subject, transducer=tr)
    else:
        r = esfield.dose_field(args.current, target_depth_cm=args.depth)
        v = dose_features.from_electric(r, subject_id=args.subject)
    print(dose_features.render(v))
    if args.as_json:
        print(json.dumps(v.as_dict(), ensure_ascii=False, indent=2, default=str))
    return 0


def _l3bench(args: argparse.Namespace) -> int:
    """L3 剂量-效应基准：种子案例 / 独立验证 / 提交评测。"""
    if args.action == "cases":
        for c in l3bench.CASES:
            print(l3bench.render_case(c))
            print()
        return 0
    if args.action == "validate":
        rep = l3_validation.build_report()
        print(l3_validation.render_report(rep))
        if args.as_json:
            print(json.dumps(rep, ensure_ascii=False, indent=2, default=str))
        return 0 if rep["overall_passed"] else 1
    if args.action == "submit":
        if not args.file:
            print("需 --file 提供提交 JSON（case_id -> 各指标数值）")
            return 1
        subs = json.loads(Path(args.file).read_text(encoding="utf-8"))
        r = l3bench.run_benchmark(subs)
        print(f"评测：passed={r['summary']['passed']} failed={r['summary']['failed']} "
              f"insufficient={r['summary']['insufficient_evidence']}")
        for c in r["cases"]:
            s = c["score"]
            verdict = {True: "PASS", False: "FAIL", None: "N/A(证据不足)"}[s["passed"]]
            print(f"  [{verdict:14s}] {s['case_id']}: {s['note']}")
        print(f"  聚合总分：{r['aggregate_score']}（{r['aggregate_refused_reason']}）")
        if args.as_json:
            print(json.dumps(r, ensure_ascii=False, indent=2, default=str))
        return 0
    return 1


def _harvest(args: argparse.Namespace) -> int:
    """P1 证据产线：候选裁判 / 独立验证。永不自动入库。"""
    if args.validate:
        rep = harvest_validation.build_report()
        print(harvest_validation.render_report(rep))
        if args.as_json:
            print(json.dumps(rep, ensure_ascii=False, indent=2, default=str))
        return 0 if rep["overall_passed"] else 1
    if not args.input:
        print("需 --input 提供候选 JSON 文件，或用 --validate 跑产线独立验证")
        return 1
    candidates = json.loads(Path(args.input).read_text(encoding="utf-8"))
    if not isinstance(candidates, list):
        print("候选文件须为 JSON 数组")
        return 1
    rep = harvest.process_candidates(candidates, offline=args.offline)
    print(harvest.render(rep))
    if args.save:
        p = harvest.save_review_queue(rep, args.save)
        print(f"人审队列已写入：{p}")
    if args.as_json:
        print(json.dumps(rep, ensure_ascii=False, indent=2, default=str))
    return 0


def _l1clauses(args: argparse.Namespace) -> int:
    """L1 规范可引用条文注册表：列出 / 引用单条 / 核对 / 重建。"""
    if args.save:
        path = l1_clauses.save_registry()
        ok, _ = l1_clauses.check_registry(path)
        print(f"注册表已重建：{path}")
        print(f"sha256={l1_clauses.fingerprint(path)}")
        return 0 if ok else 1
    if args.check:
        ok, msg = l1_clauses.check_registry()
        print(("PASS " if ok else "FAIL ") + msg)
        if ok:
            print(f"sha256={l1_clauses.fingerprint()}")
        return 0 if ok else 1
    if args.cite:
        item = l1_clauses.find_clause(args.cite)
        if item is None:
            print(f"未找到条文 ID：{args.cite}（用 --list 查看全部）")
            return 2
        if args.as_json:
            print(json.dumps(item, ensure_ascii=False, indent=2))
        else:
            print(l1_clauses.format_clause(item))
        return 0
    data = json.loads(l1_clauses.REGISTRY_PATH.read_text(encoding="utf-8"))
    if args.list:
        if args.as_json:
            print(json.dumps(data, ensure_ascii=False, indent=2))
        else:
            for item in data["clauses"] + data["rules"]:
                print(l1_clauses.format_clause(item))
            print(f"共 {len(data['clauses'])} 条条文 + {len(data['rules'])} 条判定规则；"
                  f"sha256={l1_clauses.fingerprint()}")
        return 0
    print(f"L1 条文注册表：{data['counts']['clauses']} 条条文 + "
          f"{data['counts']['rules']} 条判定规则")
    print(f"规范版本 {data['spec_version']}；代码锚点 {data['citation']['version_pin']}")
    print(f"sha256={l1_clauses.fingerprint()}")
    print("用 --list / --cite <ID> / --check / --save 查看或维护；引用格式见 CITING_L1.md")
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

    df = sub.add_parser("dosefield", help="L2 光剂量场：设备/组织参数 → 靶区剂量（含逆设计）")
    df.add_argument("--preset", help="组织预设，如 skin_nir / fat_nir / muscle_nir")
    df.add_argument("--mu-a", dest="mu_a", help="吸收系数 /cm（与 --mu-sp 一起用）")
    df.add_argument("--mu-sp", dest="mu_sp", help="约化散射系数 /cm")
    df.add_argument("--g", default="0.90", help="各向异性因子，默认 0.90")
    df.add_argument("--n", default="1.40", help="相对折射率，默认 1.40")
    df.add_argument("--target-depth", dest="target_depth", default="0.5", help="靶区深度 cm，默认 0.5")
    df.add_argument("--max-depth", dest="max_depth", default="2.0", help="最大深度 cm，默认 2.0")
    df.add_argument("--points", default="41", help="深度采样点数，默认 41")
    df.add_argument("--depths", help="自定义深度列表，逗号分隔，如 0.01,0.1,0.5,1.0")
    df.add_argument("--geometry", choices=["plane", "point"], default="plane", help="照射几何，默认 plane")
    df.add_argument("--inverse", help="逆设计：目标靶区相对剂量，反推所需表面剂量倍数")
    df.add_argument("--as-json", dest="as_json", action="store_true", help="额外输出 JSON")
    df.set_defaults(func=_dosefield)

    mb = sub.add_parser("mcbench", help="L2 验证实验：独立 MC 基准 / P1 代理误差 / 解剖敏感性")
    mb.add_argument("--which", choices=["all", "prahl", "patterson", "anatomy"], default="all")
    mb.add_argument("--photons", default="200000", help="光子数，默认 200000")
    mb.set_defaults(func=_mcbench)

    cb = sub.add_parser("calibrate", help="L2 中阶：split conformal 剂量区间（有限样本覆盖保证）")
    cb.add_argument("--preset", help="组织预设，如 muscle_nir / skin_nir")
    cb.add_argument("--mu-a", dest="mu_a", help="吸收系数 /cm（与 --mu-sp 一起用）")
    cb.add_argument("--mu-sp", dest="mu_sp", help="约化散射系数 /cm")
    cb.add_argument("--g", default="0.90", help="各向异性因子，默认 0.90")
    cb.add_argument("--n", default="1.40", help="相对折射率，默认 1.40")
    cb.add_argument("--target-depth", dest="target_depth", default="0.5", help="靶区深度 cm")
    cb.add_argument("--alpha", default="0.10", help="1-alpha 覆盖目标，默认 0.10")
    cb.add_argument("--geometry", choices=["plane", "point"], default="plane")
    cb.add_argument("--as-json", dest="as_json", action="store_true", help="输出 JSON")
    cb.add_argument("--rebuild", action="store_true", help="重建校准报告（跑独立 MC，较慢）")
    cb.add_argument("--fit", default="8", help="rebuild: fit 段参数点数")
    cb.add_argument("--calib", default="40", help="rebuild: calib 段参数点数")
    cb.add_argument("--test", default="24", help="rebuild: test 段参数点数")
    cb.add_argument("--photons", default="80000", help="rebuild: 每点光子数")
    cb.set_defaults(func=_calibrate)

    us = sub.add_parser("usfield", help="L2 超声剂量场：换能器参数 → 靶区声压/MI（含逆设计）")
    us.add_argument("--preset", help="换能器预设：tfus_transcranial / tfus_superficial / tfus_deep")
    us.add_argument("--radius", default="3.0", help="自定义孔半径 cm（无 preset 时用）")
    us.add_argument("--focal", default="8.0", help="自定义焦距 cm")
    us.add_argument("--frequency", default="0.5", help="自定义频率 MHz")
    us.add_argument("--intensity", help="表面平均声强 W/cm²（换能器功率÷孔面积，无 preset 时必填）")
    us.add_argument("--efficiency", default="0.7", help="聚焦效率 η，默认 0.7")
    us.add_argument("--skull", help="颅骨组织预设，默认 skull_thin；可用 water/skull_thin/skull_thick")
    us.add_argument("--skull-mm", dest="skull_mm", help="自定义颅骨厚度 mm（覆盖 --skull）")
    us.add_argument("--duty", default="1.0", help="占空比 (0,1]，默认 1.0")
    us.add_argument("--duration", default="1.0", help="照射时长 s（温升粗估），默认 1.0")
    us.add_argument("--inverse", help="逆设计：目标靶区峰值声压 MPa，反推所需表面声强")
    us.add_argument("--as-json", dest="as_json", action="store_true", help="额外输出 JSON")
    us.set_defaults(func=_usfield)

    ub = sub.add_parser("usbench", help="L2 超声验证：RS 积分交叉核对 / 独立验证全套 / 颅骨敏感性")
    ub.add_argument("--which", choices=["all", "rs", "skull"], default="all")
    ub.set_defaults(func=_usbench)

    es = sub.add_parser("esfield", help="L2 电场剂量场：电流 + 电极 → 靶区电场（四球解析基准，含逆设计）")
    es.add_argument("--preset", default="4sphere", help="头模型：4sphere（四层）/ uniform（均匀球）")
    es.add_argument("--current", default="1.0", help="注入电流 mA，默认 1.0")
    es.add_argument("--target-depth", dest="target_depth", default="2.2",
                    help="靶区距头皮深度 cm，默认 2.2（落在脑层）")
    es.add_argument("--electrode-area", dest="electrode_area", default="35.0",
                    help="电极面积 cm²（估计电极下电流密度），默认 35")
    es.add_argument("--harmonics", default="32", help="球谐截断阶数，默认 32")
    es.add_argument("--inverse", help="逆设计：目标靶区场强 V/m，反推所需注入电流")
    es.add_argument("--as-json", dest="as_json", action="store_true", help="额外输出 JSON")
    es.set_defaults(func=_esfield)

    eb = sub.add_parser("esbench", help="L2 电场验证：单层闭式交叉核对 / 全套 / 颅骨敏感性")
    eb.add_argument("--which", choices=["all", "closed", "skull"], default="all")
    eb.set_defaults(func=_esbench)

    hm = sub.add_parser("headmodel", help="L2 个体化解剖：代理体模 / 群体抽样 / 许可门禁")
    hm.add_argument("action", choices=["list", "show", "sample", "license"])
    hm.add_argument("--model", default="spherical", help="体模名（show/license 用）")
    hm.add_argument("--n", type=int, default=100, help="群体抽样数量")
    hm.add_argument("--seed", type=int, default=20261119, help="群体抽样随机种子")
    hm.add_argument("--commercial", action="store_true", help="按商用用途校验许可")
    hm.add_argument("--as-json", dest="as_json", action="store_true")
    hm.set_defaults(func=_headmodel)

    dfp = sub.add_parser("dosefeatures", help="L3 供给：三模态剂量场多指标特征")
    dfp.add_argument("modality", choices=["light", "ultrasound", "electric"])
    dfp.add_argument("--subject", default="analytic")
    dfp.add_argument("--preset", default="tfus_transcranial")
    dfp.add_argument("--skull", action="store_true", help="超声：是否计入颅骨层")
    dfp.add_argument("--mu-a", dest="mu_a", type=float, default=0.3)
    dfp.add_argument("--mu-sp", dest="mu_sp", type=float, default=70.0)
    dfp.add_argument("--depth", type=float, default=2.2)
    dfp.add_argument("--current", type=float, default=1.0)
    dfp.add_argument("--as-json", dest="as_json", action="store_true")
    dfp.set_defaults(func=_dosefeatures)

    lb = sub.add_parser("l3bench", help="L3 剂量-效应基准：种子案例 / 验证 / 提交评测")
    lb.add_argument("action", choices=["cases", "validate", "submit"])
    lb.add_argument("--file", help="提交 JSON 路径（submit 用）")
    lb.add_argument("--as-json", dest="as_json", action="store_true")
    lb.set_defaults(func=_l3bench)

    hv = sub.add_parser("harvest", help="P1 证据产线：候选裁判（不自动入库）/ 独立验证")
    hv.add_argument("--input", help="候选 JSON 数组文件")
    hv.add_argument("--offline", action="store_true", help="不联网，仅做 schema/去重裁判")
    hv.add_argument("--save", help="把裁判报告写入人审队列 JSON 路径")
    hv.add_argument("--validate", action="store_true", help="跑产线独立验证（不需要 --input）")
    hv.add_argument("--as-json", dest="as_json", action="store_true")
    hv.set_defaults(func=_harvest)

    lc = sub.add_parser("l1clauses", help="L1 规范可引用条文：列出 / 引用单条 / 核对 / 重建")
    lc.add_argument("--list", action="store_true", help="打印全部条文与判定规则")
    lc.add_argument("--cite", metavar="ID", help="打印单条条文的引用卡，如 IA-L1-PBM-T2-applicator_distance")
    lc.add_argument("--check", action="store_true", help="核对注册表与权威字段源 DOSE_SPEC 未脱节")
    lc.add_argument("--save", action="store_true", help="从 DOSE_SPEC 重建注册表（权威源变更时）")
    lc.add_argument("--as-json", dest="as_json", action="store_true")
    lc.set_defaults(func=_l1clauses)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
