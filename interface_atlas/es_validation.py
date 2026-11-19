"""L2 电场剂量场的独立验证基准（球谐求解器 vs 解析闭式 + 可证伪性质）.

为什么需要它：
    electric_field.py 的球谐矩阵求解器是数值实现。用自己的输出「验证」自己
    等于循环论证。需要一个不同实现路径的 ground truth：单层均匀球点电极对
    有严格解析闭式 E(r)=I/(2πσR²)·[2(1+u²)/(1-u²)²+1/(1-u²)]——由级数求和
    独立推出，与 2N×2N 矩阵求解是完全不同的路径。二者在单层退化下必须吻合。

验证清单（全部固定参数可复现）：
    1. 单层退化：球谐求解器（单层均匀球）vs 解析闭式，在扩散区逐点比对；
       判据 r≤7cm（远离表面点电极奇异区）相对误差 <1e-3。
    2. 近场截断如实报告：r→R 处（电极正下方）闭式发散、球谐有限阶截断，
       二者偏差随阶数增加而减小——这是截断误差，不是 bug，如实记录不掩盖。
    3. 电流线性 E∝I（可证伪）：I 翻倍 E 必翻倍。
    4. 电导率整体缩放反比性：所有层 σ 同乘 k → E 变为 1/k。
       推导：固定总电流 I ⇒ 表面 Neumann 边界 -σ∂φ/∂r=J(由 I 定，与 σ 无关)，
       σ→kσ 使 φ→φ/k、E→E/k；等价于「同电流下整体电阻降 k 倍⇒电压/场降 k 倍」。
       这是非平凡检验——若边界条件误写成 Dirichlet 就会失败。（注：本条初版
       期望写错成「E 不变」，被可证伪测试抓到并纠正，物理以反比为准。）
    5. 颅骨电导率单调性：σ_skull↑ → 靶区 E↑（颅骨阻抗降低，分流减少）。
    6. 逆设计闭式往返残差 ~0。
    7. 量级锚点：4 球模型 1mA → 皮层 ~0.2 V/m（文献 1mA 常用剂量量级）。
    8. 实践窗口 advisory 标注核对（2mA 非 FDA 硬阈值，如实标）。

诚实边界：
    * 闭式与球谐共享「点电极 Neumann 边界 c_l=I(2l+1)/(2πR²)」这一物理假设，
      故交叉核对锁住的是「线性系统求解正确」，不是电极模型本身（点电极近场
      本就不适用，已如实列为边界）。
    * 不做「球谐 vs FEM 真实头模」对照——那需外部 FEM 数据与网格，本仓库没有，
      不编造。四球电导率/半径来自已核实文献，作为输入锚点而非验证结论。

免责：科研验证工具，非临床。所有结论可复现（固定参数，无随机性）。
"""
from __future__ import annotations

import math

from . import electric_field as ef


def single_sphere_vs_closed_form(*, sigma: float = 0.33, R_cm: float = 9.0,
                                 current_ma: float = 1.0,
                                 n_harmonics: int = 48,
                                 radii_cm: tuple[float, ...] = (0.5, 2.0, 4.5, 7.0)) -> dict:
    """单层均匀球：球谐求解器 vs 解析闭式，扩散区逐点交叉核对。

    判据：r≤7cm（u≤0.78，远离点电极 r→R 奇异）相对误差 <1e-3。
    """
    layers = (ef.TissueLayer(R_cm, sigma, "uniform"),)
    coeffs = ef.solve_harmonics(layers, current_ma, n_harmonics=n_harmonics)
    rows = []
    for r in radii_cm:
        e_num, _ = ef.field_and_potential_on_axis(layers, r, coeffs)
        e_closed = ef.uniform_sphere_closed_form(r, sigma, R_cm, current_ma)
        rel = (e_num - e_closed) / e_closed
        rows.append({
            "r_cm": r,
            "harmonic_v_m": round(e_num, 8),
            "closed_form_v_m": round(e_closed, 8),
            "relative_error": round(rel, 8),
        })
    worst = max(abs(x["relative_error"]) for x in rows)
    return {
        "passed": bool(worst < 1e-3),
        "criterion": "relative error <1e-3 at r<=7cm (away from point-electrode singularity)",
        "method": "2N×2N spherical-harmonic linear solve vs analytic closed form (series sum)",
        "sigma_s_m": sigma, "R_cm": R_cm, "current_ma": current_ma,
        "n_harmonics": n_harmonics,
        "worst_relative_error": round(worst, 8),
        "rows": rows,
        "note": "共享点电极 Neumann 边界假设；此项锁住线性系统求解正确，不验证电极模型本身",
    }


def near_field_truncation(*, sigma: float = 0.33, R_cm: float = 9.0,
                          current_ma: float = 1.0, r_cm: float = 8.5) -> dict:
    """近场截断：r→R 处闭式发散、有限阶球谐截断，偏差随阶数减小。

    如实报告这是截断误差，不作 pass/fail（近场点电极本就不适用）。
    """
    layers = (ef.TissueLayer(R_cm, sigma, "uniform"),)
    e_closed = ef.uniform_sphere_closed_form(r_cm, sigma, R_cm, current_ma)
    rows = []
    for n in (12, 24, 48, 96):
        coeffs = ef.solve_harmonics(layers, current_ma, n_harmonics=n)
        e_num, _ = ef.field_and_potential_on_axis(layers, r_cm, coeffs)
        rows.append({
            "n_harmonics": n,
            "harmonic_v_m": round(e_num, 6),
            "closed_form_v_m": round(e_closed, 6),
            "relative_error": round((e_num - e_closed) / e_closed, 6),
        })
    monotone = all(abs(a["relative_error"]) >= abs(b["relative_error"])
                   for a, b in zip(rows, rows[1:]))
    return {
        "role": "report truncation error (NOT a pass/fail gate; near-field point electrode invalid)",
        "r_cm": r_cm, "u": round(r_cm / R_cm, 4),
        "rows": rows,
        "error_shrinks_with_more_harmonics": bool(monotone),
        "interpretation": (
            "点电极在 r→R 有奇异，有限阶球谐必然截断；闭式同处发散。"
            "靠近头皮/电极正下方的场强本就不可信（引擎已 WARN），此处只量化截断行为"
        ),
    }


def current_linearity(*, target_depth_cm: float = 2.2) -> dict:
    """E ∝ I（可证伪）：I 翻倍，E 必翻倍，电势也翻倍。"""
    pairs = []
    for i in (0.5, 1.0, 2.0):
        r = ef.dose_field(i, target_depth_cm=target_depth_cm)
        pairs.append({"current_ma": i, "field_v_m": round(r.field_v_m, 8),
                      "potential_v": round(r.potential_v, 8)})
    e1 = ef.dose_field(1.0, target_depth_cm=target_depth_cm).field_v_m
    e2 = ef.dose_field(2.0, target_depth_cm=target_depth_cm).field_v_m
    ok = abs(e2 / e1 - 2.0) < 1e-9
    return {"cases": pairs, "passed": bool(ok),
            "field_ratio_2x": round(e2 / e1, 8),
            "note": "准静态线性方程 + 固定电流源 ⇒ E 与 I 严格线性"}


def conductivity_scale_invariance(*, target_depth_cm: float = 2.2, k: float = 2.0) -> dict:
    """电导率整体缩放反比性：所有层 σ×k → E 变为 1/k（非平凡对称性检验）。

    固定总电流 I ⇒ 表面 Neumann 边界 -σ∂φ/∂r=J（J 由 I 定、与 σ 无关）。
    σ→kσ ⇒ 方程∇·(σ∇φ)=0 中 k 约掉、边界给出 ∂φ' = ∂φ/k ⇒ φ'=φ/k ⇒ E'=E/k。
    等价电路直觉：同电流下整体电阻降 k 倍 ⇒ 电压与场都降 k 倍。
    若边界条件误写成 Dirichlet（固定电势）会得到「E 不变」的错误结论。
    """
    scaled = tuple(ef.TissueLayer(l.outer_radius_cm, l.conductivity_s_m * k, l.name)
                   for l in ef.HEAD_4SPHERE)
    e1 = ef.dose_field(1.0, layers=ef.HEAD_4SPHERE, target_depth_cm=target_depth_cm).field_v_m
    e2 = ef.dose_field(1.0, layers=scaled, target_depth_cm=target_depth_cm).field_v_m
    ratio = e2 / e1
    rel = abs(ratio - 1.0 / k) / (1.0 / k)
    return {"k": k, "field_original_v_m": round(e1, 8), "field_scaled_v_m": round(e2, 8),
            "expected_ratio": round(1.0 / k, 8), "measured_ratio": round(ratio, 8),
            "relative_error": round(rel, 10), "passed": bool(rel < 1e-9),
            "note": "固定电流源 + 准静态 ⇒ 整体 σ×k 使 E×1/k；检验 Neumann 边界实现"}


def skull_monotonicity() -> dict:
    """颅骨电导率单调性：σ_skull↑ → 靶区 E↑（颅骨阻抗降低）。"""
    s = ef.skull_conductivity_sensitivity()
    return {
        "passed": bool(s["monotone_increasing"] and s["field_max_min_ratio"] > 1.5),
        "field_max_min_ratio": s["field_max_min_ratio"],
        "monotone_increasing": s["monotone_increasing"],
        "cases": s["cases"],
        "interpretation": s["interpretation"],
        "literal_anchor": s["literal_anchor"],
    }


def inverse_roundtrip() -> dict:
    """逆设计闭式往返残差 ~0。"""
    rows = []
    for target in (0.1, 0.2, 0.4):
        inv = ef.inverse_current(target)
        rows.append({
            "target_v_m": target,
            "required_ma": round(inv["required_current_ma"], 8),
            "recomputed_v_m": round(inv["recomputed_field_v_m"], 8),
            "residual_rel": round(inv["residual_rel"], 10),
        })
    worst = max(r["residual_rel"] for r in rows)
    return {"cases": rows, "worst_residual_rel": worst, "passed": bool(worst < 1e-9),
            "note": "E∝I 线性 ⇒ 闭式逆解精确，正向重算残差应 ~0"}


def magnitude_anchor() -> dict:
    """量级锚点：四球 1mA → 皮层 ~0.2 V/m（文献 1mA 常用剂量量级）。"""
    r = ef.dose_field(1.0, target_depth_cm=2.2)
    in_range = 0.1 <= r.field_v_m <= 0.5
    return {
        "field_v_m": round(r.field_v_m, 6),
        "field_mv_mm": round(r.field_mv_mm, 6),
        "expected_literature_range_v_m": [0.1, 0.5],
        "in_range": bool(in_range),
        "passed": bool(in_range),
        "source": (
            "tDCS 1-2mA 实践窗；1mA 皮层场 ~0.2 V/m 量级（Rush-Driscoll 球模型"
            "+ 文献 tDCS 仿真）；四球电导率来自 PMC9927519（已核实 2026-11-19）"
        ),
        "note": "量级锚点非精确判据；真实头模个体差异可达 2.1×",
    }


def advisory_labeling() -> dict:
    """实践窗口 advisory 标注核对：2mA 非 FDA 硬阈值，如实标。"""
    over = ef.dose_field(3.0, target_depth_cm=2.2)
    has_advisory = any("advisory" in w or "FDA" in w for w in over.warnings)
    has_over = any("3.00mA" in w or "实践窗口" in w for w in over.warnings)
    inv_over = ef.inverse_current(1.0)  # 1 V/m 需 ~5mA > 2mA，应标 within=False
    return {
        "practice_current_max_ma": ef.PRACTICE_CURRENT_MAX_MA,
        "high_current_warns": list(over.warnings),
        "advisory_labeled": bool(has_advisory),
        "over_window_flagged": bool(has_over),
        "inverse_flags_over_window": not inv_over["within_practice_window"],
        "passed": bool(has_advisory and has_over and not inv_over["within_practice_window"]),
        "note": "tES 无 FDA 硬阈值；2mA 是二十年实践窗口，冒充硬阈值即撒谎",
    }


# ---------------------------------------------------------------------------
# 汇总报告
# ---------------------------------------------------------------------------

def build_report() -> dict:
    """跑全套验证，产出可对外引用的机器可读报告（固定参数，可复现）。"""
    checks = {
        "single_sphere_vs_closed_form": single_sphere_vs_closed_form(),
        "current_linearity": current_linearity(),
        "conductivity_scale_invariance": conductivity_scale_invariance(),
        "skull_monotonicity": skull_monotonicity(),
        "inverse_roundtrip": inverse_roundtrip(),
        "magnitude_anchor": magnitude_anchor(),
        "advisory_labeling": advisory_labeling(),
    }
    gated = list(checks.keys())   # 全部为真正判据门（与光/声的 surrogate 不同，
                                  # 电模态的近场截断单列为 report，不进门）
    passed = all(checks[k].get("passed", False) for k in gated)
    return {
        "version": "L2-electric-alpha",
        "overall_passed": bool(passed),
        "gated_checks": gated,
        "reported_only": ["near_field_truncation"],
        "checks": checks,
        "near_field_truncation": near_field_truncation(),
        "disclaimer": ef.DISCLAIMER,
    }


def save_report(path: str = "data/es_validation_report.json") -> str:
    """写报告到 JSON（不 round 到丢失精度，固定内容可复现）。"""
    import json
    from pathlib import Path
    report = build_report()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return str(p)


def render_report(report: dict | None = None) -> str:
    """人类可读 Markdown。"""
    r = report or build_report()
    c = r["checks"]
    labels = {
        "single_sphere_vs_closed_form": "单层球谐 vs 解析闭式",
        "current_linearity": "电流线性 E∝I",
        "conductivity_scale_invariance": "电导率缩放反比性 (E×1/k)",
        "skull_monotonicity": "颅骨电导率单调性",
        "inverse_roundtrip": "逆设计往返残差",
        "magnitude_anchor": "量级锚点 ~0.2V/m",
        "advisory_labeling": "实践窗口 advisory 标注",
    }
    lines = [
        "## L2 电场剂量场验证报告（可复现）",
        "",
        f"**总体：{'PASS' if r['overall_passed'] else 'FAIL'}**"
        f"（判据门 {len(r['gated_checks'])} 项）",
        "",
        "| 核对项 | 方法 | 结果 |",
        "|---|---|---|",
    ]
    for k in r["gated_checks"]:
        lines.append(f"| {labels[k]} | 独立路径交叉核对/可证伪性质 | "
                     f"{'PASS' if c[k].get('passed') else 'FAIL'} |")
    sph = c["single_sphere_vs_closed_form"]
    lines += [
        "",
        f"单层球谐 vs 闭式：最坏相对误差 {sph['worst_relative_error']}"
        f"（判据 r≤7cm <1e-3）。球心场解析值 "
        f"{ef.uniform_sphere_closed_form(0.01, 0.33, 9.0, 1.0):.4f} V/m。",
        f"颅骨电导率 0.0042→0.03 S/m：靶区场 "
        f"{c['skull_monotonicity']['cases'][0]['field_v_m']:.4f}→"
        f"{c['skull_monotonicity']['cases'][-1]['field_v_m']:.4f} V/m，"
        f"变化 {c['skull_monotonicity']['field_max_min_ratio']}×"
        f"（单调={c['skull_monotonicity']['monotone_increasing']}）。",
        "",
        "> 近场截断（r→R）单列为 report 不进门：点电极近表面奇异，有限阶球谐必然截断，"
        "闭式同处发散——引擎已 WARN 该区域不可信，不掩盖也不调参。",
        "",
        "> tES 无 FDA 硬阈值；2mA 标为实践窗口（advisory），不冒充监管线。",
        "",
        f"> {r['disclaimer']}",
    ]
    return "\n".join(lines)
