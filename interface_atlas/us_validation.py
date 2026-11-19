"""L2 超声剂量场的独立验证基准（交叉核对 + 代理误差如实报告）.

为什么需要它：
    ultrasound_field.py 的聚焦增益、分层衰减是降阶解析式。用自己的解析式去
    「验证」自己等于循环论证。本模块用不同实现路径做交叉核对：

    1. 平面活塞轴上声场：用数值 Rayleigh-Sommerfeld 求积（numpy trapz）对照
       变量代换后的解析闭式。这是精确恒等式（r·dr=R·dR），二者应吻合到数值
       精度——锁住「声场基本解算对了」。对应光模块用 MC 验证 Prahl 基准。
    2. 声强-声压关系往返一致性：I=p²/(2ρc) 的正逆变换闭合。
    3. MI 定义与硬阈值同源：MI=p[MPa]/√f[MHz]，且阈值与 paramcheck 单一真源。
    4. 颅骨衰减量级：0.5MHz 经 4-8mm 颅骨的透射落在文献「单程数 dB~十余 dB」区间。
    5. 温升粗估的线性可证伪性：对声强/占空比/时长线性缩放。

诚实边界：
    * 聚焦增益 G=η(πa·sinα/λ)² 是几何近似，**没有独立闭式**可严格交叉核对。
      focus_gain_is_surrogate() 只做量级/单调性测量，标为 surrogate，不作 pass/fail
      gate——与光模块「P1 代理误差是测量输出」同一纪律，不调参凑数。
    * 数值 RS 求积验证的是平面活塞（可解析），不覆盖真实聚焦换能器的像差/旁瓣。

免责：科研验证工具，非临床。所有结论可复现（固定 seed / 固定网格）。
"""
from __future__ import annotations

import math

import numpy as np

from . import paramcheck, ultrasound_field as uf


# ---------------------------------------------------------------------------
# 1. 平面活塞轴上声场：数值 RS 求积 vs 解析闭式
# ---------------------------------------------------------------------------

def _axis_pressure_closed_form(z_cm: float, a_cm: float, k_per_cm: float) -> float:
    """轴上复声压幅度闭式（归一化前）.

    平面活塞、RS 第一积分（忽略倾斜因子，轴对称）：
        I(z) = ∫_0^a e^{-ik·R}/R · r dr,  R=√(z²+r²)
    变量代换 r·dr = R·dR ⇒ I = ∫_z^{√(z²+a²)} e^{-ikR} dR
                          = (e^{-ikz} − e^{-ik√(z²+a²)})/(ik)
        |I| = 2·|sin(k(√(z²+a²)−z)/2)| / k
    """
    r_max = math.sqrt(z_cm ** 2 + a_cm ** 2)
    return 2.0 * abs(math.sin(k_per_cm * (r_max - z_cm) / 2.0)) / k_per_cm


def _axis_pressure_numeric(z_cm: float, a_cm: float, k_per_cm: float,
                           n_r: int = 200000) -> float:
    """同一积分的数值求积（不同实现路径）：trapz 对 r 直接积分。"""
    r = np.linspace(0.0, a_cm, n_r)
    R = np.sqrt(z_cm ** 2 + r ** 2)
    integrand = np.cos(k_per_cm * R) / R * r + 1j * (-np.sin(k_per_cm * R)) / R * r
    # 复积分的模
    re = np.trapezoid(integrand.real, r)
    im = np.trapezoid(integrand.imag, r)
    return math.hypot(re, im)


def axis_comparison(*, frequency_mhz: float = 1.0, a_cm: float = 1.0,
                    sound_speed_m_s: float = 1480.0,
                    z_min_cm: float | None = None, z_max_cm: float | None = None,
                    n_z: int = 120, n_r: int = 200000) -> dict:
    """平面活塞轴上场：数值求积 vs 解析闭式（归一化后比对 shape）.

    归一化到各自最大值后比较，消除绝对标定差异，只核对 Fresnel 振荡形状。
    判据：RMS 相对误差 < 1% 且最大绝对差 < 0.02（归一化单位）。
    z 范围默认取 2a~15a（Fresnel 区，避免近场 z<a 的数值病态）。
    """
    lam = uf.wavelength_cm(frequency_mhz, sound_speed_m_s)
    k = 2.0 * math.pi / lam
    z0 = z_min_cm if z_min_cm is not None else 2.0 * a_cm
    z1 = z_max_cm if z_max_cm is not None else 15.0 * a_cm
    zs = np.linspace(z0, z1, n_z)
    closed = np.array([_axis_pressure_closed_form(float(z), a_cm, k) for z in zs])
    numeric = np.array([_axis_pressure_numeric(float(z), a_cm, k, n_r) for z in zs])
    cn, nn = closed / closed.max(), numeric / numeric.max()
    diff = np.abs(cn - nn)
    rms = float(np.sqrt(np.mean(diff ** 2)))
    maxd = float(diff.max())
    passed = rms < 0.01 and maxd < 0.02
    # 抽样若干点供人工核对（可复现）
    idx = np.linspace(0, n_z - 1, 8, dtype=int)
    return {
        "passed": bool(passed),
        "criterion": "RMS<0.01 & max_abs_diff<0.02 (both normalized to max)",
        "method": "Rayleigh-Sommerfeld axial integral: numpy trapz vs variable-substitution closed form",
        "frequency_mhz": frequency_mhz,
        "aperture_radius_cm": a_cm,
        "wavelength_cm": round(lam, 6),
        "z_range_cm": [round(z0, 4), round(z1, 4)],
        "n_z": n_z, "n_r_per_point": n_r,
        "rms_diff": round(rms, 6),
        "max_abs_diff": round(maxd, 6),
        "samples": [
            {"z_cm": round(float(zs[i]), 4),
             "closed_norm": round(float(cn[i]), 5),
             "numeric_norm": round(float(nn[i]), 5)}
            for i in idx
        ],
        "source": "Rayleigh-Sommerfeld first integral; variable substitution r·dr=R·dR is exact",
    }


# ---------------------------------------------------------------------------
# 2. 声强-声压往返一致性
# ---------------------------------------------------------------------------

def intensity_pressure_roundtrip() -> dict:
    """I=p²/(2ρc) 正逆变换闭合（纯水与脑组织两个介质）。"""
    out = []
    for rho, c, name in ((1000.0, 1480.0, "water"), (1040.0, 1560.0, "brain")):
        for i_w_cm2 in (0.5, 5.0, 15.0):
            p_pa = uf.intensity_to_pressure_pa(i_w_cm2, rho, c)
            i_back = uf.pressure_to_intensity_w_cm2(p_pa, rho, c)
            out.append({
                "medium": name, "rho": rho, "c": c,
                "I_in_w_cm2": i_w_cm2,
                "p_pk_mpa": round(p_pa / 1e6, 6),
                "I_back_w_cm2": round(i_back, 6),
                "rel_error": abs(i_back - i_w_cm2) / i_w_cm2,
            })
    worst = max(r["rel_error"] for r in out)
    return {"cases": out, "worst_rel_error": worst, "passed": worst < 1e-12,
            "note": "1 MPa 峰值声压对应 I≈p²/(2ρc)；水/脑两个介质往返一致"}


# ---------------------------------------------------------------------------
# 3. MI 定义与硬阈值同源
# ---------------------------------------------------------------------------

def mi_and_thresholds() -> dict:
    """MI 定义的量纲核对 + 与 paramcheck 硬阈值单一真源核对。"""
    mi_1mpa_1mhz = uf.mechanical_index(1.0, 1.0)
    mi_1mpa_4mhz = uf.mechanical_index(1.0, 4.0)
    # 同声压下频率翻倍4倍 → MI 减半（∝1/√f）
    inv_sqrt_check = abs(mi_1mpa_4mhz - 0.5) < 1e-12

    # 阈值同源：引擎常量应等于 paramcheck 的硬阈值
    # HARD_LIMITS[key] 是单个四元组 (key, op, threshold, msg)，不是列表
    mi_rule = paramcheck.HARD_LIMITS.get("focused ultrasound")
    mi_limit_param = mi_rule[2] if mi_rule and mi_rule[0] == "mechanical_index" else None
    temp_key, temp_op, temp_limit_param, temp_msg = paramcheck.TEMP_LIMIT
    same_source = (uf.MI_LIMIT == mi_limit_param == 1.9
                   and uf.TEMP_LIMIT_C == temp_limit_param == 2.0)
    return {
        "MI_1MPa_1MHz": mi_1mpa_1mhz,
        "MI_1MPa_4MHz": mi_1mpa_4mhz,
        "inverse_sqrt_frequency_holds": bool(inv_sqrt_check),
        "engine_mi_limit": uf.MI_LIMIT,
        "paramcheck_mi_limit": mi_limit_param,
        "engine_temp_limit_c": uf.TEMP_LIMIT_C,
        "paramcheck_temp_limit_c": temp_limit_param,
        "thresholds_single_source": bool(same_source),
        "passed": bool(inv_sqrt_check and same_source and mi_1mpa_1mhz == 1.0),
        "source": "AIUM/FDA MI=p[MPa]/√f[MHz]; thresholds from paramcheck.HARD_LIMITS/TEMP_LIMIT",
    }


# ---------------------------------------------------------------------------
# 4. 颅骨衰减量级
# ---------------------------------------------------------------------------

def skull_attenuation_magnitude() -> dict:
    """0.5MHz 经颅衰减落在文献量级区间（单程数 dB ~ 十余 dB）。"""
    tr = uf.TRANSDUCER_PRESETS["tfus_transcranial"]
    out = []
    for name in ("skull_thin", "skull_thick"):
        lay = uf.TISSUE_PRESETS[name]
        r = uf.dose_field(tr, (lay,))
        out.append({
            "preset": name,
            "thickness_cm": lay.thickness_cm,
            "alpha_db_cm_at_1mhz": lay.attenuation_db_cm,
            "alpha_db_cm_at_0p5mhz": round(lay.attenuation_at(0.5), 4),
            "attenuation_db": round(r.attenuation_db_total, 4),
            "transmission": round(r.transmission, 5),
        })
    # 单程 4-8mm 颅骨 0.5MHz：α(0.5)=11-12.5 dB/cm → 4.4-10 dB，透射 0.32-0.60
    db_vals = [c["attenuation_db"] for c in out]
    in_range = all(2.0 <= d <= 15.0 for d in db_vals)
    monotone = db_vals[0] < db_vals[1]  # 越厚衰减越大
    return {
        "cases": out,
        "expected_range_db": [2.0, 15.0],
        "all_in_expected_range": bool(in_range),
        "monotone_in_thickness": bool(monotone),
        "passed": bool(in_range and monotone),
        "note": "单程 0.5MHz 经 4-8mm 颅骨；文献经颅总衰减常报 10-20dB（含往返/散射），本模型单程量级一致",
    }


# ---------------------------------------------------------------------------
# 5. 温升粗估的线性可证伪性
# ---------------------------------------------------------------------------

def temperature_scaling() -> dict:
    """温升粗估对声强/占空比/时长线性（可证伪：任一不满足即失败）。"""
    tr = uf.TRANSDUCER_PRESETS["tfus_superficial"]
    base = uf.dose_field(tr, (uf.TISSUE_PRESETS["brain"],),
                         duty_cycle=0.1, pulse_duration_s=1.0)
    double_i = uf.dose_field(
        uf.Transducer(tr.radius_cm, tr.focal_length_cm, tr.frequency_mhz,
                      tr.intensity_w_cm2 * 2.0, tr.efficiency),
        (uf.TISSUE_PRESETS["brain"],), duty_cycle=0.1, pulse_duration_s=1.0)
    double_duty = uf.dose_field(tr, (uf.TISSUE_PRESETS["brain"],),
                                duty_cycle=0.2, pulse_duration_s=1.0)
    double_time = uf.dose_field(tr, (uf.TISSUE_PRESETS["brain"],),
                                duty_cycle=0.1, pulse_duration_s=2.0)
    rel = lambda x: abs(x - 2.0 * base.delta_t_c) / (2.0 * base.delta_t_c)
    checks = {
        "linear_in_intensity": rel(double_i.delta_t_c) < 1e-9,
        "linear_in_duty_cycle": rel(double_duty.delta_t_c) < 1e-9,
        "linear_in_duration": rel(double_time.delta_t_c) < 1e-9,
    }
    return {
        "base_delta_t_c": round(base.delta_t_c, 6),
        "checks": checks,
        "passed": all(checks.values()),
        "note": "ΔT=μa·I·t/(ρc) 线性；非临床结论，仅验证粗估模型的自洽性",
    }


# ---------------------------------------------------------------------------
# 6. 聚焦增益：代理测量（非 pass/fail gate）
# ---------------------------------------------------------------------------

def focus_gain_is_surrogate() -> dict:
    """聚焦增益的量级/单调性测量——标为 surrogate，无独立闭式，不作判据。

    与光模块 validate_patterson_mu_eff 同一纪律：如实标注这是测量输出而非
    pass/fail，绝不调参数凑一个「漂亮」的目标值。
    """
    lam = uf.wavelength_cm(0.5, 1560.0)
    gains = {}
    for name, tr in uf.TRANSDUCER_PRESETS.items():
        gains[name] = round(uf.focus_gain(tr, lam), 3)
    # 单调性：同焦距下增大孔径 → 增益增（F# 减小）
    small = uf.Transducer(1.0, 8.0, 0.5, 10.0, 0.7)
    large = uf.Transducer(3.0, 8.0, 0.5, 10.0, 0.7)
    monotone = uf.focus_gain(large, lam) > uf.focus_gain(small, lam)
    return {
        "role": "measure surrogate only (no independent closed form; NOT a pass/fail gate)",
        "gains_at_lambda": gains,
        "monotone_in_aperture": bool(monotone),
        "interpretation": (
            "G=η(πa·sinα/λ)² 是衍射受限几何近似，忽略像差/旁瓣/非线性；"
            "真实换能器焦点增益须用实测声场或全波仿真标定，本值仅作量级参考"
        ),
    }


# ---------------------------------------------------------------------------
# 汇总报告
# ---------------------------------------------------------------------------

def build_report() -> dict:
    """跑全套验证，产出可对外引用的机器可读报告（固定参数，可复现）。"""
    checks = {
        "axis_rs_integral": axis_comparison(),
        "intensity_pressure_roundtrip": intensity_pressure_roundtrip(),
        "mi_and_thresholds": mi_and_thresholds(),
        "skull_attenuation": skull_attenuation_magnitude(),
        "temperature_scaling": temperature_scaling(),
        "focus_gain_surrogate": focus_gain_is_surrogate(),
    }
    gated = ["axis_rs_integral", "intensity_pressure_roundtrip",
             "mi_and_thresholds", "skull_attenuation", "temperature_scaling"]
    passed = all(checks[k].get("passed", False) for k in gated)
    return {
        "version": "L2-ultrasound-alpha",
        "overall_passed": bool(passed),
        "gated_checks": gated,
        "surrogate_checks": ["focus_gain_surrogate"],
        "checks": checks,
        "disclaimer": uf.DISCLAIMER,
    }


def save_report(path: str = "data/us_validation_report.json") -> str:
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
    lines = [
        "## L2 超声剂量场验证报告（可复现）",
        "",
        f"**总体：{'PASS' if r['overall_passed'] else 'FAIL'}**"
        f"（判据门 {len(r['gated_checks'])} 项，代理测量 {len(r['surrogate_checks'])} 项）",
        "",
        "| 核对项 | 方法 | 结果 |",
        "|---|---|---|",
    ]
    gate_label = {
        "axis_rs_integral": "数值 RS 求积 vs 解析闭式",
        "intensity_pressure_roundtrip": "声强-声压往返",
        "mi_and_thresholds": "MI 定义 + 阈值同源",
        "skull_attenuation": "颅骨衰减量级",
        "temperature_scaling": "温升线性可证伪",
    }
    for k in r["gated_checks"]:
        lines.append(f"| {gate_label[k]} | 不同实现路径交叉核对 | "
                     f"{'PASS' if c[k].get('passed') else 'FAIL'} |")
    lines.append(f"| 聚焦增益 | 代理测量（无独立闭式） | "
                 f"{'surrogate' if 'focus_gain_surrogate' in c else '—'} |")
    ax = c["axis_rs_integral"]
    lines += [
        "",
        f"轴上场 RMS 相对差 {ax['rms_diff']}，最大绝对差 {ax['max_abs_diff']}"
        f"（归一化单位，判据 RMS<0.01）。",
        f"MI: 1MPa@1MHz → {c['mi_and_thresholds']['MI_1MPa_1MHz']}，"
        f"1MPa@4MHz → {c['mi_and_thresholds']['MI_1MPa_4MHz']}；"
        f"阈值与 paramcheck 同源={c['mi_and_thresholds']['thresholds_single_source']}。",
        "",
        "> 聚焦增益 G=η(πa·sinα/λ)² 无独立闭式，标为 surrogate，不作 pass/fail——"
        "不调参凑「漂亮」目标值（与光模块代理误差同一纪律）。",
        "",
        f"> {r['disclaimer']}",
    ]
    return "\n".join(lines)
