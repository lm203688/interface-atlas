"""L2 中阶：把 P1 代理点估计升级为带覆盖保证的剂量区间（split conformal）.

为什么不是「拟合误差曲面」：
    实测（见 data/calibration_report.json）显示 P1 扩散解相对独立 MC 的 μeff 比
    r = μeff_MC / μeff_P1 ≈ 1.156 是稳定的**系统偏差**，而随机散布仅约 2%，随
    (μa, μs', g) 无显著趋势。把这点噪声拟合成「随参数变化的曲面」是自欺——它
    拟合的是 MC 采样噪声，换个 seed 就崩。正确做法是把系统偏差单独估计，再给
    偏差残差一个**有限样本有效**的覆盖区间。

方法（split conformal prediction，无分布假设）:
    1. fit 子集：估全局乘性偏差 b = median(r)（对数空间对称化）。
    2. calib 子集：score_i = |ln(r_i / b)|，取分位数 q_{1-α} = s_(ceil((n+1)(1-α)))。
       理论保证：对 exchangeable 的 (fit, calib, test) 样本，
       P( r_new ∈ [b·e^-q, b·e^q] ) >= 1-α（有限样本，不依赖 r 的分布）。
    3. test 子集：独立实测覆盖率，验证保证成立（可证伪）。
    三段参数点、MC seed 全部分离，杜绝数据泄漏。

从 μeff 区间到剂量区间（精确，非一阶近似）:
    μeff 对 (μa, μs') 同尺度缩放是齐次的：μeff(s·μa, s·μs') = s·μeff(μa, μs')。
    因此把组织光学参数整体乘 r 即可令 μeff 乘 r，且保持输运反照率 a' 与 g 不变
    （与校准域一致）。剂量区间 = 用 r_lo / r_hi 缩放后重算的剂量场，表面仍归一
    为 1。r>1 表示真实衰减更快 → 深部剂量比 P1 点估计更低。

诚实边界:
    * 覆盖保证只对**扩散近似适用域内**的参数成立（弱吸收、扩散区）；域外由
      light_field.dose_field 的 WARN 拦截，本模块不外推保证。
    * 区间宽度反映 MC 采样噪声 + 参数域内残差，不含系统偏差估计误差的尾部；
      b 用稳健中位数，n_fit 已在报告中给出。
    * 科研降阶模型，非临床处方。

来源: split conformal 的标准有限样本结论见 Angelopoulos & Bates 2023
      "Conformal Prediction: A Gentle Introduction", Now Publishers（教程）。
      生物组织光学 MC 金标准见 Wang, Jacques & Zheng 1995 MCML。
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

from .light_field import TissueOptics, dose_field

REPORT_PATH = Path(__file__).resolve().parent.parent / "data" / "calibration_report.json"

# 校准域（覆盖组织光学参数的真实量级，见 TISSUE_PRESETS 与文献 NIR 典型值）
# 实测 r 在 μa 0.2–3.0 / μs' 40–100 稳定在 1.13–1.16，与低吸收域一致——系统偏差
# 与参数弱相关，因此采样域必须覆盖真实用例，覆盖保证才有意义。
DOMAIN = {
    "mu_a_range": [0.10, 4.00],       # 1/cm
    "mu_s_prime_range": [30.0, 130.0], # 1/cm
    "g_range": [0.60, 0.95],
    "note": "弱吸收扩散区（μa<<μs'、a'≈1）；采样域覆盖 NIR 组织光学参数真实量级",
}

ALPHAS = (0.20, 0.10, 0.05)
DEFAULT_ALPHA = 0.10


# ---------------------------------------------------------------------------
# ground-truth 观测（独立 MC）
# ---------------------------------------------------------------------------

def mu_eff_ratio(
    mu_a: float,
    mu_s_prime: float,
    g: float,
    *,
    n_photons: int = 80_000,
    seed: int = 20261118,
    max_depth_cm: float = 4.0,
) -> float | None:
    """用独立 MC 拟合 μeff，返回 r = μeff_MC / μeff_P1。

    拟合窗口取 5×mfp' 到 2/μeff_P1 的扩散区，避免近源曲率与深部噪声。
    可拟合点不足（强吸收）时返回 None，由 build_calibration 丢弃。
    """
    from .mc_light import run_mc, Layer  # 延迟导入，避免循环依赖

    if mu_a <= 0 or mu_s_prime <= 0 or not (0 <= g < 1):
        raise ValueError("需 mu_a>0, mu_s_prime>0, 0<=g<1")
    mu_s = mu_s_prime / (1.0 - g)
    layers = (Layer(float(mu_a), float(mu_s), float(g), 20.0),)
    mu_p1 = math.sqrt(3.0 * mu_a * (mu_a + mu_s_prime))
    # 高 μs' 穿透浅：max_depth 随 μeff 自适应，避免深部 bin 全空导致拟合点不足
    md = min(max_depth_cm, 4.0 / mu_p1)
    res = run_mc(layers, n_photons=n_photons, seed=seed, n_bins=400,
                 max_depth_cm=md)
    edges = np.asarray(res.bin_edges_cm)
    depth = 0.5 * (edges[:-1] + edges[1:])
    absd = np.asarray(res.depth_absorption)
    # 拟合窗口 3·mfp' 到 2.5/μeff（扩散区），信号阈值放宽到峰值 1e-5
    z0, z1 = 3.0 / mu_s_prime, min(2.5 / mu_p1, md)
    mask = (depth >= z0) & (depth <= z1) & (absd > absd.max() * 1e-5)
    if mask.sum() < 8:
        return None
    coef = np.polyfit(depth[mask], np.log(absd[mask]), 1)
    return float(-coef[0] / mu_p1)


def _sample_params(n: int, seed: int) -> list[tuple[float, float, float]]:
    """在对数光学参数域内均匀采样参数点（exchangeable）。"""
    rng = np.random.default_rng(seed)
    lo_a, hi_a = DOMAIN["mu_a_range"]
    lo_s, hi_s = DOMAIN["mu_s_prime_range"]
    lo_g, hi_g = DOMAIN["g_range"]
    mu_a = 10 ** rng.uniform(math.log10(lo_a), math.log10(hi_a), n)
    mu_sp = 10 ** rng.uniform(math.log10(lo_s), math.log10(hi_s), n)
    g = rng.uniform(lo_g, hi_g, n)
    return [(float(a), float(s), float(gg)) for a, s, gg in zip(mu_a, mu_sp, g)]


def _conformal_q(scores: list[float], alpha: float) -> float:
    """分位数 q_{1-alpha} = s_(ceil((n+1)(1-alpha)))；超界返回 inf。"""
    s = sorted(float(x) for x in scores)
    n = len(s)
    if not n:
        return float("inf")
    k = math.ceil((n + 1) * (1.0 - alpha))
    return s[k - 1] if k <= n else float("inf")


def _coverage(ratios: list[float], b: float, q: float) -> float:
    if not ratios or not math.isfinite(q):
        return 0.0
    return float(np.mean([abs(math.log(r / b)) <= q for r in ratios]))


# ---------------------------------------------------------------------------
# 报告构建（一次性，较慢；结果固化到 data/calibration_report.json）
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CalibrationReport:
    bias: float                       # b，全局乘性偏差
    q_by_alpha: dict[str, float]      # alpha -> q
    coverage_by_alpha: dict[str, float]  # 独立 test 集实测覆盖率
    falsification: dict               # 反证结果
    n_fit: int
    n_calib: int
    n_test: int
    n_photons: int
    r_fit: tuple[float, ...]
    r_calib: tuple[float, ...]
    r_test: tuple[float, ...]
    domain: dict
    seeds: dict

    def as_dict(self) -> dict:
        return {
            "version": "1.0",
            "method": "split conformal (fit/calib/test), finite-sample marginal coverage",
            "bias_b": round(self.bias, 6),
            "bias_interpretation": (
                "mu_eff_MC = b * mu_eff_P1；b>1 表示 P1 低估衰减、高估深部剂量，"
                "真实靶区剂量比 P1 点估计更低"
            ),
            "q_by_alpha": {k: round(v, 6) for k, v in self.q_by_alpha.items()},
            "coverage_by_alpha": self.coverage_by_alpha,
            "falsification": self.falsification,
            "sample_sizes": {"fit": self.n_fit, "calib": self.n_calib,
                             "test": self.n_test, "photons_each": self.n_photons},
            "r_fit": [round(x, 6) for x in self.r_fit],
            "r_calib": [round(x, 6) for x in self.r_calib],
            "r_test": [round(x, 6) for x in self.r_test],
            "domain": self.domain,
            "seeds": self.seeds,
        }


def build_calibration(
    *,
    n_fit: int = 8,
    n_calib: int = 40,
    n_test: int = 24,
    n_photons: int = 80_000,
    seed: int = 20261118,
) -> CalibrationReport:
    """跑独立 MC 构建校准报告。三段样本参数点与 seed 全分离。"""
    if min(n_fit, n_calib, n_test) < 4:
        raise ValueError("每段至少 4 个参数点")
    # conformal 分位数 s_(ceil((n+1)(1-α))) 需 ceil(...)<=n，否则区间为 inf。
    # 对最小 alpha，要求 (n_calib+1)(1-alpha) <= n_calib ⇒ n_calib >= (1-alpha)/alpha。
    min_alpha = min(ALPHAS)
    if (n_calib + 1) * (1.0 - min_alpha) > n_calib:
        raise ValueError(
            f"n_calib={n_calib} 不足以支持 alpha={min_alpha}；"
            f"至少需要 {math.ceil((1-min_alpha)/min_alpha)} 个有效 calib 点"
        )
    # 三段用不同采样 seed 与不同 MC seed 起点，保证 exchangeable 且无泄漏
    params_fit = _sample_params(n_fit, seed + 1)
    params_calib = _sample_params(n_calib, seed + 2)
    params_test = _sample_params(n_test, seed + 3)

    def _run(params, offset):
        out = []
        for i, p in enumerate(params):
            r = mu_eff_ratio(*p, n_photons=n_photons, seed=seed + offset + i)
            if r is not None:
                out.append(r)
        return out

    r_fit = _run(params_fit, 100)
    r_calib = _run(params_calib, 300)
    r_test = _run(params_test, 600)
    if min(len(r_fit), len(r_calib), len(r_test)) < 4:
        raise RuntimeError("有效样本不足，无法构建校准报告")

    b = float(np.median(r_fit))
    scores_fit_calib = [abs(math.log(r / b)) for r in r_calib]
    q_by_alpha: dict[str, float] = {}
    cov_by_alpha: dict[str, float] = {}
    for a in ALPHAS:
        q = _conformal_q(scores_fit_calib, a)
        q_by_alpha[str(a)] = q
        cov_by_alpha[str(a)] = _coverage(r_test, b, q)

    # 反证：把 alpha=0.10 的区间宽度砍半，覆盖率必须显著下降
    q10 = q_by_alpha["0.1"]
    cov_full = _coverage(r_test, b, q10)
    cov_half = _coverage(r_test, b, q10 / 2.0)
    falsification = {
        "alpha": 0.10,
        "coverage_full_interval": round(cov_full, 4),
        "coverage_half_interval": round(cov_half, 4),
        "delta": round(cov_full - cov_half, 4),
        "passed": bool(cov_full - cov_half >= 0.10),
        "meaning": "区间砍半后覆盖率显著下降 ⇒ 区间宽度确实在起作用，非恒真工具",
    }
    return CalibrationReport(
        bias=b, q_by_alpha=q_by_alpha, coverage_by_alpha=cov_by_alpha,
        falsification=falsification, n_fit=len(r_fit), n_calib=len(r_calib),
        n_test=len(r_test), n_photons=n_photons,
        r_fit=tuple(r_fit), r_calib=tuple(r_calib), r_test=tuple(r_test),
        domain=dict(DOMAIN),
        seeds={"base": seed, "sampling": [seed + 1, seed + 2, seed + 3],
               "mc_starts": [seed + 100, seed + 300, seed + 600]},
    )


def save_calibration(report: CalibrationReport, path: str | Path = REPORT_PATH) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(report.as_dict(), ensure_ascii=False, indent=2), encoding="utf-8")
    return p


def load_calibration(path: str | Path = REPORT_PATH) -> dict:
    """读固化报告；不存在时给出明确错误，不静默生成。"""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(
            f"未找到校准报告 {p}；请运行 `interface-atlas calibrate --rebuild` 生成"
        )
    return json.loads(p.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 推理：给定组织光学参数 → 带覆盖保证的剂量区间
# ---------------------------------------------------------------------------

def calibrated_interval(
    tissue: TissueOptics,
    *,
    target_depth_cm: float = 0.5,
    alpha: float = DEFAULT_ALPHA,
    geometry: str = "plane",
    report: dict | None = None,
    max_depth_cm: float | None = None,
) -> dict:
    """P1 点估计 + (1-α) 覆盖的 μeff 与靶区剂量区间。

    剂量区间通过「光学参数同尺度缩放 → μeff 齐次缩放」精确得到，保持 a'、g 不变。
    """
    rep = report if report is not None else load_calibration()
    key = str(alpha)
    if key not in rep["q_by_alpha"]:
        raise ValueError(f"alpha={alpha} 不在报告中；可用：{list(rep['q_by_alpha'])}")
    b = float(rep["bias_b"])
    q = float(rep["q_by_alpha"][key])
    if not math.isfinite(q):
        raise RuntimeError("区间分位数为 inf（样本不足），无法给出有限区间")

    r_lo, r_hi = b * math.exp(-q), b * math.exp(q)
    mu_p1 = math.sqrt(3.0 * tissue.mu_a * (tissue.mu_a + tissue.mu_s_prime))

    # 用缩放后的组织重算剂量场，取靶区相对剂量（表面归一=1）
    def _target_fluence(scale: float) -> float:
        scaled = TissueOptics(tissue.mu_a * scale, tissue.mu_s_prime * scale,
                              tissue.g, tissue.n, tissue.name + f" (x{scale:.3f})")
        res = dose_field(scaled, target_depth_cm=target_depth_cm, geometry=geometry,
                         depths_cm=(1e-6, target_depth_cm), max_depth_cm=max_depth_cm or 2.0)
        return res.points[1].fluence_rel

    # 去偏点估计：用 scale=b 重算（最接近真实值，且落在 conformal 区间内）。
    # 原始 P1 值（scale=1）系统性偏高，单列供对比，不再冒充"点估计"。
    point = _target_fluence(b)
    p1_raw = _target_fluence(1.0)
    # r 大 ⇒ 衰减快 ⇒ 深部剂量小：剂量下界用 r_hi，上界用 r_lo
    dose_lo, dose_hi = _target_fluence(r_hi), _target_fluence(r_lo)
    return {
        "alpha": alpha,
        "coverage_target": round(1.0 - alpha, 4),
        "coverage_measured_test": rep["coverage_by_alpha"].get(key),
        "bias_b": b,
        "mu_eff_p1_per_cm": round(mu_p1, 6),
        "mu_eff_interval_per_cm": [round(r_lo * mu_p1, 6), round(r_hi * mu_p1, 6)],
        "r_interval": [round(r_lo, 6), round(r_hi, 6)],
        "target_depth_cm": target_depth_cm,
        # 不 round：深靶区剂量可低至 1e-20 量级，round(8) 会抹成 0
        "target_fluence_point": point,        # 去偏后点估计
        "target_fluence_p1_raw": p1_raw,      # 原始 P1（偏乐观）
        "target_fluence_interval": [dose_lo, dose_hi],
        "geometry": geometry,
        "in_domain": in_domain(tissue),
        "coverage_statement": (
            f"若组织参数在扩散适用域内且与校准样本 exchangeable，"
            f"真实 μeff 落入所示区间的概率 >= {1-alpha:.0%}（split conformal 有限样本保证）。"
            "剂量区间由 μeff 区间经光学参数同尺度缩放精确传播，非直接对剂量做 conformal。"
        ),
        "disclaimer": "科研降阶模型，非临床处方；域外不适用覆盖保证，见 in_domain 与 WARN",
    }


def in_domain(tissue: TissueOptics) -> dict:
    """检查组织参数是否落在校准/扩散适用域内。"""
    a, s, g = tissue.mu_a, tissue.mu_s_prime, tissue.g
    la, ha = DOMAIN["mu_a_range"]
    ls, hs = DOMAIN["mu_s_prime_range"]
    lg, hg = DOMAIN["g_range"]
    checks = {
        "mu_a": la <= a <= ha,
        "mu_s_prime": ls <= s <= hs,
        "g": lg <= g <= hg,
        "weak_absorption": a < s / 5.0,
    }
    return {"ok": all(checks.values()), "checks": checks, "domain": DOMAIN}


def _fmt(x: float) -> str:
    """小剂量用科学计数法，避免 round 后显示成 0.00000 误导用户。"""
    if x == 0 or 1e-4 <= abs(x) < 1e6:
        return f"{x:.5f}"
    return f"{x:.3e}"


def render(interval: dict) -> str:
    lo, hi = interval["target_fluence_interval"]
    mlo, mhi = interval["mu_eff_interval_per_cm"]
    lines = [
        "## L2 校准剂量区间（split conformal）",
        "",
        f"- 置信水平：覆盖率目标 ≥ {interval['coverage_target']:.0%}"
        f"（独立 test 实测 {interval['coverage_measured_test']:.0%}）",
        f"- 全局偏差因子 b = {interval['bias_b']:.4f}（>1 ⇒ P1 高估深部剂量）",
        f"- μeff：P1 {interval['mu_eff_p1_per_cm']:.4f} /cm "
        f"→ 区间 [{mlo:.4f}, {mhi:.4f}] /cm",
        f"- 靶区 @{interval['target_depth_cm']:.3f} cm 相对剂量（表面=1）：",
        f"    - 去偏点估计 {_fmt(interval['target_fluence_point'])}，"
        f"覆盖区间 [{_fmt(lo)}, {_fmt(hi)}]",
        f"    - 原始 P1 {_fmt(interval['target_fluence_p1_raw'])}（偏高，仅供对比）",
        f"- 适用域：{'内' if interval['in_domain']['ok'] else '外（不声明覆盖保证）'}",
        "",
        f"> {interval['coverage_statement']}",
        f"> {interval['disclaimer']}",
    ]
    return "\n".join(lines)
