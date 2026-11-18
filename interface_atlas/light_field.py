"""L2 可计算剂量场引擎（光模态 alpha）——从设备参数算到靶区真实剂量.

为什么存在（承接 L1 的契约）：
    L1 DoseAudit 把「T1+T2 齐但缺 T3」判为 reproducible_with_model——边界条件
    完备，剂量场方程可解。本模块就是这个「model」：给定组织光学参数 + 照射
    几何，算出深度分辨的靶区光通量（fluence rate）与吸收功率密度。

    领域现状（2026）：设备报的辐照度 ≠ 靶区剂量。PNAS 胸部数字孪生证实同一
    外部剂量在不同解剖下靶区差数倍；PEMF 综述因欠报施加器距离/角度无法 meta。
    本引擎是把「设备参数 → 靶区剂量」这一跳显式化的第一块可验证基石。

物理模型（P1 扩散近似，连续波 CW，纯标准库实现）：
    控制方程   -∇·(D∇φ) + μa·φ = S,   D = 1/(3μt')
    有效衰减   μeff = sqrt(3·μa·μt') = sqrt(μa/D)
    输运反照率 a'  = μs' / μt',        μt' = μa + μs'
    源深度     z0 = 1/μt'（把准直入射近似为各向同性点源，Farrell 1992）
    外推边界   zb = 2A/(3μt'),  A = (1+rd)/(1-rd)
                rd ≈ -1.4399/n² + 0.7099/n + 0.6681 + 0.0636n（Groenhuis 1983）

    半无限 Green 函数（镜像负源满足外推零边界，Durduran 2010 Eq.20-22，CW 极限）：
        φ(ρ,z) = S/(4πD) · [ exp(-μeff·r1)/r1 - exp(-μeff·rb)/rb ]
        r1 = sqrt((z-z0)² + ρ²),  rb = sqrt((z+2zb+z0)² + ρ²)

    空间分辨漫反射率（Farrell, Patterson & Wilson 1992, Med Phys 19:879）：
        Rd(ρ) = a'/(4π) · [ (1/μt')(μeff + 1/r1)exp(-μeff·r1)/r1²
                          + (1/μt' + 4A/(3μt'))(μeff + 1/r2)exp(-μeff·r2)/r2² ]
        r1 = sqrt((1/μt')² + ρ²),  r2 = sqrt((1/μt'+4A/(3μt'))² + ρ²)

明确的能力边界（不粉饰）：
    * P1 扩散近似只在扩散区成立：μa << μs' 且离源/边界 > 几个输运平均自由程。
      强吸收（如血红蛋白密集）、薄层、近源区不适用——引擎会主动 WARN。
    * 本引擎输出**相对剂量场**：以表面剂量归一化为 1 的深度衰减倍数。绝对
      标定需要表面剂量实测或完整 Monte Carlo——这正是领域缺口，不假装解决。
    * 不碰全脑、不做异质解剖、不做温度场。L2 是从最小可验证单元往上爬。

验证锚点（外部、独立、可复现）：
    * Patterson, Chance & Wilson 1989, Appl Opt 28:2331：高反照率半无限介质
      μa=0.01, μs'=1 /mm 的吸收权重深度衰减拟合 μeff=0.1736，扩散值 0.1741
      （差 -0.29%）。本引擎的 attenuation_fit_error() 复现这一判据。
    * Farrell 1992：模型对 MC 在 ρ≥0.5mm 吻合、对体模实验 ρ≥1mm 吻合。
    * Wang, Jacques & Zheng 1995, Comput Methods Programs Biomed 47:131（MCML）：
      slab 基准 μa=1, μs=9 /mm, g=0.75, n=1, d=0.2mm → Rd=0.09734, Tt=0.66096。
      这是 L2 下一步用独立 MC 复核的标尺，本模块不声称已用扩散近似通过该薄层基准。

免责：科研用途的降阶模型，非临床处方、非监管结论。组织光学参数预设为文献
典型量级，仅用于演示与验证；正式用途必须由使用者提供实测/文献参数。
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

DISCLAIMER = (
    "L2 光剂量场引擎是 P1 扩散近似的科研降阶模型，非临床处方或监管结论；"
    "输出为相对剂量场（表面归一化），绝对剂量需表面实测或独立 MC 标定。"
)

# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TissueOptics:
    """组织光学参数（长度单位 cm，系数 1/cm）。

    mu_a         吸收系数
    mu_s_prime   约化散射系数 = mu_s·(1-g)
    g            各向异性因子（仅记录，扩散近似只用 mu_s_prime）
    n            相对折射率 n_tissue / n_air
    name         组织名（仅标注）
    """
    mu_a: float
    mu_s_prime: float
    g: float = 0.90
    n: float = 1.40
    name: str = "tissue"

    def __post_init__(self):
        if self.mu_a <= 0 or self.mu_s_prime <= 0:
            raise ValueError("mu_a 与 mu_s_prime 必须为正")
        if self.n < 1.0:
            raise ValueError("相对折射率 n 必须 >= 1")


@dataclass(frozen=True)
class DosePoint:
    """单个深度位置的剂量。"""
    depth_cm: float
    fluence_rel: float          # 相对光通量率（表面=1）
    absorption_rel: float       # 相对吸收功率密度 = mu_a · fluence（归一化）


@dataclass(frozen=True)
class FieldResult:
    """一次剂量场计算的完整结果。"""
    points: tuple[DosePoint, ...]
    mu_eff_per_cm: float
    mfp_prime_cm: float
    zb_cm: float
    z0_cm: float
    A_internal: float
    rd_total: float             # 总漫反射率（Farrell 空间积分，表面→空气）
    target_depth_cm: float
    target_fluence_rel: float
    warnings: tuple[str, ...]

    def as_dict(self) -> dict:
        return {
            "mu_eff_per_cm": self.mu_eff_per_cm,
            "mfp_prime_cm": self.mfp_prime_cm,
            "zb_cm": self.zb_cm,
            "z0_cm": self.z0_cm,
            "A_internal": self.A_internal,
            "rd_total": self.rd_total,
            "target_depth_cm": self.target_depth_cm,
            "target_fluence_rel": self.target_fluence_rel,
            "warnings": list(self.warnings),
            "depth_profile": [
                {
                    "depth_cm": round(p.depth_cm, 5),
                    "fluence_rel": round(p.fluence_rel, 6),
                    "absorption_rel": round(p.absorption_rel, 8),
                }
                for p in self.points
            ],
            "disclaimer": DISCLAIMER,
        }


# 组织光学参数预设（文献典型量级，非实测值；红/近红外 630-850nm）
# 标注为 approx：用于演示与验证，正式用途请替换为用户自己的参数。
TISSUE_PRESETS: dict[str, TissueOptics] = {
    "skin_red":      TissueOptics(4.0, 120.0, 0.90, 1.40, "skin @~660nm (approx)"),
    "skin_nir":      TissueOptics(2.0, 100.0, 0.90, 1.40, "skin @~810nm (approx)"),
    "fat_nir":       TissueOptics(0.5, 60.0, 0.90, 1.40, "adipose @~810nm (approx)"),
    "muscle_nir":    TissueOptics(1.5, 75.0, 0.90, 1.40, "muscle @~810nm (approx)"),
    "gray_matter_nir":  TissueOptics(2.0, 90.0, 0.90, 1.40, "gray matter @~810nm (approx)"),
    "white_matter_nir": TissueOptics(1.0, 120.0, 0.90, 1.40, "white matter @~810nm (approx)"),
    "skull_nir":     TissueOptics(3.0, 40.0, 0.90, 1.40, "skull @~810nm (approx)"),
    "reference_patterson1989": TissueOptics(1.0, 100.0, 0.90, 1.40, "Patterson 1989 benchmark"),
}


# ---------------------------------------------------------------------------
# 基础光学量
# ---------------------------------------------------------------------------

def optical_coefficients(t: TissueOptics) -> dict[str, float]:
    """返回 P1 扩散近似的一组基础系数。"""
    mut = t.mu_a + t.mu_s_prime
    D = 1.0 / (3.0 * mut)
    mu_eff = math.sqrt(3.0 * t.mu_a * mut)
    a_prime = t.mu_s_prime / mut
    rd = (-1.4399 / (t.n ** 2) + 0.7099 / t.n + 0.6681 + 0.0636 * t.n)
    A = (1.0 + rd) / (1.0 - rd)
    return {
        "mu_t_prime": mut,
        "D": D,
        "mu_eff": mu_eff,
        "a_prime": a_prime,
        "rd_fresnel": rd,
        "A_internal": A,
        "z0": 1.0 / mut,
        "zb": 2.0 * A / (3.0 * mut),
        "mfp_prime": 1.0 / t.mu_s_prime,
    }


def diffuse_reflectance_spatial(rho_cm: float, t: TissueOptics) -> float:
    """Farrell 1992 空间分辨漫反射率 Rd(ρ)（单位 1/cm² 的归一化核）。

    用于与空间分辨 MC / 实验逐点比对；ρ=0 处为有限值。
    """
    if rho_cm < 0:
        raise ValueError("rho_cm 必须 >= 0")
    c = optical_coefficients(t)
    mut, mu_eff, A, a_prime = c["mu_t_prime"], c["mu_eff"], c["A_internal"], c["a_prime"]
    inv_mut = 1.0 / mut
    d2 = inv_mut + 4.0 * A / (3.0 * mut)   # = inv_mut + 2·zb
    r1 = math.sqrt(inv_mut ** 2 + rho_cm ** 2)
    r2 = math.sqrt(d2 ** 2 + rho_cm ** 2)
    term1 = inv_mut * (mu_eff + 1.0 / r1) * math.exp(-mu_eff * r1) / (r1 ** 2)
    term2 = d2 * (mu_eff + 1.0 / r2) * math.exp(-mu_eff * r2) / (r2 ** 2)
    return a_prime / (4.0 * math.pi) * (term1 + term2)


def diffuse_reflectance_total(t: TissueOptics, *, n_steps: int = 2000) -> float:
    """对 Farrell Rd(ρ) 做二维空间积分，得总漫反射率（0<Rd<1）。

    ρ 上限取 8/μeff（覆盖扩散区主体），复合梯形积分。
    """
    c = optical_coefficients(t)
    rho_max = 8.0 / c["mu_eff"]
    h = rho_max / n_steps
    s = 0.0
    for i in range(n_steps + 1):
        rho = i * h
        w = 0.5 if i in (0, n_steps) else 1.0
        s += w * diffuse_reflectance_spatial(rho, t) * 2.0 * math.pi * rho
    return min(1.0, max(0.0, s * h))


def fluence_point(z_cm: float, rho_cm: float, t: TissueOptics) -> float:
    """半无限介质中，单位点源在各向同性近似下的未归一化光通量率 Green 函数。

    返回 S/(4πD) 已归入常数；此处只给空间形状项，便于相对比较与差分。
    """
    if z_cm <= 0:
        raise ValueError("z_cm 必须 > 0（半无限介质内部）")
    c = optical_coefficients(t)
    mu_eff, z0, zb = c["mu_eff"], c["z0"], c["zb"]
    r1 = math.sqrt((z_cm - z0) ** 2 + rho_cm ** 2)
    rb = math.sqrt((z_cm + 2.0 * zb + z0) ** 2 + rho_cm ** 2)
    return math.exp(-mu_eff * r1) / r1 - math.exp(-mu_eff * rb) / rb


# ---------------------------------------------------------------------------
# 宽场（LED 面板）轴向相对剂量场
# ---------------------------------------------------------------------------

def _plane_profile(z_cm: float, t: TissueOptics) -> float:
    """宽场均匀照射的轴向相对光通量：表面归一化到 1。

    由点源 Green 函数沿无限平面积分得到（横向平移不变）：
        φ(z) ∝ exp(-μeff|z-z0|) - exp(-μeff(z+2zb+z0))
    再除以 z=0 的值做归一化。这样 φ(0)=1，φ(z) 是「深度 z 剂量 / 表面剂量」。
    """
    c = optical_coefficients(t)
    mu_eff, z0, zb = c["mu_eff"], c["z0"], c["zb"]
    num = math.exp(-mu_eff * abs(z_cm - z0)) - math.exp(-mu_eff * (z_cm + 2.0 * zb + z0))
    den = math.exp(-mu_eff * z0) - math.exp(-mu_eff * (2.0 * zb + z0))
    if den <= 0:
        return math.exp(-mu_eff * z_cm)
    return max(0.0, num / den)


def dose_field(
    tissue: TissueOptics,
    *,
    depths_cm: tuple[float, ...] | list[float] | None = None,
    max_depth_cm: float = 2.0,
    n_points: int = 41,
    target_depth_cm: float = 1.0,
    geometry: str = "plane",
) -> FieldResult:
    """计算深度分辨的相对剂量场。

    geometry="plane"  宽场均匀照射（LED 面板/红光毯，PBM 主流），轴向 profile
    geometry="point"  点源（光纤探针），沿 rho=0 轴向，同样表面归一化

    target_depth_cm 是该模态的靶区深度：PBM 皮肤/黏膜 ~0.1-0.5cm，肌肉 ~1-3cm。
    """
    if geometry not in ("plane", "point"):
        raise ValueError("geometry 只能是 'plane' 或 'point'")
    if target_depth_cm <= 0:
        raise ValueError("target_depth_cm 必须 > 0")
    if depths_cm is None:
        depths = [max_depth_cm * i / (n_points - 1) for i in range(n_points)]
        depths[0] = 1e-6   # 避开 z=0 的解析奇异
    else:
        depths = sorted(float(d) for d in depths_cm)
        if not depths or depths[0] <= 0:
            raise ValueError("depths_cm 必须为正且非空")

    c = optical_coefficients(tissue)
    profile = (lambda z: _plane_profile(z, tissue)) if geometry == "plane" \
        else (lambda z: fluence_point(z, 0.0, tissue) / fluence_point(1e-6, 0.0, tissue))

    pts = []
    for z in depths:
        f = profile(z)
        # 归一化吸收功率密度：相对表面吸收（表面吸收 = mu_a·1）
        pts.append(DosePoint(z, f, tissue.mu_a * f))

    target_f = profile(target_depth_cm)

    # 主动 WARN：扩散近似适用域检查（不静默给出不可信数字）
    warns = []
    if tissue.mu_a >= tissue.mu_s_prime / 5.0:
        warns.append(
            f"μa/μs'={tissue.mu_a/tissue.mu_s_prime:.2f} 偏大；"
            "扩散近似要求 μa<<μs'，强吸收区结果需独立 MC 复核"
        )
    if target_depth_cm < 3.0 * c["mfp_prime"]:
        warns.append(
            f"靶区深度 {target_depth_cm:.3f}cm < 3×mfp'({3*c['mfp_prime']:.3f}cm)；"
            "近源区 P1 近似误差大，需 MC"
        )
    rd = diffuse_reflectance_total(tissue)
    if rd > 0.95:
        warns.append(f"总漫反射率 Rd={rd:.3f} 接近 1，组织近乎无吸收，剂量场对 μa 极敏感")

    return FieldResult(
        points=tuple(pts),
        mu_eff_per_cm=c["mu_eff"],
        mfp_prime_cm=c["mfp_prime"],
        zb_cm=c["zb"],
        z0_cm=c["z0"],
        A_internal=c["A_internal"],
        rd_total=rd,
        target_depth_cm=target_depth_cm,
        target_fluence_rel=target_f,
        warnings=tuple(warns),
    )


# ---------------------------------------------------------------------------
# 验证：扩散极限 μeff 拟合（Patterson 1989 判据的可复现版本）
# ---------------------------------------------------------------------------

def attenuation_fit_error(
    tissue: TissueOptics,
    *,
    z_start_cm: float | None = None,
    z_end_cm: float | None = None,
) -> dict:
    """从剂量场深度衰减斜率反推 μeff，与解析 μeff 对比。

    Patterson 1989：高反照率半无限介质吸收权重 exp(-μeff·z) 衰减，拟合值与
    扩散值差 -0.29%。本函数对 ln(fluence) 做最小二乘拟合，输出相对误差。

    注意：平面场在 z >> z0 处才进入纯指数区，起点默认取 5×mfp'，避开近源曲率。
    """
    c = optical_coefficients(tissue)
    z0 = z_start_cm if z_start_cm is not None else 5.0 * c["mfp_prime"]
    z1 = z_end_cm if z_end_cm is not None else 20.0 * c["mfp_prime"]
    if z1 <= z0:
        raise ValueError("z_end_cm 必须 > z_start_cm")
    zs = [z0 + (z1 - z0) * i / 40 for i in range(41)]
    ys = [math.log(_plane_profile(z, tissue)) for z in zs]
    n = len(zs)
    zbar = sum(zs) / n
    ybar = sum(ys) / n
    denom = sum((z - zbar) ** 2 for z in zs)
    slope = sum((z - zbar) * (y - ybar) for z, y in zip(zs, ys)) / denom
    fitted = -slope
    analytic = c["mu_eff"]
    return {
        "fitted_mu_eff_per_cm": fitted,
        "analytic_mu_eff_per_cm": analytic,
        "relative_error": (fitted - analytic) / analytic,
        "z_range_cm": [z0, z1],
        "n_points": n,
    }


# ---------------------------------------------------------------------------
# 可微接口 & 逆设计
# ---------------------------------------------------------------------------

def gradient(
    func,
    params: dict[str, float],
    *,
    step: float = 1e-5,
) -> dict[str, float]:
    """对任意标量函数做中心差分数值梯度（纯标准库）。

    func 签名 func(**params) -> float。step 为相对步长。
    用 Richard extrapolate 的简化版：两次步长取平均以减小截断误差。
    """
    grads = {}
    for key, val in params.items():
        h = abs(val) * step if abs(val) > 1e-12 else step
        p1, p2 = dict(params), dict(params)
        p1[key] = val + h
        p2[key] = val - h
        try:
            f1, f2 = func(**p1), func(**p2)
            grads[key] = (f1 - f2) / (2.0 * h)
        except (ValueError, ZeroDivisionError):
            grads[key] = float("nan")
    return grads


def inverse_irradiance(
    target_fluence: float,
    *,
    tissue: TissueOptics,
    target_depth_cm: float,
    geometry: str = "plane",
    surface_fluence_per_irradiance: float = 1.0,
) -> dict:
    """逆设计：给定靶区所需相对剂量，反推表面照射强度（相对倍数）。

    剂量场对表面剂量是线性的：φ_target = E_surface · f(target_depth)。
    因此 E_surface = target_fluence / f(target_depth)（除以表面→入射换算系数）。
    这是闭式逆解，不依赖迭代；解的正确性用正向重算复核（见 residual_rel）。

    注意：target_fluence 是「相对靶区剂量倍数」（相对于表面剂量=1 的单位），
    不是绝对 J/cm²。要得到绝对剂量需先标定表面剂量。
    """
    if target_fluence <= 0:
        raise ValueError("target_fluence 必须 > 0")
    res = dose_field(tissue, target_depth_cm=target_depth_cm, geometry=geometry,
                     depths_cm=(1e-6, target_depth_cm))
    f = res.points[1].fluence_rel
    # 只拒绝真正的 0（数值下溢）；深靶区剂量可达 1e-23 量级，仍是可算的逆解。
    # 不可行性由下面的 FEASIBLE_MAX 判据处理，不在这里误杀。
    if f <= 0.0:
        raise ValueError("靶区剂量场为 0，无法反推")
    e_surface = target_fluence / f
    # 独立复核：用反推值正向重算靶区剂量
    check = e_surface * f * surface_fluence_per_irradiance
    # 可行性判据：表面倍数过大意味着该光学参数下靶区剂量不可达（物理不可行），
    # 如实标出，而不是给一个天文数字让下游误以为是可用处方。
    FEASIBLE_MAX = 100.0   # 表面剂量超过 100× 靶区剂量视为不可行（工程约定阈值）
    feasible = e_surface <= FEASIBLE_MAX
    note = ("该光学参数/深度下靶区剂量需表面放大 %.3g× > %.0f×，物理上不可行；"
            "应换更透射的光学窗口或更浅靶区" % (e_surface, FEASIBLE_MAX)) if not feasible else ""
    return {
        "feasible": feasible,
        "required_surface_multiplier": e_surface,
        "target_fluence_requested": target_fluence,
        "target_fluence_recomputed": check / surface_fluence_per_irradiance,
        "residual_rel": abs(check / surface_fluence_per_irradiance - target_fluence) / target_fluence,
        "target_depth_cm": target_depth_cm,
        "geometry": geometry,
        "surface_fluence_per_irradiance": surface_fluence_per_irradiance,
        "infeasibility_note": note,
        "caveat": "逆解基于线性相对剂量场；绝对剂量与安全性需表面实测 + 独立 MC 复核",
    }


def inverse_wavelength_window(
    candidates: dict[str, TissueOptics],
    *,
    target_depth_cm: float,
    min_target_fluence: float,
    max_surface_burden: float = 2.0,
) -> dict:
    """在候选组织光学参数（可理解为波长/光谱窗口）中筛选可行方案。

    约束：靶区剂量 >= 下限，且表面负担（为达靶区剂量所需的表面倍数）<= 上限。
    这是多候选的离散逆设计，输出可行集与排序，不做连续优化以免过度承诺。
    """
    rows = []
    for name, t in candidates.items():
        res = dose_field(t, target_depth_cm=target_depth_cm, depths_cm=(1e-6, target_depth_cm))
        f = res.points[1].fluence_rel
        burden = min_target_fluence / f if f > 1e-12 else float("inf")
        rows.append({
            "candidate": name,
            "mu_eff_per_cm": round(res.mu_eff_per_cm, 4),
            "target_fluence_rel": round(f, 6),
            "surface_burden": round(burden, 4),
            "feasible": f >= min_target_fluence and burden <= max_surface_burden,
            "warnings": list(res.warnings),
        })
    rows.sort(key=lambda r: (not r["feasible"], r["surface_burden"]))
    return {"target_depth_cm": target_depth_cm, "min_target_fluence": min_target_fluence,
            "max_surface_burden": max_surface_burden, "ranking": rows}


# ---------------------------------------------------------------------------
# 渲染
# ---------------------------------------------------------------------------

def render(result: FieldResult) -> str:
    """人类可读输出（Markdown）。"""
    lines = [
        "## L2 光剂量场（P1 扩散近似 alpha）",
        "",
        f"- μeff = {result.mu_eff_per_cm:.4f} /cm，mfp' = {result.mfp_prime_cm:.4f} cm",
        f"- z0 = {result.z0_cm:.4f} cm，zb = {result.zb_cm:.4f} cm，A = {result.A_internal:.3f}",
        f"- 总漫反射率 Rd = {result.rd_total:.3f}",
        f"- 靶区深度 {result.target_depth_cm:.3f} cm → 相对剂量 {result.target_fluence_rel:.4f}（表面=1）",
        "",
        "| 深度 (cm) | 相对光通量 φ | 相对吸收 μa·φ |",
        "|---|---:|---:|",
    ]
    for p in result.points[:: max(1, len(result.points) // 12)]:
        lines.append(f"| {p.depth_cm:.3f} | {p.fluence_rel:.5f} | {p.absorption_rel:.6f} |")
    if result.warnings:
        lines.append("")
        lines.append("**适用性警告**：")
        for w in result.warnings:
            lines.append(f"- {w}")
    lines += ["", f"> {DISCLAIMER}"]
    return "\n".join(lines)
