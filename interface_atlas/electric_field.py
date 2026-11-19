"""L2 可计算剂量场引擎（电/tES-tDCS 模态 alpha）——从电极参数算到靶区电场.

为什么存在（承接 L1 契约）：
    L1 DoseAudit 把经颅电刺激报告普遍判为 not_reproducible——文献常只报
    「1 mA 阳极 M1、阴极对侧眶上」，不报电极面积、间距、个体颅骨电导率。
    本模块把「电流 + 电极几何 + 分层组织电导率 → 靶区电场 E(V/m)」这一跳
    显式化：与光/声模态平行，量化「设备报的电流 ≠ 靶区电场剂量」。

物理模型（四层同心球头模型，Rush & Driscoll 1969 经典解析基准，至今仍是
TES 正问题的标准 benchmark；准静态近似 ∇·(σ∇φ)=0）：
    每层电势球谐展开（轴对称） φ_j(r,θ)=Σ_l [A_j(u/u_j)^l + B_j(u_j/u)^{l+1}] P_l(cosθ)，
    u=r/R，u_j=层外半径/R。各阶 l 解耦，界面条件 φ 连续 + σ∂φ/∂r 连续，
    表面点电极对（阳极 θ=0 注入 +I，阴极 θ=π 流出 -I）给 Neumann 边界：
        c_l = I(2l+1)/(2πR²)  （奇数 l；偶数 l 因反对称消去）
    每层 2 个系数组成 2N×2N 线性系统，numpy.linalg.solve 逐阶求解。
    电场 E=-∇φ，轴上 θ=0 处 E_θ=0，E_r = -∂φ/∂r。

事实锚点（外部、已核实 2026-11-19，非编造）：
    * 四球电导率 S/m：头皮 0.332、颅骨 0.0083、CSF 1.79、脑 0.332；
      半径 cm：脑 8.0 / CSF 8.5 / 颅骨 9.2 / 头皮 10.0（PMC9927519；
      综述 skull 0.006-0.015、scalp 0.43、brain 0.33 同量级）。
    * 均匀球点电极对闭式：E(r)=I/(2πσR²)·[2(1+u²)/(1-u²)² + 1/(1-u²)]；
      代入 σ=0.33,R=9cm,I=1mA → 中心 E≈0.18 V/m，与文献 1mA 皮层
      ~0.2 V/m 量级一致。es_validation 用此闭式交叉核对球谐求解器（单层退化）。
    * 颅骨高阻抗：真实头模中仅约 25% 注入电流入脑（尸体量化研究，见 biomedeng
      综述）——该 25% 来自非对称头形/CSF 侧分流，本同心球轴对称模型无此效应，
      故不输出「入脑比例」（给了就是编造），只对场强负责。
    * 个体差异：Berger 团队 HCP 100 例群体头模（2025）显示组织电导率/解剖变异
      使皮层电场波动达 2.1 倍——「同电流不同靶区剂量」的定量证据。

安全口径（诚实，不编造监管硬线）：
    与超声 MI（FDA/AIUM 有明确定义）不同，tES **没有** FDA 硬阈值。
    1-2 mA 是二十年研究实践窗口，>2 mA 属研究性用法（非监管禁止）。
    本引擎把 2 mA 标为 PRACTICE_CURRENT_MAX_MA（实践上限，advisory），
    把电极下电流密度标为实践参考，均明确 NOTE 为 advisory，冒充硬阈值即撒谎。

明确的能力边界（不粉饰）：
    * 同心球、点电极：不含真实头形、脑沟回、电极-头皮接触阻抗、各向异性
      颅骨与白质。真实头模需 FEM（SimNIBS 等），本引擎是降阶可微代理。
    * 点电极在表面有奇异，靠近头皮/电极正下方场强不可信（会 WARN）；
      靶区应取脑内（默认电极下皮层 r≈7.8cm）。
    * 电场是即时场，不含 tDCS 后效应（NMDA 可塑性、膜电位时间常数）。
    * 绝对场强依赖组织电导率真值（个体差异达 2.1×），需个体化参数。

免责：科研降阶模型，非临床处方、非治疗计划系统。不构成医疗建议。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

DISCLAIMER = (
    "L2 电场剂量场引擎是四层同心球准静态降阶模型（Rush-Driscoll 解析基准），"
    "非临床处方或治疗计划系统；不含真实头形/脑沟回/各向异性/接触阻抗，"
    "点电极近表面场奇异；tES 无 FDA 硬阈值，2mA 为实践窗口（advisory）。"
)

# 实践窗口（非监管硬阈值，诚实标注 advisory）
PRACTICE_CURRENT_MAX_MA = 2.0        # tDCS 二十年研究实践上限；>2mA 为研究性用法
PRACTICE_ELECTRODE_J_MA_CM2 = 0.1     # 35cm² 电极 @2mA ≈ 0.057 mA/cm²；0.1 为实践参考
NOTE_ADVISORY = "advisory（研究实践参考，非 FDA/监管硬阈值）"


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TissueLayer:
    """同心球的一层（从内到外排列）。

    outer_radius_cm   该层外半径
    conductivity_s_m  电导率 S/m
    name              标注
    """
    outer_radius_cm: float
    conductivity_s_m: float
    name: str = "layer"

    def __post_init__(self):
        if self.outer_radius_cm <= 0:
            raise ValueError("外半径必须为正")
        if self.conductivity_s_m <= 0:
            raise ValueError("电导率必须为正")


@dataclass(frozen=True)
class ElectricDoseResult:
    """一次电场计算的完整结果。"""
    field_v_m: float                    # 靶区电场强度 |E| (V/m)
    field_mv_mm: float                  # 同上 (mV/mm)，神经元兴奋性常用单位
    potential_v: float                  # 靶区电势（相对无穷远/平均参考）
    current_ma: float                   # 注入电流
    target_depth_cm: float              # 靶区距头皮表面深度
    target_radius_cm: float             # 靶区距球心半径
    target_layer: str                   # 靶区所在组织层
    electrode_j_ma_cm2: float           # 电极下电流密度（点电极→用有限电极面积估计）
    n_harmonics: int                    # 球谐截断阶数
    center_field_v_m: float             # 球心场强（量级校验量）
    warnings: tuple[str, ...]

    def as_dict(self) -> dict:
        return {
            "field_v_m": round(self.field_v_m, 6),
            "field_mv_mm": round(self.field_mv_mm, 6),
            "potential_v": round(self.potential_v, 6),
            "current_ma": round(self.current_ma, 4),
            "target_depth_cm": round(self.target_depth_cm, 4),
            "target_radius_cm": round(self.target_radius_cm, 4),
            "target_layer": self.target_layer,
            "electrode_j_ma_cm2": round(self.electrode_j_ma_cm2, 5),
            "n_harmonics": self.n_harmonics,
            "center_field_v_m": round(self.center_field_v_m, 6),
            "warnings": list(self.warnings),
            "practice_current_max_ma": PRACTICE_CURRENT_MAX_MA,
            "practice_note": NOTE_ADVISORY,
            "disclaimer": DISCLAIMER,
        }


# 四层同心球头模型（事实锚点见模块 docstring）
HEAD_4SPHERE: tuple[TissueLayer, ...] = (
    TissueLayer(8.0, 0.332, "brain"),
    TissueLayer(8.5, 1.79, "csf"),
    TissueLayer(9.2, 0.0083, "skull"),
    TissueLayer(10.0, 0.332, "scalp"),
)
UNIFORM_SPHERE: tuple[TissueLayer, ...] = (TissueLayer(9.0, 0.33, "uniform"),)


# ---------------------------------------------------------------------------
# 球谐求解器（核心）
# ---------------------------------------------------------------------------

def _odd_l(n_harmonics: int) -> list[int]:
    """点电极反对称配置只需奇数阶：l=1,3,...,2n_harmonics-1。"""
    return [2 * k + 1 for k in range(n_harmonics)]


def solve_harmonics(
    layers: tuple[TissueLayer, ...],
    current_ma: float,
    *,
    n_harmonics: int = 32,
) -> dict[int, np.ndarray]:
    """逐阶求解球谐系数。

    返回 {l: coeff array}，每层两个系数 [A_j, B_j] 顺序排列，共 2N 个。
    基函数 φ_j(u)=A_j·(u/u_j)^l + B_j·(u_j/u)^{l+1}，u=r/R，u_j=层外半径/R。
    该基在界面处取值稳定（B_j 直接是界面电势贡献），避免 u^{-(l+1)} 爆炸。
    """
    N = len(layers)
    if N < 1:
        raise ValueError("至少一层")
    R = layers[-1].outer_radius_cm / 100.0        # 外半径 m
    u_j = np.array([lay.outer_radius_cm / layers[-1].outer_radius_cm for lay in layers])
    sigma = np.array([lay.conductivity_s_m for lay in layers])
    I = current_ma / 1000.0                       # mA → A
    coeffs: dict[int, np.ndarray] = {}

    for l in _odd_l(n_harmonics):
        M = np.zeros((2 * N, 2 * N))
        b = np.zeros(2 * N)

        def ai(j): return 2 * j
        def bi(j): return 2 * j + 1

        # 中心正则性：最内层无奇异项 B_0 = 0
        M[0, bi(0)] = 1.0
        row = 1

        # 层间界面 j=0..N-2（第 j 层外表面 = 第 j+1 层内表面，半径 u_j）
        for j in range(N - 1):
            r = u_j[j]
            ratio = r / u_j[j + 1]
            # 1) φ 连续：A_j + B_j = A_{j+1}·ratio^l + B_{j+1}/ratio^{l+1}
            M[row, ai(j)] = 1.0
            M[row, bi(j)] = 1.0
            M[row, ai(j + 1)] = -(ratio ** l)
            M[row, bi(j + 1)] = -(ratio ** (-(l + 1)))
            row += 1
            # 2) σ∂φ/∂r 连续（比例因子 1/(R·r) 两边相同，约去）
            #    第 j 层：σ_j[l·A_j - (l+1)·B_j]
            #    第 j+1 层：σ_{j+1}[l·A_{j+1}·ratio^{l-1} - (l+1)·B_{j+1}·ratio^{-(l+2)}]
            M[row, ai(j)] = sigma[j] * l
            M[row, bi(j)] = -sigma[j] * (l + 1)
            M[row, ai(j + 1)] = -sigma[j + 1] * l * (ratio ** (l - 1))
            M[row, bi(j + 1)] = sigma[j + 1] * (l + 1) * (ratio ** (-(l + 2)))
            row += 1

        # 外表面 u=1：点电流源 Neumann 边界
        # -σ[R dφ/dr] = c_l = I(2l+1)/(2πR²)  ⇒  σ[l·A_N - (l+1)·B_N] = -I(2l+1)/(2πR)
        M[row, ai(N - 1)] = sigma[N - 1] * l
        M[row, bi(N - 1)] = -sigma[N - 1] * (l + 1)
        b[row] = -I * (2 * l + 1) / (2 * math.pi * R)
        row += 1

        assert row == 2 * N
        sol = np.linalg.solve(M, b)
        sol[bi(0)] = 0.0   # 中心正则性是硬约束；求解器数值残差 ~1e-17 须显式归零，
                           # 否则 u→0 时 (u_j/u)^(l+1) 与该残差相乘仍会溢出
        coeffs[l] = sol

    return coeffs


def _layer_index(radius_cm: float, layers: tuple[TissueLayer, ...]) -> int:
    """半径（cm）所在层号（clamp 到最内层）。"""
    for j, lay in enumerate(layers):
        if radius_cm <= lay.outer_radius_cm:
            return j
    return len(layers) - 1


def field_and_potential_on_axis(
    layers: tuple[TissueLayer, ...],
    radius_cm: float,
    coeffs: dict[int, np.ndarray],
) -> tuple[float, float]:
    """轴上 (θ=0) 电场径向分量 E_r 与电势 φ。

    E_θ(θ=0)=0（轴对称），故 |E|=|E_r|。返回 (E_r [V/m], φ [V])。
    """
    R = layers[-1].outer_radius_cm / 100.0
    j = _layer_index(radius_cm, layers)
    u = (radius_cm / 100.0) / R
    u_j = layers[j].outer_radius_cm / layers[-1].outer_radius_cm

    e_r = 0.0
    phi = 0.0
    for l, c in coeffs.items():
        A, B = c[2 * j], c[2 * j + 1]
        # 最内层 B≡0（中心正则性）。u→0 时 (u_j/u)^(l+1) 会溢出，
        # 故 B 为零时直接跳过奇异项（物理上该项本就为 0，非数值近似）。
        phi += A * (u / u_j) ** l
        # dφ/dr = (1/R)[l·A·u^{l-1}/u_j^l - (l+1)·B·u_j^{l+1}/u^{l+2}]
        dedr = (1.0 / R) * (l * A * (u ** (l - 1)) / (u_j ** l))
        if abs(B) > 1e-300:
            phi += B * (u_j / u) ** (l + 1)
            dedr -= (1.0 / R) * (l + 1) * B * (u_j ** (l + 1)) / (u ** (l + 2))
        e_r -= dedr                      # E = -∇φ
    return e_r, phi


def uniform_sphere_closed_form(radius_cm: float, sigma: float, R_cm: float,
                               current_ma: float) -> float:
    """单层均匀球点电极对的轴上场强闭式（独立验证用，纯 math）。

    E(r) = I/(2πσR²)·[2(1+u²)/(1-u²)² + 1/(1-u²)]，u=r/R。
    推导：c_l=I(2l+1)/(2πR²)（奇数 l），E=Σ I(2l+1)/(2πσR²)·u^{l-1}，
    再用 Σ_{l odd}(2l+1)u^{l-1}=2(1+u²)/(1-u²)²+1/(1-u²) 求和闭式。
    验证：u→0 时 E=3I/(2πσR²)（与 l=1 主导一致）。
    """
    u = radius_cm / R_cm
    if u >= 1.0:
        raise ValueError("r 必须 < R")
    I = current_ma / 1000.0
    R = R_cm / 100.0
    s = 2.0 * (1.0 + u ** 2) / (1.0 - u ** 2) ** 2 + 1.0 / (1.0 - u ** 2)
    return I / (2.0 * math.pi * sigma * R ** 2) * s


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def dose_field(
    current_ma: float,
    *,
    layers: tuple[TissueLayer, ...] = HEAD_4SPHERE,
    target_depth_cm: float = 2.2,
    electrode_area_cm2: float = 35.0,
    n_harmonics: int = 32,
) -> ElectricDoseResult:
    """计算靶区电场 + 实践窗口核对。

    target_depth_cm   靶区距头皮外表面的深度；默认 2.2cm → r=7.8cm，落在脑层内
                      （头皮外半径 10cm，脑外半径 8cm）
    electrode_area_cm2 有限电极面积，仅用于估计电极下电流密度（点电极奇异）
    """
    if current_ma <= 0:
        raise ValueError("注入电流必须为正")
    R_cm = layers[-1].outer_radius_cm
    target_r = R_cm - target_depth_cm
    if target_r <= 0:
        raise ValueError("靶区深度超过头半径")
    if electrode_area_cm2 <= 0:
        raise ValueError("电极面积必须为正")

    coeffs = solve_harmonics(layers, current_ma, n_harmonics=n_harmonics)
    e_r, phi = field_and_potential_on_axis(layers, target_r, coeffs)
    center_e, _ = field_and_potential_on_axis(layers, 0.01, coeffs)  # 0.01cm 避开 u→0 数值奇异

    j_layer = _layer_index(target_r, layers)
    tname = layers[j_layer].name

    # 电极下电流密度（用有限电极面积估计，点电极本身奇异）
    electrode_j = current_ma / electrode_area_cm2

    # 注：不输出「入脑电流比例」。同心球轴对称模型下径向电流必穿过每层，
    # 无法定义有意义的「入脑比例」；文献 25% 来自非对称真实头模的 CSF 侧向
    # 分流，本降阶模型给该数字等于编造。诚实留空，见 L2_ELECTRIC.md。
    warns: list[str] = []
    if current_ma > PRACTICE_CURRENT_MAX_MA:
        warns.append(
            f"电流 {current_ma:.2f}mA > 实践窗口 {PRACTICE_CURRENT_MAX_MA}mA"
            f"（{NOTE_ADVISORY}）；>2mA 属研究性用法，需 IRB/安全评估"
        )
    if electrode_j > PRACTICE_ELECTRODE_J_MA_CM2:
        warns.append(
            f"电极下电流密度 {electrode_j:.3f}mA/cm² > 实践参考 "
            f"{PRACTICE_ELECTRODE_J_MA_CM2}mA/cm²（{NOTE_ADVISORY}）；"
            "过大有皮肤灼伤/电化学风险，宜增大电极面积"
        )
    if tname != "brain":
        warns.append(f"靶区落在 {tname} 层而非脑内；点电极近表面场奇异，结果不可信")
    if target_depth_cm < (R_cm - layers[-3].outer_radius_cm):
        warns.append("靶区距头皮过近，点电极近似在近场失效，需有限电极模型")
    warns.append(
        f"tES 无 FDA 硬阈值；实践窗口标注为 {NOTE_ADVISORY}，"
        "个体差异可使皮层场波动约 2.1×（HCP 群体头模）"
    )

    return ElectricDoseResult(
        field_v_m=abs(e_r),
        field_mv_mm=abs(e_r),                  # 1 V/m = 1 mV/mm（单位恒等）
        potential_v=phi,
        current_ma=current_ma,
        target_depth_cm=target_depth_cm,
        target_radius_cm=target_r,
        target_layer=tname,
        electrode_j_ma_cm2=electrode_j,
        n_harmonics=n_harmonics,
        center_field_v_m=abs(center_e),
        warnings=tuple(warns),
    )


# ---------------------------------------------------------------------------
# 逆设计与敏感性
# ---------------------------------------------------------------------------

def inverse_current(
    target_field_v_m: float,
    *,
    layers: tuple[TissueLayer, ...] = HEAD_4SPHERE,
    target_depth_cm: float = 2.2,
    electrode_area_cm2: float = 35.0,
    n_harmonics: int = 32,
) -> dict:
    """逆设计：给定靶区目标场强，反推所需注入电流（闭式，因 E 对 I 严格线性）。

    球谐系统对 I 线性 ⇒ E(I)=k·I。用参考电流 I_ref 求 k，再 I_req=E_target/k。
    正向重算复核残差。同时标出是否越过实践窗口（advisory，非硬阈值）。
    """
    if target_field_v_m <= 0:
        raise ValueError("目标场强必须为正")
    I_ref = 1.0
    ref = dose_field(I_ref, layers=layers, target_depth_cm=target_depth_cm,
                     electrode_area_cm2=electrode_area_cm2, n_harmonics=n_harmonics)
    k = ref.field_v_m / I_ref
    if k <= 0:
        raise ValueError("参考场强为 0，无法反推")
    i_req = target_field_v_m / k
    check = dose_field(i_req, layers=layers, target_depth_cm=target_depth_cm,
                       electrode_area_cm2=electrode_area_cm2, n_harmonics=n_harmonics)
    residual = abs(check.field_v_m - target_field_v_m) / target_field_v_m
    within = bool(i_req <= PRACTICE_CURRENT_MAX_MA)
    return {
        "within_practice_window": within,
        "required_current_ma": i_req,
        "target_field_v_m": target_field_v_m,
        "recomputed_field_v_m": check.field_v_m,
        "residual_rel": residual,
        "field_per_ma_v_m": k,
        "practice_current_max_ma": PRACTICE_CURRENT_MAX_MA,
        "note": (
            f"所需电流 {i_req:.3f}mA 超过实践窗口 {PRACTICE_CURRENT_MAX_MA}mA"
            f"（{NOTE_ADVISORY}）；属研究性用法，需安全评估"
        ) if not within else NOTE_ADVISORY,
        "caveat": "逆解基于线性同心球模型；绝对处方需个体头模/电导率实值复核",
    }


def skull_conductivity_sensitivity(
    *,
    skull_sigmas_s_m: tuple[float, ...] = (0.0042, 0.0083, 0.015, 0.03),
    current_ma: float = 1.0,
    target_depth_cm: float = 2.2,
    n_harmonics: int = 24,
) -> dict:
    """颅骨电导率敏感性：同电流/电极，颅骨电导率变 → 靶区场变（个体差异论点）。

    文献颅骨电导率跨度 0.002-0.03 S/m（Berger 群体头模 Beta 分布），
    本函数扫描该区间，输出靶区场强变化倍数。
    """
    cases = []
    for sk in skull_sigmas_s_m:
        layers = tuple(
            TissueLayer(lay.outer_radius_cm, sk if lay.name == "skull" else lay.conductivity_s_m, lay.name)
            for lay in HEAD_4SPHERE
        )
        r = dose_field(current_ma, layers=layers, target_depth_cm=target_depth_cm,
                       n_harmonics=n_harmonics)
        cases.append({
            "skull_conductivity_s_m": sk,
            "field_v_m": round(r.field_v_m, 6),
        })
    fields = [c["field_v_m"] for c in cases if c["field_v_m"] > 0]
    ratio = round(max(fields) / min(fields), 3) if len(fields) >= 2 else None
    return {
        "cases": cases,
        "field_max_min_ratio": ratio,
        "monotone_increasing": all(a["field_v_m"] <= b["field_v_m"]
                                   for a, b in zip(cases, cases[1:])),
        "interpretation": (
            "同一注入电流、仅颅骨电导率不同 → 靶区电场按颅骨阻抗变化；"
            "与 HCP 群体头模「组织电导率变异致皮层场波动 2.1×」同向，"
            "是「设备参数 ≠ 靶区剂量」的电模态定量证据"
        ),
        "literal_anchor": "Berger 等 HCP 100 例群体头模（2025）：组织电导率变异使皮层电场波动达 2.1×",
    }


# ---------------------------------------------------------------------------
# 渲染
# ---------------------------------------------------------------------------

def render(result: ElectricDoseResult) -> str:
    """人类可读输出（Markdown）。"""
    lines = [
        "## L2 电场剂量场（四层同心球，Rush-Driscoll 解析基准 alpha）",
        "",
        f"- 靶区（{result.target_layer}, 距头皮 {result.target_depth_cm:.1f}cm）"
        f"电场 **{result.field_v_m:.4f} V/m**（{result.field_mv_mm:.4f} mV/mm）",
        f"- 靶区电势 {result.potential_v:.4f} V；球心场 {result.center_field_v_m:.4f} V/m（量级校验）",
        f"- 注入电流 {result.current_ma:.2f} mA；电极下电流密度 {result.electrode_j_ma_cm2:.3f} mA/cm²",
        f"- 球谐截断 {result.n_harmonics} 阶；实践窗口 ≤{PRACTICE_CURRENT_MAX_MA}mA（{NOTE_ADVISORY}）",
    ]
    if result.warnings:
        lines += ["", "**适用性说明**："]
        lines += [f"- {w}" for w in result.warnings]
    lines += ["", f"> {DISCLAIMER}"]
    return "\n".join(lines)
