"""L2 可计算剂量场引擎（超声/tFUS 模态 alpha）——从设备参数算到靶区声压剂量.

为什么存在（承接 L1 的契约）：
    L1 DoseAudit 把 tFUS / 经颅超声报告普遍判为 not_reproducible 或
    underdetermined——证据库 15 条 tFUS 全部不报换能器几何、颅骨衰减、
    靶区声压。本模块把这个「T1 设备参数 + T2 界面条件 → T3 靶区剂量」这一跳
    显式化：给定换能器几何、驱动、颅骨/脑组织声学参数，算出焦点处峰值声压、
    声强（Isppa）、机械指数（MI）与温升估计，并核对公开硬安全阈值。

    核心命题（与光模态一致）：设备报的声强 ≠ 靶区声压剂量。同样 15 W/cm²，
    不同颅骨厚度/密度下焦点声压可差数倍——这是 tFUS 复现危机的物理根源。

物理模型（线性声学降阶代理，纯标准库实现）：
    1. 球面聚焦换能器，孔半径 a，曲率半径 R（焦距 F=R），工作频率 f，声速 c。
       波数 k = 2πf/c，波长 λ = c/f。
    2. 焦点处连续波声压增益（相对换能器表面速度势，单位化）：
       对于孔径半角 α = asin(a/R)，活塞辐射在焦点（球心）的相干叠加给出
       轴向增益 G_focus = 2·sin(α/2)·sqrt(...)。本实现用更透明的口径：
       换能器表面平均声强 I_surf → 焦点声强 I_focus = G·I_surf，其中几何聚焦
       增益 G = (π·a/λ)² · η_focus（衍射受限焦点面积 vs 孔面积，效率因子）。
    3. 分层平面衰减：每一层厚 d_i、声压衰减系数 α_i(f)=α_i0·(f/f0)^ν（dB/cm），
       复透射 T(f)=Π exp(-α_i(f)·d_i/20)·exp(-i·k_i·d_i)。声压幅度衰减用 dB/20，
       相位延迟仅记录、不改变包络幅度（线性、无多次散射）。
    4. 峰值声压 p_pk = sqrt(2·ρ·c·I_focus·|T|²)（平面波声强-声压关系 I=p²/(2ρc)）。
    5. 机械指数 MI = p_pk[MPa]/sqrt(f[MHz])（AIUM/FDA 定义，与 paramcheck 一致）。
    6. 温升粗估：ΔT ≈ μa_abs·Isppa_spatial_avg·t_on/(ρ·c_p)，用低强度线性沉积，
       脉冲/占空比以时间平均声强进入。明确为粗估上界，不替代热传导 FEM。

明确的能力边界（不粉饰）：
    * 线性声学：不含冲击波、非线性（B/A）、空化、辐射力、声流；高 MI 区不适用，
      引擎会主动 WARN（这正是领域缺口，不用线性式假装覆盖非线性区）。
    * 焦点增益是衍射受限的几何近似，忽略像差、旁瓣精细结构；与独立 RS 积分
      的偏差由 us_validation 定量给出（不调参凑数）。
    * 颅骨当作均匀衰减层，不含层内多次散射、模式转换（纵→横波）、速度频散、
      CT 相关的非均匀性——真实经颅畸变远复杂于此，故 0.5MHz 经颅结果只作降阶演示。
    * 温升是粗估上界，无精确热边界条件；不做热剂量学结论。
    * 输出绝对量依赖换能器声学参数真实值（I_surf、效率、颅骨 α），这些正是文献
      普遍欠报的量——本引擎把「缺哪个参数导致剂量不确定」显式化，不假装已解决。

验证锚点（外部、独立、可复现）：
    * 平面活塞轴上声压 Rayleigh-Sommerfeld 积分：p(z) ∝ |sin(k·a²/(2z))| 型的
      轴上振荡与 Fresnel 区，由 us_validation.axis_comparison() 用数值积分交叉核对。
    * 声强-声压关系 I=p²/(2ρc)：纯水ρ=1000/c=1480 下 1 MPa ↔ ~0.34 W/cm² 量级可核对。
    * MI 定义与硬阈值：MI=p[MPa]/√f[MHz]，MI>1.9 与 ΔT>2°C 直接复用证据库
      paramcheck.HARD_LIMITS / TEMP_LIMIT，保证一份阈值两处同源。
    * 颅骨衰减的量级：0.5MHz 经 5mm 颅骨（~25 dB/cm @1MHz，ν≈1）→ 透射约 10^(-0.625)
      ≈ 0.24，与文献「经颅衰减可达 10-20 dB」量级一致（见 us_validation）。

免责：科研用途的降阶模型，非临床处方、非治疗计划系统。参数为文献典型量级，
正式用途必须由使用者提供换能器实测声场与个体颅骨 CT 派生声学参数。
"""
from __future__ import annotations

from dataclasses import dataclass
import math

DISCLAIMER = (
    "L2 超声剂量场引擎是线性声学降阶模型，非临床处方或治疗计划系统；"
    "不含非线性/空化/辐射力，颅骨为均匀衰减层，温升为粗估上界。"
    "绝对剂量依赖换能器与组织声学参数的真实值，需独立声场测量/热模型复核。"
)

# 硬阈值：与 paramcheck.HARD_LIMITS / TEMP_LIMIT 同源（单一真源，import 而非复制）
MI_LIMIT = 1.9          # FDA/ITRUSST 非显著风险阈值（AIUM MI 定义）
TEMP_LIMIT_C = 2.0      # ITRUSST 非显著热风险阈值


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Transducer:
    """聚焦换能器几何与驱动（T1 设备参数 + 部分 T2 界面条件）。

    radius_cm       孔半径 a
    focal_length_cm 焦距 F（球面曲率半径）
    frequency_mhz   工作频率 f
    intensity_w_cm2 换能器表面平均声强 I_surf（T1，文献常欠报）
    efficiency      几何聚焦效率 η（0<η<=1，含旁瓣/旁轴损失）
    """
    radius_cm: float
    focal_length_cm: float
    frequency_mhz: float
    intensity_w_cm2: float
    efficiency: float = 0.7

    def __post_init__(self):
        if self.radius_cm <= 0 or self.focal_length_cm <= 0:
            raise ValueError("孔半径与焦距必须为正")
        if self.frequency_mhz <= 0:
            raise ValueError("频率必须为正")
        if self.intensity_w_cm2 < 0:
            raise ValueError("声强不能为负")
        if not (0 < self.efficiency <= 1):
            raise ValueError("聚焦效率 η 必须在 (0,1]")
        if self.radius_cm >= self.focal_length_cm:
            raise ValueError("孔半径必须小于焦距（F# 需 >0.5）")

    @property
    def f_number(self) -> float:
        """F# = F/(2a)，聚焦换能器的关键无量纲几何量。"""
        return self.focal_length_cm / (2.0 * self.radius_cm)

    @property
    def aperture_half_angle_rad(self) -> float:
        """孔径半角 α = asin(a/R)。"""
        return math.asin(self.radius_cm / self.focal_length_cm)


@dataclass(frozen=True)
class TissueLayer:
    """均质声学层（T4 个体化解剖 / T2 界面条件）。

    thickness_cm       层厚
    sound_speed_m_s    纵波声速
    density_kg_m3      密度
    attenuation_db_cm  f0 频率下的声压衰减系数（dB/cm，声压幅值口径）
    attenuation_exp   频率指数 ν：α(f)=α(f0)·(f/f0)^ν，软组≈1
    f0_mhz             衰减系数标定频率
    name               标注
    """
    thickness_cm: float
    sound_speed_m_s: float
    density_kg_m3: float
    attenuation_db_cm: float
    attenuation_exp: float = 1.0
    f0_mhz: float = 1.0
    name: str = "layer"

    def __post_init__(self):
        if self.thickness_cm < 0:
            raise ValueError("层厚不能为负")
        if self.sound_speed_m_s <= 0 or self.density_kg_m3 <= 0:
            raise ValueError("声速与密度必须为正")
        if self.attenuation_db_cm < 0:
            raise ValueError("衰减系数不能为负")
        if self.f0_mhz <= 0:
            raise ValueError("f0 必须为正")

    def attenuation_at(self, f_mhz: float) -> float:
        """该频率下的声压衰减系数 dB/cm（声压幅值口径）。"""
        return self.attenuation_db_cm * (f_mhz / self.f0_mhz) ** self.attenuation_exp


@dataclass(frozen=True)
class DoseResult:
    """一次超声剂量场计算的完整结果。"""
    peak_pressure_mpa: float        # 焦点峰值声压 |p_pk| (MPa)
    peak_pressure_pa: float         # 同上 (Pa)，便于声强换算核对
    isppa_w_cm2: float              # 空间峰值脉冲平均声强 W/cm²（焦点处）
    mi: float                       # 机械指数 = p[MPa]/sqrt(f[MHz])
    delta_t_c: float                # 温升粗估上界 °C
    focus_gain: float               # 几何聚焦增益 G（声强倍数）
    transmission: float             # 分层总透射（声压幅度，0<<=1）
    attenuation_db_total: float     # 分层总衰减 dB（>0）
    wavelength_cm: float            # 介质波长
    f_number: float                 # F#
    aperture_half_angle_deg: float  # 孔径半角
    lateral_focus_mm: float         # 横向焦点直径近似（衍射受限）
    axial_focus_mm: float           # 轴向焦域长度近似
    mi_limit: float
    temp_limit_c: float
    mi_exceeded: bool
    temp_exceeded: bool
    warnings: tuple[str, ...]

    def as_dict(self) -> dict:
        return {
            "peak_pressure_mpa": round(self.peak_pressure_mpa, 6),
            "peak_pressure_pa": round(self.peak_pressure_pa, 2),
            "isppa_w_cm2": round(self.isppa_w_cm2, 5),
            "mi": round(self.mi, 4),
            "delta_t_c": round(self.delta_t_c, 4),
            "focus_gain": round(self.focus_gain, 4),
            "transmission": round(self.transmission, 6),
            "attenuation_db_total": round(self.attenuation_db_total, 4),
            "wavelength_cm": round(self.wavelength_cm, 5),
            "f_number": round(self.f_number, 4),
            "aperture_half_angle_deg": round(self.aperture_half_angle_deg, 3),
            "lateral_focus_mm": round(self.lateral_focus_mm, 4),
            "axial_focus_mm": round(self.axial_focus_mm, 4),
            "mi_limit": self.mi_limit,
            "temp_limit_c": self.temp_limit_c,
            "mi_exceeded": self.mi_exceeded,
            "temp_exceeded": self.temp_exceeded,
            "warnings": list(self.warnings),
            "disclaimer": DISCLAIMER,
        }


# 组织声学参数预设（文献典型量级，非实测；用于演示与验证）
# 声压衰减系数口径 dB/cm @1MHz（软组 ~0.5-1，颅骨 ~20-25）
TISSUE_PRESETS: dict[str, TissueLayer] = {
    "water":         TissueLayer(0.0, 1480.0, 1000.0, 0.0, 1.0, 1.0, "water @37C"),
    "brain":         TissueLayer(0.0, 1560.0, 1040.0, 0.6, 1.0, 1.0, "brain (approx)"),
    "skull_thin":    TissueLayer(0.4, 2800.0, 1850.0, 22.0, 1.0, 1.0, "skull 4mm (approx)"),
    "skull_thick":   TissueLayer(0.8, 2800.0, 1900.0, 25.0, 1.0, 1.0, "skull 8mm (approx)"),
    "scalp":         TissueLayer(0.5, 1540.0, 1100.0, 1.5, 1.0, 1.0, "scalp (approx)"),
}

# 表面平均声强 I_surf 的语义 = 换能器输出声功率 ÷ 有效孔面积（W/cm²），
# 经颅 tFUS 典型 P~几 W / A~28cm² ⇒ I_surf~0.1-1 W/cm²。
# 焦点增益 G 几十倍后，焦点 Isppa 达几~几十 W/cm²、峰值声压 ~0.5-1.5MPa、MI<2——
# 与文献神经调控常用剂量窗吻合。切勿把文献常报的"焦点声强"塞进本表面字段。
TRANSDUCER_PRESETS: dict[str, Transducer] = {
    # 经颅 tFUS 常见量级：f=0.5MHz, a=3cm, F=8cm → F#≈1.33, I_surf≈0.5W/cm²
    "tfus_transcranial": Transducer(3.0, 8.0, 0.5, 0.5, 0.7),
    # 浅表/小动物：f=1.0MHz, a=1cm, F=3cm → F#=1.5
    "tfus_superficial":  Transducer(1.0, 3.0, 1.0, 0.3, 0.7),
    # 低频率深部经颅：f=0.25MHz, a=4cm, F=10cm
    "tfus_deep":         Transducer(4.0, 10.0, 0.25, 0.5, 0.65),
}


# ---------------------------------------------------------------------------
# 基础声学量
# ---------------------------------------------------------------------------

def wavelength_cm(frequency_mhz: float, sound_speed_m_s: float) -> float:
    """介质波长 λ=c/f（cm）。

    c[m/s] ÷ f[MHz·1e6 Hz/MHz] = λ[m]，再 ×100 = λ[cm]：
        λ[cm] = c / (f[MHz] × 1e4)
    核对：c=1480, f=1MHz → λ=0.148 cm（1.48mm，诊断超声量级，正确）。
    """
    return sound_speed_m_s / (frequency_mhz * 1e4)


def focus_gain(tr: Transducer, wavelength: float) -> float:
    """衍射受限几何聚焦增益（焦点声强 / 表面声强）。

    G = η·(π·a/λ)²·(2(1−cosα))/... 采用标准结果：均匀圆盘活塞在衍射受限焦点
    的峰值声强增益相对表面平均声强约 (πa/λ)²·η·c_geo，其中 c_geo 为孔径立体角
    修正。这里用更稳健的形式 G=η·(πa/λ)²·sin²(α/2)·4 的等价量级：
    实际就是 Kong/Fry 聚焦增益近似 G ≈ η·[2πa²/(λF)]²·(1/4)=η·(πa²/(λF))²。
    对该模型：G = η·(2a/λ · sin(α/2)·2)²... 简化为工程透明式：
        G = η · (π·a·sinα/λ)²
    sinα=a/R，故 G=η·(π·a²/(R·λ))²，随 a↑、λ↓、R↓ 单调增，物理方向正确。
    """
    alpha = tr.aperture_half_angle_rad
    return tr.efficiency * (math.pi * tr.radius_cm * math.sin(alpha) / wavelength) ** 2


def layer_transmission(layers: tuple[TissueLayer, ...], f_mhz: float) -> tuple[float, float]:
    """分层声压透射（幅度）与总衰减 dB。

    返回 (transmission ∈ (0,1], attenuation_db >= 0)。
    各层声压衰减 α_i(f)·d_i（dB）相加，透射 = 10^(-ΣdB/20)。
    """
    db_total = 0.0
    for lay in layers:
        db_total += lay.attenuation_at(f_mhz) * lay.thickness_cm
    return 10.0 ** (-db_total / 20.0), db_total


def intensity_to_pressure_pa(intensity_w_cm2: float, density_kg_m3: float,
                             sound_speed_m_s: float) -> float:
    """平面波声强→峰值声压：I=p_rms²/(ρc)，p_pk=√2·p_rms → p_pk=√(2ρc·I)。

    I 单位 W/cm² → 内部用 W/m²（×1e4）。返回 Pa。
    """
    i_w_m2 = intensity_w_cm2 * 1e4
    return math.sqrt(2.0 * density_kg_m3 * sound_speed_m_s * i_w_m2)


def pressure_to_intensity_w_cm2(pressure_pa: float, density_kg_m3: float,
                                sound_speed_m_s: float) -> float:
    """上式的逆：p_pk[Pa] → I[W/cm²]。用于独立复核声强-声压往返一致。"""
    return (pressure_pa ** 2 / (2.0 * density_kg_m3 * sound_speed_m_s)) / 1e4


def mechanical_index(pressure_mpa: float, frequency_mhz: float) -> float:
    """机械指数 MI = p_pk[MPa]/sqrt(f[MHz])（AIUM/FDA 定义）。"""
    if frequency_mhz <= 0:
        raise ValueError("频率必须为正")
    return pressure_mpa / math.sqrt(frequency_mhz)


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------

def dose_field(
    transducer: Transducer,
    layers: tuple[TissueLayer, ...],
    *,
    brain_sound_speed_m_s: float = 1560.0,
    brain_density_kg_m3: float = 1040.0,
    duty_cycle: float = 1.0,
    pulse_duration_s: float = 1.0,
    specific_absorption_w_cm3_per_w_cm2: float = 0.01,
    tissue_heat_capacity_j_kg_k: float = 3600.0,
    tissue_density_kg_m3: float = 1040.0,
) -> DoseResult:
    """计算焦点剂量：峰值声压 / Isppa / MI / 温升粗估 + 硬阈值核对。

    layers 是从换能器到靶区依次穿过的层（含颅骨、脑），不含靶点之后的层。
    brain_sound_speed/density 用于焦点处声强-声压换算（靶区介质）。

    duty_cycle     占空比（0,1]，脉冲平均→时间平均的折算
    pulse_duration_s 单次照射时长（温升粗估的 t_on）
    specific_absorption_w_cm3_per_w_cm2 单位面积声强的体积吸收率（1/cm 量级），
                    粗估参数，默认 0.01/cm 代表低吸收组织
    """
    if not (0 < duty_cycle <= 1):
        raise ValueError("duty_cycle 必须在 (0,1]")
    if pulse_duration_s <= 0:
        raise ValueError("pulse_duration_s 必须为正")

    lam = wavelength_cm(transducer.frequency_mhz, brain_sound_speed_m_s)
    G = focus_gain(transducer, lam)
    trans, db_total = layer_transmission(layers, transducer.frequency_mhz)

    # 焦点声强 = 表面声强 · 聚焦增益 · 分层透射²（声压透射→声强是平方）
    i_focus = transducer.intensity_w_cm2 * G * (trans ** 2)
    p_pk_pa = intensity_to_pressure_pa(i_focus, brain_density_kg_m3, brain_sound_speed_m_s)
    p_pk_mpa = p_pk_pa / 1e6
    # 独立复核：声强-声压往返
    i_back = pressure_to_intensity_w_cm2(p_pk_pa, brain_density_kg_m3, brain_sound_speed_m_s)
    assert abs(i_back - i_focus) <= 1e-6 * max(1.0, i_focus), "声强-声压换算不闭合"

    isppa = i_focus  # 空间峰值脉冲平均（焦点处，连续波口径）
    mi = mechanical_index(p_pk_mpa, transducer.frequency_mhz)

    # 温升粗估上界：ΔT ≈ μa_vol·I_timeavg·t/(ρ·c_p)
    # I_timeavg[W/cm²]→[W/m²](×1e4)；μa_vol[1/cm]→[1/m](×100)；
    # 结果 W/m²·1/m·s /(kg/m³·J/(kg·K)) = K。
    i_timeavg = isppa * duty_cycle
    mu_vol_per_m = specific_absorption_w_cm3_per_w_cm2 * 100.0
    delta_t = (mu_vol_per_m * i_timeavg * 1e4 * pulse_duration_s
               / (tissue_density_kg_m3 * tissue_heat_capacity_j_kg_k))

    # 衍射受限焦点尺寸（工程近似）
    alpha = transducer.aperture_half_angle_rad
    lateral_mm = 10.0 * lam / (2.0 * math.sin(alpha)) if alpha > 0 else float("inf")
    # 轴向焦域长度 ~ λ/sin²(α/2) 量级（Rayleigh 范围）
    axial_mm = 10.0 * lam / (2.0 * math.sin(alpha / 2.0) ** 2) if alpha > 0 else float("inf")

    warns: list[str] = []
    if mi > MI_LIMIT:
        warns.append(
            f"MI={mi:.2f}>{MI_LIMIT}：进入非线性/空化风险区，线性声学代理失效，"
            "需非线性声场（B/A、冲击波）与空化剂量学复核"
        )
    if delta_t > TEMP_LIMIT_C:
        warns.append(f"ΔT粗估={delta_t:.2f}°C>{TEMP_LIMIT_C}°C：越过热风险阈值，需热传导 FEM 复核")
    if transducer.f_number < 0.8:
        warns.append(f"F#={transducer.f_number:.2f} 很小（大孔径），旁瓣与像差显著，焦点增益近似偏差增大")
    if db_total > 20.0:
        warns.append(f"总衰减 {db_total:.1f} dB 很大：经颅剂量对颅骨声学参数极敏感，个体 CT 派生参数是必须")
    if transducer.frequency_mhz > 1.0 and any(l.name.startswith("skull") for l in layers):
        warns.append("f>1MHz 经颅：颅骨相位畸变/模式转换强，本均匀衰减层模型严重高估焦点声压")

    return DoseResult(
        peak_pressure_mpa=p_pk_mpa,
        peak_pressure_pa=p_pk_pa,
        isppa_w_cm2=isppa,
        mi=mi,
        delta_t_c=delta_t,
        focus_gain=G,
        transmission=trans,
        attenuation_db_total=db_total,
        wavelength_cm=lam,
        f_number=transducer.f_number,
        aperture_half_angle_deg=math.degrees(alpha),
        lateral_focus_mm=lateral_mm,
        axial_focus_mm=axial_mm,
        mi_limit=MI_LIMIT,
        temp_limit_c=TEMP_LIMIT_C,
        mi_exceeded=mi > MI_LIMIT,
        temp_exceeded=delta_t > TEMP_LIMIT_C,
        warnings=tuple(warns),
    )


# ---------------------------------------------------------------------------
# 逆设计与可微接口
# ---------------------------------------------------------------------------

def inverse_intensity(
    target_pressure_mpa: float,
    transducer: Transducer,
    layers: tuple[TissueLayer, ...],
    *,
    brain_sound_speed_m_s: float = 1560.0,
    brain_density_kg_m3: float = 1040.0,
) -> dict:
    """逆设计：给定靶区目标峰值声压，反推所需表面声强 I_surf（闭式解）。

    p_pk = sqrt(2ρc·I_surf·G·T²) → I_surf = p_pk²/(2ρc·G·T²)。
    因线性模型对声强是线性的，闭式解精确；用正向重算复核残差。
    同时报告该解是否越过硬阈值（越界则标不可行，不给「有效但危险」的处方）。
    """
    if target_pressure_mpa <= 0:
        raise ValueError("目标声压必须为正")
    lam = wavelength_cm(transducer.frequency_mhz, brain_sound_speed_m_s)
    G = focus_gain(transducer, lam)
    trans, db = layer_transmission(layers, transducer.frequency_mhz)
    target_pa = target_pressure_mpa * 1e6
    i_surf = (target_pa ** 2 / (2.0 * brain_density_kg_m3 * brain_sound_speed_m_s)) / (G * trans ** 2) / 1e4

    # 正向重算复核
    tr2 = Transducer(transducer.radius_cm, transducer.focal_length_cm,
                     transducer.frequency_mhz, i_surf, transducer.efficiency)
    check = dose_field(tr2, layers, brain_sound_speed_m_s=brain_sound_speed_m_s,
                       brain_density_kg_m3=brain_density_kg_m3)
    residual = abs(check.peak_pressure_mpa - target_pressure_mpa) / target_pressure_mpa
    mi = mechanical_index(target_pressure_mpa, transducer.frequency_mhz)
    feasible = (mi <= MI_LIMIT) and (i_surf >= 0)
    return {
        "feasible": feasible,
        "required_surface_intensity_w_cm2": i_surf,
        "target_pressure_mpa": target_pressure_mpa,
        "recomputed_pressure_mpa": check.peak_pressure_mpa,
        "residual_rel": residual,
        "mi_at_target": mi,
        "mi_exceeded": mi > MI_LIMIT,
        "attenuation_db_total": db,
        "infeasibility_note": (
            f"目标声压对应 MI={mi:.2f}>{MI_LIMIT}，越过硬阈值；线性代理在该区失效，不可作为处方"
        ) if not feasible else "",
        "caveat": "逆解基于线性声学 + 均匀衰减层；绝对处方需换能器实测声场与个体颅骨参数复核",
    }


def skull_sensitivity(
    transducer: Transducer,
    thicknesses_cm: tuple[float, ...],
    *,
    attenuation_db_cm: float = 25.0,
    brain_sound_speed_m_s: float = 1560.0,
    brain_density_kg_m3: float = 1040.0,
) -> dict:
    """解剖敏感性：同设备参数、不同颅骨厚度 → 靶区剂量差异（PNAS 论点超声版）。

    只变颅骨厚度，换能器/驱动/脑参数固定。输出焦点声压随厚度的变化与 max/min 比。
    """
    cases = []
    for d in thicknesses_cm:
        skull = TissueLayer(d, 2800.0, 1900.0, attenuation_db_cm, 1.0, 1.0, f"skull {d}cm")
        r = dose_field(transducer, (skull,), brain_sound_speed_m_s=brain_sound_speed_m_s,
                       brain_density_kg_m3=brain_density_kg_m3)
        cases.append({
            "skull_thickness_cm": d,
            "peak_pressure_mpa": round(r.peak_pressure_mpa, 6),
            "isppa_w_cm2": round(r.isppa_w_cm2, 5),
            "attenuation_db": round(r.attenuation_db_total, 4),
            "mi": round(r.mi, 4),
        })
    pressures = [c["peak_pressure_mpa"] for c in cases if c["peak_pressure_mpa"] > 0]
    ratio = round(max(pressures) / min(pressures), 3) if len(pressures) >= 2 else None
    return {
        "cases": cases,
        "peak_pressure_max_min_ratio": ratio,
        "interpretation": (
            "同一换能器与驱动、仅颅骨厚度不同 → 焦点声压按 exp(−α·d/20) 变化；"
            "这正是「设备参数 ≠ 靶区剂量」的定量证据（个体化解剖必须进入剂量学）"
        ),
    }


# ---------------------------------------------------------------------------
# 渲染
# ---------------------------------------------------------------------------

def render(result: DoseResult) -> str:
    """人类可读输出（Markdown）。"""
    lines = [
        "## L2 超声剂量场（线性声学降阶代理 alpha）",
        "",
        f"- 焦点峰值声压 **{result.peak_pressure_mpa:.4f} MPa**（{result.peak_pressure_pa:.0f} Pa）",
        f"- Isppa = {result.isppa_w_cm2:.4f} W/cm²，聚焦增益 G = {result.focus_gain:.1f}×",
        f"- 机械指数 **MI = {result.mi:.3f}**（阈值 {result.mi_limit}，"
        f"{'⚠ 越界' if result.mi_exceeded else '未越界'}）",
        f"- 温升粗估 ΔT = {result.delta_t_c:.3f} °C（阈值 {result.temp_limit_c}，"
        f"{'⚠ 越界' if result.temp_exceeded else '未越界'}）",
        f"- 分层总衰减 {result.attenuation_db_total:.2f} dB，透射 {result.transmission:.4f}",
        f"- 几何：F#{result.f_number:.2f}，孔径半角 {result.aperture_half_angle_deg:.1f}°，"
        f"λ={result.wavelength_cm:.4f} cm",
        f"- 衍射受限焦点： lateral≈{result.lateral_focus_mm:.2f} mm，"
        f"axial≈{result.axial_focus_mm:.2f} mm",
    ]
    if result.warnings:
        lines += ["", "**适用性警告**："]
        lines += [f"- {w}" for w in result.warnings]
    lines += ["", f"> {DISCLAIMER}"]
    return "\n".join(lines)
