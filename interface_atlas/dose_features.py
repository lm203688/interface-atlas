"""L3 剂量-效应基准的供给层：三模态统一剂量场特征向量。

为什么不是单一标量
------------------
VCC2025 的教训是单一 MAE leaderboard 会把模型引向"预测均值、丢掉分布"。
剂量-效应领域更严重：Ghaempanah et al. 2025 实测证明，换一套 CSF 分割
管线后电场**峰值几乎不变**（0.236->0.224 V/m），但**聚焦性体积掉了 ~39%**
（10.9->6.43 cm3）。若只用峰值做特征，"剂量相同"的结论会掩盖靶区外
暴露的巨大差异。因此本层强制输出**多指标**：峰值 + 空间衰减 + 穿透
深度 + 聚焦性 + 分位数覆盖，并对每个指标标注可信度来源。

诚实边界
--------
* 三模态解析引擎（P1 扩散 / 轴上 RS / 同心球球谐）都**不返回完整 3D 场**，
  聚焦性与覆盖体积在解析解下只能是 proxy（用焦域尺寸或轴上分位数粗估），
  逐项标 ``analytic | proxy | unavailable``，绝不冒充数值场体积分。
* 要拿到可信的聚焦性/靶区外体积，须接数值场（MC / FEM / SimNIBS），
  本模块提供 :func:`from_numeric_field` 钩子，由用户接入真实场数据。
* 跨模态峰值单位不同（V/m、MPa、相对光通量），**不做虚假归一化排名**；
  跨模态比较只在 L3 基准里按"各模态自己的效应锚点"分别做。

MIT License.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Sequence

# ---------------------------------------------------------------------------

DISCLAIMER = (
    "DoseFeatureVector 是多指标剂量场描述子，不是疗效预测；"
    "标注为 proxy 的指标来自解析粗估或自动选轴剖面，不能等同严格数值场体积分；"
    "标注为 analytic 的体积指标来自真解体积分，其准确性仍受来源网格/边界条件约束；"
    "跨模态峰值单位不同，禁止直接横向排名。"
)

RELIABILITY_LEVELS = ("analytic", "proxy", "unavailable")
MODALITIES = ("photobiomodulation", "focused ultrasound", "transcranial electric")


@dataclass(frozen=True)
class DoseFeatureVector:
    """一次剂量场计算的多指标特征。

    modality            模态规范名
    subject_id          头模/样本 ID（来自 anatomy.HeadModel）
    anatomy_source      phantom | subject_imaging | population
    peak                峰值剂量（模态原生单位）
    peak_unit           单位字符串
    peak_reference      峰值参考标度（MI 限 / 实践上限 / 表面=1），仅解释用
    penetration_depth_cm  达到峰值 1/e 的深度（光/声）或靶区深度（电）
    decay_per_cm        轴向衰减率（ln 峰值的负梯度，1/cm）
    fwhm_cm             半高全宽（焦域轴向长度 / 光穿透半高深度）
    coverage_depths_cm  达到峰值各分位数的深度 {分位数: cm}
    focality_proxy      聚焦性代理（越小越聚焦；0=点，1=全靶区）
    focality_basis      聚焦性代理的推导来源说明
    reliability         逐指标来源：analytic | proxy | unavailable
    warnings            适用性与诚实告警
    provenance          溯源（引擎版本、体模 DOI、seed 等）
    """
    modality: str
    subject_id: str
    anatomy_source: str
    peak: float
    peak_unit: str
    peak_reference: str = ""
    penetration_depth_cm: float = float("nan")
    decay_per_cm: float = float("nan")
    fwhm_cm: float = float("nan")
    coverage_depths_cm: dict[float, float] = field(default_factory=dict)
    coverage_volumes_cm3: dict[float, float] = field(default_factory=dict)
    focality_proxy: float = float("nan")
    focality_basis: str = ""
    reliability: dict[str, str] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "modality": self.modality,
            "subject_id": self.subject_id,
            "anatomy_source": self.anatomy_source,
            "peak": round(self.peak, 6) if math.isfinite(self.peak) else None,
            "peak_unit": self.peak_unit,
            "peak_reference": self.peak_reference,
            "penetration_depth_cm": _r(self.penetration_depth_cm),
            "decay_per_cm": _r(self.decay_per_cm),
            "fwhm_cm": _r(self.fwhm_cm),
            "coverage_depths_cm": {str(k): _r(v) for k, v in self.coverage_depths_cm.items()},
            "coverage_volumes_cm3": {str(k): _r(v) for k, v in self.coverage_volumes_cm3.items()},
            "focality_proxy": _r(self.focality_proxy),
            "focality_basis": self.focality_basis,
            "reliability": dict(self.reliability),
            "warnings": list(self.warnings),
            "provenance": self.provenance,
            "disclaimer": DISCLAIMER,
        }


def _r(x: float) -> float | None:
    return round(x, 6) if math.isfinite(x) else None


def _crossing_depths(points: Sequence, frac: float, peak: float) -> float | None:
    """轴上曲线首次降到 peak*frac 的深度（线性插值）；无交点返回 None。"""
    if peak <= 0:
        return None
    target = peak * frac
    for a, b in zip(points, points[1:]):
        va = getattr(a, "fluence_rel", None)
        vb = getattr(b, "fluence_rel", None)
        if va is None or vb is None:
            continue
        if va >= target >= vb:
            if va == vb:
                return a.depth_cm
            t = (va - target) / (va - vb)
            return a.depth_cm + t * (b.depth_cm - a.depth_cm)
    return None


# ---------------------------------------------------------------------------
# 三模态抽取器（纯函数，不修改上游结果）
# ---------------------------------------------------------------------------

def from_light_field(result: Any, *, subject_id: str = "analytic",
                     anatomy_source: str = "phantom") -> DoseFeatureVector:
    """从 light_field.FieldResult 抽特征。

    光模态的 fluence_rel 是归一化纵深剖面（表面=1），峰值/衰减/穿透深度都
    有轴上解析依据；聚焦性在 P1 扩散下无可靠定义，标 unavailable。
    """
    pts = tuple(result.points)
    peak = max(p.fluence_rel for p in pts)
    cov = {}
    for frac in (0.9, 1.0 / math.e, 0.5, 0.1):
        d = _crossing_depths(pts, frac, peak)
        if d is not None:
            cov[frac] = d
    d50 = cov.get(0.5, float("nan"))
    penetr = cov.get(1.0 / math.e, float("nan"))  # 1/e 交点已进 cov
    penetr = float("nan") if penetr is None else penetr
    if math.isnan(penetr):
        penetr = result.target_depth_cm * 2.0  # 粗上界，标 proxy
    # 衰减率：用 0.9->0.1 的双对数斜率，比单点更稳健
    if 0.9 in cov and 0.1 in cov and cov[0.1] > cov[0.9] > 0:
        decay = math.log(0.9 / 0.1) / (cov[0.1] - cov[0.9])
    else:
        decay = result.mu_eff_per_cm  # 解析回退
    warns = list(result.warnings)
    return DoseFeatureVector(
        modality="photobiomodulation",
        subject_id=subject_id, anatomy_source=anatomy_source,
        peak=peak, peak_unit="relative fluence (surface=1)",
        peak_reference="surface fluence = 1",
        penetration_depth_cm=penetr, decay_per_cm=decay,
        fwhm_cm=d50 if not math.isnan(d50) else float("nan"),
        coverage_depths_cm=cov,
        focality_proxy=float("nan"),
        focality_basis="P1 扩散无横向场，聚焦性不可用（需 MC/数值场）",
        reliability={
            "peak": "analytic", "penetration_depth_cm": "proxy" if cov.get(1.0/math.e) is None else "analytic",
            "decay_per_cm": "analytic", "fwhm_cm": "proxy",
            "focality_proxy": "unavailable",
        },
        warnings=tuple(warns + ["聚焦性不可用：光场聚焦需数值场，解析代理不提供横向分布"]),
        provenance={"mu_eff_per_cm": result.mu_eff_per_cm, "n_points": len(pts)},
    )


def from_ultrasound(result: Any, *, subject_id: str = "analytic",
                    anatomy_source: str = "phantom",
                    transducer: Any = None) -> DoseFeatureVector:
    """从 ultrasound_field.DoseResult 抽特征。

    超声有关轴/焦域尺寸（lateral_focus_mm / axial_focus_mm）。聚焦性 proxy =
    轴向焦域长度 / 焦距 F（无量纲，越小越聚焦）；**需要传入 ``transducer``**
    才能拿到 F，因为 DoseResult 不携带换能器几何。缺 transducer 时该项如实
    标 unavailable，不用 λ·f# 冒充孔径/焦距（那是焦斑尺度，概念错误）。
    MI/温升是 FDA/ITRUSST 阈值，作为 peak_reference。
    """
    lat = result.lateral_focus_mm / 10.0
    ax = result.axial_focus_mm / 10.0
    fwhm = ax
    focal_length = getattr(transducer, "focal_length_cm", None) if transducer is not None else None
    has_focal = isinstance(focal_length, (int, float)) and focal_length > 0
    foc = ax / focal_length if has_focal else float("nan")
    path_cm = float(focal_length) if has_focal else float("nan")
    # 衰减率：总衰减 dB -> 幅度比取 ln。缺焦距时用波长×f# 量级作 fallback 并标 proxy。
    denom = path_cm if has_focal else 2.0 * result.wavelength_cm * max(result.f_number, 1.0)
    decay = result.attenuation_db_total * math.log(10) / 20.0 / max(denom, 1e-9)
    return DoseFeatureVector(
        modality="focused ultrasound",
        subject_id=subject_id, anatomy_source=anatomy_source,
        peak=result.peak_pressure_mpa, peak_unit="MPa (peak pressure)",
        peak_reference=f"MI limit {result.mi_limit}; MI={result.mi:.3f}",
        penetration_depth_cm=focal_length if has_focal else float("nan"),
        decay_per_cm=decay, fwhm_cm=fwhm,
        coverage_depths_cm={0.5: fwhm},
        focality_proxy=foc,
        focality_basis=("轴向焦域长度/焦距 F（无量纲，越小越聚焦；proxy）" if has_focal
                        else "缺 transducer 几何（焦距），聚焦性不可用"),
        reliability={
            "peak": "analytic", "penetration_depth_cm": "analytic" if has_focal else "unavailable",
            "decay_per_cm": "proxy", "fwhm_cm": "proxy",
            "focality_proxy": "proxy" if has_focal else "unavailable",
        },
        warnings=tuple(list(result.warnings) + [
            "焦域为衍射受限粗估，非 k-Wave/声场数值积分结果",
            "颅骨折射与相位像差未建模，焦点位置可能偏移（Kop 2026 印证模拟剂量不可靠）",
        ]),
        provenance={"focus_gain": result.focus_gain, "mi": result.mi,
                    "attenuation_db_total": result.attenuation_db_total,
                    "focal_length_cm": focal_length if has_focal else None},
    )


def from_electric(result: Any, *, subject_id: str = "analytic",
                  anatomy_source: str = "phantom",
                  head_radius_cm: float = 10.0) -> DoseFeatureVector:
    """从 electric_field.ElectricDoseResult 抽特征。

    同心球球谐轴上解给靶区场强与球心场，可算"内外比"聚焦 proxy；横向
    分布与皮层聚焦性无解析解，标 unavailable。电流实践上限 2mA 作参考。
    """
    ratio = result.center_field_v_m / result.field_v_m if result.field_v_m > 0 else float("nan")
    # 内外比越小 -> 越聚焦在表层；这里只作粗代理，不解释成聚焦性高低
    foc = ratio if math.isfinite(ratio) else float("nan")
    path = max(result.target_depth_cm, 1e-9)
    decay = math.log(max(result.field_v_m, 1e-12) / max(result.center_field_v_m, 1e-12)) \
        / max(head_radius_cm - path, 1e-9) if path < head_radius_cm else float("nan")
    return DoseFeatureVector(
        modality="transcranial electric",
        subject_id=subject_id, anatomy_source=anatomy_source,
        peak=result.field_v_m, peak_unit="V/m (target field)",
        peak_reference=f"practice current <=2 mA; I={result.current_ma:.3f} mA (advisory)",
        penetration_depth_cm=result.target_depth_cm,
        decay_per_cm=decay, fwhm_cm=float("nan"),
        coverage_depths_cm={},
        focality_proxy=foc,
        focality_basis="球心场/靶区场比（粗代理；非皮层聚焦性体积分）",
        reliability={
            "peak": "analytic", "penetration_depth_cm": "analytic",
            "decay_per_cm": "proxy", "fwhm_cm": "unavailable", "focality_proxy": "proxy",
        },
        warnings=tuple(list(result.warnings) + [
            "同心球无皮层横向场，聚焦性不可用（需 FEM/SimNIBS 数值场）",
            "电导率个体差异可使靶区场变异 ~2.1x（Berger 2025），点估计不代表群体",
        ]),
        provenance={"target_layer": result.target_layer,
                    "n_harmonics": result.n_harmonics, "current_ma": result.current_ma},
    )


def from_numeric_field(*, modality: str, subject_id: str,
                       peak: float, peak_unit: str,
                       field_grid: Any, voxel_size_cm: float,
                       threshold_frac: float = 0.5,
                       anatomy_source: str = "subject_imaging") -> DoseFeatureVector:
    """从完整 3D 数值场抽取真正可信的聚焦性/覆盖体积。

    ``field_grid`` 为 ndarray（3D 剂量幅值），``voxel_size_cm`` 为各向同性
    体素边长。此处才计算体积分——由用户接入 MC/FEM/k-Wave/SimNIBS 场后调用。
    需要 numpy；未安装抛 ImportError。
    """
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError("from_numeric_field 需要 numpy") from exc
    grid = np.asarray(field_grid, dtype=float)
    if grid.ndim != 3 or voxel_size_cm <= 0:
        raise ValueError("field_grid 须为 3D，voxel_size_cm 须为正")
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        raise ValueError("场数据全为非有限值")
    p = float(finite.max())
    vox = voxel_size_cm ** 3
    n_half = int((finite >= threshold_frac * p).sum())
    cov_vol = n_half * vox
    total_vol = finite.size * vox
    foc = cov_vol / total_vol if total_vol > 0 else float("nan")
    return DoseFeatureVector(
        modality=modality, subject_id=subject_id, anatomy_source=anatomy_source,
        peak=p, peak_unit=peak_unit,
        penetration_depth_cm=float("nan"), decay_per_cm=float("nan"),
        fwhm_cm=float("nan"),
        coverage_depths_cm={threshold_frac: cov_vol},  # 此处语义为体积 cm3
        focality_proxy=foc,
        focality_basis=f"数值场 |field|>={threshold_frac}·peak 体积/总体积（真实体积分）",
        reliability={"peak": "analytic", "focality_proxy": "analytic",
                     "coverage_depths_cm": "analytic",
                     "penetration_depth_cm": "unavailable", "decay_per_cm": "unavailable"},
        warnings=("coverage_depths_cm 在本抽取器里语义为体积 cm3（非深度），"
                  "消费方须按 provenance.field_source 区分",),
        provenance={"field_source": "numeric", "voxel_size_cm": voxel_size_cm,
                    "coverage_volume_cm3": cov_vol, "total_volume_cm3": total_vol},
    )


def from_numeric_field_full(field: Any, *, modality: str, subject_id: str,
                            axis: int | None = None,
                            volume_fracs: Sequence[float] = (0.9, 1.0 / math.e, 0.5, 0.1),
                            anatomy_source: str = "subject_imaging") -> DoseFeatureVector:
    """全指标数值场抽取器：从真解算峰值 + 多分位体积 + 沿轴穿透/衰减/FWHM.

    与 :func:`from_numeric_field` 的区别
    -----------------------------------
    旧钩子只算峰值与半高体积，穿透/衰减/FWHM 全 unavailable。本抽取器补齐全
    指标——但**诚实按轴**：3D 场没有天然"深度轴"，穿透深度/衰减/FWHM 只是在
    某条过峰值的 1D 剖面上定义：

    * ``axis`` 由调用方显式给定（0/1/2）→ 这三项标 ``analytic``（剖面上严格
      计算），并在 basis 里写明是哪根轴。
    * ``axis=None`` → 自动选半高体积跨度最长的轴，三项标 ``proxy``（自动选轴
      不是剂量学轴，仅供粗看），basis 里写明自动选择。
    * 峰值在轴边界或剖面无交点 → 该项 ``unavailable``，不编造。

    体积分（覆盖体积/聚焦性）不依赖轴，始终按各向异性体素 ``sx*sy*sz`` 严格
    计算，标 ``analytic``。多分位体积写入 ``coverage_volumes_cm3``（新字段），
    不再污染 `coverage_depths_cm` 的深度语义。

    ``field`` 为 :class:`field_import.NumericField`（鸭子类型：需有 grid /
    unit / voxel_size_cm / source_tool / source_ref / warnings / provenance）。
    """
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError("from_numeric_field_full 需要 numpy") from exc

    grid = np.asarray(field.grid, dtype=float)
    if grid.ndim != 3:
        raise ValueError("NumericField.grid 须为 3D")
    vox = tuple(float(v) for v in field.voxel_size_cm)
    if len(vox) != 3 or any(v <= 0 for v in vox):
        raise ValueError("voxel_size_cm 须为三元正数组")
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        raise ValueError("场数据全为非有限值")

    peak = float(finite.max())
    flat_idx = int(np.argmax(grid))
    peak_idx = tuple(int(x) for x in np.unravel_index(flat_idx, grid.shape))
    vox_vol = vox[0] * vox[1] * vox[2]
    total_vol = float(grid.size) * vox_vol

    coverage_vol: dict[float, float] = {}
    for frac in volume_fracs:
        n = int((finite >= frac * peak).sum())
        coverage_vol[round(float(frac), 6)] = n * vox_vol

    half_frac = 0.5
    half_vol = coverage_vol.get(half_frac)
    if half_vol is None:
        half_vol = float((finite >= half_frac * peak).sum()) * vox_vol
    foc = half_vol / total_vol if total_vol > 0 else float("nan")
    axis_basis = "半高体积 cm3 / 总体积 cm3（各向异性体素体积分）"

    def axis_profile(ax: int) -> Any:
        """过峰值沿 ax 轴的 1D 剖面：固定另两轴索引于峰值。"""
        idx = list(peak_idx)
        idx[ax] = slice(None)
        return np.asarray(grid[tuple(idx)], dtype=float)

    # 选轴：显式给定 → analytic；自动选跨度最长轴 → proxy
    auto_axis = False
    if axis is None:
        spans = []
        for ax in range(3):
            prof_ax = axis_profile(ax)
            above = np.asarray(prof_ax >= half_frac * peak, dtype=bool)
            spans.append(int(above.sum()))
        axis = int(np.argmax(spans))
        auto_axis = True
    if axis not in (0, 1, 2):
        raise ValueError("axis 须为 0 / 1 / 2 或 None")
    profile_level = "proxy" if auto_axis else "analytic"
    cadence = vox[axis]

    # 过峰值的 1D 剖面
    prof = axis_profile(axis)
    pk_i = peak_idx[axis]

    def cross(side: int, frac: float) -> float | None:
        """从峰值往 side 方向找首次降到 peak*frac 的距离 cm（线性插值）。"""
        target = peak * frac
        i = pk_i
        while 0 <= i + side < prof.size:
            j = i + side
            a, b = prof[i], prof[j]
            if not (math.isfinite(a) and math.isfinite(b)):
                return None
            if a >= target >= b:
                if a == b:
                    return abs(j - pk_i) * cadence
                t = (a - target) / (a - b)
                cross_pos = i + side * t  # 沿 side 方向插值，负向不能写成 i+t
                return abs(cross_pos - pk_i) * cadence
            i = j
        return None

    d_lo = cross(-1, 1.0 / math.e)
    d_hi = cross(+1, 1.0 / math.e)
    penetr = max((d for d in (d_lo, d_hi) if d is not None), default=float("nan"))

    f_lo = cross(-1, half_frac)
    f_hi = cross(+1, half_frac)
    fwhm = (f_lo + f_hi) if (f_lo is not None and f_hi is not None) else float("nan")

    # 衰减率：用深侧 0.9→1/e 两点对数斜率（比单点稳健）
    deep_side = +1 if (d_hi or 0.0) >= (d_lo or 0.0) else -1
    c09 = cross(deep_side, 0.9)
    c1e = cross(deep_side, 1.0 / math.e)
    if c09 is not None and c1e is not None and c1e > c09 > 0:
        decay = math.log(0.9 * math.e) / (c1e - c09)
    else:
        decay = float("nan")

    axis_name = "xyz"[axis]
    axis_basis_profile = (
        f"过峰值沿 axis={axis}({axis_name}) 的 1D 剖面"
        + ("（自动选轴，非剂量学轴）" if auto_axis else "（调用方显式声明）")
    )
    reliability: dict[str, str] = {
        "peak": "analytic",
        "coverage_volumes_cm3": "analytic",
        "focality_proxy": "analytic",
        "penetration_depth_cm": profile_level if math.isfinite(penetr) else "unavailable",
        "decay_per_cm": profile_level if math.isfinite(decay) else "unavailable",
        "fwhm_cm": profile_level if math.isfinite(fwhm) else "unavailable",
    }
    warns = list(field.warnings) + [
        "穿透深度/衰减/FWHM 为过峰值 1D 剖面指标，非 3D 体积量；"
        + ("轴为自动选择（proxy），请显式传 axis 以冻结口径" if auto_axis
           else f"轴={axis_name} 由调用方声明"),
    ]

    prov = dict(field.provenance)
    prov.update({
        "field_source": field.source_tool,
        "source_ref": field.source_ref,
        "field_key": field.field_key,
        "field_unit": field.unit,
        "scale_applied": field.scale_applied,
        "voxel_size_cm": vox,
        "voxel_volume_cm3": vox_vol,
        "total_volume_cm3": total_vol,
        "half_max_volume_cm3": half_vol,
        "peak_index": peak_idx,
        "profile_axis": axis,
        "profile_axis_auto": auto_axis,
        "anatomy_source": anatomy_source,
    })
    return DoseFeatureVector(
        modality=modality, subject_id=subject_id, anatomy_source=anatomy_source,
        peak=peak, peak_unit=field.unit,
        peak_reference="外部数值场真解（网格/边界条件依赖）",
        penetration_depth_cm=penetr, decay_per_cm=decay, fwhm_cm=fwhm,
        coverage_depths_cm={},  # 真解无"深度分位"，体积走 coverage_volumes_cm3
        coverage_volumes_cm3=coverage_vol,
        focality_proxy=foc, focality_basis=axis_basis,
        reliability=reliability,
        warnings=tuple(warns),
        provenance=prov,
    )


# ---------------------------------------------------------------------------
# 群体特征分布（为 L3 供多主体输入，而非单一样本）
# ---------------------------------------------------------------------------

def cohort_summary(vectors: Iterable[DoseFeatureVector]) -> dict[str, Any]:
    """汇总一组特征向量的多指标分布（均值/标准差/分位数/可信度计数）。

    明确不输出单一"群体均值"作为结论——只给分布，消费方（L3 基准）
    必须逐指标使用，避免 VCC 单一 MAE 覆辙。
    """
    vs = tuple(vectors)
    if not vs:
        raise ValueError("vectors 不能为空")
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError("cohort_summary 需要 numpy") from exc

    def dist(getter: Callable[[DoseFeatureVector], float]) -> dict[str, Any]:
        vals = np.array([getter(v) for v in vs], dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return {"n": 0}
        return {
            "n": int(vals.size),
            "mean": round(float(vals.mean()), 6),
            "sd": round(float(vals.std(ddof=1)), 6) if vals.size > 1 else 0.0,
            "q05": round(float(np.percentile(vals, 5)), 6),
            "q50": round(float(np.percentile(vals, 50)), 6),
            "q95": round(float(np.percentile(vals, 95)), 6),
            "max_min_ratio": round(float(vals.max() / vals.min()), 4)
            if vals.min() > 0 else None,
        }

    rel_count: dict[str, int] = {}
    for v in vs:
        for k, lvl in v.reliability.items():
            key = f"{k}:{lvl}"
            rel_count[key] = rel_count.get(key, 0) + 1

    return {
        "n_subjects": len(vs),
        "modalities": sorted(set(v.modality for v in vs)),
        "peak": dist(lambda v: v.peak),
        "penetration_depth_cm": dist(lambda v: v.penetration_depth_cm),
        "decay_per_cm": dist(lambda v: v.decay_per_cm),
        "fwhm_cm": dist(lambda v: v.fwhm_cm),
        "focality_proxy": dist(lambda v: v.focality_proxy),
        "reliability_counts": rel_count,
        "single_metric_leaderboard_forbidden": True,
        "note": "多指标分布；禁止聚合成单一分数排名（VCC2025 教训）。",
    }


def render(v: DoseFeatureVector) -> str:
    lines = [
        f"剂量场特征：{v.modality}  subject={v.subject_id}  anatomy={v.anatomy_source}",
        f"  峰值 = {v.peak:.4g} {v.peak_unit}" + (f"  ({v.peak_reference})" if v.peak_reference else ""),
        f"  穿透深度 = {_r(v.penetration_depth_cm)} cm   轴向衰减 = {_r(v.decay_per_cm)} 1/cm",
        f"  半高宽 = {_r(v.fwhm_cm)} cm   聚焦性 proxy = {_r(v.focality_proxy)}",
        f"  聚焦性依据：{v.focality_basis or '—'}",
        "  分层覆盖（分位数→深度/体积）：",
    ]
    for frac, val in sorted(v.coverage_depths_cm.items()):
        lines.append(f"    ≥{frac:.0%} 深度 : {_r(val)} cm")
    for frac, val in sorted(v.coverage_volumes_cm3.items()):
        lines.append(f"    ≥{frac:.0%} 体积 : {_r(val)} cm3")
    lines.append("  可信度：" + ", ".join(f"{k}={lv}" for k, lv in v.reliability.items()))
    for w in v.warnings:
        lines.append(f"  ⚠ {w}")
    lines.append(f"  注：{DISCLAIMER}")
    return "\n".join(lines)
