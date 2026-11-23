"""L2/L3 数值真解导入层：外部电场/声场 → 本包标准 DoseFeatureVector.

为什么存在
----------
P2 验收要求"电场/热场真解导入，SimNIBS / k-Wave 输出可读"。此前
``dose_features.from_numeric_field`` 只是一个 ndarray 钩子——调用方得自己把
SimNIBS 的 NIfTI、k-Wave 的 MAT/HDF5 读成数组、自己猜体素间距与单位。这正是
"外部工具不愿采用我们的口径"的原因：接入口太原始。

本模块把接入口标准化为一条明确契约：

    SimNIBS NIfTI / k-Wave MAT(HDF5) / 开放 NPZ
        → NumericField（显式单位 + 各向异性体素 + provenance）
        → dose_features.from_numeric_field_full（全指标）
        → DoseFeatureVector

诚实边界（硬性）
----------------
* **不猜单位**：场单位必须由调用方显式传入 ``field_unit``。SimNIBS E 场常为
  V/m、k-Wave 压力常为 Pa，但不同导出脚本会预缩放——猜单位就是伪造剂量。
* **不复制外部数据**：只读入文件路径，计算完即弃数组；不把受许可约束的真解
  写进本仓。provenance 记录来源路径与许可声明，由调用方先对头模做许可核验。
* **可选依赖优雅降级**：nibabel / h5py 未装时，NIfTI/HDF5 适配器抛带安装指引的
  ImportError，**判定内核与 CLI 其余命令不受影响**；scipy.io 可读 v5 MAT；
  开放 NPZ 契约零依赖。
* **不宣称真解等于体内剂量**：数值场是某套网格/边界条件/组织参数下的解，
  与活体剂量的差距不是本模块能消除的，告警如实保留。
* 各向异性体素按 (sx, sy, sz) 分别处理，体积分用 sx*sy*sz，不假设立方体素。

MIT License.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Any

DISCLAIMER = (
    "NumericField 是外部数值真解的只读包装；单位由调用方显式声明，本包不猜测。"
    "数值场是特定网格/边界条件/组织参数下的解，不等于体内实测剂量；"
    "不复制、不入库受许可约束的外部数据。"
)


@dataclass(frozen=True)
class NumericField:
    """标准数值场中间表示（只读、计算后即弃）。

    grid             3D numpy 数组，剂量幅值（标量场，如 |E| 或 |p|）
    unit             场单位字符串，由调用方显式声明（如 'V/m' / 'MPa'）
    voxel_size_cm    (sx, sy, sz) 各向异性体素边长 cm；grid 轴序与之对应
    source_tool      来源工具标识（simnibs / kwave / npz / array）
    source_ref       来源文件路径或标识（不复制内容）
    field_key        来源中的场变量名（MAT 变量 / NIfTI 说明 / NPZ key）
    scale_applied    读入时施加的线性缩放（如 mm→cm 不涉及；Pa→MPa=1e-6）
    provenance       溯源字典（许可声明、affine、shape 等）
    warnings         适用性与诚实告警
    """
    grid: Any
    unit: str
    voxel_size_cm: tuple[float, float, float]
    source_tool: str
    source_ref: str
    field_key: str = ""
    scale_applied: float = 1.0
    provenance: dict[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(self.grid.shape)

    @property
    def voxel_volume_cm3(self) -> float:
        sx, sy, sz = self.voxel_size_cm
        return sx * sy * sz

    def describe(self) -> str:
        lines = [
            f"数值场：source={self.source_tool}  key={self.field_key or '—'}",
            f"  shape={self.shape}  voxel(cm)={tuple(round(v, 6) for v in self.voxel_size_cm)}"
            f"  体素体积={self.voxel_volume_cm3:.6g} cm3",
            f"  单位={self.unit}（调用方声明）  缩放={self.scale_applied:g}",
            f"  来源引用={self.source_ref}（只读，未复制入库）",
        ]
        for k, v in self.provenance.items():
            lines.append(f"  溯源 {k}={v}")
        for w in self.warnings:
            lines.append(f"  ⚠ {w}")
        lines.append(f"  注：{DISCLAIMER}")
        return "\n".join(lines)


def _require_numpy() -> Any:
    try:
        import numpy as np  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError("field_import 需要 numpy（pip install numpy）") from exc
    return np


def _validate_grid(grid: Any) -> Any:
    np = _require_numpy()
    arr = np.asarray(grid)
    if arr.ndim != 3:
        raise ValueError(f"场数据须为 3D，实际 {arr.ndim}D shape={arr.shape}")
    if arr.size == 0:
        raise ValueError("场数据为空")
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise ValueError("场数据全为非有限值（NaN/inf）")
    return arr


def _norm_voxel(voxel_size_cm: Any) -> tuple[float, float, float]:
    """接受标量（各向同性）或三元组（各向异性），统一成三元组。"""
    if isinstance(voxel_size_cm, (int, float)):
        v = float(voxel_size_cm)
        if v <= 0:
            raise ValueError("voxel_size_cm 须为正")
        return (v, v, v)
    vals = tuple(float(x) for x in voxel_size_cm)
    if len(vals) != 3 or any(v <= 0 for v in vals):
        raise ValueError("voxel_size_cm 须为三元正数组 (sx, sy, sz)")
    return vals


def _check_unit(field_unit: str) -> str:
    if not isinstance(field_unit, str) or not field_unit.strip():
        raise ValueError("field_unit 必须显式声明（如 'V/m' / 'MPa'）；本包不猜单位")
    return field_unit.strip()


# ---------------------------------------------------------------------------
# 开放 NPZ 契约（零依赖，推荐外部工具优先导出这种格式）
# ---------------------------------------------------------------------------

def from_npz(path: str | Path, *, field_unit: str,
             field_key: str = "field") -> NumericField:
    """读取开放 NPZ 契约，零第三方依赖（numpy 已有）。

    NPZ 必须含：
      field           3D 幅值数组
      voxel_size_cm   标量或三元数组，单位 cm
    可选：
      unit            来源声明的单位（仍要求调用方用 field_unit 显式确认）
      source_tool     来源工具字符串
      provenance      JSON 字符串（许可、affine 等）
    """
    np = _require_numpy()
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"NPZ 不存在：{p}")
    unit = _check_unit(field_unit)
    with np.load(p, allow_pickle=False) as z:
        if field_key not in z:
            raise KeyError(f"NPZ 缺字段 {field_key!r}；现有：{list(z.files)}")
        grid = _validate_grid(z[field_key])
        if "voxel_size_cm" not in z:
            raise KeyError("NPZ 缺 voxel_size_cm（标量或三元，单位 cm）")
        vox = _norm_voxel(z["voxel_size_cm"])
        embedded_unit = str(z["unit"]) if "unit" in z else ""
        source_tool = str(z["source_tool"]) if "source_tool" in z else "npz"
        provenance: dict[str, Any] = {}
        if "provenance" in z:
            import json
            provenance = json.loads(str(z["provenance"]))

    warns: list[str] = []
    if embedded_unit and embedded_unit != unit:
        warns.append(f"NPZ 内声明单位 {embedded_unit!r} 与调用方 field_unit {unit!r} 不一致，以调用方为准")
    return NumericField(
        grid=grid, unit=unit, voxel_size_cm=vox,
        source_tool=source_tool, source_ref=str(p),
        field_key=field_key, scale_applied=1.0,
        provenance={**provenance, "npz_unit_declared": embedded_unit or None,
                    "shape": tuple(grid.shape)},
        warnings=tuple(warns),
    )


def save_npz(field: NumericField, path: str | Path, *,
             provenance_json: str = "") -> Path:
    """把 NumericField 写成开放 NPZ 契约（供外部工具消费/回灌）。

    注意：写出的是调用方已读入的数组——若来源受许可约束，调用方须自行确认
    有权转存。本函数不判断许可，只在 provenance 里保留许可声明字段。
    """
    np = _require_numpy()
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        p, field=field.grid,
        voxel_size_cm=np.asarray(field.voxel_size_cm, dtype=float),
        unit=field.unit, source_tool=field.source_tool,
        provenance=provenance_json or __import__("json").dumps(field.provenance, default=str),
    )
    return p


# ---------------------------------------------------------------------------
# SimNIBS NIfTI（nibabel 可选）
# ---------------------------------------------------------------------------

def from_simnibs_nifti(path: str | Path, *, field_unit: str,
                       license_declared: str = "",
                       subject_id: str = "") -> NumericField:
    """读 SimNIBS 导出的 NIfTI 标量场（如 |E| 幅值图，常见 V/m）。

    体素间距取 ``header.get_zooms()[:3]``，单位 mm → cm（÷10），**各向异性
    分别保留**。affine 原样写入 provenance，供消费方核对空间朝向——本包不
    自动重采样到某个"标准空间"，那会悄悄改变剂量场。

    需要 nibabel；未装时抛带指引的 ImportError，不影响其他命令。
    """
    p = Path(path)
    unit = _check_unit(field_unit)
    if not p.exists():
        raise FileNotFoundError(f"NIfTI 不存在：{p}")
    try:
        import nibabel as nib  # type: ignore
    except ImportError as exc:  # pragma: no cover - 依赖可选
        raise ImportError(
            "from_simnibs_nifti 需要 nibabel（pip install nibabel）；"
            "或先用外部脚本导出为开放 NPZ（field + voxel_size_cm）再用 from_npz 读取。"
        ) from exc

    img = nib.load(str(p))
    grid = _validate_grid(img.get_fdata())
    zooms_mm = tuple(float(z) for z in img.header.get_zooms()[:3])
    vox_cm = tuple(z / 10.0 for z in zooms_mm)
    vox_cm = _norm_voxel(vox_cm)  # 校验正数
    aff = getattr(img, "affine", None)

    warns = [
        "NIfTI 空间朝向未自动重采样；axis 顺序须由调用方核对 affine",
        "SimNIBS 组织电导率/边界条件不在本文件内，数值场误差无法由本包评估",
        "数值场 ≠ 体内实测剂量；个体差异（电导率 ~2.1x）不在此文件中消除",
    ]
    if not license_declared:
        warns.append("未声明来源数据许可；受 DUA/CC BY-NC 约束的数据须先做许可核验")

    return NumericField(
        grid=grid, unit=unit, voxel_size_cm=vox_cm,
        source_tool="simnibs", source_ref=str(p),
        field_key="nifti_scalar_field", scale_applied=1.0,
        provenance={
            "subject_id": subject_id or None,
            "license_declared": license_declared or None,
            "zooms_mm": zooms_mm,
            "affine": (aff.tolist() if aff is not None else None),
            "data_dtype": str(img.get_data_dtype()),
            "shape": tuple(grid.shape),
        },
        warnings=tuple(warns),
    )


# ---------------------------------------------------------------------------
# k-Wave MAT(v5 via scipy.io) / HDF5(v7.3 via h5py，可选)
# ---------------------------------------------------------------------------

def from_kwave_mat(path: str | Path, *, field_unit: str,
                   variable: str = "p", scale: float = 1.0,
                   voxel_size_cm: Any,
                   license_declared: str = "",
                   subject_id: str = "") -> NumericField:
    """读 k-Wave 输出场（MAT v5 用 scipy.io，v7.3 用 h5py 可选）。

    ``variable`` 为场数组变量名（k-Wave 常用 ``p`` / ``p_final`` / ``sensor_data``）。
    ``voxel_size_cm`` 必须由调用方给出：k-Wave 的 dx/dy/dz 通常在仿真脚本里，
    不在 .mat 文件元数据里——本包不猜网格间距。
    ``scale`` 为显式线性缩放（如源数组单位 Pa 要转 MPa，传 1e-6）；默认 1.0。

    需要 scipy（v5）或 h5py（v7.3）；两者缺一时给明确指引。
    """
    p = Path(path)
    unit = _check_unit(field_unit)
    vox_cm = _norm_voxel(voxel_size_cm)
    if not p.exists():
        raise FileNotFoundError(f"k-Wave 文件不存在：{p}")
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("scale 须为正有限数")

    raw: Any = None
    backend = ""
    # 先试 scipy.io（MAT v5）；失败再试 h5py（v7.3）
    try:
        from scipy.io import loadmat  # type: ignore
        mat = loadmat(str(p), variable_names=[variable], squeeze_me=False)
        if variable not in mat:
            raise KeyError(f"MAT 缺变量 {variable!r}")
        raw = mat[variable]
        backend = "scipy.io (MAT v5)"
    except ImportError:
        raw = None
    except Exception:
        raw = None  # 可能是 v7.3，落到 h5py

    if raw is None:
        try:
            import h5py  # type: ignore
        except ImportError as exc:  # pragma: no cover - 依赖可选
            raise ImportError(
                "读取 k-Wave MAT 需要 scipy（MAT v5）或 h5py（MAT v7.3）；"
                "或导出为开放 NPZ 后用 from_npz。pip install scipy h5py"
            ) from exc
        with h5py.File(str(p), "r") as fh:
            if variable not in fh:
                raise KeyError(f"HDF5 缺数据集 {variable!r}；现有：{list(fh.keys())[:10]}")
            raw = fh[variable][...]
        backend = "h5py (HDF5 v7.3)"

    grid = _validate_grid(raw)
    if scale != 1.0:
        grid = grid * scale
        grid = _validate_grid(grid)

    warns = [
        f"k-Wave 网格间距由调用方提供（文件不含 dx/dy/dz），已按 {vox_cm} cm 处理",
        "k-Wave 声场为仿真解；换能器几何/颅骨参数/热-声耦合不在本文件内，剂量误差无法由本包评估",
        "数值场 ≠ 体内实测声压；数值孔径与网格分辨率须由调用方核对",
    ]
    if not license_declared:
        warns.append("未声明来源数据许可；受 DUA/机构约束的数据须先做许可核验")

    return NumericField(
        grid=grid, unit=unit, voxel_size_cm=vox_cm,
        source_tool="kwave", source_ref=str(p),
        field_key=variable, scale_applied=scale,
        provenance={
            "subject_id": subject_id or None,
            "license_declared": license_declared or None,
            "read_backend": backend,
            "shape": tuple(grid.shape),
        },
        warnings=tuple(warns),
    )


def from_array(grid: Any, *, field_unit: str, voxel_size_cm: Any,
               source_tool: str = "array", source_ref: str = "in-memory",
               field_key: str = "", scale: float = 1.0) -> NumericField:
    """从内存 3D 数组构造 NumericField（测试/外部已解析好的场）。"""
    unit = _check_unit(field_unit)
    vox_cm = _norm_voxel(voxel_size_cm)
    arr = _validate_grid(grid)
    if scale != 1.0:
        arr = _validate_grid(arr * scale)
    return NumericField(
        grid=arr, unit=unit, voxel_size_cm=vox_cm,
        source_tool=source_tool, source_ref=source_ref,
        field_key=field_key, scale_applied=scale,
        provenance={"shape": tuple(arr.shape), "source": "in-memory array"},
        warnings=("内存数组来源：单位/网格由调用方声明，本包未核验来源工具输出语义",),
    )
