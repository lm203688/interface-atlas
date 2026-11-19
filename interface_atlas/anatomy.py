"""L2 个体化解剖层：公开头模/体模的数据契约 + 可复现代理体模 + 导入适配器。

为什么需要这一层
----------------
L1 剂量报告规范 v1（`dosimetry.py`）把 ``anatomy_source`` 定义为 tier-4
（subject，个体化解剖）字段，取值 ``phantom | subject_imaging | population``；
但 L2 三个剂量场引擎此前都跑在解析球/分层球上——设备参数到靶区剂量的
换算里，**解剖个体差异**这一维一直是空的。已有实测证据说明它不是小事：

* Berger et al. 2025 (Sci Data, DOI 10.1038/s41597-025-04886-0) 发布 100 个
  HCP 头模 + 变异组织电导率：同一 montage 下皮层场强个体差异约 2.1x。
* Ghaempanah et al. 2025 (BeNe Brain Stimulation) 只改 CSF 分割管线，
  电场峰值几乎不变（0.236 -> 0.224 V/m），但聚焦性体积掉了 ~39%。
  即：单一标量剂量无法描述个体差异，必须带解剖特征与场形指标。

本模块做什么 / 不做什么
----------------------
做：
* 定义头模数据契约 (:class:`HeadModel`, :class:`TissueLayerGeometry`)，
  带 DOI / license / commercial_ok / 坐标空间 / source 溯源字段。
* 内置三档**可复现代理体模**：``spherical``（解析同心球，与
  :data:`electric_field.HEAD_4SPHERE` 同构）、``mni152_style``（分层椭球代理）、
  ``hcp_population``（群体统计代理）。它们是用来跑通链路的 proxy，
  **不是**真实 MNI152 图像或真实 HCP 头模，绝不冒充。
* 提供真实数据导入适配器：NIfTI 分割掩码 -> 层厚/体积特征（nibabel 可选）；
  lead field 矩阵（HDF5，可选）。真实数据由用户自己提供。
* 许可门禁 :func:`check_license`：CC BY-NC / HCP DUA 数据不得进入商用管线，
  调用方必须显式声明 ``commercial=False``，否则 fail-closed。
* :func:`population_sample`：群体解剖抽样，给蒙特卡洛剂量区间供数。

不做（诚实边界）：
* 不下载、不打包 Berger/cuNIBS 的 HCP 头模（9.9 GB，CC BY 4.0 + HCP DUA）。
* 不打包 SimNIBS ``ernie_extended``（CC BY-NC 4.0）——许可与本项目 MIT 冲突。
* 群体分布默认参数是 **illustrative（示例）**，不是从真实队列提取的经验值；
  报告里如实标注，用户可用自己的队列统计替换。

MIT License. 数据许可见各体模 ``license`` 字段。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Any, Iterable

# ---------------------------------------------------------------------------
# 事实锚点常量（公开来源，勿凭记忆篡改）
# ---------------------------------------------------------------------------

SOURCES: dict[str, dict[str, str]] = {
    "analytic": {
        "doi": "N/A (解析构造)",
        "license": "MIT",
        "citation": "interface-atlas 解析代理体模",
        "note": "同心球分层取自 electric_field.HEAD_4SPHERE 的公开电导率值。",
    },
    "mni152_style": {
        "doi": "10.1016/j.neuroimage.2020.117044",
        "license": "MIT (本代理) / 原始分割法 CC BY 4.0",
        "citation": "Puonti et al. 2020 NeuroImage (CHARM 分割法)",
        "note": "本体模是分层椭球解析代理，不是真实 MNI152 图像；"
                "仅借用其组织分层与尺度量级。",
    },
    "hcp_population": {
        "doi": "10.1038/s41597-025-04886-0",
        "license": "CC BY 4.0 + HCP Open Access DUA",
        "citation": "Berger et al. 2025 Scientific Data (100 HCP head models)",
        "note": "本体模是群体统计代理，默认参数为 illustrative；"
                "真实 100 头模见 Zenodo 10.5281/zenodo.13259679，受 HCP DUA 约束。",
    },
    "siminibs_ernie": {
        "doi": "N/A",
        "license": "见 SimNIBS example dataset",
        "citation": "SimNIBS datasets (ernie, MNI152, spherical)",
        "note": "https://simnibs.github.io/simnibs/build/html/dataset.html",
    },
    "siminibs_ernie_extended": {
        "doi": "10.1162/imag_a_00379",
        "license": "CC BY-NC 4.0",
        "citation": "Van Hoornweder et al. 2024 Imaging Neuroscience",
        "note": "非商用许可；不得打包进 MIT 仓库或用于商用剂量管线。",
    },
    "cunibs_hcp": {
        "doi": "10.1038/s41597-025-04886-0",
        "license": "CC BY 4.0 + HCP Open Access DUA",
        "citation": "cuNIBS / CHARM 100 HCP meshes, Zenodo record 23511725",
        "note": "9.9 GB HDF5；本包只定义读取契约，不含数据。",
    },
}

# 四层组织默认物性。电导率来自四球头模型公开值（见 electric_field 模块）；
# 声衰减为量级参考，标注为 illustrative，需真实数据替换。
DEFAULT_CONDUCTIVITY_S_M = {
    "brain": 0.332, "csf": 1.79, "skull": 0.0083, "scalp": 0.332,
}
DEFAULT_ATTENUATION_DB_CM_MHZ = {
    "brain": 0.50, "csf": 0.05, "skull": 20.0, "scalp": 0.60,
}
ANATOMY_SOURCE_CHOICES = ("phantom", "subject_imaging", "population")
COMMERCIAL_ALLOWED_LICENSES = ("MIT", "CC BY 4.0", "Apache-2.0", "BSD-3-Clause")
ILLUSTRATIVE_TAG = "illustrative（示例分布，非真实队列经验值；用真实数据替换）"

DISCLAIMER = (
    "anatomy 层提供的是头模数据契约与可复现代理体模，不是临床患者模型；"
    "内置代理体模不能用于临床剂量决策。真实头模需用户自行接入并遵守其数据许可。"
)


# ---------------------------------------------------------------------------
# 数据契约
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TissueLayerGeometry:
    """一层组织的几何 + 物性。

    outer_radius_cm   球/椭球外半径（cm）；椭球代理下为长半轴。
    shape_ratio       椭球另两轴 / 长半轴（1.0 = 正球）。
    conductivity_s_m  电导率（S/m），喂 electric_field。
    attenuation_db_cm_mhz  声衰减（dB/cm/MHz），喂 ultrasound_field。
    name              组织名：brain / csf / skull / scalp。
    """
    name: str
    outer_radius_cm: float
    conductivity_s_m: float
    attenuation_db_cm_mhz: float
    shape_ratio: float = 1.0

    def __post_init__(self):
        if self.outer_radius_cm <= 0:
            raise ValueError("outer_radius_cm 必须为正")
        if not 0 < self.shape_ratio <= 1:
            raise ValueError("shape_ratio 必须落在 (0, 1]")
        if self.conductivity_s_m <= 0:
            raise ValueError("conductivity_s_m 必须为正")
        if self.attenuation_db_cm_mhz < 0:
            raise ValueError("attenuation_db_cm_mhz 不能为负")

    @property
    def thickness_cm(self) -> float:
        """由外层半径差得到本层厚度由 HeadModel 计算；单层时返回半径。"""
        return self.outer_radius_cm  # 由 HeadModel.layer_thicknesses() 覆盖

    def as_dict(self) -> dict:
        return {
            "name": self.name,
            "outer_radius_cm": round(self.outer_radius_cm, 4),
            "shape_ratio": round(self.shape_ratio, 4),
            "conductivity_s_m": round(self.conductivity_s_m, 6),
            "attenuation_db_cm_mhz": round(self.attenuation_db_cm_mhz, 4),
        }


@dataclass(frozen=True)
class HeadModel:
    """一个头模的完整溯源 + 几何契约。

    source          analytic | mni152_style | hcp_population | 用户自定义 tag
    anatomy_source  dosimetry v1 tier-4 取值：phantom | subject_imaging | population
    license         数据许可原文（不是代码的 MIT）
    commercial_ok   是否允许商用管线；由 check_license 校验
    layers          由内到外排列的组织层（必须递增半径）
    illustrative    物性/统计是否为示例值（非真实队列）
    provenance      额外溯源（DOI、抓取日期、脚本版本等）
    """
    subject_id: str
    source: str
    license: str
    layers: tuple[TissueLayerGeometry, ...]
    anatomy_source: str = "phantom"
    coordinate_space: str = "subject_ras_cm"
    doi: str = "N/A"
    citation: str = ""
    commercial_ok: bool = True
    illustrative: bool = False
    provenance: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.anatomy_source not in ANATOMY_SOURCE_CHOICES:
            raise ValueError(f"anatomy_source 必须是 {ANATOMY_SOURCE_CHOICES} 之一")
        rs = [lay.outer_radius_cm for lay in self.layers]
        if any(b <= a for a, b in zip(rs, rs[1:])):
            raise ValueError("layers 必须由内到外严格递增外半径")
        if not self.license:
            raise ValueError("license 不能为空——数据许可必须显式声明")

    @property
    def outer_radius_cm(self) -> float:
        return self.layers[-1].outer_radius_cm

    def layer_thicknesses(self) -> dict[str, float]:
        """各层厚度（cm）= 本层外半径 − 内层外半径；最内层从球心算起。"""
        out: dict[str, float] = {}
        for i, lay in enumerate(self.layers):
            if i == 0:
                out[lay.name] = lay.outer_radius_cm
            else:
                out[lay.name] = lay.outer_radius_cm - self.layers[i - 1].outer_radius_cm
        return out

    def skull_thickness_mm(self) -> float:
        for lay, thick in self.layer_thicknesses().items():
            if lay == "skull":
                return thick * 10.0
        return float("nan")

    def as_dict(self) -> dict:
        return {
            "subject_id": self.subject_id,
            "source": self.source,
            "anatomy_source": self.anatomy_source,
            "coordinate_space": self.coordinate_space,
            "license": self.license,
            "commercial_ok": self.commercial_ok,
            "doi": self.doi,
            "citation": self.citation,
            "illustrative": self.illustrative,
            "outer_radius_cm": round(self.outer_radius_cm, 4),
            "layer_thickness_cm": {k: round(v, 4) for k, v in self.layer_thicknesses().items()},
            "skull_thickness_mm": round(self.skull_thickness_mm(), 3),
            "layers": [lay.as_dict() for lay in self.layers],
            "provenance": self.provenance,
            "disclaimer": DISCLAIMER,
        }


# ---------------------------------------------------------------------------
# 许可门禁（fail-closed）
# ---------------------------------------------------------------------------

def check_license(model: HeadModel, *, commercial: bool = False) -> tuple[bool, str]:
    """校验头模许可是否允许当前用途。

    commercial=True 时，只有 :data:`COMMERCIAL_ALLOWED_LICENSES` 白名单通过；
    CC BY-NC / HCP DUA / 未识别许可一律拒绝（fail-closed，不猜测）。
    返回 ``(ok, reason)``，不抛异常，调用方决定是否 halt。

    >>> check_license(_MNI152_STYLE, commercial=True)[0]
    True
    >>> check_license(_ERNIE_EXTENDED, commercial=True)[0]
    False
    """
    if not commercial:
        return True, "非商用/研究用途：按数据原文许可署名与遵守 DUA 即可"
    lic = model.license.strip()
    if "NC" in lic or "DUA" in lic or "非商用" in lic:
        return False, f"许可含非商用/使用条款约束：{lic!r}；不得进入商用管线"
    allowed = [a for a in COMMERCIAL_ALLOWED_LICENSES if lic.startswith(a)]
    if not allowed:
        return False, f"许可未列入商用白名单 {COMMERCIAL_ALLOWED_LICENSES}：{lic!r}；需人工复核"
    if not model.commercial_ok:
        return False, "commercial_ok=False：上游明确禁止商用"
    return True, f"许可 {lic!r} 允许商用（仍需按原文署名）"


def require_license(model: HeadModel, *, commercial: bool = False) -> None:
    """check_license 的 fail-fast 版本，违规直接抛 ValueError。"""
    ok, reason = check_license(model, commercial=commercial)
    if not ok:
        raise ValueError(f"许可门禁拒绝：{reason}")


# ---------------------------------------------------------------------------
# 内置可复现代理体模（不是真实影像数据）
# ---------------------------------------------------------------------------

def _layers_from_radii(radii_cm: Iterable[float],
                       names: Iterable[str] = ("brain", "csf", "skull", "scalp"),
                       shape_ratio: float = 1.0,
                       *,
                       conductivity: dict[str, float] | None = None,
                       attenuation: dict[str, float] | None = None) -> tuple[TissueLayerGeometry, ...]:
    cond = conductivity or DEFAULT_CONDUCTIVITY_S_M
    att = attenuation or DEFAULT_ATTENUATION_DB_CM_MHZ
    return tuple(
        TissueLayerGeometry(name, float(r), cond[name], att[name], shape_ratio)
        for name, r in zip(names, radii_cm)
    )


SPHERICAL_PHANTOM = HeadModel(
    subject_id="analytic_4sphere",
    source="analytic",
    license=SOURCES["analytic"]["license"],
    layers=_layers_from_radii((8.0, 8.5, 9.2, 10.0)),
    anatomy_source="phantom",
    doi=SOURCES["analytic"]["doi"],
    citation=SOURCES["analytic"]["citation"],
    commercial_ok=True,
    provenance={"note": SOURCES["analytic"]["note"]},
)

# 分层椭球代理：外半径 10cm，三轴比 1 : 0.88 : 0.92（量级参考，illustrative）
_MNI152_STYLE = HeadModel(
    subject_id="mni152_style_proxy",
    source="mni152_style",
    license=SOURCES["mni152_style"]["license"],
    layers=_layers_from_radii((7.4, 7.9, 8.6, 9.5), shape_ratio=0.88),
    anatomy_source="phantom",
    doi=SOURCES["mni152_style"]["doi"],
    citation=SOURCES["mni152_style"]["citation"],
    commercial_ok=True,
    illustrative=True,
    provenance={"note": SOURCES["mni152_style"]["note"],
                "shape_ratio_secondary": 0.92},
)

# 群体统计代理：颅骨 6.5mm、脑外半径 8.0cm（illustrative，非 Berger 实测分布）
_HCP_POPULATION = HeadModel(
    subject_id="hcp_population_proxy",
    source="hcp_population",
    license=SOURCES["hcp_population"]["license"],
    layers=_layers_from_radii((8.0, 8.5, 9.15, 10.0)),
    anatomy_source="population",
    doi=SOURCES["hcp_population"]["doi"],
    citation=SOURCES["hcp_population"]["citation"],
    commercial_ok=False,   # HCP DUA：保守起见默认禁止商用管线
    illustrative=True,
    provenance={"note": SOURCES["hcp_population"]["note"],
                "distribution_note": ILLUSTRATIVE_TAG},
)

# CC BY-NC 体模示例：只用于演示许可门禁，不打包任何真实数据
_ERNIE_EXTENDED = HeadModel(
    subject_id="ernie_extended_metadata_only",
    source="siminibs_ernie_extended",
    license=SOURCES["siminibs_ernie_extended"]["license"],
    layers=_layers_from_radii((7.5, 8.0, 8.7, 9.6)),
    anatomy_source="subject_imaging",
    doi=SOURCES["siminibs_ernie_extended"]["doi"],
    citation=SOURCES["siminibs_ernie_extended"]["citation"],
    commercial_ok=False,
    provenance={"note": SOURCES["siminibs_ernie_extended"]["note"] + "；本对象仅含元数据，无真实数据"},
)

BUILTIN_MODELS: dict[str, HeadModel] = {
    "spherical": SPHERICAL_PHANTOM,
    "mni152_style": _MNI152_STYLE,
    "hcp_population": _HCP_POPULATION,
    "ernie_extended_metadata": _ERNIE_EXTENDED,
}


def get_model(name: str) -> HeadModel:
    """取内置代理体模。name 见 :data:`BUILTIN_MODELS`。"""
    if name not in BUILTIN_MODELS:
        raise KeyError(f"未知体模 {name!r}；可选：{sorted(BUILTIN_MODELS)}")
    return BUILTIN_MODELS[name]


# ---------------------------------------------------------------------------
# 群体解剖抽样（给蒙特卡洛剂量区间供数）
# ---------------------------------------------------------------------------

# 示例分布参数——不是真实队列经验值，报告必须带 illustrative 标签。
# 用户可用 load_population_stats() 从自己的 CSV 替换。
EXAMPLE_POPULATION_STATS: dict[str, dict[str, float]] = {
    "brain_outer_radius_cm": {"mean": 8.0, "sd": 0.35},
    "csf_thickness_mm": {"mean": 5.0, "sd": 1.2},
    "skull_thickness_mm": {"mean": 6.5, "sd": 1.4},
    "scalp_thickness_mm": {"mean": 8.0, "sd": 1.5},
    "skull_conductivity_mult": {"mean": 1.0, "sd": 0.35},  # Berger 2025 电导率变异量级
}


def population_sample(*, n: int = 100, seed: int = 20261119,
                      source: str = "hcp_population",
                      stats: dict[str, dict[str, float]] | None = None) -> tuple[HeadModel, ...]:
    """从群体统计代理抽样 n 个头模（固定 seed 可复现）。

    ``stats`` 默认 :data:`EXAMPLE_POPULATION_STATS`（illustrative 示例分布）；
    传入 :func:`load_population_stats` 的返回值即可用真实队列统计替换。
    """
    if n < 1:
        raise ValueError("n 必须 >= 1")
    import numpy as np

    st = stats if stats is not None else EXAMPLE_POPULATION_STATS
    missing = set(EXAMPLE_POPULATION_STATS) - set(st)
    if missing:
        raise ValueError(f"群体统计缺参数：{sorted(missing)}")

    def draw(param: str, lo: float) -> float:
        return max(lo, float(rng.normal(st[param]["mean"], st[param]["sd"])))

    rng = np.random.default_rng(seed)
    base = get_model(source)
    out: list[HeadModel] = []
    for i in range(n):
        brain_r = draw("brain_outer_radius_cm", 5.0)
        csf = draw("csf_thickness_mm", 1.0) / 10.0
        skull = draw("skull_thickness_mm", 1.5) / 10.0
        scalp = draw("scalp_thickness_mm", 2.0) / 10.0
        kmult = draw("skull_conductivity_mult", 0.2)
        radii = (brain_r, brain_r + csf, brain_r + csf + skull,
                 brain_r + csf + skull + scalp)
        cond = dict(DEFAULT_CONDUCTIVITY_S_M)
        cond["skull"] *= kmult
        layers = _layers_from_radii(radii, conductivity=cond)
        out.append(replace(
            base,
            subject_id=f"{source}_s{i:04d}",
            layers=layers,
            illustrative=stats is None,
            provenance={**base.provenance, "sample_index": i, "seed": seed,
                        "distribution_note": ILLUSTRATIVE_TAG if stats is None
                                             else "user_supplied_stats"},
        ))
    return tuple(out)


def load_population_stats(csv_path: str) -> dict[str, dict[str, float]]:
    """从 CSV（列：parameter,mean,sd）加载真实群体统计，替换示例分布。

    不内置任何真实队列数字——这是给用户接入自己数据的钩子。
    """
    import csv
    stats: dict[str, dict[str, float]] = {}
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            stats[row["parameter"].strip()] = {
                "mean": float(row["mean"]), "sd": float(row["sd"]),
            }
    required = set(EXAMPLE_POPULATION_STATS)
    missing = required - set(stats)
    if missing:
        raise ValueError(f"群体统计 CSV 缺参数：{sorted(missing)}")
    return stats


# ---------------------------------------------------------------------------
# 真实数据导入适配器（可选依赖，缺失时给出明确指引，不静默失败）
# ---------------------------------------------------------------------------

def import_nifti_segmentation(nifti_path: str, *, subject_id: str,
                              tissue_labels: dict[str, int] | None = None,
                              license: str = "CC BY 4.0",
                              commercial: bool = False,
                              source: str = "subject_imaging") -> HeadModel:
    """从 NIfTI 分割掩码提取组织体积/等效半径，构造 HeadModel。

    tissue_labels 默认 {brain:1, csf:2, skull:3, scalp:4}；体素间距来自 affine。
    计算每组织体素数 * 体素体积 = 体积，再取等效球半径 (3V/4π)^(1/3)，
    由内到外累加成层。这是**几何粗化**：真实mesh的精细结构不保留，
    报告里以 ``coarsening='equivalent_sphere'`` 标注，不冒充精细个体化。

    需要 nibabel；未安装时抛 ImportError 并给出安装指引。
    """
    require_license_head = HeadModel(  # 仅用于许可门禁校验（layers 临时构造）
        subject_id=subject_id, source=source, license=license,
        layers=_layers_from_radii((1.0, 2.0, 3.0, 4.0)),
        anatomy_source="subject_imaging" if source == "subject_imaging" else "population",
        commercial_ok=not commercial,
    )
    if commercial:
        require_license(require_license_head, commercial=True)

    try:
        import nibabel as nib  # type: ignore
    except ImportError as exc:  # pragma: no cover - 依赖可选
        raise ImportError(
            "import_nifti_segmentation 需要 nibabel；pip install nibabel。"
            "若只需跑通链路，请改用内置代理体模 anatomy.get_model('mni152_style')。"
        ) from exc

    labels = tissue_labels or {"brain": 1, "csf": 2, "skull": 3, "scalp": 4}
    img = nib.load(nifti_path)
    data = img.get_fdata()
    vox_mm = img.header.get_zooms()[:3]
    voxel_volume_cm3 = float(vox_mm[0] * vox_mm[1] * vox_mm[2]) / 1000.0

    volumes_cm3 = {name: int((data == lab).sum()) * voxel_volume_cm3
                   for name, lab in labels.items()}
    if any(v <= 0 for v in volumes_cm3.values()):
        raise ValueError(f"分割掩码缺组织或体素数不足：{volumes_cm3}")

    # 由内到外累加体积 -> 等效球半径
    radii: list[float] = []
    acc = 0.0
    for name in ("brain", "csf", "skull", "scalp"):
        acc += volumes_cm3[name]
        radii.append((3.0 * acc / (4.0 * math.pi)) ** (1.0 / 3.0))
    if any(b <= a for a, b in zip(radii, radii[1:])):
        raise ValueError("等效半径未严格递增，分割标签可能错位")

    return HeadModel(
        subject_id=subject_id, source=source, license=license,
        layers=_layers_from_radii(tuple(radii)),
        anatomy_source="subject_imaging",
        coordinate_space=str(img.affine.dtype),  # 占位：真实空间由 affine 定义
        commercial_ok=not commercial,
        illustrative=False,
        provenance={
            "nifti_path": nifti_path,
            "tissue_labels": labels,
            "voxel_volume_cm3": round(voxel_volume_cm3, 6),
            "tissue_volumes_cm3": {k: round(v, 3) for k, v in volumes_cm3.items()},
            "coarsening": "equivalent_sphere",
            "note": "几何粗化：等效球半径，不保留精细解剖；"
                    "coords 字段仅记录 affine dtype，需用户自行核对空间朝向。",
        },
    )


def load_lead_field(hdf5_path: str, *, subject_id: str,
                    dataset: str = "lead_field") -> dict[str, Any]:
    """读取 lead field 矩阵（HDF5，如 cuNIBS 格式），返回元数据 + 矩阵。

    需要 h5py（可选）。本函数不复制数据、不做许可判断——调用方须先对
    对应 HeadModel 调 :func:`require_license`，因为 HCP 来源数据受 DUA 约束。
    """
    try:
        import h5py  # type: ignore
    except ImportError as exc:  # pragma: no cover - 依赖可选
        raise ImportError("load_lead_field 需要 h5py；pip install h5py") from exc

    with h5py.File(hdf5_path, "r") as fh:
        if dataset not in fh:
            raise KeyError(f"HDF5 中找不到 dataset {dataset!r}；现有：{list(fh.keys())[:10]}")
        matrix = fh[dataset][...]
    return {
        "subject_id": subject_id,
        "dataset": dataset,
        "shape": tuple(matrix.shape),
        "matrix": matrix,
        "note": "lead field 由外部数据提供；本包不打包、不复制、不绕过其数据许可。",
    }


# ---------------------------------------------------------------------------
# 展示
# ---------------------------------------------------------------------------

def render(model: HeadModel) -> str:
    lines = [
        f"头模：{model.subject_id}",
        f"  来源={model.source}  anatomy_source={model.anatomy_source}",
        f"  许可={model.license}  商用={'允许' if model.commercial_ok else '禁止'}",
        f"  DOI={model.doi}  坐标空间={model.coordinate_space}",
        f"  外半径={model.outer_radius_cm:.2f} cm  颅骨={model.skull_thickness_mm():.2f} mm",
        "  分层（内→外，厚度 cm）：",
    ]
    for name, thick in model.layer_thicknesses().items():
        lines.append(f"    - {name:6s} {thick:.3f}")
    if model.illustrative:
        lines.append(f"  ⚠ {ILLUSTRATIVE_TAG}")
    lines.append(f"  注：{DISCLAIMER}")
    return "\n".join(lines)
