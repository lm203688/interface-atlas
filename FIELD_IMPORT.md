# 数值真解导入契约（P2）

> 让外部真解（SimNIBS 电场 / k-Wave 声场）进入本包口径：读得进、单位明确、指标可核对、数据不入库。

## 1. 为什么需要导入层

本包的 L2 解析引擎（光 P1、超声轴上 RS、电四球球谐）**不返回完整 3D 场**，因此聚焦性/覆盖体积只能是 proxy。真正的聚焦性需要数值真解。但"有真解"和"能用统一口径比较"是两件事——不同工具的场文件格式、体素间距、单位、缩放各不相同。没有导入契约，外部研究者不会采用我们的字段定义。

本层把接入标准化：

```text
SimNIBS NIfTI / k-Wave MAT(HDF5) / 开放 NPZ
  → NumericField（显式单位 + 各向异性体素 + provenance）
  → from_numeric_field_full（全指标）
  → DoseFeatureVector
```

## 2. 三条入口与依赖

| 入口 | 函数 | 依赖 | 说明 |
|---|---|---|---|
| 开放 NPZ | `from_npz` | numpy（已有） | **推荐**。格式见 §3，零新增依赖 |
| SimNIBS NIfTI | `from_simnibs_nifti` | nibabel（可选） | 读 `|E|` 幅值图；体素间距取 `header.get_zooms()[:3]`，mm→cm 各向异性保留 |
| k-Wave MAT | `from_kwave_mat` | scipy.io（v5）/ h5py（v7.3，可选） | 场变量名与体素间距**必须由调用方提供** |
| 内存数组 | `from_array` | numpy | 已解析好的 3D 数组 |

依赖缺失时，对应适配器抛 `ImportError` 并给出安装指引或"导出 NPZ 再读"的替代路径；**判定内核、CLI 其余命令、MCP/Web 入口完全不受影响**。

## 3. 开放 NPZ 契约（推荐外部工具优先导出）

```
field          3D 幅值数组（必需），如 |E| 或 |p|
voxel_size_cm  标量或三元 (sx,sy,sz)，单位 cm（必需）
unit           来源单位字符串（可选；调用方仍须显式确认）
source_tool    来源标识（可选）
provenance     JSON 字符串（可选；许可、affine 等）
```

`from_npz(path, field_unit=...)` 的 `field_unit` **必须显式传**，且与内嵌 `unit` 不一致时会告警（以调用方为准）。`save_npz(field, path)` 可把 `NumericField` 写回该契约供其他工具消费。

## 4. 三条铁律（诚实边界）

1. **不猜单位**。`field_unit` 必填。SimNIBS 常为 V/m、k-Wave 常为 Pa，但导出脚本可能预缩放——猜单位等于伪造剂量。需要缩放用 `scale`（如 Pa→MPa 传 `1e-6`），并记入 provenance。
2. **不复制外部数据**。适配器只读路径、计算即用完弃数组，不把受许可约束的真解写进本仓。调用方须先对头模做许可核验（`anatomy.require_license`），并把许可声明传入 `license_declared`。
3. **真解 ≠ 体内剂量**。数值场是特定网格/边界条件/组织参数下的解；电导率个体差异（~2.1×）、颅骨建模误差、网格分辨率都不是本层能消除的。告警如实保留。

## 5. 全指标的诚实分级

`from_numeric_field_full` 从真解算：峰值、多分位覆盖体积、聚焦性（半高体积/总体积）、穿透深度、轴向衰减、FWHM。分级：

| 指标 | 来源 | 标注 |
|---|---|---|
| 峰值 / 覆盖体积 / 聚焦性 | 3D 体积分，各向异性 `sx·sy·sz` | `analytic` |
| 穿透深度 / 衰减 / FWHM | 过峰值的 **1D 剖面** | 显式 `--axis` → `analytic`；自动选轴 → `proxy`；无交点 → `unavailable` |

**为什么轴要分级**：3D 场没有天然"深度轴"，穿透深度只在某条剖面上定义。自动选跨度最长的轴只是粗看，不是剂量学轴——要冻结口径必须显式给定 `axis`（0/1/2）。覆盖体积写入新字段 `coverage_volumes_cm3`（cm³），不再占用 `coverage_depths_cm` 的深度语义。

## 6. CLI 用法

```bash
# 先看导入元数据（单位、体素、缩放、许可、告警），不抽特征
python -m interface_atlas.cli fieldimport npz field.npz --unit V/m --describe

# 抽全指标：显式轴 2（冻结口径），输出 JSON
python -m interface_atlas.cli fieldimport simnibs E_mag.nii.gz \
  --unit V/m --subject S01 --license "CC BY 4.0" \
  --axis 2 --modality "transcranial electric" --as-json

# k-Wave MAT v5：变量 p，Pa→MPa，体素各向异性
python -m interface_atlas.cli fieldimport kwave p_final.mat \
  --unit MPa --variable p --scale 1e-6 --voxel 0.1,0.1,0.2 \
  --axis 2 --modality "focused ultrasound"
```

## 7. 与 L1 条文的关系

导入保留来源场单位，不做跨模态归一化。要按 L1 规范审计某条报告是否报全了界面条件（T2）/靶区剂量（T3），用 `doseaudit`；本层只负责把真解变成可比较的多指标特征。跨模态峰值单位不同，禁止直接横向排名——这是 L3 基准的硬规则。

## 8. 免责

本层是数据导入与特征抽取工具，不输出处方、不判定疗效或安全。`analytic` 只表示"从给定场严格计算"，不代表该场准确反映体内真实剂量。
