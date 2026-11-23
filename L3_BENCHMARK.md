# L3 剂量-效应公开基准（种子版）与个体化解剖层

> 版本 v0.9.0 · 2026-11-19 · 社区草案，非临床指南
>
> 本版做了两件事：(1) 把 L2 从解析球推进到**个体化解剖层**（接公开头模契约 + 群体抽样 + 许可门禁）；(2) 建立 **L3 剂量-效应种子基准**——它的目的不是排名模型，而是统一"剂量该怎么报、效应该怎么比"的口径，并让每条结论可被证伪。

---

## 零、为什么这一层必须存在：三层因果断裂

物理干预领域有个被反复观察到的尴尬：**设备报的参数一样，靶区剂量差几倍；靶区剂量一样，生物效应又不一致**。2025–2026 年的一手研究把这条链路上的断裂点钉死了：

| 断裂 | 一手实证锚点 | 对本项目的含义 |
|---|---|---|
| **① 模拟剂量 ≠ 体内剂量** | Kop et al. 2026 (bioRxiv, `10.64898/2026.09.27.754841`)：24 人双盲，**MR-ARFI 实测组织位移预测了 LGN/V1 的 BOLD 效应，而模拟 in-situ 声强与位移、效应都不相关** | L2 输出的声压/场强点估计必须配体内剂量锚点（MR-ARFI / ECAP / 实测），否则"剂量-效应"结论不可信。颅骨与声耦合差异未被模拟 |
| **② 体内剂量 ≠ 生物效应** | Spivak et al. 2026 (*Neuromodulation*, `10.1016/j.neurom.2026.03.015`)：tFUS 升到 10.08 W/cm²（14× FDA 诊断限）仍无不良事件，但剂量-效应 **mixed and inconclusive**（n=10） | 剂量上去效应没跟上真实存在。基准必须容纳**非单调**剂量-效应，禁止预设线性 |
| **③ 暴露度量本身多维** | Ghaempanah et al. 2025 (BeNe Brain Stimulation)：换 CSF 分割管线，电场峰值 0.236→0.224 V/m（~5%，p=0.071），但聚焦体积 10.9→6.43 cm³（**~39%**，p=3.6e−14, r=0.87） | 单一标量"剂量"会漏掉最重要的差异。基准强制多指标：峰值+聚焦+覆盖+分布 |

**结论**：任何声称"X 剂量 → Y 效应"的工具，如果只给一个标量、不报体内剂量锚点、不报聚焦性、不允许非单调，就与该领域已知事实矛盾。本基准按此设门禁。

---

## 一、L2 个体化解剖层（`anatomy.py`）

### 1.1 数据契约

`dosimetry` 剂量报告规范 v1 把 `anatomy_source` 定义为 tier-4 字段（`phantom | subject_imaging | population`），但此前没有模块消费它。本层补上：

- `TissueLayerGeometry`：组织层的半径、电导率（喂电场）、声衰减（喂超声）、椭球轴比。
- `HeadModel`：subject_id / source / anatomy_source / coordinate_space / DOI / citation / **license / commercial_ok** / illustrative / provenance。

**头模必须显式声明数据许可**，未声明即构造失败（fail-closed）。

### 1.2 内置代理体模（不是真实影像，绝不冒充）

| 体模 | 来源锚点 | 颅骨 | 许可 | 商用 |
|---|---|---|---|---|
| `spherical` | 与 `electric_field.HEAD_4SPHERE` 同构的解析同心球 | 7.0 mm | MIT | 允许 |
| `mni152_style` | 借 Puonti 2020 (`10.1016/j.neuroimage.2020.117044`) 分层尺度的椭球代理 | 7.0 mm | MIT（代理）/ CC BY 4.0（原法） | 允许 |
| `hcp_population` | Berger 2025 (*Sci Data*, `10.1038/s41597-025-04886-0`) 100 HCP 头模的群体统计代理 | 6.5 mm | CC BY 4.0 + **HCP DUA** | 禁止 |
| `ernie_extended_metadata` | 仅元数据占位（Van Hoornweder 2024, `10.1162/imag_a_00379`） | — | **CC BY-NC 4.0** | 禁止 |

后两者明确标注 `illustrative`——群体分布参数是**示例值，不是从真实队列提取的经验值**；真实 100 头模（9.9 GB）由用户自行从 Zenodo `10.5281/zenodo.13259679` 获取并遵守 HCP DUA。

### 1.3 导入适配器（真实数据由用户接入）

- `import_nifti_segmentation()`：从 NIfTI 分割掩码按体素体积算等效球半径，构造个体化 `HeadModel`。**几何粗化**（标 `coarsening=equivalent_sphere`），不冒充精细头模。需 `nibabel`（可选依赖）。
- `load_lead_field()`：读取 cuNIBS 风格 HDF5 lead-field 矩阵，只返回元数据+矩阵，不复制数据、不绕过许可。
- `load_population_stats()`：从 CSV 加载真实群体统计替换示例分布。
- `check_license(model, commercial=...)`：**fail-closed** 许可门禁——CC BY-NC / DUA / 未识别许可一律拒绝商用。四个内置体模的门禁行为均已实测。

### 1.4 群体抽样

`population_sample(n, seed)` 从群体分布抽 n 个头模（固定 seed 可复现），为蒙特卡洛剂量区间供数——接上 v0.6 的 conformal 思路，让"点估计+区间"扩展到个体化维度。抽样实测：颅骨厚度跨度 >1 mm、落在生理区间 2–15 mm，同 seed 完全一致。

---

## 二、L3 剂量场特征层（`dose_features.py`）

### 2.1 多指标向量（拒绝单一 MAE）

`DoseFeatureVector` 每次输出都含：**峰值 + 轴向衰减 + 穿透深度 + 半高宽 + 分层覆盖 + 聚焦性 proxy + 逐项可信度标注**。设计直接回应断裂③——Ghaempanah 2025 的 39% 聚焦性差异，在只报峰值的体系里完全消失。

| 模态 | 峰值单位 | 聚焦性来源 | 诚实的不可用项 |
|---|---|---|---|
| 光 P1 | 相对光通量（表面=1） | — | 聚焦性 unavailable（无横向场，需 MC） |
| 超声 | MPa 峰值声压 | 轴向焦域/焦距（需传 transducer，否则 unavailable） | 颅骨折射/相位像差未建模 |
| 电 | V/m 靶区场 | 球心场/靶区场比（proxy） | 皮层横向聚焦 unavailable（需 FEM/SimNIBS） |

**关键诚实点**：
- 每个指标标 `analytic | proxy | unavailable`；解析引擎不返回完整 3D 场，聚焦性只能是 proxy，绝不冒充数值场体积分。
- 跨模态峰值单位不同（V/m、MPa、相对光通量），**不做虚假归一化排名**。
- `from_numeric_field_full()` + `field_import`：用户把 SimNIBS / k-Wave / NPZ 真解接入后，才算真正可信的聚焦体积分（已用单峰值体素与高斯椭球解析交叉核对）。覆盖体积进 `coverage_volumes_cm3`（cm³），穿透/衰减/FWHM 按是否显式选轴分级——自动选轴只标 proxy，不冒充剂量学轴。导入契约见 `FIELD_IMPORT.md`。

### 2.2 群体汇总

`cohort_summary()` 给每个指标的均值/标准差/分位数/极值比，明确标 `single_metric_leaderboard_forbidden=True`——只给分布，不聚合成一个分数。

---

## 三、L3 种子基准（`l3bench.py`）

### 3.1 四个种子案例

| 案例 | 模态 | 剂量锚点 | 效应发现 | 证据级别 |
|---|---|---|---|---|
| **B1 THIELSCHER-2026** | tES | 0.2 V/m（43 人×8 靶区） | 统一场强框架 + 跨样本泛化，已用于 ~1000 数据集 | peer-reviewed (`10.1016/j.neuroimage.2026.121882`) |
| **B2 KOP-2026** | tFUS | 500 kHz, 50 W/cm², MR-ARFI 位移 | 位移预测 BOLD，模拟强度不预测（24 人双盲） | preprint (`10.64898/2026.09.27.754841`) |
| **B3 SPIVAK-2026** | tFUS | 0.72–10.08 W/cm²（14× FDA 限） | 安全但剂量-效应 inconclusive（n=10） | peer-reviewed (`10.1016/j.neurom.2026.03.015`) |
| **B4 GHAEMPANAH-2025** | tES | 1 mA 额部 montage（47 人） | CSF 分割使聚焦体积 -39% 而峰值 ~不变 | conference poster |

### 3.2 准入硬门

每条案例**必须**有：非空的测量手段、涉及哪条断裂、应报指标、来源 DOI（同行评议必须真实 DOI；预印本/会议可显式标 N/A 并注明场合）、以及**可证伪推翻条件**（≥20 字，明确"什么新证据会推翻它"）。缺推翻条件的案例拒绝收录——这是本基准与"摘抄文献做宣传"的根本区别。

### 3.3 评测协议（多指标，拒绝刷分）

- 提交方按 `peak_accuracy / focality_agreement / coverage_calibration / dose_response_shape / anatomy_sensitivity` 五类口径报数。
- **缺任一必需指标 → `passed=None`（证据不足，不猜、不聚合）**。
- 逐项过阈值，**不做加权求和**（避免单一分数掩盖维度差异）。
- `aggregate_score` 恒为 `None`——明确拒绝"谁最准"的排行榜。

---

## 四、本轮验证与门禁

独立验证 `l3_validation.py` 用**与被测代码不同的实现路径**交叉核对，8 项全 PASS：

| 检验 | 方法（可证伪点） | 结果 |
|---|---|---|
| `geometry_conservation` | 独立重算层厚和 == 外半径（几何恒等式） | PASS |
| `license_gate` | 构造 CC BY-NC / DUA 数据，验证商用被拒、研究放行 | PASS |
| `population_sampling` | 同 seed 复现 + 不同 seed 分化 + 生理区间 | PASS |
| `feature_multimetric` | 多指标齐出 + 逐项可信度 + 跨模态不归一化 | PASS |
| `numeric_field_focality` | 单峰值体素体积占比的解析值 vs 体素求和 | PASS |
| `benchmark_protocol` | 空提交全 N/A + 越界拒过 + 无聚合分 | PASS |
| `three_gaps_anchored` | 三条断裂必须有实证锚点与 DOI（Kop/Spivak 必在） | PASS |
| `l3_benchmark_selfcheck` | 基准自检报告整体通过 | PASS |

全套测试 **123 passed**（原 97 + 新 26）；凭证扫描 CLEAN。机器报告：`data/l3_validation_report.json`。

---

## 五、边界声明（不可省略）

1. 本基准**不是临床指南或疗效结论**，任何条目不得解读为治疗建议。
2. 内置体模是**解析/统计代理**，不是真实患者模型；真实头模需用户自行接入并遵守数据许可。
3. 群体分布默认参数是 **illustrative 示例值**，不是经验分布；真实数据用 `load_population_stats()` 替换。
4. 解析引擎的聚焦性是 proxy，只有接入数值场才算真实体积分。
5. 案例 B2 是**预印本**、B4 是**会议摘要**，均已在 caveat 中标注，不作为确定结论。
6. 三层断裂意味着：**从"模拟剂量"直接跳到"疗效"的任何推断都不被本基准支持**。

---

## 六、可复现命令

```bash
# 内置体模 / 群体抽样 / 许可门禁
python -m interface_atlas.cli headmodel list
python -m interface_atlas.cli headmodel sample --n 100 --seed 42
python -m interface_atlas.cli headmodel license --model ernie_extended_metadata --commercial

# 三模态多指标特征
python -m interface_atlas.cli dosefeatures light --depth 1.0
python -m interface_atlas.cli dosefeatures ultrasound --skull
python -m interface_atlas.cli dosefeatures electric --current 1.0 --depth 2.2

# L3 基准
python -m interface_atlas.cli l3bench cases
python -m interface_atlas.cli l3bench validate
```

MIT License. 案例文字为对一手研究的事实转述，版权归原作者；本包仅索引 DOI，不复制全文。
