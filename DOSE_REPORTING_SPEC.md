# 物理干预剂量报告规范 v1（社区草案）

> Physical Intervention Dose Reporting Spec v1 — interface-atlas
>
> 状态：**社区草案 v1**，非 ICRU / FDA / ISO 官方标准。欢迎引用、复现、挑错。
> 机器可读版本：`schema/dose_report_v1.schema.json`（JSON Schema draft 2020-12）
> 实证依据：`data/dose_gap_map.json`（由 `python -m interface_atlas.cli dosegap` 生成）

## 0. 为什么需要这份规范

放疗成为精确医学靠的不是更强的加速器，而是三件事：**ICRU 剂量报告标准**（统一"剂量"是什么意思）、**DICOM-RT**（机器可读的剂量流转）、**治疗计划系统 TPS**（从设备参数算到靶区剂量）。

物理干预（光 / 电 / 声 / 机械）今天处在放疗 1970 年代的位置：有设备、有疗效信号，但没有剂量学标准、没有 TPS、没有可复现的处方。后果是可测量的：

- PBM 伞形综述（15 个 meta、9000+ 人）：35 个终点**零个达高确定性**，65.7% 无显著效应，剂量学异质性是主因。
- PEMF/EMTT 系统综述（Biomedicines 2026）：多数试验不报线圈几何、脉冲形状、施加器距离与角度，**无法做 meta 分析**，作者喊出 "urgent call to action for standardized reporting"。
- 神经调控闭环综述（2026）：连 "neural dose" 的术语、归一化 ECAP 幅值、剂量精度都无共识。

本规范要解决的就是这个断点。它不宣称疗效，只规定**要让一份物理干预报告可被复现、可被跨研究比较，至少要报什么**。

## 1. 核心创新：四层剂量本体（TIER）

"设备报的参数 ≠ 靶区真实剂量"。同一台设备、同一组旋钮，在不同解剖下靶区剂量可以差数倍（PNAS 胸部数字孪生已定量证实）。因此剂量报告必须分四层：

| 层 | 名称 | 定义 | 谁提供 |
|---|---|---|---|
| **T1** | 设备参数 device | 设备面板/说明书直接给出的量：波长、频率、功率、脉宽、幅度 | 设备 |
| **T2** | 界面条件 interface | 从设备到靶区的传递条件：施加器距离、入射角、线圈几何、颅骨衰减、电极位置、组织参数 | 实验者（**文献最常漏报**） |
| **T3** | 靶区剂量 target | 真正到达靶区的剂量场：靶区 fluence / 声压 / 电场、靶区深度、组织光学/声学/电学参数 | 测量或计算 |
| **T4** | 个体化解剖 subject | 用群体平均值、体模，还是个体影像 | 实验者 |

**这四层构成一个方程系统**：T1 + T2 是剂量场方程的完备边界条件，T3 是这个方程的解。由此得出可复现性的四档判定：

| 判定 | 含义 | 条件 |
|---|---|---|
| `not_reproducible` | 不可复现，等同宣称 | T1 有字段完全缺失 |
| `underdetermined` | 方程欠定 | T1 齐、T2 缺（界面条件不足） |
| `reproducible_with_model` | 可由剂量场引擎求解 | T1+T2 齐、T3 缺 |
| `reproducible_from_report` | 报告自带靶区剂量，直接可复现 | T1/T2/T3 全齐 |

这个判定的意义：它把"能不能补算"形式化了。T1+T2 齐但 T3 缺，不是"数据烂"，而是"边界条件完备、可以交给计算引擎求解"——这正是治疗计划系统的定位，也是本项目 L2 的工作。

## 2. 各模态必报字段

单位允许等价写法（如 `660 nm` / `wavelength=660nm`）。完整字段与别名见 `interface_atlas/dosimetry.py` 的 `DOSE_SPEC`。

### 2.1 光生物调节 photobiomodulation

| 层 | 字段 | 单位 |
|---|---|---|
| T1 | 波长 wavelength | nm |
| T1 | 辐照度/功率密度 power_density | mW/cm² |
| T1 | 能量密度 fluence | J/cm² |
| T1 | 照射时长 exposure_time | s |
| T1 | 光斑面积 spot_area | cm² |
| T2 | 施加器-皮肤距离 applicator_distance | cm |
| T2 | 入射角 applicator_angle | deg |
| T2 | 光束分布 beam_profile | uniform / gaussian |
| T3 | 靶区实际光通量 target_fluence | J/cm² |
| T3 | 靶区深度 target_depth | cm |
| T3 | 组织光学参数 tissue_optics | μa / μs (1/cm) |
| T4 | 解剖来源 anatomy_source | phantom / subject_imaging / population |

### 2.2 迷走神经刺激 vagus nerve stimulation

T1：频率(Hz)、脉宽(µs)、电流(mA)、单次时长(min)、部位(cervical/auricular/abdominal/other)。
T2：电极几何与触点配置、靶纤维束、接口阻抗(kΩ)。
T3：神经剂量（ECAP 归一化，mV）、靶器官 engagement（如脾神经激活）。
T4：解剖来源（population / subject_imaging）。

注：VNS 的 "neural dose" 目前学界无共识，本规范采用可观测代理量（ECAP 幅值）+ 靶器官功能读数双报，不强行定义单一数值。

### 2.3 聚焦超声 focused ultrasound

T1：频率(MHz)、机械指数(MI)、声强(W/cm² 或 mW/cm²，注明 ISPTA/ISPPA)、脉冲宽度、占空比(%)、温升(°C)。
T2：换能器几何、颅骨衰减(dB)、声窗。
T3：靶区声压(MPa)、靶区温度(°C)、组织声学参数(衰减 dB/cm/MHz + 声速)。
T4：解剖来源。

硬阈值（非显著风险边界，越界改判）：MI ≤ 1.9，温升 ≤ 2°C（FDA / ITRUSST）。**注意：这是改判条件，不是疗效保证。**

### 2.4 机械载荷 mechanical loading

T1：应变/剪切应力(%/Pa/dyne·cm⁻²)、加载频率(Hz)、单次时长(min)、基质刚度(kPa)。
T2：加载波形(static/cyclic/impact)、接触几何。
T3：靶区实际应变(%)、组织粘弹性参数。
T4：解剖来源。

### 2.5 附录：代谢 / 微生物（不纳入主评分）

这类模态无传统剂量单位，用暴露量替代：摄入量(g/day 或 CFU 或 mmol/L)、暴露时长、化合物/菌株身份、靶部位浓度。v1 只定义字段，不纳入 completeness_score 主评分——因为"剂量"语义与物理能量不同，强行合并会产生假可比性。

## 3. 完整性分数（completeness_score）

每个字段的状态与得分：

| 状态 | 含义 | 分 |
|---|---|---|
| PRESENT | 带数值和可核对单位 | 1.0 |
| NO_UNIT | 有数字无单位（跨研究不可比） | 0.5 |
| RELATIVE | 只有"低剂量/适量"等相对描述 | 0.25 |
| MISSING | 完全缺失 | 0.0 |

层权重：T1=0.40、T2=0.25、T3=0.25、T4=0.10。总分 = Σ(层均分 × 层权重)，范围 0–1。

**分数只度量报告完整性，不度量疗效或安全性。** 一份完整度 1.0 的报告仍可能无效；一份 0.1 的报告也不能推出"无效"——只能推出"不可复现"。

## 4. 实测基线与改进目标

对现有 103 条证据库逐条审计（`data/dose_gap_map.json`）：

| 指标 | 实测值 |
|---|---|
| 平均完整度 | **0.111** |
| 最高完整度 | 0.40 |
| `not_reproducible` | 99 / 103 |
| `underdetermined` | 4 / 103 |
| `reproducible_with_model` / `from_report` | **0 / 103** |

缺失最严重的层：**T2 界面条件与 T3 靶区剂量几乎全缺**——光生物调节 37 条全部不报施加器距离；聚焦超声 15 条全部不报换能器几何、颅骨衰减、靶区声压与温度；全部 85 条主表记录不报解剖来源。

这不是"数据不够"，是**报告标准缺位**的直接量化。规范 v1 的目标不是让旧文献变好，而是让**新报告从一开始就能达到 `reproducible_with_model`**（即 T1+T2 齐），从而可进入计算引擎求解。

## 5. 如何使用

```bash
# 审计单份报告（人类可读）
python -m interface_atlas.cli doseaudit --modality photobiomodulation \
  --set wavelength=660nm --set fluence=6J/cm2 --set exposure_time=240s

# 输出 JSON（供机器消费）
python -m interface_atlas.cli doseaudit --modality focused ultrasound --json ...

# 批量审计证据库 → 缺失图谱
python -m interface_atlas.cli dosegap

# Python API
from interface_atlas import dosimetry
rep = dosimetry.audit({"id": "X", "modality": "photobiomodulation", "wavelength": "660 nm"})
print(rep.reproducibility, rep.completeness_score, rep.missing_params)
```

## 6. 引用与版本

- 版本：v1.0.0（2026-11-18）
- 仓库：https://github.com/lm203688/interface-atlas
- 许可证：MIT
- 变更规则：字段增减属于 minor；层结构或判定语义变更属于 major。任何变更必须同步更新 JSON Schema、本规范和 `dose_gap_map.json`。

## 7. 免责

本规范是社区草案，不构成监管、临床或安全结论。参数落在窗口内不代表有效，越界不代表必然有害。使用者应自行做 IRB / 设备安全评估。
