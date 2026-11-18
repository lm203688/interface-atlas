# interface-atlas

非药物「调用接口层」的两个可计算地基：**A1 接收器能力图谱** + **A2 物理干预证据与剂量图谱**。

核心命题：物理信号本身不认识细胞。所谓「细胞特异性」来自五层外部约束（空间定位、能量沉积、接收器表达、网络拓扑、剂量窗）。本仓库把其中可计算的 **L3 接收器表达** 与 **证据/剂量/副作用** 做成结构化、可审计、可复现的数据与工具。

## 两个模块

| 模块 | 回答什么 | 数据源 | 产物 |
|---|---|---|---|
| A1 `receptors` | 哪些细胞表达哪些物理信号接收器 | Human Protein Atlas v25.1 单细胞 RNA（155 种细胞类型、20162 基因） | 13×155 表达长表、中文热图、阈值敏感性审计 |
| A2 `evidence` | 什么模态、什么剂量、对什么终点、证据多强、副作用是什么 | **103 条证据**（PBM/VNS/tFUS/机械/代谢，全部附一手 DOI/PMID 且逐条核验） | schema v1.2 CSV、查询器、剂量分层审计、内容指纹 |
| `link` | 文献说的细胞，在数据里真有接收器吗 | A1 + A2 交叉核验 + `aliases.csv` | 匹配状态表，匹配不上就标 UNMATCHED，附证据库指纹 |
| `claimcheck` | 一句设备宣称，证据到底支不支持 | A2 证据库 | 四级判定（支持/反驳/超出/无证据）+ 命中证据 ID + 推翻条件 |
| `paramcheck` | 一份参数报告写得够不够资格谈复现 | 物理参数规格（光/电/声/机械） | 五级判定（完整/不足/越界/无证据/部分）+ 硬阈值核对 + 推翻条件 |
| `dosimetry` | 一份剂量报告缺什么、能不能复现、能不能交给计算引擎补算 | **剂量报告规范 v1**（四层剂量本体，光/电/声/机械 + 代谢附录） | completeness_score + 四档可复现性判定 + `schema/dose_report_v1.schema.json` + 字段缺失图谱 |
| `light_field` | 设备/组织参数 → 靶区真实剂量场（光模态） | P1 扩散近似 + Farrell 1992 偶极子 + 外推边界 | 深度分辨光通量/吸收场、闭式逆设计、可微接口、适用性告警 |
| `mc_light` | 独立金标准校验剂量场、量化代理误差 | 1D 多层 Monte Carlo（numpy 向量化） | Prahl 基准验证、P1 代理误差测量、解剖敏感性实验 |

> **v0.4 战略转向**：项目从「证据审计工具」升级为「物理干预剂量学操作系统」的 L1 层。放疗靠 ICRU 报告标准 + DICOM-RT + TPS 成为精确医学；物理干预仍缺这一层。详见 [`DOSE_REPORTING_SPEC.md`](DOSE_REPORTING_SPEC.md) 与项目战略文档。`paramcheck`（判决）与 `dosimetry`（完整性诊断）正交并存，前者不放宽。
>
> **v0.5 L2 落地**：L1 把「T1+T2 齐缺 T3」判为 `reproducible_with_model`——`light_field` 就是这个可计算的 model。承接「设备参数 ≠ 靶区剂量」这一核心难题，用独立 MC（`mc_light`）作金标准、P1 扩散解作可微代理，误差如实量化。详见 [`L2_VALIDATION.md`](L2_VALIDATION.md)。

## 安装与快速开始

```bash
pip install -r requirements.txt

# 下载 HPA 单细胞 RNA 表（16.4 MB，不入库）
curl -L -o data/rna_single_cell_type.tsv.zip \
  https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip

# A1：重建表达矩阵并审计
python -m interface_atlas.cli receptors --rebuild
python -m interface_atlas.cli receptors --audit
python -m interface_atlas.cli receptors --heatmap data/receptor_heatmap.png

# A2：证据库查询与审计
python -m interface_atlas.cli evidence --init
python -m interface_atlas.cli evidence --audit
python -m interface_atlas.cli evidence --gene CHRNA7 --grade A

# 交叉核验
python -m interface_atlas.cli link

# 宣称-证据判定（v0.2）
python -m interface_atlas.cli claim "红光面板能抗衰老延寿" --markdown
python -m interface_atlas.cli claim "植入式迷走刺激能治疗类风湿关节炎"

# 参数报告预检（v0.3）
python -m interface_atlas.cli paramcheck --modality focused_ultrasound \
  --set frequency=0.5MHz --set mechanical_index=2.4 --set intensity=15W/cm2 \
  --set pulse_duration=0.5ms --set duty_cycle=70% --set temperature_rise=1.2C \
  --cell-type "cortical neurons" --endpoint "motor modulation"
python -m interface_atlas.cli paramcheck --batch   # 批量预检 103 条证据

# 剂量报告完整性审计（v0.4，规范 v1）
python -m interface_atlas.cli doseaudit --modality photobiomodulation \
  --set wavelength=660nm --set power_density=25mW/cm2 --set fluence=6J/cm2 \
  --set exposure_time=240s --set spot_area=1cm2 \
  --set applicator_distance=2cm --set applicator_angle=90deg --set beam_profile=gaussian
python -m interface_atlas.cli doseaudit --batch   # 批量审计：平均完整度 0.114
python -m interface_atlas.cli dosegap             # 重建字段缺失图谱 data/dose_gap_map.json

# L2 剂量场：参数 → 靶区剂量（v0.5，光模态）
python -m interface_atlas.cli dosefield --preset skin_nir --target-depth 0.3
python -m interface_atlas.cli dosefield --preset muscle_nir --target-depth 1.0 --inverse 0.01
python -m interface_atlas.cli mcbench --which all --photons 200000   # 独立验证三件套
```

`data/receptor_expression.csv` 与 `data/evidence_matrix.csv` 已随仓库提供，不下载原始 zip 也能查询；`--rebuild` 才需要原始 zip。

## 接收器面板（13 个基因）

| 模态 | 基因 | 接收器 |
|---|---|---|
| 光生物调节 | COX4I1 | 细胞色素 c 氧化酶亚基 IV |
| 光/昼夜 | CRY1, CRY2 | 隐花色素 |
| 迷走/胆碱能 | CHRNA7 | α7 烟碱型乙酰胆碱受体 |
| 机械力 | PIEZO1, TRPV4 | 机械敏感离子通道 |
| 热 / 冷 | TRPV1 / TRPM8 | 温度敏感通道 |
| 微生物代谢物 | FFAR2, FFAR3 | GPR43 / GPR41（短链脂肪酸受体） |
| 内源 GLP-1 | GLP1R | GLP-1 受体 |
| 胆汁酸 | GPBAR1, NR1H4 | TGR5 / FXR |

## v0.2：从「能查」到「能判」

v0.2 补上路线图的价值拐点——**ClaimCheck 推断层**。输入一句宣称，输出可追溯、可证伪的判定，不输出恒真的"不可信"裁决。

| 判定 | 触发条件 | 示例 |
|---|---|---|
| `SUPPORTED` 证据支持 | 同模态+同终点+同方向，有 A/B 级人体证据 | 「植入式迷走刺激治疗类风湿关节炎」→ VNS-003 |
| `REFUTED` 证据反驳 | 同模态同终点存在 A/B 级人体「无效应」证据 | 「红光面板抗衰老」→ PBM-009 |
| `EXCEEDED` 证据超出 | 仅机制/动物证据，或宣称终点超出实测范围 | 把胶原合成升级成「抗衰延寿」 |
| `NO_EVIDENCE` 无证据 | 库中无对应记录 | 「特定频率电磁波修复细胞」 |

每条判定都带：命中/反证证据 ID、触发规则、**推翻条件**、边界声明、证据库内容指纹。含「抗衰/延寿」的宣称即使命中胶原等替代终点，也不允许升级为 SUPPORTED——这正是最常见的越界宣称。

同版本另外三项升级：
- **schema v1.1**：A 级证据缺一手来源 = 审计硬失败；剂量分四层（数值/相对/物理参数/缺失），来源未报告的剂量标警告而**不编造**；SHA256 内容指纹保证跨版本可复现。
- **三条 A 级来源闭环**：VNS-003 SetPoint FDA PMA P240039、VNS-004 Schiweck 2024（DOI 10.1016/j.bbi.2023.12.008）、VNS-005 Wolf 2021（DOI 10.1111/psyp.13933）全部补一手链接。顺带更正 VNS-005 研究数：正式发表版 16 项（此前 17 项为预印本数字），结论不变。
- **别名外置 + 可核验率**：`aliases.csv` 每条带 basis 依据，替代硬编码；可核验率 32% → 39.3%（11/28）。剩余 17 条多为患者群体/物理屏障/HPA 无对应细胞分类，**不硬凑**。

## v0.3：证据库扩到 103 条 + ParamCheck

**证据库从 28 条扩到 103 条**，覆盖五类模态：

| 模态 | 条目范围 | 重点 |
|---|---|---|
| PBM 光生物调节 | PBM-001~037 | 双相剂量窗、视网膜安全阈值、动物→人穿透率、口腔黏膜炎 Meta |
| VNS 迷走刺激 | VNS-001~018 | 植入式癫痫/抑郁/RA/克罗恩、耳夹消费级、HRV 反证、副作用 |
| tFUS 聚焦超声 | FUS-001~015 | 神经调控、BBB 开放、MI≤1.9 / 温升≤2°C 硬阈值、颅骨衰减 |
| 机械力 | MECH-001~015 | 成骨/软骨/肌腱/内皮、剪切力数值、过载损伤 |
| 代谢接口 | MET-001~018 | SCFA→L 细胞 GLP-1、胆汁酸 TGR5/FXR、内源 GLP-1 |

**来源纪律**：75 条新增证据的 DOI/PMID **逐条经 CrossRef 与 NCBI eutils 核验为真**，虚构链接零容忍。等级分布 A15/B42/C46；42 条带数值剂量，8 条物理参数，1 条相对剂量；52 条来源本身未报告剂量——**如实标缺失，不编造**。schema 升 v1.2.0；内容指纹 `99b2fdd535ca`。交叉核验 14/103 可对齐（13.6%，新条目多为患者群体/组织/生物标志物，不是细胞类型，**不硬凑**）。

**ParamCheck：参数报告预检**。同一个物理干预，参数写法决定它是可复现报告还是pure宣称：

| 判定 | 触发 | 例子 |
|---|---|---|
| `COMPLETE` | 关键参数齐 + 剂量带数值单位 + 未越界 | 660 nm / 25 mW/cm² / 6 J/cm² / 240 s / 1 cm² |
| `PARTIAL` | 参数齐但有字段缺单位等非致命问题 | 剂量有数字无单位 |
| `INSUFFICIENT` | 关键参数不全或剂量只有相对描述 | 「低剂量红光照射」 |
| `OUT_OF_WINDOW` | 越过公开硬安全阈值 | tFUS MI=2.4（>1.9）、温升 3.5°C（>2°C） |
| `NO_EVIDENCE` | 无任何可核对参数 | 「62–72 Hz 频率修复细胞」 |

判定**依被测报告而变**（改一个参数判定就变，有可证伪测试），带推翻条件；硬阈值越界时任何"有效"叙事都不能抵消风险改判。免责：本地自查线索，非监管/临床结论。

## v0.4：剂量报告规范 v1 —— 从审计工具到生态入口

战略转向：不做更大的证据库，做**物理干预剂量学的那一层**。放疗靠 ICRU 报告标准 + DICOM-RT + 治疗计划系统成为精确医学，物理干预尚缺这一层（PEMF 综述因欠报剂量无法做 meta；VNS 连 "neural dose" 都无共识）。

**核心贡献：四层剂量本体**。设备报的参数 ≠ 靶区真实剂量（PNAS 胸部数字孪生证实同一外部剂量在不同解剖下靶区差数倍）：

| 层 | 内容 | 作用 |
|---|---|---|
| T1 设备参数 | 波长/频率/功率/脉宽/幅度 | 方程输入 |
| T2 界面条件 | 施加器距离/角度、线圈几何、颅骨衰减、电极位置 | 方程边界条件（文献最常漏报） |
| T3 靶区剂量 | 靶区 fluence/声压/电场、组织光学声学参数 | 方程的解 |
| T4 个体化解剖 | 群体值/体模/个体影像 | 可否个体化 |

由此得到四档判定：`not_reproducible`（T1 缺）→ `underdetermined`（T1 齐 T2 缺，方程欠定）→ `reproducible_with_model`（T1+T2 齐，可交 L2 计算引擎求解）→ `reproducible_from_report`（T1/T2/T3 全齐）。**判定依报告逐层变化，有可证伪测试。**

**实测基线**（`dosegap` 生成，103 条）：平均完整度 **0.114**，最高 0.40；99 条 `not_reproducible`、4 条 `underdetermined`、`reproducible_*` **0 条**。缺失最严重的是 T2 界面条件与 T3 靶区剂量——PBM 37 条全不报施加器距离，tFUS 15 条全不报换能器几何/颅骨衰减/靶区声压，85 条主表记录全不报解剖来源。**这不是数据不够，是报告标准缺位的量化证据。**

产物：[`DOSE_REPORTING_SPEC.md`](DOSE_REPORTING_SPEC.md)（人类可读规范）、`schema/dose_report_v1.schema.json`（机器可读，JSON Schema draft 2020-12）、`data/dose_gap_map.json`（实证图谱）。规范是社区草案 v1，非 ICRU/FDA/ISO 官方标准；分数度量完整性，不度量疗效或安全。`paramcheck` 保留原有五级判决不放宽，与本模块正交并存。

## 首批发现（可证伪，非结论）

1. **广布接收器不能宣称细胞特异**：在 ≥1 nCPM 阈值下，COX4I1、CRY1、CRY2、PIEZO1、TRPV1 在超过 80% 的细胞类型中检出。`audit()` 会显式预警。
2. **α7nAChR 的宣称与 RNA 数据有张力**：文献常写「巨噬细胞 α7nAChR 介导抗炎」（A2 记录 VNS-001），但 HPA 中巨噬细胞 CHRNA7 仅 **0.2 nCPM**；最高的是脉络丛上皮细胞（68.9）与肾上腺髓质细胞（54.4）。这可能是 RNA 检测灵敏度、激活状态依赖或蛋白-RNA 不一致，**不等于文献错误**，但足以说明不能把通路叙事当成细胞事实。
3. **SCFA 受体不在肠道独占**：FFAR2 最高表达于中性粒细胞（236.2），肠道内分泌细胞仅排第三（33.1）；「纤维→L 细胞→GLP-1」的叙事只是其中一条通路。
4. **交叉核验可核验率低**：v0.2 时 11/28（39%）；v0.3 扩库后为 14/103（13.6%）——新条目大量是患者群体、组织与生物标志物，本就不是 HPA 细胞类型。工具不硬凑，UNMATCHED 一律标出；低可核验率本身就是「物理干预宣称难以落到具体细胞」的量化证据。

## 统计口径与边界

- `nCPM` 是归一化每百万计数，**只在同一基因内跨细胞类型可比**，不能跨基因比绝对高低。
- 阈值 1 nCPM 是探索性约定，不是生物学真理；`audit()` 报告 0.5/1/5/10/20 五档的阳性细胞数，阈值敏感性完全披露。
- `specificity_tau` 是类 tau 指数 `(Σ max(x−mean,0)) / Σx`，取值 0–1，越大越集中于少数细胞类型。
- RNA ≠ 膜蛋白丰度；免疫受体表达高度依赖激活状态。本图谱是**假设生成器**，不是疗效证据。
- v0.2 曾缺来源的 3 条 A 级记录（VNS-003/004/005）已补一手 DOI/FDA 链接；v0.3 全部 103 条来源均经 CrossRef/NCBI 核验，A 级缺源硬失败仍为审计退出门。
- v0.3 的 52 条缺剂量是「来源本身未报告」，不是遗漏；ParamCheck 会把这类报告判为 INSUFFICIENT/NO_EVIDENCE 而非编造数值。

## 证据等级

| 等级 | 含义 |
|---|---|
| A | 系统评价/Meta、监管批准、大型队列一致结果 |
| B | 单项 RCT 或中等确定性系统评价 |
| C | 动物/体外机制清楚，人体终点不可靠 |
| D | 概念、假说或营销宣称 |

## 测试

```bash
PYTHONPATH=. python -m pytest tests -q
```

## 许可与引用

代码 MIT。HPA 数据为 **CC BY 4.0**，允许研究与商业使用但须署名（2026-11-16 核实原文）。
使用本仓库的 A1 数据须引用：Karlsson M et al., *A single-cell type transcriptomics map of human tissues*, Sci Adv (2021), DOI 10.1126/sciadv.abh2169，并注明 HPA v25.1 与下载 URL。

本仓库不构成医疗建议。
