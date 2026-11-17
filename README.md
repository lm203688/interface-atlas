# interface-atlas

非药物「调用接口层」的两个可计算地基：**A1 接收器能力图谱** + **A2 物理干预证据与剂量图谱**。

核心命题：物理信号本身不认识细胞。所谓「细胞特异性」来自五层外部约束（空间定位、能量沉积、接收器表达、网络拓扑、剂量窗）。本仓库把其中可计算的 **L3 接收器表达** 与 **证据/剂量/副作用** 做成结构化、可审计、可复现的数据与工具。

## 两个模块

| 模块 | 回答什么 | 数据源 | 产物 |
|---|---|---|---|
| A1 `receptors` | 哪些细胞表达哪些物理信号接收器 | Human Protein Atlas v25.1 单细胞 RNA（155 种细胞类型、20162 基因） | 13×155 表达长表、中文热图、阈值敏感性审计 |
| A2 `evidence` | 什么模态、什么剂量、对什么终点、证据多强、副作用是什么 | 28 条种子证据（已附一手来源） | schema v1.1 CSV、查询器、剂量分层审计、内容指纹 |
| `link` | 文献说的细胞，在数据里真有接收器吗 | A1 + A2 交叉核验 + `aliases.csv` | 匹配状态表，匹配不上就标 UNMATCHED，附证据库指纹 |
| `claimcheck` | 一句设备宣称，证据到底支不支持 | A2 证据库 | 四级判定（支持/反驳/超出/无证据）+ 命中证据 ID + 推翻条件 |

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

## 首批发现（可证伪，非结论）

1. **广布接收器不能宣称细胞特异**：在 ≥1 nCPM 阈值下，COX4I1、CRY1、CRY2、PIEZO1、TRPV1 在超过 80% 的细胞类型中检出。`audit()` 会显式预警。
2. **α7nAChR 的宣称与 RNA 数据有张力**：文献常写「巨噬细胞 α7nAChR 介导抗炎」（A2 记录 VNS-001），但 HPA 中巨噬细胞 CHRNA7 仅 **0.2 nCPM**；最高的是脉络丛上皮细胞（68.9）与肾上腺髓质细胞（54.4）。这可能是 RNA 检测灵敏度、激活状态依赖或蛋白-RNA 不一致，**不等于文献错误**，但足以说明不能把通路叙事当成细胞事实。
3. **SCFA 受体不在肠道独占**：FFAR2 最高表达于中性粒细胞（236.2），肠道内分泌细胞仅排第三（33.1）；「纤维→L 细胞→GLP-1」的叙事只是其中一条通路。
4. **交叉核验通过率仅 39%**：28 条证据中 11 条能对齐到 HPA 细胞类型，17 条因非单一细胞类型、命名歧义或 HPA 无对应细胞而无法核验。工具不硬凑，全部标出。

## 统计口径与边界

- `nCPM` 是归一化每百万计数，**只在同一基因内跨细胞类型可比**，不能跨基因比绝对高低。
- 阈值 1 nCPM 是探索性约定，不是生物学真理；`audit()` 报告 0.5/1/5/10/20 五档的阳性细胞数，阈值敏感性完全披露。
- `specificity_tau` 是类 tau 指数 `(Σ max(x−mean,0)) / Σx`，取值 0–1，越大越集中于少数细胞类型。
- RNA ≠ 膜蛋白丰度；免疫受体表达高度依赖激活状态。本图谱是**假设生成器**，不是疗效证据。
- A2 中缺来源的 3 条 A 级记录（VNS-003/004/005）被显式列为 unresolved，不删除、不静默补链。

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
