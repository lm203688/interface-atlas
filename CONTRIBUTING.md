# 贡献指南：30 分钟提交第一条证据

欢迎外部贡献。这个项目最需要的是**真实、可追溯的证据条目**——不是更多观点。下面是一条最小路径，从零到可合并 PR。

---

## 0. 先懂红线

- 每条证据必须能追溯到原文（DOI / PMID / PMCID，或明确标注机构页/指南页）。
- **AI 抽取的候选不算证据**。机器只做 schema、去重、来源可达性核验；科学真伪由你核对原文后签字。
- 不接受"营销宣称 + 推断"当证据；反证（无效应/负结果）与正结果同样欢迎，而且更稀缺。
- 不提供临床建议，不替任何产品背书。

---

## 1. 环境（2 分钟）

```bash
git clone https://github.com/lm203688/interface-atlas.git
cd interface-atlas
python -m pip install -e . pytest
python -m pytest tests -q    # 应先全绿
```

---

## 2. 写一条候选（10 分钟）

新建 `candidates.json`，字段对齐证据库 `data/evidence_matrix.csv`：

```json
[
  {
    "id": "PBM-013",
    "modality": "photobiomodulation",
    "signal": "635 nm red",
    "receiver_gene": "COX4I1",
    "receiver": "CCO",
    "cell_type": "human fibroblast",
    "endpoint": "collagen synthesis",
    "effect_direction": "increase",
    "effect_size": "moderate (SMD 0.5)",
    "dose": "5",
    "dose_unit": "J/cm2",
    "model_system": "human RCT",
    "human_evidence": "yes",
    "evidence_grade": "B",
    "adverse_event": "",
    "safety_boundary": "biphasic dose response",
    "source_url": "https://doi.org/10.xxxx/yyyy",
    "notes": "原文核对记录：表 2，n=42，2026-11-21"
  }
]
```

字段约束：
- `effect_direction` ∈ `increase | decrease | physical_constraint | mixed | none | biphasic`
- `evidence_grade` ∈ `A | B | C | D`
- `human_evidence` ∈ `yes | no | mixed`
- `source_url` 必须是 http(s) 链接，优先 DOI / PubMed / PMC。

**不要直接编一个 DOI 让它看起来像真的**——产线会查，维护者也会查。

---

## 3. 跑机器裁判（2 分钟）

```bash
# 先离线：只查 schema 与去重，不发网络请求
python -m interface_atlas.cli harvest --input candidates.json --offline

# 再在线：核验 DOI / PMID / PMCID 是否真实存在
python -m interface_atlas.cli harvest --input candidates.json --save data/harvest_review_queue.json
```

三态含义：

| 结果 | 含义 | 你要做什么 |
|---|---|---|
| `ready_for_review` | 来源存在、无重复 | **仍需你核对原文**才能提交 |
| `needs_manual_verification` | 来源不受公开 API 支持或限流/离线 | 亲自打开链接确认，在 notes 写核验记录 |
| `rejected` | schema 错、重复、或来源明确 404 | 修正或放弃；不要换个体裁重新塞 |

注意：CrossRef 经常 429 限流，这时候是 `unverified` 而不是"伪造"。别把限流当通过，也别当 rejection。

---

## 4. 人工核对原文（10 分钟，不可跳过）

打开原文，确认三件事：

1. 条目里的 endpoint / 方向 / 效应量确实出现在原文（**不是摘要、不是二手转述**）；
2. 剂量参数与单位来自原文，不是你的换算；
3. `notes` 写清楚定位信息（表号 / 图号 / 段名 / 核对日期）。

这一步是项目最值钱的质量闸门。AI 不能代签。

---

## 5. 提交 PR（4 分钟）

把核对过的条目追加到 `data/evidence_matrix.csv`，然后：

```bash
python -m pytest tests -q
python scripts/generate_sbom.py --check
git checkout -b evidence/PBM-013
git add data/evidence_matrix.csv
git commit -m "evidence: add PBM-013 (DOI verified 2026-11-21)"
git push origin evidence/PBM-013
```

PR 描述里请写：来源链接、原文定位、机器裁判输出摘要、为什么这条值得入库（正结果/反证/剂量学缺口）。

---

## 代码贡献

- 改动判定逻辑必须同时更新测试与 `CHANGELOG.md`，并说明证据指纹是否变化。
- 新增模块要有独立验证（用不同实现路径交叉核对），不能只测自己的函数返回值。
- 不要提交凭证、密钥、带 PII 的数据；不要在 `.md` 里复制凭证扫描正则（会自命中）。
- 许可按"代码 / 权重 / 数据"三层核验，动手前读 LICENSE 原文；CC BY-NC / HCP DUA 数据不进 MIT 仓。

---

## 审核流程

维护者按 `GOVERNANCE.md` 的五道门审核。合并后更新 `CITATION.cff` 版本与证据指纹。你的名字会进 CHANGELOG 致谢（如希望匿名请说明）。
