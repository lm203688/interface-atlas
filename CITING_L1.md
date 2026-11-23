# 如何引用 L1 剂量报告规范

> 面向外部研究者、系统综述作者、规范制定者。目标：让"设备参数 → 靶区剂量"的报告口径有稳定、可核对、可版本冻结的引用锚点。
> 规范正文：`DOSE_REPORTING_SPEC.md`　|　机器可读条文：`schema/l1_clauses_v1.json`　|　唯一权威字段源：`interface_atlas/dosimetry.py::DOSE_SPEC`

## 1. 为什么需要一份引用指南

规范若只能被笼统引用（"按 interface-atlas 规范报告波长"），就无法精确到条、无法版本冻结、无法机器核对，也就无法积累"领域依赖度"。本指南把引用分成两级：**整份规范**与**单条条文**。每条条文都有稳定 ID，ID 以语义字段名为锚——插入新字段不会让已发布 ID 漂移。

## 2. 引用整份规范

```text
interface-atlas contributors. Physical Intervention Dose Reporting Spec v1
(interface-atlas L1). Spec version 1.0.0, package v1.4.0.
https://github.com/lm203688/interface-atlas
```

- 规范版本：**v1.0.0**（2026-11-18 制定；条文注册表 v1，2026-11-21 编译）
- 代码版本锚点：见下方"版本冻结"
- BibTeX：

```bibtex
@misc{interfaceatlas2026l1,
  title  = {Physical Intervention Dose Reporting Spec v1 (interface-atlas L1)},
  author = {{interface-atlas contributors}},
  year   = {2026},
  version = {1.0.0},
  url    = {https://github.com/lm203688/interface-atlas},
  note   = {Community draft; not an ICRU/FDA/ISO standard}
}
```

## 3. 引用单条条文

条文 ID 格式：`IA-L1-<模态短码>-<层>-<字段key>`

```text
... 施加器-皮肤距离按 interface-atlas L1 条文 IA-L1-PBM-T2-applicator_distance
(规范 v1.0.0 / 包 v1.4.0) 报告。
```

机器读取：

```bash
python -m interface_atlas.cli l1clauses --cite IA-L1-PBM-T2-applicator_distance
```

模态短码（一经发布不得更改）：

| 短码 | 模态 | 条文数 | 备注 |
|---|---|---|---|
| PBM | photobiomodulation | 12 | 光生物调节 |
| VNS | vagus nerve stimulation | 11 | 迷走神经刺激 |
| FUS | focused ultrasound | 12 | 聚焦超声 |
| MECH | mechanical loading | 9 | 机械载荷 |
| MET | metabolic / microbiota | 4 | 附录，不纳入主评分 |

层：`T1` 设备参数 / `T2` 界面条件 / `T3` 靶区剂量 / `T4` 个体化解剖。层权重 T1=0.40、T2=0.25、T3=0.25、T4=0.10。

判定规则也可逐条引用：`IA-L1-RULE-001`（not_reproducible）至 `IA-L1-RULE-005`（completeness_score）。

可直接引用的示例条文（原文，非 paraphrase）：

- `IA-L1-PBM-T2-applicator_distance`：报告须给出施加器-皮肤距离（cm），作为界面条件（T2）。须描述从设备到靶区的传递条件；该项是跨研究比较的关键，也是文献最常漏报项。
- `IA-L1-FUS-T2-skull_attenuation`：报告须给出颅骨衰减（dB or text），作为界面条件（T2）。须描述从设备到靶区的传递条件；该项是跨研究比较的关键，也是文献最常漏报项。
- `IA-L1-VNS-T3-neural_dose`：报告须给出神经剂量(ECAP归一化)（mV or normalized ECAP），作为靶区剂量（T3）。须为测量值，或明确标注计算来源的估计值；模拟量不得冒充体内实测。

完整 48 条 + 5 条判定规则以注册表为准：

```bash
python -m interface_atlas.cli l1clauses --list        # 人类可读全文
python -m interface_atlas.cli l1clauses --save        # 重建注册表（权威源变更时）
python -m interface_atlas.cli l1clauses --check       # 核对注册表未与权威源脱节
```

## 4. 版本冻结与变更规则

| 变更类型 | 版本动作 | 对已发布 ID 的影响 |
|---|---|---|
| 字段新增 / 删除 | minor（规范 v1.x） | 新增字段获得新 ID；删除字段的 ID 进入 deprecated 名单，**不复用** |
| 别名、单位等价写法增补 | patch | 不变 |
| 层结构、`T` 语义、判定规则、层权重变更 | major（规范 v2） | 必须重新编号，提供 v1→v2 映射表 |
| 模态短码变更 | major | 禁止；短码一经发布永久固定 |

注册表 `spec_version` 为规范版本，`citation.version_pin` 为代码版本锚点。引用时应同时写两者。

版本冻结实现：注册表每次由 `DOSE_SPEC` 确定性编译（键排序 + 固定缩进），`sha256(schema/l1_clauses_v1.json)` 即该版条文的指纹——**引用时同时记录该 sha256，即可精确冻结内容**，不依赖尚未生成的 commit SHA。`--check` 可在任何时刻验证注册库与权威源未脱节。当前指纹：

```text
sha256(l1_clauses_v1.json) = c20b491ebddf4a65134003722897c54ab384e04fb1c3f7fe8fb4a6b8d575d581
```

维护者签名 tag（见 `RELEASING.md`）后，可改以 `v1.4.0` tag 作为额外锚点。指纹来自权威源编译，条文内容或版本锚点任一变化都会改变指纹。

## 5. 实证出处的诚实边界

注册表每条挂 `evidence_refs`，来自本仓 103 条真实证据库（`data/evidence_matrix.csv`），并标注来源类型：`pmc` / `doi` / `institutional`。这些出处是**该模态的实证背景**——说明"这个字段为什么必要"，例如 PEMF/EMTT 综述因欠报线圈几何而无法 meta 分析。

必须明确：

- 出处**不代表**该字段被这些文献验证；规范字段是报告要求，不是疗效结论。
- 出处**不构成**疗效或安全性证据；机构页与 DOI 一视同仁地标记为背景，不冒充金标准。
- 限流（如 CrossRef 429）或离线导致的无法核验，标 `unverified`，不判伪也不冒充通过（见 `GOVERNANCE.md` 与 `harvest.py` 的三态核验）。

## 6. 免责

本规范是社区草案，**非 ICRU / FDA / ISO 官方标准**。参数落在窗口内不代表有效，越界不代表必然有害。`completeness_score` 只度量报告完整性，不度量疗效或安全性。使用者应自行做 IRB 与设备安全评估。
