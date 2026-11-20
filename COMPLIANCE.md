# 合规自查：EU AI Act 义务映射

版本：1.0　|　日期：2026-11-20　|　性质：**工程自查文档，不构成法律意见**

本文把 interface-atlas 已有的技术设计，逐条映射到欧盟 AI Act 中**当前已生效**的义务。目的有二：一是让企业采购审查有第一份可出示的材料；二是明确哪些不是我们能满足的、哪些需要法务判断。

法规时间线依据：`人体编程_生产级提升与前沿提升建议报告_2026-11-20.md` §2.3 的公开来源调研。**正式合规结论须由法务复核，本文不替代。**

---

## 1. 适用时间线（据公开调研，需法务确认）

| 义务 | 生效/适用日 | 本项目相关性 |
|---|---|---|
| 第 5 条禁止性规定 | 2025-02-02 起 | 项目伦理红线与之方向一致：不做健康个体增强声明、不帮助规避医疗监管 |
| 第 4 条 AI 素养义务 | 2025-02-02 起（Omnibus 将措辞从"确保"弱化为"采取措施支持"，义务本身保留） | **核心相关**：工具必须让使用者知道"何时不该信" |
| 第 50 条透明度义务与罚则 | 2026-08-02 起 | **核心相关**：输出必须可识别为 AI 生成、带出处、可追溯 |
| Annex III 高风险义务 | 推迟至 2027-12-02 | 本项目当前不进入该分类（非就业/教育/信贷/执法等场景） |
| Annex I 医疗器械相关 | 推迟至 2028-08-02 | 本项目明确声明非医疗器械、非临床决策工具 |

罚则上限：€3500 万或全球营业额 7%（适用于相关违法行为）。这条数字的意义不是吓人，而是说明**大客户会认真审**——所以合规资产要提前备好。

---

## 2. 技术映射：设计 → 义务

| 义务要求 | 本项目的实现方式 | 可出示的证据 |
|---|---|---|
| **Art.50**：AI 生成内容可识别、可追溯 | 每个服务响应强制带 `meta`：`library` / `version` / `evidence_fingerprint` / `evidence_count` / `boundary`。判定结果不是裸文本，永远带来源指纹 | `interface_atlas/service.py: meta(), _wrap()`；`tests/test_service.py` 断言每次响应都带 meta |
| **Art.50**：不冒充人类/权威结论 | `BOUNDARY` 常量明确写"科研降阶模型与社区草案规范，非临床处方、非监管结论、不构成医疗建议"；三入口输出一致 | `service.py: BOUNDARY`；`SECURITY.md` §2 三入口表 |
| **Art.50**：引用可查 | 每条判定带命中证据 ID（`matched_ids` / `counter_ids`），103 条证据均有 DOI 且经 CrossRef/NCBI 核验；证据库指纹 `99b2fdd535ca` 可用于版本比对 | `evidence.py`；`data/evidence_matrix.csv`；`CITATION.cff` 附 HPA 数据源引用 |
| **Art.4**：使用者理解能力边界 | 每条判定附**推翻条件**——即"什么情况下不该信这个结论"；这是把"素养"工程化的直接手段，比写免责声明强 | `claimcheck.py` 四级判定输出；`examples/agent_demo.md`；`L3_BENCHMARK.md` 三层因果断裂 |
| **Art.4**：Agent 使用者同样被覆盖 | MCP 10 个工具的 description 均写明适用域与"不可用于临床决策 / 非监管结论 / 非真实患者影像"；输出 `structuredContent` 同源带指纹 | `mcp_server.py: _TOOLS`；`MCP_SECURITY_REVIEW.md` §2 逐工具能力声明 |
| **Art.5**：不做禁止性用途 | 伦理红线：不得用于增强健康个体的声明、不得规避医疗监管。参数越界只报 `OUT_OF_WINDOW`，不提供"绕过阈值"的调参建议 | 项目战略文档边界声明；`paramcheck.py` 硬阈值单一真源 |
| 高风险分类规避 | 项目定位为科研证据审计与降阶剂量代理，不输出诊断/处方，不接入治疗决策回路 | `service.py: BOUNDARY`；`l3bench.py` 拒绝"模拟剂量→疗效"推断 |

---

## 3. 可出示的证据包

企业审查或客户尽调时，一次给出以下材料，不需要临时整理：

```bash
# 1. 版本与证据指纹（证明结论可追溯到具体库状态）
python -m interface_atlas.cli --version
python -c "from interface_atlas import service; print(service.meta())"

# 2. 一次判定样例：输出自带边界、证据 ID、推翻条件
python -m interface_atlas.cli claim "红光面板能抗衰老" --as-json

# 3. 供应链清单（CycloneDX 1.5）
python scripts/generate_sbom.py --check && cat sbom.json

# 4. 安全政策与威胁模型
cat SECURITY.md MCP_SECURITY_REVIEW.md

# 5. 测试与凭证门禁证据
python -m pytest tests -q
```

指标可用性：CI 在 Ubuntu/Windows × Python 3.11–3.14 八组矩阵跑全部测试；延迟与确定性基线见 `benchmarks/baseline.json`。

---

## 4. 诚实边界

- 本文**不是**合规认证，也不声称"符合 EU AI Act"。它只说明哪些设计要求可以对应到哪条义务。
- 法规适用分类（尤其是否触及 Annex I/III）取决于实际部署方式与 jurisdiction，须由法务判断。项目当前**没有**医疗器械定位。
- 第 50 条针对"AI 生成/操纵内容"的可识别性要求是通用义务；本项目输出的是**带证据的判定**而非生成式内容，映射方式是"让来源和边界可见"，不是"添加 AI 水印"。
- 时间线来源为本轮项目调研报告，非官方法规文本；正式用途请核对 EUR-Lex 原文。
- 罚则数字仅用于说明审查严肃性，不构成对本项目风险敞口的量化评估。
