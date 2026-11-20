# CHANGELOG

遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 与语义化版本。

## [1.1.0] - 2026-11-20

按《生产级提升与前沿提升建议报告》P0 六项落地：工程外壳硬化，**判定内核零改动**。

### 新增

- **安全件**：`SECURITY.md`（威胁模型：CLI/MCP/Web 三入口攻击面 + 数据分类 + 供应链现状）；
  `MCP_SECURITY_REVIEW.md`（10 个工具逐项对照 NSA 指南与 OWASP MCP Top 10 检查单，含 4 项如实标注"未做"）。
- **供应链件**：`scripts/generate_sbom.py` + `sbom.json`（CycloneDX 1.5，确定性生成）。
  **纠正报告措辞**：判定内核与三入口零第三方依赖，但完整包声明 numpy/pandas/matplotlib 三个直接依赖，SBOM 如实列出，不谎称整包零依赖。
- **合规件**：`COMPLIANCE.md`（EU AI Act Art.50 透明度 / Art.4 素养 / Art.5 禁止性规定的逐条工程映射，
  附可出示证据包；明确非法律意见、时间线需法务核对）。
- **可靠性件**：`benchmarks/latency.py` + `baseline.json`（实测最差 P95 4.34ms，上限 200ms；
  8 类响应 100% 确定性），CI 新增 SBOM 一致性与 SLO 基准两步。
- **注入边界**：`service.py` 统一夹取 `n ≤ 1000`、`limit ≤ 500`、文本 ≤ 4000 字符，
  CLI/MCP/Web 三入口同时生效；MCP 10 个 schema 全部 `additionalProperties: false`，
  服务端真正拒绝未知键（不只是文档声明）。
- **发布流程**：`RELEASING.md`（tag → 签名 → 发布说明；**签名需本人 GPG/sigstore 密钥，AI 不代持**）。
- **测试 +16 项**：SBOM/依赖漂移、输入边界、schema 闭合、三入口同源、确定性、文档存在性，累计 172 passed。

### 未做（如实）

- 签名 tag、Zenodo DOI、PyPI 发布（前两项需本人账号/密钥；PyPI 属 P1 决策）。
- MCP 整行请求字节上限、速率限制、输出 schema 声明、会话审计日志（见 MCP 过检清单）。

## [1.0.0] - 2026-11-19

产品级发布：内核不变，新增三个对外入口与一套工程信誉件。全部判定逻辑与 v0.9.0 完全同源。

### 新增

- **MCP 服务**（`interface_atlas/mcp_server.py`）：10 个 MCP 工具（claim_check / paramcheck /
  dose_audit / dose_field / evidence_query / headmodel_list / headmodel_sample /
  headmodel_license / l3_cases / atlas_report），stdio JSON-RPC、零依赖，
  供 AI 代理直接调用；每个工具输出自带证据指纹、库版本与边界声明。
- **Web 服务**（`interface_atlas/webapp.py`）：纯标准库 REST API（`/api/claim`、
  `/api/paramcheck`、`/api/doseaudit`、`/api/dosefield`、`/api/evidence`、
  `/api/report`、`/api/headmodels`、`/api/l3cases`、`/api/meta`），仅绑定 127.0.0.1。
- **交互站**（`web/index.html`）：五步评估旅程（宣称→证据→参数→剂量→基准）一屏跑通。
- **共享服务层**（`interface_atlas/service.py`）：CLI / API / MCP 三入口的统一内核门面，
  保证结果逐字节一致。
- **工程信誉件**：GitHub Actions（3.11–3.14 × Ubuntu/Windows 矩阵 + 凭证扫描）、
  CITATION.cff、本 CHANGELOG。
- **测试 27 项**：test_service（11）/ test_mcp（10）/ test_webapp（10，部分重叠计数），
  含三入口同源性、防恒真（判定随输入变化）、fail-closed 许可门禁、
  stdio 端到端会话等可证伪断言。

### 不变

- 判定内核（claimcheck / paramcheck / dosimetry / 三模态剂量场 / L3 基准）零改动。
- 五条质量门禁（可证伪 / 可复现 / 不编造 / 许可闭环 / 边界声明）原样保留。

### 边界

- MCP / Web 层只做搬运与呈现，不构成医疗建议；剂量学输出仍为科研降阶模型。
