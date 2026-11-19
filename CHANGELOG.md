# CHANGELOG

遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/) 与语义化版本。

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
