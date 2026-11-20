# MCP 工具安全过检清单

日期：2026-11-20　|　对象：`interface_atlas/mcp_server.py`（10 个工具，stdio 传输）　|　依据：NSA MCP 安全指南（CSI_MCP_SECURITY.PDF, 2026-05）与 OWASP MCP Top 10 的检查维度

**性质说明**：这是工程自查，不是第三方安全认证。每条状态都给出代码证据，未做的一律标"未做"，不用"已考虑"冒充"已实现"。

---

## 1. 检查单

| 检查维度 | 状态 | 证据 / 说明 |
|---|---|---|
| 传输面最小化 | ✅ 已做 | 仅 stdio，newline-delimited JSON-RPC；无网络监听、无 HTTP 出站 |
| 工具只读 | ✅ 已做 | 10 个工具全部经 `service.py` 读数据后计算；无写文件、无 shell、无环境变量读取 |
| 能力声明明确 | ✅ 已做 | 每个工具 description 写明适用域与"不可用于临床决策 / 非监管结论 / 非真实患者"等边界 |
| 输入键白名单 | ✅ v1.1.0 起 | 10 个 `inputSchema` 全部 `additionalProperties: false`；服务端 `_validate_args()` 真正拒绝未知键，不是只写在文档里 |
| 资源耗尽防护 | ✅ v1.1.0 起 | `service.py` 统一夹取：`n ≤ 1000`、`limit ≤ 500`、文本 ≤ 4000 字符；三入口共享 |
| 错误不泄露栈 | ⚠️ 部分 | 参数错误/未知工具返回精简文案；未预期异常仍向 **stderr** 打印 traceback（不进响应，但 stderr 可能被客户端采集） |
| 工具名枚举固定 | ✅ 已做 | `_DISPATCH` 白名单；未知工具返回 `isError: true`，不做模糊匹配 |
| 请求大小上限（整行） | ❌ 未做 | 单行 JSON 无字节上限；超大合法文本由 `MAX_TEXT_CHARS` 挡住，但畸形超长 JSON 仍会被解析 |
| 速率限制 | ❌ 未做 | 本地 stdio 工具，无多租户场景；如被代理高频调用，由宿主资源限制兜底 |
| 会话级审计日志 | ❌ 未做 | 响应自带证据指纹可追溯结论，但不记录调用历史（有意不写盘，避免留存用户文本） |
| 工具输出 schema 声明 | ❌ 未做 | 目前返回 `content` + `structuredContent`，未在 MCP 层声明 `outputSchema` |
| 影子服务器防护（OWASP MCP09） | ⚠️ 依赖部署方 | 服务器身份为 `interface-atlas`、版本来自内核；客户端应固定 `python -m interface_atlas.mcp_server` 启动命令，不接入未知同名服务器 |

**有意的设计取舍**：本地 stdio 工具没有也不需要有 OAuth / API key / 多租户隔离。加上这些会把"零配置本地调用"的优势抵消掉。真正的风险不是鉴权，而是**代理被提示注入后把工具当武器**——因此防护重点放在"只读 + 输入边界 + 输出自带边界声明"三处。

---

## 2. 逐工具能力声明

| 工具 | 能做到 | 明确不做 | 输入边界 |
|---|---|---|---|
| `claim_check` | 四级宣称判定 + 证据 ID + 推翻条件 | 不做临床决策、不查库外证据 | `text` ≤ 4000 字符，非空 |
| `paramcheck` | 五级参数预检 + 公开硬阈值越界告警 | 非监管结论 | `payload` 须含 modality |
| `dose_audit` | T1–T4 完整性 + 可复现性四档 | 不替报告背书 | `payload` 须含 modality |
| `dose_field` | 三模态降阶剂量场 + 逆设计 | 非临床处方、非数值场真解 | modality 枚举三选一 |
| `evidence_query` | 103 条证据多维过滤 | 不生成新证据、不联网补库 | `limit` 夹到 1–500 |
| `headmodel_list` | 内置体模清单 + 许可门禁结论 | 不提供真实患者影像 | 无参数 |
| `headmodel_sample` | 固定 seed 群体抽样统计 | 非经验分布、非真实患者 | `n` 夹到 1–1000 |
| `headmodel_license` | fail-closed 商用许可判定 | 不替代法律顾问 | `name` 必填 |
| `l3_cases` | 4 个基准案例 + 推翻条件 | 不聚合总分、不宣称疗效 | 无参数 |
| `atlas_report` | 宣称→证据一次跑通 | 非临床决策 | `text` ≤ 4000 字符 |

---

## 3. 可证伪测试

`tests/test_mcp.py` 与 `tests/test_production_hardening.py` 覆盖：

- 未知工具、未知参数、空文本 → `isError: true`，不崩服务；
- 10 个 schema 全部 `additionalProperties: false`（逐工具断言）；
- `n=10⁹` → 实际抽样 1000；`limit=10⁹` → 不超过 500；5000 字符文本 → `ValueError`；
- 同一宣称经 MCP 与经内核调用，判定与证据指纹逐字节一致（零分叉）；
- 完整 stdio 会话：initialize → notification → tools/call，响应 id 顺序正确。

这些是"能推翻本清单结论"的测试：任何一条退化都会让对应的"已做"变成假话。
