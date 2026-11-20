# 安全政策与威胁模型

适用版本：interface-atlas ≥ 1.1.0　|　最近更新：2026-11-20

本项目是**本地科研工具**，不是面向公网的多租户服务。本文如实声明攻击面、已做缓解与尚未达标项，供企业安全团队评估，而不是宣称"绝对安全"。

---

## 1. 漏洞披露

请用 GitHub 私有漏洞报告入口提交，不要在公开 Issue 里贴可利用细节：

> https://github.com/lm203688/interface-atlas/security/advisories/new

处理流程：收到后 5 个工作日内确认；确认后 30 天内给出修复或书面说明；修复随下一个补丁版本发布并在 CHANGELOG 记录。报告者署名默认公开，可要求匿名。

---

## 2. 威胁模型：三个入口，一个内核

| 入口 | 传输/边界 | 能做什么 | 不能做什么 | 主要风险 |
|---|---|---|---|---|
| CLI | 本机进程，用户自有参数 | 全部判定与剂量场计算 | 无网络出站 | 恶意参数（超大输入）导致本地资源占用 |
| MCP Server | stdio，newline-delimited JSON-RPC 2.0，由 AI 代理在本机拉起 | 10 个只读工具，全部经 `service.py` 调内核 | 无文件写入、无网络出站、无 shell 执行、无环境变量读取 | 代理被提示注入后传入超大 `n` / `limit` / 长文本 |
| Web API | `ThreadingHTTPServer`，**默认仅绑定 127.0.0.1:7100** | 10 个只读端点（6 POST / 4 GET） | 无写操作端点、无认证态、不持久化用户数据 | `--host 0.0.0.0` 时无认证、无 TLS，等同把只读工具裸奔到局域网 |

**共同缓解**：所有外部输入的长度/数量边界统一在 `service.py` 实施，CLI / MCP / Web 三入口共享，避免"修了一个门、另两个还开着"。v1.1.0 起：

- `claim_check` / `report` 文本 ≤ 4000 字符；
- `evidence_query` 的 `limit` 被夹到 1–500；
- `headmodel_sample` 的 `n` 被夹到 1–1000；
- 所有 MCP `inputSchema` 标注 `additionalProperties: false`，拒绝未知键。

**边界说明**：这些是资源耗尽防护，不是"输入语义校验"。语义层面的错误（如不存在的 preset）仍由内核 fail-fast 抛 `ValueError`，MCP 转成 `isError: true`，Web 转成结构化错误响应——不静默接受、不返回猜测值。

---

## 3. 数据分类

| 数据 | 来源 | 许可 | 是否含个人信息 |
|---|---|---|---|
| 证据库 103 条 | 公开文献，DOI 经 CrossRef/NCBI 核验 | 条目为事实性摘录；原论文版权归各自出版方 | 否 |
| HPA 单细胞表达 | Human Protein Atlas v25.1 | CC BY 4.0（须署名，见 CITATION.cff） | 否（聚合表达数据） |
| 内置头模 | 解析/统计代理体模，明确标 `illustrative` | 本项目 MIT | 否，非真实患者 |
| 用户请求体 | 运行时传入，不落盘、不写日志 | — | 取决于用户输入；服务端不做留存 |

许可门禁为 fail-closed：标为 CC BY-NC 或受 HCP DUA 约束的外部体模，商用请求一律拒绝，且这类数据不打包进 MIT 仓库。

---

## 4. 依赖与供应链

**必须澄清一处常见误述**：项目**判定内核与三个入口**是零第三方依赖（仅 Python 标准库）；但**完整 Python 包不是零依赖**。`pyproject.toml` 声明 3 个直接运行时依赖：

| 依赖 | 声明约束 | 实际用途 | 是否在主判定链路上 |
|---|---|---|---|
| numpy | ≥1.24 | 光 MC 金标准、电场球谐求解、剂量场数值运算 | 是（L2 数值引擎） |
| pandas | ≥2.0 | 受体表达 CSV 读取 | 否（受体面板） |
| matplotlib | ≥3.7 | 受体热图生成 | 否（可选绘图） |

可选延迟导入：`nibabel`（NIfTI 体模导入，未安装时该功能不可用，不影响其余能力）。

机器可读清单见 `sbom.json`（CycloneDX 1.5，由 `scripts/generate_sbom.py` 生成）。**未达标项，如实列出**：

- 直接依赖只有版本下界，未锁定传递依赖，无 lockfile；
- release tag 尚未签名（需维护者本人的 GPG/sigstore 密钥；流程见 `RELEASING.md`）；
- 依赖漏洞自动扫描尚未接入（当前 CI 只跑测试与凭证扫描）。

---

## 5. 已落实的防护

- CI（GitHub Actions）在 Ubuntu/Windows × Python 3.11–3.14 八组矩阵跑全部测试，并在每组做凭证入库扫描；
- 三入口输出统一携带证据指纹、库版本与边界声明，代理拿到的结论可追溯；
- MCP 对未知工具、参数错误、内部异常分别返回结构化 `isError`，单次调用不会杀死服务进程；
- 无 `subprocess` / `os.system` / `eval` / `exec` 调用，无运行时文件写入（唯一写函数 `evidence.save_seed()` 是开发期重建种子数据的工具，不在判定链路）。

---

## 6. 部署建议

1. 默认本机使用即可，不要把 Web 服务暴露到公网。
2. 确需局域网访问时，用 `python -m interface_atlas.webapp --host <addr>` 之前，先在外层加认证代理与 TLS；本工具自身不提供这两项。
3. 把 MCP 客户端的 cwd 指向只读工作区；stdio 服务器继承启动它的用户权限，不要以高权限账户拉起。
4. 升级前核对 CHANGELOG；判定内核变更会同时改变证据指纹，可用于确认版本一致性。
