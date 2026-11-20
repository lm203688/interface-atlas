# 发布流程（tag → 签名 → 发布说明）

版本策略：语义化版本（MAJOR.MINOR.PATCH）。判定内核或输出 schema 变化 → MINOR；仅文档/工程件 → PATCH；破坏性接口变化 → MAJOR。

---

## 1. 发布前门禁（AI 可自主完成）

```bash
cd interface-atlas
PYTHONPATH=. python -m pytest tests -q          # 全部测试
PYTHONPATH=. python benchmarks/latency.py --check  # SLO：P95 < 200ms + 100% 确定性
python scripts/generate_sbom.py --check         # SBOM 与 pyproject 一致
# 凭证扫描：CI 在每组矩阵自动执行（见 .github/workflows/test.yml 同名步骤）。
# 注意：不要把扫描正则复制进任何 .md —— 扫描器会命中文档里的正则字面量，
# 制造"仓库有凭证"的假阳性（本文件初版就踩过这个坑）。
PYTHONPATH=. python -m interface_atlas.cli --version
```

门禁通过后：更新 `CHANGELOG.md`（新增条目 + 日期）、`pyproject.toml` / `interface_atlas/__init__.py` / `CITATION.cff` 三处版本号，重新生成 `sbom.json` 与 `benchmarks/baseline.json`，通过 GitHub Contents API 推送 main 并逐文件 SHA 复验。

---

## 2. 签名 tag（**需维护者本人执行，AI 不代持密钥**）

本仓库采用 GitHub Contents API 直推，**本地没有 git 仓库对象**，因此不能在本机直接 `git tag -s`。维护者本人按以下任一路径操作：

**路径 A（推荐，本地 clone 后签名）**
```bash
git clone https://github.com/lm203688/interface-atlas.git
cd interface-atlas
git tag -s v1.1.0 -m "v1.1.0: 生产级硬化（安全/合规/SBOM/SLO）"
git push origin v1.1.0
```

**路径 B（GitHub Web）**
在 Releases 页新建 release 并勾选 "Sign with GitHub"（sigstore 签名，无需本地 GPG）。

前置条件：路径 A 需要本地 GPG key 且已配置 `git config user.signingkey`；路径 B 需要仓库写权限与 GitHub 账号。**这两项都是账号/密钥硬阻塞，AI 不代持、不绕过。**

---

## 3. 发布说明

签名 tag 后在 GitHub Releases 写发布说明，至少包含：

1. 本版本相对上一版的**可验证变化**（不是形容词，是命令或测试结果）；
2. 证据指纹（`python -c "from interface_atlas import service; print(service.meta()['evidence_fingerprint'])"`）；
3. 三入口一致性声明：CLI / MCP / Web 共享 `service.py`，判定内核零改动；
4. 已知边界：科研降阶模型、非临床处方、非监管结论；
5. 附件：`sbom.json`、`benchmarks/baseline.json`。

---

## 4. 发布后（P1，需本人账号）

- **Zenodo DOI**：归档该 tag，把 DOI 回填 `CITATION.cff`（使引用从"形式正确"变成"可解析到持久记录"）。
- **PyPI 决策**：`pip install interface-atlas` 当前**未发布**；要么发布，要么在 README 写明刻意不发布及理由——静默缺席才是问题。
- **采用度量**：每月记录一次 GitHub traffic 与外部首次引用，写入 MAINTAINERS 回顾（该文件尚未建立，属 P1 治理件）。

---

## 5. 当前状态

| 项 | 状态 |
|---|---|
| v1.0.0 已推送 main | ✅ `2026-11-19T07:18:59Z`，69 文件 |
| v1.0.0 签名 tag | ❌ 未签名（需本人密钥） |
| v1.1.0 门禁 | 见本轮推送复验结果 |
| Zenodo DOI / PyPI | ❌ 未做（P1） |
