# MCP 自然语言驱动示例：一次完整的评估旅程

> 前提：MCP 客户端（Claude Desktop / Cursor / Kimi Work 等）中配置命令
> `python -m interface_atlas.mcp_server`（在仓库根目录、已安装依赖的环境中）。
> 以下对话展示 AI 代理如何用 10 个工具跑通「宣称 → 证据 → 剂量完整性 → 靶区剂量 → 基准口径」五步旅程。
> 所有结论来自 interface-atlas 判定内核，每条带证据 ID、推翻条件与证据库指纹（`99b2fdd535ca`，v1.0.0）。

---

**用户**：帮我审一下这句营销宣称——"红光面板照一照就能抗衰老、延缓衰老"。

**代理**（调用 `claim_check`）：

```
verdict: REFUTED
reason: 存在 1 条 A/B 级人体『无效应』证据：PBM-009(A)（Son 2025 伞形综述：
        15 个 meta / >9000 人 / 35 终点，0 个高确定性，23 个终点无显著效应）
counter_ids: [PBM-009]
falsification: 需要针对该人群/剂量窗的新 RCT 推翻现有 meta 结论；
               注意反证针对群体平均效应，不排除亚组获益
```

**用户**：那他们的剂量参数写得怎么样？能复现吗？

**代理**（调用 `dose_audit`，以典型红光面板说明书字段为例）：

```
payload: {modality: "photobiomodulation", wavelength: "660nm",
          power_density: "25mW/cm2", fluence: "6J/cm2", exposure_time: "240s"}
completeness_score: 0.35
reproducibility: not_reproducible
missing_params: [applicator_distance, applicator_angle, beam_profile,
                 target_fluence, target_depth, tissue_optics, anatomy_source]
reason: T2 界面条件全缺 —— 施加器距离/角度决定到达皮肤的真实功率密度，
        不报就无法跨设备比较（这正是 PBM 领域 35 终点零高确定性的剂量学根因）
```

**用户**：如果补上界面条件，真实到达深层组织的剂量怎么算？

**代理**（调用 `dose_field`，light 模态）：

```
params: {preset: "muscle_nir", target_depth_cm: 1.0}
结果要点：表面归一相对光通量场；μeff 决定深度衰减；
          适用域外（强吸收/近源）引擎会主动 WARN，不静默给不可信数字。
注意：这是相对剂量场；绝对剂量需表面实测 + 独立 MC 标定。
```

**用户**：颅骨对聚焦超声的影响呢？厂商说"穿颅无忧"。

**代理**（调用 `dose_field`，ultrasound 模态）：

```
params: {preset: "tfus_transcranial", skull: "skull_thin"}
结果要点：焦点峰值声压随颅骨厚度 2mm→8mm 变化约 2.37× ——
          设备报的同一声强，靶区剂量差两倍以上。
          MI≤1.9 与温升≤2°C 为硬阈值，越界即 OUT_OF_WINDOW，与疗效叙事无关。
```

**用户**：这些结论在学术上处于什么位置？能进临床试验设计吗？

**代理**（调用 `l3_cases`）：

```
B1 THIELSCHER-2026：tES 统一场强框架（0.2 V/m，peer-reviewed）
B2 KOP-2026：MR-ARFI 实测位移预测效应，模拟强度不预测（preprint）
B3 SPIVAK-2026：tFUS 剂量升至 14× FDA 限仍 inconclusive（peer-reviewed）
B4 GHAEMPANAH-2025：分割管线使聚焦体积 -39% 而峰值不变（conference）
```

**用户**：这套工具本身的边界是什么？

**代理**（调用 `meta` / 各响应自带）：

```
版本 v1.0.0 · 证据库 103 条（A15/B42/C46）· 指纹 99b2fdd535ca
边界：科研降阶模型与社区草案规范，非临床处方、非监管结论、不构成医疗建议；
      自发报告≠发生率；RNA≠膜蛋白；设备参数≠靶区剂量。
```

---

## 工具清单

| 工具 | 一句话 |
|---|---|
| `claim_check` | 一句宣称 → 四级判定 + 推翻条件 |
| `paramcheck` | 参数报告 → 五级预检（越硬阈值即 OUT_OF_WINDOW） |
| `dose_audit` | 剂量报告 → T1–T4 完整性 + 可复现性四档 |
| `dose_field` | 设备/组织参数 → 靶区真实剂量（light/ultrasound/electric） |
| `evidence_query` | 103 条证据库查询（模态/基因/细胞/等级/关键词） |
| `headmodel_list` / `headmodel_sample` / `headmodel_license` | 个体化解剖层：体模清单 / 群体抽样 / fail-closed 许可门禁 |
| `l3_cases` | L3 种子基准（4 案例，每条带推翻条件） |
| `atlas_report` | 一次调用跑通「宣称→证据」链 |
