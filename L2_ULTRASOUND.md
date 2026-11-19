# L2 可计算剂量场引擎——超声/tFUS 模态（alpha，v0.7.0）

> 承接 L1「剂量报告规范 v1」。L1 把 tFUS/经颅超声报告普遍判为 `not_reproducible`——
> 证据库 15 条 tFUS **全部**不报换能器几何、颅骨衰减、靶区声压（平均完整度 0.114，
> `reproducible_*` 0 条）。本模块是「设备参数 + 界面条件 → 靶区声压剂量」这一跳的
> 可计算 model，与光模态 `light_field` 平行，把超声的剂量黑箱显式化。

## 一句话

给定换能器几何（孔半径、焦距、频率、表面声强）+ 颅骨/脑声学参数，算出焦点处
**峰值声压、声强 Isppa、机械指数 MI、温升粗估**，核对公开硬阈值，并能逆设计
（给靶区目标声压反推所需表面声强）。

## 物理模型（线性声学降阶代理）

| 量 | 公式 | 说明 |
|---|---|---|
| 波长 | λ = c / (f·1e4) cm | c[m/s], f[MHz]；1MHz 水 → 0.148cm |
| 孔径半角 | α = asin(a/R) | a 孔半径，R 焦距 |
| 聚焦增益 | G = η·(π·a·sinα / λ)² | 衍射受限焦点声强 / 表面声强 |
| 分层衰减 | α_i(f)=α_i0·(f/f0)^ν；T=10^(−Σαd/20) | dB/20 为声压幅度口径 |
| 焦点声强 | I_focus = I_surf · G · T² | 声压透射 T → 声强是 T² |
| 峰值声压 | p_pk = √(2ρc·I_focus) | 平面波 I=p²/(2ρc) |
| 机械指数 | MI = p_pk[MPa] / √f[MHz] | AIUM/FDA 定义 |
| 温升粗估 | ΔT = μ_vol·I_timeavg·t / (ρ·c_p) | 线性沉积上界，非热剂量学 |

**硬阈值与证据库单一真源**：`MI_LIMIT=1.9`、`TEMP_LIMIT_C=2.0` 直接对齐
`paramcheck.HARD_LIMITS["focused ultrasound"]` 与 `TEMP_LIMIT`（测试锁定同源，
改一处另一处必须跟着变）。越界即报 `exceeded=True` + WARN，不静默给「有效剂量」。

## 实测剂量落点（固定参数，可复现）

| 预设 | 颅骨 | 焦点声压 | Isppa | MI | ΔT 粗估 | 越界 |
|---|---|---:|---:|---:|---:|---|
| tfus_transcranial (0.5MHz, a3cm, F8cm, I_surf 0.5) | 4mm | 0.727 MPa | 16.3 W/cm² | 1.029 | 0.044°C | 否 |
| 同上 | 8mm | 0.382 MPa | 4.5 W/cm² | 0.540 | 0.012°C | 否 |
| tfus_deep (0.25MHz, a4cm, F10cm) | 4mm | 0.642 MPa | 12.7 W/cm² | 1.284 | 0.034°C | 否 |
| 高声强场景 (I_surf 5 W/cm²) | 4mm | — | — | 3.25 | — | **是（MI>1.9 报 WARN）** |

**表面声强语义**：I_surf = 换能器输出声功率 ÷ 有效孔面积。经颅 tFUS 典型
P~几 W / A~28cm² ⇒ I_surf~0.1–1 W/cm²。切勿把文献常报的「焦点声强」塞进表面字段
（开发中实测发现：误填 15 W/cm² 会让焦点 MI=9.35 越界——引擎正确报警，但预设不真实）。

**解剖敏感性**（PNAS 论点超声版）：同换能器/驱动，仅颅骨厚度 2→8mm，
焦点声压按 exp(−α·d/20) 变化，max/min 比 ≈ **2.37×**。这就是「设备参数 ≠
靶区剂量」的定量证据——个体化解剖必须进入剂量学。

## 独立验证（`us_validation.py`，机器可读报告 `data/us_validation_report.json`）

| 核对项 | 方法 | 结果 |
|---|---|---|
| 平面活塞轴上声场 | 数值 RS 求积（numpy trapz 20万点）vs 变量代换闭式 | PASS，RMS=0，max=0 |
| 声强-声压往返 | I=p²/(2ρc) 正逆变换，水/脑两介质 | PASS，误差<1e-12 |
| MI 定义 + 阈值同源 | 频率×4→MI 减半；阈值与 paramcheck 一致 | PASS |
| 颅骨衰减量级 | 0.5MHz 经 4/8mm 颅骨 → 4.4/10.0 dB（透射 0.60/0.32） | PASS，落文献区间 |
| 温升线性 | 对声强/占空比/时长线性 | PASS |
| 聚焦增益 | **代理测量（无独立闭式）** | surrogate，不作 pass/fail |

**关于 RMS=0（务请注意，防止误读）**：轴上场核对用的是 RS 积分
`∫₀ᵃ e^{−ikR}/R · r dr`，变量代换 `r·dr = R·dR` 是**精确恒等式**（非近似），
故解析闭式与 20万点 trapz 数值积分必然逐点吻合到浮点精度。这一项验证的是
**RS 积分的实现写对了**（对应光模块用 MC 复现 Prahl 基准的 Tt），**不是**验证
聚焦换能器代理。聚焦增益 `G=η(πa·sinα/λ)²` 没有独立闭式可严格核对，已如实标为
surrogate、不调参凑「漂亮」目标值——与光模块「P1 代理误差是测量输出」同一纪律。

## 能力边界（不粉饰）

- **线性声学**：不含冲击波、非线性 B/A、空化、辐射力、声流；高 MI 区代理失效并 WARN。
- **均匀衰减层颅骨**：不含层内多次散射、纵→横波模式转换、速度频散、CT 非均匀性；
  真实经颅畸变远复杂，0.5MHz 结果只作降阶演示。f>1MHz 经颅会额外 WARN。
- **温升粗估上界**：无精确热边界条件，不做热剂量学结论。
- **绝对剂量依赖真值**：换能器实测声场、个体颅骨 CT 派生声学参数——这些正是文献普遍
  欠报的量。本引擎把「缺哪个参数导致剂量不确定」显式化，不假装已解决。

## CLI

```bash
# 正向剂量场（预设 + 颅骨层）
python -m interface_atlas.cli usfield --preset tfus_transcranial --skull skull_thin
# 自定义表面声强 / 颅骨厚度
python -m interface_atlas.cli usfield --preset tfus_deep --intensity 1.0 --skull-mm 0.6
# 逆设计：目标靶区声压 → 所需表面声强
python -m interface_atlas.cli usfield --preset tfus_transcranial --inverse 0.5 --skull skull_thin
# 解剖敏感性（颅骨厚度扫描）
python -m interface_atlas.cli usbench --which skull
# 全套独立验证报告
python -m interface_atlas.cli usbench --which all
```

## 测试

```bash
PYTHONPATH=. python -m pytest tests -q     # 81 passed
```

## 免责

科研降阶模型，非临床处方、非治疗计划系统。不构成医疗建议。代码 MIT。
