"""interface-atlas: 非药物接口层的接收器图谱与证据剂量图谱.

子模块：
  receptors  —— A1 接收器能力图谱（HPA 单细胞 RNA）
  evidence   —— A2 物理干预证据与剂量图谱（schema v1.2，103 条）
  claimcheck —— 宣称-证据对齐判定（可证伪、可追溯）
  link       —— A1↔A2 交叉核验
  paramcheck   —— 物理干预参数报告预检（本地自查，五级判决）
  dosimetry    —— 剂量报告规范 v1：四层剂量本体 + 完整性审计（L1 生态入口）
  light_field  —— L2 光剂量场：P1 扩散近似正向解 + 逆设计（可微、CPU 可跑）
  mc_light     —— L2 独立 Monte Carlo 求解器（验证金标准 + 解剖敏感性）
  calibration  —— L2 中阶：split conformal 剂量区间（有限样本覆盖保证）
  ultrasound_field —— L2 超声/tFUS 剂量场：换能器参数 + 颅骨衰减 → 靶区声压/MI/温升
  us_validation —— L2 超声独立验证：RS 积分交叉核对 + 代理误差如实报告
  electric_field —— L2 电/tES-tDCS 剂量场：四球解析解 + 靶区电场 + 逆设计
  es_validation —— L2 电场独立验证：单层球谐 vs 解析闭式 + 可证伪性质
  anatomy      —— L2 个体化解剖：头模数据契约 + 代理体模 + 导入适配器 + 许可门禁
  dose_features —— L3 供给：三模态多指标剂量场特征（拒绝单一标量）
  l3bench      —— L3 剂量-效应种子基准：三层断裂 + 多指标 + 可证伪
  l3_validation —— L2 个体化 + L3 基准的独立验证
  service      —— 产品化服务层：CLI / Web API / MCP 三入口共享内核（v1.0）
  webapp       —— 零依赖本地 Web 服务（REST API + 交互站，仅 127.0.0.1）
  mcp_server   —— MCP 服务：判定内核以 MCP 工具形式暴露给 AI 代理

MIT License.
"""

__version__ = "1.1.0"

from . import calibration, dosimetry, light_field, mc_light, paramcheck  # noqa: F401
from . import ultrasound_field, us_validation  # noqa: F401
from . import electric_field, es_validation  # noqa: F401
from . import anatomy, dose_features, l3bench, l3_validation  # noqa: F401
