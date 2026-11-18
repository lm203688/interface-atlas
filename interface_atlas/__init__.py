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

MIT License.
"""

__version__ = "0.5.0"

from . import dosimetry, light_field, mc_light, paramcheck  # noqa: F401
