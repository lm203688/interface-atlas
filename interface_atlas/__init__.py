"""interface-atlas: 非药物接口层的接收器图谱与证据剂量图谱.

四个子模块：
  receptors  —— A1 接收器能力图谱（HPA 单细胞 RNA）
  evidence   —— A2 物理干预证据与剂量图谱（schema v1.2，103 条）
  claimcheck —— 宣称-证据对齐判定（可证伪、可追溯）
  link       —— A1↔A2 交叉核验
  paramcheck —— 物理干预参数报告预检（本地自查，分级判定）

MIT License.
"""

__version__ = "0.3.0"

from . import paramcheck  # noqa: F401
