"""interface-atlas MCP 服务：把判定内核以 MCP 工具形式暴露给 AI 代理.

零依赖（仅标准库），stdio 传输，newline-delimited JSON-RPC 2.0。
实现 MCP 2025-11-25 协议的 initialize / ping / tools/list / tools/call 与
notifications/initialized（静默处理）。

设计原则（继承项目质量门禁）：
  1. 与 CLI / Web API 完全同源 —— 全部调用 service.py，判定内核零分叉；
  2. 每个工具输出自带证据指纹、库版本与边界声明，代理拿到的结论可审计；
  3. 工具描述写明适用域与不可做的事（不输出疗效承诺、不越出证据库）。

启动：
    python -m interface_atlas.mcp_server
或在 MCP 客户端配置中指向该命令。
"""
from __future__ import annotations

import json
import sys
import traceback

from . import service

PROTOCOL_VERSION = "2025-11-25"

_TOOLS: list[dict] = [
    {
        "name": "claim_check",
        "description": (
            "判定一句物理干预/非药物疗法宣称是否有证据支持。返回四级判定 "
            "SUPPORTED（A/B 级人体证据支持）/ REFUTED（有人体反证）/ EXCEEDED（宣称超出实测终点）"
            "/ NO_EVIDENCE（库中无对应），每条带命中证据 ID、推翻条件与证据库指纹。"
            "适用：科研选题、营销文案审计、合规初筛。不可用于临床决策。"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {"text": {"type": "string", "description": "一句宣称，如『红光面板能抗衰老』"}},
            "required": ["text"],
        },
    },
    {
        "name": "paramcheck",
        "description": (
            "物理干预参数报告预检（五级判决）：参数齐 + 剂量带单位 = COMPLETE；"
            "越过公开硬安全阈值（如 tFUS MI>1.9、温升>2°C）= OUT_OF_WINDOW。"
            "本地自查线索，非监管结论。"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "payload": {"type": "object", "description": "含 modality 与参数键值，如 {'modality':'photobiomodulation','wavelength':'660nm','fluence':'6J/cm2'}"},
                "source_id": {"type": "string", "description": "可选，报告标识"},
            },
            "required": ["payload"],
        },
    },
    {
        "name": "dose_audit",
        "description": (
            "按《物理干预剂量报告规范 v1》（社区草案）审计一份剂量报告的完整性："
            "四层剂量本体 T1 设备/T2 界面/T3 靶区/T4 解剖，输出 completeness_score 与"
            "可复现性四档（not_reproducible / underdetermined / reproducible_with_model /"
            " reproducible_from_report）。"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "payload": {"type": "object", "description": "含 modality 与已报字段键值"},
                "source_id": {"type": "string"},
            },
            "required": ["payload"],
        },
    },
    {
        "name": "dose_field",
        "description": (
            "L2 剂量场引擎（科研降阶代理）：给定设备/组织参数计算靶区真实剂量。"
            "light=光 P1 扩散；ultrasound=超声聚焦+颅骨衰减；electric=四层球谐电场。"
            "设备参数≠靶区剂量：同一外部参数在不同解剖下靶区可差数倍。"
        ),
        "inputSchema": {
            "type": "object",
            "properties": {
                "modality": {"type": "string", "enum": ["light", "ultrasound", "electric"]},
                "params": {"type": "object", "description": "如 {'preset':'muscle_nir','target_depth_cm':1.0}"},
            },
            "required": ["modality", "params"],
        },
    },
    {
        "name": "evidence_query",
        "description": "查询物理干预证据库（103 条，A/B/C/D 分级，来源经 CrossRef/NCBI 核验）。支持模态/基因/细胞/等级/人体证据/关键词过滤。",
        "inputSchema": {
            "type": "object",
            "properties": {
                "modality": {"type": "string"}, "gene": {"type": "string"},
                "cell_type": {"type": "string"}, "grade_min": {"type": "string", "enum": ["A", "B", "C", "D"]},
                "human_only": {"type": "boolean"}, "keyword": {"type": "string"},
                "limit": {"type": "integer", "default": 50},
            },
        },
    },
    {
        "name": "headmodel_list",
        "description": "列出内置解剖代理体模及许可门禁结论（fail-closed：CC BY-NC/DUA 一律拒绝商用）。",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "headmodel_sample",
        "description": "群体解剖抽样（固定 seed 可复现；illustrative 示例分布，非经验值），为个体化剂量区间供数。",
        "inputSchema": {
            "type": "object",
            "properties": {
                "n": {"type": "integer", "default": 100},
                "seed": {"type": "integer", "default": 20261119},
            },
        },
    },
    {
        "name": "headmodel_license",
        "description": "按 fail-closed 许可门禁校验指定体模是否允许（商）用。",
        "inputSchema": {
            "type": "object",
            "properties": {
                "name": {"type": "string"}, "commercial": {"type": "boolean", "default": False},
            },
            "required": ["name"],
        },
    },
    {
        "name": "l3_cases",
        "description": "L3 剂量-效应种子基准案例（4 个一手案例，每条带可证伪推翻条件；多指标协议，拒绝聚合总分）。",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "atlas_report",
        "description": "完整评估旅程：一句宣称 → 四级判定 + 相关证据（一次调用跑通『宣称→证据』链）。",
        "inputSchema": {
            "type": "object",
            "properties": {"text": {"type": "string"}},
            "required": ["text"],
        },
    },
]

_DISPATCH = {
    "claim_check": lambda a: service.claim_check(a["text"]),
    "paramcheck": lambda a: service.param_check(a["payload"], a.get("source_id", "")),
    "dose_audit": lambda a: service.dose_audit(a["payload"], a.get("source_id", "")),
    "dose_field": lambda a: service.dose_field(a["modality"], a.get("params", {})),
    "evidence_query": lambda a: service.evidence_query(
        modality=a.get("modality"), gene=a.get("gene"), cell_type=a.get("cell_type"),
        grade_min=a.get("grade_min"), human_only=bool(a.get("human_only", False)),
        keyword=a.get("keyword"), limit=int(a.get("limit", 50))),
    "headmodel_list": lambda a: service.headmodel_list(),
    "headmodel_sample": lambda a: service.headmodel_sample(
        n=int(a.get("n", 100)), seed=int(a.get("seed", 20261119))),
    "headmodel_license": lambda a: service.headmodel_license(
        a["name"], commercial=bool(a.get("commercial", False))),
    "l3_cases": lambda a: service.l3_cases(),
    "atlas_report": lambda a: service.report(a["text"]),
}


def _result(payload: dict) -> dict:
    return {
        "content": [{"type": "text", "text": json.dumps(payload, ensure_ascii=False, indent=2)}],
        "structuredContent": payload,
        "isError": False,
    }


def _error(code: int, message: str, req_id) -> dict:
    return {"jsonrpc": "2.0", "id": req_id, "error": {"code": code, "message": message}}


def handle_message(msg: dict) -> dict | None:
    """处理一条 JSON-RPC 消息；notification 返回 None."""
    method = msg.get("method", "")
    req_id = msg.get("id")
    if method.startswith("notifications/"):
        return None
    if method == "initialize":
        return {"jsonrpc": "2.0", "id": req_id, "result": {
            "protocolVersion": PROTOCOL_VERSION,
            "capabilities": {"tools": {"listChanged": False}},
            "serverInfo": {"name": "interface-atlas", "version": service.meta()["version"]},
        }}
    if method == "ping":
        return {"jsonrpc": "2.0", "id": req_id, "result": {}}
    if method == "tools/list":
        return {"jsonrpc": "2.0", "id": req_id, "result": {"tools": _TOOLS}}
    if method == "tools/call":
        params = msg.get("params", {})
        name = params.get("name", "")
        fn = _DISPATCH.get(name)
        if fn is None:
            return {"jsonrpc": "2.0", "id": req_id, "result": {
                "content": [{"type": "text",
                             "text": f"未知工具 {name}；可用：{', '.join(_DISPATCH)}"}],
                "isError": True,
            }}
        try:
            return {"jsonrpc": "2.0", "id": req_id, "result": _result(fn(params.get("arguments", {})))}
        except (ValueError, KeyError) as exc:
            return {"jsonrpc": "2.0", "id": req_id, "result": {
                "content": [{"type": "text", "text": f"参数错误：{exc}"}],
                "isError": True,
            }}
        except Exception as exc:  # 服务不得因单次调用崩溃
            traceback.print_exc(file=sys.stderr)
            return {"jsonrpc": "2.0", "id": req_id, "result": {
                "content": [{"type": "text", "text": f"内部错误：{type(exc).__name__}: {exc}"}],
                "isError": True,
            }}
    if req_id is None:
        return None
    return _error(-32601, f"未知方法 {method}", req_id)


def serve_stdio(in_stream=None, out_stream=None) -> int:
    """stdio 主循环：逐行读写 newline-delimited JSON."""
    in_stream = in_stream or sys.stdin
    out_stream = out_stream or sys.stdout
    for line in in_stream:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            out_stream.write(json.dumps(_error(-32700, "Parse error", None)) + "\n")
            out_stream.flush()
            continue
        resp = handle_message(msg)
        if resp is not None:
            out_stream.write(json.dumps(resp, ensure_ascii=False) + "\n")
            out_stream.flush()
    return 0


def main() -> int:
    return serve_stdio()


if __name__ == "__main__":
    raise SystemExit(main())
