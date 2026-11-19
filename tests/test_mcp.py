"""mcp_server 的协议测试：initialize / tools/list / tools/call 与错误路径."""
from __future__ import annotations

import json

from interface_atlas import mcp_server, service


def _call(server_fn, method, params=None, req_id=1):
    msg = {"jsonrpc": "2.0", "id": req_id, "method": method}
    if params is not None:
        msg["params"] = params
    return server_fn(msg)


def test_initialize_negotiates_protocol():
    resp = _call(mcp_server.handle_message, "initialize",
                 {"protocolVersion": "2025-11-25", "capabilities": {}, "clientInfo": {"name": "test"}})
    assert resp["result"]["protocolVersion"] == "2025-11-25"
    assert resp["result"]["serverInfo"]["name"] == "interface-atlas"
    assert resp["result"]["serverInfo"]["version"] == service.meta()["version"]


def test_ping():
    assert _call(mcp_server.handle_message, "ping")["result"] == {}


def test_tools_list_covers_journey():
    resp = _call(mcp_server.handle_message, "tools/list")
    names = {t["name"] for t in resp["result"]["tools"]}
    expected = {"claim_check", "paramcheck", "dose_audit", "dose_field",
                "evidence_query", "headmodel_list", "headmodel_sample",
                "headmodel_license", "l3_cases", "atlas_report"}
    assert expected <= names
    for t in resp["result"]["tools"]:
        assert t["inputSchema"]["type"] == "object"
    # 边界声明写进工具描述（至少核心判定工具必须带）
    desc = {t["name"]: t["description"] for t in resp["result"]["tools"]}
    assert "不可" in desc["claim_check"] or "非" in desc["claim_check"]


def test_tools_call_claim_check_same_as_kernel():
    resp = _call(mcp_server.handle_message, "tools/call",
                 {"name": "claim_check",
                  "arguments": {"text": "植入式迷走刺激能治疗类风湿关节炎"}})
    assert resp["result"]["isError"] is False
    payload = resp["result"]["structuredContent"]
    assert payload["result"]["verdict"] == "SUPPORTED"
    assert payload["result"]["evidence_fingerprint"] == service.meta()["evidence_fingerprint"]
    # 文本通道与结构化通道一致
    assert json.loads(resp["result"]["content"][0]["text"]) == payload


def test_tools_call_verdict_varies_with_input():
    """防恒真：同一工具对不同输入给出不同判定."""
    v1 = _call(mcp_server.handle_message, "tools/call",
               {"name": "claim_check", "arguments": {"text": "红光面板能抗衰老"}}
               )["result"]["structuredContent"]["result"]["verdict"]
    v2 = _call(mcp_server.handle_message, "tools/call",
               {"name": "claim_check", "arguments": {"text": "植入式迷走刺激能治疗类风湿关节炎"}}
               )["result"]["structuredContent"]["result"]["verdict"]
    assert v1 != v2


def test_tools_call_unknown_tool_is_error():
    resp = _call(mcp_server.handle_message, "tools/call",
                 {"name": "make_claim_true", "arguments": {}})
    assert resp["result"]["isError"] is True
    assert "未知工具" in resp["result"]["content"][0]["text"]


def test_tools_call_bad_args_returns_iserror_not_crash():
    resp = _call(mcp_server.handle_message, "tools/call",
                 {"name": "dose_field", "arguments": {"modality": "light", "params": {"preset": "nope"}}})
    assert resp["result"]["isError"] is True
    resp2 = _call(mcp_server.handle_message, "tools/call",
                  {"name": "claim_check", "arguments": {"text": "  "}})
    assert resp2["result"]["isError"] is True


def test_notifications_are_silent():
    assert mcp_server.handle_message({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None


def test_unknown_method_error():
    resp = _call(mcp_server.handle_message, "resources/list")
    assert resp["error"]["code"] == -32601


def test_stdio_session_end_to_end():
    """模拟一次完整 stdio 会话：initialize → tools/list → tools/call."""
    lines = [
        json.dumps({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                    "params": {"protocolVersion": "2025-11-25", "capabilities": {}}}),
        json.dumps({"jsonrpc": "2.0", "method": "notifications/initialized"}),
        json.dumps({"jsonrpc": "2.0", "id": 2, "method": "tools/call",
                    "params": {"name": "claim_check", "arguments": {"text": "红光面板能抗衰老"}}}),
    ]
    import io
    inp, out = io.StringIO("\n".join(lines) + "\n"), io.StringIO()
    assert mcp_server.serve_stdio(inp, out) == 0
    responses = [json.loads(l) for l in out.getvalue().splitlines() if l.strip()]
    assert [r["id"] for r in responses] == [1, 2]  # notification 无响应
    assert responses[1]["result"]["structuredContent"]["result"]["verdict"] == "REFUTED"
