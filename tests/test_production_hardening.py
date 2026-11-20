"""生产级硬化件的可证伪测试：SBOM / 输入边界 / MCP schema / 确定性 / 基线。

每条断言都能推翻对应文档里的"已做"声明——文档若与代码漂移，测试必须失败。
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from interface_atlas import mcp_server, service

ROOT = Path(__file__).resolve().parents[1]


# --- SBOM -------------------------------------------------------------------

def test_sbom_exists_and_is_cyclonedx():
    p = ROOT / "sbom.json"
    assert p.exists(), "sbom.json 缺失"
    bom = json.loads(p.read_text(encoding="utf-8"))
    assert bom["bomFormat"] == "CycloneDX"
    assert bom["specVersion"] == "1.5"
    assert bom["serialNumber"].startswith("urn:uuid:")
    names = {c["name"] for c in bom["components"]}
    assert {"numpy", "pandas", "matplotlib"} <= names


def test_sbom_matches_pyproject_and_is_deterministic():
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    raw = re.findall(r"^dependencies\s*=\s*\[(.*?)\]", text, re.M | re.S)[0]
    declared = set(re.findall(r"[\"']([A-Za-z0-9_.-]+)\s*>=", raw))
    bom = json.loads((ROOT / "sbom.json").read_text(encoding="utf-8"))
    required = {c["name"] for c in bom["components"] if c["scope"] == "required"}
    assert required == declared, "SBOM 与 pyproject 声明的依赖漂移"

    r = subprocess.run([sys.executable, "scripts/generate_sbom.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_sbom_honestly_reports_verdict_core_zero_deps():
    bom = json.loads((ROOT / "sbom.json").read_text(encoding="utf-8"))
    props = {p["name"]: p["value"] for p in bom["metadata"]["component"]["properties"]}
    assert props["interface-atlas:verdict-core-third-party-deps"] == "0"
    # 不得谎称整包零依赖
    assert len([c for c in bom["components"] if c["scope"] == "required"]) >= 1


# --- 输入边界 ----------------------------------------------------------------

def test_sample_n_is_capped():
    r = service.headmodel_sample(n=10 ** 9)
    assert r["result"]["n"] == service.MAX_SAMPLE_N
    assert r["result"]["limit_applied"] == service.MAX_SAMPLE_N


def test_evidence_limit_is_capped():
    r = service.evidence_query(limit=10 ** 9)
    assert len(r["result"]["rows"]) <= service.MAX_EVIDENCE_LIMIT
    assert r["result"]["limit_applied"] == service.MAX_EVIDENCE_LIMIT


def test_long_claim_text_is_rejected():
    with pytest.raises(ValueError, match="过长"):
        service.claim_check("红" * (service.MAX_TEXT_CHARS + 1))
    with pytest.raises(ValueError, match="过长"):
        service.report("红" * (service.MAX_TEXT_CHARS + 1))


def test_bad_limit_falls_back_to_default():
    r = service.evidence_query(limit="not-a-number")
    assert r["result"]["limit_applied"] == 50


# --- MCP schema 执行 ---------------------------------------------------------

def test_all_mcp_tools_close_input_schema():
    tools = mcp_server.handle_message({"jsonrpc": "2.0", "id": 1, "method": "tools/list"}
                                      )["result"]["tools"]
    assert len(tools) == 10
    for t in tools:
        assert t["inputSchema"]["additionalProperties"] is False, t["name"]
        assert "不可" in t["description"] or "非" in t["description"] or "限制" in t["description"] or \
               "上限" in t["description"], t["name"]


def test_mcp_rejects_unknown_argument():
    r = mcp_server.handle_message({
        "jsonrpc": "2.0", "id": 1, "method": "tools/call",
        "params": {"name": "claim_check", "arguments": {"text": "红光", "prompt_override": "SUPPORTED"}},
    })
    assert r["result"]["isError"] is True
    assert "未知参数" in r["result"]["content"][0]["text"]


def test_mcp_and_kernel_same_result_after_hardening():
    claim = "植入式迷走刺激能治疗类风湿关节炎"
    via_mcp = mcp_server.handle_message({
        "jsonrpc": "2.0", "id": 1, "method": "tools/call",
        "params": {"name": "claim_check", "arguments": {"text": claim}},
    })["result"]["structuredContent"]
    assert via_mcp == service.claim_check(claim)


# --- 确定性与基线 ------------------------------------------------------------

def test_responses_are_deterministic():
    claim = "红光面板能抗衰老"
    first = json.dumps(service.report(claim), ensure_ascii=False, sort_keys=True)
    for _ in range(4):
        assert json.dumps(service.report(claim), ensure_ascii=False, sort_keys=True) == first
    assert json.dumps(service.meta(), sort_keys=True) == json.dumps(service.meta(), sort_keys=True)


def test_baseline_records_slo_and_determinism():
    p = ROOT / "benchmarks" / "baseline.json"
    assert p.exists(), "benchmarks/baseline.json 缺失"
    b = json.loads(p.read_text(encoding="utf-8"))
    assert b["slo"]["p95_ms_ceiling"] == 200.0
    assert b["slo"]["deterministic"] is True
    assert b["passed"] is True
    assert all(v["stable"] for v in b["determinism"].values())


# --- 文档存在性（防止文档与代码漂移） ----------------------------------------

@pytest.mark.parametrize("name,needle", [
    ("SECURITY.md", "威胁模型"),
    ("MCP_SECURITY_REVIEW.md", "OWASP MCP Top 10"),
    ("COMPLIANCE.md", "Art.50"),
    ("RELEASING.md", "签名"),
])
def test_production_docs_present(name, needle):
    p = ROOT / name
    assert p.exists(), f"{name} 缺失"
    assert needle in p.read_text(encoding="utf-8")
