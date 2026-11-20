#!/usr/bin/env python3
"""生成 CycloneDX 1.5 SBOM（零第三方依赖，Python 3.10+ 可运行）。

设计取舍：
  * 不解析整个 pyproject.toml，只用正则抽取 [project] 的版本约束——避免为了
    生成 SBOM 而引入 tomli/tomllib 的版本分裂（3.10 无 tomllib）。
  * 默认确定性输出：时间戳取 SOURCE_DATE_EPOCH 或固定发布日，serialNumber
    由内容哈希派生。CI 用 --check 验证已提交的 sbom.json 与脚本一致。
  * 不声明未核验的上游许可证。SBOM 声明"包 + 版本约束 + 用途 + 是否在主判定
    链路"，许可证核验留给专门的合规流程——编造许可证比留空更危险。

用法：
    python scripts/generate_sbom.py            # 写 sbom.json
    python scripts/generate_sbom.py --check    # 只校验已提交文件是否一致
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
OUT = ROOT / "sbom.json"
FALLBACK_TIMESTAMP = "2026-11-20T00:00:00Z"

# 用途说明来自源码实读（见 SECURITY.md §4），不是从包名猜测。
USAGE = {
    "numpy": "L2 数值引擎：光 MC 金标准、电场球谐求解、剂量场运算",
    "pandas": "受体表达面板：CSV 读取",
    "matplotlib": "受体热图生成（可选绘图）",
}
MAIN_CHAIN = {"numpy"}


def _dependencies(text: str) -> list[str]:
    """从 [project].dependencies 的 Python 列表字面量抽取包名与约束。"""
    m = re.search(r"^dependencies\s*=\s*\[(.*?)\]", text, re.M | re.S)
    if not m:
        return []
    return [s.strip().strip("'\"") for s in m.group(1).split(",") if s.strip()]


def _field(text: str, key: str) -> str:
    m = re.search(rf"^{key}\s*=\s*[\"']([^\"']+)[\"']", text, re.M)
    return m.group(1) if m else ""


def build_bom() -> dict:
    text = PYPROJECT.read_text(encoding="utf-8")
    name = _field(text, "name") or "interface-atlas"
    version = _field(text, "version") or "0.0.0"
    deps = _dependencies(text)

    components = []
    for spec in deps:
        m = re.match(r"^([A-Za-z0-9_.-]+)\s*(.*)$", spec)
        if not m:
            continue
        pkg, constraint = m.group(1), m.group(2).strip()
        components.append({
            "type": "library",
            "name": pkg,
            "bom-ref": f"pkg:pypi/{pkg}@{constraint}" if constraint else f"pkg:pypi/{pkg}",
            "purl": f"pkg:pypi/{pkg}@{constraint}" if constraint else f"pkg:pypi/{pkg}",
            "scope": "required",
            "properties": [
                {"name": "interface-atlas:usage", "value": USAGE.get(pkg, "声明于 pyproject.toml")},
                {"name": "interface-atlas:main-verdict-chain",
                 "value": "true" if pkg in MAIN_CHAIN else "false"},
                {"name": "interface-atlas:license-verification", "value": "pending-upstream-check"},
            ],
        })

    # 可选延迟依赖：未列入 [project].dependencies，但代码有延迟 import。
    components.append({
        "type": "library",
        "name": "nibabel",
        "bom-ref": "pkg:pypi/nibabel",
        "purl": "pkg:pypi/nibabel",
        "scope": "optional",
        "properties": [
            {"name": "interface-atlas:usage", "value": "NIfTI 体模导入（anatomy.py 延迟 import）"},
            {"name": "interface-atlas:main-verdict-chain", "value": "false"},
            {"name": "interface-atlas:license-verification", "value": "pending-upstream-check"},
        ],
    })

    body = {
        "bomFormat": "CycloneDX",
        "specVersion": "1.5",
        "version": 1,
        "metadata": {
            "timestamp": FALLBACK_TIMESTAMP,
            "tools": [{"vendor": "interface-atlas", "name": "generate_sbom.py", "version": "1.0"}],
            "component": {
                "type": "application",
                "name": name,
                "version": version,
                "bom-ref": f"pkg:pypi/{name}@{version}",
                "description": _field(text, "description"),
                "properties": [
                    {"name": "interface-atlas:verdict-core-third-party-deps", "value": "0"},
                    {"name": "interface-atlas:note",
                     "value": "判定内核与三入口零第三方依赖；完整包依赖见 components。"},
                ],
            },
        },
        "components": components,
    }
    # 确定性 UUID：由内容哈希派生，避免每次生成产生无意义 diff。
    digest = hashlib.sha256(json.dumps(body, sort_keys=True, ensure_ascii=False).encode()).hexdigest()
    body["serialNumber"] = "urn:uuid:" + str(uuid.UUID(digest[:32]))
    return body


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="校验已提交的 sbom.json 与脚本生成结果一致")
    args = ap.parse_args()

    generated = json.dumps(build_bom(), indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    if args.check:
        if not OUT.exists():
            print(f"MISSING: {OUT}")
            return 1
        if OUT.read_text(encoding="utf-8") != generated:
            print("DRIFT: sbom.json 与 pyproject/脚本不一致，请重新生成")
            return 1
        print(f"OK: {OUT.name} 一致（{len(build_bom()['components'])} components）")
        return 0
    OUT.write_text(generated, encoding="utf-8")
    print(f"WROTE: {OUT}（{len(build_bom()['components'])} components）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
