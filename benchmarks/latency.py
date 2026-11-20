#!/usr/bin/env python3
"""判定延迟基准 + 确定性检查（零第三方依赖）。

SLO（来自生产级提升报告 §3.1）：
  * 单机判定 P95 < 200 ms
  * 同一版本内同输入 100% 确定性（结果逐字节一致）

用法：
    python benchmarks/latency.py            # 实测并写 baseline.json
    python benchmarks/latency.py --check    # 只校验 SLO，不写文件（CI 用）
"""
from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

from interface_atlas import service

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "benchmarks" / "baseline.json"
SLO_P95_MS = 200.0
REPEAT = 200
WARMUP = 5

CASES = {
    "meta": lambda: service.meta(),
    "claim_check": lambda: service.claim_check("植入式迷走刺激能治疗类风湿关节炎"),
    "paramcheck": lambda: service.param_check(
        {"modality": "focused ultrasound", "mechanical_index": 1.2}),
    "dose_audit": lambda: service.dose_audit({"modality": "photobiomodulation"}),
    "evidence_query": lambda: service.evidence_query(limit=50),
    "dose_field_light": lambda: service.dose_field("light", {"preset": "muscle_nir", "target_depth_cm": 1.0}),
    "l3_cases": lambda: service.l3_cases(),
    "report": lambda: service.report("红光面板能抗衰老"),
}


def _percentile(xs: list[float], q: float) -> float:
    s = sorted(xs)
    i = min(len(s) - 1, max(0, int(round(q * (len(s) - 1)))))
    return s[i]


def measure(cases: dict, repeat: int = REPEAT) -> dict:
    out = {}
    for name, fn in cases.items():
        for _ in range(WARMUP):
            fn()
        samples = []
        for _ in range(repeat):
            t0 = time.perf_counter()
            fn()
            samples.append((time.perf_counter() - t0) * 1000.0)
        out[name] = {
            "n": repeat,
            "p50_ms": round(_percentile(samples, 0.50), 4),
            "p95_ms": round(_percentile(samples, 0.95), 4),
            "p99_ms": round(_percentile(samples, 0.99), 4),
            "max_ms": round(max(samples), 4),
        }
    return out


def determinism_check(repeat: int = 25) -> dict:
    """同输入重复调用，JSON 序列化结果必须逐字节一致。"""
    out = {}
    for name, fn in CASES.items():
        first = json.dumps(fn(), ensure_ascii=False, sort_keys=True)
        stable = all(json.dumps(fn(), ensure_ascii=False, sort_keys=True) == first
                     for _ in range(repeat - 1))
        out[name] = {"stable": bool(stable), "bytes": len(first.encode("utf-8"))}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true", help="只校验 SLO，不写 baseline.json")
    ap.add_argument("--repeat", type=int, default=REPEAT)
    args = ap.parse_args()

    measured = measure(CASES, args.repeat)
    det = determinism_check()
    worst = max((v["p95_ms"], k) for k, v in measured.items())
    all_stable = all(v["stable"] for v in det.values())
    slo_ok = worst[0] < SLO_P95_MS

    report = {
        "slo": {"p95_ms_ceiling": SLO_P95_MS, "deterministic": True},
        "environment": {"note": "本机单次实测；跨机器只比阈值，不与基线数值硬比对"},
        "repeat": args.repeat,
        "latency_ms": measured,
        "determinism": det,
        "worst_case": {"case": worst[1], "p95_ms": worst[0]},
        "passed": bool(slo_ok and all_stable),
    }

    print(f"worst P95 = {worst[0]:.3f} ms ({worst[1]})  ceiling = {SLO_P95_MS:.0f} ms")
    print(f"determinism = {'PASS' if all_stable else 'FAIL'} ({len(det)} cases)")
    print(f"SLO = {'PASS' if report['passed'] else 'FAIL'}")

    if args.check:
        return 0 if report["passed"] else 1

    BASELINE.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"WROTE: {BASELINE}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
