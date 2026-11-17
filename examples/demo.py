"""最小演示：查接收器、查证据、做交叉核验、判定宣称."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from interface_atlas import claimcheck, evidence, link, receptors  # noqa: E402

print("=== A1：α7nAChR 高表达细胞 ===")
print(receptors.top_cells(receptors.load_matrix(), gene="CHRNA7", n=5).to_string(index=False))

print("\n=== A2：A 级证据 ===")
for r in evidence.query(grade_min="A"):
    print(f"{r['id']} | {r['cell_type']} | {r['endpoint']} | {r['effect_direction']}")

print("\n=== link：无法核验的记录数 ===")
rows = link.crosslink()
print(link.summarize(rows))

print("\n=== ClaimCheck：宣称-证据判定 ===")
for claim in ["红光面板能抗衰老延寿", "植入式迷走刺激能治疗类风湿关节炎",
              "特定频率电磁波能修复细胞"]:
    r = claimcheck.check(claim)
    print(f"[{r['verdict']}] {claim} → {r['reason'][:60]}")
