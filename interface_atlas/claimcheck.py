"""ClaimCheck：把一句干预宣称，判定为可追溯、可证伪的结论.

设计原则（继承 vec-trust-kit v0.1 的教训）：
  - 不输出恒真的"不可信"裁决；判定必须取决于被测宣称与证据库内容；
  - 每条判定带命中证据 ID、触发规则、推翻条件、边界声明；
  - 反证单独列出——有 A 级"无效应"证据时输出 REFUTED，而不是含糊的"证据不足"；
  - 宣称超出证据范围（如把体外机制说成人体抗衰）时明确标 EXCEEDED；
  - 无法解析或库中无对应时输出 NO_EVIDENCE，不猜测、不编造。

用法：
  from interface_atlas.claimcheck import check
  result = check("红光面板能抗衰老")
  print(result["verdict"], result["reason"], result["falsification"])
"""
from __future__ import annotations

from . import evidence

VERDICTS = {
    "SUPPORTED": "证据支持：同模态、同终点、同方向，且有人体 A/B 级证据",
    "REFUTED": "证据反驳：同模态同终点存在 A/B 级人体『无效应/相反』证据",
    "EXCEEDED": "证据超出：仅有机制/动物/低等级证据，或宣称终点超出实际测量范围",
    "NO_EVIDENCE": "无证据：库中无对应模态或终点的记录，无法判定",
}

# 归一化词典：所有关键词命中都记录，便于 reason 追溯
MODALITY_KW: dict[str, tuple[str, ...]] = {
    "photobiomodulation": ("光", "红光", "近红外", "nir", "pbm", "photobiomodulation", "激光", "led"),
    "vagus nerve stimulation": ("迷走", "耳夹", "耳刺激", "vns", "tavns", "vagus"),
    "focused ultrasound": ("超声", "聚焦超声", "tfus", "ultrasound"),
    "mechanical loading": ("机械", "运动", "力学", "加载", "剪切力", "mechanical", "exercise"),
    "metabolic / microbiota": ("饮食", "纤维", "scfa", "菌群", "益生菌", "胆汁酸", "代谢", "diet", "fiber"),
}

# 终点归一化 → 与证据库 endpoint 字段匹配的关键词
ENDPOINT_KW: dict[str, tuple[str, ...]] = {
    "cytokine": ("炎症", "抗炎", "tnf", "il-6", "il6", "crp", "inflammation", "cytokine"),
    "cognitive": ("认知", "记忆", "脑功能", "cognitive", "memory"),
    "osteoarthritis": ("骨关节炎", "osteoarthritis", "womac"),
    "rheumatoid": ("类风湿", "rheumatoid"),
    "hair density": ("脱发", "生发", "头发", "hair"),
    "fatigue": ("疲劳", "纤维肌痛", "fatigue", "fibromyalgia"),
    "collagen": ("胶原", "皮肤", "collagen", "skin", "皱纹"),
    "GLP-1": ("血糖", "胰岛素", "glp-1", "glp1", "减重", "饱腹", "insulin", "satiety"),
    "ATP production": ("线粒体", "atp", "能量代谢"),
    "HRV": ("hrv", "心率变异"),
    "anti-aging": ("抗衰", "衰老", "寿命", "健康寿命", "anti-aging", "longevity", "长寿"),
}

INCREASE_KW = ("提高", "增加", "促进", "改善", "增强", "生发", "抗衰", "减少皱纹",
               "increase", "improve", "boost", "enhance", "stimulate")
DECREASE_KW = ("降低", "减少", "抑制", "消炎", "抗炎", "缓解", "decrease", "reduce",
               "lower", "suppress", "inhibit")
CONSUMER_KW = ("消费级", "家用", "耳夹", "贴片", "consumer", "home device")
IMPLANT_KW = ("植入", "implant", "手术")


def _hits(text: str, table: dict[str, tuple[str, ...]]) -> list[str]:
    low = text.lower()
    return [name for name, kws in table.items() if any(kw in low for kw in kws)]


def _parse_claim(text: str) -> dict:
    low = text.lower()
    modalities = _hits(low, MODALITY_KW)
    endpoints = _hits(low, ENDPOINT_KW)
    direction = "increase"
    if any(k in low for k in DECREASE_KW) and not any(k in low for k in INCREASE_KW):
        direction = "decrease"
    elif any(k in low for k in INCREASE_KW):
        direction = "increase"
    return {
        "modalities": modalities,
        "endpoints": endpoints,
        "direction": direction,
        "consumer": any(k in low for k in CONSUMER_KW),
        "implant": any(k in low for k in IMPLANT_KW),
    }


def _endpoint_match(row: dict, target: str) -> bool:
    """终点是否相关。证据库的终点信息分散在 endpoint 与 cell_type 两列，须合并匹配。

    例如 PBM-005 的终点写在 cell_type='knee osteoarthritis'，
    VNS-003 的终点写在 cell_type='rheumatoid arthritis patients'，
    只看 endpoint 列会漏掉这两条关键记录。
    """
    e = f"{row.get('endpoint', '')} {row.get('cell_type', '')}".lower()
    t = target.lower()
    if t == "glp-1":
        return any(k in e for k in ("glp-1", "glp1", "insulin", "satiety", "分泌"))
    if t == "cytokine":
        return any(k in e for k in ("tnf", "il-6", "crp", "cytokine", "inflammatory", "炎症"))
    if t == "anti-aging":
        return False  # 库中没有任何死亡率/健康寿命终点，刻意不映射
    return t in e or e in t


def check(claim: str, rows: list[dict] | None = None) -> dict:
    """判定一条宣称。

    返回 dict：verdict / reason / matched_ids / counter_ids / falsification /
    boundaries / parsed / evidence_fingerprint。
    """
    rows = rows if rows is not None else evidence.load()
    parsed = _parse_claim(claim)
    base = {
        "claim": claim,
        "parsed": parsed,
        "evidence_fingerprint": evidence.fingerprint(rows),
        "boundaries": [
            "证据等级 A/B/C/D 见 README；C 级仅机制，不等于人体疗效",
            "nCPM 是 RNA 水平，不等于膜蛋白丰度",
            "本判定基于证据库当前版本，不代表医疗建议",
        ],
    }
    if not parsed["modalities"] or not parsed["endpoints"]:
        return {**base, "verdict": "NO_EVIDENCE",
                "reason": f"无法解析宣称的模态或终点：modality={parsed['modalities']}, endpoint={parsed['endpoints']}",
                "matched_ids": [], "counter_ids": [],
                "falsification": "把宣称改写为明确的『干预方式 + 生物终点』后重新判定"}

    # anti-aging 是"终点超出"信号：库中无死亡率/健康寿命终点
    anti_aging = "anti-aging" in parsed["endpoints"]
    real_endpoints = [e for e in parsed["endpoints"] if e != "anti-aging"]

    candidates: list[dict] = []
    for r in rows:
        if not any(m in r["modality"].lower() for m in parsed["modalities"]):
            continue
        if real_endpoints and not any(_endpoint_match(r, e) for e in real_endpoints):
            continue
        candidates.append(r)

    if not candidates:
        return {**base, "verdict": "NO_EVIDENCE",
                "reason": f"库中无对应记录：modality={parsed['modalities']}, endpoint={parsed['endpoints']}",
                "matched_ids": [], "counter_ids": [],
                "falsification": "提供该模态/终点的一手研究后重新判定；当前不是『无效』，是『未知』"}

    strong = [r for r in candidates if r["evidence_grade"] in {"A", "B"}]
    human = [r for r in strong if r["human_evidence"] == "yes"]

    # 反证优先：同模态同终点，A/B 级人体证据明确 none
    refuting = [r for r in human if r["effect_direction"] == "none"]
    if refuting:
        return {**base, "verdict": "REFUTED",
                "reason": f"存在 {len(refuting)} 条 A/B 级人体『无效应』证据：" +
                          "; ".join(f"{r['id']}({r['evidence_grade']})" for r in refuting),
                "matched_ids": [], "counter_ids": [r["id"] for r in refuting],
                "falsification": "需要针对该人群/剂量窗的新 RCT 推翻现有 meta 结论；"
                                 "注意反证针对的是群体平均效应，不排除亚组获益"}

    # 抗衰/延寿是库中没有的直接终点。即使同时命中胶原等替代终点，
    # 也不允许把替代终点升级为"抗衰有效"——这是最常见的越界宣称。
    if anti_aging:
        return {**base, "verdict": "EXCEEDED",
                "reason": "宣称含『抗衰/延寿』，但库中无任何死亡率/健康寿命/衰老终点；"
                          "胶原合成等替代终点不能外推为抗衰疗效",
                "matched_ids": [r["id"] for r in candidates], "counter_ids": [],
                "falsification": "需要以死亡率为终点的长期人体 RCT；替代终点（胶原、炎症因子）不足以升级判定"}

    # 方向一致的人体强证据
    agreeing = [r for r in human if r["effect_direction"] == parsed["direction"]]
    if agreeing:
        ids = [r["id"] for r in agreeing]
        return {**base, "verdict": "SUPPORTED",
                "reason": f"命中 {len(agreeing)} 条 A/B 级人体同向证据：" + ", ".join(ids),
                "matched_ids": ids, "counter_ids": [],
                "falsification": "若新证据显示效应量不显著或人群不匹配，判定应降级；"
                                 "注意证据不支持『抗衰/延寿』这类未测量的远端终点"}

    # 消费级宣称但证据来自植入式：模态相同但干预层级不同
    if parsed["consumer"] and any(r["signal"].lower().find("implant") >= 0 or
                                  "implant" in r["notes"].lower() for r in candidates):
        return {**base, "verdict": "EXCEEDED",
                "reason": "宣称是消费级设备，但同模态强证据来自植入式设备；"
                          "消费级耳夹可能未激活目标神经",
                "matched_ids": [r["id"] for r in candidates], "counter_ids": [],
                "falsification": "提供消费级设备针对同一终点的假刺激对照 RCT 可升级判定"}

    # 其余：仅有机制/低等级证据，或宣称超出测量终点
    weak = [r["id"] for r in candidates if r["evidence_grade"] in {"C", "D"}]
    if anti_aging:
        reason = "宣称含『抗衰/延寿』，但库中没有任何死亡率或健康寿命终点；" \
                 "现有证据最多到疾病症状或生物标志物层面"
    else:
        reason = f"仅命中 {len(weak)} 条 C/D 级或机制证据，无 A/B 级人体同向证据"
    return {**base, "verdict": "EXCEEDED",
            "reason": reason,
            "matched_ids": [r["id"] for r in candidates], "counter_ids": [],
            "falsification": "需要同模态、同终点、同方向的 A/B 级人体证据；"
                             "仅机制清楚不足以支持人体宣称"}


def to_markdown(result: dict) -> str:
    lines = [
        f"## ClaimCheck：{result['claim']}",
        "",
        f"**判定**：{result['verdict']} — {VERDICTS[result['verdict']]}",
        "",
        f"**依据**：{result['reason']}",
        "",
        f"**命中证据**：{', '.join(result['matched_ids']) or '无'}",
        f"**反证证据**：{', '.join(result['counter_ids']) or '无'}",
        f"**证据库指纹**：`{result['evidence_fingerprint']}`",
        "",
        f"**推翻条件**：{result['falsification']}",
        "",
        "**边界**：",
    ]
    lines += [f"- {b}" for b in result["boundaries"]]
    return "\n".join(lines)
