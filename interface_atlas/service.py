"""产品化服务层：CLI / Web API / MCP 三个入口共享的同一内核.

设计原则（继承项目质量门禁）：
  1. 本模块只做"搬运与呈现" —— 所有判定逻辑仍在 claimcheck/paramcheck/dosimetry 等内核模块；
  2. 每个响应都带证据指纹、库版本与边界声明，可被外部审计；
  3. 所有输出 JSON 可序列化，不给不可信的跨进程对象。
"""
from __future__ import annotations

from . import __version__, anatomy, claimcheck, dosimetry, evidence, l3bench, paramcheck
from . import light_field, ultrasound_field as usfield, electric_field as esfield

BOUNDARY = (
    "科研降阶模型与社区草案规范，非临床处方、非监管结论、不构成医疗建议；"
    "自发报告≠发生率；RNA≠膜蛋白；设备参数≠靶区剂量。"
)


def meta() -> dict:
    """服务元信息：版本 + 证据指纹 + 边界声明（每次响应必带）."""
    rows = evidence.load()
    return {
        "library": "interface-atlas",
        "version": __version__,
        "evidence_fingerprint": evidence.fingerprint(rows),
        "evidence_count": len(rows),
        "boundary": BOUNDARY,
    }


def _wrap(payload: dict) -> dict:
    return {"meta": meta(), "result": payload}


def claim_check(text: str) -> dict:
    """宣称 → 四级判定（SUPPORTED/REFUTED/EXCEEDED/NO_EVIDENCE）+ 推翻条件."""
    if not text or not text.strip():
        raise ValueError("宣称文本不能为空")
    return _wrap(claimcheck.check(text.strip()))


def param_check(payload: dict, source_id: str = "") -> dict:
    """参数报告 → 五级预检（COMPLETE/PARTIAL/INSUFFICIENT/OUT_OF_WINDOW/NO_EVIDENCE）."""
    if not isinstance(payload, dict) or not payload.get("modality"):
        raise ValueError("payload 需为含 modality 字段的对象")
    rep = paramcheck.check(payload, source_id=source_id or str(payload.get("id", "")))
    return _wrap(__import__("json").loads(paramcheck.to_json(rep)))


def dose_audit(payload: dict, source_id: str = "") -> dict:
    """剂量报告 → T1–T4 完整性审计 + 可复现性四档判定."""
    if not isinstance(payload, dict) or not payload.get("modality"):
        raise ValueError("payload 需为含 modality 字段的对象")
    rep = dosimetry.audit(payload, source_id=source_id or str(payload.get("id", "")))
    return _wrap(__import__("json").loads(dosimetry.to_json(rep)))


def evidence_query(modality: str | None = None, gene: str | None = None,
                   cell_type: str | None = None, grade_min: str | None = None,
                   human_only: bool = False, keyword: str | None = None,
                   limit: int = 50) -> dict:
    """证据库查询（只读）."""
    rows = evidence.load()
    rows = evidence.query(rows, modality=modality, gene=gene, cell_type=cell_type,
                          grade_min=grade_min, human_only=human_only, keyword=keyword)
    return _wrap({"count": len(rows), "rows": rows[: int(limit)]})


def dose_field(modality: str, params: dict) -> dict:
    """三模态剂量场：设备/组织参数 → 靶区剂量（降阶代理，如实标注边界）.

    modality: light | ultrasound | electric
    """
    params = dict(params or {})
    if modality == "light":
        preset = params.get("preset")
        if preset:
            if preset not in light_field.TISSUE_PRESETS:
                raise ValueError(f"未知组织预设 {preset}；可用：{', '.join(light_field.TISSUE_PRESETS)}")
            tissue = light_field.TISSUE_PRESETS[preset]
        else:
            tissue = light_field.TissueOptics(
                float(params.get("mu_a", 0.3)), float(params.get("mu_sp", 70.0)),
                float(params.get("g", 0.90)), float(params.get("n", 1.40)), "custom",
            )
        res = light_field.dose_field(
            tissue,
            depths_cm=params.get("depths_cm"),
            max_depth_cm=float(params.get("max_depth_cm", 2.0)),
            n_points=int(params.get("n_points", 41)),
            target_depth_cm=float(params.get("target_depth_cm", 0.5)),
            geometry=params.get("geometry", "plane"),
        )
        return _wrap(res.as_dict())
    if modality == "ultrasound":
        preset = params.get("preset", "tfus_transcranial")
        if preset not in usfield.TRANSDUCER_PRESETS:
            raise ValueError(f"未知换能器预设 {preset}；可用：{', '.join(usfield.TRANSDUCER_PRESETS)}")
        tr = usfield.TRANSDUCER_PRESETS[preset]
        skull_name = params.get("skull", "skull_thin")
        if skull_name not in usfield.TISSUE_PRESETS:
            raise ValueError(f"未知组织预设 {skull_name}；可用：{', '.join(usfield.TISSUE_PRESETS)}")
        layers = (usfield.TISSUE_PRESETS[skull_name],)
        res = usfield.dose_field(
            tr, layers,
            duty_cycle=float(params.get("duty_cycle", 1.0)),
            pulse_duration_s=float(params.get("pulse_duration_s", 1.0)),
        )
        return _wrap(res.as_dict())
    if modality == "electric":
        layers = esfield.UNIFORM_SPHERE if params.get("preset") == "uniform" else esfield.HEAD_4SPHERE
        res = esfield.dose_field(
            float(params.get("current_ma", 1.0)),
            layers=layers,
            target_depth_cm=float(params.get("target_depth_cm", 2.2)),
            electrode_area_cm2=float(params.get("electrode_area_cm2", 35.0)),
            n_harmonics=int(params.get("n_harmonics", 32)),
        )
        return _wrap(res.as_dict())
    raise ValueError(f"未知模态 {modality}；可用：light / ultrasound / electric")


def headmodel_list() -> dict:
    """内置代理体模清单（含许可门禁结论）."""
    out = []
    for name, m in anatomy.BUILTIN_MODELS.items():
        ok, reason = anatomy.check_license(m, commercial=True)
        d = m.as_dict()
        d["name"] = name
        d["commercial_ok"] = ok
        d["license_note"] = reason
        out.append(d)
    return _wrap({"models": out, "note": "内置体模为解析/统计代理，非真实患者影像；真实头模由用户接入并遵守数据许可。"})


def headmodel_sample(n: int = 100, seed: int = 20261119) -> dict:
    """群体解剖抽样（illustrative 示例分布，固定 seed 可复现）."""
    models = anatomy.population_sample(n=int(n), seed=int(seed))
    skulls = [m.skull_thickness_mm() for m in models]
    return _wrap({
        "n": len(models), "seed": int(seed),
        "skull_thickness_mm": {"min": min(skulls), "max": max(skulls),
                               "mean": sum(skulls) / len(skulls)},
        "note": "示例分布非经验值；真实群体统计用 load_population_stats() 替换。",
    })


def headmodel_license(name: str, commercial: bool = False) -> dict:
    """fail-closed 许可门禁."""
    m = anatomy.get_model(name)
    ok, reason = anatomy.check_license(m, commercial=bool(commercial))
    return _wrap({"model": name, "commercial": bool(commercial), "allowed": ok, "reason": reason})


def l3_cases() -> dict:
    """L3 基准种子案例（每条带推翻条件）."""
    return _wrap({"cases": [c.as_dict() for c in l3bench.CASES],
                  "protocol_note": "多指标协议，拒绝单一 MAE 与聚合总分；预印本/会议条目已标注。"})


def report(text: str) -> dict:
    """完整评估旅程：宣称 → 判定 → 剂量完整性 → 相关证据（一次请求跑通）."""
    claim = claimcheck.check(text.strip())
    matched = claim.get("matched_ids", []) + claim.get("counter_ids", [])
    rows = evidence.load()
    related = [r for r in rows if r.get("id") in matched][:10]
    return _wrap({
        "claim": text.strip(),
        "verdict": claim,
        "related_evidence": related,
        "journey": ["claim", "evidence", "dose_audit", "dose_field", "l3_bench"],
        "note": "剂量场补算与基准比对请用 dose_field / l3_cases 端点按模态细化。",
    })
