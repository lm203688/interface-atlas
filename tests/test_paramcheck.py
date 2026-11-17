from interface_atlas import paramcheck, evidence


def _pbm_complete():
    return {
        "modality": "photobiomodulation",
        "wavelength": "660 nm",
        "power_density": "25 mW/cm2",
        "fluence": "6 J/cm2",
        "exposure_time": "240 s",
        "spot_area": "1 cm2",
        "cell_type": "fibroblasts",
        "endpoint": "collagen synthesis",
        "source_url": "https://doi.org/10.1089/pho.2018.4455",
    }


def test_complete_report_passes():
    rep = paramcheck.check(_pbm_complete())
    assert rep.verdict == "COMPLETE"
    assert rep.passed
    assert not any(c.status == "FAIL" for c in rep.checks)


def test_missing_key_parameter_changes_verdict():
    rep = paramcheck.check(_pbm_complete())
    partial = dict(_pbm_complete())
    partial.pop("power_density")
    rep2 = paramcheck.check(partial)
    assert rep2.verdict == "INSUFFICIENT"      # 判定随报告改变，不是恒真
    assert rep2.verdict != rep.verdict
    assert any(c.name == "key_parameters" and c.status == "FAIL" for c in rep2.checks)


def test_relative_dose_is_rejected_not_invented():
    partial = dict(_pbm_complete())
    partial["fluence"] = "low dose"
    rep = paramcheck.check(partial)
    dose_check = next(c for c in rep.checks if c.name == "dose")
    assert dose_check.status == "FAIL"
    assert rep.verdict == "INSUFFICIENT"


def test_tfus_mi_over_limit_is_out_of_window():
    fus = {
        "modality": "focused ultrasound",
        "frequency": "0.5 MHz",
        "mechanical_index": "2.4",
        "intensity": "15 W/cm2",
        "pulse_duration": "0.5 ms",
        "duty_cycle": "70%",
        "temperature_rise": "1.2 C",
        "cell_type": "cortical neurons",
        "endpoint": "motor cortex modulation",
        "source_url": "https://doi.org/10.1038/s41598-020-62265-8",
    }
    rep = paramcheck.check(fus)
    assert rep.verdict == "OUT_OF_WINDOW"
    assert not rep.passed
    assert any(c.name == "hard_limits" and "1.9" in c.detail for c in rep.checks)
    # 只把 MI 改回窗内，判定必须变化
    fus["mechanical_index"] = "1.2"
    rep2 = paramcheck.check(fus)
    assert rep2.verdict != "OUT_OF_WINDOW"


def test_temperature_rise_hard_limit():
    fus = {
        "modality": "focused ultrasound",
        "frequency": "1 MHz", "mechanical_index": "1.0",
        "intensity": "10 W/cm2", "pulse_duration": "1 ms",
        "duty_cycle": "50%", "temperature_rise": "3.5 C",
        "source_url": "https://doi.org/10.1016/j.brs.2025.10.007",
    }
    rep = paramcheck.check(fus)
    assert rep.verdict == "OUT_OF_WINDOW"
    assert any("温升" in c.detail for c in rep.checks if c.name == "hard_limits")


def test_unknown_modality_is_no_evidence():
    rep = paramcheck.check({"modality": "quantum frequency healing", "dose": "62-72 Hz"})
    assert rep.verdict == "NO_EVIDENCE"
    assert any(c.status == "FAIL" for c in rep.checks)


def test_batch_over_evidence_library_is_falsifiable():
    rows = evidence.load()
    reports = [paramcheck.check_evidence_row(r) for r in rows]
    verdicts = {r.verdict for r in reports}
    # 103 条里必须出现多种判定，不能全部 PASS 或全部 FAIL（防止恒真/硬编码）
    assert len(verdicts) >= 3
    assert "INSUFFICIENT" in verdicts  # 文献摘要普遍报不全参数 → 如实判不足
    assert "NO_EVIDENCE" in verdicts   # 无参数的纯宣称 → 判无证据
    # 判定必须与可核对参数量相关，不能随机：NO_EVIDENCE 的应确实没有参数
    for rep in reports:
        if rep.verdict == "NO_EVIDENCE":
            assert all(c.status == "FAIL" for c in rep.checks if c.name == "dose")
    assert all(rep.falsification for rep in reports)


def test_output_includes_disclaimer_and_falsification():
    rep = paramcheck.check(_pbm_complete())
    text = paramcheck.render(rep)
    assert "本地自查" in text and "推翻条件" in text
    assert rep.to_dict()["verdict"] == "COMPLETE"
