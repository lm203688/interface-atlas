"""webapp 的端到端测试：起真实 HTTP 服务，逐端点核对与 CLI 内核一致."""
from __future__ import annotations

import json
import threading
import urllib.error
import urllib.request

import pytest

from interface_atlas import service, webapp


@pytest.fixture(scope="module")
def server():
    srv = webapp.make_server("127.0.0.1", 0)  # 随机空闲端口
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()
    srv.server_close()
    t.join(timeout=5)


def _get(base, path):
    try:
        with urllib.request.urlopen(base + path, timeout=15) as r:
            return r.status, json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read().decode("utf-8"))


def _post(base, path, payload):
    req = urllib.request.Request(base + path, data=json.dumps(payload).encode("utf-8"),
                                 headers={"Content-Type": "application/json"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            return r.status, json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read().decode("utf-8"))


def test_index_page_served(server):
    with urllib.request.urlopen(server + "/", timeout=15) as r:
        assert r.status == 200
        assert "interface-atlas" in r.read().decode("utf-8")


def test_meta_endpoint(server):
    status, d = _get(server, "/api/meta")
    assert status == 200
    assert d["version"] == service.meta()["version"]
    assert d["evidence_count"] == 103


def test_claim_endpoint_matches_kernel(server):
    """API 与 CLI 内核同源性：HTTP 返回值必须等于本地内核返回值."""
    status, d = _post(server, "/api/claim", {"text": "植入式迷走刺激能治疗类风湿关节炎"})
    assert status == 200
    kernel = service.claim_check("植入式迷走刺激能治疗类风湿关节炎")["result"]
    assert d["result"]["verdict"] == kernel["verdict"] == "SUPPORTED"
    assert d["result"]["matched_ids"] == kernel["matched_ids"]


def test_paramcheck_endpoint_out_of_window(server):
    status, d = _post(server, "/api/paramcheck",
                      {"payload": {"modality": "focused ultrasound",
                                   "mechanical_index": "2.4", "intensity": "15W/cm2"}})
    assert status == 200
    assert d["result"]["verdict"] == "OUT_OF_WINDOW"


def test_doseaudit_endpoint(server):
    status, d = _post(server, "/api/doseaudit",
                      {"payload": {"modality": "photobiomodulation", "wavelength": "660nm"}})
    assert status == 200
    assert d["result"]["reproducibility"] in ("not_reproducible", "underdetermined")


def test_dosefield_endpoint_three_modalities(server):
    for modality, params in [("light", {"preset": "muscle_nir"}),
                             ("ultrasound", {"preset": "tfus_transcranial"}),
                             ("electric", {"current_ma": 1.0})]:
        status, d = _post(server, "/api/dosefield", {"modality": modality, "params": params})
        assert status == 200, d
        assert d["result"]
    status, d = _post(server, "/api/dosefield", {"modality": "bogus", "params": {}})
    assert status == 400
    assert "参数错误" in d["error"]


def test_evidence_endpoint_filter(server):
    status, d = _get(server, "/api/evidence?modality=photobiomodulation")
    assert status == 200
    assert 0 < d["result"]["count"] < 103
    status, d = _get(server, "/api/evidence")
    assert d["result"]["count"] == 103


def test_headmodels_and_l3_endpoints(server):
    status, d = _get(server, "/api/headmodels")
    assert status == 200 and len(d["result"]["models"]) == 4
    status, d = _get(server, "/api/l3cases")
    assert status == 200 and len(d["result"]["cases"]) == 4


def test_report_endpoint(server):
    status, d = _post(server, "/api/report", {"text": "红光面板能抗衰老"})
    assert status == 200
    assert d["result"]["verdict"]["verdict"] == "REFUTED"
    assert d["result"]["related_evidence"]


def test_unknown_endpoint_404(server):
    status, d = _get(server, "/api/nope")
    assert status == 404
    assert "endpoints" in d


def test_bad_json_400(server):
    req = urllib.request.Request(server + "/api/claim", data=b"not-json",
                                 headers={"Content-Type": "application/json"}, method="POST")
    with pytest.raises(urllib.error.HTTPError) as e:
        urllib.request.urlopen(req, timeout=15)
    assert e.value.code == 400
