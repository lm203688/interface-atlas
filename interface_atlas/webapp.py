"""interface-atlas 零依赖 Web 服务：REST API + 单页交互站.

纯标准库 http.server，无第三方依赖；与 CLI / MCP 完全同源（都走 service.py）。
仅绑定 127.0.0.1（本地工具，不做公网暴露；公网部署请自行加认证与 TLS）。

启动：
    python -m interface_atlas.webapp --port 7100
"""
from __future__ import annotations

import argparse
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from . import service

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
}


def _json_bytes(payload: dict, status: int = 200) -> tuple[int, bytes, str]:
    return status, json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8"), "application/json; charset=utf-8"


def _bad_request(exc: ValueError) -> tuple[int, bytes, str]:
    return _json_bytes({"error": f"参数错误：{exc}", "meta": service.meta()}, status=400)


def route(method: str, path: str, query: dict, body: dict):
    """路由 → (status, body, content_type)；只读 GET + 计算型 POST，无写操作."""
    if method == "GET" and path in ("/", "/index.html"):
        index = WEB_DIR / "index.html"
        if index.exists():
            return 200, index.read_bytes(), "text/html; charset=utf-8"
        return _json_bytes({"error": "web/index.html 缺失"}, status=404)

    if method == "POST" and path == "/api/claim":
        return _json_bytes(service.claim_check(body.get("text", "")))
    if method == "POST" and path == "/api/paramcheck":
        return _json_bytes(service.param_check(body.get("payload", {}), body.get("source_id", "")))
    if method == "POST" and path == "/api/doseaudit":
        return _json_bytes(service.dose_audit(body.get("payload", {}), body.get("source_id", "")))
    if method == "POST" and path == "/api/dosefield":
        return _json_bytes(service.dose_field(body.get("modality", ""), body.get("params", {})))
    if method == "POST" and path == "/api/report":
        return _json_bytes(service.report(body.get("text", "")))

    if method == "GET" and path == "/api/evidence":
        q = {k: v[0] for k, v in query.items() if v and v[0]}
        return _json_bytes(service.evidence_query(
            modality=q.get("modality"), gene=q.get("gene"), cell_type=q.get("cell_type"),
            grade_min=q.get("grade_min"), human_only=q.get("human_only") == "true",
            keyword=q.get("keyword"), limit=int(q.get("limit", 50))))
    if method == "GET" and path == "/api/headmodels":
        return _json_bytes(service.headmodel_list())
    if method == "GET" and path == "/api/headmodel/sample":
        return _json_bytes(service.headmodel_sample(
            n=int(query.get("n", ["100"])[0]), seed=int(query.get("seed", ["20261119"])[0])))
    if method == "GET" and path == "/api/l3cases":
        return _json_bytes(service.l3_cases())
    if method == "GET" and path == "/api/meta":
        return _json_bytes(service.meta())

    return _json_bytes({"error": f"未知端点 {method} {path}",
                        "endpoints": ["/api/claim", "/api/paramcheck", "/api/doseaudit",
                                      "/api/dosefield", "/api/evidence", "/api/report",
                                      "/api/headmodels", "/api/l3cases", "/api/meta"]},
                       status=404)


class Handler(BaseHTTPRequestHandler):
    server_version = "interface-atlas-web"

    def log_message(self, fmt, *args):  # 静默访问日志，保持测试输出干净
        pass

    def _send(self, status: int, body: bytes, content_type: str):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        try:
            status, body, ctype = route("GET", parsed.path, parse_qs(parsed.query), {})
        except ValueError as exc:
            status, body, ctype = _bad_request(exc)
        except Exception as exc:  # 服务不得因单次请求崩溃
            status, body, ctype = _json_bytes({"error": f"内部错误：{type(exc).__name__}: {exc}"}, status=500)
        self._send(status, body, ctype)

    def do_POST(self):
        parsed = urlparse(self.path)
        try:
            length = int(self.headers.get("Content-Length", 0))
            body = json.loads(self.rfile.read(length).decode("utf-8")) if length else {}
        except (json.JSONDecodeError, UnicodeDecodeError):
            self._send(400, json.dumps({"error": "请求体需为 JSON"}, ensure_ascii=False).encode("utf-8"),
                       "application/json; charset=utf-8")
            return
        try:
            status, payload, ctype = route("POST", parsed.path, parse_qs(parsed.query), body)
        except ValueError as exc:
            status, payload, ctype = _bad_request(exc)
        except Exception as exc:
            status, payload, ctype = _json_bytes({"error": f"内部错误：{type(exc).__name__}: {exc}"}, status=500)
        self._send(status, payload, ctype)


def make_server(host: str = "127.0.0.1", port: int = 7100) -> ThreadingHTTPServer:
    return ThreadingHTTPServer((host, port), Handler)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="interface-atlas-web", description="interface-atlas 本地 Web 服务（仅 127.0.0.1）")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=7100)
    args = p.parse_args(argv)
    srv = make_server(args.host, args.port)
    print(f"interface-atlas Web 服务已启动：http://{args.host}:{args.port}/ （Ctrl+C 停止）")
    print("边界：" + service.BOUNDARY)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
