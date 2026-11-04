#!/usr/bin/env python3
"""Local server for DPDP review pages with in-page decisions.

Serves <data>/<run>/review/*.html at http://127.0.0.1:<port>/<run>/review/<step>.html and accepts the analyst's decision from the
page's form (POST /api/decision). The decision is written to <data>/<run>/review/pending_decision.json; Claude waits for that
file (pipeline/wait_decision.py), applies the decision with the step script, and the step script archives it (common.R).

Safety: binds 127.0.0.1 only; writes only pending_decision.json inside an existing run folder under <data>; a random token made at
start is put into served pages and must come back in the X-DPDP-Token header (other websites can't read it, and the custom header
blocks simple cross-site posts). No other files are written.

Usage: python3 pipeline/review_server.py [--data /home/user/Claude_Shiny/data] [--port 8765]
"""
import argparse, json, os, re, secrets, sys, time
from http.server import ThreadingHTTPServer, BaseHTTPRequestHandler
from urllib.parse import urlparse, parse_qs, unquote

STEPS = {"01_load", "03_parameters", "04_noise", "05_filter", "06_normalize", "07_impute", "08_rollup", "09_qc", "10_stats"}
RUN_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FILE_RE = re.compile(r"^[A-Za-z0-9_.-]+\.html$")
MAX_BODY = 20000


def make_handler(data_root, token):
    class H(BaseHTTPRequestHandler):
        server_version = "dpdp-review/1"

        def log_message(self, fmt, *args):
            sys.stderr.write("%s %s\n" % (time.strftime("%H:%M:%S"), fmt % args))

        def send(self, code, body, ctype="application/json"):
            b = body.encode() if isinstance(body, str) else body
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(b)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(b)

        def run_dir(self, run):
            if not run or not RUN_RE.match(run):
                return None
            d = os.path.join(data_root, run)
            return d if os.path.isdir(os.path.join(d, "review")) else None

        def do_GET(self):
            u = urlparse(self.path)
            parts = [unquote(p) for p in u.path.split("/") if p]
            if u.path == "/api/status":
                q = parse_qs(u.query)
                d = self.run_dir((q.get("run") or [""])[0]); step = (q.get("step") or [""])[0]
                if not d or step not in STEPS:
                    return self.send(400, json.dumps({"error": "bad run or step"}))
                pend = os.path.join(d, "review", "pending_decision.json")
                pending = None
                if os.path.exists(pend):
                    try:
                        pending = json.load(open(pend)).get("step")
                    except Exception:
                        pending = "?"
                page = os.path.join(d, "review", step + ".html")
                return self.send(200, json.dumps({"pending_step": pending,
                                                  "page_mtime": os.path.getmtime(page) if os.path.exists(page) else 0}))
            if not parts:   # index of runs
                runs = sorted(r for r in os.listdir(data_root) if os.path.isfile(os.path.join(data_root, r, "review", "index.html")))
                items = "".join('<li><a href="/%s/review/index.html">%s</a></li>' % (r, r) for r in runs)
                return self.send(200, "<!doctype html><meta charset=utf-8><title>DPDP runs</title><h1>DPDP runs</h1><ul>%s</ul>" % items, "text/html; charset=utf-8")
            if len(parts) == 3 and parts[1] == "review" and FILE_RE.match(parts[2]):
                d = self.run_dir(parts[0])
                f = os.path.join(d, "review", parts[2]) if d else None
                if f and os.path.isfile(f):
                    html = open(f, "rb").read().replace(b"__DPDP_TOKEN__", token.encode())
                    return self.send(200, html, "text/html; charset=utf-8")
            return self.send(404, json.dumps({"error": "not found"}))

        def do_POST(self):
            if urlparse(self.path).path != "/api/decision":
                return self.send(404, json.dumps({"error": "not found"}))
            if self.headers.get("X-DPDP-Token") != token:
                return self.send(403, json.dumps({"error": "bad token: reload the page from this server"}))
            n = int(self.headers.get("Content-Length") or 0)
            if n <= 0 or n > MAX_BODY:
                return self.send(400, json.dumps({"error": "bad body size"}))
            try:
                dec = json.loads(self.rfile.read(n))
            except Exception:
                return self.send(400, json.dumps({"error": "not JSON"}))
            d = self.run_dir(str(dec.get("run", ""))); step = str(dec.get("step", ""))
            if not d or step not in STEPS or not isinstance(dec.get("choice"), str):
                return self.send(400, json.dumps({"error": "bad run, step or choice"}))
            rec = {"run": os.path.basename(d), "step": step, "choice": dec["choice"][:200],
                   "fields": {str(k)[:60]: str(v)[:500] for k, v in (dec.get("fields") or {}).items()},
                   "note": str(dec.get("note", ""))[:500], "analyst": str(dec.get("analyst", ""))[:100],
                   "submitted": time.strftime("%Y-%m-%d %H:%M:%S"), "source": "review page"}
            path = os.path.join(d, "review", "pending_decision.json")
            tmp = path + ".tmp"
            with open(tmp, "w") as f:
                json.dump(rec, f, indent=1)
            os.replace(tmp, path)
            self.log_message("decision %s %s: %s", rec["run"], step, rec["choice"])
            return self.send(200, json.dumps({"ok": True, "submitted": rec["submitted"]}))
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.environ.get("DPDP_DATA", "/home/user/Claude_Shiny/data"))
    ap.add_argument("--port", type=int, default=int(os.environ.get("DPDP_PORT", "8765")))
    a = ap.parse_args()
    data_root = os.path.realpath(a.data)
    token = secrets.token_urlsafe(24)
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), make_handler(data_root, token))
    print("DPDP review server: http://127.0.0.1:%d/  (data: %s)" % (a.port, data_root), flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
