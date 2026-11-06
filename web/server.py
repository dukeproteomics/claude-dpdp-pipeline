#!/usr/bin/env python3
"""DPDP web app: the pipeline in a browser tab, at https://<server>/dpdp/ (behind Apache).

One page per user: the run list and review pages on the left, and on the right a terminal running the official Claude Code
under the user's own Linux account and their own Claude login (the same thing RStudio's Terminal tab gives, in an app layout).

- Login: Linux username + password (PAM service "dpdp-web"), members of DPDP_WEB_GROUP only (default rstudio-users).
- Terminal: one session per user, started as that user (needs root, like RStudio Server). It keeps running when the tab is
  closed; reopening the page reattaches and replays the recent output.
- Review pages: <DPDP_DATA>/<run>/review/*.html, served with a per-login token. "Send to Claude" writes
  <run>/review/pending_decision.json (same format as pipeline/review_server.py), which Claude picks up with wait_decision.py.
- Upload: report + sample sheet into <DPDP_DATA>/input/<name>/.

Dev mode (not root): only the account running the server can use it, and login is by a one-time link printed at start
(no password), e.g.  DPDP_DATA=... python3 web/server.py --dev
"""
import argparse, asyncio, fcntl, grp, hmac, html as htmlmod, json, os, pwd, re, secrets, shutil, signal, struct, subprocess, sys, termios, time
from aiohttp import web, WSMsgType

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.geteuid() == 0
STEPS = {"01_load", "03_parameters", "04_noise", "05_filter", "06_normalize", "06b_normalize_compare", "07_impute", "07b_impute_compare", "08_rollup",
         "09_qc", "10_stats"}
RUN_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
FILE_RE = re.compile(r"^[A-Za-z0-9_.-]+\.html$")
USER_RE = re.compile(r"^[a-z_][a-z0-9_.-]{0,31}$")
SESSION_HOURS = 12
SCROLLBACK = 400_000


def log(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- config
def config():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dev", action="store_true", help="not root: one-time login link for the current user, no password")
    a = ap.parse_args()
    c = dict(base=os.environ.get("DPDP_WEB_BASE", "/dpdp").rstrip("/"),
             host=os.environ.get("DPDP_WEB_HOST", "127.0.0.1"),
             port=int(os.environ.get("DPDP_WEB_PORT", "8770")),
             data=os.path.realpath(os.environ.get("DPDP_DATA", "/srv/dpdp/data")),
             app=os.environ.get("DPDP_APP", "/srv/dpdp/shiny_DPDP"),
             plugin=os.environ.get("DPDP_PLUGIN", "/srv/dpdp/dpdp-marketplace"),
             group=os.environ.get("DPDP_WEB_GROUP", "rstudio-users"),
             workdir=os.environ.get("DPDP_WEB_WORKDIR", "dpdp"),   # Claude's folder, inside each user's home
             dev=a.dev)
    if not ROOT and not c["dev"]:
        sys.exit("Run as root (it starts terminals as each user, like RStudio Server), or use --dev to test as yourself.")
    if ROOT and c["dev"]:
        sys.exit("--dev is only for testing as a normal user.")
    if c["dev"] and c["host"] not in ("127.0.0.1", "localhost"):
        sys.exit("--dev listens on 127.0.0.1 only.")
    os.makedirs(os.path.join(c["data"], "input"), exist_ok=True)
    return c


# ---------------------------------------------------------------- login
class Auth:
    def __init__(self, cfg):
        self.cfg, self.sessions, self.fails = cfg, {}, {}
        self.dev_token = secrets.token_urlsafe(24) if cfg["dev"] else None
        self.me = pwd.getpwuid(os.geteuid()).pw_name

    def allowed(self, user):
        if not USER_RE.match(user or ""):
            return False
        try:
            pw = pwd.getpwnam(user)
        except KeyError:
            return False
        if self.cfg["dev"]:
            return user == self.me
        if pw.pw_uid == 0:
            return False
        try:
            g = grp.getgrnam(self.cfg["group"])
        except KeyError:
            return False
        return user in g.gr_mem or pw.pw_gid == g.gr_gid

    def password_ok(self, user, password):
        import pam   # python-pam
        return pam.pam().authenticate(user, password, service="dpdp-web")

    def throttled(self, key):
        now = time.time()
        f = [t for t in self.fails.get(key, []) if now - t < 900]
        self.fails[key] = f
        return len(f) >= 8

    def new_session(self, user):
        tok = secrets.token_urlsafe(32)
        self.sessions[tok] = {"user": user, "csrf": secrets.token_urlsafe(24), "t": time.time()}
        return tok

    def get(self, request):
        s = self.sessions.get(request.cookies.get("dpdp_session", ""))
        if not s:
            return None
        if time.time() - s["t"] > SESSION_HOURS * 3600:
            self.sessions.pop(request.cookies.get("dpdp_session"), None)
            return None
        s["t"] = time.time()
        return s


# ---------------------------------------------------------------- terminals
class Term:
    """One Claude Code session per user: a PTY running web/session_start.sh as that user."""

    def __init__(self, user, cfg, on_exit):
        self.user, self.clients, self.buf, self.alive = user, set(), bytearray(), True
        pw = pwd.getpwnam(user)
        work = os.path.join(pw.pw_dir, cfg["workdir"])
        env = {"HOME": pw.pw_dir, "USER": user, "LOGNAME": user, "SHELL": pw.pw_shell or "/bin/bash",
               "PATH": os.path.join(pw.pw_dir, ".local/bin") + ":/usr/local/bin:/usr/bin:/bin",
               "TERM": "xterm-256color", "COLORTERM": "truecolor", "LANG": os.environ.get("LANG", "C.UTF-8"),
               "DPDP_DATA": cfg["data"], "DPDP_APP": cfg["app"], "DPDP_PLUGIN": cfg["plugin"], "DPDP_WEB": "1",
               "DPDP_WEB_WORK": work}
        self.master, slave = os.openpty()
        kw = {}
        if ROOT:
            groups = [g.gr_gid for g in grp.getgrall() if user in g.gr_mem]
            kw = dict(user=pw.pw_uid, group=pw.pw_gid, extra_groups=groups)
        # setsid -c: new session with the PTY as controlling terminal (so Ctrl-C and resizing work)
        self.proc = subprocess.Popen(["setsid", "-c", "/bin/bash", os.path.join(HERE, "session_start.sh")],
                                     stdin=slave, stdout=slave, stderr=slave, env=env, cwd=pw.pw_dir, umask=0o002,
                                     close_fds=True, **kw)
        os.close(slave)
        os.set_blocking(self.master, False)
        self.on_exit = on_exit
        asyncio.get_running_loop().add_reader(self.master, self._read)
        log("terminal started for", user, "pid", self.proc.pid)

    def _read(self):
        try:
            data = os.read(self.master, 65536)
        except BlockingIOError:
            return
        except OSError:
            data = b""
        if not data:
            return self._closed()
        self.buf += data
        if len(self.buf) > SCROLLBACK:
            del self.buf[:len(self.buf) - SCROLLBACK]
        for ws in list(self.clients):
            asyncio.ensure_future(self._send(ws, data))

    async def _send(self, ws, data):
        try:
            await ws.send_bytes(data)
        except Exception:
            self.clients.discard(ws)

    def _closed(self):
        if not self.alive:
            return
        self.alive = False
        loop = asyncio.get_running_loop()
        loop.remove_reader(self.master)
        try:
            os.close(self.master)
        except OSError:
            pass
        try:
            self.proc.wait(timeout=2)
        except Exception:
            pass
        log("terminal ended for", self.user)
        for ws in list(self.clients):
            asyncio.ensure_future(self._send(ws, b"\r\n\x1b[33m[Claude session ended. Click 'Restart Claude' to start a new one.]\x1b[0m\r\n"))
        self.on_exit(self)

    def write(self, data):
        if self.alive:
            try:
                os.write(self.master, data)
            except OSError:
                pass

    def resize(self, cols, rows):
        if self.alive and 0 < cols < 1000 and 0 < rows < 500:
            fcntl.ioctl(self.master, termios.TIOCSWINSZ, struct.pack("HHHH", rows, cols, 0, 0))

    def stop(self):
        if self.alive:
            for sig in (signal.SIGHUP, signal.SIGTERM):
                try:
                    os.killpg(self.proc.pid, sig)
                except Exception:
                    pass
                try:
                    self.proc.wait(timeout=1.5)
                    break
                except Exception:
                    pass
            self._closed()


class Terms:
    def __init__(self, cfg):
        self.cfg, self.by_user = cfg, {}

    def get(self, user):
        t = self.by_user.get(user)
        if t is None or not t.alive:
            t = Term(user, self.cfg, lambda x: self.by_user.pop(x.user, None) if self.by_user.get(x.user) is x else None)
            self.by_user[user] = t
        return t

    def restart(self, user, old):
        old.stop()
        new = self.get(user)
        if new is not old:
            new.clients |= old.clients; old.clients.clear()
        return new


# ---------------------------------------------------------------- app
def make_app(cfg):
    A, T, B = Auth(cfg), Terms(cfg), cfg["base"]
    static = os.path.join(HERE, "static")
    routes = web.RouteTableDef()

    def secure(request):
        return request.headers.get("X-Forwarded-Proto") == "https"

    def page(fname, **subs):
        html = open(os.path.join(static, fname), encoding="utf-8").read()
        for k, v in subs.items():
            html = html.replace("__%s__" % k.upper(), htmlmod.escape(v, quote=True))
        return web.Response(text=html, content_type="text/html", headers={"Cache-Control": "no-store"})

    def need(request, csrf=False):
        s = A.get(request)
        if not s:
            raise web.HTTPUnauthorized(text=json.dumps({"error": "not logged in"}), content_type="application/json")
        if csrf and not hmac.compare_digest(request.headers.get("X-DPDP-Token", ""), s["csrf"]):
            raise web.HTTPForbidden(text=json.dumps({"error": "bad token: reload the page"}), content_type="application/json")
        return s

    def run_dir(run):
        if not run or not RUN_RE.match(run) or run in (".", ".."):
            return None
        d = os.path.join(cfg["data"], run)
        return d if os.path.isdir(os.path.join(d, "review")) else None

    def chown_to(path, user):
        if ROOT:
            pw = pwd.getpwnam(user)
            try:
                os.chown(path, pw.pw_uid, -1)
            except OSError:
                pass

    # -- http to https (Apache sends X-Forwarded-Proto on the 443 site only)
    @web.middleware
    async def https_only(request, handler):
        if ROOT and not secure(request) and request.headers.get("X-Forwarded-For"):
            host = request.headers.get("Host", "").split(":")[0]
            raise web.HTTPFound("https://%s%s" % (host, request.rel_url))
        resp = await handler(request)
        resp.headers.setdefault("X-Frame-Options", "SAMEORIGIN")
        resp.headers.setdefault("X-Content-Type-Options", "nosniff")
        return resp

    @routes.get(B)
    async def to_slash(request):
        raise web.HTTPFound(B + "/")

    @routes.get(B + "/")
    async def home(request):
        s = A.get(request)
        if not s:
            return page("login.html", base=B, msg="")
        pw = pwd.getpwnam(s["user"])
        return page("app.html", base=B, user=s["user"], name=(pw.pw_gecos.split(",")[0] or s["user"]), csrf=s["csrf"])

    @routes.post(B + "/login")
    async def login(request):
        form = await request.post()
        user, password = str(form.get("username", "")).strip().lower(), str(form.get("password", ""))
        key = request.headers.get("X-Forwarded-For", request.remote or "?").split(",")[0].strip()
        if A.throttled(key) or A.throttled(user):
            return page("login.html", base=B, msg="Too many failed logins. Wait 15 minutes.")
        ok = A.allowed(user) and not cfg["dev"] and await asyncio.get_running_loop().run_in_executor(None, A.password_ok, user, password)
        if not ok:
            A.fails.setdefault(key, []).append(time.time()); A.fails.setdefault(user, []).append(time.time())
            await asyncio.sleep(1.5)
            log("login failed:", user, "from", key)
            return page("login.html", base=B, msg="Login failed. Use your server (RStudio) username and password. Ask the admin to add you to the '%s' group if needed." % cfg["group"])
        log("login:", user, "from", key)
        r = web.HTTPFound(B + "/")
        r.set_cookie("dpdp_session", A.new_session(user), path=B + "/", httponly=True, secure=secure(request), samesite="Strict",
                     max_age=SESSION_HOURS * 3600)
        raise r

    @routes.get(B + "/dev-login")
    async def dev_login(request):
        if not cfg["dev"] or not A.dev_token or not hmac.compare_digest(request.query.get("t", ""), A.dev_token):
            raise web.HTTPForbidden(text="no")
        r = web.HTTPFound(B + "/")
        r.set_cookie("dpdp_session", A.new_session(A.me), path=B + "/", httponly=True, samesite="Strict")
        raise r

    @routes.post(B + "/logout")
    async def logout(request):
        need(request, csrf=True)
        A.sessions.pop(request.cookies.get("dpdp_session", ""), None)
        r = web.json_response({"ok": True}); r.del_cookie("dpdp_session", path=B + "/")
        return r

    @routes.get(B + "/static/{f}")
    async def static_file(request):
        f = request.match_info["f"]
        if not re.match(r"^[A-Za-z0-9_.-]+\.(js|css|svg|png)$", f):
            raise web.HTTPNotFound()
        p = os.path.join(static, f)
        if not os.path.isfile(p):
            raise web.HTTPNotFound()
        return web.FileResponse(p, headers={"Cache-Control": "max-age=3600"})

    # -- runs and review pages
    @routes.get(B + "/api/runs")
    async def runs(request):
        need(request)
        out = []
        for r in sorted(os.listdir(cfg["data"])):
            rv = os.path.join(cfg["data"], r, "review")
            if not RUN_RE.match(r) or not os.path.isdir(rv):
                continue
            pages = []
            for f in sorted(os.listdir(rv)):
                m = re.match(r"^(\d\d[a-z]?)_([a-z_]+)\.html$", f)
                if m:
                    pages.append({"file": f, "step": f[:-5], "label": m.group(2).replace("_", " ").title(),
                                  "mtime": os.path.getmtime(os.path.join(rv, f))})
            if not pages:
                continue
            pend = os.path.join(rv, "pending_decision.json")
            try:
                owner = pwd.getpwuid(os.stat(os.path.join(cfg["data"], r)).st_uid).pw_name
            except KeyError:
                owner = "?"
            out.append({"run": r, "owner": owner, "pages": pages, "pending": os.path.exists(pend),
                        "mtime": max([p["mtime"] for p in pages] or [os.path.getmtime(rv)])})
        out.sort(key=lambda x: -x["mtime"])
        return web.json_response({"runs": out, "data": cfg["data"]})

    @routes.get(B + "/runs/{run}/review/{f}")
    async def review_page(request):
        s = need(request)
        d, f = run_dir(request.match_info["run"]), request.match_info["f"]
        if not d or not FILE_RE.match(f) or not os.path.isfile(os.path.join(d, "review", f)):
            raise web.HTTPNotFound(text="Page not found")
        html = open(os.path.join(d, "review", f), "rb").read().replace(b"__DPDP_TOKEN__", s["csrf"].encode())
        # pages built before the API paths were relative
        api = (B + "/runs/api/").encode()
        html = html.replace(b'fetch("/api/', b'fetch("' + api)
        # the analyst is the logged-in user
        pw = pwd.getpwnam(s["user"])
        who = json.dumps(pw.pw_gecos.split(",")[0] or s["user"]).encode()
        html = html.replace(b"</head>", b"<script>window.DPDP_WEB_ANALYST=" + who + b";</script></head>", 1)
        return web.Response(body=html, content_type="text/html", headers={"Cache-Control": "no-store"})

    @routes.get(B + "/runs/api/status")
    async def status(request):
        need(request)
        d, step = run_dir(request.query.get("run", "")), request.query.get("step", "")
        if not d or step not in STEPS:
            return web.json_response({"error": "bad run or step"}, status=400)
        pend, pending = os.path.join(d, "review", "pending_decision.json"), None
        if os.path.exists(pend):
            try:
                pending = json.load(open(pend)).get("step")
            except Exception:
                pending = "?"
        pg = os.path.join(d, "review", step + ".html")
        return web.json_response({"pending_step": pending, "page_mtime": os.path.getmtime(pg) if os.path.exists(pg) else 0})

    @routes.post(B + "/runs/api/decision")
    async def decision(request):
        s = need(request, csrf=True)
        if (request.content_length or 0) > 20000:
            return web.json_response({"error": "bad body size"}, status=400)
        try:
            dec = await request.json()
        except Exception:
            return web.json_response({"error": "not JSON"}, status=400)
        d, step = run_dir(str(dec.get("run", ""))), str(dec.get("step", ""))
        if not d or step not in STEPS or not isinstance(dec.get("choice"), str):
            return web.json_response({"error": "bad run, step or choice"}, status=400)
        pw = pwd.getpwnam(s["user"])
        rec = {"run": os.path.basename(d), "step": step, "choice": dec["choice"][:200],
               "fields": {str(k)[:60]: str(v)[:500] for k, v in (dec.get("fields") or {}).items()},
               "note": str(dec.get("note", ""))[:500], "analyst": (pw.pw_gecos.split(",")[0] or s["user"])[:100],
               "user": s["user"], "submitted": time.strftime("%Y-%m-%d %H:%M:%S"), "source": "web app"}
        path = os.path.join(d, "review", "pending_decision.json")
        with open(path + ".tmp", "w") as fh:
            json.dump(rec, fh, indent=1)
        os.chmod(path + ".tmp", 0o664); chown_to(path + ".tmp", s["user"])
        os.replace(path + ".tmp", path)
        log("decision", s["user"], rec["run"], step, rec["choice"])
        return web.json_response({"ok": True, "submitted": rec["submitted"]})

    # -- upload: report + sample sheet into <data>/input/<name>/
    @routes.post(B + "/api/upload")
    async def upload(request):
        s = need(request, csrf=True)
        name = request.query.get("name", "")
        if not RUN_RE.match(name) or name in (".", ".."):
            return web.json_response({"error": "Name: letters, digits, _ . - only"}, status=400)
        dest = os.path.join(cfg["data"], "input", name)
        os.makedirs(dest, exist_ok=True); os.chmod(dest, 0o2775); chown_to(dest, s["user"])
        saved, reader = [], await request.multipart()
        while True:
            part = await reader.next()
            if part is None:
                break
            fn = re.sub(r"[^A-Za-z0-9_.+-]", "_", os.path.basename(part.filename or "")).lstrip(".")
            if not re.match(r"^.+\.(tsv|txt|csv|xlsx|xls|zip)$", fn, re.I):
                return web.json_response({"error": "File type not allowed: %s (tsv, txt, csv, xlsx, xls, zip)" % fn}, status=400)
            tmp, n = os.path.join(dest, "." + fn + ".part"), 0
            with open(tmp, "wb") as fh:
                while True:
                    chunk = await part.read_chunk(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk); n += len(chunk)
            os.chmod(tmp, 0o664); chown_to(tmp, s["user"])
            os.replace(tmp, os.path.join(dest, fn))
            saved.append({"file": fn, "bytes": n})
        log("upload", s["user"], dest, [x["file"] for x in saved])
        return web.json_response({"ok": True, "folder": dest, "files": saved})

    # -- terminal
    @routes.get(B + "/ws/term")
    async def term_ws(request):
        s = need(request)
        origin = request.headers.get("Origin", "")
        host = request.headers.get("Host", "")
        if not origin or origin.split("://", 1)[-1] != host:
            raise web.HTTPForbidden(text="bad origin")
        ws = web.WebSocketResponse(heartbeat=30, max_msg_size=1 << 20)
        await ws.prepare(request)
        t = T.get(s["user"])
        if t.buf:
            await ws.send_bytes(bytes(t.buf))
        t.clients.add(ws)
        try:
            async for msg in ws:
                if msg.type != WSMsgType.TEXT:
                    continue
                try:
                    m = json.loads(msg.data)
                except Exception:
                    continue
                t = T.by_user.get(s["user"]) or t
                if m.get("t") == "i":
                    t.write(str(m.get("d", "")).encode())
                elif m.get("t") == "r":
                    t.resize(int(m.get("c", 0)), int(m.get("r", 0)))
                elif m.get("t") == "restart":
                    log("restart requested by", s["user"])
                    t = T.restart(s["user"], t)
                    for c in list(t.clients):
                        asyncio.ensure_future(t._send(c, b"\x1bc"))
                    t.resize(int(m.get("c", 0)) or 120, int(m.get("r", 0)) or 40)
        finally:
            for x in list(T.by_user.values()):
                x.clients.discard(ws)
        return ws

    app = web.Application(middlewares=[https_only], client_max_size=1 << 20)
    app.add_routes(routes)
    app["auth"] = A
    return app


def main():
    cfg = config()
    app = make_app(cfg)
    if cfg["dev"]:
        print("DEV login link (this account only, valid while the server runs):\n  http://127.0.0.1:%d%s/dev-login?t=%s"
              % (cfg["port"], cfg["base"], app["auth"].dev_token), flush=True)
    log("DPDP web app on %s:%d%s/  data %s  group %s" % (cfg["host"], cfg["port"], cfg["base"], cfg["data"], cfg["group"]))
    web.run_app(app, host=cfg["host"], port=cfg["port"], print=None, access_log=None)


if __name__ == "__main__":
    main()
