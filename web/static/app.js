(function () {
  "use strict";
  const meta = n => document.querySelector("meta[name=" + n + "]").content;
  const BASE = meta("dpdp-base"), CSRF = meta("dpdp-csrf");
  const $ = id => document.getElementById(id);

  // ------------------------------------------------------------ terminal (official Claude Code, as this user)
  const term = new Terminal({ cursorBlink: true, fontSize: 14, scrollback: 5000, allowProposedApi: true, lineHeight: 1.15,
    fontFamily: '"DejaVu Sans Mono", Menlo, Consolas, monospace', theme: { background: "#1e1e1e" } });
  const fit = new FitAddon.FitAddon();
  term.loadAddon(fit);
  term.loadAddon(new WebLinksAddon.WebLinksAddon((e, uri) => window.open(uri, "_blank", "noopener")));
  term.open($("term"));
  let ws = null, retry = 0;
  const send = m => { if (ws && ws.readyState === 1) ws.send(JSON.stringify(m)); };
  const size = () => { try { fit.fit(); } catch (e) {} send({ t: "r", c: term.cols, r: term.rows }); };
  function connect() {
    const proto = location.protocol === "https:" ? "wss://" : "ws://";
    ws = new WebSocket(proto + location.host + BASE + "/ws/term");
    ws.binaryType = "arraybuffer";
    $("tstat").textContent = "connecting...";
    ws.onopen = () => { retry = 0; $("tstat").textContent = "connected"; term.reset(); size(); setTimeout(size, 300); };
    ws.onmessage = e => term.write(new Uint8Array(e.data));
    ws.onclose = () => {
      $("tstat").textContent = "disconnected, reconnecting...";
      fetch(BASE + "/api/runs").then(r => { if (r.status === 401) location.reload(); }).catch(() => {});
      setTimeout(connect, Math.min(10000, 1000 * Math.pow(2, retry++)));
    };
  }
  term.onData(d => send({ t: "i", d: d }));
  window.addEventListener("resize", size);
  connect();
  $("restart").onclick = () => {
    if (confirm("End this Claude session and start a new one? (Claude can continue the last conversation.)")) send({ t: "restart", c: term.cols, r: term.rows });
    term.focus();
  };
  const typeIntoClaude = text => { send({ t: "i", d: text }); term.focus(); };

  // ------------------------------------------------------------ runs and review pages
  let runs = [], cur = localStorage.getItem("dpdp_run") || "", curFile = "", seen = {}, firstLoad = true;
  const sel = $("run"), frame = $("page");
  function showPage(run, file) {
    if (!run || !file) { frame.hidden = true; $("empty").hidden = false; curFile = ""; return; }
    $("empty").hidden = true; frame.hidden = false;
    if (run !== cur || file !== curFile) frame.src = BASE + "/runs/" + encodeURIComponent(run) + "/review/" + encodeURIComponent(file);
    cur = run; curFile = file; localStorage.setItem("dpdp_run", run);
    drawSteps();
  }
  function drawSteps() {
    const r = runs.find(x => x.run === cur), nav = $("steps");
    nav.innerHTML = "";
    if (!r) return;
    r.pages.forEach(p => {
      const a = document.createElement("a");
      a.href = "#"; a.textContent = p.step.split("_")[0] + " " + p.label.replace(/\bQc\b/, "QC");
      if (p.file === curFile) a.className = "on";
      a.onclick = e => { e.preventDefault(); $("follow").checked = false; showPage(cur, p.file); };
      nav.appendChild(a);
    });
    $("pending").hidden = !r.pending;
  }
  const newest = r => r.pages.reduce((a, p) => (!a || p.mtime > a.mtime ? p : a), null);
  function refresh() {
    fetch(BASE + "/api/runs").then(r => { if (r.status === 401) { location.reload(); throw 0; } return r.json(); }).then(j => {
      const before = Object.assign({}, seen);
      runs = j.runs;
      const opts = runs.map(r => r.run + "|" + r.owner).join(",");
      if (sel.dataset.opts !== opts) {
        sel.innerHTML = "";
        runs.forEach(r => { const o = document.createElement("option"); o.value = r.run; o.textContent = r.run + "  (" + r.owner + ")"; sel.appendChild(o); });
        sel.dataset.opts = opts;
      }
      if (!runs.length) { showPage("", ""); return; }
      if (!runs.find(r => r.run === cur)) cur = runs[0].run;
      // follow: a run that changed since the last poll (Claude is working on it) and its newest page
      if ($("follow").checked && !firstLoad) {
        const moved = runs.filter(r => before[r.run] !== undefined ? r.mtime > before[r.run] : true);
        if (moved.length) { const r = moved[0], p = newest(r); if (p) { showPage(r.run, p.file); } }
      }
      runs.forEach(r => { seen[r.run] = r.mtime; });
      sel.value = cur;
      const r = runs.find(x => x.run === cur);
      if (firstLoad || !curFile) { const p = r && newest(r); showPage(cur, p ? p.file : ""); }
      firstLoad = false;
      drawSteps();
    }).catch(() => {});
  }
  sel.onchange = () => { const r = runs.find(x => x.run === sel.value), p = r && newest(r); showPage(sel.value, p ? p.file : ""); };
  refresh(); setInterval(refresh, 4000);

  // ------------------------------------------------------------ new project upload
  const dlg = $("upload");
  $("newproj").onclick = () => { $("upmsg").textContent = ""; $("upprog").hidden = true; $("upgo").disabled = false; dlg.showModal(); };
  $("upcancel").onclick = () => dlg.close();
  $("upform").onsubmit = e => {
    e.preventDefault();
    const name = $("upname").value.trim(), fd = new FormData();
    if (!/^[A-Za-z0-9_.-]+$/.test(name)) { $("upmsg").textContent = "Name: letters, digits, _ . - only"; return; }
    fd.append("report", $("uprep").files[0]); fd.append("samples", $("upss").files[0]);
    const x = new XMLHttpRequest();
    x.open("POST", BASE + "/api/upload?name=" + encodeURIComponent(name));
    x.setRequestHeader("X-DPDP-Token", CSRF);
    $("upprog").hidden = false; $("upgo").disabled = true; $("upmsg").textContent = "Uploading...";
    x.upload.onprogress = ev => { if (ev.lengthComputable) $("upprog").value = 100 * ev.loaded / ev.total; };
    x.onload = () => {
      let j = {}; try { j = JSON.parse(x.responseText); } catch (err) {}
      if (x.status !== 200) { $("upmsg").textContent = "Upload failed: " + (j.error || x.status); $("upgo").disabled = false; return; }
      const f = j.files.map(f => j.folder + "/" + f.file);
      dlg.close();
      typeIntoClaude("Process project " + name + " from " + f[0] + " and " + f[1]);
    };
    x.onerror = () => { $("upmsg").textContent = "Upload failed (connection)."; $("upgo").disabled = false; };
    x.send(fd);
  };

  // ------------------------------------------------------------ layout and logout
  const split = $("split"), left = $("left");
  split.onmousedown = e => {
    e.preventDefault(); split.classList.add("drag"); frame.style.pointerEvents = "none";
    const move = ev => { left.style.flexBasis = Math.max(280, Math.min(window.innerWidth - 340, ev.clientX)) + "px"; };
    const up = () => { split.classList.remove("drag"); frame.style.pointerEvents = ""; document.removeEventListener("mousemove", move);
      document.removeEventListener("mouseup", up); localStorage.setItem("dpdp_split", left.style.flexBasis); size(); };
    document.addEventListener("mousemove", move); document.addEventListener("mouseup", up);
  };
  if (localStorage.getItem("dpdp_split")) left.style.flexBasis = localStorage.getItem("dpdp_split");
  $("logout").onclick = () => fetch(BASE + "/logout", { method: "POST", headers: { "X-DPDP-Token": CSRF } }).then(() => location.reload());
})();
