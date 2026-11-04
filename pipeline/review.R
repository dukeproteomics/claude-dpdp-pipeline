# Review pages: one self-contained HTML file per step (graphs embedded, interactive parts inline),
# plus an index of all steps and recorded decisions. Written to <run_dir>/review/.
# The pages are read-only: the analyst decides in the Claude session, and the decision is logged to decisions.jsonl.

suppressPackageStartupMessages({ library(jsonlite); library(htmltools) })

review_dir <- function(run_dir) { d <- file.path(run_dir, "review"); dir.create(d, showWarnings = FALSE); d }

esc <- function(x) htmltools::htmlEscape(as.character(x))
fmt_n <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
fmt_pct <- function(x, d = 2) paste0(formatC(as.numeric(x), format = "f", digits = d), "%")

embed_png <- function(path, alt, caption = NULL) {
  path <- gsub("/+", "/", path)
  if (!file.exists(path)) return(sprintf('<p class="missing">Plot not found: %s</p>', esc(path)))
  sprintf('<figure><img alt="%s" src="data:image/png;base64,%s"/><figcaption>%s<span class="file">%s</span></figcaption></figure>',
          esc(alt), base64enc::base64encode(path), if (is.null(caption)) "" else paste0(esc(caption), "<br/>"), esc(path))
}

html_table <- function(df, numeric_cols = integer(0), class = "") {
  head <- paste0("<tr>", paste0("<th", ifelse(seq_along(df) %in% numeric_cols, ' class="num"', ""), ">", esc(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td", ifelse(seq_along(df) %in% numeric_cols, ' class="num"', ""), ">", vapply(df[i, ], function(v) as.character(v), ""), "</td>", collapse = ""), "</tr>"), "")
  sprintf('<table class="%s"><thead>%s</thead><tbody>%s</tbody></table>', class, head, paste(rows, collapse = ""))
}

read_decisions <- function(run_dir, step = NULL) {
  f <- file.path(run_dir, "decisions.jsonl")
  if (!file.exists(f)) return(list())
  recs <- lapply(readLines(f, warn = FALSE), function(l) tryCatch(fromJSON(l, simplifyVector = FALSE), error = function(e) NULL))
  recs <- Filter(Negate(is.null), recs)
  if (!is.null(step)) recs <- Filter(function(r) identical(r$step, step), recs)
  recs
}

plotly_js <- function() {
  f <- list.files(system.file("htmlwidgets/lib/plotlyjs", package = "plotly"), pattern = "\\.min\\.js$", full.names = TRUE)[1]
  if (is.na(f)) stop("plotly.js not found in the plotly R package")
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

REVIEW_CSS <- "
.viz-root { color-scheme: light; --surface-0:#f4f4f2; --surface-1:#fcfcfb; --text-primary:#0b0b0b; --text-secondary:#52514e;
  --text-muted:#6f6e69; --rule:#d9d8d3; --grid:#ebeae6; --series-1:#2a78d6; --ref-dynamic:#eb6834; --ref-fixed:#008300;
  --accent:#2a78d6; --warn-bg:#fff4d6; --warn-ink:#6b4a00; --ok-bg:#e3f4ea; --ok-ink:#0f5132; }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme=\"light\"])) .viz-root { color-scheme: dark; --surface-0:#111110;
  --surface-1:#1a1a19; --text-primary:#f0efec; --text-secondary:#c3c2b7; --text-muted:#9b9a92; --rule:#3a3a37; --grid:#2a2a28;
  --series-1:#3987e5; --ref-dynamic:#d95926; --ref-fixed:#2f9e2f; --accent:#3987e5; --warn-bg:#3a2e0c; --warn-ink:#f5d58a;
  --ok-bg:#12301f; --ok-ink:#a8e0bf; } }
body { margin:0; font: 15px/1.5 system-ui, -apple-system, 'Segoe UI', sans-serif; }
.viz-root { background:var(--surface-0); color:var(--text-primary); min-height:100vh; }
header { background:var(--surface-1); border-bottom:1px solid var(--rule); padding:14px 28px; display:flex; gap:24px; align-items:baseline; flex-wrap:wrap; }
header h1 { font-size:20px; margin:0; } header .meta { color:var(--text-secondary); font-size:13px; }
header a { color:var(--accent); }
main { max-width:1180px; margin:0 auto; padding:20px 28px 60px; }
section { background:var(--surface-1); border:1px solid var(--rule); border-radius:8px; padding:16px 20px; margin:16px 0; }
section h2 { font-size:16px; margin:0 0 10px; }
.status { display:inline-block; padding:2px 10px; border-radius:999px; font-size:13px; font-weight:600; }
.status.waiting { background:var(--warn-bg); color:var(--warn-ink); } .status.done { background:var(--ok-bg); color:var(--ok-ink); }
.tiles { display:flex; gap:12px; flex-wrap:wrap; }
.tile { border:1px solid var(--rule); border-radius:8px; padding:10px 14px; min-width:150px; }
.tile .v { font-size:22px; font-weight:600; font-variant-numeric: tabular-nums; } .tile .l { color:var(--text-secondary); font-size:13px; }
table { border-collapse:collapse; font-size:14px; margin:6px 0; max-width:100%; } td { overflow-wrap:anywhere; } td:nth-child(-n+3) { white-space:nowrap; overflow-wrap:normal; } th, td { padding:4px 10px; border-bottom:1px solid var(--grid); text-align:left; }
th { color:var(--text-secondary); font-weight:600; } td.num, th.num { text-align:right; font-variant-numeric: tabular-nums; }
figure { margin:8px 0 16px; } figure img { max-width:min(100%, 760px); border:1px solid var(--rule); border-radius:4px; background:#fff; }
figcaption { color:var(--text-secondary); font-size:13px; } .file { color:var(--text-muted); font-size:12px; }
.grid2 { display:grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap:16px; }
.options li { margin:6px 0; } code { background:var(--grid); padding:1px 5px; border-radius:4px; font-size:13px; }
pre { background:var(--grid); padding:10px; border-radius:6px; overflow-x:auto; font-size:13px; }
.note { color:var(--text-secondary); font-size:13px; } .missing { color:var(--warn-ink); }
.control { display:flex; gap:14px; align-items:center; flex-wrap:wrap; margin:6px 0 10px; }
.control input[type=range] { width:420px; } .control .val { font-weight:600; font-variant-numeric: tabular-nums; min-width:200px; }
.nextbox { border:2px solid var(--warn-ink); background:var(--warn-bg); } .nextbox h2 { color:var(--warn-ink); }
.dchoice { margin:6px 0; } .dfield { margin-left:18px; color:var(--text-secondary); font-size:14px; } .dfield input, .dfield select, .dbar input { font:inherit; font-size:14px; padding:2px 6px; }
.dbar { margin-top:10px; display:flex; gap:8px; align-items:center; flex-wrap:wrap; } .dsend { font:inherit; font-weight:600; padding:6px 16px; border-radius:6px; border:1px solid var(--accent); background:var(--accent); color:#fff; cursor:pointer; } .dsend:disabled { opacity:0.5; cursor:default; }
.dmsg { font-weight:600; color:var(--text-primary); }
.nextbox.done { border-color:var(--ok-ink); background:var(--ok-bg); } .nextbox.done h2 { color:var(--ok-ink); }
.control button { font:inherit; font-size:13px; padding:3px 10px; border:1px solid var(--rule); border-radius:6px; background:var(--surface-1); color:var(--text-primary); cursor:pointer; }
"

# What the analyst does next, shown at the top of every page: the decision form (sent to Claude through pipeline/review_server.py)
# while the step waits, or what comes next when it is done. Choices match the AskUserQuestion options in each step's skill.
ch <- function(value, label, fields = list()) list(value = value, label = label, fields = fields)
fld <- function(name, label, type = "text", default = "", options = NULL) list(name = name, label = label, type = type, default = default, options = options)

decision_spec <- function(step, variant, p) {
  g <- function(k, d = "") { v <- p[[k]]; if (is.null(v) || length(v) == 0 || is.na(v[1])) d else as.character(v[1]) }
  switch(paste0(step, if (nzchar(variant)) paste0(":", variant) else ""),
    "01_load" = list(q = "Does the design match the data?", nxt = "Step 3 - Parameters", choices = list(
      ch("accept", "Yes, continue"), ch("fix", "No, fix the sample sheet", list(fld("problem", "What is wrong?"))))),
    "03_parameters" = list(q = "Keep these parameters?", nxt = "Step 4 - Noise", choices = list(
      ch("accept", "Yes, continue"), ch("change", "Change settings", list(fld("change", "Which setting and value?"))))),
    "04_noise" = list(q = "Which noise setting?", nxt = "Step 5 - Filter", choices = list(
      ch("none", "none - keep every value"),
      ch("dynamic", paste0("dynamic - remove values below the inflection (", g("noise_inflection"), ")"), list(fld("noise_keep", "noise_keep (put values back for precursors that pass the filter)", "checkbox", "false"))),
      ch("fixed", "fixed - remove values below a baseline", list(fld("noise_baseline_value", "baseline intensity", "number", g("noise_baseline_value", "1")),
                                                              fld("noise_keep", "noise_keep (put values back for precursors that pass the filter)", "checkbox", "false"))))),
    "05_filter" = list(q = "Apply these filter settings?", nxt = "Step 6 - Normalize", choices = list(
      ch("apply", "Apply as shown"), ch("change", "Change settings", list(fld("change", "Which filter and value? (e.g. misaligned_cutoff=60)"))))),
    "06_normalize" = list(q = "Which proteins to leave out of the factors, and which method?", nxt = "Step 7 - Impute",
      global = list(fld("norm_type", "Method", "select", "sltmm", c("sltmm", "sl", "tmm", "quantile", "loess", "vsn", "lr"))),
      choices = list(ch("default_exclude", "Exclude: pipeline default (added material only)"), ch("current_exclude", "Exclude: keep the current rule"),
                     ch("custom_exclude", "Exclude: my own regex", list(fld("norm_exclude_grep", "description regex", "text", g("norm_exclude_grep")))))),
    "07_impute" = list(q = "How should the data be imputed?", nxt = "Step 8 - Rollup", choices = list(
      ch("duke", "duke (near the group mean), bottom-X for the rest", list(fld("missing_cutoff", "missing_cutoff %", "number", g("missing_cutoff", "50")), fld("bottom_x", "bottom_x %", "number", g("bottom_x", "2")))),
      ch("bottomx", "bottomx only", list(fld("bottom_x", "bottom_x %", "number", g("bottom_x", "2")))))),
    "08_rollup" = list(q = "Which rollup?", nxt = "Step 9 - QC", choices = list(
      ch("iq_maxlfq_scaled", "iq_maxlfq, scaled"), ch("iq_maxlfq", "iq_maxlfq, not scaled"), ch("sum", "sum"))),
    "09_qc" = list(q = "Accept the QC and continue to Stats?", nxt = "Step 10 - Stats", choices = list(
      ch("accept_flag", "Accept and flag samples", list(fld("flag_samples", "samples to flag (comma separated)", "text", ""))),
      ch("accept", "Accept"), ch("stop", "Stop here"))),
    "10_stats" = list(q = "Run stats with these settings?", nxt = "Stats results", choices = list(
      ch("run", "Run as shown"), ch("run_adjp", "Run, with Up/Down calls by adjusted p"),
      ch("change", "Change settings or comparisons", list(fld("change", "What to change?"))))),
    "10_stats:results" = list(q = "Accept the stats results?", nxt = "Finished (Pathway and Report are not built yet)", choices = list(
      ch("accept", "Accept"), ch("change", "Change cutoffs and re-run", list(fld("change", "New p type / p cutoff / FC cutoff"))),
      ch("plots", "Add protein plots", list(fld("proteins", "accessions (comma separated)"))))),
    NULL)
}

next_step_box <- function(run_dir, step, status_class, variant = "") {
  db <- file.path(run_dir, paste0(basename(run_dir), ".db"))
  if (!file.exists(db)) db <- list.files(run_dir, pattern = "\\.db$", full.names = TRUE)[1]
  p <- if (!is.na(db) && file.exists(db)) tryCatch(get_params(db), error = function(e) list()) else list()
  spec <- decision_spec(step, variant, p)
  if (is.null(spec)) return("")
  if (!identical(status_class, "waiting"))
    return(sprintf('<section class="nextbox done"><h2>This step is done</h2><p>Next: <b>%s</b>. Claude continues in the session.</p></section>', esc(spec$nxt)))
  sprintf('<section class="nextbox" id="decision"><h2>Your decision: %s</h2><div id="dform"></div>
<p class="note" id="dhelp">Pick an option and press <b>Send to Claude</b>, or answer in the Claude window instead.</p></section>
<script>window.DPDP_DECISION = %s;</script>', esc(spec$q),
    toJSON(list(run = basename(run_dir), step = step, spec = spec), auto_unbox = TRUE, null = "null"))
}

DECISION_JS <- '
(function () {
  const S = window.DPDP_DECISION; if (!S) return;
  const box = document.getElementById("dform"), help = document.getElementById("dhelp");
  const token = document.querySelector("meta[name=dpdp-token]").content, served = location.protocol.startsWith("http") && token !== "__DPDP_" + "TOKEN__";
  const el = (t, a = {}, h = "") => { const e = document.createElement(t); Object.entries(a).forEach(([k, v]) => e.setAttribute(k, v)); e.innerHTML = h; return e; };
  const field = (f, pre) => { const w = el("label", { class: "dfield" }, f.label + " ");
    let inp;
    if (f.type === "select") { inp = el("select"); (f.options || []).forEach(o => inp.appendChild(el("option", o === f.default ? { selected: "" } : {}, o))); }
    else if (f.type === "checkbox") { inp = el("input", { type: "checkbox" }); inp.checked = f.default === "true"; }
    else { inp = el("input", { type: f.type === "number" ? "number" : "text", step: "any", size: "40" }); inp.value = f.default || ""; }
    inp.dataset.name = f.name; inp.dataset.pre = pre; w.appendChild(inp); return w; };
  S.spec.choices.forEach((c, i) => {
    const row = el("div", { class: "dchoice" }), lab = el("label");
    const r = el("input", { type: "radio", name: "dchoice", value: c.value }); if (i === 0) r.checked = true;
    lab.appendChild(r); lab.appendChild(document.createTextNode(" " + c.label)); row.appendChild(lab);
    (c.fields || []).forEach(f => row.appendChild(field(f, c.value))); box.appendChild(row); });
  (S.spec.global || []).forEach(f => box.appendChild(field(f, "")));
  const note = el("input", { type: "text", size: "60", placeholder: "optional: why (goes in the decision log)" });
  const who = el("input", { type: "text", size: "20", placeholder: "your name" }); who.value = localStorage.getItem("dpdp_analyst") || "";
  const bar = el("div", { class: "dbar" }); bar.appendChild(el("span", {}, "Note ")); bar.appendChild(note); bar.appendChild(el("span", {}, " Analyst ")); bar.appendChild(who);
  const btn = el("button", { type: "button", class: "dsend" }, "Send to Claude"); bar.appendChild(btn); box.appendChild(bar);
  const msg = el("p", { class: "dmsg" }); box.appendChild(msg);
  window.dpdpSetField = (name, value) => { box.querySelectorAll("[data-name=\\"" + name + "\\"]").forEach(i => { if (!i.dataset.edited) i.value = value; }); };
  Object.entries(window.DPDP_PREFILL || {}).forEach(([k, v]) => window.dpdpSetField(k, v));
  box.addEventListener("input", e => { if (e.target.dataset) e.target.dataset.edited = "1"; });
  if (!served) { btn.disabled = true; msg.innerHTML = "Buttons need the review server. Ask Claude to start it, then open <code>http://127.0.0.1:8765/" + S.run + "/review/" + S.step + ".html</code>. Or answer in the Claude window."; return; }
  let mtime0 = null;
  const poll = () => fetch("/api/status?run=" + encodeURIComponent(S.run) + "&step=" + S.step).then(r => r.json()).then(st => {
    if (mtime0 === null) mtime0 = st.page_mtime; else if (st.page_mtime > mtime0) { location.reload(); return; }
    if (st.pending_step === S.step && !btn.disabled) { btn.disabled = true; msg.textContent = "A decision for this step is waiting for Claude."; }
  }).catch(() => {});
  poll(); setInterval(poll, 2000);
  btn.addEventListener("click", () => {
    const choice = box.querySelector("input[name=dchoice]:checked").value, fields = {};
    box.querySelectorAll("[data-name]").forEach(i => { if (i.dataset.pre === "" || i.dataset.pre === choice) fields[i.dataset.name] = i.type === "checkbox" ? String(i.checked) : i.value; });
    localStorage.setItem("dpdp_analyst", who.value);
    btn.disabled = true; msg.textContent = "Sending...";
    fetch("/api/decision", { method: "POST", headers: { "Content-Type": "application/json", "X-DPDP-Token": token },
      body: JSON.stringify({ run: S.run, step: S.step, choice: choice, fields: fields, note: note.value, analyst: who.value }) })
      .then(r => r.json().then(j => ({ ok: r.ok, j }))).then(({ ok, j }) => {
        if (!ok) { btn.disabled = false; msg.textContent = "Not sent: " + (j.error || "error"); return; }
        msg.textContent = "Sent to Claude at " + j.submitted.slice(11, 16) + ". This page reloads when Claude has applied it.";
        help.textContent = ""; })
      .catch(() => { btn.disabled = false; msg.textContent = "Not sent: is the review server running?"; });
  });
})();
'

write_review_page <- function(run_dir, step, title, status, status_class, body, extra_head = "", decision_variant = "") {
  body <- paste0(next_step_box(run_dir, step, status_class, decision_variant), body, '<script>', DECISION_JS, '</script>')
  extra_head <- paste0('<meta name="dpdp-token" content="__DPDP_TOKEN__"/>', extra_head)
  p_prefix <- basename(run_dir)
  html <- sprintf('<!doctype html><html lang="en"><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>%s - %s</title><style>%s</style>%s</head><body><div class="viz-root">
<header><h1>%s</h1><span class="status %s">%s</span><span class="meta">Run <b>%s</b> &middot; page made %s &middot; <a href="index.html">all steps</a></span></header>
<main>%s</main></div></body></html>',
    esc(title), esc(p_prefix), REVIEW_CSS, extra_head, esc(title), status_class, esc(status), esc(p_prefix),
    format(Sys.time(), "%Y-%m-%d %H:%M"), body)
  f <- file.path(review_dir(run_dir), paste0(step, ".html"))
  writeLines(html, f, useBytes = TRUE)
  write_review_index(run_dir)
  f
}

# index.html: every step page, plus every decision recorded in decisions.jsonl.
write_review_index <- function(run_dir) {
  d <- review_dir(run_dir)
  pages <- sort(setdiff(list.files(d, pattern = "\\.html$"), "index.html"))
  dec <- read_decisions(run_dir)
  dec_rows <- if (length(dec)) data.frame(
    Time = vapply(dec, function(r) r$time %||% "", ""), Step = vapply(dec, function(r) r$step %||% "", ""),
    Analyst = vapply(dec, function(r) r$analyst %||% "", ""),
    Settings = vapply(dec, function(r) esc(paste(names(r$settings), vapply(r$settings, function(v) paste(unlist(v), collapse = ","), ""), sep = "=", collapse = "; ")), ""),
    Note = vapply(dec, function(r) esc(r$note %||% ""), ""), stringsAsFactors = FALSE) else NULL
  links <- if (length(pages)) paste0("<ul>", paste0('<li><a href="', pages, '">', esc(sub("\\.html$", "", pages)), "</a></li>", collapse = ""), "</ul>") else "<p>No step pages yet.</p>"
  body <- paste0('<section><h2>Step review pages</h2>', links, '</section><section><h2>Decisions recorded</h2>',
                 if (is.null(dec_rows)) "<p>None yet.</p>" else html_table(dec_rows), '</section>')
  html <- sprintf('<!doctype html><html lang="en"><head><meta charset="utf-8"/><title>%s - review</title><style>%s</style></head><body><div class="viz-root">
<header><h1>Run %s</h1><span class="meta">updated %s</span></header><main>%s</main></div></body></html>',
    esc(basename(run_dir)), REVIEW_CSS, esc(basename(run_dir)), format(Sys.time(), "%Y-%m-%d %H:%M"), body)
  writeLines(html, file.path(d, "index.html"), useBytes = TRUE)
}
