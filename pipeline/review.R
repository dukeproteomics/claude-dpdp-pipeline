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

# The app's static plots, kept for reference behind a collapsed panel (the interactive charts above are the ones to read).
app_png_panel <- function(files, captions = NULL, title = "The app&apos;s original plots") {
  files <- files[file.exists(gsub("/+", "/", files))]
  if (!length(files)) return("")
  if (is.null(captions)) captions <- sub("\\.png$", "", basename(files))
  paste0('<details><summary>', title, '</summary><div class="grid2">',
         paste(mapply(function(f, cp) embed_png(f, basename(f), cp), files, captions), collapse = ""), '</div></details>')
}

# Per-sample bar chart, coloured by group (fixed group order = sample_groups order). `series`: named list of numeric vectors
# (one per sample); with more than one series a select switches between them (one y-axis at a time, never two).
sample_bar <- function(id, labels, groups, series, ytitle = names(series), group_order = unique(groups), height = 300) {
  if (length(ytitle) == 1 && length(series) > 1) ytitle <- rep(ytitle, length(series))
  d <- list(samples = labels, groups = groups, order = group_order, series = lapply(series, function(v) unname(as.numeric(v))),
            ytitle = as.list(setNames(ytitle, names(series))))
  sprintf('%s<div id="%s_sb" style="height:%dpx" role="img" aria-label="Per-sample bar chart"></div>
<details><summary>Table view</summary><div id="%s_sbt"></div></details>
<script>window.DPDP_BARS = (window.DPDP_BARS || []).concat([{ id: "%s", d: %s }]);</script>',
    if (length(series) > 1) sprintf('<div class="control"><label for="%s_sbs">Show</label><select id="%s_sbs">%s</select></div>', id, id,
      paste0('<option>', esc(names(series)), '</option>', collapse = "")) else "",
    id, height, id, id, toJSON(d, auto_unbox = TRUE, digits = NA))
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
  --text-muted:#6f6e69; --rule:#d9d8d3; --grid:#ebeae6; --series-1:#2a78d6; --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100;
  --series-5:#e87ba4; --series-6:#008300; --series-7:#6250d6; --series-8:#e34948; --div-neg:#2a78d6; --div-mid:#f0efec; --div-pos:#e34948;
  --seq-lo:#cde2fb; --seq-hi:#104281; --ref-dynamic:#eb6834; --ref-fixed:#008300;
  --accent:#2a78d6; --warn-bg:#fff4d6; --warn-ink:#6b4a00; --ok-bg:#e3f4ea; --ok-ink:#0f5132; }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme=\"light\"])) .viz-root { color-scheme: dark; --surface-0:#111110;
  --surface-1:#1a1a19; --text-primary:#f0efec; --text-secondary:#c3c2b7; --text-muted:#9b9a92; --rule:#3a3a37; --grid:#2a2a28;
  --series-1:#3987e5; --series-2:#d95926; --series-3:#199e70; --series-4:#c98500; --series-5:#d55181; --series-6:#008300;
  --series-7:#9085e9; --series-8:#e66767; --div-neg:#3987e5; --div-mid:#383835; --div-pos:#e66767; --seq-lo:#1c2a3a; --seq-hi:#86b6ef;
  --ref-dynamic:#d95926; --ref-fixed:#2f9e2f; --accent:#3987e5; --warn-bg:#3a2e0c; --warn-ink:#f5d58a;
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

# Impute form: the run's default first (it is preselected) and marked. User decision 2026-11-04: limpa for raw, perseus for background.
impute_choices_ordered <- function(background, choices) {
  def <- if (background) "perseus" else "limpa"
  k <- which(vapply(choices, `[[`, "", "value") == def)
  choices[[k]]$label <- paste0(choices[[k]]$label, " - default for ", if (background) "background" else "raw", " exports")
  c(choices[k], choices[-k])
}

decision_spec <- function(step, variant, p, run_dir = NULL) {
  g <- function(k, d = "") { v <- p[[k]]; if (is.null(v) || length(v) == 0 || is.na(v[1])) d else as.character(v[1]) }
  spec <- switch(paste0(step, if (nzchar(variant)) paste0(":", variant) else ""),
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
      global = list(fld("norm_type", "Method", "select", "sltmm", if (exists("norm_method_names")) norm_method_names() else
        c("sltmm", "sl", "tmm", "median", "medianratio", "rlr", "loess", "rt", "rt_loess", "vsn", "quantile", "none"))),
      choices = list(ch("default_exclude", if (identical(g("ptm"), "1") || toupper(g("ptm")) == "TRUE") paste0("Pipeline default: ", g("ptm_grep"), " peptides only, added material excluded") else "Exclude: pipeline default (added material only)"),
                     if (identical(g("ptm"), "1") || toupper(g("ptm")) == "TRUE") ch("all_peptides", paste0("All peptides (not only ", g("ptm_grep"), "), added material excluded")),
                     ch("current_exclude", "Exclude: keep the current rule"),
                     ch("custom_exclude", "Exclude: my own regex", list(fld("norm_exclude_grep", "description regex", "text", g("norm_exclude_grep")))),
                     ch("preset_adh", "Factors from the yeast ADH spike only (designs where everything changes)"),
                     ch("preset_carboxylases", "Factors from the endogenous biotinylated carboxylases only (proximity labeling / streptavidin pulldown)"),
                     ch("custom_include", "Factors from my own reference proteins (e.g. the bait)", list(fld("norm_include_grep", "accession regex", "text", ""))))),
    "07_impute" = list(q = "How should the data be imputed?", nxt = "Step 8 - Rollup", choices = impute_choices_ordered(identical(g("sp_fill_mode"), "background"), list(
      ch("duke", "duke (near the group mean), bottom-X for the rest", list(fld("missing_cutoff", "missing_cutoff %", "number", g("missing_cutoff", "50")), fld("bottom_x", "bottom_x %", "number", g("bottom_x", "2")))),
      ch("bottomx", "bottomx only", list(fld("bottom_x", "bottom_x %", "number", g("bottom_x", "2")))),
      ch("perseus", "perseus (down-shifted normal per sample)", list(fld("perseus_shift", "shift (SD)", "number", g("perseus_shift", "1.8")), fld("perseus_width", "width (SD)", "number", g("perseus_width", "0.3")))),
      ch("minprob", "minprob", list(fld("minprob_q", "quantile", "number", g("minprob_q", "0.01")))),
      ch("qrilc", "qrilc"),
      ch("limpa", "limpa: no imputation (DPC-Quant rollup and dpcDE stats follow)")))),
    "08_rollup" = list(q = "Which rollup?", nxt = "Step 9 - QC", choices = if (!is.null(run_dir) && identical(applied_impute_type(file.path(run_dir, "x.db")), "limpa"))
      list(ch("limpa", "limpa DPC-Quant (Impute was limpa, so this is the only rollup)", list(fld("limpa_background", "Spectronaut background values (background exports)", "select", g("limpa_background", "missing"), c("missing", "use")))))
      else if (identical(g("data_output"), "Peptide")) list(ch("sum", "sum of the charge states (Peptide output: the only rollup without limpa)"))
      else list(ch("iq_maxlfq_scaled", "iq_maxlfq, scaled"), ch("iq_maxlfq", "iq_maxlfq, not scaled"), ch("sum", "sum"))),
    "09_qc" = list(q = "Accept the QC and continue to Stats?", nxt = "Step 10 - Stats", choices = list(
      ch("accept_flag", "Accept and flag samples", list(fld("flag_samples", "samples to flag (comma separated)", "text", ""))),
      ch("accept", "Accept"), ch("stop", "Stop here"))),
    "10_stats" = list(q = "Run stats with these settings?", nxt = "Stats results", choices = list(
      ch("run", "Run as shown"), ch("run_adjp", "Run, with Up/Down calls by adjusted p"),
      ch("change", "Change settings or comparisons", list(fld("change", "What to change?"))))),
    "10_stats:results" = list(q = "Accept the stats results?", nxt = "Finished (Pathway and Report are not built yet)", choices = list(
      ch("accept", "Accept"), ch("change", "Change cutoffs and re-run", list(fld("change", "New p type / p cutoff / FC cutoff"))),
      ch("plots", if (identical(g("data_output"), "Peptide")) "Add peptide plots (the peptides of these proteins)" else "Add protein plots", list(fld("proteins", "protein accessions (comma separated)"))))),
    NULL)
  if (!is.null(spec)) spec$choices <- Filter(Negate(is.null), spec$choices)
  spec
}

next_step_box <- function(run_dir, step, status_class, variant = "") {
  db <- file.path(run_dir, paste0(basename(run_dir), ".db"))
  if (!file.exists(db)) db <- list.files(run_dir, pattern = "\\.db$", full.names = TRUE)[1]
  p <- if (!is.na(db) && file.exists(db)) tryCatch(get_params(db), error = function(e) list()) else list()
  spec <- decision_spec(step, variant, p, run_dir)
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
  const API = location.pathname.replace(/[^\\/]+\\/review\\/[^\\/]*$/, "") + "api/";   // review_server: /api/; web app: <base>/runs/api/
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
  const who = el("input", { type: "text", size: "20", placeholder: "your name" }); who.value = window.DPDP_WEB_ANALYST || localStorage.getItem("dpdp_analyst") || ""; if (window.DPDP_WEB_ANALYST) who.readOnly = true;
  const bar = el("div", { class: "dbar" }); bar.appendChild(el("span", {}, "Note ")); bar.appendChild(note); bar.appendChild(el("span", {}, " Analyst ")); bar.appendChild(who);
  const btn = el("button", { type: "button", class: "dsend" }, "Send to Claude"); bar.appendChild(btn); box.appendChild(bar);
  const msg = el("p", { class: "dmsg" }); box.appendChild(msg);
  window.dpdpSetField = (name, value, force) => { box.querySelectorAll("[data-name=\\"" + name + "\\"]").forEach(i => { if (force || !i.dataset.edited) { i.value = value; if (force) i.dataset.edited = "1"; } }); };
  window.dpdpChoose = value => { const r = box.querySelector("input[name=dchoice][value=\\"" + value + "\\"]"); if (r) r.checked = true; };
  Object.entries(window.DPDP_PREFILL || {}).forEach(([k, v]) => window.dpdpSetField(k, v));
  box.addEventListener("input", e => { if (e.target.dataset) e.target.dataset.edited = "1"; });
  if (!served) { btn.disabled = true; msg.innerHTML = "Buttons need the review server. Ask Claude to start it, then open <code>http://127.0.0.1:8765/" + S.run + "/review/" + S.step + ".html</code>. Or answer in the Claude window."; return; }
  let mtime0 = null;
  const poll = () => fetch(API + "status?run=" + encodeURIComponent(S.run) + "&step=" + S.step).then(r => r.json()).then(st => {
    if (mtime0 === null) mtime0 = st.page_mtime; else if (st.page_mtime > mtime0) { location.reload(); return; }
    if (st.pending_step === S.step && !btn.disabled) { btn.disabled = true; msg.textContent = "A decision for this step is waiting for Claude."; }
  }).catch(() => {});
  poll(); setInterval(poll, 2000);
  btn.addEventListener("click", () => {
    const choice = box.querySelector("input[name=dchoice]:checked").value, fields = {};
    box.querySelectorAll("[data-name]").forEach(i => { if (i.dataset.pre === "" || i.dataset.pre === choice) fields[i.dataset.name] = i.type === "checkbox" ? String(i.checked) : i.value; });
    localStorage.setItem("dpdp_analyst", who.value);
    btn.disabled = true; msg.textContent = "Sending...";
    fetch(API + "decision", { method: "POST", headers: { "Content-Type": "application/json", "X-DPDP-Token": token },
      body: JSON.stringify({ run: S.run, step: S.step, choice: choice, fields: fields, note: note.value, analyst: who.value }) })
      .then(r => r.json().then(j => ({ ok: r.ok, j }))).then(({ ok, j }) => {
        if (!ok) { btn.disabled = false; msg.textContent = "Not sent: " + (j.error || "error"); return; }
        msg.textContent = "Sent to Claude at " + j.submitted.slice(11, 16) + ". This page reloads when Claude has applied it.";
        help.textContent = ""; })
      .catch(() => { btn.disabled = false; msg.textContent = "Not sent: is the review server running? (web app: reload the page)"; });
  });
})();
'

write_review_page <- function(run_dir, step, title, status, status_class, body, extra_head = "", decision_variant = "") {
  body <- paste0('<script>', COMMON_JS, '</script>', next_step_box(run_dir, step, status_class, decision_variant), body, '<script>', DECISION_JS, '</script>')
  need_plotly <- function() if (!grepl("Plotly.newPlot", body, fixed = TRUE) && !grepl("plotly", extra_head, fixed = TRUE)) body <<- paste0(body, '<script>', plotly_js(), '</script>')
  if (grepl("DPDP_CHARTS", body, fixed = TRUE)) { need_plotly(); body <- paste0(body, '<script>', SAMPLE_CHARTS_JS, '</script>') }
  if (grepl("DPDP_BARS", body, fixed = TRUE)) { need_plotly(); body <- paste0(body, '<script>', SAMPLE_BAR_JS, '</script>') }
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

# --- Shared interactive sample charts (replace the app's static bar / box / histogram PNGs) ---------------------------------
# sample_charts(m, labels, groups, id, ref) -> HTML + data for: total intensity per sample (bars by group), log2 intensity per
# sample (box plots from precomputed quantiles), and the log2 histogram of all values with reference lines.
# m: matrix of raw intensities (samples in columns). ref: named list of vertical lines for the histogram (log2 values).
sample_charts <- function(m, labels, groups, id, ref = list(), title_suffix = "", group_order = unique(groups)) {
  m[!is.na(m) & m <= 0] <- NA
  lx <- log2(m)
  q <- t(apply(lx, 2, function(v) as.numeric(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE))))
  edges <- seq(floor(min(lx, na.rm = TRUE) * 5) / 5, ceiling(max(lx, na.rm = TRUE) * 5) / 5, by = 0.2)
  h <- tabulate(findInterval(lx[!is.na(lx)], edges, rightmost.closed = TRUE), nbins = length(edges) - 1)
  d <- list(samples = labels, groups = groups, total = as.numeric(colSums(m, na.rm = TRUE)), measured = as.integer(colSums(!is.na(m))),
            q = unname(split(round(q, 3), row(q))), edges = edges, hist = h, ref = ref, order = group_order)
  sprintf('<div class="grid2"><div><h3>Total intensity per sample%s</h3><div id="%s_bar" style="height:320px" role="img" aria-label="Total intensity per sample"></div></div>
<div><h3>log2 intensity per sample%s</h3><div id="%s_box" style="height:320px" role="img" aria-label="log2 intensity distribution per sample"></div></div></div>
<h3>log2 intensity histogram, all values%s</h3><div id="%s_hist" style="height:340px" role="img" aria-label="Histogram of log2 intensities"></div>
<details><summary>Table view</summary><div id="%s_table"></div></details>
<script>window.DPDP_CHARTS = (window.DPDP_CHARTS || []).concat([{ id: "%s", d: %s }]);</script>',
    title_suffix, id, title_suffix, id, title_suffix, id, id, id, toJSON(d, auto_unbox = TRUE, digits = NA))
}

# Group colour = categorical slot by the group's position in the run's group order (entity, not rank; 8 slots, validated).
COMMON_JS <- '
window.DPDP_GCOL = function (g, order) { const i = Math.max(0, (order || []).indexOf(g));
  return getComputedStyle(document.querySelector(".viz-root")).getPropertyValue("--series-" + ((i % 8) + 1)).trim(); };
'

SAMPLE_BAR_JS <- '
(function () {
  const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
  const f3 = v => v === null ? "-" : (Math.abs(v) >= 1000 ? Math.round(v).toLocaleString("en-US") : (+v).toPrecision(3));
  (window.DPDP_BARS || []).forEach(({ id, d }) => {
    const names = Object.keys(d.series), sel = document.getElementById(id + "_sbs"), gs = d.order;
    function draw() {
      const k = sel ? sel.value : names[0], y = d.series[k];
      Plotly.react(id + "_sb", gs.filter(g => d.groups.includes(g)).map(g => { const ix = d.groups.map((x, i) => x === g ? i : -1).filter(i => i >= 0);
        return { type: "bar", name: g, x: ix.map(i => d.samples[i]), y: ix.map(i => y[i]), marker: { color: DPDP_GCOL(g, gs) },
          hovertemplate: "%{x}<br>" + k + ": %{y:,.4~g}<extra>" + g + "</extra>" }; }),
        Object.assign(base(), { margin: { l: 70, r: 10, t: 10, b: 60 }, legend: { orientation: "h", y: 1.12 }, bargap: 0.25,
          xaxis: { categoryorder: "array", categoryarray: d.samples }, yaxis: { title: d.ytitle[k], gridcolor: C("--grid"), rangemode: "tozero" } }),
        { displaylogo: false, responsive: true });
    }
    if (sel) sel.addEventListener("change", draw); draw();
    document.getElementById(id + "_sbt").innerHTML = "<table><thead><tr><th>Sample</th><th>Group</th>" + names.map(n => "<th class=num>" + n + "</th>").join("") +
      "</tr></thead><tbody>" + d.samples.map((s, i) => "<tr><td>" + s + "</td><td>" + d.groups[i] + "</td>" + names.map(n => "<td class=num>" + f3(d.series[n][i]) + "</td>").join("") + "</tr>").join("") + "</tbody></table>";
  });
})();
'

SAMPLE_CHARTS_JS <- '
(function () {
  const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const GC = ["--series-1", "--series-2", "--series-3", "--series-4"], fmt = n => Math.round(n).toLocaleString("en-US");
  const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
  const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
  (window.DPDP_CHARTS || []).forEach(({ id, d }) => {
    const gs = [...new Set(d.groups)], col = g => DPDP_GCOL(g, d.order || gs), idx = g => d.groups.map((x, i) => x === g ? i : -1).filter(i => i >= 0);
    const med = d.total.slice().sort((a, b) => a - b)[Math.floor(d.total.length / 2)];
    Plotly.newPlot(id + "_bar", gs.map(g => ({ type: "bar", name: g, x: idx(g).map(i => d.samples[i]), y: idx(g).map(i => d.total[i]), marker: { color: col(g) },
      customdata: idx(g).map(i => (d.total[i] / med).toFixed(3)), hovertemplate: "%{x}<br>total %{y:.3s}<br>%{customdata} x median<extra></extra>" })),
      Object.assign(base(), { margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: "h", y: 1.12 }, bargap: 0.25,
        xaxis: { categoryorder: "array", categoryarray: d.samples }, yaxis: Object.assign(ax("total intensity"), { rangemode: "tozero" }) }), { displaylogo: false, responsive: true });
    Plotly.newPlot(id + "_box", gs.map(g => ({ type: "box", name: g, x: idx(g).map(i => d.samples[i]),
      lowerfence: idx(g).map(i => d.q[i][0]), q1: idx(g).map(i => d.q[i][1]), median: idx(g).map(i => d.q[i][2]), q3: idx(g).map(i => d.q[i][3]), upperfence: idx(g).map(i => d.q[i][4]),
      marker: { color: col(g) }, fillcolor: col(g) + "55", line: { width: 1.5 }, hoverinfo: "y+x" })),
      Object.assign(base(), { margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: "h", y: 1.12 },
        xaxis: { categoryorder: "array", categoryarray: d.samples }, yaxis: ax("log2 intensity (box 25-75%, whiskers 5-95%)") }), { displaylogo: false, responsive: true });
    const centers = d.hist.map((_, i) => (d.edges[i] + d.edges[i + 1]) / 2), refs = Object.entries(d.ref || {});
    const top = Math.max(...d.hist);
    Plotly.newPlot(id + "_hist", [{ type: "bar", x: centers, y: d.hist, width: 0.2, marker: { color: C("--series-1") },
      customdata: centers.map(c => Math.pow(2, c).toPrecision(3)), hovertemplate: "log2 %{x:.1f} (intensity ~%{customdata}): %{y:,} values<extra></extra>" }],
      Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 50 }, showlegend: false, bargap: 0.02,
        shapes: refs.map(([k, r]) => ({ type: "line", x0: r.x, x1: r.x, y0: 0, y1: 1, yref: "paper", line: { color: C(r.color || "--text-muted"), width: 2, dash: r.dash || "dot" } })),
        annotations: refs.map(([k, r], j) => ({ x: r.x, y: 1 - j * 0.08, yref: "paper", xanchor: "left", xshift: 4, showarrow: false,
          text: k + " " + r.x.toFixed(1) + " (" + fmt(Math.pow(2, r.x)) + ")", font: { size: 11, color: C("--text-primary") } })),
        xaxis: ax("log2 intensity"), yaxis: Object.assign(ax("values"), { rangemode: "tozero" }) }), { displaylogo: false, responsive: true });
    document.getElementById(id + "_table").innerHTML = "<table><thead><tr><th>Sample</th><th>Group</th><th class=num>total</th><th class=num>x median</th><th class=num>values</th><th class=num>median log2</th></tr></thead><tbody>" +
      d.samples.map((s, i) => `<tr><td>${s}</td><td>${d.groups[i]}</td><td class=num>${d.total[i].toPrecision(4)}</td><td class=num>${(d.total[i] / med).toFixed(3)}</td><td class=num>${fmt(d.measured[i])}</td><td class=num>${d.q[i][2].toFixed(2)}</td></tr>`).join("") + "</tbody></table>";
  });
})();
'
