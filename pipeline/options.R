# Step options: one list per step (pipeline/options/NN_<step>.R) that both the step script and its review page use.
#
# - The step script takes its defaults and checks from the list (opts_spec, opts_validate), so the page can't offer a value
#   the script would refuse.
# - The review page turns the list into the editable Options panel (options_form): every option, like the app's tab, with
#   Preview / Apply / Reset to defaults. Editing a field turns Apply off until a preview has rebuilt the page with those
#   settings (user decision 2026-11-09: apply only after preview).
# - The page sends all values with the decision (choice "preview" or "apply"); the step script reads them with --from_page
#   (page_args), so Claude doesn't retype them. Apply must match the last preview (apply_settings).

opt <- function(name, label, type = c("number", "checkbox", "select", "text"), default, meaning = "", choices = NULL,
                depends = NULL, min = NULL, max = NULL, step = NULL, app = TRUE, disabled = NULL, info = NULL) {
  type <- match.arg(type)
  list(name = name, label = label, type = type, default = default, meaning = meaning, choices = choices, depends = depends,
       min = min, max = max, step = step, app = app, disabled = disabled, info = info)
}
opt_section <- function(title, options, note = NULL) list(title = title, note = note, options = options)

opts_flat <- function(sections) do.call(c, lapply(sections, `[[`, "options"))
opts_spec <- function(sections) { o <- opts_flat(sections); setNames(lapply(o, `[[`, "default"), vapply(o, `[[`, "", "name")) }

# Stop on a value the step can't use (select outside its choices, number outside min/max).
opts_validate <- function(sections, settings) {
  for (o in opts_flat(sections)) {
    if (!is.null(o$depends) && !as_flag(settings[[o$depends]])) next   # option not in use
    v <- settings[[o$name]]
    if (o$type == "select" && !is.null(o$choices) && is.null(o$disabled) && !as.character(v) %in% o$choices)
      stop(sprintf("--%s must be one of: %s (got %s)", o$name, paste(o$choices, collapse = ", "), v), call. = FALSE)
    if (o$type == "number") {
      if (is.na(v)) stop(sprintf("--%s must be a number", o$name), call. = FALSE)
      if (!is.null(o$min) && v < o$min) stop(sprintf("--%s must be >= %s (got %s)", o$name, o$min, v), call. = FALSE)
      if (!is.null(o$max) && v > o$max) stop(sprintf("--%s must be <= %s (got %s)", o$name, o$max, v), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# --from_page: take the option values the analyst sent from the Options panel (review/pending_decision.json) as arguments.
# Arguments given on the command line win. A preview archives the page decision (it is not a decision for the log); an apply
# leaves it for log_decision(), which archives it and records the page choice.
page_args <- function(args, run_dir, step, sections, action) {
  if (!isTRUE(args$from_page) && !identical(args$from_page, "TRUE")) return(args)
  pend <- file.path(run_dir, "review", "pending_decision.json")
  if (!file.exists(pend)) stop("--from_page: no decision from the review page (", pend, ")", call. = FALSE)
  pd <- fromJSON(pend, simplifyVector = FALSE)
  if (!identical(pd$step, step)) stop("--from_page: the page decision is for step ", pd$step, ", not ", step, call. = FALSE)
  if (!identical(pd$choice, action)) stop(sprintf("--from_page: the page asked for '%s', not '%s'", pd$choice, action), call. = FALSE)
  for (k in intersect(vapply(opts_flat(sections), `[[`, "", "name"), names(pd$fields)))
    if (is.null(args[[k]])) args[[k]] <- as.character(pd$fields[[k]])
  if (is.null(args$note) && nzchar(pd$note %||% "")) args$note <- pd$note
  if (is.null(args$analyst) && nzchar(pd$analyst %||% "")) args$analyst <- pd$analyst
  if (action == "preview") {
    done <- file.path(run_dir, "review", "decisions_done"); dir.create(done, showWarnings = FALSE)
    file.rename(pend, file.path(done, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", step, "_preview.json")))
  }
  args
}

# Apply = the settings of the last preview. Option arguments (from --from_page or the command line) must equal them;
# a different value stops the step: preview first.
apply_settings <- function(run_dir, step, spec, args, settings) {
  f <- file.path(run_dir, "steps", step, "summary.json")
  prev <- if (file.exists(f)) fromJSON(f, simplifyVector = FALSE) else NULL
  if (is.null(prev) || !identical(prev$action, "preview"))
    stop("Apply needs a preview first: run --action=preview with these settings and show the page.", call. = FALSE)
  conv <- function(k, v) if (is.logical(spec[[k]])) as_flag(v) else if (is.numeric(spec[[k]])) as.numeric(v) else as.character(v)
  out <- settings
  for (k in names(spec)) if (!is.null(prev$settings[[k]])) out[[k]] <- conv(k, prev$settings[[k]])
  given <- intersect(names(spec), names(args))
  diff <- Filter(function(k) !isTRUE(all.equal(conv(k, args[[k]]), out[[k]])), given)
  if (length(diff)) stop("These settings differ from the last preview: ", paste(sprintf("%s (preview %s, now %s)", diff,
    vapply(diff, function(k) as.character(out[[k]]), ""), vapply(diff, function(k) as.character(args[[k]]), "")), collapse = "; "),
    ". Preview them first.", call. = FALSE)
  out
}

# The Options panel. values: what the form starts with; previewed: the settings the page shows (NULL if none); applied: read-only.
options_form <- function(step, sections, values, previewed = NULL, applied = FALSE, run_dir = NULL, title = "Options") {
  secs <- lapply(sections, function(s) {
    s$options <- lapply(s$options, function(o) { o$value <- values[[o$name]] %||% o$default; o })
    s })
  d <- list(step = step, run = if (!is.null(run_dir)) basename(run_dir) else "", applied = applied,
            previewed = if (is.null(previewed)) NULL else previewed[intersect(names(previewed), vapply(opts_flat(sections), `[[`, "", "name"))],
            sections = secs)
  sprintf('<section id="options"><h2>%s</h2>%s<div id="optform"></div></section><style>%s</style><script>window.DPDP_OPTIONS = %s;</script>',
    esc(title),
    if (applied) '<p class="note">These are the settings that were applied (read-only).</p>'
    else '<p class="note">Every option of this step, as in the app. Change any of them, press <b>Preview with these settings</b> to see the effect on this page, then <b>Apply</b>. Changed options are marked; <b>pipeline</b> marks options the app does not have.</p>',
    OPTIONS_CSS, toJSON(d, auto_unbox = TRUE, null = "null", digits = NA))
}

OPTIONS_CSS <- "
#options table.opts { width:100%; border-collapse:collapse; margin:4px 0 14px; table-layout:fixed; }
#options table.opts td { padding:5px 8px; border-bottom:1px solid var(--grid); vertical-align:top; white-space:normal; overflow-wrap:anywhere; }
#options td.olab { width:30%; font-weight:500; } #options td.oval { width:17%; white-space:nowrap; } #options td.omean { width:53%; color:var(--text-secondary); font-size:13px; }
#options td.oval input[type=number], #options td.oval input[type=text], #options td.oval select { font:inherit; font-size:14px; padding:2px 6px; width:9em; }
#options tr.changed td.olab::before { content:'\\25CF  '; color:var(--warn-ink); }
#options tr.changed td.oval input, #options tr.changed td.oval select { border:2px solid var(--warn-ink); border-radius:4px; }
#options tr.off td { opacity:0.45; } #options .odef { display:block; font-size:12px; color:var(--warn-ink); margin-top:2px; }
#options .pipe { font-size:11px; padding:0 5px; border:1px solid var(--rule); border-radius:8px; color:var(--text-secondary); margin-left:4px; }
#options h3 { margin:14px 0 2px; }
#optbar { position:sticky; bottom:0; background:var(--surface-1); border-top:2px solid var(--accent); padding:10px 4px; display:flex; gap:8px; align-items:center; flex-wrap:wrap; z-index:5; }
#optbar button { font:inherit; font-weight:600; padding:6px 14px; border-radius:6px; border:1px solid var(--accent); background:var(--accent); color:#fff; cursor:pointer; }
#optbar button.quiet { background:transparent; color:var(--accent); } #optbar button:disabled { opacity:0.45; cursor:default; }
#optbar input { font:inherit; font-size:14px; padding:3px 6px; } #optmsg { flex-basis:100%; font-weight:600; }
"

OPTIONS_JS <- '
(function () {
  const O = window.DPDP_OPTIONS; if (!O) return;
  const root = document.getElementById("optform");
  const API = location.pathname.replace(/[^\\/]+\\/review\\/[^\\/]*$/, "") + "api/";
  const tokEl = document.querySelector("meta[name=dpdp-token]"), token = tokEl ? tokEl.content : "";
  const served = location.protocol.startsWith("http") && token && token !== "__DPDP_" + "TOKEN__";
  const el = (t, a = {}, h = "") => { const e = document.createElement(t); Object.entries(a).forEach(([k, v]) => e.setAttribute(k, v)); e.innerHTML = h; return e; };
  const isTrue = v => v === true || String(v).toLowerCase() === "true" || String(v) === "1";
  const eq = (o, a, b) => { if (a === undefined || a === null || b === undefined || b === null) return a == b;
    if (o.type === "checkbox") return isTrue(a) === isTrue(b);
    if (o.type === "number") return Number(a) === Number(b);
    return String(a) === String(b); };
  const show = (o, v) => o.type === "checkbox" ? (isTrue(v) ? "on" : "off") : String(v);
  const X = {};
  O.sections.forEach(sec => {
    root.appendChild(el("h3", {}, sec.title));
    if (sec.note) root.appendChild(el("p", { class: "note" }, sec.note));
    const tb = el("table", { class: "opts" });
    sec.options.forEach(o => {
      const tr = el("tr"), td = el("td", { class: "oval" });
      tr.appendChild(el("td", { class: "olab" }, o.label + (o.app === false ? \'<span class="pipe">pipeline</span>\' : "")));
      let inp;
      if (o.type === "checkbox") { inp = el("input", { type: "checkbox" }); inp.checked = isTrue(o.value); }
      else if (o.type === "select") { inp = el("select"); (o.choices || []).forEach(c => inp.appendChild(el("option", { value: c }, c))); inp.value = String(o.value); }
      else { inp = el("input", { type: o.type === "number" ? "number" : "text", step: o.step || "any" });
        if (o.min !== null && o.min !== undefined) inp.min = o.min; if (o.max !== null && o.max !== undefined) inp.max = o.max; inp.value = o.value; }
      inp.setAttribute("aria-label", o.label.replace(/<[^>]+>/g, ""));
      td.appendChild(inp); const def = el("span", { class: "odef" }); td.appendChild(def); tr.appendChild(td);
      tr.appendChild(el("td", { class: "omean" }, (o.meaning || "") + (o.info ? " <b>" + o.info + "</b>" : "") + (o.disabled ? \' <span class="missing">\' + o.disabled + "</span>" : "")));
      tb.appendChild(tr); X[o.name] = { o: o, inp: inp, tr: tr, def: def };
      inp.addEventListener("input", update); inp.addEventListener("change", update);
    });
    root.appendChild(tb);
  });
  const val = x => x.inp.type === "checkbox" ? x.inp.checked : x.inp.value;
  const bar = el("div", { id: "optbar" });
  const reset = el("button", { type: "button", class: "quiet" }, "Reset to defaults");
  const prev = el("button", { type: "button" }, "Preview with these settings");
  const apply = el("button", { type: "button" }, "Apply these settings");
  const note = el("input", { type: "text", size: "40", placeholder: "optional: why (goes in the decision log)" });
  const who = el("input", { type: "text", size: "16", placeholder: "your name" });
  who.value = window.DPDP_WEB_ANALYST || localStorage.getItem("dpdp_analyst") || ""; if (window.DPDP_WEB_ANALYST) who.readOnly = true;
  const msg = el("span", { id: "optmsg" });
  [reset, prev, apply].forEach(b => bar.appendChild(b)); bar.appendChild(el("span", {}, "Note")); bar.appendChild(note);
  bar.appendChild(el("span", {}, "Analyst")); bar.appendChild(who); bar.appendChild(msg);
  if (!O.applied) root.appendChild(bar);
  let pending = false, sent = false;
  function update() {
    const dirty = [];
    Object.values(X).forEach(x => {
      const o = x.o, dep = o.depends && X[o.depends], on = !dep || dep.inp.checked;
      x.inp.disabled = O.applied || !!o.disabled || !on;
      x.tr.classList.toggle("off", !on || !!o.disabled);
      const ch = !eq(o, val(x), o.default);
      x.tr.classList.toggle("changed", ch); x.def.textContent = ch ? "default: " + show(o, o.default) : "";
      if (O.previewed && !eq(o, val(x), O.previewed[o.name])) dirty.push(o.label.replace(/<[^>]+>/g, ""));
    });
    if (O.applied) return;
    const busy = pending || sent || !served;
    prev.disabled = busy; reset.disabled = busy;
    apply.disabled = busy || !O.previewed || dirty.length > 0;
    if (!served) msg.textContent = "Buttons need the review server (or the DPDP web app). You can also tell Claude in the chat.";
    else if (sent || pending) msg.textContent = "Sent to Claude. This page reloads when Claude has run it.";
    else if (!O.previewed) msg.textContent = "Preview these settings first; Apply turns on when the page shows them.";
    else if (dirty.length) msg.textContent = "Changed since the preview: " + dirty.join(", ") + ". Preview to see the effect before applying.";
    else msg.textContent = "The page shows these settings. Apply them, or change options and preview again.";
  }
  reset.onclick = () => { Object.values(X).forEach(x => { const d = x.o.default;
    if (x.inp.type === "checkbox") x.inp.checked = isTrue(d); else x.inp.value = String(d); }); update(); };
  const send = choice => {
    const fields = {}; Object.values(X).forEach(x => { fields[x.o.name] = String(val(x)); });
    localStorage.setItem("dpdp_analyst", who.value); sent = true; update(); msg.textContent = "Sending...";
    fetch(API + "decision", { method: "POST", headers: { "Content-Type": "application/json", "X-DPDP-Token": token },
      body: JSON.stringify({ run: O.run, step: O.step, choice: choice, fields: fields, note: note.value, analyst: who.value }) })
      .then(r => r.json().then(j => ({ ok: r.ok, j }))).then(({ ok, j }) => {
        if (!ok) { sent = false; update(); msg.textContent = "Not sent: " + (j.error || "error"); return; }
        msg.textContent = (choice === "apply" ? "Apply" : "Preview") + " sent to Claude at " + j.submitted.slice(11, 16) + ". This page reloads when Claude has run it."; })
      .catch(() => { sent = false; update(); msg.textContent = "Not sent: is the server running? Reload the page."; });
  };
  prev.onclick = () => send("preview");
  apply.onclick = () => send("apply");
  update();
  if (served && !O.applied) {
    let mtime0 = null;
    const poll = () => fetch(API + "status?run=" + encodeURIComponent(O.run) + "&step=" + O.step).then(r => r.json()).then(st => {
      if (mtime0 === null) mtime0 = st.page_mtime; else if (st.page_mtime > mtime0) { location.reload(); return; }
      const p = st.pending_step === O.step; if (p !== pending) { pending = p; update(); } }).catch(() => {});
    poll(); setInterval(poll, 2000);
  }
})();
'
