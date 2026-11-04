# Review page for Step 10 - Stats (Setup, Comparisons, Graphs, Protein Plots, Excel).
# Decisions: stats settings + comparisons (before the calculation); then accept the results or change cutoffs.
# Shows: settings and comparisons; significance with raw p vs adjusted p side by side; an interactive volcano where the
# p type and cutoffs can be changed (counts and the significant list update; same Up/Down rule as the app's Stats column);
# raw p-value histogram; PCA; the static graphs from 10b; the protein plots from 10c; output files.

review_stats <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  sd_dir <- file.path(run_dir, "steps", "10_stats")
  rd <- function(f) if (file.exists(file.path(sd_dir, f))) fromJSON(file.path(sd_dir, f), simplifyVector = FALSE) else NULL
  summ <- rd("summary.json"); graphs <- rd("graphs.json"); outp <- rd("output.json")
  applied <- !is.null(summ) && !identical(summ$action, "preview")
  design <- read_table("design", db_path)
  sc <- if (table_exists("stats_comp", db_path)) read_table("stats_comp", db_path) else NULL
  s <- summ$settings %||% list()
  g <- function(k) s[[k]] %||% p[[k]]
  pc <- as.numeric(g("pvalue_cutoff")); fcc <- as.numeric(g("foldchange_cutoff")); use_adj <- as_flag(g("checkbox_filter_adjpval"))
  locs <- function(x) as.integer(trimws(strsplit(as.character(x), ",")[[1]]))
  sp <- qc_spike_settings(run_dir)

  # --- per-comparison data for the volcano ---
  comps <- list()
  if (applied && !is.null(sc)) for (i in seq_len(nrow(sc))) {
    nm <- sc$Name[i]; tb <- sc$Final_Table_Name[i]
    if (!table_exists(tb, db_path)) next
    d <- read_table(tb, db_path)
    fc <- d[[paste0(nm, "_FC")]]; pv <- d[[paste0(nm, "_pval")]]; ap <- d[[paste0(nm, "_adjpval")]] %||% rep(NA_real_, nrow(d))
    nN <- length(locs(sc$N_loc[i])); nD <- length(locs(sc$D_loc[i]))
    di <- strsplit(as.character(d$Detected_Imputed), ".", fixed = TRUE)
    onoff <- vapply(di, function(v) { if (length(v) < nN + nD) return(""); a <- all(v[1:nN] == "-"); b <- all(v[(nN + 1):(nN + nD)] == "-")
      if (b && !a) paste0("only in ", sc$FactorsN[i]) else if (a && !b) paste0("only in ", sc$FactorsD[i]) else "" }, "")
    comps[[nm]] <- list(name = nm, N = sc$FactorsN[i], D = sc$FactorsD[i], acc = d$Accession, gene = d$Genes, fc = signif(fc, 8),
                        p = signif(pv, 8), adj = signif(ap, 8), mf = d[[paste0(nm, "_mf")]], onoff = onoff, stats = d$Stats,
                        prec = d$Precursors, flag = flag_contaminants(d$Accession, d$Description, sp$spike, sp$constant))
  }
  cres <- if (!is.null(summ$comparisons)) setNames(summ$comparisons, vapply(summ$comparisons, function(x) x$comparison, "")) else list()

  status <- if (applied) "Stats calculated: review the results" else "Waiting for approval of the settings"
  status_class <- "waiting"
  dec <- read_decisions(run_dir, "10_stats")
  acc_dec <- Filter(function(r) identical(r$settings$accepted, TRUE), dec); last_acc <- if (length(acc_dec)) acc_dec[[length(acc_dec)]] else NULL
  if (!is.null(last_acc) && !is.null(summ$generated) && last_acc$time < summ$generated) last_acc <- NULL
  if (!is.null(last_acc)) { status <- "Results accepted"; status_class <- "done" }

  set_df <- data.frame(Setting = c("Normalization tested", "Pooled QC carried along", "p-value cutoff", "Fold change cutoff",
                                   "Up/Down calls use", "Adjusted p column", "Missing factor", "Peptide refilter"),
    Value = c(esc(summ$stat_norm %||% p$stat_norm), esc(summ$comp_spqc %||% p$comp_spqc), esc(pc), esc(fcc),
              if (use_adj) sprintf("<b>adjusted p</b> (%s)", esc(g("padjust_options"))) else "<b>raw p</b> (checkbox_filter_adjpval off)",
              if (as_flag(g("checkbox_adjpval"))) esc(g("padjust_options")) else "off", esc(g("missing_factor")), if (as_flag(g("peptide_refilter"))) "on" else "off"),
    Meaning = c("Which normalized data the tests use.", "Pool shown in tables and plots; not part of the test.", "",
                "Absolute fold change needed for Up/Down.", "Raw p gives more calls but no false-discovery control; adjusted p controls FDR.",
                "Added to the tables either way.", "Proteins with a smaller share of measured (not imputed) values in the comparison are removed.",
                "Re-filters precursors for the comparison's samples and rolls up again."), stringsAsFactors = FALSE)
  comp_df <- if (!is.null(sc) && nrow(sc)) data.frame(Comparison = esc(sc$Name), Numerator = esc(sc$FactorsN), Denominator = esc(sc$FactorsD),
    Samples = sprintf("%s vs %s (+%s pool)", sc$N, sc$D, sc$SPQC), stringsAsFactors = FALSE) else NULL

  decision_html <- paste0('<section><h2>', if (applied) "Settings used" else "Decision needed: approve the settings and comparisons", '</h2>',
    html_table(set_df), if (!is.null(comp_df)) paste0('<h3>Comparisons</h3>', html_table(comp_df)) else "",
    if (!is.null(last_acc)) sprintf('<p class="note">Results accepted %s%s. %s</p>', esc(last_acc$time), if (!is.null(last_acc$analyst)) paste0(" by ", esc(last_acc$analyst)) else "", esc(last_acc$note %||% "")) else "",
    if (applied && is.null(last_acc)) '<p class="note">Review the results below, then accept them in the Claude session or change cutoffs or comparisons and re-run.</p>' else "",
    '</section>')

  res_html <- ""
  if (length(cres)) res_html <- paste0('<section><h2>Significance</h2>', paste(vapply(cres, function(r) {
    pd <- r$pvalue_distribution
    warn <- if (!use_adj && (r$significant_if_adjusted_p %||% 0) == 0 && (r$significant %||% 0) > 0)
      sprintf('<p class="missing">The %s Up/Down calls use <b>raw p</b>. <b>None pass FDR</b> (smallest adjusted p = %s). With %s proteins tested, about %s would have raw p &lt; %s by chance alone.</p>',
              fmt_n(r$significant), esc(pd$min_adjpval), fmt_n(r$proteins), fmt_n(round(r$proteins * pc)), esc(pc)) else ""
    paste0('<h3>', esc(r$comparison), ' (', esc(r$N), ' vs ', esc(r$D), ')</h3><div class="tiles">',
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">proteins tested (%s removed by missing factor)</div></div>', fmt_n(r$proteins), fmt_n(r$removed_by_missing_factor)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">Up / Down by %s p: %s / %s</div></div>', fmt_n(r$significant), if (use_adj) "adjusted" else "raw", fmt_n(r$up), fmt_n(r$down)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">would pass with adjusted p &le; %s</div></div>', fmt_n(r$significant_if_adjusted_p), esc(pc)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">smallest adjusted p</div></div>', esc(pd$min_adjpval)),
      sprintf('<div class="tile"><div class="v">%s / %s / %s</div><div class="l">raw p &lt; 0.05 / 0.01 / 0.001</div></div>', fmt_n(pd$raw_p_lt_0.05), fmt_n(pd$raw_p_lt_0.01), fmt_n(pd$raw_p_lt_0.001)),
      '</div>', warn)
  }, ""), collapse = ""), '</section>')

  volc_html <- if (length(comps)) paste0('<section><h2>Volcano explorer</h2>',
    '<p class="note">Change the p type and cutoffs to see how the Up/Down calls change (same rule as the app&apos;s Stats column). ',
    'Ringed points were measured in only one group (the other group is all imputed); their large fold changes come from imputation. ',
    'Diamonds mark <b>flagged</b> proteins (digestion trypsin, human keratins, QC spike-ins): kept in the analysis, but likely not biology.</p>',
    '<div class="control">', if (length(comps) > 1) paste0('<label for="vcomp">Comparison</label><select id="vcomp">', paste0('<option>', esc(names(comps)), '</option>', collapse = ""), '</select>') else "",
    '<label for="vpt">p type</label><select id="vpt"><option value="p"', if (!use_adj) " selected" else "", '>raw p</option><option value="adj"', if (use_adj) " selected" else "", '>adjusted p</option></select>',
    '<label for="vpc">p cutoff</label><select id="vpc">', paste0('<option', ifelse(c(0.001, 0.01, 0.05, 0.1, 0.2) == pc, " selected", ""), '>', c(0.001, 0.01, 0.05, 0.1, 0.2), '</option>', collapse = ""), '</select>',
    '<label for="vfc">FC cutoff</label><input type="range" id="vfc" min="1" max="4" step="0.1" value="', fcc, '"/><span class="val" id="vfcv"></span>',
    '<label><input type="checkbox" id="vonoff"/> leave out proteins measured in one group only</label></div>',
    '<div class="tiles" id="vtiles"></div><div id="volcano" style="height:520px" role="img" aria-label="Volcano plot"></div>',
    '<details open><summary>Significant proteins at these settings (table view)</summary><div class="control"><input type="search" id="vq" placeholder="filter by gene or accession" style="font:inherit;font-size:13px;padding:3px 8px"/></div><div id="vtable" style="max-height:420px;overflow:auto"></div></details>',
    '<h2>Raw p-value distribution</h2><p class="note">Under no real differences the bars are flat. A peak at low p means real signal; how far the first bar rises above the flat level is a rough count of true differences.</p>',
    '<div id="phist" style="height:260px" role="img" aria-label="Histogram of raw p-values"></div><p class="note" id="phistnote"></p></section>') else ""

  pca_html <- ""
  if (!is.null(graphs) && length(graphs$info)) {
    sets <- Filter(function(x) !is.null(x$pca), graphs$info)
    if (length(sets)) pca_html <- paste0('<section><h2>PCA</h2><div class="control"><label for="pset">Samples</label><select id="pset">',
      paste0('<option>', esc(names(sets)), '</option>', collapse = ""), '</select></div><div id="pca" style="height:440px" role="img" aria-label="PCA of samples"></div>',
      '<p class="note" id="pcanote"></p></section>')
  }

  stat_plots <- unlist(graphs$plots); prot_plots <- unlist(outp$protein_plots$plots %||% outp$plots)
  xl <- list.files(file.path(run_dir, summ$stat_norm %||% p$stat_norm %||% ""), pattern = "\\.xlsx$", full.names = TRUE)
  graphs_html <- paste0(
    if (length(stat_plots)) paste0('<section><h2>Stats graphs</h2><div class="grid2">', paste(vapply(stat_plots, function(f) embed_png(f, basename(f), sub("\\.png$", "", basename(f))), ""), collapse = ""), '</div></section>') else "",
    if (length(prot_plots)) paste0('<section><h2>Protein plots</h2><div class="grid2">', paste(vapply(prot_plots, function(f) embed_png(f, basename(f), sub("\\.png$", "", basename(f))), ""), collapse = ""), '</div></section>') else "",
    if (length(xl)) paste0('<section><h2>Output files</h2><ul>', paste0('<li><code>', esc(xl), '</code></li>', collapse = ""), '</ul>',
      if (!use_adj) '<p class="note">The Excel Up/Down column uses raw p (checkbox_filter_adjpval off). Table1 Sample Info is created empty by the app.</p>' else "", '</section>') else "")

  js_data <- list(comps = comps, pc = pc, fcc = fcc, pca = if (!is.null(graphs)) lapply(Filter(function(x) !is.null(x$pca), graphs$info), function(x) list(pca = x$pca, nearest = x$nearest_sample)) else NULL)
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA, na = "null"), ';', STATS_JS, '</script>')
  extra_css <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100; --up:#e34948; --down:#2a78d6; } td.flag { color:var(--warn-ink); font-weight:600; } section h3 { font-size:15px; margin:14px 0 4px; }
.control select { font:inherit; font-size:13px; } .control input[type=range] { width:160px; } .control .val { min-width:60px; }</style>'
  body <- paste0(decision_html, res_html, volc_html, pca_html, graphs_html, js)
  write_review_page(run_dir, "10_stats", "Step 10 - Stats", status, status_class, body, extra_head = extra_css,
                    decision_variant = if (applied) "results" else "")
}

STATS_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
const cfg = { displaylogo: false, responsive: true };
const l2 = f => f === null ? null : (f < 0 ? -Math.log2(-f) : Math.log2(f));
const names = Object.keys(D.comps);
function vdraw() {
  if (!names.length) return;
  const cs = document.getElementById("vcomp"), X = D.comps[cs ? cs.value : names[0]];
  const pt = document.getElementById("vpt").value, pc = +document.getElementById("vpc").value, fc = +document.getElementById("vfc").value, hide = document.getElementById("vonoff").checked;
  document.getElementById("vfcv").textContent = fc.toFixed(1) + "x";
  const P = X[pt], n = X.acc.length, call = new Array(n);
  let up = 0, down = 0, upOn = 0, downOn = 0;
  for (let i = 0; i < n; i++) { const q = P[i]; let c = "";
    if (q !== null && q <= pc && X.fc[i] !== null) { if (X.fc[i] >= fc) c = "Up"; else if (X.fc[i] <= -fc) c = "Down"; }
    if (hide && X.onoff[i]) c = "";
    call[i] = c; if (c === "Up") { up++; if (X.onoff[i]) upOn++; } if (c === "Down") { down++; if (X.onoff[i]) downOn++; } }
  document.getElementById("vtiles").innerHTML = `<div class="tile"><div class="v">${fmt(up + down)}</div><div class="l">significant (${pt === "p" ? "raw" : "adjusted"} p ≤ ${pc}, |FC| ≥ ${fc.toFixed(1)})</div></div>` +
    `<div class="tile"><div class="v">${fmt(up)} / ${fmt(down)}</div><div class="l">Up / Down (${X.N} vs ${X.D})</div></div>` +
    `<div class="tile"><div class="v">${fmt(upOn + downOn)}</div><div class="l">of these measured in one group only</div></div>` +
    `<div class="tile"><div class="v">${fmt(call.filter((c, i) => (c === "Up" || c === "Down") && X.flag[i]).length)}</div><div class="l">of these flagged (contaminant / spike-in)</div></div>`;
  const groups = { "": [], Up: [], Down: [] }, ring = [];
  for (let i = 0; i < n; i++) { if (P[i] === null) continue; if (hide && X.onoff[i]) continue; groups[call[i]].push(i); if (X.onoff[i]) ring.push(i); }
  const tr = (idx, name, col, size) => ({ type: "scattergl", mode: "markers", name: name, x: idx.map(i => l2(X.fc[i])), y: idx.map(i => -Math.log10(P[i])),
    text: idx.map(i => `${X.gene[i] || ""} ${X.acc[i]}<br>FC ${(+X.fc[i]).toPrecision(4)}  p ${(+X.p[i]).toPrecision(3)}  adj p ${X.adj[i] === null ? "-" : (+X.adj[i]).toPrecision(3)}<br>${X.prec[i]} precursors${X.onoff[i] ? "<br>" + X.onoff[i] : ""}`),
    hovertemplate: "%{text}<extra></extra>", marker: { size: size, color: C(col), opacity: name === "not significant" ? 0.35 : 0.85 } });
  const traces = [tr(groups[""], "not significant", "--text-muted", 4), tr(groups.Up, "Up", "--up", 6), tr(groups.Down, "Down", "--down", 6)];
  traces.push({ type: "scattergl", mode: "markers", name: "one group only", x: ring.map(i => l2(X.fc[i])), y: ring.map(i => -Math.log10(P[i])), hoverinfo: "skip",
    marker: { size: 10, color: "rgba(0,0,0,0)", line: { width: 1.5, color: C("--text-primary") } } });
  const fl = []; for (let i = 0; i < n; i++) if (X.flag[i] && P[i] !== null && !(hide && X.onoff[i])) fl.push(i);
  traces.push({ type: "scattergl", mode: "markers", name: "flagged", x: fl.map(i => l2(X.fc[i])), y: fl.map(i => -Math.log10(P[i])),
    text: fl.map(i => `${X.gene[i] || ""} ${X.acc[i]}<br>${X.flag[i]}`), hovertemplate: "%{text}<extra></extra>",
    marker: { symbol: "diamond-open", size: 12, color: C("--series-4"), line: { width: 2, color: C("--series-4") } } });
  const lf = Math.log2(fc), ly = -Math.log10(pc);
  const xs = X.fc.filter(v => v !== null).map(v => Math.abs(l2(v))), xm = Math.max(...xs) * 1.05;
  Plotly.react("volcano", traces, Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 50 }, showlegend: true, legend: { orientation: "h", y: 1.07 },
    shapes: [{ type: "line", x0: lf, x1: lf, y0: 0, y1: 1, yref: "paper", line: { dash: "dot", color: C("--text-muted") } }, { type: "line", x0: -lf, x1: -lf, y0: 0, y1: 1, yref: "paper", line: { dash: "dot", color: C("--text-muted") } },
      { type: "line", x0: 0, x1: 1, xref: "paper", y0: ly, y1: ly, line: { dash: "dot", color: C("--text-muted") } }],
    xaxis: Object.assign(ax("log2 fold change (" + X.N + " / " + X.D + ")"), { range: [-xm, xm] }), yaxis: ax(pt === "p" ? "-log10 raw p" : "-log10 adjusted p") }), cfg);
  const q = (document.getElementById("vq").value || "").toLowerCase();
  const sig = []; for (let i = 0; i < n; i++) { const c = call[i]; if ((c === "Up" || c === "Down") && (!q || (X.gene[i] || "").toLowerCase().includes(q) || X.acc[i].toLowerCase().includes(q))) sig.push(i); }
  sig.sort((a, b) => P[a] - P[b]);
  document.getElementById("vtable").innerHTML = `<table><thead><tr><th>Accession</th><th>Gene</th><th>Call</th><th class=num>FC</th><th class=num>raw p</th><th class=num>adj p</th><th class=num>precursors</th><th>Measured</th><th>Flag</th></tr></thead><tbody>` +
    sig.slice(0, 500).map(i => `<tr><td>${X.acc[i]}</td><td>${X.gene[i] || ""}</td><td>${call[i]}</td><td class=num>${(+X.fc[i]).toPrecision(4)}</td><td class=num>${(+X.p[i]).toPrecision(3)}</td><td class=num>${X.adj[i] === null ? "-" : (+X.adj[i]).toPrecision(3)}</td><td class=num>${X.prec[i]}</td><td>${X.onoff[i]}</td><td class="${X.flag[i] ? "flag" : ""}">${X.flag[i]}</td></tr>`).join("") + "</tbody></table>" +
    (sig.length > 500 ? `<p class="note">First 500 of ${sig.length}.</p>` : "");
  // p histogram
  const bins = new Array(20).fill(0); X.p.forEach(v => { if (v !== null) bins[Math.min(19, Math.floor(v * 20))]++; });
  const flat = bins.slice(10).reduce((s, v) => s + v, 0) / 10, excess = Math.max(0, Math.round(bins[0] - flat));
  Plotly.react("phist", [{ type: "bar", x: bins.map((_, i) => (i + 0.5) / 20), y: bins, width: 0.045, marker: { color: C("--series-1") }, hovertemplate: "p %{x:.3f}: %{y:,} proteins<extra></extra>" }],
    Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 45 }, showlegend: false,
      shapes: [{ type: "line", x0: 0, x1: 1, y0: flat, y1: flat, line: { dash: "dot", color: C("--text-muted") } }],
      annotations: [{ x: 0.98, y: flat, text: "flat level, estimated from p > 0.5", showarrow: false, yanchor: "bottom", xanchor: "right", font: { size: 12, color: C("--text-primary") } }],
      xaxis: ax("raw p-value"), yaxis: Object.assign(ax("proteins"), { rangemode: "tozero" }) }), cfg);
  document.getElementById("phistnote").textContent = `First bar (p < 0.05): ${fmt(bins[0])} proteins against a flat level of about ${fmt(flat)}, so roughly ${fmt(excess)} more than chance alone would give.`;
}
["vcomp", "vpt", "vpc", "vonoff"].forEach(id => { const e = document.getElementById(id); if (e) e.addEventListener("change", vdraw); });
["vfc", "vq"].forEach(id => { const e = document.getElementById(id); if (e) e.addEventListener("input", vdraw); });
vdraw();
if (D.pca && Object.keys(D.pca).length) {
  const ps = document.getElementById("pset");
  function pdraw() {
    const s = D.pca[ps.value], sc = s.pca.scores, gr = [...new Set(sc.map(r => r.group))], cols = ["--series-1", "--series-2", "--series-3", "--up"];
    Plotly.react("pca", gr.map((g, k) => { const r = sc.filter(z => z.group === g); return { type: "scatter", mode: "markers+text", name: g, x: r.map(z => z.PC1), y: r.map(z => z.PC2),
      text: r.map(z => z.Sample), textposition: "top center", marker: { size: 12, color: C(cols[k % 4]), line: { width: 2, color: C("--surface-1") } }, hovertemplate: "%{text}<br>PC1 %{x:.1f}<br>PC2 %{y:.1f}<extra></extra>" }; }),
      Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 50 }, showlegend: true, legend: { orientation: "h", y: 1.08 },
        xaxis: ax("PC1 (" + s.pca.variance_pct.PC1 + "%)"), yaxis: ax("PC2 (" + s.pca.variance_pct.PC2 + "%)") }), cfg);
    document.getElementById("pcanote").textContent = s.nearest ? "Nearest sample (sample distance): " + Object.entries(s.nearest).map(([a, b]) => a + " → " + b).join(", ") : "";
  }
  ps.addEventListener("change", pdraw); pdraw();
}
'
