# Review page for Step 7b - imputation comparison. Informational: it doesn't record a decision; the Impute step does.

review_impute_compare <- function(db_path) {
  run_dir <- dirname(db_path)
  f <- file.path(run_dir, "steps", "07_impute", "compare.json")
  if (!file.exists(f)) stop("Run 07b_impute_compare.R first.")
  C <- fromJSON(f, simplifyVector = FALSE)
  M <- C$methods; ids <- names(M)
  g <- function(m, k) { v <- M[[m]][[k]]; if (is.null(v)) NA else v }
  gm <- function(m, k) { v <- M[[m]]$masked[[k]]; if (is.null(v)) NA else v }
  cols <- list(
    list(k = "rmse", lab = "Masking error (RMSE, log2)", src = "masked", low = TRUE),
    list(k = "rmse_mar", lab = "... gaps likely missed (MAR)", src = "masked", low = TRUE),
    list(k = "rmse_mnar", lab = "... gaps likely low (MNAR)", src = "masked", low = TRUE),
    list(k = "bias", lab = "Masking bias (log2)", src = "masked", low = NA),
    list(k = "null_pct_p01", lab = "No-difference splits: % p < 0.01 (expect 1)", src = "m", low = TRUE),
    list(k = "null_calls", lab = "No-difference splits: calls", src = "m", low = TRUE),
    list(k = "spqc_median_cv", lab = "SPQC median protein CV %", src = "m", low = TRUE),
    list(k = "adh1_cv", lab = "ADH1 CV % (9 samples)", src = "m", low = TRUE),
    list(k = "calls", lab = paste0(C$comparison, " calls (raw p, FC)"), src = "m", low = NA),
    list(k = "calls_onegroup", lab = "... measured in one group only", src = "m", low = NA),
    list(k = "onegroup_median_abs_log2fc", lab = "... their median |log2 FC|", src = "m", low = NA),
    list(k = "calls_adj", lab = "calls by adjusted p", src = "m", low = NA),
    list(k = "pi0", lab = "estimated share of unchanged proteins", src = "m", low = NA))
  val <- function(m, c) if (c$src == "masked") gm(m, c$k) else g(m, c$k)
  best <- lapply(cols, function(c) { v <- sapply(ids, function(m) suppressWarnings(as.numeric(val(m, c))))
    if (isTRUE(c$low) && any(!is.na(v))) ids[which(v == min(v, na.rm = TRUE))] else character(0) })
  head <- paste0("<tr><th>Method</th>", paste0("<th class=num>", esc(sapply(cols, `[[`, "lab")), "</th>", collapse = ""), "</tr>")
  rows <- vapply(ids, function(m) paste0("<tr><td>", esc(M[[m]]$label), "</td>", paste0(vapply(seq_along(cols), function(j) {
    v <- val(m, cols[[j]]); s <- if (is.null(v) || is.na(v)) "-" else format(v, big.mark = ",")
    sprintf('<td class="num%s">%s</td>', if (m %in% best[[j]]) " best" else "", s) }, ""), collapse = ""), "</tr>"), "")
  tbl <- sprintf("<div style='overflow-x:auto'><table class='cmp'><thead>%s</thead><tbody>%s</tbody></table></div>", head, paste(rows, collapse = ""))
  body <- paste0('<section><h2>What this compares</h2><p>Every method starts from the same filtered precursors (', fmt_n(C$start$precursors), ', ',
    fmt_n(C$start$gaps), ' gaps), scaled by SL factors from measured values. Each is rolled up with MaxLFQ (limpa uses its own DPC-Quant) and tested ',
    'with limma, using the same code. Lower is better for errors, false positives and CVs; the lowest value in those columns is marked. ',
    'There is no single score: read the columns together.</p><ul class="options">',
    sprintf('<li><b>Masking test:</b> %s measured values in %s complete precursors were hidden the way this dataset loses values (more often at low intensity), refilled, and compared with the truth. Not possible for Spectronaut background, no imputation or limpa.</li>', fmt_n(C$masking$hidden), fmt_n(C$masking$complete_precursors)),
    '<li><b>No-difference splits:</b> three comparisons where each side has one A, one B and one pool sample, so nothing should differ. Group-aware imputations (duke, mixed) are re-run with the two sides as their groups, as they would be in a real comparison. About 1% of proteins should have p &lt; 0.01 by chance; more means false positives.</li>',
    '<li><b>One group only:</b> calls where one group has no measured precursor; their fold change depends on how the gaps were filled.</li></ul></section>',
    '<section><h2>Scores</h2>', tbl, '</section>',
    '<section><h2>Charts</h2><div class="grid2"><div><h3>Masking error (lower is better)</h3><div id="c_rmse" style="height:300px"></div></div>',
    '<div><h3>False positives in no-difference splits (1% expected)</h3><div id="c_null" style="height:300px"></div></div>',
    '<div><h3>SPQC median protein CV (lower is better)</h3><div id="c_cv" style="height:300px"></div></div>',
    '<div><h3>', esc(C$comparison), ' calls and how many are one-group-only</h3><div id="c_calls" style="height:300px"></div></div></div>',
    '<h3>Raw p-value histograms, ', esc(C$comparison), '</h3><div id="c_ph" style="height:320px"></div></section>',
    sprintf('<p class="note">Run took %s s. Data: %s.</p>', C$seconds, esc(C$db_path)),
    '<script>', plotly_js(), '</script><script>const M = ', toJSON(M, auto_unbox = TRUE, digits = NA, null = "null", na = "null"), ';', CMP_JS, '</script>')
  extra <- '<style>.cmp td.best { font-weight:700; background:var(--ok-bg); color:var(--ok-ink); } .cmp th { font-size:12px; }</style>'
  write_review_page(run_dir, "07b_impute_compare", "Step 7b - Imputation comparison", "Information for the Impute decision", "done", body, extra_head = extra)
}

CMP_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
const ids = Object.keys(M), lab = ids.map(m => M[m].label);
const base = (yt) => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 11 },
  margin: { l: 240, r: 20, t: 10, b: 40 }, showlegend: false, yaxis: { automargin: true, autorange: "reversed" }, xaxis: { title: yt, gridcolor: C("--grid"), zeroline: false } });
const bar = (id, vals, xt, ref) => { const L = base(xt); if (ref !== undefined) L.shapes = [{ type: "line", x0: ref, x1: ref, y0: 0, y1: 1, yref: "paper", line: { dash: "dot", color: C("--text-muted") } }];
  Plotly.newPlot(id, [{ type: "bar", orientation: "h", y: lab, x: vals, marker: { color: C("--series-1") }, hovertemplate: "%{y}: %{x}<extra></extra>" }], L, { displayModeBar: false, responsive: true }); };
bar("c_rmse", ids.map(m => M[m].masked ? M[m].masked.rmse : null), "RMSE (log2)");
bar("c_null", ids.map(m => M[m].null_pct_p01), "% of proteins with p < 0.01", 1);
bar("c_cv", ids.map(m => M[m].spqc_median_cv), "median CV %");
Plotly.newPlot("c_calls", [
  { type: "bar", orientation: "h", name: "measured in both groups", y: lab, x: ids.map(m => M[m].calls - M[m].calls_onegroup), marker: { color: C("--series-1") } },
  { type: "bar", orientation: "h", name: "one group only", y: lab, x: ids.map(m => M[m].calls_onegroup), marker: { color: "#eb6834" } }],
  Object.assign(base("calls"), { barmode: "stack", showlegend: true, legend: { orientation: "h", y: 1.12 } }), { displayModeBar: false, responsive: true });
Plotly.newPlot("c_ph", ids.map((m, i) => ({ type: "scatter", mode: "lines", name: M[m].label, x: M[m].phist.map((_, k) => (k + 0.5) / 20), y: M[m].phist, line: { width: 2 } })),
  { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 11 }, margin: { l: 60, r: 20, t: 10, b: 45 },
    legend: { orientation: "h", y: -0.25 }, xaxis: { title: "raw p-value (bins of 0.05)", gridcolor: C("--grid") }, yaxis: { title: "proteins", gridcolor: C("--grid"), rangemode: "tozero" } },
  { displaylogo: false, responsive: true });
'
