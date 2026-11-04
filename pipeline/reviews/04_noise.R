# Review page for Step 4 - Noise (human decision step).
# Shows: the decision to make and what each option does (with counts), the app's inflection and start histogram plots,
# an interactive cutoff explorer (values / precursors below any cutoff, per sample), and the applied result once decided.

review_noise <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  start <- read_table("precursor_start", db_path)
  cols <- (ncol(start) - n_s + 1):ncol(start)
  design <- read_table("design", db_path)
  labels <- if (!is.null(design$Label) && length(design$Label) == n_s) design$Label else names(start)[cols]
  m <- as.matrix(start[, cols]); m[m <= 0] <- NA
  lg <- log2(m)
  sd_dir <- file.path(run_dir, "steps", "04_noise")
  applied <- if (file.exists(file.path(sd_dir, "summary.json"))) fromJSON(file.path(sd_dir, "summary.json")) else NULL
  if (!is.null(applied) && !identical(applied$action, "apply")) applied <- NULL
  decisions <- read_decisions(run_dir, "04_noise")
  last_dec <- if (length(decisions)) decisions[[length(decisions)]] else NULL

  # --- explorer data: binned log2 intensities, cumulative counts below each bin edge ---
  step <- 0.1
  lo <- floor(min(lg, na.rm = TRUE) / step) * step; hi <- ceiling(max(lg, na.rm = TRUE) / step) * step
  edges <- round(seq(lo, hi, by = step), 4)
  cum_below <- function(v) { v <- v[!is.na(v)]; c(0, cumsum(tabulate(findInterval(v, edges, rightmost.closed = TRUE), nbins = length(edges) - 1))) }
  per_sample <- lapply(seq_len(n_s), function(j) cum_below(lg[, j]))
  hist_total <- diff(Reduce(`+`, per_sample))
  k_min <- as.integer(p$filter_min_measured_all %||% 2); if (is.na(k_min) || k_min < 1) k_min <- 2
  measured <- rowSums(!is.na(lg))
  kth <- apply(lg, 1, function(r) { r <- sort(r[!is.na(r)], decreasing = TRUE); if (length(r) >= k_min) r[k_min] else NA })
  rmin <- suppressWarnings(apply(lg, 1, min, na.rm = TRUE)); rmax <- suppressWarnings(apply(lg, 1, max, na.rm = TRUE))
  ok <- measured > 0
  explorer <- list(
    edges = edges, hist = hist_total, samples = labels, per_sample = per_sample,
    sample_values = as.integer(colSums(!is.na(lg))),
    rows_total = sum(ok), rows_short = sum(ok & measured < k_min), k_min = k_min,
    cum_rowmin = cum_below(rmin[ok]), cum_rowmax = cum_below(rmax[ok]), cum_kth = cum_below(kth[ok & measured >= k_min]),
    inflection = as.numeric(p$noise_inflection), fixed = as.numeric(p$noise_baseline_value))

  infl <- as.numeric(p$noise_inflection); fixed <- as.numeric(p$noise_baseline_value)
  total <- as.numeric(p$noise_total)
  dyn_n <- as.numeric(p$dynamic_noise_count); fix_n <- as.numeric(p$custom_noise_count)
  have_infl <- !is.na(infl) && infl > 0

  # --- decision section ---
  status <- if (is.null(applied)) "Waiting for the analyst's decision" else sprintf("Decided: %s", applied$settings$noise_type)
  status_class <- if (is.null(applied)) "waiting" else "done"
  options_html <- sprintf('<ol class="options">
<li><b>none</b> - keep every value. Values removed: 0.</li>
<li><b>fixed</b> - remove values below a set intensity (<code>noise_baseline_value</code>, now %s). Values removed at %s: <b>%s</b> (%s of %s).</li>
<li><b>dynamic</b> - remove values below the inflection point of the intensity curve (%s, log2 %.2f). Values removed: <b>%s</b> (%s).</li>
</ol>
<p><b>noise_keep</b> (with fixed or dynamic): noise values are removed so the filter treats them as missing; precursors that are
mostly noise then fail the filter. For every precursor that <i>passes</i> the filter, its measured values are put back, so they are not imputed.
Without noise_keep, removed values stay missing and are imputed later.</p>',
    esc(fixed), esc(fixed), fmt_n(fix_n), fmt_pct(fix_n / total * 100), fmt_n(total),
    esc(infl), if (have_infl) log2(infl) else NA, fmt_n(dyn_n), fmt_pct(dyn_n / total * 100))
  decision_html <- if (is.null(applied)) paste0(
    '<section><h2>Decision needed</h2>', options_html,
    '<p class="note">Look at the graphs and the explorer below, then give your choice in the Claude session. ',
    'Claude records it in <code>decisions.jsonl</code> and only then applies it.</p></section>')
  else {
    s <- applied$settings
    rem <- applied$values_set_to_na_per_sample
    rem_df <- data.frame(Sample = labels, Removed = fmt_n(unlist(rem)), Percent = fmt_pct(unlist(rem) / explorer$sample_values * 100))
    paste0('<section><h2>Decision</h2><div class="tiles">',
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">noise_type</div></div>', esc(s$noise_type)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">baseline used</div></div>', esc(applied$baseline_used)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">values removed</div></div>', fmt_n(applied$values_set_to_na)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">noise_keep</div></div>', if (isTRUE(s$noise_keep)) "on" else "off"),
      '</div>',
      if (!is.null(last_dec)) sprintf('<p class="note">Recorded %s%s: %s</p>', esc(last_dec$time),
                                      if (!is.null(last_dec$analyst)) paste0(" by ", esc(last_dec$analyst)) else "", esc(last_dec$note %||% "")) else "",
      if (applied$values_set_to_na > 0) paste0("<h2>Values removed per sample</h2>", html_table(rem_df, 2:3)) else "",
      '<details><summary>Options that were available</summary>', options_html, '</details></section>')
  }

  tiles <- paste0('<section><h2>Key numbers</h2><div class="tiles">',
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">precursors</div></div>', fmt_n(nrow(start))),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">measured data points</div></div>', fmt_n(total)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">inflection point (log2 %.2f)</div></div>', esc(infl), if (have_infl) log2(infl) else NA),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">below inflection (%s)</div></div>', fmt_n(dyn_n), fmt_pct(dyn_n / total * 100)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">below fixed baseline %s</div></div>', fmt_n(fix_n), esc(fixed)),
    '</div></section>')

  qc <- p$qc_path
  # Interactive version of the app's inflection plot: every data point sorted from highest to lowest log2 intensity
  # (every 50th point drawn), with the inflection (dynamic) and fixed baseline marked.
  sv <- sort(lg[!is.na(lg)], decreasing = TRUE); keep <- unique(c(seq(1, length(sv), by = 50), length(sv)))
  curve <- list(x = keep, y = round(sv[keep], 3), n = length(sv),
                infl_log2 = if (have_infl) log2(infl) else NA, infl_rank = if (have_infl) sum(sv >= log2(infl)) else NA,
                fixed_log2 = if (fixed > 0) log2(fixed) else NA)
  graphs <- paste0('<section><h2>Graphs</h2><h3>Sorted intensities (inflection plot)</h3>',
    '<p class="note">Every measured value, sorted from highest to lowest. The curve drops steeply at the low end; the inflection is where that drop starts. ',
    'Values right of the orange line are what <b>dynamic</b> removes; zoom into the tail to see it.</p>',
    '<div id="inflplot" style="height:380px" role="img" aria-label="Sorted log2 intensities with the inflection point"></div>',
    sprintf('<script>window.DPDP_CURVE = %s;</script>', toJSON(curve, auto_unbox = TRUE, digits = NA, na = "null")),
    '<details><summary>The app&apos;s original plots</summary><div class="grid2">',
    embed_png(file.path(qc, "Inflection_Point.png"), "Inflection point plot", "App plot: inflection point."),
    embed_png(file.path(qc, "Precursor_Start_Histogram.png"), "Start histogram", "App plot: start histogram."),
    if (!is.null(applied)) embed_png(file.path(qc, "Precursor_NoiseFiltered_Histogram.png"), "Noise filtered histogram", "App plot: histogram after the noise step.") else "",
    '</div></details></section>')

  explorer_html <- sprintf('<section><h2>Cutoff explorer</h2>
<p class="note">Move the slider to see what any cutoff would remove. It is a guide for choosing; the exact counts for fixed and
dynamic are in the options above. Precursors with fewer than %d values left would fail the filter&apos;s minimum-measured rule
(<code>filter_min_measured_all</code> = %d).</p>
<div class="control"><label for="cut">Cutoff (log2)</label><input id="cut" type="range" min="0" max="%d" step="1" value="0"/>
<span class="val" id="cutval"></span>
<button type="button" data-jump="inflection">Go to inflection</button><button type="button" data-jump="fixed">Go to fixed baseline</button></div>
<div class="tiles" id="cuttiles"></div>
<div id="histplot" style="height:380px" role="img" aria-label="Histogram of log2 intensities of all data points, with reference lines"></div>
<details open><summary>Per sample at this cutoff (table view)</summary><div id="cuttable"></div></details>
</section>', k_min, k_min, length(edges) - 1)

  js <- sprintf('<script>%s</script><script>
const D = %s;
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const pct = (a, b) => b ? (100 * a / b).toFixed(2) + "%%" : "-";
const centers = D.hist.map((_, i) => (D.edges[i] + D.edges[i + 1]) / 2);
const idxFor = v => Math.min(D.edges.length - 1, Math.max(0, Math.round((Math.log2(v) - D.edges[0]) / (D.edges[1] - D.edges[0]))));
function refLine(x, color, dash, text, yref) {
  return { shape: { type: "line", x0: x, x1: x, y0: 0, y1: 1, yref: "paper", line: { color: color, width: 2, dash: dash } },
           ann: { x: x, y: yref, yref: "paper", text: text, showarrow: false, xanchor: "left", xshift: 4, font: { size: 12, color: C("--text-primary") } } };
}
const fixedLog = D.fixed > 0 ? Math.log2(D.fixed) : null, inflLog = D.inflection > 0 ? Math.log2(D.inflection) : null;
const refs = [];
if (inflLog !== null) refs.push(refLine(inflLog, C("--ref-dynamic"), "solid", "dynamic: inflection " + D.inflection, 0.98));
if (fixedLog !== null && fixedLog >= D.edges[0]) refs.push(refLine(fixedLog, C("--ref-fixed"), "dot", "fixed baseline " + D.fixed, 0.88));
const layout = { margin: { l: 70, r: 20, t: 10, b: 50 }, paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"),
  font: { color: C("--text-secondary"), size: 12 }, bargap: 0.05, showlegend: false,
  xaxis: { title: "log2 intensity", gridcolor: C("--grid"), zeroline: false },
  yaxis: { title: "data points", gridcolor: C("--grid"), zeroline: false, rangemode: "tozero" } };
const trace = { type: "bar", x: centers, y: D.hist, width: 0.1, marker: { color: C("--series-1") },
  customdata: D.hist.map((_, i) => [Math.pow(2, D.edges[i]).toPrecision(4), Math.pow(2, D.edges[i + 1]).toPrecision(4)]),
  hovertemplate: "log2 %%{x:.2f}<br>intensity %%{customdata[0]} to %%{customdata[1]}<br>%%{y:,} data points<extra></extra>" };
Plotly.newPlot("histplot", [trace], layout, { displaylogo: false, responsive: true });
const slider = document.getElementById("cut");
function update() {
  const k = +slider.value, e = D.edges[k], raw = Math.pow(2, e);
  document.getElementById("cutval").textContent = "log2 " + e.toFixed(1) + "  =  intensity " + raw.toPrecision(4);
  const removed = D.per_sample.map(c => c[k]), tot = removed.reduce((a, b) => a + b, 0), allv = D.sample_values.reduce((a, b) => a + b, 0);
  const touched = D.cum_rowmin[k], gone = D.cum_rowmax[k], fail = D.cum_kth[k] + D.rows_short;
  document.getElementById("cuttiles").innerHTML =
    `<div class="tile"><div class="v">${fmt(tot)}</div><div class="l">values below (${pct(tot, allv)})</div></div>` +
    `<div class="tile"><div class="v">${fmt(touched)}</div><div class="l">precursors with a value below</div></div>` +
    `<div class="tile"><div class="v">${fmt(gone)}</div><div class="l">precursors with every value below</div></div>` +
    `<div class="tile"><div class="v">${fmt(fail)}</div><div class="l">precursors with &lt; ${D.k_min} values above (incl. ${fmt(D.rows_short)} already short)</div></div>`;
  document.getElementById("cuttable").innerHTML = "<table><thead><tr><th>Sample</th><th class=num>values below</th><th class=num>%% of sample</th><th class=num>values measured</th></tr></thead><tbody>" +
    D.samples.map((s, i) => `<tr><td>${s}</td><td class=num>${fmt(removed[i])}</td><td class=num>${pct(removed[i], D.sample_values[i])}</td><td class=num>${fmt(D.sample_values[i])}</td></tr>`).join("") + "</tbody></table>";
  const cut = refLine(e, C("--text-primary"), "dash", "explore " + raw.toPrecision(4), 0.78);
  const all = refs.concat([cut]);
  Plotly.relayout("histplot", { shapes: all.map(r => r.shape), annotations: all.map(r => r.ann) });
}
slider.addEventListener("input", update);
document.querySelectorAll("[data-jump]").forEach(b => b.addEventListener("click", () => {
  const v = b.dataset.jump === "inflection" ? D.inflection : D.fixed; if (v > 0) { slider.value = idxFor(v); update(); } }));
slider.value = D.inflection > 0 ? idxFor(D.inflection) : 0; update();
(function () { const K = window.DPDP_CURVE; if (!K) return;
  const sh = [], an = [];
  if (K.infl_rank !== null) { sh.push({ type: "line", x0: K.infl_rank, x1: K.infl_rank, y0: 0, y1: 1, yref: "paper", line: { color: C("--ref-dynamic"), width: 2 } });
    an.push({ x: K.infl_rank, y: 0.95, yref: "paper", xanchor: "right", xshift: -4, showarrow: false, text: "dynamic: " + fmt(K.n - K.infl_rank) + " values below " + D.inflection, font: { size: 12, color: C("--text-primary") } }); }
  if (K.fixed_log2 !== null) { sh.push({ type: "line", x0: 0, x1: 1, xref: "paper", y0: K.fixed_log2, y1: K.fixed_log2, line: { color: C("--ref-fixed"), width: 2, dash: "dot" } });
    an.push({ x: 0.02, xref: "paper", y: K.fixed_log2, yanchor: "bottom", showarrow: false, text: "fixed baseline " + D.fixed, font: { size: 12, color: C("--text-primary") } }); }
  Plotly.newPlot("inflplot", [{ type: "scattergl", mode: "lines", x: K.x, y: K.y, line: { width: 2, color: C("--series-1") },
    hovertemplate: "rank %%{x:,}<br>log2 %%{y:.2f}<extra></extra>" }],
    { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, margin: { l: 60, r: 20, t: 10, b: 50 },
      showlegend: false, shapes: sh, annotations: an, xaxis: { title: "rank (highest to lowest)", gridcolor: C("--grid"), zeroline: false },
      yaxis: { title: "log2 intensity", gridcolor: C("--grid"), zeroline: false } }, { displaylogo: false, responsive: true });
})();
</script>', plotly_js(), toJSON(explorer, auto_unbox = TRUE, digits = NA))

  body <- paste0(decision_html, tiles, graphs, explorer_html, js)
  write_review_page(run_dir, "04_noise", "Step 4 - Noise", status, status_class, body)
}
