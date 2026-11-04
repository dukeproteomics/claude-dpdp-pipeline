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
  groups_tbl <- read_table("sample_groups", db_path)
  gcols_n <- setNames(lapply(seq_len(nrow(groups_tbl)), function(i) as.integer(groups_tbl$start[i]):as.integer(groups_tbl$end[i])), groups_tbl$Group)
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

  # --- replicate reproducibility: does a value at this intensity repeat in the other replicates? ---
  pool <- if (!is.null(p$comp_spqc) && p$comp_spqc %in% groups_tbl$Group) p$comp_spqc else if ("SPQC" %in% groups_tbl$Group) "SPQC" else NA
  rep_sets <- if (!is.na(pool)) list(gcols_n[[pool]]) else Filter(function(cc) length(cc) >= 3, gcols_n)
  rep_src <- if (!is.na(pool)) sprintf("the %d %s injections (one pooled sample run %d times: pure technical replicates)", length(gcols_n[[pool]]), pool, length(gcols_n[[pool]]))
             else "replicates within each group (biological + technical variation)"
  rb <- seq(floor(min(lg, na.rm = TRUE) * 2) / 2, ceiling(max(lg, na.rm = TRUE) * 2) / 2, by = 0.5)
  acc <- list(n = 0, all = 0, none = 0); cvl <- list(); mnl <- list()
  for (cc in rep_sets) {
    sub_lg <- lg[, cc, drop = FALSE]; k <- length(cc)
    for (j in seq_len(k)) { v <- sub_lg[, j]; ok <- !is.na(v); oth <- rowSums(!is.na(sub_lg[ok, -j, drop = FALSE]))
      b <- findInterval(v[ok], rb, rightmost.closed = TRUE)
      acc$n <- acc$n + tabulate(b, length(rb) - 1); acc$all <- acc$all + tabulate(b[oth == k - 1], length(rb) - 1); acc$none <- acc$none + tabulate(b[oth == 0], length(rb) - 1) }
    full <- rowSums(!is.na(sub_lg)) == k; raw <- 2^sub_lg[full, , drop = FALSE]
    cvl[[length(cvl) + 1]] <- apply(raw, 1, stats::sd) / rowMeans(raw) * 100; mnl[[length(mnl) + 1]] <- rowMeans(sub_lg[full, , drop = FALSE])
  }
  cvv <- unlist(cvl); mnv <- unlist(mnl); cb <- findInterval(mnv, rb, rightmost.closed = TRUE)
  med_cv <- vapply(seq_len(length(rb) - 1), function(i) { x <- cvv[cb == i]; if (length(x) >= 20) round(stats::median(x), 1) else NA_real_ }, 0)
  ok_bin <- acc$n >= 50
  redet <- ifelse(ok_bin, round(100 * acc$all / pmax(acc$n, 1), 1), NA); never <- ifelse(ok_bin, round(100 * acc$none / pmax(acc$n, 1), 1), NA)
  cross <- function(y, bad) { i <- which(!is.na(y) & bad(y)); if (!length(i)) NA else rb[max(i) + 1] }   # lowest log2 above which every bin is fine
  repro <- list(centers = (head(rb, -1) + tail(rb, -1)) / 2, n = acc$n, redet = redet, never = never, cv = med_cv, src = rep_src,
                cv50 = cross(med_cv, function(y) y >= 50), redet70 = cross(redet, function(y) y < 70),
                infl = if (have_infl) log2(infl) else NA)
  repro_html <- paste0('<section id="repro"><h2>Replicate reproducibility</h2>',
    sprintf('<p class="note">Noise is signal that does not repeat. Using %s: for each intensity, the share of values that are also measured in <b>all</b> the other replicates, and the replicate CV. ', esc(rep_src)),
    'A noise floor shows as the point where re-detection falls and CV climbs. Move the slider to pick a cutoff; <b>Use as fixed baseline</b> puts it in the decision form above.</p>',
    '<div class="control"><label for="rcut">Cutoff (log2)</label><input id="rcut" type="range" min="0" max="', length(rb) * 4, '" step="1"/><span class="val" id="rcutv"></span>',
    '<button type="button" id="ruse">Use as fixed baseline</button><button type="button" data-rjump="cv50">CV 50% point</button><button type="button" data-rjump="redet70">re-detected 70% point</button><button type="button" data-rjump="infl">current inflection</button></div>',
    '<div class="tiles" id="rtiles"></div><div class="grid2"><div><h3>Re-detected in all other replicates</h3><div id="rredet" style="height:300px" role="img" aria-label="Share of values re-detected by intensity"></div></div>',
    '<div><h3>Median replicate CV</h3><div id="rcv" style="height:300px" role="img" aria-label="Replicate CV by intensity"></div></div></div>',
    '<details><summary>Table view</summary>', html_table(data.frame(`log2 intensity` = sprintf("%.1f-%.1f", head(rb, -1), tail(rb, -1)), Values = fmt_n(acc$n),
      `Re-detected in all others` = ifelse(is.na(redet), "-", paste0(redet, "%")), `In no other` = ifelse(is.na(never), "-", paste0(never, "%")),
      `Median CV` = ifelse(is.na(med_cv), "-", paste0(med_cv, "%")), check.names = FALSE)[acc$n > 0, ], 2:5), '</details>',
    sprintf('<script>window.DPDP_REPRO = %s;</script>', toJSON(repro, auto_unbox = TRUE, digits = NA, na = "null")), '</section>')

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

  body <- paste0(decision_html, tiles, graphs, repro_html, explorer_html, js, '<script>', REPRO_JS, '</script>')
  write_review_page(run_dir, "04_noise", "Step 4 - Noise", status, status_class, body)
}

REPRO_JS <- '
(function () {
  const R = window.DPDP_REPRO; if (!R) return;
  const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const fmt = n => Math.round(n).toLocaleString("en-US");
  const lo = R.centers[0] - 0.25, step = 0.125, sl = document.getElementById("rcut");
  const xOf = k => lo + k * step, kOf = x => Math.max(0, Math.min(+sl.max, Math.round((x - lo) / step)));
  const at = (arr, x) => { let i = R.centers.findIndex(c => x < c + 0.25); if (i < 0) i = R.centers.length - 1; return arr[i]; };
  const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, margin: { l: 55, r: 15, t: 10, b: 45 }, showlegend: false });
  const refs = [["CV 50%", R.cv50, "--series-3", "dash"], ["re-detected 70%", R.redet70, "--ref-fixed", "dash"], ["inflection", R.infl, "--ref-dynamic", "solid"]].filter(r => r[1] !== null);
  function draw() {
    const x = xOf(+sl.value), raw = Math.pow(2, x);
    document.getElementById("rcutv").textContent = "log2 " + x.toFixed(2) + " = intensity " + raw.toPrecision(3);
    const shapes = refs.map(r => ({ type: "line", x0: r[1], x1: r[1], y0: 0, y1: 1, yref: "paper", line: { color: C(r[2]), width: 2, dash: r[3] } }))
      .concat([{ type: "line", x0: x, x1: x, y0: 0, y1: 1, yref: "paper", line: { color: C("--text-primary"), width: 2, dash: "dot" } }]);
    const anns = refs.map((r, j) => ({ x: r[1], y: 0.98 - j * 0.1, yref: "paper", xanchor: "left", xshift: 3, showarrow: false, text: r[0], font: { size: 11, color: C("--text-primary") } }));
    const tr = (y, name, fmtY) => [{ type: "scatter", mode: "lines+markers", x: R.centers, y: y, line: { width: 2, color: C("--series-1") }, marker: { size: 6, color: C("--series-1") },
      customdata: R.n, hovertemplate: "log2 %{x:.2f}: " + fmtY + "<br>%{customdata:,} values<extra></extra>" }];
    Plotly.react("rredet", tr(R.redet, "re-detected", "%{y:.1f}% re-detected"), Object.assign(base(), { shapes, annotations: anns,
      xaxis: { title: "log2 intensity", gridcolor: C("--grid"), zeroline: false, range: [lo, 16] }, yaxis: { title: "% of values", range: [0, 100], gridcolor: C("--grid") } }), { displayModeBar: false, responsive: true });
    Plotly.react("rcv", tr(R.cv, "CV", "median CV %{y:.1f}%"), Object.assign(base(), { shapes, annotations: anns,
      xaxis: { title: "log2 intensity", gridcolor: C("--grid"), zeroline: false, range: [lo, 16] }, yaxis: { title: "median CV %", rangemode: "tozero", gridcolor: C("--grid") } }), { displayModeBar: false, responsive: true });
    let below = 0; if (typeof D !== "undefined") { const ki = D.edges.findIndex(e => e >= x - 1e-9); const kk = ki < 0 ? D.edges.length - 1 : ki; below = D.per_sample.reduce((s, c) => s + c[kk], 0); }
    const tot = (typeof D !== "undefined") ? D.sample_values.reduce((a, b) => a + b, 0) : 0;
    const rd = at(R.redet, x), cv = at(R.cv, x);
    document.getElementById("rtiles").innerHTML =
      `<div class="tile"><div class="v">${rd === null ? "-" : rd + "%"}</div><div class="l">re-detected in all other replicates at this intensity</div></div>` +
      `<div class="tile"><div class="v">${cv === null ? "-" : cv + "%"}</div><div class="l">median replicate CV at this intensity</div></div>` +
      `<div class="tile"><div class="v">${fmt(below)}</div><div class="l">values below this cutoff (${tot ? (100 * below / tot).toFixed(2) : "-"}% of all)</div></div>`;
  }
  sl.addEventListener("input", draw);
  document.querySelectorAll("[data-rjump]").forEach(b => { const v = R[b.dataset.rjump]; if (v === null) { b.disabled = true; return; } b.addEventListener("click", () => { sl.value = kOf(v); draw(); }); });
  document.getElementById("ruse").addEventListener("click", () => {
    const raw = Math.round(Math.pow(2, xOf(+sl.value)) * 100) / 100;
    if (window.dpdpChoose) window.dpdpChoose("fixed"); if (window.dpdpSetField) window.dpdpSetField("noise_baseline_value", String(raw), true);
    const d = document.getElementById("decision"); if (d) d.scrollIntoView({ behavior: "smooth" }); });
  sl.value = kOf(R.cv50 !== null ? R.cv50 : (R.infl !== null ? R.infl : R.centers[0])); draw();
})();
'
