# Review page for Step 5 - Filter.
# Preview state: the settings to approve, the noise-filtered histogram, and "what-if" explorers that show how many precursors
#   each filter would remove on its own (from precursor_noise, the filter's input) across its range of settings.
# Applied state: adds the filter cascade, counts before/after, missing values per sample, and the app's bar/box plots
#   (start vs filtered).
# The what-if counts follow the app's rules (minimum_filter, filter_min_group, misaligned test_alignment, filter_percentCV,
# precursor_spqc_quality). Precursor quality has no what-if (it is a per-protein iterative test); its count comes from the run.

review_filter <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  groups <- read_table("sample_groups", db_path)
  design <- read_table("design", db_path)
  df <- read_table("precursor_noise", db_path)
  info <- ncol(df) - n_s
  m <- as.matrix(df[, (info + 1):ncol(df)]); m[!is.na(m) & m == 0] <- NA
  labels <- if (!is.null(design$Label) && length(design$Label) == n_s) design$Label else colnames(m)
  n_rows <- nrow(m)
  grp_of_all <- rep(groups$Group, as.integer(groups$Count))
  sd_dir <- file.path(run_dir, "steps", "05_filter")
  summ <- if (file.exists(file.path(sd_dir, "summary.json"))) fromJSON(file.path(sd_dir, "summary.json"), simplifyVector = FALSE) else NULL
  applied <- !is.null(summ) && identical(summ$action, "apply")
  s <- if (!is.null(summ)) summ$settings else list()
  getp <- function(k) { v <- s[[k]] %||% p[[k]]; v }
  flag <- function(k) as_flag(getp(k))
  num <- function(k) as.numeric(getp(k))
  cum_ge <- function(x, grid) vapply(grid, function(g) sum(x >= g), 0)   # count of rows with x >= g

  # --- what-if 1: minimum measured across all samples (remove rows with measured < k) ---
  measured <- rowSums(!is.na(m))
  wi_min <- list(x = 1:n_s, series = list(all = vapply(1:n_s, function(k) sum(measured < k), 0)), current = num("filter_min_measured_all"))

  # --- group helpers ---
  gcols <- lapply(seq_len(nrow(groups)), function(i) as.integer(groups$start[i]):as.integer(groups$end[i]))
  names(gcols) <- groups$Group
  frac <- sapply(gcols, function(cc) rowSums(!is.na(m[, cc, drop = FALSE])) / length(cc))
  if (is.null(dim(frac))) frac <- matrix(frac, ncol = length(gcols))

  # --- what-if 2: % measured in at least one group (remove rows whose best group is below X%) ---
  maxfrac <- apply(frac, 1, max)
  xg <- 0:100
  wi_pct <- list(x = xg, series = list(all = vapply(xg, function(x) sum(maxfrac < x / 100 - 1e-9), 0)), current = num("filter_x_percent_value"))

  # --- what-if 3: misaligned (a group with missing% > cutoff, < 100, and mean >= intensity cutoff) ---
  lmean <- as.numeric(p$intensity_mean); lsd <- as.numeric(p$intensity_sd)
  sds <- seq(0, 2, by = 0.25)
  gmean <- sapply(gcols, function(cc) rowMeans(m[, cc, drop = FALSE], na.rm = TRUE))
  missp <- (1 - frac) * 100
  mis_series <- lapply(sds, function(k) {
    cut <- 2^(lmean + k * lsd)
    q <- missp; q[!(q < 100 & !is.na(gmean) & gmean >= cut)] <- -1
    best <- apply(q, 1, max)
    vapply(xg, function(x) sum(best > x), 0)
  })
  names(mis_series) <- sprintf("mean + %.2f SD (%s)", sds, formatC(round(2^(lmean + sds * lsd)), format = "d", big.mark = ","))
  cur_sd <- num("intensity_cutoff_sd")
  wi_mis <- list(x = xg, series = mis_series, current = num("misaligned_cutoff"),
                 current_series = names(mis_series)[which.min(abs(sds - cur_sd))])

  # --- what-if 4: CV of one group (remove rows with CV >= value; CV NA counts as 0, as in the app) ---
  cvx <- seq(0, 200, by = 1)
  cv_series <- lapply(gcols, function(cc) {
    sub <- m[, cc, drop = FALSE]; av <- rowMeans(sub, na.rm = TRUE); sdv <- apply(sub, 1, stats::sd, na.rm = TRUE)
    cv <- signif(100 * sdv / av, 3); cv[is.na(cv)] <- 0
    cum_ge(cv, cvx)
  })
  wi_cv <- list(x = cvx, series = cv_series, current = num("filter_cv_value"), current_series = getp("filter_cv_group"))

  # --- what-if 5: SPQC ratio (SPQC mean / sample mean off by more than +-accuracy %, SPQC mean > intensity) ---
  wi_spqc <- NULL
  spqc_cols <- grep("SPQC", colnames(m))
  if (length(spqc_cols) && length(spqc_cols) < n_s) {
    sample_avg <- rowMeans(m[, -spqc_cols, drop = FALSE], na.rm = TRUE); spqc_avg <- rowMeans(m[, spqc_cols, drop = FALSE], na.rm = TRUE)
    acc <- spqc_avg / sample_avg * 100; dev <- abs(acc - 100)
    ints <- sort(unique(c(0, 100, 250, 500, 1000, 2000, 5000, num("precursor_spqc_intensity"))))
    ax <- seq(5, 200, by = 5)
    sp_series <- lapply(ints, function(it) { ok <- !is.na(dev) & !is.na(spqc_avg) & spqc_avg > it; vapply(ax, function(a) sum(ok & dev > a), 0) })
    names(sp_series) <- paste0("SPQC mean > ", formatC(ints, format = "d", big.mark = ","))
    wi_spqc <- list(x = ax, series = sp_series, current = num("precursor_spqc_accuracy"),
                    current_series = paste0("SPQC mean > ", formatC(num("precursor_spqc_intensity"), format = "d", big.mark = ",")))
  }

  # --- settings table ---
  onoff <- function(k) if (flag(k)) "<b>on</b>" else "off"
  set_rows <- list(
    c("Minimum measured, all samples", "always on", sprintf("at least <b>%s</b> of %d samples", esc(getp("filter_min_measured_all")), n_s),
      "Removes precursors measured in too few samples overall."),
    c("Measured in X% of a group", onoff("filter_x_percent"), sprintf("X = <b>%s%%</b>", esc(getp("filter_x_percent_value"))),
      "Keeps a precursor only if at least one group has it in X% of its samples."),
    c("Misaligned", onoff("checkbox_misaligned"),
      sprintf("missing > <b>%s%%</b> in a group, group mean &ge; mean + <b>%s</b> SD (%s); target <b>%s</b>", esc(getp("misaligned_cutoff")),
              esc(getp("intensity_cutoff_sd")), fmt_n(p$intensity_cutoff), esc(getp("misaligned_target"))),
      "A precursor that is high where measured but missing in much of a group is likely misaligned. target=dataset (default) removes the precursor; target=group blanks that group, and limpa is then not available (blanked groups read as low)."),
    c("CV of one group", onoff("filter_cv"), sprintf("group <b>%s</b>, CV &lt; <b>%s%%</b>", esc(getp("filter_cv_group")), esc(getp("filter_cv_value"))),
      "Removes precursors that vary too much within the chosen group (usually SPQC)."),
    c("Precursor quality", onoff("precursor_quality"),
      sprintf("%s: CV &gt; <b>%s%%</b>, intensity &ge; %s, values &ge; %s", if (identical(getp("precursor_quality_method"), "profile")) paste0("<b>profile</b> (residual to the ", feat(p, FALSE), " profile)") else "app (share of the protein signal)",
              esc(getp("precursor_quality_sd")), esc(getp("precursor_quality_intensity")), esc(getp("precursor_quality_min"))),
      if (is_peptide_output(p)) "Peptide output: within each peptide (3+ charge states), removes charge states that don't follow the peptide (profile method; the app's method does nothing for peptides). Phosphopeptides of one protein are never compared with each other."
      else "Within each protein (3+ precursors), removes precursors that don't follow the protein (profile) or whose share of the protein signal varies too much (app)."),
    c("SPQC ratio", onoff("precursor_spqc_ratio"),
      sprintf("SPQC/sample mean outside 100 &plusmn; <b>%s%%</b>, SPQC mean &gt; %s", esc(getp("precursor_spqc_accuracy")), esc(getp("precursor_spqc_intensity"))),
      "The pool should sit near the sample average; precursors where it doesn't are removed."),
    c("PTM only", onoff("filter_ptm"), esc(p$ptm_grep %||% ""), "Keeps only sequences matching the PTM pattern (PTM datasets)."))
  set_df <- as.data.frame(do.call(rbind, set_rows), stringsAsFactors = FALSE); names(set_df) <- c("Filter", "Status", "Setting", "What it does")

  # --- status + decision ---
  status <- if (applied) "Applied" else "Waiting for approval of the settings"
  status_class <- if (applied) "done" else "waiting"
  dec <- read_decisions(run_dir, "05_filter"); last_dec <- if (length(dec)) dec[[length(dec)]] else NULL
  decision_html <- paste0('<section><h2>', if (applied) "Settings applied" else "Decision needed: approve these settings", '</h2>',
    html_table(set_df),
    { mn <- unlist(summ$cutoff$misaligned_notes); if (length(mn)) paste0('<p class="missing">', paste(esc(mn), collapse = "<br/>"), '</p>') else "" },
    if (applied && !is.null(last_dec)) sprintf('<p class="note">Recorded %s%s. %s</p>', esc(last_dec$time),
       if (!is.null(last_dec$analyst)) paste0(" by ", esc(last_dec$analyst)) else "", esc(last_dec$note %||% "")) else "",
    if (!applied) '<p class="note">Use the explorers below to see what each filter would remove at other settings, then approve or change the settings in the Claude session.</p>' else "",
    '</section>')

  # --- applied results ---
  results_html <- ""; result_data <- NULL
  if (applied) {
    cb <- summ$counts_before; cf <- summ$counts_final
    casc <- summ$cascade
    stage_lab <- c(ptm_isolate = "PTM only", min_measured_all = "Minimum measured", min_group_percent = "X% of a group",
                   misaligned = "Misaligned", reapply_minimum = "Re-apply minimum", cv = "CV of group",
                   precursor_quality = "Precursor quality", spqc_ratio = "SPQC ratio", remove_duplicates = "Duplicates")
    casc_df <- data.frame(Stage = vapply(casc, function(x) stage_lab[[x$stage]] %||% x$stage, ""),
                          Precursors = vapply(casc, function(x) as.numeric(x$precursors_removed), 0),
                          Peptides = vapply(casc, function(x) as.numeric(x$peptides_removed), 0),
                          Proteins = vapply(casc, function(x) as.numeric(x$proteins_removed), 0), stringsAsFactors = FALSE)
    miss <- unlist(summ$missing_after_filter_per_sample); missp <- unlist(summ$missing_after_filter_pct)
    grp_of <- rep(groups$Group, as.integer(groups$Count))
    result_data <- list(stages = casc_df$Stage, prec = casc_df$Precursors, pep = casc_df$Peptides, prot = casc_df$Proteins,
                        samples = labels, miss = as.numeric(miss), missp = as.numeric(missp), group = grp_of, groups = groups$Group)
    nk <- summ$noise_keep
    warn <- unlist(summ$warnings)
    results_html <- paste0('<section><h2>Result</h2><div class="tiles">',
      sprintf('<div class="tile"><div class="v">%s &rarr; %s</div><div class="l">precursors (-%s)</div></div>', fmt_n(cb$precursors), fmt_n(cf$precursors), fmt_n(cb$precursors - cf$precursors)),
      sprintf('<div class="tile"><div class="v">%s &rarr; %s</div><div class="l">peptides</div></div>', fmt_n(cb$peptides), fmt_n(cf$peptides)),
      sprintf('<div class="tile"><div class="v">%s &rarr; %s</div><div class="l">proteins</div></div>', fmt_n(cb$proteins), fmt_n(cf$proteins)),
      if (!is.null(nk)) sprintf('<div class="tile"><div class="v">%s</div><div class="l">noise values put back (noise_keep), in %s precursors</div></div>', fmt_n(nk$values_restored), fmt_n(nk$precursors_with_restored_values)) else "",
      '</div>',
      if (length(warn)) paste0('<p class="missing">', paste(esc(warn), collapse = "<br/>"), '</p>') else "",
      '<h2>Removed at each stage</h2><div id="cascade" style="height:340px" role="img" aria-label="Precursors removed at each filter stage"></div>',
      '<details><summary>Table view</summary>', html_table(transform(casc_df, Precursors = fmt_n(Precursors), Peptides = fmt_n(Peptides), Proteins = fmt_n(Proteins)), 2:4), '</details>',
      '<h2>Missing values per sample after filtering</h2><p class="note">These are what the Impute step fills.</p>',
      '<div id="missing" style="height:320px" role="img" aria-label="Percent missing values per sample after filtering"></div>',
      '<details><summary>Table view</summary>', html_table(data.frame(Sample = labels, Group = grp_of, Missing = fmt_n(miss), Percent = fmt_pct(missp)), 3:4), '</details>',
      '</section>')
  }

  qc <- p$qc_path
  lmean <- as.numeric(p$intensity_mean); lsd <- as.numeric(p$intensity_sd)
  ref <- list()
  if (is.finite(lmean) && is.finite(lsd)) ref <- list(`mean` = list(x = lmean, color = "--text-primary", dash = "solid"), `-1 SD` = list(x = lmean - lsd), `+1 SD` = list(x = lmean + lsd))
  if (is.finite(as.numeric(p$intensity_cutoff))) ref$`misaligned cutoff` <- list(x = log2(as.numeric(p$intensity_cutoff)), color = "--series-3", dash = "dash")
  nb <- suppressWarnings(as.numeric(p$noise_baseline_value)); if (identical(p$noise_type, "fixed") && is.finite(nb) && nb > 1) ref$`noise baseline` <- list(x = log2(nb), color = "--series-2", dash = "dash")
  if (identical(p$noise_type, "dynamic") && is.finite(as.numeric(p$noise_inflection))) ref$`noise inflection` <- list(x = log2(as.numeric(p$noise_inflection)), color = "--series-2", dash = "dash")
  charts <- if (applied) {
    fm <- as.matrix(read_table("precursor_filter", db_path)[, (info + 1):(info + n_s)])
    paste0('<p class="note">After filtering', if (!is.null(summ$noise_keep)) " (noise values put back for kept precursors, noise_keep)" else "", '.</p>', sample_charts(fm, labels, grp_of_all, "filt", ref))
  } else paste0('<p class="note">The filter&apos;s input (after the noise step). The misaligned cutoff line is the intensity the misaligned filter uses.</p>', sample_charts(m, labels, grp_of_all, "filt", ref))
  graphs <- paste0('<section><h2>Graphs</h2>', charts,
    '<details><summary>The app&apos;s original plots</summary><div class="grid2">',
    embed_png(file.path(qc, "Precursor_NoiseFiltered_Histogram.png"), "Noise filtered histogram", "App plot: histogram of the filter's input."),
    if (applied) paste0(
      embed_png(file.path(qc, "Precursor_Start_barplot.png"), "Start bar plot", "Start: total intensity per sample."),
      embed_png(file.path(qc, "Precursor_Filter_barplot.png"), "Filter bar plot", "Filtered: total intensity per sample."),
      embed_png(file.path(qc, "Precursor_Start_boxplot.png"), "Start box plot", "Start: log2 intensity per sample."),
      embed_png(file.path(qc, "Precursor_Filter_boxplot.png"), "Filter box plot", "Filtered: log2 intensity per sample.")) else "",
    '</div></details></section>')

  explorers <- list(
    list(id = "wi_min", title = "Minimum measured, all samples", xlab = "minimum values measured", unit = "",
         note = "Precursors measured in fewer samples than this are removed.", d = wi_min, on = TRUE),
    list(id = "wi_pct", title = "Measured in X% of a group", xlab = "X (% of the best group's samples)", unit = "%",
         note = "Precursors whose best group has fewer than X% measured are removed.", d = wi_pct, on = flag("filter_x_percent")),
    list(id = "wi_mis", title = "Misaligned", xlab = "missing % in a group above", unit = "%",
         note = sprintf("Precursors with a group that is more than this %% missing (but not 100%%) and has a group mean at or above the intensity cutoff. Counts are for target=dataset (whole precursor removed). Choose the intensity cutoff in the list."),
         d = wi_mis, on = flag("checkbox_misaligned")),
    list(id = "wi_cv", title = "CV of one group", xlab = "CV cutoff (%)", unit = "%",
         note = "Precursors with CV at or above this in the chosen group are removed.", d = wi_cv, on = flag("filter_cv")))
  if (!is.null(wi_spqc)) explorers[[length(explorers) + 1]] <- list(id = "wi_spqc", title = "SPQC ratio", xlab = "allowed deviation from 100% (+-%)", unit = "%",
    note = "Precursors whose SPQC mean differs from the sample mean by more than this are removed (when the SPQC mean is above the intensity in the list).",
    d = wi_spqc, on = flag("precursor_spqc_ratio"))
  ex_html <- paste0('<section><h2>What-if explorers</h2><p class="note">Each panel counts the precursors one filter would remove <b>on its own</b>, ',
    'applied to the filter&apos;s input (', fmt_n(n_rows), ' precursors). In the real run the filters run in sequence, so the stage counts above ',
    'are smaller. Orange dashed line = current setting; the black dot follows the slider. Precursor quality has no explorer.</p><div class="grid2">',
    paste(vapply(explorers, function(e) sprintf(
      '<div class="panel"><h3>%s <span class="status %s">%s</span></h3><p class="note">%s</p><div class="control">%s<input type="range" id="%s_r" min="0" max="%d" step="1"/><span class="val" id="%s_v"></span></div><div id="%s_p" style="height:240px" role="img" aria-label="%s"></div></div>',
      esc(e$title), if (e$on) "done" else "off", if (e$on) "on" else "off", esc(e$note),
      if (length(e$d$series) > 1) sprintf('<select id="%s_s">%s</select>', e$id, paste0('<option', ifelse(names(e$d$series) == (e$d$current_series %||% ""), " selected", ""), '>', esc(names(e$d$series)), '</option>', collapse = "")) else "",
      e$id, length(e$d$x) - 1, e$id, e$id, esc(e$title)), ""), collapse = ""),
    '</div></section>')

  js_data <- list(n_rows = n_rows, explorers = lapply(explorers, function(e) list(id = e$id, xlab = e$xlab, unit = e$unit, x = e$d$x, series = e$d$series, current = e$d$current)),
                  result = result_data, group_colors = c("--series-1", "--series-2", "--series-3", "--series-4", "--series-5", "--series-6"))
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA), ';', FILTER_JS, '</script>')
  extra_css <- '<style>.panel { border:1px solid var(--rule); border-radius:8px; padding:10px 14px; } .panel h3 { font-size:15px; margin:0 0 4px; }
.status.off { background:var(--grid); color:var(--text-secondary); } .viz-root { --series-5:#e87ba4; --series-6:#008300; }
.control select { font:inherit; font-size:13px; } .control input[type=range] { width:220px; }</style>'
  body <- paste0(decision_html, results_html, graphs, ex_html, js)
  write_review_page(run_dir, "05_filter", "Step 5 - Filter", status, status_class, body, extra_head = extra_css)
}

FILTER_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const base = () => ({ margin: { l: 60, r: 16, t: 8, b: 44 }, paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"),
  font: { color: C("--text-secondary"), size: 12 }, showlegend: false });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
D.explorers.forEach(e => {
  const r = document.getElementById(e.id + "_r"), sel = document.getElementById(e.id + "_s"), v = document.getElementById(e.id + "_v");
  const key = () => sel ? sel.value : Object.keys(e.series)[0];
  const ci = Math.max(0, e.x.findIndex(x => x >= e.current)); r.value = ci;
  function draw() {
    const y = e.series[key()], i = +r.value, xi = e.x[i];
    const exp = { x: [xi], y: [y[i]], mode: "markers", marker: { size: 9, color: C("--text-primary"), line: { width: 2, color: C("--surface-1") } },
      hovertemplate: "explore " + xi + e.unit + "<br>%{y:,} removed<extra></extra>" };
    const curShape = { type: "line", x0: e.current, x1: e.current, y0: 0, y1: 1, yref: "paper", line: { color: C("--ref-dynamic"), width: 2, dash: "dash" } };
    const curAnn = { x: e.current, y: 1, yref: "paper", text: "current " + e.current + e.unit + ": " + fmt(y[ci]), showarrow: false,
      xanchor: e.current > (e.x[0] + e.x[e.x.length - 1]) / 2 ? "right" : "left", xshift: e.current > (e.x[0] + e.x[e.x.length - 1]) / 2 ? -4 : 4,
      yanchor: "bottom", font: { size: 12, color: C("--text-primary") } };
    const line = { x: e.x, y: y, mode: "lines", line: { shape: "hv", width: 2, color: C("--series-1") },
      hovertemplate: "%{x}" + e.unit + ": %{y:,} precursors removed<extra></extra>" };
    Plotly.react(e.id + "_p", [line, exp], Object.assign(base(), { margin: { l: 60, r: 16, t: 24, b: 44 }, shapes: [curShape], annotations: [curAnn], xaxis: ax(e.xlab), yaxis: Object.assign(ax("precursors removed"), { rangemode: "tozero" }) }), { displayModeBar: false, responsive: true });
    v.textContent = xi + e.unit + ": " + fmt(y[i]) + " removed (" + (100 * y[i] / D.n_rows).toFixed(2) + "%)";
  }
  r.addEventListener("input", draw); if (sel) sel.addEventListener("change", draw); draw();
});
if (D.result) {
  const R = D.result;
  Plotly.newPlot("cascade", [{ type: "bar", orientation: "h", y: R.stages, x: R.prec, marker: { color: C("--series-1") },
    customdata: R.stages.map((_, i) => [R.pep[i], R.prot[i]]),
    hovertemplate: "%{y}<br>%{x:,} precursors<br>%{customdata[0]:,} peptides<br>%{customdata[1]:,} proteins<extra></extra>",
    text: R.prec.map(fmt), textposition: "outside", cliponaxis: false }],
    Object.assign(base(), { margin: { l: 130, r: 60, t: 8, b: 40 }, xaxis: ax("precursors removed"), yaxis: { autorange: "reversed" }, bargap: 0.35 }), { displaylogo: false, responsive: true });
  const gc = {}; R.groups.forEach((g, i) => gc[g] = C(D.group_colors[i % D.group_colors.length]));
  const traces = R.groups.map(g => { const idx = R.group.map((x, i) => x === g ? i : -1).filter(i => i >= 0);
    return { type: "bar", name: g, x: idx.map(i => R.samples[i]), y: idx.map(i => R.missp[i]), marker: { color: gc[g] },
      customdata: idx.map(i => R.miss[i]), hovertemplate: "%{x}<br>%{y:.2f}% missing (%{customdata:,} values)<extra>" + g + "</extra>" }; });
  Plotly.newPlot("missing", traces, Object.assign(base(), { showlegend: true, legend: { orientation: "h", y: 1.12 }, bargap: 0.3,
    xaxis: { categoryorder: "array", categoryarray: R.samples }, yaxis: Object.assign(ax("% missing"), { rangemode: "tozero" }) }), { displaylogo: false, responsive: true });
}
'
