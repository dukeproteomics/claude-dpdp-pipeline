# Review page for Step 9 - QC (spike-ins). A review step: it doesn't change the data. The analyst accepts it (optionally
# noting samples) before Stats.
# Lab rules (user, 2026-11-04):
#   casein (spike at design levels): the sample sheet is correct; a sample that doesn't match its level is a sample problem.
#   ADH (constant, added by volume before injection on the EvoSep tip): judge it WITHOUT normalization (impute table).
# Shows: spike measured vs expected per sample with a tolerance slider and the pool check (SPQC vs the mean of the study
# samples), measured precursors per sample, constant spike per sample (un-normalized by default) with CV and max/min, app plots.

review_qc <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  sd_dir <- file.path(run_dir, "steps", "09_qc")
  if (!file.exists(file.path(sd_dir, "summary.json"))) stop("Run 09_qc.R first.")
  summ <- fromJSON(file.path(sd_dir, "summary.json"), simplifyVector = FALSE)
  design <- read_table("design", db_path); design <- design[order(design$Raw_Order), ]
  groups <- read_table("sample_groups", db_path)
  grp_of <- rep(groups$Group, as.integer(groups$Count))
  labels <- design$Label
  pool <- if (!is.null(p$comp_spqc) && p$comp_spqc %in% groups$Group) p$comp_spqc else if ("SPQC" %in% groups$Group) "SPQC" else NA
  dec <- read_decisions(run_dir, "09_qc")
  accepted <- Filter(function(r) identical(r$settings$accepted, TRUE), dec)
  last_acc <- if (length(accepted)) accepted[[length(accepted)]] else NULL
  # accepted only counts if it came after the latest QC run
  if (!is.null(last_acc) && !is.null(summ$generated) && last_acc$time < summ$generated) last_acc <- NULL

  # measured precursors per sample for the spike and constant accessions (from the filtered table, before imputation)
  filt <- read_table("precursor_filter", db_path)
  fm <- as.matrix(filt[, (ncol(filt) - n_s + 1):ncol(filt)])
  prec_counts <- function(accs) { r <- filt$Accession %in% accs; list(total = sum(r), measured = as.integer(colSums(!is.na(fm[r, , drop = FALSE])))) }

  st <- summ$settings
  res <- setNames(summ$results, vapply(summ$results, function(r) r$norm, ""))
  spike <- lapply(res, function(r) if (length(r$spike$accessions)) list(
    intensity = vapply(r$spike$per_sample, function(x) as.numeric(x$intensity), 0),
    level = vapply(r$spike$per_sample, function(x) as.numeric(x$level), 0),
    by_level = r$spike$by_level) else NULL)
  constant <- lapply(res, function(r) lapply(r$constant, function(c) list(acc = c$accession, name = c$name,
    values = as.numeric(unlist(c$values)), cv = c$cv_pct_all_samples, mm = c$max_over_min)))
  norms <- names(res)
  default_spike_norm <- norms[1]; default_const_norm <- if ("impute" %in% norms) "impute" else norms[1]

  status <- if (!is.null(last_acc)) "Accepted" else "Waiting for the analyst's review"
  status_class <- if (!is.null(last_acc)) "done" else "waiting"
  decision_html <- paste0('<section><h2>', if (!is.null(last_acc)) "Review recorded" else "Review needed", '</h2>',
    '<p>QC does not change the data. Check the spike-ins below, then accept the QC in the Claude session (optionally noting samples) before Stats.</p>',
    '<ul class="options"><li><b>Casein</b> (', esc(paste(unlist(st$spike_accession), collapse = ", ")), ') is spiked at the levels in the sample sheet (',
    esc(grep("Level", names(design), value = TRUE)[1]), '). The sheet is taken as correct: a sample that doesn&apos;t match its level is a <b>sample</b> problem.</li>',
    '<li><b>ADH</b> (', esc(paste(unlist(st$constant_accession), collapse = ", ")), ') is added by volume before injection (EvoSep tip). Judge it <b>without normalization</b>; ',
    'normalization can make it less flat, which is expected.</li></ul>',
    if (!is.null(last_acc)) sprintf('<p class="note">Accepted %s%s.%s %s</p>', esc(last_acc$time), if (!is.null(last_acc$analyst)) paste0(" by ", esc(last_acc$analyst)) else "",
      if (length(unlist(last_acc$settings$flagged_samples))) paste0(" Flagged samples: <b>", esc(paste(unlist(last_acc$settings$flagged_samples), collapse = ", ")), "</b>.") else "",
      esc(last_acc$note %||% "")) else "",
    '</section>')

  spike_html <- if (length(Filter(Negate(is.null), spike))) {
    pc <- prec_counts(unlist(st$spike_accession))
    paste0('<section><h2>Casein spike: measured vs design level</h2>',
      '<p class="note">Bars: summed spike intensity per sample. Diamonds: expected from the design level, scaled to the mean of the lowest-level samples. ',
      'Pool samples (', esc(pool %||% ""), ') are compared with the mean of the study samples instead of a level. A sample is flagged when measured / expected is outside the tolerance.</p>',
      '<div class="control"><label for="snorm">Normalization</label><select id="snorm">', paste0('<option', ifelse(norms == default_spike_norm, " selected", ""), '>', esc(norms), '</option>', collapse = ""), '</select>',
      '<label for="tol">Tolerance</label><input type="range" id="tol" min="10" max="200" step="5" value="50"/><span class="val" id="tolv"></span></div>',
      '<div class="tiles" id="spiketiles"></div>',
      '<div id="spikeplot" style="height:360px" role="img" aria-label="Spike intensity per sample with expected levels"></div>',
      '<details open><summary>Per sample (table view)</summary><div id="spiketable"></div></details>',
      sprintf('<p class="note">%d spike precursors after filtering. Few measured precursors means the spike is near the detection limit there and imputation fills much of it.</p>', pc$total),
      '</section>')
  } else '<section><h2>Casein spike</h2><p>Spike accessions not found.</p></section>'

  const_html <- if (length(constant[[1]])) paste0('<section><h2>ADH: constant spike per sample</h2>',
    '<p class="note">Should be flat across samples. Default view is the un-normalized (impute) table, the right view for a spike added by volume.</p>',
    '<div class="control"><label for="cnorm">Normalization</label><select id="cnorm">', paste0('<option', ifelse(norms == default_const_norm, " selected", ""), '>', esc(norms), '</option>', collapse = ""), '</select></div>',
    '<div class="tiles" id="consttiles"></div><div class="grid2" id="constgrid">',
    paste0('<div id="const_', seq_along(constant[[1]]), '" style="height:300px" role="img" aria-label="Constant spike per sample"></div>', collapse = ""), '</div></section>') else ""

  qc <- p$qc_path
  app_plots <- unlist(summ$plots)
  graphs <- paste0('<section><h2>App plots</h2><div class="grid2">', paste(vapply(app_plots, function(f) embed_png(f, basename(f), sub("_barplot\\.png$", "", basename(f))), ""), collapse = ""), '</div></section>')

  pcs <- prec_counts(unlist(st$spike_accession))
  pcc <- lapply(unlist(st$constant_accession), function(a) prec_counts(a)); names(pcc) <- unlist(st$constant_accession)
  js_data <- list(samples = labels, groups = grp_of, pool = pool, spike = spike, constant = constant, spike_prec = pcs, const_prec = pcc)
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA, null = "null"), ';', QC_JS, '</script>')
  extra_css <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; } .control select { font:inherit; font-size:13px; }
td.flag { color:var(--warn-ink); font-weight:600; } .control input[type=range] { width:200px; }</style>'
  body <- paste0(decision_html, spike_html, const_html, graphs, js)
  write_review_page(run_dir, "09_qc", "Step 9 - QC", status, status_class, body, extra_head = extra_css)
}

QC_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
const cfg = { displaylogo: false, responsive: true };
const GC = {}; [...new Set(D.groups)].forEach((g, i) => GC[g] = ["--series-1", "--series-2", "--series-3"][i % 3]);
function spikeDraw() {
  const n = document.getElementById("snorm").value, S = D.spike[n], tol = +document.getElementById("tol").value;
  document.getElementById("tolv").textContent = "±" + tol + "%";
  const minL = Math.min(...S.level), lowIdx = S.level.map((l, i) => l === minL ? i : -1).filter(i => i >= 0);
  const ref = lowIdx.reduce((s, i) => s + S.intensity[i], 0) / lowIdx.length;
  const isPool = D.groups.map(g => g === D.pool), stIdx = D.groups.map((g, i) => g !== D.pool ? i : -1).filter(i => i >= 0);
  const studyMean = stIdx.reduce((s, i) => s + S.intensity[i], 0) / Math.max(1, stIdx.length);
  // study samples: expected from the design level; pool samples: expected = mean of the study samples (what pooling gives)
  const exp = S.level.map((l, i) => isPool[i] && D.pool ? studyMean : ref * l / minL), ratio = S.intensity.map((v, i) => v / exp[i]);
  const flag = ratio.map(r => r > 1 + tol / 100 || r < 1 / (1 + tol / 100));
  const bars = [...new Set(D.groups)].map(g => { const idx = D.groups.map((x, i) => x === g ? i : -1).filter(i => i >= 0);
    return { type: "bar", name: g, x: idx.map(i => D.samples[i]), y: idx.map(i => S.intensity[i]), marker: { color: C(GC[g]) },
      customdata: idx.map(i => [S.level[i], fmt(exp[i]), ratio[i].toFixed(2)]),
      hovertemplate: "%{x}<br>measured %{y:,.0f}<br>level %{customdata[0]}, expected %{customdata[1]}<br>measured/expected %{customdata[2]}<extra></extra>" }; });
  bars.push({ type: "scatter", mode: "markers", name: "expected (level; pool = study mean)", x: D.samples, y: exp, marker: { symbol: "diamond", size: 11, color: C("--text-primary"), line: { width: 2, color: C("--surface-1") } },
    hovertemplate: "%{x}: expected %{y:,.0f}<extra></extra>" });
  Plotly.react("spikeplot", bars, Object.assign(base(), { margin: { l: 70, r: 10, t: 10, b: 50 }, showlegend: true, legend: { orientation: "h", y: 1.12 },
    xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: Object.assign(ax("summed spike intensity (log scale)"), { type: "log" }) }), cfg);
  // pool check
  let tiles = "";
  if (D.pool) { const st = D.groups.map((g, i) => g !== D.pool ? i : -1).filter(i => i >= 0), pl = D.groups.map((g, i) => g === D.pool ? i : -1).filter(i => i >= 0);
    const mS = st.reduce((s, i) => s + S.intensity[i], 0) / st.length, mP = pl.reduce((s, i) => s + S.intensity[i], 0) / pl.length;
    tiles += `<div class="tile"><div class="v">${fmt(mP)}</div><div class="l">${D.pool} mean</div></div><div class="tile"><div class="v">${fmt(mS)}</div><div class="l">mean of the study samples (what a pool should show)</div></div><div class="tile"><div class="v">${(mP / mS).toFixed(2)}</div><div class="l">${D.pool} / study mean</div></div>`; }
  tiles += `<div class="tile"><div class="v">${flag.filter(f => f).length} of ${flag.length}</div><div class="l">samples outside ±${tol}%</div></div>`;
  document.getElementById("spiketiles").innerHTML = tiles;
  document.getElementById("spiketable").innerHTML = "<table><thead><tr><th>Sample</th><th>Group</th><th class=num>level</th><th class=num>measured</th><th class=num>expected</th><th class=num>measured / expected</th><th class=num>precursors measured</th><th></th></tr></thead><tbody>" +
    D.samples.map((s, i) => `<tr><td>${s}</td><td>${D.groups[i]}</td><td class=num>${isPool[i] && D.pool ? "pool" : S.level[i]}</td><td class=num>${fmt(S.intensity[i])}</td><td class=num>${fmt(exp[i])}</td><td class=num>${ratio[i].toFixed(2)}</td><td class=num>${D.spike_prec.measured[i]} of ${D.spike_prec.total}</td><td class="${flag[i] ? "flag" : ""}">${flag[i] ? (ratio[i] > 1 ? "high" : "low") + (isPool[i] && D.pool ? " vs study mean" : " - check sample") : ""}</td></tr>`).join("") + "</tbody></table>";
}
if (document.getElementById("snorm")) { document.getElementById("snorm").addEventListener("change", spikeDraw); document.getElementById("tol").addEventListener("input", spikeDraw); spikeDraw(); }
function constDraw() {
  const n = document.getElementById("cnorm").value, L = D.constant[n];
  let tiles = "";
  L.forEach((c, k) => {
    const mean = c.values.reduce((s, v) => s + v, 0) / c.values.length, pr = D.const_prec[c.acc];
    tiles += `<div class="tile"><div class="v">${c.cv}%</div><div class="l">${c.acc} ${c.name} CV (max/min ${c.mm})</div></div>`;
    const tr = [...new Set(D.groups)].map(g => { const idx = D.groups.map((x, i) => x === g ? i : -1).filter(i => i >= 0);
      return { type: "bar", name: g, x: idx.map(i => D.samples[i]), y: idx.map(i => c.values[i]), marker: { color: C(GC[g]) },
        customdata: idx.map(i => [(100 * c.values[i] / mean - 100).toFixed(1), pr ? pr.measured[i] + " of " + pr.total : ""]),
        hovertemplate: "%{x}: %{y:,.0f}<br>%{customdata[0]}% vs mean<br>precursors measured %{customdata[1]}<extra></extra>" }; });
    Plotly.react("const_" + (k + 1), tr, Object.assign(base(), { title: { text: c.acc + " " + c.name + " (" + n + ")", font: { size: 13, color: C("--text-primary") } },
      margin: { l: 70, r: 10, t: 36, b: 50 }, showlegend: k === 0, legend: { orientation: "h", y: -0.2 },
      shapes: [{ type: "line", x0: 0, x1: 1, xref: "paper", y0: mean, y1: mean, line: { color: C("--text-muted"), width: 1, dash: "dot" } }],
      xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: Object.assign(ax("protein intensity"), { rangemode: "tozero" }) }), { displayModeBar: false, responsive: true });
  });
  document.getElementById("consttiles").innerHTML = tiles;
}
if (document.getElementById("cnorm")) { document.getElementById("cnorm").addEventListener("change", constDraw); constDraw(); }
'
