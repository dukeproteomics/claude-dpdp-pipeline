# Review page for Step 8 - Rollup (precursor -> protein).
# Decision: rollup_method (only sum and iq_maxlfq work in the app) and maxlfq_scale.
# Shows: precursors per protein, a sum vs iq_maxlfq comparison computed here (protein CV per group, same CV formula as
# the app's percentCV_gw on the raw scale), a protein explorer (precursor profiles with both rollups), and after apply the
# app's CV bar plot, protein bar plots, and the CV summary.

ROLLUP_WORKING <- c("sum", "iq_maxlfq")

review_rollup <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  design <- read_table("design", db_path)
  groups <- read_table("sample_groups", db_path)
  labels <- design$Label
  grp_of <- rep(groups$Group, as.integer(groups$Count))
  gcols <- lapply(seq_len(nrow(groups)), function(i) as.integer(groups$start[i]):as.integer(groups$end[i]))
  norms <- trimws(strsplit(p$norm_type, ",")[[1]])
  sd_dir <- file.path(run_dir, "steps", "08_rollup")
  summ <- if (file.exists(file.path(sd_dir, "summary.json"))) fromJSON(file.path(sd_dir, "summary.json"), simplifyVector = FALSE) else NULL
  applied <- !is.null(summ) && identical(summ$action %||% "apply", "apply")
  s <- summ$settings %||% list()
  method <- s$rollup_method %||% p$rollup_method; scale <- as_flag(s$maxlfq_scale %||% p$maxlfq_scale)

  cvs <- function(m) sapply(gcols, function(cc) { x <- m[, cc, drop = FALSE]; signif(100 * apply(x, 1, stats::sd) / rowMeans(x), 3) })
  ecdf_pts <- function(v) { v <- sort(v[is.finite(v)]); g <- c(seq(0, 50, by = 0.5), seq(51, 100, by = 1)); list(x = g, y = round(vapply(g, function(t) mean(v <= t), 0) * 100, 2)) }

  comp <- list(); explorer <- NULL; ppp <- NULL
  for (n in norms) {
    tab <- paste0("precursor_impute_", n)
    if (!table_exists(tab, db_path)) next
    df <- read_table(tab, db_path)
    x <- as.matrix(df[, (ncol(df) - n_s + 1):ncol(df)])
    acc <- df$Accession
    summed <- rowsum(x, acc, reorder = TRUE, na.rm = TRUE)
    lx <- log2(x); lx[!is.finite(lx)] <- NA
    long <- data.frame(protein_list = rep(acc, n_s), sample_list = rep(colnames(x), each = nrow(x)), quant = as.vector(lx), id = rep(seq_len(nrow(x)), n_s))
    long <- long[!is.na(long$quant), ]
    est <- suppressMessages(utils::capture.output(fl <- iq::fast_MaxLFQ(long), type = "output")); lfq <- 2^fl$estimate
    lfq <- lfq[rownames(summed), , drop = FALSE]
    lfq_s <- lfq * (rowMeans(summed) / rowMeans(lfq, na.rm = TRUE))
    cv_sum <- cvs(summed); cv_lfq <- cvs(lfq)
    comp[[n]] <- list(
      sum = lapply(seq_along(gcols), function(i) ecdf_pts(cv_sum[, i])), iq_maxlfq = lapply(seq_along(gcols), function(i) ecdf_pts(cv_lfq[, i])),
      median = list(sum = round(apply(cv_sum, 2, median, na.rm = TRUE), 1), iq_maxlfq = round(apply(cv_lfq, 2, median, na.rm = TRUE), 1)),
      mean = list(sum = round(colMeans(cv_sum, na.rm = TRUE), 1), iq_maxlfq = round(colMeans(cv_lfq, na.rm = TRUE), 1)),
      under20 = list(sum = round(colMeans(cv_sum < 20, na.rm = TRUE) * 100, 1), iq_maxlfq = round(colMeans(cv_lfq < 20, na.rm = TRUE) * 100, 1)))
    if (is.null(explorer)) {
      npre <- table(acc)
      ppp <- list(k = 1:20, n = vapply(1:20, function(k) if (k < 20) sum(npre == k) else sum(npre >= 20), 0), proteins = length(npre))
      spikes <- intersect(c("P00330", "P00331", "P02662", "P02663", "P02666"), names(npre))
      set.seed(1); mid <- names(npre)[npre >= 2 & npre <= 15]
      pick <- unique(c(spikes, names(sort(npre, decreasing = TRUE))[1:8], sample(mid, min(25, length(mid)))))
      genes <- setNames(df$Genes[match(pick, acc)], pick); desc <- setNames(df$Description[match(pick, acc)], pick)
      explorer <- list(norm = n, proteins = lapply(pick, function(a) {
        r <- which(acc == a)
        list(acc = a, label = sprintf("%s %s (%d precursors)%s", a, genes[[a]] %||% "", length(r), if (a %in% spikes) " - spike-in" else ""),
             desc = desc[[a]], prec = unname(lapply(r, function(i) round(lx[i, ], 3))),
             sum = round(log2(summed[a, ]), 3), lfq = round(log2(if (scale) lfq_s[a, ] else lfq[a, ]), 3))
      }))
    }
  }

  # --- status + decision ---
  status <- if (applied) sprintf("Applied: %s", method) else "Waiting for the analyst's decision"
  status_class <- if (applied) "done" else "waiting"
  dec <- read_decisions(run_dir, "08_rollup"); last_dec <- if (length(dec)) dec[[length(dec)]] else NULL
  options_html <- '<ul class="options">
<li><b>iq_maxlfq</b> - MaxLFQ (iq package): protein value from the median of precursor ratios between samples. Robust to a precursor that is off
in some samples. <b>maxlfq_scale</b> on: each protein is rescaled to the level of its summed intensity (does not change CVs or fold changes).</li>
<li><b>sum</b> - add the precursor intensities. Dominated by the most intense precursors.</li>
<li>Not available: <b>median</b>, <b>mean</b>, <b>median_polish</b>, <b>topn</b> - the app&apos;s functions fail (tested on 11488).</li></ul>'
  decision_html <- paste0('<section><h2>', if (applied) "Decision" else "Decision needed", '</h2>',
    if (applied) paste0('<div class="tiles">',
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">rollup_method</div></div>', esc(method)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">maxlfq_scale</div></div>', if (scale) "on" else "off"),
      if (!is.null(summ$tables)) paste(vapply(summ$tables, function(t) sprintf('<div class="tile"><div class="v">%s</div><div class="l">proteins (%s)</div></div>', fmt_n(t$proteins), esc(t$norm)), ""), collapse = "") else "",
      '</div>', if (!is.null(last_dec)) sprintf('<p class="note">Recorded %s%s. %s</p>', esc(last_dec$time), if (!is.null(last_dec$analyst)) paste0(" by ", esc(last_dec$analyst)) else "", esc(last_dec$note %||% "")) else "",
      '<details><summary>Options</summary>', options_html, '</details>') else options_html,
    if (!applied) '<p class="note">Compare the methods below, then give the choice in the Claude session.</p>' else "",
    '</section>')

  applied_html <- ""
  if (applied && !is.null(summ$summary_cv)) {
    scv <- summ$summary_cv
    rows <- lapply(summ$tables, function(t) {
      data.frame(Normalization = t$norm, Group = groups$Group,
                 `Mean CV (app)` = vapply(groups$Group, function(g) { r <- Filter(function(z) z$Group == g, scv)[[1]]; sprintf("%.1f%%", as.numeric(r[[t$norm]])) }, ""),
                 `Median CV` = sprintf("%.1f%%", unlist(t$median_cv_by_group)), `Proteins CV < 20%` = sprintf("%.1f%%", unlist(t$proteins_cv_under_20_pct)), check.names = FALSE)
    })
    applied_html <- paste0('<section><h2>Result: protein CV by group</h2><p class="note">The app reports the mean CV, which a high-CV tail pulls up; the median is the typical protein. ',
      'The pooled QC (SPQC) should be lowest: it measures technical variation only.</p>', html_table(do.call(rbind, rows), 3:5), '</section>')
  }

  ppp_html <- if (!is.null(ppp)) sprintf('<section><h2>Precursors per protein</h2><div class="tiles"><div class="tile"><div class="v">%s</div><div class="l">proteins</div></div>
<div class="tile"><div class="v">%s</div><div class="l">with 1 precursor (%s)</div></div><div class="tile"><div class="v">%s</div><div class="l">with 1 or 2 precursors (%s)</div></div></div>
<p class="note">With one precursor, every method gives the same protein value; differences between methods come from proteins with several precursors.</p>
<div id="pppplot" style="height:260px" role="img" aria-label="Number of proteins by precursor count"></div></section>',
    fmt_n(ppp$proteins), fmt_n(ppp$n[1]), fmt_pct(ppp$n[1] / ppp$proteins * 100, 1), fmt_n(ppp$n[1] + ppp$n[2]), fmt_pct((ppp$n[1] + ppp$n[2]) / ppp$proteins * 100, 1)) else ""

  comp_html <- paste0('<section><h2>sum vs iq_maxlfq: protein CV within each group</h2>',
    '<p class="note">Both rollups calculated here from the imputed precursor table. Each curve is the share of proteins with a CV at or below x. ',
    'A curve further up and left means tighter replicates. Same CV formula as the app.</p>',
    if (length(comp) > 1) paste0('<div class="control"><label for="cnorm">Normalization</label><select id="cnorm">', paste0('<option>', esc(names(comp)), '</option>', collapse = ""), '</select></div>') else "",
    '<div class="grid3" id="cvgrid">', paste0('<div><div id="cv_', seq_along(gcols), '" style="height:300px" role="img" aria-label="CV curves group ', esc(groups$Group), '"></div></div>', collapse = ""), '</div>',
    '<details open><summary>Table view</summary><div id="cvtable"></div></details></section>')

  expl_html <- if (!is.null(explorer)) sprintf('<section><h2>Protein explorer</h2>
<p class="note">Grey lines: the protein&apos;s precursors (log2, %s). Coloured: the protein value from sum and from iq_maxlfq%s. Proteins shown: spike-ins, the
proteins with the most precursors, and a random set with 2 to 15 precursors.</p>
<div class="control"><label for="prot">Protein</label><select id="prot">%s</select>
<label for="pview">Show</label><select id="pview"><option value="rel">relative: each line centred on its mean</option><option value="abs">absolute log2</option></select></div>
<div id="protplot" style="height:380px" role="img" aria-label="Precursor and protein values across samples"></div></section>',
    esc(explorer$norm), if (scale) " (scaled to the sum level)" else "",
    paste0('<option value="', seq_along(explorer$proteins) - 1, '">', vapply(explorer$proteins, function(z) esc(z$label), ""), '</option>', collapse = "")) else ""

  qc <- p$qc_path
  graphs <- if (applied) paste0('<section><h2>Graphs</h2><div class="grid2">',
    embed_png(file.path(qc, "CV_barplot.png"), "CV bar plot", "App plot: average protein CV% per group for each normalization."),
    paste(vapply(norms, function(n) embed_png(file.path(qc, paste0("protein_", n, "_barplot.png")), paste("Protein", n), sprintf("App plot: protein total intensity per sample (%s).", n)), ""), collapse = ""),
    '</div></section>') else ""

  js_data <- list(groups = groups$Group, samples = labels, comp = comp, ppp = ppp, explorer = explorer, method = method)
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA), ';', ROLLUP_JS, '</script>')
  extra_css <- '<style>.viz-root { --series-2:#eb6834; } .grid3 { display:grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap:12px; }
.control select { font:inherit; font-size:13px; max-width:560px; }</style>'
  body <- paste0(decision_html, applied_html, comp_html, ppp_html, expl_html, graphs, js)
  write_review_page(run_dir, "08_rollup", "Step 8 - Rollup", status, status_class, body, extra_head = extra_css)
}

ROLLUP_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
const cfg = { displaylogo: false, responsive: true };
const MC = { sum: "--series-1", iq_maxlfq: "--series-2" };
function drawCV() {
  const sel = document.getElementById("cnorm"), n = sel ? sel.value : Object.keys(D.comp)[0], c = D.comp[n];
  D.groups.forEach((g, i) => {
    const tr = ["sum", "iq_maxlfq"].map(m => ({ type: "scatter", mode: "lines", name: m + (m === D.method ? " (current setting)" : ""), x: c[m][i].x, y: c[m][i].y,
      line: { width: 2, color: C(MC[m]), dash: m === "sum" ? "dot" : "solid" }, hovertemplate: m + ": %{y:.1f}% of proteins with CV ≤ %{x}%<extra></extra>" }));
    Plotly.react("cv_" + (i + 1), tr, Object.assign(base(), { title: { text: g, font: { size: 14, color: C("--text-primary") } }, margin: { l: 55, r: 10, t: 36, b: 45 },
      showlegend: i === 0, legend: { x: 0.45, y: 0.1 }, xaxis: Object.assign(ax("protein CV %"), { range: [0, 60] }), yaxis: Object.assign(ax("% of proteins"), { range: [0, 100] }) }), { displayModeBar: false, responsive: true });
  });
  document.getElementById("cvtable").innerHTML = "<table><thead><tr><th>Group</th><th class=num>median CV sum</th><th class=num>median CV iq_maxlfq</th><th class=num>mean CV sum</th><th class=num>mean CV iq_maxlfq</th><th class=num>CV &lt; 20% sum</th><th class=num>CV &lt; 20% iq_maxlfq</th></tr></thead><tbody>" +
    D.groups.map((g, i) => `<tr><td>${g}</td><td class=num>${c.median.sum[i]}%</td><td class=num>${c.median.iq_maxlfq[i]}%</td><td class=num>${c.mean.sum[i]}%</td><td class=num>${c.mean.iq_maxlfq[i]}%</td><td class=num>${c.under20.sum[i]}%</td><td class=num>${c.under20.iq_maxlfq[i]}%</td></tr>`).join("") + "</tbody></table>";
}
const cn = document.getElementById("cnorm"); if (cn) cn.addEventListener("change", drawCV); drawCV();
if (D.ppp) Plotly.newPlot("pppplot", [{ type: "bar", x: D.ppp.k.map(k => k === 20 ? "20+" : String(k)), y: D.ppp.n, marker: { color: C("--series-1") },
  hovertemplate: "%{x} precursors: %{y:,} proteins<extra></extra>" }], Object.assign(base(), { margin: { l: 60, r: 10, t: 10, b: 45 }, bargap: 0.2,
  xaxis: { title: "precursors per protein", type: "category" }, yaxis: Object.assign(ax("proteins"), { rangemode: "tozero" }) }), cfg);
if (D.explorer) {
  const ps = document.getElementById("prot");
  function drawP() {
    const z = D.explorer.proteins[+ps.value], rel = document.getElementById("pview").value === "rel";
    const ctr = a => { if (!rel) return a; const v = a.filter(q => q !== null && isFinite(q)), m = v.reduce((s, q) => s + q, 0) / Math.max(1, v.length); return a.map(q => q === null ? null : q - m); };
    const tr = z.prec.map(ctr).map((y, j) => ({ type: "scatter", mode: "lines", x: D.samples, y: y, line: { width: 1, color: C("--text-muted") }, opacity: 0.6,
      name: "precursors", legendgroup: "p", showlegend: j === 0, hovertemplate: "precursor " + (j + 1) + "<br>%{x}: %{y:.2f}<extra></extra>" }));
    tr.push({ type: "scatter", mode: "lines+markers", name: "sum", x: D.samples, y: ctr(z.sum), line: { width: 2, color: C("--series-1"), dash: "dot" },
      marker: { size: 8, color: C("--series-1") }, hovertemplate: "sum %{x}: %{y:.2f}<extra></extra>" });
    tr.push({ type: "scatter", mode: "lines+markers", name: "iq_maxlfq", x: D.samples, y: ctr(z.lfq), line: { width: 2, color: C("--series-2") },
      marker: { size: 8, color: C("--series-2") }, hovertemplate: "iq_maxlfq %{x}: %{y:.2f}<extra></extra>" });
    Plotly.react("protplot", tr, Object.assign(base(), { margin: { l: 60, r: 10, t: 30, b: 80 }, showlegend: true, legend: { orientation: "h", y: -0.2 },
      title: { text: z.desc, font: { size: 13, color: C("--text-primary") } }, xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: ax(rel ? "log2 intensity minus line mean" : "log2 intensity") }), cfg);
  }
  ps.addEventListener("change", drawP); document.getElementById("pview").addEventListener("change", drawP); drawP();
}
'
