# Review page for Step 6 - Normalize.
# Two decisions: (1) which proteins to leave out when calculating the factors (exclude/include regex),
# (2) the normalization method.
# Shows: the excluded proteins (current regex vs the pipeline default rule), per-sample scaling factors from several
# methods (calculated here from the normalization subset, same formulas as the app), the app's NormData and normalized
# bar plots, and per-sample log2 distributions before vs after the applied method.

NORM_DEFAULT_EXCLUDE <- "trypsin, org=(sus scrofa|bos taurus)|keratin, type [^,]*, org=homo sapiens|(alpha-s1|alpha-s2|beta)-casein, org=bos taurus"

review_normalize <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  design <- read_table("design", db_path)
  groups <- read_table("sample_groups", db_path)
  filt <- read_table("precursor_filter", db_path)
  info <- ncol(filt) - n_s
  labels <- if (!is.null(design$Label) && length(design$Label) == n_s) design$Label else names(filt)[(info + 1):ncol(filt)]
  grp_of <- rep(groups$Group, as.integer(groups$Count))
  sd_dir <- file.path(run_dir, "steps", "06_normalize")
  summ <- if (file.exists(file.path(sd_dir, "summary.json"))) fromJSON(file.path(sd_dir, "summary.json"), simplifyVector = FALSE) else NULL
  applied <- !is.null(summ) && identical(summ$action, "apply")
  s <- summ$settings %||% list()
  excl_mode <- s$norm_exclude %||% p$norm_exclude; excl_grep <- s$norm_exclude_grep %||% p$norm_exclude_grep
  incl_mode <- s$norm_include %||% p$norm_include; incl_grep <- s$norm_include_grep %||% p$norm_include_grep
  normdata <- if (table_exists("precursor_normdata", db_path)) read_table("precursor_normdata", db_path) else filt
  smat <- function(d) { x <- as.matrix(d[, (ncol(d) - n_s + 1):ncol(d)]); x[!is.na(x) & x == 0] <- NA; x }

  # --- excluded proteins: current subset vs the pipeline default rule ---
  in_cur <- !filt$PrecursorId %in% normdata$PrecursorId
  # pipeline default rule: leave out added material; on PTM runs also everything that isn't a PTM peptide
  ptm <- is_ptm_run(p)
  in_def <- grepl(NORM_DEFAULT_EXCLUDE, filt$Description, ignore.case = TRUE) | (ptm & !is_ptm_row(filt, p))
  in_def_mat <- grepl(NORM_DEFAULT_EXCLUDE, filt$Description, ignore.case = TRUE)
  key <- paste(filt$Accession, filt$Description, sep = "\r")
  # the excluded-proteins list shows proteins left out by description (added material); the non-PTM peptides are counted, not listed
  keep_rows <- !ptm | is_ptm_row(filt, p)
  ex <- data.frame(key = key, cur = in_cur & keep_rows, def = in_def_mat & keep_rows, stringsAsFactors = FALSE)
  ex <- ex[ex$cur | ex$def, ]
  agg <- if (nrow(ex)) do.call(rbind, lapply(split(ex, ex$key), function(d) data.frame(key = d$key[1], n = nrow(d), cur = any(d$cur), def = any(d$def)))) else NULL
  same_rule <- identical(excl_mode, "description") && identical(tolower(excl_grep), tolower(NORM_DEFAULT_EXCLUDE)) && (if (ptm) identical(incl_mode, "sequence") && identical(tolower(incl_grep), tolower(p$ptm_grep)) else identical(incl_mode, "not_used"))
  ex_table <- ""
  if (!is.null(agg)) {
    agg$Accession <- sub("\r.*", "", agg$key); agg$Description <- sub(".*\r", "", agg$key)
    agg$Rule <- ifelse(agg$cur & agg$def, "both", ifelse(agg$cur, "current only", "default only"))
    agg <- agg[order(match(agg$Rule, c("current only", "default only", "both")), -agg$n), ]
    show <- data.frame(Accession = esc(agg$Accession), Description = esc(substr(agg$Description, 1, 90)), Precursors = fmt_n(agg$n),
                       `Excluded by` = agg$Rule, check.names = FALSE)
    ex_table <- html_table(show, 3)
  }
  n_cur_prot <- if (is.null(agg)) 0 else sum(agg$cur); n_def_prot <- if (is.null(agg)) 0 else sum(agg$def)
  n_cur_prec <- sum(in_cur); n_def_prec <- sum(in_def)

  # --- scaling factors (multiply the sample by this), current subset and default-rule subset ---
  factor_set <- function(x) {
    tot <- colSums(x, na.rm = TRUE); med <- matrixStats::colMedians(x, na.rm = TRUE)
    sl <- mean(tot) / tot
    x0 <- x; x0[is.na(x0)] <- 0
    tmm <- 1 / edgeR::calcNormFactors(x0, method = "TMM", sumTrim = 0.1)
    xs <- sweep(x0, 2, sl, "*"); sltmm <- sl / edgeR::calcNormFactors(xs, method = "TMM", sumTrim = 0.1)
    list(SL = round(sl, 4), `SL + TMM (estimate)` = round(sltmm, 4), TMM = round(tmm, 4), Median = round(mean(med) / med, 4))
  }
  fac <- list(current = factor_set(smat(normdata)))
  if (!same_rule) fac$default <- factor_set(smat(filt[!in_def, ]))

  # --- distributions before vs after (log2 quantiles per sample) ---
  qtab <- function(x) { lx <- log2(x); t(apply(lx, 2, function(v) as.numeric(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)))) }
  dist <- list(before = qtab(smat(filt)))
  method <- if (applied) summ$norm_type else NULL
  if (applied && table_exists(summ$table, db_path)) dist$after <- qtab(smat(read_table(summ$table, db_path)))

  # --- status + decision ---
  status <- if (applied) sprintf("Applied: %s", method) else "Waiting for the analyst's decision"
  status_class <- if (applied) "done" else "waiting"
  dec <- read_decisions(run_dir, "06_normalize"); last_dec <- if (length(dec)) dec[[length(dec)]] else NULL
  methods_html <- '<ul class="options">
<li><b>sltmm</b> - SL now, then TMM after imputation (the lab default).</li>
<li><b>sl</b> - sample loading: scale every sample to the same total intensity of the normalization subset.</li>
<li><b>tmm</b> - trimmed mean of M-values (edgeR): robust to a minority of proteins that change.</li>
<li><b>quantile</b>, <b>loess</b> (cyclic loess), <b>vsn</b>, <b>lr</b> (linear regression): reshape each sample&apos;s distribution, not just scale it.</li>
<li><b>directlfq</b> - needs the python setup (not yet wired in the pipeline).</li>
<li>Not available: <b>mti</b> (no function in the app), <b>ai</b> and <b>mi</b> (the app calls them with the wrong arguments).</li></ul>'
  rule_html <- sprintf('<p>Current subset: exclude by <b>%s</b> <code>%s</code>%s. It leaves out <b>%s proteins%s / %s precursors</b> of %s.</p>',
    esc(excl_mode), esc(excl_grep), if (!identical(incl_mode, "not_used")) sprintf(', include only %s <code>%s</code>', esc(incl_mode), esc(incl_grep)) else "",
    fmt_n(n_cur_prot), if (ptm) " (by the exclude rule)" else "", fmt_n(n_cur_prec), fmt_n(nrow(filt)))
  if (ptm) rule_html <- paste0(rule_html, sprintf(
    '<p>PTM run (<code>%s</code>): the samples are enriched, so the pipeline default calculates the factors from the <b>%s PTM precursors</b> only; the <b>%s</b> others are not supposed to be there. They stay in the data, are normalized with the same factors, and are left out of Stats.</p>',
    esc(p$ptm_grep), fmt_n(sum(is_ptm_row(filt, p))), fmt_n(sum(!is_ptm_row(filt, p)))))
  if (!same_rule) rule_html <- paste0(rule_html, sprintf(
    '<p class="missing">This is not the pipeline default rule (only added material: digestion trypsin, human contaminant keratins, bovine alpha/beta casein%s). The default <code>%s</code> would leave out <b>%s precursors</b>. See the list and the factors for both below.</p>',
    if (ptm) sprintf("; and on PTM runs, every precursor without <code>%s</code>", esc(p$ptm_grep)) else "", esc(NORM_DEFAULT_EXCLUDE), fmt_n(n_def_prec)))
  decision_html <- paste0('<section><h2>', if (applied) "Decision" else "Decisions needed", '</h2>',
    '<h3>1. Proteins left out of the factor calculation</h3>', rule_html,
    '<h3>2. Normalization method</h3>',
    if (applied) sprintf('<div class="tiles"><div class="tile"><div class="v">%s</div><div class="l">method</div></div><div class="tile"><div class="v">%s%% &rarr; %s%%</div><div class="l">CV of sample totals, before &rarr; after</div></div></div>%s',
                         esc(method), esc(summ$total_cv_pct_before), esc(summ$total_cv_pct_after), if (!is.null(summ$note)) paste0('<p class="note">', esc(summ$note), '</p>') else "") else "",
    if (applied) paste0('<details><summary>Methods that were available</summary>', methods_html, '</details>') else methods_html,
    if (applied && !is.null(last_dec)) sprintf('<p class="note">Recorded %s%s. %s</p>', esc(last_dec$time),
       if (!is.null(last_dec$analyst)) paste0(" by ", esc(last_dec$analyst)) else "", esc(last_dec$note %||% "")) else "",
    if (!applied) '<p class="note">Review the excluded list, the factors and the plots below, then give both choices in the Claude session.</p>' else "",
    '</section>')

  factors_html <- paste0('<section><h2>Scaling factors per sample</h2>',
    '<p class="note">What each method would multiply each sample by (1 = unchanged), calculated here from the normalization subset with the app&apos;s formulas. ',
    'SL + TMM is an estimate: the real TMM step runs after imputation. Distribution methods (quantile, loess, vsn, lr) have no single factor. ',
    'Factors far from 1 mean that sample&apos;s loading differs; methods that agree give confidence.</p>',
    if (!same_rule) '<div class="control"><label for="facset">Subset</label><select id="facset"><option value="current">current regex</option><option value="default">pipeline default rule</option></select></div>' else "",
    '<div id="facplot" style="height:380px" role="img" aria-label="Scaling factor per sample for each method"></div>',
    '<details><summary>Table view</summary><div id="factable"></div></details></section>')

  excl_html <- paste0('<section><h2>Proteins left out of the factor calculation</h2>',
    if (nzchar(ex_table)) paste0('<details', if (!same_rule) ' open' else '', '><summary>', fmt_n(nrow(agg)), ' proteins (current rule and pipeline default)</summary>', ex_table, '</details>') else '<p>None.</p>',
    '</section>')

  qc <- p$qc_path
  totals <- list(`filtered` = colSums(smat(filt), na.rm = TRUE), `normalization subset` = colSums(smat(normdata), na.rm = TRUE))
  if (applied && table_exists(summ$table, db_path)) totals[[paste0("after ", method, if (identical(method, "sltmm")) " (SL part)" else "")]] <- colSums(smat(read_table(summ$table, db_path)), na.rm = TRUE)
  graphs <- paste0('<section><h2>Graphs</h2><h3>Total intensity per sample</h3>',
    '<p class="note">Filtered data, the normalization subset (what the factors are calculated from), and after normalization. Hover for each total relative to the mean.</p>',
    '<div id="totplot" style="height:340px" role="img" aria-label="Total intensity per sample before and after normalization"></div>',
    sprintf('<script>window.DPDP_TOTALS = %s;</script>', toJSON(lapply(totals, unname), digits = NA)),
    '<details><summary>The app&apos;s original plots</summary><div class="grid2">',
    embed_png(file.path(qc, "Precursor_Filter_barplot.png"), "Filtered bar plot", "App plot: filtered data."),
    embed_png(file.path(qc, "Precursor_NormData_barplot.png"), "NormData bar plot", "App plot: normalization subset."),
    if (applied) embed_png(file.path(qc, paste0("Precursor_", method, "_barplot.png")), "Normalized bar plot", sprintf("App plot: after %s.", method)) else "",
    '</div></details><h2>log2 intensity per sample', if (applied) ', before vs after' else '', '</h2>',
    '<p class="note">Box = 25th to 75th percentile, line = median, whiskers = 5th to 95th percentile. After normalization the medians should line up.</p>',
    '<div id="boxplot" style="height:420px" role="img" aria-label="log2 intensity distribution per sample"></div></section>')

  js_data <- list(samples = labels, groups = grp_of, factors = fac, dist = lapply(dist, function(q) unname(split(q, row(q)))),
                  method = method %||% "")
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA), ';', NORM_JS, '</script>')
  extra_css <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100; } section h3 { font-size:15px; margin:14px 0 4px; }
.control select { font:inherit; font-size:13px; }</style>'
  body <- paste0(decision_html, factors_html, graphs, excl_html, js)
  write_review_page(run_dir, "06_normalize", "Step 6 - Normalize", status, status_class, body, extra_head = extra_css)
}

NORM_JS <- '
(function () { const T = window.DPDP_TOTALS; if (!T) return;
  const css0 = getComputedStyle(document.querySelector(".viz-root")), C0 = n => css0.getPropertyValue(n).trim(), cols = ["--text-muted", "--series-1", "--series-2"];
  Plotly.newPlot("totplot", Object.keys(T).map((k, i) => { const mean = T[k].reduce((a, b) => a + b, 0) / T[k].length;
    return { type: "bar", name: k, x: D.samples, y: T[k], marker: { color: C0(cols[i % 3]) }, customdata: T[k].map(v => (v / mean).toFixed(3)),
      hovertemplate: "%{x}<br>" + k + ": %{y:.3s}<br>%{customdata} x mean<extra></extra>" }; }),
    { barmode: "group", paper_bgcolor: C0("--surface-1"), plot_bgcolor: C0("--surface-1"), font: { color: C0("--text-secondary"), size: 12 },
      margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: "h", y: 1.12 }, bargap: 0.2, bargroupgap: 0.05,
      xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: { title: "total intensity", gridcolor: C0("--grid"), rangemode: "tozero" } },
    { displaylogo: false, responsive: true });
})();
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const SER = ["--series-1", "--series-2", "--series-3", "--series-4"];
const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
function drawFactors() {
  const sel = document.getElementById("facset"), set = D.factors[sel ? sel.value : "current"];
  const names = Object.keys(set);
  const traces = names.map((m, i) => ({ type: "scatter", mode: "markers", name: m, y: D.samples, x: set[m],
    marker: { size: 10, color: C(SER[i % SER.length]), line: { width: 2, color: C("--surface-1") } },
    hovertemplate: "%{y}<br>" + m + ": %{x:.3f}<extra></extra>" }));
  Plotly.react("facplot", traces, Object.assign(base(), { margin: { l: 80, r: 20, t: 10, b: 50 }, showlegend: true, legend: { orientation: "h", y: 1.1 },
    shapes: [{ type: "line", x0: 1, x1: 1, y0: 0, y1: 1, yref: "paper", line: { color: C("--text-muted"), width: 1, dash: "dot" } }],
    xaxis: ax("scaling factor (1 = unchanged)"), yaxis: { autorange: "reversed", gridcolor: C("--grid") } }), { displaylogo: false, responsive: true });
  document.getElementById("factable").innerHTML = "<table><thead><tr><th>Sample</th>" + names.map(m => "<th class=num>" + m + "</th>").join("") +
    "</tr></thead><tbody>" + D.samples.map((s, j) => "<tr><td>" + s + "</td>" + names.map(m => "<td class=num>" + set[m][j].toFixed(3) + "</td>").join("") + "</tr>").join("") + "</tbody></table>";
}
const fs = document.getElementById("facset"); if (fs) fs.addEventListener("change", drawFactors); drawFactors();
const sets = Object.keys(D.dist);
const label = { before: "before (filtered)", after: D.method === "sltmm" ? "after sltmm (SL part; TMM is applied after imputation)" : "after " + D.method };
const btr = sets.map((k, i) => ({ type: "box", name: label[k], x: D.samples, offsetgroup: k,
  q1: D.dist[k].map(q => q[1]), median: D.dist[k].map(q => q[2]), q3: D.dist[k].map(q => q[3]),
  lowerfence: D.dist[k].map(q => q[0]), upperfence: D.dist[k].map(q => q[4]),
  marker: { color: C(SER[i]) }, line: { width: 1.5 }, fillcolor: C(SER[i]) + "55", hoverinfo: "y+name" }));
Plotly.newPlot("boxplot", btr, Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 60 }, boxmode: "group", showlegend: true,
  legend: { orientation: "h", y: 1.08 }, xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: ax("log2 intensity") }), { displaylogo: false, responsive: true });
'
