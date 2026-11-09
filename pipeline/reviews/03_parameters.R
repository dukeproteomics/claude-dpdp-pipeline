# Review page for Step 3 - Parameters. The Options panel (pipeline/options/03_parameters.R: every option that applies to this data
# source; Run with these settings / Accept), start counts, per-sample start numbers, and the app's start bar/box plots and histogram.
# Accept with: Rscript pipeline/accept.R --db=<db> --step=03_parameters --decided_by=user [--from_page]

review_parameters <- function(db_path) {
  run_dir <- dirname(db_path)
  summ <- fromJSON(file.path(run_dir, "steps", "03_parameters", "summary.json"), simplifyVector = FALSE)
  a <- review_accept_status(run_dir, "03_parameters", summ$generated)
  s <- summ$settings
  ls <- file.path(run_dir, "steps", "01_load", "summary.json")
  opts <- parameters_options(get_params(db_path), read_table("design", db_path),
                             if (file.exists(ls)) fromJSON(ls, simplifyVector = FALSE)$detected else NULL)
  form <- options_form("03_parameters", opts, s, previewed = s, applied = !is.null(a), run_dir = run_dir, title = "Parameters",
    words = list(preview = "Run", apply = "Accept", previewing = "the last run", applied = "accepted",
                 effect = "(re-runs this step, about a minute, and rebuilds this page)"))
  ps <- summ$per_sample
  lab <- read_table("design", db_path)$Label
  grp <- unlist(lapply(summ$sample_groups, function(g) rep(g$Group, g$Count)))
  pdat <- list(s = if (length(lab) == length(ps)) lab else vapply(ps, function(x) x$sample, ""), g = grp,
               rel = vapply(ps, function(x) as.numeric(x$total_vs_median), 0), miss = vapply(ps, function(x) as.numeric(x$missing), 0),
               med = vapply(ps, function(x) as.numeric(x$median_log2), 0))
  pp <- get_params(db_path); qc <- pp$qc_path
  st <- read_table("precursor_start", db_path); sm <- as.matrix(st[, (ncol(st) - length(pdat$s) + 1):ncol(st)])
  lmean <- as.numeric(pp$intensity_mean); lsd <- as.numeric(pp$intensity_sd)
  b2 <- as.numeric(quantile(log2(sm[!is.na(sm) & sm > 0]), 0.02))
  ref <- list()
  if (is.finite(lmean) && is.finite(lsd)) ref <- list(`mean` = list(x = lmean, color = "--text-primary", dash = "solid"),
    `-1 SD` = list(x = lmean - lsd), `+1 SD` = list(x = lmean + lsd),
    `misaligned cutoff` = list(x = log2(as.numeric(pp$intensity_cutoff)), color = "--series-3", dash = "dash"))
  ref$`bottom 2%` <- list(x = b2, color = "--series-2", dash = "dash")
  start_charts <- sample_charts(sm, pdat$s, grp, "start", ref)
  body <- paste0(if (!is.null(a)) paste0('<section><h2>Review recorded</h2>', accept_note(a), '</section>') else "", form,
    '<section><h2>Start data</h2><div class="tiles">',
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">precursors</div></div>', fmt_n(summ$counts_start$precursors)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">peptides</div></div>', fmt_n(summ$counts_start$peptides)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">proteins</div></div>', fmt_n(summ$counts_start$proteins)),
    paste(vapply(summ$sample_groups, function(g) sprintf('<div class="tile"><div class="v">%s</div><div class="l">samples in %s</div></div>', g$Count, esc(g$Group)), ""), collapse = ""),
    '</div><h2>Total intensity per sample, relative to the median sample</h2><p class="note">1 = the median sample. Far from 1 means that sample was loaded heavier or lighter; normalization corrects this.</p>',
    '<div id="rel" style="height:280px" role="img" aria-label="Total intensity relative to the median per sample"></div>',
    '<details><summary>Table view</summary>', html_table(data.frame(Sample = pdat$s, Group = grp, `Total vs median` = sprintf("%.3f", pdat$rel),
      `Median log2` = sprintf("%.2f", pdat$med), Missing = fmt_n(pdat$miss), check.names = FALSE), 3:5), '</details></section>',
    '<section><h2>Graphs</h2>', start_charts,
    '<p class="note">Histogram lines as in the app: mean, &plusmn; 1 SD, the misaligned intensity cutoff (mean + ', esc(pp$intensity_cutoff_sd %||% "0.5"),
    ' SD, used in Filter) and the bottom 2% (the bottom-X fill range in Impute). A noise population would show as a shoulder or second peak at the low end.</p>',
    '<details><summary>The app&apos;s original plots</summary><div class="grid2">',
    embed_png(file.path(qc, "Precursor_Start_barplot.png"), "Start bar plot", "App plot: total intensity per sample."),
    embed_png(file.path(qc, "Precursor_Start_boxplot.png"), "Start box plot", "App plot: log2 intensity distribution per sample."),
    embed_png(file.path(qc, "Precursor_Start_Histogram.png"), "Start histogram", "App plot: log2 intensity histogram."),
    '</div></details></section>',
    '<script>', plotly_js(), '</script><script>const D = ', toJSON(pdat, digits = NA), ';
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim(), cols = ["--series-1", "--series-2", "--series-3"];
Plotly.newPlot("rel", [...new Set(D.g)].map((g, k) => { const i = D.g.map((x, j) => x === g ? j : -1).filter(j => j >= 0);
  return { type: "bar", name: g, x: i.map(j => D.s[j]), y: i.map(j => D.rel[j]), marker: { color: DPDP_GCOL(g, [...new Set(D.g)]) }, hovertemplate: "%{x}: %{y:.3f} x median<extra></extra>" }; }),
  { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, margin: { l: 60, r: 10, t: 10, b: 50 },
    legend: { orientation: "h", y: 1.12 }, shapes: [{ type: "line", x0: 0, x1: 1, xref: "paper", y0: 1, y1: 1, line: { dash: "dot", color: C("--text-muted") } }],
    xaxis: { categoryorder: "array", categoryarray: D.s }, yaxis: { title: "total / median sample", gridcolor: C("--grid"), rangemode: "tozero" } },
  { displaylogo: false, responsive: true });</script>')
  extra <- '<style>.viz-root { }</style>'
  write_review_page(run_dir, "03_parameters", "Step 3 - Parameters", if (is.null(a)) "Waiting for the analyst's review" else "Accepted",
                    if (is.null(a)) "waiting" else "done", body, extra_head = extra)
}
