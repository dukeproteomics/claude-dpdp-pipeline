# Review page for Step 3 - Parameters. Shows the settings (primary group, output level, PTM, peptide selection), sample groups,
# per-sample start numbers, and the app's start bar/box plots and histogram. Accept with:
#   Rscript pipeline/accept.R --db=<db> --step=03_parameters --decided_by=user

review_parameters <- function(db_path) {
  run_dir <- dirname(db_path)
  summ <- fromJSON(file.path(run_dir, "steps", "03_parameters", "summary.json"), simplifyVector = FALSE)
  a <- review_accept_status(run_dir, "03_parameters", summ$generated)
  s <- summ$settings
  set_df <- data.frame(Setting = c("primary_group", "data_output", "ptm", "peptide_select", "use_isoform", "multi_tmt"),
    Value = c(esc(s$primary_group), esc(s$data_output), if (isTRUE(s$ptm)) paste0("yes (", esc(s$ptm_grep), ", localization ", esc(s$ptm_local), ")") else "no",
              esc(s$peptide_select), if (isTRUE(s$use_isoform)) "yes" else "no", if (isTRUE(s$multi_tmt)) "yes" else "no"),
    Meaning = c("Full = use every protein accession in a group; Primary = only the first.", "Report proteins (rollup in Step 8) or peptides.",
                "Phospho/PTM analysis.", "Which peptides count for a protein: Razor, Unique or Shared.", "Keep protein isoforms apart.", "Several TMT sets."))
  ps <- summ$per_sample
  lab <- read_table("design", db_path)$Label
  grp <- unlist(lapply(summ$sample_groups, function(g) rep(g$Group, g$Count)))
  pdat <- list(s = if (length(lab) == length(ps)) lab else vapply(ps, function(x) x$sample, ""), g = grp,
               rel = vapply(ps, function(x) as.numeric(x$total_vs_median), 0), miss = vapply(ps, function(x) as.numeric(x$missing), 0),
               med = vapply(ps, function(x) as.numeric(x$median_log2), 0))
  qc <- get_params(db_path)$qc_path
  body <- paste0('<section><h2>', if (is.null(a)) "Review needed" else "Review recorded", '</h2>', html_table(set_df), accept_note(a), '</section>',
    '<section><h2>Start data</h2><div class="tiles">',
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">precursors</div></div>', fmt_n(summ$counts_start$precursors)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">peptides</div></div>', fmt_n(summ$counts_start$peptides)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">proteins</div></div>', fmt_n(summ$counts_start$proteins)),
    paste(vapply(summ$sample_groups, function(g) sprintf('<div class="tile"><div class="v">%s</div><div class="l">samples in %s</div></div>', g$Count, esc(g$Group)), ""), collapse = ""),
    '</div><h2>Total intensity per sample, relative to the median sample</h2><p class="note">1 = the median sample. Far from 1 means that sample was loaded heavier or lighter; normalization corrects this.</p>',
    '<div id="rel" style="height:280px" role="img" aria-label="Total intensity relative to the median per sample"></div>',
    '<details><summary>Table view</summary>', html_table(data.frame(Sample = pdat$s, Group = grp, `Total vs median` = sprintf("%.3f", pdat$rel),
      `Median log2` = sprintf("%.2f", pdat$med), Missing = fmt_n(pdat$miss), check.names = FALSE), 3:5), '</details></section>',
    '<section><h2>Graphs</h2><div class="grid2">',
    embed_png(file.path(qc, "Precursor_Start_barplot.png"), "Start bar plot", "App plot: total intensity per sample."),
    embed_png(file.path(qc, "Precursor_Start_boxplot.png"), "Start box plot", "App plot: log2 intensity distribution per sample."),
    embed_png(file.path(qc, "Precursor_Start_Histogram.png"), "Start histogram", "App plot: log2 intensity histogram (a low-end shoulder would suggest noise)."),
    '</div></section>',
    '<script>', plotly_js(), '</script><script>const D = ', toJSON(pdat, digits = NA), ';
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim(), cols = ["--series-1", "--series-2", "--series-3"];
Plotly.newPlot("rel", [...new Set(D.g)].map((g, k) => { const i = D.g.map((x, j) => x === g ? j : -1).filter(j => j >= 0);
  return { type: "bar", name: g, x: i.map(j => D.s[j]), y: i.map(j => D.rel[j]), marker: { color: C(cols[k % 3]) }, hovertemplate: "%{x}: %{y:.3f} x median<extra></extra>" }; }),
  { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, margin: { l: 60, r: 10, t: 10, b: 50 },
    legend: { orientation: "h", y: 1.12 }, shapes: [{ type: "line", x0: 0, x1: 1, xref: "paper", y0: 1, y1: 1, line: { dash: "dot", color: C("--text-muted") } }],
    xaxis: { categoryorder: "array", categoryarray: D.s }, yaxis: { title: "total / median sample", gridcolor: C("--grid"), rangemode: "tozero" } },
  { displaylogo: false, responsive: true });</script>')
  extra <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; }</style>'
  write_review_page(run_dir, "03_parameters", "Step 3 - Parameters", if (is.null(a)) "Waiting for the analyst's review" else "Accepted",
                    if (is.null(a)) "waiting" else "done", body, extra_head = extra)
}
