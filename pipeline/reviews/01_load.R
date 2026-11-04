# Review page for Step 1 - Load. Review step: the analyst confirms the design matches the data (sample IDs, groups,
# detected format) before Parameters. Accept with: Rscript pipeline/accept.R --db=<db> --step=01_load --decided_by=user

review_accept_status <- function(run_dir, step, generated) {
  acc <- Filter(function(r) identical(r$settings$accepted, TRUE), read_decisions(run_dir, step))
  a <- if (length(acc)) acc[[length(acc)]] else NULL
  if (!is.null(a) && !is.null(generated) && a$time < generated) a <- NULL
  a
}
accept_note <- function(a) if (is.null(a)) "" else sprintf('<p class="note">Accepted %s%s. %s</p>', esc(a$time),
  if (!is.null(a$analyst)) paste0(" by ", esc(a$analyst)) else "", esc(a$note %||% ""))

review_load <- function(db_path) {
  run_dir <- dirname(db_path)
  summ <- fromJSON(file.path(run_dir, "steps", "01_load", "summary.json"), simplifyVector = FALSE)
  a <- review_accept_status(run_dir, "01_load", summ$generated)
  det <- summ$detected; cn <- summ$counts
  ids <- summ$sample_id_check
  id_df <- data.frame(Order = vapply(ids, function(x) as.character(x$Raw_Order), ""), Label = vapply(ids, function(x) esc(x$Label), ""),
    Group = vapply(ids, function(x) esc(x$Group), ""), ID = vapply(ids, function(x) esc(x$ID), ""),
    `Data column` = vapply(ids, function(x) esc(x$data_column), ""),
    Match = vapply(ids, function(x) if (isTRUE(x$app_check) && isTRUE(x$full_id_check)) "yes" else "<b class=\"bad\">NO</b>", ""), check.names = FALSE)
  all_ok <- all(vapply(ids, function(x) isTRUE(x$app_check) && isTRUE(x$full_id_check), TRUE))
  org <- summ$organisms_top
  org_df <- data.frame(Organism = esc(names(org)), Proteins = fmt_n(unlist(org)), check.names = FALSE)
  miss <- unlist(summ$missing_per_sample); missf <- unlist(summ$missing_fraction_per_sample)
  grp <- vapply(ids, function(x) x$Group, "")
  body <- paste0('<section><h2>', if (is.null(a)) "Review needed: does the design match the data?" else "Review recorded", '</h2>',
    '<p>Check that every sample in the sample sheet found its data column, the groups are right, and the detected format is what you expect.</p>',
    if (all_ok) '<p><b>All sample IDs matched their data columns.</b></p>' else '<p class="missing"><b>Some sample IDs did not match a data column.</b> Fix the sample sheet before going on.</p>',
    if (length(unlist(summ$warnings))) paste0('<p class="missing"><b>Warnings:</b> ', paste(esc(unlist(summ$warnings)), collapse = '; '), '</p>') else '',
    accept_note(a), '</section>',
    '<section><h2>Detected format</h2><div class="tiles">',
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">source (SP = Spectronaut)</div></div>', esc(det$data_source)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">data format</div></div>', esc(det$raw_data_format)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">output</div></div>', esc(det$data_output)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">PTM data</div></div>', if (isTRUE(det$ptm)) "yes" else "no"),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">precursors</div></div>', fmt_n(cn$precursors)),
    sprintf('<div class="tile"><div class="v">%s</div><div class="l">proteins</div></div>', fmt_n(cn$proteins)),
    sprintf('<div class="tile"><div class="v">%s / %s</div><div class="l">samples in data / in sheet</div></div>', cn$sample_columns, cn$design_samples),
    '</div></section>',
    '<section><h2>Samples</h2>', html_table(id_df), '</section>',
    '<section><h2>Missing values per sample (raw data)</h2><div id="miss" style="height:300px" role="img" aria-label="Missing values per sample"></div>',
    '<details><summary>Table view</summary>', html_table(data.frame(Sample = names(miss), Group = grp, Missing = fmt_n(miss), Percent = fmt_pct(missf * 100)), 3:4), '</details></section>',
    '<section><h2>Organisms</h2><p class="note">Proteins per organism in the search results. Spike-ins and contaminants show up here (e.g. yeast ADH, bovine casein, human keratins, pig trypsin).</p>',
    html_table(org_df, 2), '</section>',
    '<script>', plotly_js(), '</script><script>const D = ', toJSON(list(s = names(miss), g = grp, m = as.numeric(missf * 100), n = as.numeric(miss)), digits = NA), ';
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim(), cols = ["--series-1", "--series-2", "--series-3"];
Plotly.newPlot("miss", [...new Set(D.g)].map((g, k) => { const i = D.g.map((x, j) => x === g ? j : -1).filter(j => j >= 0);
  return { type: "bar", name: g, x: i.map(j => D.s[j]), y: i.map(j => D.m[j]), customdata: i.map(j => D.n[j]), marker: { color: C(cols[k % 3]) },
    hovertemplate: "%{x}: %{y:.2f}% missing (%{customdata:,})<extra></extra>" }; }),
  { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, margin: { l: 60, r: 10, t: 10, b: 50 },
    legend: { orientation: "h", y: 1.12 }, xaxis: { categoryorder: "array", categoryarray: D.s }, yaxis: { title: "% missing", gridcolor: C("--grid"), rangemode: "tozero" } },
  { displaylogo: false, responsive: true });</script>')
  extra <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; } .bad { color:var(--warn-ink); }</style>'
  write_review_page(run_dir, "01_load", "Step 1 - Load", if (is.null(a)) "Waiting for the analyst's review" else "Accepted",
                    if (is.null(a)) "waiting" else "done", body, extra_head = extra)
}
