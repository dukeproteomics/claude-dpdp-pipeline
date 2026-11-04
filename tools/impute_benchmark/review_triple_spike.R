# Review page for the triple-spike imputation benchmark (informational; no decision).
review_triple_spike <- function(out_dir) {
  B <- fromJSON(file.path(out_dir, "benchmark.json"), simplifyVector = FALSE)
  M <- B$methods; ids <- names(M)
  cols <- list(
    list(k = "mask_rmse", lab = "Masking error (RMSE, log2)", low = TRUE),
    list(k = "proteins", lab = "Proteins", low = NA),
    list(k = "tp_total", lab = "True calls (6 comparisons)", low = FALSE),
    list(k = "fp_total", lab = "False calls", low = TRUE),
    list(k = "fdr_mean", lab = "Observed FDR %, mean (target 5)", low = TRUE),
    list(k = "fdr_max", lab = "Observed FDR %, worst comparison", low = TRUE),
    list(k = "sensitivity_mean", lab = "Sensitivity % (spike proteins with |FC| >= 1.5)", low = FALSE),
    list(k = "human_p01_mean", lab = "Human (unchanged): % p < 0.01 (expect 1)", low = TRUE),
    list(k = "abs_error_mean", lab = "Spike log2 FC error, |median| vs expected", low = TRUE),
    list(k = "iqr_mean", lab = "Spike log2 FC spread (IQR)", low = TRUE),
    list(k = "human_iqr_mean", lab = "Human log2 FC spread (IQR)", low = TRUE))
  val <- function(m, c) { v <- if (c$k == "mask_rmse") M[[m]]$masked$rmse else M[[m]][[c$k]]; if (is.null(v)) NA else v }
  best <- lapply(cols, function(c) { v <- sapply(ids, function(m) suppressWarnings(as.numeric(val(m, c))))
    if (is.na(c$low) || all(is.na(v))) character(0) else if (c$low) ids[which(v == min(v, na.rm = TRUE))] else ids[which(v == max(v, na.rm = TRUE))] })
  head <- paste0("<tr><th>Method</th>", paste0("<th class=num>", esc(sapply(cols, `[[`, "lab")), "</th>", collapse = ""), "</tr>")
  rows <- vapply(ids, function(m) paste0("<tr><td>", esc(M[[m]]$label), "</td>", paste0(vapply(seq_along(cols), function(j) {
    v <- val(m, cols[[j]]); s <- if (is.na(v)) "-" else format(v, big.mark = ",")
    sprintf('<td class="num%s">%s</td>', if (m %in% best[[j]]) " best" else "", s) }, ""), collapse = ""), "</tr>"), "")
  tbl <- sprintf("<div style='overflow-x:auto'><table class='cmp'><thead>%s</thead><tbody>%s</tbody></table></div>", head, paste(rows, collapse = ""))
  comp_rows <- do.call(rbind, lapply(B$comparisons, function(cp) data.frame(Comparison = cp$name,
    `E. coli` = cp$expected$ecoli, Human = cp$expected$human, Yeast = cp$expected$yeast, check.names = FALSE)))
  gaps <- unlist(B$gaps_per_sample); med_g <- median(gaps)
  odd <- names(gaps)[gaps > 2 * med_g]
  opts <- paste0(vapply(B$comparisons, function(cp) sprintf('<option value="%s">%s</option>', esc(cp$name), esc(cp$name)), ""), collapse = "")
  body <- paste0(
    '<section><h2>What this tests</h2><p>Three proteomes mixed at known levels: human constant (log2 FC 0), E. coli and yeast at different levels per group. ',
    'Every method starts from the same ', fmt_n(B$start$precursors), ' precursors (', fmt_n(B$start$gaps), ' gaps; ',
    paste(names(B$start$species), fmt_n(unlist(B$start$species)), sep = " ", collapse = ", "), ' proteins), ',
    'normalized so the human precursors measured in every sample have equal totals. Then MaxLFQ (limpa: DPC-Quant) and limma (limpa: dpcDE). ',
    'Calls are adj p &le; 0.05 and |FC| &ge; 1.5. A true call is an E. coli or yeast protein with an expected change of at least 1.5-fold, in the right direction; ',
    'a false call is a human protein, or a spike protein in the wrong direction.</p>',
    html_table(comp_rows, 2:4), '<p class="note">Expected log2 FC per comparison (numerator vs denominator).</p>',
    if (length(odd)) sprintf('<p class="missing">Sample %s has %s gaps, more than twice the median sample (%s). It is kept, so comparisons with its group are harder for every method.</p>',
                            esc(paste(odd, collapse = ", ")), fmt_n(max(gaps[odd])), fmt_n(med_g)) else "",
    '</section>',
    '<section><h2>Scores (all 6 comparisons)</h2>', tbl, '<p class="note">The best value in each scored column is marked. Masking error is not possible for background, no imputation or limpa.</p></section>',
    '<section><h2>By comparison</h2><div class="control"><label>Comparison <select id="cmpSel">', opts, '</select></label></div>',
    '<div class="grid2"><div><h3>True calls (higher is better)</h3><div id="c_tp" style="height:340px"></div></div>',
    '<div><h3>Observed FDR % (target 5)</h3><div id="c_fdr" style="height:340px"></div></div></div>',
    '<h3>log2 FC per species (box: 25-75%, whiskers: 5-95%); dotted lines are the expected values</h3><div id="c_fc" style="height:560px"></div>',
    '<div id="t_cmp"></div></section>',
    sprintf('<p class="note">Run took %s s. Raw: %s. Background: %s.</p>', B$seconds, esc(B$raw), esc(B$background %||% "none")),
    '<script>', plotly_js(), '</script><script>const M = ', toJSON(M, auto_unbox = TRUE, digits = NA, null = "null", na = "null"),
    ', COMPS = ', toJSON(B$comparisons, auto_unbox = TRUE, digits = NA), ';', TS_JS, '</script>')
  extra <- '<style>.cmp td.best { font-weight:700; background:var(--ok-bg); color:var(--ok-ink); } .cmp th { font-size:12px; } .cmp td.num { white-space:nowrap; }
.viz-root { --sp-human:#2a78d6; --sp-ecoli:#1baf7a; --sp-yeast:#eb6834; }
@media (prefers-color-scheme: dark) { .viz-root { --sp-human:#3987e5; --sp-ecoli:#199e70; --sp-yeast:#d95926; } }</style>'
  write_review_page(out_dir, "triple_spike_benchmark", "Imputation benchmark - three-species spike-in", "Information", "done", body, extra_head = extra)
}

TS_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
const ids = Object.keys(M), lab = ids.map(m => M[m].label);
const SP = [["human", "Human"], ["ecoli", "E. coli"], ["yeast", "Yeast"]];
const base = (xt) => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 11 },
  margin: { l: 300, r: 20, t: 10, b: 40 }, showlegend: false, yaxis: { automargin: true, autorange: "reversed" }, xaxis: { title: xt, gridcolor: C("--grid"), zeroline: false } });
const cfg = { displayModeBar: false, responsive: true };
function draw() {
  const c = document.getElementById("cmpSel").value, cp = COMPS.find(x => x.name === c);
  const per = m => M[m].per[c];
  Plotly.react("c_tp", [{ type: "bar", orientation: "h", y: lab, x: ids.map(m => per(m).tp), marker: { color: C("--series-1") },
    customdata: ids.map(m => [per(m).possible, per(m).sensitivity]), hovertemplate: "%{y}<br>%{x} true calls of %{customdata[0]} possible (%{customdata[1]}%)<extra></extra>" }], base("true calls"), cfg);
  const L2 = base("observed FDR %"); L2.shapes = [{ type: "line", x0: 5, x1: 5, y0: 0, y1: 1, yref: "paper", line: { dash: "dot", color: C("--text-muted") } }];
  L2.annotations = [{ x: 5, y: 1, yref: "paper", text: "5% target", showarrow: false, xanchor: "left", yanchor: "bottom", font: { color: C("--text-muted") } }];
  Plotly.react("c_fdr", [{ type: "bar", orientation: "h", y: lab, x: ids.map(m => per(m).fdr), marker: { color: C("--series-1") },
    customdata: ids.map(m => [per(m).fp_human, per(m).fp_direction]), hovertemplate: "%{y}<br>FDR %{x}%: %{customdata[0]} human, %{customdata[1]} wrong direction<extra></extra>" }], L2, cfg);
  const tr = SP.map(([k, nm]) => { const q = ids.map(m => per(m).species[k].q || [null, null, null, null, null]);
    return { type: "box", orientation: "h", name: nm, y: lab, lowerfence: q.map(v => v[0]), q1: q.map(v => v[1]), median: q.map(v => v[2]), q3: q.map(v => v[3]), upperfence: q.map(v => v[4]),
      marker: { color: C("--sp-" + k) }, line: { width: 1.5 }, hoverinfo: "x+name" }; });
  const L3 = base("log2 fold change"); L3.boxmode = "group"; L3.showlegend = true; L3.legend = { orientation: "h", y: 1.06 };
  L3.shapes = SP.map(([k]) => ({ type: "line", x0: cp.expected[k], x1: cp.expected[k], y0: 0, y1: 1, yref: "paper", line: { dash: "dot", width: 2, color: C("--sp-" + k) } }));
  Plotly.react("c_fc", tr, L3, cfg);
  let h = "<table class=cmp><thead><tr><th>Method</th><th class=num>true calls</th><th class=num>possible</th><th class=num>human called</th><th class=num>wrong direction</th><th class=num>FDR %</th>" +
    SP.map(([k, nm]) => "<th class=num>" + nm + " n</th><th class=num>" + nm + " median error</th><th class=num>" + nm + " IQR</th>").join("") + "</tr></thead><tbody>";
  ids.forEach(m => { const r = per(m); h += "<tr><td>" + M[m].label + "</td><td class=num>" + r.tp + "</td><td class=num>" + r.possible + "</td><td class=num>" + r.fp_human +
    "</td><td class=num>" + r.fp_direction + "</td><td class=num>" + r.fdr + "</td>" + SP.map(([k]) => { const s = r.species[k];
      return "<td class=num>" + s.n + "</td><td class=num>" + s.median_error + "</td><td class=num>" + s.iqr + "</td>"; }).join("") + "</tr>"; });
  document.getElementById("t_cmp").innerHTML = "<h3>Table view</h3><div style=\\"overflow-x:auto\\">" + h + "</tbody></table></div>";
}
document.getElementById("cmpSel").addEventListener("change", draw); draw();
'
