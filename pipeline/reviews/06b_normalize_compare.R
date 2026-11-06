# Review page for Step 6b - normalization comparison. Informational: the Normalize step records the decision.

review_normalize_compare <- function(db_path) {
  run_dir <- dirname(db_path)
  f <- file.path(run_dir, "steps", "06_normalize", "compare.json")
  if (!file.exists(f)) stop("Run 06b_normalize_compare.R first.")
  C <- fromJSON(f, simplifyVector = FALSE)
  M <- C$methods; ids <- names(M)
  cols <- list(
    list(k = "cv_within", lab = "Median CV % within groups", low = TRUE),
    list(k = "cv_pool", lab = paste0("Median CV % in ", C$pool %||% "pool"), low = TRUE),
    list(k = "adh_cv", lab = "ADH spike CV % (information)", low = NA),
    list(k = "int_trend", lab = "Sample spread by intensity (log2 SD)", low = TRUE),
    list(k = "rt_trend", lab = "Sample spread by RT (log2 SD)", low = TRUE),
    list(k = "fc_spread", lab = "Fold-change spread vs SL (1 = same)", low = NA),
    list(k = "truth_error", lab = "Known-truth error (mean |log2|)", low = TRUE))
  if (!isTRUE(C$has_rt)) cols <- Filter(function(c) c$k != "rt_trend", cols)
  if (!isTRUE(C$truth)) cols <- Filter(function(c) c$k != "truth_error", cols)
  if (isTRUE(C$adh_circular)) cols <- lapply(cols, function(c) { if (c$k == "adh_cv") { c$low <- NA; c$lab <- "ADH CV % (circular: subset is ADH)" }; c })
  if (all(vapply(ids, function(m) is.null(M[[m]]$adh_cv), TRUE))) cols <- Filter(function(c) c$k != "adh_cv", cols)
  if (is.null(C$pool) || is.na(C$pool)) cols <- Filter(function(c) c$k != "cv_pool", cols)
  num <- function(v) if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
  best <- lapply(cols, function(c) { v <- sapply(ids, function(m) num(M[[m]][[c$k]]))
    if (isTRUE(c$low) && any(!is.na(v))) ids[which(v <= min(v, na.rm = TRUE) + 1e-9)] else character(0) })
  warn <- function(m, c) c$k == "fc_spread" && !is.na(num(M[[m]]$fc_spread)) && abs(num(M[[m]]$fc_spread) - 1) > 0.1
  head <- paste0("<tr><th>Method</th><th>Corrects</th>", paste0("<th class=num>", esc(sapply(cols, `[[`, "lab")), "</th>", collapse = ""), "</tr>")
  fam <- c(scaling = "one factor per sample", intensity = "intensity bias", rt = "drift along the gradient", distribution = "whole distribution", none = "nothing")
  rows <- vapply(ids, function(m) paste0("<tr><td><b>", esc(m), "</b> <span class=note>", esc(M[[m]]$label), "</span></td><td>", esc(fam[[M[[m]]$family]]), "</td>",
    paste0(vapply(seq_along(cols), function(j) { v <- num(M[[m]][[cols[[j]]$k]]); s <- if (is.na(v)) "-" else format(v, big.mark = ",")
      sprintf('<td class="num%s%s">%s</td>', if (m %in% best[[j]]) " best" else "", if (warn(m, cols[[j]])) " warnc" else "", s) }, ""), collapse = ""), "</tr>"), "")
  tbl <- sprintf("<div style='overflow-x:auto'><table class='cmp'><thead>%s</thead><tbody>%s</tbody></table></div>", head, paste(rows, collapse = ""))
  fcm <- if (length(C$pairs)) paste0('<details><summary>Median log2 fold change of each comparison (what the method treats as unchanged)</summary>',
    html_table(cbind(data.frame(Method = ids), setNames(as.data.frame(do.call(rbind, lapply(ids, function(m) unlist(M[[m]]$fc_median)))), unlist(C$pairs))),
               seq_len(length(C$pairs)) + 1), '</details>') else ""
  body <- paste0('<section><h2>What this compares</h2>',
    sprintf('<p>Every method is fitted on the same normalization subset (%s of %s precursors; exclude <code>%s</code>; include %s) and applied to all values. ',
      fmt_n(C$subset$used_for_factors), fmt_n(C$subset$precursors), esc(C$subset$rule_exclude %||% ""), esc(C$subset$rule_include %||% "")),
    'Scores use measured values only. The subset decides what is assumed constant; change it on the Step 6 page, then re-run this comparison. ',
    'There is no single score: read the columns together.</p><ul class="options">',
    '<li><b>CVs</b>: lower means replicates agree better. Methods that shrink differences (fold-change spread below 1) also lower CVs, so a low CV alone is not proof.</li>',
    '<li><b>ADH spike CV</b> (information, not a ranking): ADH is added at the same amount to every sample, so when protein loads differ, a correct loading normalization makes ADH vary; its CV is usually lowest with <b>none</b>. A big rise after a method means that method moved the samples a lot.</li>',
    '<li><b>Sample spread by intensity / RT</b>: how much the samples still differ from each other within each intensity or retention-time decile (SD of log ratios, averaged over deciles). Scaling methods cannot remove these; curve methods can.</li>',
    '<li><b>Fold-change spread vs SL</b>: scaling methods are exactly 1. Values far from 1 (marked) mean the method changes the size of the differences between groups.</li>',
    if (isTRUE(C$adh_circular)) '<li class="missing"><b>ADH CV is not a test here:</b> the normalization subset is the ADH spike itself, so the factors were made to equalize it.</li>' else "",
    if (isTRUE(C$truth)) '<li><b>Known truth</b>: group names give the species mix, so the expected log2 ratios are known (human, E. coli, yeast).</li>' else "",
    '</ul></section>',
    '<section><h2>Scores</h2>', tbl, fcm, '</section>',
    '<section><h2>Charts</h2><div class="control"><label for="msel">Method</label> <select id="msel">',
    paste0(sprintf('<option value="%s"%s>%s</option>', ids, ifelse(ids == "sl", " selected", ""), esc(ids)), collapse = ""), '</select></div>',
    '<h3>Each sample&apos;s log ratio to the reference profile by intensity, centred across samples (all lines on 0 = no bias between samples)</h3><div id="c_int" style="height:360px"></div>',
    if (isTRUE(C$has_rt)) '<h3>The same by retention time</h3><div id="c_rt" style="height:360px"></div>' else "",
    '<div class="grid2"><div><h3>Median CV % within groups</h3><div id="c_cv" style="height:320px"></div></div>',
    '<div><h3>Fold-change spread vs SL</h3><div id="c_fc" style="height:320px"></div></div></div></section>',
    sprintf('<p class="note">Run took %s s. Data: %s.</p>', C$seconds, esc(C$db_path)),
    '<script>', plotly_js(), '</script><script>const CMP = ', toJSON(C, auto_unbox = TRUE, digits = NA, null = "null", na = "null"), ';', NCMP_JS, '</script>')
  extra <- '<style>.cmp td.best { font-weight:700; background:var(--ok-bg); color:var(--ok-ink); } .cmp td.warnc { background:var(--warn-bg, #fff3cd); } .cmp th { font-size:12px; }</style>'
  write_review_page(run_dir, "06b_normalize_compare", "Step 6b - Normalization comparison", "Information for the Normalize decision", "done", body, extra_head = extra)
}

NCMP_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
const M = CMP.methods, ids = Object.keys(M), S = CMP.samples, G = CMP.groups, ug = [...new Set(G)];
const gcol = g => C("--series-" + ((ug.indexOf(g) % 6) + 1)) || "#888";
const lay = (xt, yt) => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 11 },
  margin: { l: 60, r: 20, t: 10, b: 45 }, legend: { orientation: "h", y: -0.2 }, xaxis: { title: xt, gridcolor: C("--grid"), zeroline: false },
  yaxis: { title: yt, gridcolor: C("--grid"), zeroline: true, zerolinecolor: C("--text-muted") } });
function curves(id, mid, cur) { if (!document.getElementById(id) || !cur) return;
  Plotly.react(id, S.map((s, j) => ({ type: "scatter", mode: "lines+markers", name: s, legendgroup: G[j], x: mid, y: cur[j],
    line: { width: 1.5, color: gcol(G[j]) }, marker: { size: 4 }, hovertemplate: s + " (" + G[j] + ")<br>%{y:.3f}<extra></extra>" })),
    lay(id === "c_rt" ? "retention time (decile medians)" : "log2 intensity (decile medians)", "log2 ratio to reference"), { displaylogo: false, responsive: true }); }
function draw() { const m = document.getElementById("msel").value; curves("c_int", M[m].int_mid, M[m].int_curves); curves("c_rt", M[m].rt_mid, M[m].rt_curves); }
document.getElementById("msel").addEventListener("change", draw); draw();
const hb = (id, vals, xt, ref) => { const L = lay(xt, ""); L.margin.l = 110; L.yaxis = { automargin: true, autorange: "reversed" }; L.showlegend = false;
  if (ref !== undefined) L.shapes = [{ type: "line", x0: ref, x1: ref, y0: 0, y1: 1, yref: "paper", line: { dash: "dot", color: C("--text-muted") } }];
  Plotly.newPlot(id, [{ type: "bar", orientation: "h", y: ids, x: vals, marker: { color: C("--series-1") }, hovertemplate: "%{y}: %{x}<extra></extra>" }], L, { displayModeBar: false, responsive: true }); };
hb("c_cv", ids.map(m => M[m].cv_within), "median CV %");
hb("c_fc", ids.map(m => M[m].fc_spread), "IQR of log2 FC relative to SL", 1);
'
