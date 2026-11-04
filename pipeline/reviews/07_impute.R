# Review page for Step 7 - Impute.
# Decisions: impute_type (duke / bottomx), missing_cutoff (which rows get duke), bottom_x (bottom-X fill range).
# Rules (app: impute_duke, impute_bottomx, impute_meta_data_bg):
#   duke     per group, rows with 0 < missing <= Count * (100 - missing_cutoff)/100 get group mean (log2) + U(-1,1) * SD,
#            SD from 20 intensity bins of the group (built from precursor_filter). All other missing values -> bottom-X.
#   bottom-X uniform between the lowest log2 value and the X-th percentile of all values of the table being imputed.
# Shows: fill plan explorer (missing_cutoff), bottom-X range explorer on the intensity histogram, duke SD by intensity,
# the app's missing value plots; after apply: measured vs duke vs bottom-X histogram, imputed % per sample, the TMM factors
# applied here for sltmm, and the imputed bar plots.

review_impute <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  design <- read_table("design", db_path)
  groups <- read_table("sample_groups", db_path)
  labels <- if (!is.null(design$Label) && length(design$Label) == n_s) design$Label else NULL
  grp_of <- rep(groups$Group, as.integer(groups$Count))
  norms <- trimws(strsplit(p$norm_type, ",")[[1]])
  sd_dir <- file.path(run_dir, "steps", "07_impute")
  summ <- if (file.exists(file.path(sd_dir, "summary.json"))) fromJSON(file.path(sd_dir, "summary.json"), simplifyVector = FALSE) else NULL
  applied <- !is.null(summ) && identical(summ$action, "apply")
  s <- summ$settings %||% list()
  itype <- s$impute_type %||% p$impute_type; bx <- as.numeric(s$bottom_x %||% p$bottom_x); mc <- as.numeric(s$missing_cutoff %||% p$missing_cutoff)
  smat <- function(d) { x <- as.matrix(d[, (ncol(d) - n_s + 1):ncol(d)]); x[!is.na(x) & x <= 0] <- NA; x }
  filt <- smat(read_table("precursor_filter", db_path))
  if (is.null(labels)) labels <- colnames(filt)
  na_mask <- is.na(filt)
  gcols <- lapply(seq_len(nrow(groups)), function(i) as.integer(groups$start[i]):as.integer(groups$end[i]))

  # rows by number missing, per group (drives the duke / bottom-X split for any missing_cutoff)
  rows_by_k <- lapply(gcols, function(cc) { k <- rowSums(na_mask[, cc, drop = FALSE]); as.integer(tabulate(k[k > 0], nbins = length(cc))) })

  # bottom-X range for X = 1..25, from the first normalized table (the app imputes each norm table separately)
  tab1 <- paste0("precursor_norm_", norms[1])
  lv <- sort(log2(smat(read_table(tab1, db_path)))); lv <- lv[lv > 0]; nn <- length(lv)
  bx_grid <- 1:25
  bx_max <- vapply(bx_grid, function(x) lv[ceiling(x * nn / 100)], 0); bx_min <- lv[1]
  step <- 0.1; edges <- seq(floor(min(lv) / step) * step, ceiling(max(lv) / step) * step, by = step)
  hist_meas <- tabulate(findInterval(lv, edges, rightmost.closed = TRUE), nbins = length(edges) - 1)

  # duke SD by intensity bin, per group
  bins <- lapply(groups$Group, function(g) { b <- read_table(paste0("impute_bin_", g), db_path); list(lo = b$min2, hi = b$max, sd = b$sd, mid = (b$min + pmin(b$max, max(lv))) / 2) })
  names(bins) <- groups$Group

  # --- after apply: measured / duke / bottom-X distributions (impute-only table, before TMM) and per-sample % ---
  after <- NULL
  if (applied && table_exists("precursor_impute_impute", db_path)) {
    imp <- smat(read_table("precursor_impute_impute", db_path))
    duke_mask <- matrix(FALSE, nrow(filt), ncol(filt))
    if (identical(itype, "duke")) for (i in seq_along(gcols)) {
      cc <- gcols[[i]]; k <- rowSums(na_mask[, cc, drop = FALSE]); ok <- k > 0 & k <= length(cc) * (100 - mc) / 100
      duke_mask[ok, cc] <- na_mask[ok, cc]
    }
    bot_mask <- na_mask & !duke_mask
    li <- log2(imp)
    hb <- function(v) tabulate(findInterval(v[!is.na(v) & v > 0], edges, rightmost.closed = TRUE), nbins = length(edges) - 1)
    after <- list(meas = hb(li[!na_mask]), duke = hb(li[duke_mask]), bot = hb(li[bot_mask]),
                  duke_pct = round(colSums(duke_mask) / nrow(filt) * 100, 2), bot_pct = round(colSums(bot_mask) / nrow(filt) * 100, 2),
                  duke_n = colSums(duke_mask), bot_n = colSums(bot_mask))
    if ("sltmm" %in% norms && table_exists("precursor_impute_sltmm", db_path) && table_exists("precursor_norm_sltmm", db_path)) {
      a <- smat(read_table("precursor_impute_sltmm", db_path)); b <- smat(read_table("precursor_norm_sltmm", db_path))
      after$tmm <- round(apply(a / b, 2, median, na.rm = TRUE), 4)
    }
  }

  # --- status + decision ---
  status <- if (applied) sprintf("Applied: %s", itype) else "Waiting for the analyst's decision"
  status_class <- if (applied) "done" else "waiting"
  dec <- read_decisions(run_dir, "07_impute"); last_dec <- if (length(dec)) dec[[length(dec)]] else NULL
  tot_missing <- sum(na_mask)
  options_html <- sprintf('<ul class="options">
<li><b>duke</b> - a precursor missing in only a few samples of a group is filled near that group&apos;s mean (mean &plusmn; random &times; the SD of
precursors at that intensity). Everything else gets bottom-X. <code>missing_cutoff</code> = %s%% decides &quot;only a few&quot;: duke is used when the group has
at least %s%% of its samples measured.</li>
<li><b>bottomx</b> - every missing value is filled with a random value between the lowest intensity and the bottom <b>%s%%</b> of all data
(<code>bottom_x</code>), i.e. treated as below detection.</li></ul>
<p class="note">Other methods in the app&apos;s list are not implemented in <code>impute_apply_bg2()</code>. The random numbers are seeded
(set.seed(123)), so a re-run gives the same values.</p>', esc(mc), esc(mc), esc(bx))
  decision_html <- paste0('<section><h2>', if (applied) "Decision" else "Decision needed", '</h2>',
    if (applied) paste0('<div class="tiles">',
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">impute_type</div></div>', esc(itype)),
      sprintf('<div class="tile"><div class="v">%s%%</div><div class="l">missing_cutoff</div></div>', esc(mc)),
      sprintf('<div class="tile"><div class="v">%s%%</div><div class="l">bottom_x</div></div>', esc(bx)),
      sprintf('<div class="tile"><div class="v">%s</div><div class="l">values filled (%s of all)</div></div>', fmt_n(tot_missing), fmt_pct(tot_missing / length(na_mask) * 100)),
      '</div>', if (!is.null(summ$note)) paste0('<p class="note">', esc(summ$note), '</p>') else "",
      if (!is.null(last_dec)) sprintf('<p class="note">Recorded %s%s. %s</p>', esc(last_dec$time), if (!is.null(last_dec$analyst)) paste0(" by ", esc(last_dec$analyst)) else "", esc(last_dec$note %||% "")) else "",
      '<details><summary>Options</summary>', options_html, '</details>') else options_html,
    if (!applied) '<p class="note">Use the explorers below, then give the choice in the Claude session.</p>' else "",
    '</section>')

  plan_html <- '<section><h2>Fill plan: duke vs bottom-X by group</h2>
<p class="note">How the missing values of each group would be filled at the chosen <code>missing_cutoff</code>. A group with a large bottom-X share has many
precursors missing in most of its samples; those get low values, which can create large fold changes.</p>
<div class="control"><label for="mc">missing_cutoff</label><input type="range" id="mc" min="0" max="100" step="1"/><span class="val" id="mcv"></span>
<label for="itype">impute_type</label><select id="itype"><option>duke</option><option>bottomx</option></select></div>
<div id="planplot" style="height:260px" role="img" aria-label="Missing values per group, filled by duke or bottom-X"></div>
<details><summary>Table view</summary><div id="plantable"></div></details></section>'

  bx_html <- sprintf('<section><h2>Bottom-X fill range</h2>
<p class="note">Bottom-X values are drawn evenly between the lowest measured intensity and the X-th percentile of <code>%s</code> (shaded band).</p>
<div class="control"><label for="bx">bottom_x</label><input type="range" id="bx" min="0" max="%d" step="1"/><span class="val" id="bxv"></span></div>
<div id="bxplot" style="height:340px" role="img" aria-label="Intensity histogram with the bottom-X fill range"></div></section>', esc(tab1), length(bx_grid) - 1)

  duke_html <- '<section><h2>Duke fill spread by intensity</h2>
<p class="note">Duke fills are the group mean &plusmn; up to this SD (log2). The SD comes from the measured precursors in the same intensity bin (20 bins per group).</p>
<div id="dukeplot" style="height:300px" role="img" aria-label="Duke SD by intensity bin per group"></div></section>'

  after_html <- if (!is.null(after)) paste0('<section><h2>Result: where the filled values sit</h2>',
    '<p class="note">log2 intensities of the impute-only table (before TMM): measured values vs duke fills vs bottom-X fills. Bottom-X should sit at the low end;
duke fills should look like the measured distribution.</p><div id="afterhist" style="height:340px" role="img" aria-label="Histogram of measured, duke and bottom-X values"></div>',
    '<h2>Filled values per sample</h2><div id="afterbar" style="height:300px" role="img" aria-label="Percent of values filled per sample"></div>',
    '<details><summary>Table view</summary>', html_table(data.frame(Sample = labels, Group = grp_of, Duke = fmt_n(after$duke_n), `Duke %` = fmt_pct(after$duke_pct),
      `Bottom-X` = fmt_n(after$bot_n), `Bottom-X %` = fmt_pct(after$bot_pct), check.names = FALSE), 3:6), '</details>',
    if (!is.null(after$tmm)) paste0('<h2>TMM factors applied here (sltmm)</h2><p class="note">For sltmm the app applies TMM after imputation. Median ratio of <code>precursor_impute_sltmm</code> to <code>precursor_norm_sltmm</code> per sample.</p>',
      html_table(data.frame(Sample = labels, Group = grp_of, `TMM factor` = formatC(after$tmm, format = "f", digits = 3), check.names = FALSE), 3)) else "",
    '</section>') else ""

  qc <- p$qc_path
  graphs <- paste0('<section><h2>Graphs</h2><div class="grid2">',
    embed_png(file.path(qc, "missing_bar_plot.png"), "Missing values bar plot", "App plot: missing values per sample."),
    embed_png(file.path(qc, "missing_percent_plot.png"), "Missing values percent plot", "App plot: percent missing per sample."),
    if (applied) paste(vapply(norms, function(n) embed_png(file.path(qc, paste0("Precursor_impute_", n, "_barplot.png")), paste("Imputed", n),
      sprintf("Pipeline plot: total intensity per sample after imputation (%s).%s", n, if (n == "sltmm") " Totals need not be equal after TMM." else "")), ""), collapse = "") else "",
    '</div></section>')

  js_data <- list(groups = groups$Group, counts = as.integer(groups$Count), rows_by_k = rows_by_k, mc = mc, itype = itype,
                  edges = edges, hist = hist_meas, bx_grid = bx_grid, bx_max = bx_max, bx_min = bx_min, bx = bx,
                  bins = bins, samples = labels, sample_group = grp_of, after = after)
  js <- paste0('<script>', plotly_js(), '</script><script>const D = ', toJSON(js_data, auto_unbox = TRUE, digits = NA), ';', IMPUTE_JS, '</script>')
  extra_css <- '<style>.viz-root { --series-2:#eb6834; --series-3:#1baf7a; --series-4:#eda100; } .control select { font:inherit; font-size:13px; }
.control input[type=range] { width:260px; }</style>'
  body <- paste0(decision_html, after_html, plan_html, bx_html, duke_html, graphs, js)
  write_review_page(run_dir, "07_impute", "Step 7 - Impute", status, status_class, body, extra_head = extra_css)
}

IMPUTE_JS <- '
const css = getComputedStyle(document.querySelector(".viz-root"));
const C = n => css.getPropertyValue(n).trim();
const fmt = n => Math.round(n).toLocaleString("en-US");
const SER = ["--series-1", "--series-2", "--series-3", "--series-4"];
const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
const cfg = { displaylogo: false, responsive: true };
const centers = D.hist.map((_, i) => (D.edges[i] + D.edges[i + 1]) / 2);
// fill plan
const mc = document.getElementById("mc"), it = document.getElementById("itype"); mc.value = D.mc; it.value = D.itype;
function plan() {
  const cut = +mc.value, duke = [], bot = [];
  D.groups.forEach((g, i) => { const maxm = D.counts[i] * (100 - cut) / 100; let d = 0, b = 0;
    D.rows_by_k[i].forEach((rows, j) => { const k = j + 1, v = k * rows; if (it.value === "duke" && k <= maxm) d += v; else b += v; }); duke.push(d); bot.push(b); });
  document.getElementById("mcv").textContent = cut + "%: duke when at least " + (100 - cut) + "% of a group is measured" + (it.value === "bottomx" ? " (not used for bottomx)" : "");
  const tr = (name, x, col) => ({ type: "bar", orientation: "h", name: name, y: D.groups, x: x, marker: { color: C(col), line: { width: 2, color: C("--surface-1") } },
    hovertemplate: "%{y}: %{x:,} " + name + "<extra></extra>" });
  Plotly.react("planplot", [tr("duke", duke, "--series-3"), tr("bottom-X", bot, "--series-2")], Object.assign(base(), { barmode: "stack", showlegend: true,
    legend: { orientation: "h", y: 1.15 }, margin: { l: 70, r: 20, t: 10, b: 44 }, xaxis: ax("missing values"), yaxis: { autorange: "reversed" } }), cfg);
  document.getElementById("plantable").innerHTML = "<table><thead><tr><th>Group</th><th class=num>missing</th><th class=num>duke</th><th class=num>bottom-X</th><th class=num>bottom-X %</th></tr></thead><tbody>" +
    D.groups.map((g, i) => `<tr><td>${g}</td><td class=num>${fmt(duke[i] + bot[i])}</td><td class=num>${fmt(duke[i])}</td><td class=num>${fmt(bot[i])}</td><td class=num>${(100 * bot[i] / Math.max(1, duke[i] + bot[i])).toFixed(1)}%</td></tr>`).join("") + "</tbody></table>";
}
mc.addEventListener("input", plan); it.addEventListener("change", plan); plan();
// bottom-X range
const bxr = document.getElementById("bx"); bxr.value = Math.max(0, D.bx_grid.indexOf(D.bx));
function bxdraw() {
  const i = +bxr.value, x = D.bx_grid[i], hi = D.bx_max[i];
  document.getElementById("bxv").textContent = x + "%: fill between log2 " + D.bx_min.toFixed(2) + " and " + hi.toFixed(2) + " (intensity " + Math.pow(2, D.bx_min).toPrecision(3) + " to " + Math.pow(2, hi).toPrecision(4) + ")";
  Plotly.react("bxplot", [{ type: "bar", x: centers, y: D.hist, width: 0.1, marker: { color: C("--series-1") }, hovertemplate: "log2 %{x:.2f}: %{y:,} values<extra></extra>" }],
    Object.assign(base(), { margin: { l: 70, r: 20, t: 10, b: 50 }, showlegend: false, bargap: 0.05,
      shapes: [{ type: "rect", x0: D.bx_min, x1: hi, y0: 0, y1: 1, yref: "paper", fillcolor: C("--series-2"), opacity: 0.25, line: { width: 0 } }],
      annotations: [{ x: hi, y: 0.95, yref: "paper", text: "bottom " + x + "%", showarrow: false, xanchor: "left", xshift: 4, font: { color: C("--text-primary") } }],
      xaxis: ax("log2 intensity"), yaxis: Object.assign(ax("measured values"), { rangemode: "tozero" }) }), cfg);
}
bxr.addEventListener("input", bxdraw); bxdraw();
// duke SD
Plotly.newPlot("dukeplot", D.groups.map((g, i) => ({ type: "scatter", mode: "lines+markers", name: g, x: D.bins[g].mid, y: D.bins[g].sd,
  line: { width: 2, color: C(SER[i % 4]) }, marker: { size: 8, color: C(SER[i % 4]), line: { width: 2, color: C("--surface-1") } },
  hovertemplate: g + "<br>group mean log2 %{x:.2f}<br>SD %{y:.3f}<extra></extra>" })),
  Object.assign(base(), { margin: { l: 60, r: 20, t: 10, b: 50 }, showlegend: true, legend: { orientation: "h", y: 1.12 },
    xaxis: ax("group mean log2 intensity (bin)"), yaxis: Object.assign(ax("SD used for duke fills (log2)"), { rangemode: "tozero" }) }), cfg);
// after
if (D.after) {
  const A = D.after, h = (y, name, col) => ({ type: "bar", x: centers, y: y, width: 0.1, name: name, marker: { color: C(col) }, hovertemplate: name + " log2 %{x:.2f}: %{y:,}<extra></extra>" });
  Plotly.newPlot("afterhist", [h(A.meas, "measured", "--series-1"), h(A.duke, "duke fill", "--series-3"), h(A.bot, "bottom-X fill", "--series-2")],
    Object.assign(base(), { barmode: "stack", bargap: 0, showlegend: true, legend: { orientation: "h", y: 1.1 }, margin: { l: 70, r: 20, t: 10, b: 50 },
      xaxis: ax("log2 intensity"), yaxis: Object.assign(ax("values"), { rangemode: "tozero" }) }), cfg);
  const b = (y, name, col, n) => ({ type: "bar", name: name, x: D.samples, y: y, customdata: n, marker: { color: C(col), line: { width: 2, color: C("--surface-1") } },
    hovertemplate: "%{x}<br>" + name + ": %{y:.2f}% (%{customdata:,})<extra></extra>" });
  Plotly.newPlot("afterbar", [b(A.duke_pct, "duke", "--series-3", A.duke_n), b(A.bot_pct, "bottom-X", "--series-2", A.bot_n)],
    Object.assign(base(), { barmode: "stack", showlegend: true, legend: { orientation: "h", y: 1.12 }, margin: { l: 60, r: 20, t: 10, b: 50 },
      xaxis: { categoryorder: "array", categoryarray: D.samples }, yaxis: Object.assign(ax("% of values filled"), { rangemode: "tozero" }) }), cfg);
}
'
