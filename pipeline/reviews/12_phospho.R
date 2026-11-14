# Review page for Step 12 - Phospho motifs and kinases: the Options panel (Run / Accept); per comparison and direction the Kinase
# Library kinase enrichment (volcano + table, Ser/Thr or Tyr), the enrichment logo (residue x position), the motif-x motifs, and the
# sites with their windows and top kinases.

review_phospho <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  sd <- file.path(run_dir, "steps", "12_phospho")
  summ <- if (file.exists(file.path(sd, "summary.json"))) fromJSON(file.path(sd, "summary.json"), simplifyVector = FALSE) else NULL
  res <- if (file.exists(file.path(sd, "results.json"))) fromJSON(file.path(sd, "results.json"), simplifyVector = FALSE) else NULL
  dec <- Filter(function(r) isTRUE(r$settings$accepted), read_decisions(run_dir, "12_phospho"))
  acc <- if (length(dec) && !is.null(summ) && dec[[length(dec)]]$time >= summ$generated) dec[[length(dec)]] else NULL
  form <- options_form("12_phospho", phospho_options(p), summ$settings %||% list(), previewed = if (!is.null(summ)) summ$settings else NULL, applied = !is.null(acc),
                       run_dir = run_dir, title = "Phospho options",
                       words = list(preview = "Run", apply = "Accept", previewing = "the last run", applied = "accepted",
                                    effect = "(UniProt sequences, Kinase Library, logos, motif-x; a minute or two)"))
  info <- if (!is.null(summ)) paste0('<section><h2>About this analysis</h2><p class="note">',
    sprintf("Kinase Library %s (percentile rank, top %s kinases per site; Fisher&apos;s exact test, BH, Up or Down v the unchanged sites). Logos and motif-x: Up or Down v every tested site. Windows: %s of %s sites (UniProt). Excel (every table of every comparison): %s.",
            esc(summ$kinase_library), esc(summ$settings$phospho_kl_rank), fmt_n(summ$windows$with_window), fmt_n(summ$windows$sites),
            { u <- run_file_link(run_dir, summ$excel); if (nzchar(u)) sprintf('<a href="%s" download>download %s</a> (<code>%s</code>)', esc(u), esc(basename(summ$excel)), esc(summ$excel)) else sprintf("<code>%s</code>", esc(summ$excel)) }),
    if (!is.null(acc)) sprintf(" Accepted %s%s.", esc(acc$time), if (!is.null(acc$analyst)) paste0(" by ", esc(acc$analyst)) else "") else "", '</p>',
    '<div id="phabout">', html_table(data.frame(Comparison = vapply(summ$comparisons, function(x) esc(x$name), ""), Sites = vapply(summ$comparisons, function(x) fmt_n(x$with_window), ""),
      Up = vapply(summ$comparisons, function(x) fmt_n(x$up), ""), Down = vapply(summ$comparisons, function(x) fmt_n(x$down), ""),
      `Up: kinases (adj p)` = vapply(summ$comparisons, function(x) esc(paste(unlist(x$up_kinases), collapse = ", ")), ""),
      `Down: kinases (adj p)` = vapply(summ$comparisons, function(x) esc(paste(unlist(x$down_kinases), collapse = ", ")), ""), check.names = FALSE), 2:4), '</div>',
    '<p class="note">The Kinase Library describes human kinases; for other organisms read the kinases as motif classes (e.g. CDK / MAPK = proline-directed, PKA / AKT = basophilic, CK2 = acidophilic).</p></section>')
  else '<section><h2>Not run yet</h2><p>Set the options and press <b>Run with these settings</b>.</p></section>'
  body <- paste0(form, info, if (!is.null(res)) paste0(
    '<section><h2>Kinases</h2><div class="control"><label for="phcomp">Comparison</label><select id="phcomp"></select>',
    '<span class="phseg" id="phdir"><button type="button" data-v="up" class="on">Up</button><button type="button" data-v="down">Down</button></span>',
    '<span class="phseg" id="phkt"><button type="button" data-v="ser_thr" class="on">Ser/Thr</button><button type="button" data-v="tyrosine">Tyr</button></span></div>',
    '<p class="note" id="phnote"></p><div id="phvol" style="height:520px" role="img" aria-label="Kinase enrichment volcano"></div>',
    '<div id="phkin" style="max-height:380px;overflow:auto"></div></section>',
    '<section><h2>Enrichment logo</h2><p class="note" id="phlogonote"></p><div id="phlogo" style="height:560px" role="img" aria-label="Residue enrichment by position"></div></section>',
    '<section><h2>Motifs (motif-x)</h2><div id="phmotifs"></div></section>',
    '<section><h2>Sites</h2><div class="control"><input type="search" id="phq" placeholder="filter sites (gene, accession, site, kinase)" style="font:inherit;font-size:14px;padding:3px 8px;width:22em"/>',
    '<label><input type="checkbox" id="phsig" checked/> Up / Down only</label><span class="note" id="phcount"></span></div><div id="phsites" style="max-height:520px;overflow:auto"></div></section>',
    '<style>.phseg button { font:inherit; padding:3px 12px; border:1px solid var(--accent); background:transparent; color:var(--accent); cursor:pointer; } .phseg { margin-left:8px; }',
    '.phseg button.on { background:var(--accent); color:#fff; } .phwin { font-family:monospace; letter-spacing:0.5px; } .phwin b { color:#c0392b; }</style>',
    '<script>', plotly_js(), '</script><script>window.DPDP_PH = ', toJSON(list(run = basename(run_dir), fdr = summ$settings$phospho_kl_fdr %||% 0.1, comps = res,
      has_explore = file.exists(file.path(run_dir, "review", "10e_explore.html"))), auto_unbox = TRUE, digits = NA, null = "null", na = "null"),
    ';</script><script>', PHOSPHO_JS, '</script>') else "")
  write_review_page(run_dir, "12_phospho", "Step 12 - Phospho motifs and kinases", if (!is.null(acc)) "Accepted" else if (!is.null(res)) "Review the results" else "Set the options",
                    if (!is.null(acc)) "done" else "waiting", body)
}

PHOSPHO_JS <- '
(function () {
  const W = window.DPDP_PH; if (!W) return;
  const $ = id => document.getElementById(id), css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\\"": "&quot;" }[c]));
  const FAM = ["#2a78d6", "#e8743b", "#19a979", "#ed4a7b", "#945ecf", "#13a4b4", "#525df4", "#bf399e", "#6c8893", "#ee6868", "#a8a032", "#7d5c3e"];
  let ci = 0, dir = "up", kt = "ser_thr";
  W.comps.forEach((c, k) => $("phcomp").appendChild(Object.assign(document.createElement("option"), { value: k, textContent: c.name })));
  const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } });
  const cfg = { displaylogo: false, responsive: true };
  function kinases() {
    const c = W.comps[ci], K = (c.kinases || {})[kt] || {}, rows = K[dir] || [];
    $("phnote").textContent = (dir === "up" ? "Up" : "Down") + " " + (kt === "ser_thr" ? "Ser/Thr" : "Tyr") + " sites in " + c.name + ": " + (dir === "up" ? K.up_sites : K.down_sites) + ", against " + K.unchanged_sites + " unchanged sites. " +
      (K.note ? K.note + "." : K.error ? "Error: " + K.error : rows.length ? "" : "No kinase table.");
    if (!rows.length) { Plotly.purge("phvol"); $("phkin").innerHTML = ""; return; }
    const fams = [...new Set(rows.map(r => r.family || "other"))].sort(), sig = r => r.fisher_adj_pval <= W.fdr;
    const lab = rows.filter(r => sig(r) && r.log2_freq_factor > 0).slice(0, 15).map(r => r.kinase);
    const tr = fams.map((f, k) => { const s = rows.filter(r => (r.family || "other") === f);
      return { type: "scatter", mode: "markers+text", name: f, x: s.map(r => r.log2_freq_factor), y: s.map(r => -Math.log10(Math.max(r.fisher_adj_pval, 1e-300))),
        text: s.map(r => lab.includes(r.kinase) ? r.kinase : ""), textposition: "top center", textfont: { size: 10 },
        customdata: s.map(r => [r.kinase, r.fg_counts, r.fg_total, r.bg_counts, r.bg_total, r.fisher_adj_pval]), meta: { cd: ["kinase", "sites", "of sites", "unchanged sites", "of unchanged", "adj p"] },
        hovertemplate: "<b>%{customdata[0]}</b> (" + f + ")<br>%{customdata[1]} of %{customdata[2]} sites v %{customdata[3]} of %{customdata[4]}<br>log2 frequency factor %{x:.2f}, adj p %{customdata[5]:.2g}<extra></extra>",
        marker: { size: 8, color: FAM[k % FAM.length], opacity: s.map(r => sig(r) ? 0.95 : 0.25) } }; });
    Plotly.react("phvol", tr, Object.assign(base(), { margin: { l: 60, r: 10, t: 10, b: 50 }, legend: { orientation: "h", y: 1.1 },
      shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, y0: -Math.log10(W.fdr), y1: -Math.log10(W.fdr), line: { dash: "dot", color: C("--text-muted") } },
               { type: "line", yref: "paper", y0: 0, y1: 1, x0: 0, x1: 0, line: { dash: "dot", color: C("--text-muted") } }],
      xaxis: { title: "log2 frequency factor (substrate share in " + dir + " v unchanged sites)", gridcolor: C("--grid"), zeroline: false },
      yaxis: { title: "-log10 adjusted p", gridcolor: C("--grid"), zeroline: false } }), cfg);
    $("phkin").innerHTML = "<table><thead><tr><th>Kinase</th><th>Family</th><th class=num>Sites</th><th class=num>Unchanged sites</th><th class=num>log2 FF</th><th class=num>adj p</th></tr></thead><tbody>" +
      rows.filter(r => r.log2_freq_factor > 0).slice(0, 60).map(r => "<tr" + (sig(r) ? "" : " style=\\"opacity:0.55\\"") + "><td><b>" + esc(r.kinase) + "</b></td><td>" + esc(r.family) + "</td><td class=num>" + r.fg_counts + " / " + r.fg_total +
        " (" + (+r.fg_percent).toFixed(1) + "%)</td><td class=num>" + r.bg_counts + " / " + r.bg_total + "</td><td class=num>" + (+r.log2_freq_factor).toFixed(2) + "</td><td class=num>" + (+r.fisher_adj_pval).toExponential(1) + "</td></tr>").join("") + "</tbody></table>";
  }
  function logo() {
    const L = W.comps[ci][dir].logo;
    if (!L) { Plotly.purge("phlogo"); $("phlogonote").textContent = "Too few sites for a logo."; return; }
    const z = L.score, cap = Math.max(L.threshold * 2, 5), txt = z.map((row, i) => row.map(v => Math.abs(v) >= L.threshold ? L.residues[i] : ""));
    Plotly.react("phlogo", [{ type: "heatmap", z: z, x: L.positions.map(p => (p > 0 ? "+" : "") + p), y: L.residues, zmid: 0, zmin: -cap, zmax: cap,
      colorscale: [[0, "#2e6fb7"], [0.5, "#ffffff"], [1, "#c0392b"]], text: txt, texttemplate: "%{text}", textfont: { size: 12, color: "#111" },
      customdata: L.freq, hovertemplate: "%{y} at %{x}: %{z:.1f} (log10 p; + over-represented)<br>%{customdata:.1%} of the sites<extra></extra>",
      colorbar: { title: { text: "log10 p" } } }],
      Object.assign(base(), { margin: { l: 40, r: 10, t: 10, b: 50 }, xaxis: { title: "position (site = 0)", type: "category" }, yaxis: { type: "category", autorange: "reversed" } }), cfg);
    $("phlogonote").textContent = (dir === "up" ? "Up" : "Down") + " S/T sites (" + L.n_fg + ") v every tested S/T site (" + L.n_bg + "): binomial significance of each residue at each position (red over-, blue under-represented); letters mark |log10 p| >= " + L.threshold + " (Bonferroni 0.05).";
  }
  function motifs() {
    const M = W.comps[ci][dir].motifs || [];
    $("phmotifs").innerHTML = M.length ? "<table><thead><tr><th>Motif</th><th>Central</th><th class=num>Score</th><th class=num>Sites</th><th class=num>Background</th><th class=num>Fold</th></tr></thead><tbody>" +
      M.map(m => "<tr><td class=phwin>" + esc(m.motif).replace(/([a-z])/, "<b>$1</b>") + "</td><td>" + m.central + "</td><td class=num>" + m.score + "</td><td class=num>" + m.fg_matches + " / " + m.fg_size +
        "</td><td class=num>" + m.bg_matches + " / " + m.bg_size + "</td><td class=num>" + (+m.fold).toFixed(2) + "</td></tr>").join("") + "</tbody></table><p class=note>Order of discovery (motif-x removes the sites of each motif before searching again); Sites / Background = matches among the remaining sites.</p>"
      : "<p class=note>No motif (lower the minimum occurrences for small sets).</p>";
  }
  const cache = {};
  function sites() {
    const c = W.comps[ci];
    if (!cache[c.site_file]) { if (!location.protocol.startsWith("http")) { $("phcount").textContent = "The site table loads through the review server or the web app."; return; }
      $("phcount").textContent = "loading..."; const want = ci;
      fetch(c.site_file).then(r => r.json()).then(t => { cache[c.site_file] = t; if (ci === want) sites(); }).catch(() => { $("phcount").textContent = "Could not load the site table."; }); return; }
    const S = cache[c.site_file], q = $("phq").value.toLowerCase(), only = $("phsig").checked;
    const rows = S.filter(s => (!only || s.call) && (!q || [s.gene, s.acc, s.site, s.top_kinases].join(" ").toLowerCase().includes(q)));
    $("phcount").textContent = rows.length + " of " + S.length + " sites" + (rows.length > 300 ? " (first 300 shown)" : "");
    const win = w => esc(w.slice(0, 7)) + "<b>" + esc(w[7]) + "</b>" + esc(w.slice(8));
    $("phsites").innerHTML = "<table><thead><tr><th>Protein</th><th>Site</th><th>Window</th><th class=num>FC</th><th class=num>p</th><th>Call</th><th>Top kinases (percentile)</th></tr></thead><tbody>" +
      rows.slice(0, 300).map(s => "<tr><td>" + esc(s.gene || s.acc) + " " + (W.has_explore ? "<a href=\\"10e_explore.html#" + encodeURIComponent(s.acc) + "\\" target=_blank>" + esc(s.acc) + "</a>" : esc(s.acc)) + "</td><td>" + esc(s.site) +
        "</td><td class=phwin>" + win(s.window) + "</td><td class=num>" + (s.fc == null ? "" : (+s.fc).toFixed(2)) + "</td><td class=num>" + (s.p == null ? "" : (+s.p).toPrecision(2)) + "</td><td>" + esc(s.call) + "</td><td>" + esc(s.top_kinases) + "</td></tr>").join("") + "</tbody></table>";
  }
  // downloads (user, 2026-11-13): every table in full (all kinases, all filtered sites), not only the rows shown
  const D_ = () => (dir === "up" ? "Up" : "Down"), K_ = () => (kt === "ser_thr" ? "Ser-Thr" : "Tyr");
  if (window.dpdpCsvButton) {
    dpdpCsvButton("phabout", "Summary");
    dpdpCsvButton("phkin", () => "Kinases " + W.comps[ci].name + " " + D_() + " " + K_(), () => { const rows = (((W.comps[ci].kinases || {})[kt] || {})[dir]) || [];
      const cols = rows.length ? Object.keys(rows[0]) : []; return { cols: ["comparison", "direction", "kinase type"].concat(cols), rows: rows.map(r => [W.comps[ci].name, D_(), K_()].concat(cols.map(k => r[k]))) }; });
    dpdpCsvButton("phmotifs", () => "Motifs " + W.comps[ci].name + " " + D_(), () => { const M = W.comps[ci][dir].motifs || [], cols = M.length ? Object.keys(M[0]) : [];
      return { cols: ["comparison", "direction"].concat(cols), rows: M.map(m => [W.comps[ci].name, D_()].concat(cols.map(k => m[k]))) }; });
    dpdpCsvButton("phsites", () => "Sites " + W.comps[ci].name + ($("phsig").checked ? " Up-Down" : "") + ($("phq").value ? " " + $("phq").value : ""), () => {
      const S = cache[W.comps[ci].site_file] || [], q = $("phq").value.toLowerCase(), only = $("phsig").checked;
      const rows = S.filter(s => (!only || s.call) && (!q || [s.gene, s.acc, s.site, s.top_kinases].join(" ").toLowerCase().includes(q)));
      return { cols: ["comparison", "accession", "gene", "site", "window", "FC", "p", "call", "peptides", "top kinases (percentile)"],
        rows: rows.map(s => [W.comps[ci].name, s.acc, s.gene, s.site, s.window, s.fc, s.p, s.call, s.peptides, s.top_kinases]) }; });
  }
  const all = () => { kinases(); logo(); motifs(); sites(); };
  $("phcomp").addEventListener("change", () => { ci = +$("phcomp").value; all(); });
  ["phdir", "phkt"].forEach(id => $(id).querySelectorAll("button").forEach(b => b.addEventListener("click", () => {
    $(id).querySelectorAll("button").forEach(x => x.classList.toggle("on", x === b)); if (id === "phdir") dir = b.dataset.v; else kt = b.dataset.v; all(); })));
  ["phq", "phsig"].forEach(id => $(id).addEventListener("input", sites));
  all();
})();
'
