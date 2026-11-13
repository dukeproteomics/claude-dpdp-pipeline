# Review page for Step 11 - Pathway: the Options panel (Run / Accept); per comparison the enriched terms of the Up and the Down
# proteins (STRING, against the tested proteins) and the STRING network of the hits (clusters, what each is enriched for).
# Clicking a term highlights its proteins in the network and lists them; clicking a protein opens it in the Explore tab.

review_pathway <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  sd <- file.path(run_dir, "steps", "11_pathway")
  summ <- if (file.exists(file.path(sd, "summary.json"))) fromJSON(file.path(sd, "summary.json"), simplifyVector = FALSE) else NULL
  res <- if (file.exists(file.path(sd, "results.json"))) fromJSON(file.path(sd, "results.json"), simplifyVector = FALSE) else NULL
  dec <- Filter(function(r) isTRUE(r$settings$accepted), read_decisions(run_dir, "11_pathway"))
  acc <- if (length(dec) && !is.null(summ) && dec[[length(dec)]]$time >= summ$generated) dec[[length(dec)]] else NULL
  orgs <- unlist(summ$organisms_in_data %||% list())
  opts <- pathway_options(p, orgs)
  form <- options_form("11_pathway", opts, summ$settings %||% list(), previewed = if (!is.null(summ)) summ$settings else NULL, applied = !is.null(acc),
                       run_dir = run_dir, title = "Pathway options",
                       words = list(preview = "Run", apply = "Accept", previewing = "the last run", applied = "accepted",
                                    effect = "(asks string-db.org; about a minute per comparison)"))
  info <- if (!is.null(summ)) paste0('<section><h2>About this analysis</h2><p class="note">',
    sprintf("STRING %s (%s), species %s%s. Background: the proteins tested in each comparison. Up and Down tested separately. ",
            esc(summ$string$version), esc(summ$string$base), esc(summ$species), if (!is.null(summ$organism)) paste0(" (", esc(summ$organism), ")") else ""),
    sprintf('Excel: <code>%s</code>.', esc(summ$excel)), if (!is.null(acc)) sprintf(" Accepted %s%s.", esc(acc$time), if (!is.null(acc$analyst)) paste0(" by ", esc(acc$analyst)) else "") else "", '</p>',
    html_table(data.frame(Comparison = vapply(summ$comparisons, function(x) esc(x$name), ""),
      Tested = vapply(summ$comparisons, function(x) fmt_n(x$tested), ""), `In STRING` = vapply(summ$comparisons, function(x) fmt_n(x$mapped_background), ""),
      Up = vapply(summ$comparisons, function(x) fmt_n(x$up), ""), `Up terms` = vapply(summ$comparisons, function(x) fmt_n(x$up_terms), ""),
      Down = vapply(summ$comparisons, function(x) fmt_n(x$down), ""), `Down terms` = vapply(summ$comparisons, function(x) fmt_n(x$down_terms), ""),
      Network = vapply(summ$comparisons, function(x) sprintf("%s proteins, %s links", fmt_n(x$network_proteins), fmt_n(x$network_edges)), ""), check.names = FALSE), 2:7),
    '</section>') else '<section><h2>Not run yet</h2><p>Set the options and press <b>Run with these settings</b>.</p></section>'
  body <- paste0(form, info, if (!is.null(res)) paste0(
    '<section><h2>Enriched terms</h2><div class="control"><label for="pwcomp">Comparison</label><select id="pwcomp"></select>',
    '<span class="pwdir"><button type="button" data-dir="up" class="on">Up</button><button type="button" data-dir="down">Down</button></span>',
    '<label for="pwsrc">Source</label><select id="pwsrc"></select><input type="search" id="pwq" placeholder="filter terms" style="font:inherit;font-size:14px;padding:3px 8px"/></div>',
    '<p class="note" id="pwnote"></p>',
    '<div class="control"><label for="pwchart">Chart</label><select id="pwchart"><option value="bar">top terms, Up and Down (bars)</option><option value="dot">top terms of this direction (dots)</option></select>',
    '<label for="pwtopn">Terms</label><select id="pwtopn"><option>10</option><option selected>15</option><option>20</option><option>30</option></select></div>',
    '<div id="pworachart" style="height:520px" role="img" aria-label="Top enriched terms"></div>',
    '<div id="pwterms" style="max-height:460px;overflow:auto"></div></section>',
    '<section><h2>Ranked analysis (GSEA)</h2><p class="note" id="pwgnote"></p>',
    '<div class="control"><label for="pwgsrc">Source</label><select id="pwgsrc"></select><input type="search" id="pwgq" placeholder="filter gene sets" style="font:inherit;font-size:14px;padding:3px 8px"/>',
    '<label for="pwgtopn">Bars</label><select id="pwgtopn"><option>10</option><option selected>15</option><option>20</option></select></div>',
    '<div id="pwgchart" style="height:560px" role="img" aria-label="Top gene sets by NES"></div>',
    '<div id="pwgterms" style="max-height:420px;overflow:auto"></div>',
    '<h3 id="pwgest">Enrichment plot</h3><p class="note">Click a gene set in the table to draw its running enrichment score (the classic GSEA plot): the curve rises where the set&apos;s proteins sit among the highest-ranked (up) proteins, falls where they sit among the lowest.</p>',
    '<div id="pwgrun" style="height:380px" role="img" aria-label="GSEA running enrichment score"></div></section>',
    '<section><h2>Overview across comparisons</h2><p class="note">The top terms of every comparison side by side: what is shared and what is specific. Over-representation: colour = direction (Up red, Down blue) and -log10 FDR, size = hits. GSEA: colour = NES, size = -log10 FDR. Empty = not significant there.</p>',
    '<div class="control"><label for="pwomethod">Method</label><select id="pwomethod"><option value="ora">over-representation</option><option value="gsea">GSEA</option></select>',
    '<label for="pwosrc">Source</label><select id="pwosrc"></select><label for="pwotopn">Top per comparison</label><select id="pwotopn"><option>5</option><option selected>8</option><option>12</option><option>20</option></select></div>',
    '<div id="pwoverview" style="height:640px" role="img" aria-label="Terms across comparisons"></div></section>',
    '<section><h2>Network of the hits</h2><p class="note">STRING interactions among the Up (red) and Down (blue) proteins; size = fold change. Clusters: Louvain on the confidence scores. ',
    'Click a term above to highlight its proteins, a cluster below to highlight it, or a protein to open it in the Explore tab.</p>',
    '<div class="control"><label for="pwcolor">Colour</label><select id="pwcolor"><option value="call">Up / Down</option><option value="cluster">cluster</option></select>',
    '<label><input type="checkbox" id="pwlabels" checked/> names</label><label><input type="checkbox" id="pwiso"/> proteins without interactions</label><button type="button" id="pwclear">Clear highlight</button><a id="pwlink" target="_blank" rel="noopener">Open this network on string-db.org</a></div>',
    '<div id="pwnet" style="height:720px" role="img" aria-label="STRING network of the hits"></div><p class="note" id="pwnetnote"></p>',
    '<h3>Clusters</h3><div id="pwclusters"></div><h3 id="pwmemt">Proteins</h3><div id="pwmembers" style="max-height:420px;overflow:auto"></div></section>',
    '<style>.pwdir button { font:inherit; padding:3px 12px; border:1px solid var(--accent); background:transparent; color:var(--accent); cursor:pointer; }',
    '.pwdir button.on { background:var(--accent); color:#fff; } #pwterms table tbody tr, #pwclusters tbody tr { cursor:pointer; }',
    '#pwterms tr.sel td, #pwclusters tr.sel td { background:var(--warn-bg); font-weight:600; } #pwlink { margin-left:12px; }</style>',
    '<script>', plotly_js(), '</script><script>window.DPDP_PW = ', toJSON(list(run = basename(run_dir), comps = res, has_explore = file.exists(file.path(run_dir, "review", "10e_explore.html"))),
      auto_unbox = TRUE, digits = NA, null = "null", na = "null"), ';</script><script>', PATHWAY_JS, '</script>') else "")
  write_review_page(run_dir, "11_pathway", "Step 11 - Pathway (STRING)", if (!is.null(acc)) "Accepted" else if (!is.null(res)) "Review the results" else "Set the options",
                    if (!is.null(acc)) "done" else "waiting", body)
}

PATHWAY_JS <- '
(function () {
  const W = window.DPDP_PW; if (!W) return;
  const $ = id => document.getElementById(id), css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\\"": "&quot;" }[c]));
  const UP = "#c0392b", DOWN = "#2e6fb7", PAL = ["#2a78d6", "#e8743b", "#19a979", "#ed4a7b", "#945ecf", "#13a4b4", "#525df4", "#bf399e", "#6c8893", "#ee6868"];
  let ci = 0, dir = "up", selTerm = null, selCluster = null;
  W.comps.forEach((c, k) => $("pwcomp").appendChild(Object.assign(document.createElement("option"), { value: k, textContent: c.name })));
  const terms = () => (W.comps[ci][dir].terms || []);
  const nameOf = id => (W.comps[ci].names || {})[id] || id.replace(/^\\d+\\./, "");
  function fillSources() { const s = $("pwsrc"), have = [...new Set(terms().map(t => t.source))];
    s.innerHTML = "<option value=\\"\\">all</option>" + have.map(x => "<option>" + esc(x) + "</option>").join(""); }
  function drawTerms() {
    const c = W.comps[ci], d = c[dir], q = $("pwq").value.toLowerCase(), src = $("pwsrc").value;
    $("pwnote").textContent = (dir === "up" ? "Up" : "Down") + " in " + c.name + ": " + d.n + " proteins (" + d.mapped + " in STRING), against " + c.mapped_background + " tested. " + d.note;
    const rows = terms().map((t, k) => [t, k]).filter(([t]) => (!src || t.source === src) && (!q || (t.description + " " + t.term).toLowerCase().includes(q)));
    $("pwterms").innerHTML = rows.length ? "<table><thead><tr><th>Source</th><th>Term</th><th class=num>Hits</th><th class=num>In term</th><th class=num>Fold</th><th class=num>FDR</th><th>Proteins</th></tr></thead><tbody>" +
      rows.map(([t, k]) => { const ids = t.ids || [], names = ids.map(nameOf);
        return "<tr data-k=" + k + (selTerm === k ? " class=sel" : "") + "><td>" + esc(t.source) + "</td><td title=\\"" + esc(t.term) + "\\">" + esc(t.description) + "</td><td class=num>" + t.hits + "</td><td class=num>" + t.term_size +
          "</td><td class=num>" + (+t.fold).toFixed(2) + "</td><td class=num>" + (+t.fdr).toExponential(1) + "</td><td>" + esc(names.slice(0, 8).join(", ")) + (names.length > 8 ? ", ..." : "") + "</td></tr>"; }).join("") + "</tbody></table>"
      : "<p class=note>No terms.</p>";
  }
  $("pwterms").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; const k = +tr.dataset.k; selTerm = selTerm === k ? null : k; selCluster = null; drawTerms(); drawNet(); members(); });
  // network
  const log2fc = f => f == null ? 0 : Math.sign(f) * Math.log2(Math.abs(f));
  function highlighted() { const net = W.comps[ci].network; if (!net) return null;
    if (selTerm !== null) return new Set(terms()[selTerm].ids || []);
    if (selCluster !== null) return new Set(net.nodes.filter(n => n.cluster === selCluster).map(n => n.id));
    return null; }
  function drawNet() {
    const net = W.comps[ci].network;
    if (!net || !net.nodes || !net.nodes.length) { Plotly.purge("pwnet"); $("pwnetnote").textContent = "No network (fewer than 2 hits in STRING)."; $("pwclusters").innerHTML = ""; $("pwlink").hidden = true; return; }
    const N = net.nodes.filter(n => $("pwiso").checked || n.degree > 0), byId = {}; N.forEach((n, k) => { byId[n.id] = n; });
    const niso = net.nodes.length - net.nodes.filter(n => n.degree > 0).length;
    const H = highlighted(), ex = [], ey = [];
    (net.edges || []).forEach(e => { const a = byId[e.a], b = byId[e.b]; if (!a || !b) return; if (H && !(H.has(e.a) && H.has(e.b))) return; ex.push(a.x, b.x, null); ey.push(a.y, b.y, null); });
    const fade = [], ex2 = [], ey2 = [];
    if (H) (net.edges || []).forEach(e => { const a = byId[e.a], b = byId[e.b]; if (!a || !b || (H.has(e.a) && H.has(e.b))) return; ex2.push(a.x, b.x, null); ey2.push(a.y, b.y, null); });
    const byCluster = $("pwcolor").value === "cluster", labelAll = $("pwlabels").checked && N.length <= 150, labels = $("pwlabels").checked;
    const groups = byCluster ? [...new Set(N.map(n => n.cluster))].sort((a, b) => a - b).map(k => [k ? "cluster " + k : "no interactions", n => n.cluster === k, k ? PAL[(k - 1) % PAL.length] : C("--text-muted")])
      : [["Up", n => n.call === "Up", UP], ["Down", n => n.call === "Down", DOWN]];
    const traces = [{ type: "scatter", mode: "lines", x: ex2, y: ey2, line: { width: 0.6, color: C("--grid") }, hoverinfo: "skip", showlegend: false },
      { type: "scatter", mode: "lines", x: ex, y: ey, line: { width: 0.9, color: C("--text-muted") }, hoverinfo: "skip", showlegend: false }];
    groups.forEach(([name, f, col]) => { const s = N.filter(f); if (!s.length) return;
      traces.push({ type: "scatter", mode: labels ? "markers+text" : "markers", name: name, x: s.map(n => n.x), y: s.map(n => n.y),
        text: s.map(n => (labelAll || (H && H.has(n.id))) ? n.name : ""), hovertext: s.map(n => n.name), textposition: "top center", textfont: { size: 10, color: C("--text-secondary") },
        customdata: s.map(n => [n.acc, n.gene, n.fc, n.p, n.cluster, n.degree]),
        hovertemplate: "<b>%{hovertext}</b> %{customdata[0]}<br>FC %{customdata[2]:.2f}, p %{customdata[3]:.2g}<br>cluster %{customdata[4]}, %{customdata[5]} links<extra></extra>",
        marker: { size: s.map(n => 7 + 4 * Math.min(4, Math.abs(log2fc(n.fc)))), color: col, opacity: s.map(n => (!H || H.has(n.id)) ? 0.95 : 0.15),
                  line: { width: s.map(n => (H && H.has(n.id)) ? 2 : 0.5), color: C("--text-primary") } } }); });
    Plotly.react("pwnet", traces, { paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 },
      margin: { l: 10, r: 10, t: 10, b: 10 }, legend: { orientation: "h", y: 1.03 }, hovermode: "closest", dragmode: "pan",
      xaxis: { visible: false }, yaxis: { visible: false, scaleanchor: "x" } }, { displaylogo: false, responsive: true, scrollZoom: true });
    const pn = $("pwnet"); if (!pn.__click) { pn.__click = true; pn.on("plotly_click", ev => { const pt = ev.points && ev.points[0]; if (!pt || !pt.customdata) return;
      if (W.has_explore) window.open("10e_explore.html#" + encodeURIComponent(pt.customdata[0]), "_blank"); }); }
    $("pwnetnote").textContent = net.nodes.length + " of " + net.of + " hits" + (net.nodes.length < net.of ? " (the smallest p)" : "") + ", " + (net.edges || []).length + " interactions" +
      (niso ? "; " + niso + " without interactions" + ($("pwiso").checked ? " (grid below)" : " (hidden)") : "") + ". " +
      (N.length > 150 && $("pwlabels").checked ? "Names: shown up to 150 proteins, and for highlighted ones; hover for the others. " : "") + "Scroll to zoom, drag to move.";
    if (net.link) { $("pwlink").hidden = false; $("pwlink").href = net.link; } else $("pwlink").hidden = true;
    const cl = net.clusters || [];
    $("pwclusters").innerHTML = cl.length ? "<table><thead><tr><th>Cluster</th><th class=num>Proteins</th><th>Enriched for (top 3)</th></tr></thead><tbody>" +
      cl.map(c => "<tr data-c=" + c.cluster + (selCluster === c.cluster ? " class=sel" : "") + "><td>cluster " + c.cluster + "</td><td class=num>" + c.size + "</td><td>" + esc((c.terms || []).join("; ") || "-") + "</td></tr>").join("") + "</tbody></table>"
      : "<p class=note>No cluster of 5 or more proteins.</p>";
  }
  $("pwclusters").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; const k = +tr.dataset.c; selCluster = selCluster === k ? null : k; selTerm = null; drawTerms(); drawNet(); members(); });
  function members() {
    const net = W.comps[ci].network, H = highlighted();
    if (!H) { $("pwmemt").textContent = "Proteins"; $("pwmembers").innerHTML = "<p class=note>Click a term or a cluster to list its proteins.</p>"; return; }
    const byId = Object.assign({}, W.comps[ci].hits || {}); ((net && net.nodes) || []).forEach(n => { byId[n.id] = Object.assign({}, byId[n.id], n); });
    const ids = [...H].sort((a, b) => ((byId[a] || {}).p ?? 2) - ((byId[b] || {}).p ?? 2));
    $("pwmemt").textContent = (selTerm !== null ? terms()[selTerm].description : "Cluster " + selCluster) + ": " + ids.length + " proteins";
    $("pwmembers").innerHTML = "<table><thead><tr><th>Protein</th><th>Accession</th><th class=num>FC</th><th class=num>p</th><th>Call</th></tr></thead><tbody>" +
      ids.map(id => { const n = byId[id]; const nm = n ? n.name : nameOf(id);
        const a = n ? (W.has_explore ? "<a href=\\"10e_explore.html#" + encodeURIComponent(n.acc) + "\\" target=_blank>" + esc(n.acc) + "</a>" : esc(n.acc)) : "";
        return "<tr><td>" + esc(nm) + "</td><td>" + a + "</td><td class=num>" + (n && n.fc != null ? (+n.fc).toFixed(2) : "") + "</td><td class=num>" + (n && n.p != null ? (+n.p).toPrecision(2) : "") +
          "</td><td>" + (n ? esc(n.call) : "") + "</td></tr>"; }).join("") + "</tbody></table>";
  }
  $("pwcomp").addEventListener("change", () => { ci = +$("pwcomp").value; selTerm = null; selCluster = null; fillSources(); drawTerms(); drawNet(); members(); });
  document.querySelectorAll(".pwdir button").forEach(b => b.addEventListener("click", () => { dir = b.dataset.dir; document.querySelectorAll(".pwdir button").forEach(x => x.classList.toggle("on", x === b));
    selTerm = null; fillSources(); drawTerms(); drawNet(); members(); }));
  ["pwsrc", "pwq"].forEach(id => $(id).addEventListener("input", drawTerms));
  ["pwcolor", "pwlabels", "pwiso"].forEach(id => $(id).addEventListener("change", drawNet));
  $("pwclear").onclick = () => { selTerm = null; selCluster = null; drawTerms(); drawNet(); members(); };
  // --- publication charts: top over-represented terms ---
  const short = (t, n) => (t.length > n ? t.slice(0, n - 3) + "..." : t);
  const lfdr = v => -Math.log10(Math.max(+v, 1e-300));
  const bp = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 } }), cfg = { displaylogo: false, responsive: true };
  function oraChart() {
    const c = W.comps[ci], n = +$("pwtopn").value, src = $("pwsrc").value, pick = d => (c[d].terms || []).filter(t => !src || t.source === src).slice(0, n);
    if ($("pwchart").value === "bar") {
      const U = pick("up"), D = pick("down"), all = D.slice().reverse().map(t => ["Down", t]).concat(U.slice().reverse().map(t => ["Up", t]));
      if (!all.length) { Plotly.purge("pworachart"); return; }
      const y = all.map(([d, t], k) => short(t.description, 60) + "\\u200b".repeat(d === "Up" ? 1 : 2));
      Plotly.react("pworachart", ["Down", "Up"].map(d => { const s = all.map((x, k) => [x, k]).filter(([x]) => x[0] === d);
        return { type: "bar", orientation: "h", name: d, y: s.map(([x, k]) => y[k]), x: s.map(([x]) => (d === "Up" ? 1 : -1) * lfdr(x[1].fdr)),
          marker: { color: d === "Up" ? UP : DOWN }, customdata: s.map(([x]) => [x[1].source, x[1].hits, x[1].term_size, x[1].fold, x[1].fdr]),
          hovertemplate: "%{y}<br>%{customdata[0]}: %{customdata[1]} of %{customdata[2]}, fold %{customdata[3]}, FDR %{customdata[4]:.2g}<extra>" + d + "</extra>" }; }),
        Object.assign(bp(), { barmode: "overlay", margin: { l: 360, r: 20, t: 10, b: 50 }, legend: { orientation: "h", y: 1.06 },
          xaxis: { title: "-log10 FDR (Down left, Up right)", gridcolor: C("--grid"), zeroline: true, zerolinecolor: C("--text-muted") },
          yaxis: { categoryorder: "array", categoryarray: y, automargin: true, tickfont: { size: 11 } } }), cfg);
    } else {
      const T = pick(dir).slice().reverse(); if (!T.length) { Plotly.purge("pworachart"); return; }
      const nh = c[dir].mapped || 1;
      Plotly.react("pworachart", [{ type: "scatter", mode: "markers", x: T.map(t => t.hits / nh), y: T.map(t => short(t.description, 60)),
        marker: { size: T.map(t => 6 + 2.2 * Math.sqrt(t.hits)), color: T.map(t => lfdr(t.fdr)), colorscale: dir === "up" ? [[0, "#f6c9c2"], [1, UP]] : [[0, "#c8d9ef"], [1, DOWN]],
                  colorbar: { title: { text: "-log10 FDR" } }, line: { width: 0.5, color: C("--text-primary") } },
        customdata: T.map(t => [t.source, t.hits, t.term_size, t.fold, t.fdr]),
        hovertemplate: "%{y}<br>%{customdata[0]}: %{customdata[1]} of %{customdata[2]}, fold %{customdata[3]}, FDR %{customdata[4]:.2g}<extra></extra>" }],
        Object.assign(bp(), { margin: { l: 360, r: 20, t: 10, b: 50 }, xaxis: { title: "hit ratio (term hits / " + (dir === "up" ? "Up" : "Down") + " proteins)", gridcolor: C("--grid") },
          yaxis: { automargin: true, tickfont: { size: 11 } } }), cfg);
    }
  }
  // --- GSEA ---
  let gsel = null;
  const G = () => W.comps[ci].gsea;
  function gseaSources() { const g = G(), s = $("pwgsrc"), have = g ? [...new Set((g.terms || []).map(t => t.source))] : [];
    s.innerHTML = "<option value=\\"\\">all</option>" + have.map(x => "<option>" + esc(x) + "</option>").join(""); }
  function gseaRows() { const g = G(); if (!g) return []; const q = $("pwgq").value.toLowerCase(), src = $("pwgsrc").value;
    return (g.terms || []).filter(t => (!src || t.source === src) && (!q || (t.description + " " + t.term).toLowerCase().includes(q))); }
  function gseaDraw() {
    const g = G();
    if (!g) { $("pwgnote").textContent = "GSEA was not run (switch it on in the options)."; ["pwgchart", "pwgrun"].forEach(id => Plotly.purge(id)); $("pwgterms").innerHTML = ""; return; }
    $("pwgnote").textContent = W.comps[ci].name + ": " + g.n + " proteins ranked by sign(FC) x -log10(p). " + g.note;
    const R = gseaRows(), n = +$("pwgtopn").value, P = R.filter(t => t.NES > 0).slice(0, n), N = R.filter(t => t.NES < 0).slice(0, n);
    const bars = N.slice().reverse().concat(P.slice().reverse()).sort((a, b) => a.NES - b.NES);
    if (bars.length) Plotly.react("pwgchart", [{ type: "bar", orientation: "h", y: bars.map(t => short(t.description, 60)), x: bars.map(t => t.NES),
        marker: { color: bars.map(t => t.NES > 0 ? UP : DOWN), opacity: bars.map(t => 0.45 + 0.55 * Math.min(1, lfdr(t.fdr) / 6)) },
        customdata: bars.map(t => [t.source, t.size, t.fdr]), hovertemplate: "%{y}<br>%{customdata[0]}, %{customdata[1]} proteins<br>NES %{x:.2f}, FDR %{customdata[2]:.2g}<extra></extra>" }],
      Object.assign(bp(), { margin: { l: 360, r: 20, t: 10, b: 50 }, xaxis: { title: "normalized enrichment score (NES; up right, down left)", gridcolor: C("--grid"), zerolinecolor: C("--text-muted") },
        yaxis: { automargin: true, tickfont: { size: 11 } } }), cfg);
    else Plotly.purge("pwgchart");
    $("pwgterms").innerHTML = R.length ? "<table><thead><tr><th>Source</th><th>Gene set</th><th class=num>Size</th><th class=num>NES</th><th class=num>FDR</th><th>Leading edge</th></tr></thead><tbody>" +
      R.map(t => "<tr data-t=\\"" + esc(t.term) + "\\"" + (gsel === t.term ? " class=sel" : "") + "><td>" + esc(t.source) + "</td><td title=\\"" + esc(t.term) + "\\">" + esc(t.description) + "</td><td class=num>" + t.size +
        "</td><td class=num style=\\"color:" + (t.NES > 0 ? UP : DOWN) + "\\">" + (+t.NES).toFixed(2) + "</td><td class=num>" + (+t.fdr).toExponential(1) + "</td><td>" +
        esc((t.leading || []).slice(0, 10).map(nameOf).join(", ")) + ((t.leading || []).length > 10 ? ", ..." : "") + "</td></tr>").join("") + "</tbody></table>" : "<p class=note>No gene set passes the FDR cutoff.</p>";
    runPlot();
  }
  function runPlot() {
    const g = G(); if (!g || !gsel || !g.members || !g.members[gsel]) { Plotly.purge("pwgrun"); $("pwgest").textContent = "Enrichment plot"; return; }
    const st = g.ranked_stat, mem = new Set(g.members[gsel]), N = st.length, Nh = mem.size;
    let sumw = 0; mem.forEach(k => { sumw += Math.abs(st[k]); });
    const ys = [0]; let run = 0, best = 0;
    for (let k = 0; k < N; k++) { run += mem.has(k) ? Math.abs(st[k]) / sumw : -1 / (N - Nh); ys.push(run); if (Math.abs(run) > Math.abs(best)) best = run; }
    const t = (g.terms || []).concat(g.all_top || []).find(x => x.term === gsel) || {};
    $("pwgest").textContent = "Enrichment plot: " + (t.description || gsel) + " (NES " + (+t.NES).toFixed(2) + ", FDR " + (+t.fdr).toPrecision(2) + ")";
    const step = Math.max(1, Math.floor(N / 1500)), xs = [], yv = []; for (let k = 0; k <= N; k += step) { xs.push(k); yv.push(ys[k]); }
    const mk = [...mem].sort((a, b) => a - b);
    Plotly.react("pwgrun", [{ type: "scatter", mode: "lines", x: xs, y: yv, line: { color: best > 0 ? UP : DOWN, width: 2 }, name: "running score", hovertemplate: "rank %{x}: %{y:.3f}<extra></extra>" },
      { type: "scatter", mode: "markers", x: mk, y: mk.map(() => 0), yaxis: "y2", marker: { symbol: "line-ns-open", size: 14, color: C("--text-primary") }, name: "set members", hovertemplate: "rank %{x}<extra></extra>" },
      { type: "bar", x: xs, y: xs.map(k => st[Math.min(k, N - 1)]), yaxis: "y3", marker: { color: xs.map(k => st[Math.min(k, N - 1)] > 0 ? UP : DOWN) }, name: "ranking score", hoverinfo: "skip" }],
      Object.assign(bp(), { showlegend: false, margin: { l: 60, r: 20, t: 10, b: 50 },
        xaxis: { title: "rank (up to down)", gridcolor: C("--grid"), range: [0, N] }, yaxis: { domain: [0.42, 1], title: "enrichment score", gridcolor: C("--grid"), zerolinecolor: C("--text-muted") },
        yaxis2: { domain: [0.3, 0.4], visible: false, range: [-1, 1] }, yaxis3: { domain: [0, 0.26], title: "score", gridcolor: C("--grid") } }), cfg);
  }
  $("pwgterms").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; gsel = tr.dataset.t; gseaDraw(); $("pwgrun").scrollIntoView({ behavior: "smooth", block: "center" }); });
  ["pwgsrc", "pwgq", "pwgtopn"].forEach(id => $(id).addEventListener("input", gseaDraw));
  // --- overview across comparisons ---
  function ovSources() { const s = $("pwosrc"), have = new Set(); W.comps.forEach(c => { ["up", "down"].forEach(d => (c[d].terms || []).forEach(t => have.add(t.source))); ((c.gsea || {}).terms || []).forEach(t => have.add(t.source)); });
    s.innerHTML = "<option value=\\"\\">all</option>" + [...have].sort().map(x => "<option>" + esc(x) + "</option>").join(""); }
  function overview() {
    const m = $("pwomethod").value, src = $("pwosrc").value, n = +$("pwotopn").value, comps = W.comps.map(c => c.name);
    const cell = {}, pick = new Map();
    W.comps.forEach(c => {
      if (m === "ora") ["up", "down"].forEach(d => (c[d].terms || []).filter(t => !src || t.source === src).forEach((t, k) => {
          const key = t.term, v = cell[key + "|" + c.name]; const val = { s: (d === "up" ? 1 : -1) * lfdr(t.fdr), size: t.hits, t: t, d: d };
          if (!v || Math.abs(val.s) > Math.abs(v.s)) cell[key + "|" + c.name] = val; if (k < n) pick.set(key, t.description); }));
      else (((c.gsea || {}).terms) || []).filter(t => !src || t.source === src).forEach((t, k) => { cell[t.term + "|" + c.name] = { s: t.NES, size: lfdr(t.fdr), t: t };
          if (k < n) pick.set(t.term, t.description); });
    });
    const terms = [...pick.keys()];
    if (!terms.length) { Plotly.purge("pwoverview"); return; }
    // order terms by the comparison where they are strongest
    terms.sort((a, b) => { const fa = comps.findIndex(c => cell[a + "|" + c]), fb = comps.findIndex(c => cell[b + "|" + c]); return fb - fa || (pick.get(b) < pick.get(a) ? -1 : 1); });
    const X = [], Y = [], Z = [], S = [], H = [];
    terms.forEach(t => comps.forEach(c => { const v = cell[t + "|" + c]; if (!v) return; X.push(c); Y.push(short(pick.get(t), 60)); Z.push(v.s);
      S.push(m === "ora" ? 6 + 2.2 * Math.sqrt(v.size) : 6 + 3 * Math.min(8, v.size));
      H.push(m === "ora" ? (v.d === "up" ? "Up" : "Down") + ", " + v.t.hits + " of " + v.t.term_size + ", FDR " + (+v.t.fdr).toPrecision(2) : "NES " + (+v.t.NES).toFixed(2) + ", FDR " + (+v.t.fdr).toPrecision(2)); }));
    const lim = Math.max(...Z.map(Math.abs), 1);
    $("pwoverview").style.height = Math.max(360, 110 + 22 * terms.length) + "px";
    Plotly.react("pwoverview", [{ type: "scatter", mode: "markers", x: X, y: Y, text: H, hovertemplate: "%{y}<br>%{x}: %{text}<extra></extra>",
      marker: { size: S, color: Z, cmin: -lim, cmax: lim, colorscale: [[0, DOWN], [0.5, "#f2f2f2"], [1, UP]], line: { width: 0.5, color: C("--text-primary") },
                colorbar: { title: { text: m === "ora" ? "-log10 FDR (Down -, Up +)" : "NES" } } } }],
      Object.assign(bp(), { margin: { l: 360, r: 20, t: 10, b: 120 }, xaxis: { type: "category", categoryorder: "array", categoryarray: comps, tickangle: -30, gridcolor: C("--grid") },
        yaxis: { type: "category", categoryorder: "array", categoryarray: terms.map(t => short(pick.get(t), 60)).reverse(), automargin: true, tickfont: { size: 11 }, gridcolor: C("--grid") } }), cfg);
  }
  ["pwomethod", "pwosrc", "pwotopn"].forEach(id => $(id).addEventListener("input", overview));
  ["pwchart", "pwtopn"].forEach(id => $(id).addEventListener("input", oraChart));
  $("pwsrc").addEventListener("input", oraChart); $("pwq").addEventListener("input", () => {});
  $("pwcomp").addEventListener("change", () => { gsel = null; gseaSources(); gseaDraw(); oraChart(); });
  document.querySelectorAll(".pwdir button").forEach(b => b.addEventListener("click", oraChart));
  fillSources(); drawTerms(); drawNet(); members(); oraChart(); gseaSources(); gseaDraw(); ovSources(); overview();
})();
'
