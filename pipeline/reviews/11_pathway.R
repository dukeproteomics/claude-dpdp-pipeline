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
    sprintf('Excel (every table of every comparison): <a href="%s" download>download %s</a> <span class="note">(<code>%s</code>)</span>.', esc(run_file_link(run_dir, summ$excel)), esc(basename(summ$excel)), esc(summ$excel)), if (!is.null(acc)) sprintf(" Accepted %s%s.", esc(acc$time), if (!is.null(acc$analyst)) paste0(" by ", esc(acc$analyst)) else "") else "", '</p>',
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
    '<section><h2>GO with the hierarchy (topGO, grouped)</h2><p class="note" id="pwgonote"></p>',
    '<div class="control"><label for="pwgoont">Ontology</label><select id="pwgoont"><option value="BP">biological process</option><option value="MF">molecular function</option><option value="CC">cellular component</option></select>',
    '<span class="pwdir" id="pwgodir"><button type="button" data-dir="up" class="on">Up</button><button type="button" data-dir="down">Down</button></span>',
    '<label for="pwgochart">Chart</label><select id="pwgochart"><option value="bar">groups (bars)</option><option value="tree">groups and terms (treemap)</option></select></div>',
    '<div id="pwgochartdiv" style="height:480px" role="img" aria-label="GO term groups"></div>',
    '<div id="pwgoterms" style="max-height:460px;overflow:auto"></div>',
    '<p class="note">weight01 p: topGO&apos;s test that discounts a parent term&apos;s signal coming from a more specific child term (not corrected for multiple testing, as its authors advise); classic p: the plain Fisher test for comparison. Groups: terms with GOSemSim similarity above the cutoff; the first (most significant) term represents the group. Click a term to highlight its proteins in the network below.</p></section>',
    '<section><h2>GO slim profile</h2><p class="note" id="pwslnote"></p>',
    '<div class="control"><label for="pwslont">Ontology</label><select id="pwslont"><option value="BP">biological process</option><option value="MF">molecular function</option><option value="CC">cellular component</option></select>',
    '<label for="pwslorder">Order</label><select id="pwslorder"><option value="tested">share of the tested proteins</option><option value="sig">most enriched in Up or Down</option></select>',
    '<label for="pwsltopn">Categories</label><select id="pwsltopn"><option>15</option><option selected>25</option><option>40</option><option value="999">all</option></select></div>',
    '<div id="pwslchart" style="height:640px" role="img" aria-label="GO slim profile"></div>',
    '<div id="pwslterms" style="max-height:420px;overflow:auto"></div>',
    '<p class="note">Broad GO categories (the GO Consortium&apos;s generic slim; a protein counts in every category above its own terms). Bars: the share of the annotated tested, Up and Down proteins in each category; a category can hold more than one bar&apos;s worth because proteins sit in several. * = Up or Down enriched v the tested proteins (Fisher, BH FDR &le; 0.05). Click a category for its drill-down.</p></section>',
    '<section><h2>Overview across comparisons</h2><p class="note">The top terms of every comparison side by side: what is shared and what is specific. Over-representation: colour = direction (Up red, Down blue) and -log10 FDR, size = hits. GSEA: colour = NES, size = -log10 FDR. Empty = not significant there.</p>',
    '<div class="control"><label for="pwomethod">Method</label><select id="pwomethod"><option value="ora">over-representation</option><option value="gsea">GSEA</option><option value="go">GO (topGO, biological process)</option></select>',
    '<label for="pwosrc">Source</label><select id="pwosrc"></select><label for="pwotopn">Top per comparison</label><select id="pwotopn"><option>5</option><option selected>8</option><option>12</option><option>20</option></select></div>',
    '<div id="pwoverview" style="height:640px" role="img" aria-label="Terms across comparisons"></div></section>',
    '<section><h2>Network of the hits</h2><p class="note">STRING interactions among the Up (red) and Down (blue) proteins; size = fold change. Clusters: Louvain on the confidence scores. ',
    'Click a term above to highlight its proteins, a cluster below to highlight it, or a protein to open it in the Explore tab.</p>',
    '<div class="control"><label for="pwcolor">Colour</label><select id="pwcolor"><option value="call">Up / Down</option><option value="cluster">cluster</option></select>',
    '<label><input type="checkbox" id="pwlabels" checked/> names</label><label><input type="checkbox" id="pwiso"/> proteins without interactions</label><button type="button" id="pwclear">Clear highlight</button><a id="pwlink" target="_blank" rel="noopener">Open this network on string-db.org</a></div>',
    '<div id="pwnet" style="height:720px" role="img" aria-label="STRING network of the hits"></div><p class="note" id="pwnetnote"></p>',
    '<h3>Clusters</h3><div id="pwclusters"></div><h3 id="pwmemt">Proteins</h3><div id="pwmembers" style="max-height:420px;overflow:auto"></div></section>',
    '<section id="pwdrill"><h2 id="pwdrt">Term drill-down</h2><p class="note" id="pwdrnote">Click a term (over-representation, GSEA, GO or GO slim) or a network cluster: its tested proteins are shown on the volcano and in a heatmap.</p>',
    '<div class="control"><label><input type="checkbox" id="pwdrlab" checked/> names on the volcano (20 smallest p)</label><label for="pwdrscale">Heatmap</label><select id="pwdrscale"><option value="z">z-score per protein</option><option value="log2">log2 intensity</option></select>',
    '<label for="pwdrrows">Rows</label><select id="pwdrrows"><option value="hits">hits only</option><option value="50" selected>50 smallest p</option><option value="100">100 smallest p</option><option value="300">300 smallest p</option></select></div>',
    '<div id="pwdrvol" style="height:520px" role="img" aria-label="Volcano with the term proteins highlighted"></div>',
    '<div id="pwdrheat" style="height:520px" role="img" aria-label="Heatmap of the term proteins"></div>',
    '<h3 id="pwdrtt">Proteins of the term</h3><div id="pwdrtab" style="max-height:420px;overflow:auto"></div></section>',
    '<style>.pwdir button { font:inherit; padding:3px 12px; border:1px solid var(--accent); background:transparent; color:var(--accent); cursor:pointer; }',
    '.pwdir button.on { background:var(--accent); color:#fff; } .gotog { font:inherit; font-size:11px; padding:0 5px; border:1px solid var(--rule); border-radius:8px; background:transparent; color:var(--text-secondary); cursor:pointer; } #pwterms table tbody tr, #pwclusters tbody tr { cursor:pointer; }',
    '#pwterms tr.sel td, #pwclusters tr.sel td { background:var(--warn-bg); font-weight:600; } #pwlink { margin-left:12px; } #pwclusters td, #pwterms td, #pwgterms td, #pwgoterms td { white-space:normal; }</style>',
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
  $("pwterms").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; const k = +tr.dataset.k; window.__pwGoSel = null; selTerm = selTerm === k ? null : k; selCluster = null; drawTerms(); drawNet(); members();
    const t = selTerm !== null ? terms()[selTerm] : null; window.__pwDrill(t ? { term: t.term, label: t.description, ids: t.ids, src: t.source + ", " + (dir === "up" ? "Up" : "Down") } : null); });
  // network
  const log2fc = f => f == null ? 0 : Math.sign(f) * Math.log2(Math.abs(f));
  function highlighted() { const net = W.comps[ci].network; if (!net) return null;
    if (window.__pwGoSel) return new Set(window.__pwGoSel.ids || []);
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
  $("pwclusters").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; const k = +tr.dataset.c; selCluster = selCluster === k ? null : k; selTerm = null; drawTerms(); drawNet(); members();
    window.__pwDrill(selCluster !== null ? { label: "network cluster " + k, ids: W.comps[ci].network.nodes.filter(n => n.cluster === k).map(n => n.id), src: "STRING network" } : null); });
  function members() {
    const net = W.comps[ci].network, H = highlighted();
    if (!H) { $("pwmemt").textContent = "Proteins"; $("pwmembers").innerHTML = "<p class=note>Click a term or a cluster to list its proteins.</p>"; return; }
    const byId = Object.assign({}, W.comps[ci].hits || {}); ((net && net.nodes) || []).forEach(n => { byId[n.id] = Object.assign({}, byId[n.id], n); });
    const ids = [...H].sort((a, b) => ((byId[a] || {}).p ?? 2) - ((byId[b] || {}).p ?? 2));
    $("pwmemt").textContent = (window.__pwGoSel ? window.__pwGoSel.description : selTerm !== null ? terms()[selTerm].description : "Cluster " + selCluster) + ": " + ids.length + " proteins";
    $("pwmembers").innerHTML = "<table><thead><tr><th>Protein</th><th>Accession</th><th class=num>FC</th><th class=num>p</th><th>Call</th></tr></thead><tbody>" +
      ids.map(id => { const n = byId[id]; const nm = n ? n.name : nameOf(id);
        const a = n ? (W.has_explore ? "<a href=\\"10e_explore.html#" + encodeURIComponent(n.acc) + "\\" target=_blank>" + esc(n.acc) + "</a>" : esc(n.acc)) : "";
        return "<tr><td>" + esc(nm) + "</td><td>" + a + "</td><td class=num>" + (n && n.fc != null ? (+n.fc).toFixed(2) : "") + "</td><td class=num>" + (n && n.p != null ? (+n.p).toPrecision(2) : "") +
          "</td><td>" + (n ? esc(n.call) : "") + "</td></tr>"; }).join("") + "</tbody></table>";
  }
  $("pwcomp").addEventListener("change", () => { ci = +$("pwcomp").value; selTerm = null; selCluster = null; fillSources(); drawTerms(); drawNet(); members(); });
  document.querySelectorAll(".pwdir:not(#pwgodir) button").forEach(b => b.addEventListener("click", () => { dir = b.dataset.dir; document.querySelectorAll(".pwdir:not(#pwgodir) button").forEach(x => x.classList.toggle("on", x === b));
    selTerm = null; fillSources(); drawTerms(); drawNet(); members(); }));
  ["pwsrc", "pwq"].forEach(id => $(id).addEventListener("input", drawTerms));
  ["pwcolor", "pwlabels", "pwiso"].forEach(id => $(id).addEventListener("change", drawNet));
  $("pwclear").onclick = () => { window.__pwGoSel = null; selTerm = null; selCluster = null; drawTerms(); drawNet(); members(); if (window.__pwGoDraw) window.__pwGoDraw(); window.__pwDrill(null); };
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
          marker: { color: d === "Up" ? UP : DOWN }, customdata: s.map(([x]) => [x[1].source, x[1].hits, x[1].term_size, x[1].fold, x[1].fdr, x[1].term, x[1].description]), meta: { cd: ["source", "hits", "term size", "fold", "FDR", "term", "description"] },
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
        customdata: T.map(t => [t.source, t.hits, t.term_size, t.fold, t.fdr, t.term, t.description]), meta: { cd: ["source", "hits", "term size", "fold", "FDR", "term", "description"] },
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
        customdata: bars.map(t => [t.source, t.size, t.fdr, t.term, t.description]), meta: { cd: ["source", "size", "FDR", "term", "description"] }, hovertemplate: "%{y}<br>%{customdata[0]}, %{customdata[1]} proteins<br>NES %{x:.2f}, FDR %{customdata[2]:.2g}<extra></extra>" }],
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
  $("pwgterms").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; gsel = tr.dataset.t; gseaDraw(); $("pwgrun").scrollIntoView({ behavior: "smooth", block: "center" });
    const t = (G().terms || []).find(x => x.term === gsel); window.__pwDrill({ term: gsel, label: t ? t.description : gsel, ids: t ? t.leading : [], src: "GSEA" + (t ? ", NES " + (+t.NES).toFixed(2) : "") }); });
  ["pwgsrc", "pwgq", "pwgtopn"].forEach(id => $(id).addEventListener("input", gseaDraw));
  // --- GO (topGO, grouped) ---
  let godir = "up", open_ = new Set();
  function goData() { const g = W.comps[ci].go; return g ? (g[$("pwgoont").value] || {})[godir] : null; }
  function goDraw() {
    const d = goData();
    if (!W.comps[ci].go) { $("pwgonote").textContent = "GO (topGO) was not run (switch it on in the options)."; Plotly.purge("pwgochartdiv"); $("pwgoterms").innerHTML = ""; return; }
    const T = (d && d.terms) || []; $("pwgonote").textContent = W.comps[ci].name + ", " + (godir === "up" ? "Up" : "Down") + ": " + ((d && d.note) || "");
    if (!T.length) { Plotly.purge("pwgochartdiv"); $("pwgoterms").innerHTML = "<p class=note>No term.</p>"; return; }
    const groups = {}; T.forEach(t => { (groups[t.group] = groups[t.group] || []).push(t); });
    const G2 = Object.keys(groups).map(k => groups[k]).sort((a, b) => a[0].p - b[0].p), col = godir === "up" ? UP : DOWN;
    if ($("pwgochart").value === "bar") {
      const R = G2.slice(0, 25).reverse();
      Plotly.react("pwgochartdiv", [{ type: "bar", orientation: "h", y: R.map(g => short(g[0].description, 55) + (g.length > 1 ? " (+" + (g.length - 1) + ")" : "")), x: R.map(g => lfdr(g[0].p)),
        marker: { color: col }, customdata: R.map(g => [g[0].term, g[0].hits, g[0].annotated, g[0].fold, g[0].p, g.slice(1).map(t => t.description).join("; "), g[0].description]), meta: { cd: ["GO id", "hits", "annotated", "fold", "weight01 p", "similar terms in the group", "description"] },
        hovertemplate: "%{y}<br>%{customdata[0]}: %{customdata[1]} of %{customdata[2]}, fold %{customdata[3]}, weight01 p %{customdata[4]:.2g}<br>%{customdata[5]}<extra></extra>" }],
        Object.assign(bp(), { margin: { l: 360, r: 20, t: 10, b: 50 }, xaxis: { title: "-log10 p (topGO weight01)", gridcolor: C("--grid") }, yaxis: { automargin: true, tickfont: { size: 11 } } }), cfg);
    } else {
      const ids = [], labels = [], parents = [], values = [], colors = [];
      G2.forEach((g, k) => { const gid = "g" + k; ids.push(gid); labels.push(short(g[0].description, 40)); parents.push(""); values.push(0); colors.push(PAL[k % PAL.length]);
        g.forEach(t => { ids.push(gid + "/" + t.term); labels.push(short(t.description, 40)); parents.push(gid); values.push(lfdr(t.p)); colors.push(PAL[k % PAL.length]); }); });
      Plotly.react("pwgochartdiv", [{ type: "treemap", ids, labels, parents, values, branchvalues: "remainder", marker: { colors }, textinfo: "label",
        hovertemplate: "%{label}<br>-log10 p %{value:.2f}<extra></extra>", tiling: { pad: 2 } }], Object.assign(bp(), { margin: { l: 4, r: 4, t: 4, b: 4 } }), cfg);
    }
    const row = (t, rep, n) => "<tr data-term=\\"" + t.term + "\\"" + (window.__pwGoSel && window.__pwGoSel.term === t.term ? " class=sel" : "") + (rep ? "" : " style=\\"font-size:12px\\"") + "><td>" +
      (rep ? (n > 1 ? "<button type=button class=gotog data-g=\\"" + t.group + "\\">" + (open_.has(String(t.group)) ? "&minus;" : "+") + (n - 1) + "</button> " : "") + "<b>" + esc(t.description) + "</b>" : "&nbsp;&nbsp;&nbsp;&#8627; " + esc(t.description)) +
      " <small>" + t.term + "</small></td><td class=num>" + t.hits + " / " + t.annotated + "</td><td class=num>" + t.expected + "</td><td class=num>" + (+t.fold).toFixed(2) + "</td><td class=num>" + (+t.p).toExponential(1) +
      "</td><td class=num>" + (+t.p_classic).toExponential(1) + "</td><td>" + esc((t.ids || []).slice(0, 8).map(nameOf).join(", ")) + ((t.ids || []).length > 8 ? ", ..." : "") + "</td></tr>";
    $("pwgoterms").innerHTML = "<table><thead><tr><th>GO term (group)</th><th class=num>Hits / annotated</th><th class=num>Expected</th><th class=num>Fold</th><th class=num>weight01 p</th><th class=num>classic p</th><th>Proteins</th></tr></thead><tbody>" +
      G2.map(g => row(g[0], true, g.length) + (open_.has(String(g[0].group)) ? g.slice(1).map(t => row(t, false, 0)).join("") : "")).join("") + "</tbody></table>";
  }
  window.__pwGoDraw = goDraw;
  $("pwgoterms").addEventListener("click", e => { const b = e.target.closest("button.gotog"); if (b) { const g = b.dataset.g; open_.has(g) ? open_.delete(g) : open_.add(g); goDraw(); return; }
    const tr = e.target.closest("tbody tr"); if (!tr) return; const t = ((goData() || {}).terms || []).find(x => x.term === tr.dataset.term); if (!t) return;
    window.__pwGoSel = window.__pwGoSel && window.__pwGoSel.term === t.term ? null : t; selTerm = null; selCluster = null; goDraw(); drawTerms(); drawNet(); members();
    window.__pwDrill(window.__pwGoSel ? { term: t.term, label: t.description, ids: t.ids, src: "GO topGO, " + (godir === "up" ? "Up" : "Down"), go: true } : null); });
  ["pwgoont", "pwgochart"].forEach(id => $(id).addEventListener("input", () => { open_.clear(); goDraw(); }));
  document.querySelectorAll("#pwgodir button").forEach(b => b.addEventListener("click", () => { godir = b.dataset.dir; document.querySelectorAll("#pwgodir button").forEach(x => x.classList.toggle("on", x === b)); open_.clear(); window.__pwGoSel = null; goDraw(); drawNet(); members(); }));
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
      else if (m === "go") ["up", "down"].forEach(d => ((((c.go || {}).BP || {})[d] || {}).terms || []).filter(t => t.representative).forEach((t, k) => {
          const key = t.term, v = cell[key + "|" + c.name], val = { s: (d === "up" ? 1 : -1) * lfdr(t.p), size: t.hits, t: Object.assign({ fdr: t.p, term_size: t.annotated }, t), d: d };
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
      S.push(m !== "gsea" ? 6 + 2.2 * Math.sqrt(v.size) : 6 + 3 * Math.min(8, v.size));
      H.push(m !== "gsea" ? (v.d === "up" ? "Up" : "Down") + ", " + v.t.hits + " of " + v.t.term_size + ", FDR " + (+v.t.fdr).toPrecision(2) : "NES " + (+v.t.NES).toFixed(2) + ", FDR " + (+v.t.fdr).toPrecision(2)); }));
    const lim = Math.max(...Z.map(Math.abs), 1);
    $("pwoverview").style.height = Math.max(360, 110 + 22 * terms.length) + "px";
    Plotly.react("pwoverview", [{ type: "scatter", mode: "markers", x: X, y: Y, text: H, hovertemplate: "%{y}<br>%{x}: %{text}<extra></extra>",
      marker: { size: S, color: Z, cmin: -lim, cmax: lim, colorscale: [[0, DOWN], [0.5, "#f2f2f2"], [1, UP]], line: { width: 0.5, color: C("--text-primary") },
                colorbar: { title: { text: m === "ora" ? "-log10 FDR (Down -, Up +)" : m === "go" ? "-log10 p weight01 (Down -, Up +)" : "NES" } } } }],
      Object.assign(bp(), { margin: { l: 360, r: 20, t: 10, b: 120 }, xaxis: { type: "category", categoryorder: "array", categoryarray: comps, tickangle: -30, gridcolor: C("--grid") },
        yaxis: { type: "category", categoryorder: "array", categoryarray: terms.map(t => short(pick.get(t), 60)).reverse(), automargin: true, tickfont: { size: 11 }, gridcolor: C("--grid") } }), cfg);
  }
  $("pwoverview").dpdpData = () => { const g = $("pwoverview"), t = (g.data || [])[0]; if (!t) return { cols: [], rows: [] };
    const m = $("pwomethod").value, v = m === "gsea" ? "NES" : m === "go" ? "signed -log10 weight01 p (Up +, Down -)" : "signed -log10 FDR (Up +, Down -)";
    return { cols: ["comparison", "term", v, "details"], rows: Array.from(t.x).map((x, i) => [x, t.y[i], t.marker.color[i], t.text[i]]) }; };
  ["pwomethod", "pwosrc", "pwotopn"].forEach(id => $(id).addEventListener("input", overview));
  // --- term drill-down: volcano and heatmap of the tested proteins of the chosen term (data per comparison, loaded on demand) ---
  const DR = {}; let drSel = null;
  function drLoad(k) { if (DR[k]) return Promise.resolve(DR[k]); const f = W.comps[k].drill_file;
    if (!f) return Promise.reject(new Error("not made: run the Pathway step again"));
    if (!location.protocol.startsWith("http")) return Promise.reject(new Error("open the page through the review server or the web app"));
    return fetch(f).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }).then(d => { d.sidIx = {}; d.prot.sid.forEach((s, j) => { if (s) d.sidIx[s] = j; }); DR[k] = d; return d; }); }
  window.__pwDrill = sel => { drSel = sel; drill(); };
  function drill() {
    const k = ci;
    if (!drSel) { $("pwdrt").textContent = "Term drill-down"; ["pwdrvol", "pwdrheat"].forEach(id => Plotly.purge(id)); $("pwdrtab").innerHTML = "";
      $("pwdrnote").textContent = "Click a term (over-representation, GSEA, GO or GO slim) or a network cluster: its tested proteins are shown on the volcano and in a heatmap."; return; }
    $("pwdrt").textContent = "Term drill-down: " + drSel.label;
    drLoad(k).then(D => { if (k !== ci) return; drDraw(D); }).catch(e => { $("pwdrnote").textContent = "Drill-down data could not be loaded (" + e.message + ")."; });
  }
  function drDraw(D) {
    const P = D.prot, n = P.acc.length, site = D.level === "site";
    const TT = drSel.go ? (D.go_terms || {}) : (D.terms || {});
    let ix = drSel.term && TT[drSel.term], partial = !!drSel.term && (D.partial || []).includes(drSel.term) && !drSel.go;
    if (!ix) { ix = (drSel.ids || []).map(s => D.sidIx[s]).filter(j => j != null); partial = !!drSel.term; }
    ix = [...new Set(ix)];
    const inT = new Set(ix), l2 = f => f == null ? null : Math.sign(f) * Math.log2(Math.max(Math.abs(f), 1)), lp = v => v == null ? null : -Math.log10(Math.max(v, 1e-300));
    const nu = ix.filter(j => P.call[j] === "Up").length, nd = ix.filter(j => P.call[j] === "Down").length;
    const med = (a => { a = a.filter(x => x != null).sort((x, y) => x - y); return a.length ? a[Math.floor(a.length / 2)] : null; })(ix.map(j => l2(P.fc[j])));
    $("pwdrnote").textContent = ix.length + " of " + n + " tested proteins in this term" + (drSel.src ? " (" + drSel.src + ")" : "") + ": " + nu + " Up, " + nd + " Down" +
      (med != null ? "; median log2 FC " + med.toFixed(2) : "") + "." + (partial ? " Only its hits are shown: the STRING download does not list this term (KEGG) or lists it differently (e.g. STRING clusters), so its unchanged members are not known." : "") +
      (site ? " Phospho: each protein is shown by its most significant site." : "");
    // volcano
    const rest = [], mem = { Up: [], Down: [], "": [] };
    for (let j = 0; j < n; j++) { if (P.fc[j] == null || P.p[j] == null) continue; if (inT.has(j)) (mem[P.call[j]] || mem[""]).push(j); else rest.push(j); }
    const labSet = new Set($("pwdrlab").checked ? ix.filter(j => P.p[j] != null).sort((a, b) => P.p[a] - P.p[b]).slice(0, 20) : []), lab = labSet.size > 0;
    const nm = j => (P.name[j] || P.acc[j]) + (site && P.site[j] ? " " + P.site[j] : "");
    const tr = (js, name, col, size, op, text) => ({ type: "scatter", mode: text ? "markers+text" : "markers", name: name, x: js.map(j => l2(P.fc[j])), y: js.map(j => lp(P.p[j])),
      text: text ? js.map(j => labSet.has(j) ? nm(j) : "") : undefined, textposition: "top center", textfont: { size: 10, color: C("--text-secondary") }, hovertext: js.map(nm), customdata: js.map(j => [P.acc[j], P.fc[j], P.p[j]]),
      hovertemplate: "<b>%{hovertext}</b> %{customdata[0]}<br>FC %{customdata[1]:.2f}, p %{customdata[2]:.2g}<extra></extra>", marker: { size: size, color: col, opacity: op, line: { width: text === null ? 0 : 0.5, color: C("--text-primary") } } });
    Plotly.react("pwdrvol", [tr(rest, "other tested", C("--text-muted"), 4, 0.35, null), tr(mem[""], "in term, unchanged", "#7a7a7a", 8, 0.95, lab ? 1 : 0), tr(mem.Down, "in term, Down", DOWN, 9, 0.95, lab ? 1 : 0), tr(mem.Up, "in term, Up", UP, 9, 0.95, lab ? 1 : 0)],
      Object.assign(bp(), { margin: { l: 60, r: 20, t: 10, b: 50 }, legend: { orientation: "h", y: 1.07 }, hovermode: "closest",
        xaxis: { title: "log2 fold change (" + W.comps[ci].name + ")", gridcolor: C("--grid"), zerolinecolor: C("--text-muted") }, yaxis: { title: "-log10 p", gridcolor: C("--grid") } }), cfg);
    const pv = $("pwdrvol"); if (!pv.__click) { pv.__click = true; pv.on("plotly_click", ev => { const pt = ev.points && ev.points[0]; if (pt && pt.customdata && W.has_explore) window.open("10e_explore.html#" + encodeURIComponent(pt.customdata[0]), "_blank"); }); }
    // heatmap
    const rv = $("pwdrrows").value, z = $("pwdrscale").value === "z";
    let H = ix.filter(j => P.p[j] != null).sort((a, b) => P.p[a] - P.p[b]);
    H = rv === "hits" ? H.filter(j => P.call[j] === "Up" || P.call[j] === "Down") : H.slice(0, +rv);
    H.sort((a, b) => (P.fc[b] ?? 0) - (P.fc[a] ?? 0));
    if (!H.length || !D.v || !D.v.length) { Plotly.purge("pwdrheat"); }
    else {
      const row = j => { const a = (D.v[j] || []).map(x => x == null ? null : x / 100); if (!z) return a; const v = a.filter(x => x !== null), m = v.reduce((s, x) => s + x, 0) / Math.max(1, v.length),
        sd = Math.sqrt(v.reduce((s, x) => s + (x - m) ** 2, 0) / Math.max(1, v.length - 1)); return a.map(x => x === null || !sd ? null : (x - m) / sd); };
      const Z = H.map(row), ys = H.map(j => nm(j) + "​".repeat(0) + (H.filter(x => nm(x) === nm(j)).length > 1 ? " (" + P.acc[j] + ")" : ""));
      const xs = D.samples.map((s, c) => s), zmax = z ? 2.5 : null;
      const ann = D.order.map(g => { const c = D.groups.map((x, q) => x === g ? q : -1).filter(q => q >= 0); return { x: (c[0] + c[c.length - 1]) / 2, y: 1.0, xref: "x", yref: "paper", yanchor: "bottom", text: "<b>" + esc(g) + "</b>", showarrow: false, font: { color: DPDP_GCOL(g, D.order) } }; });
      const shapes = D.order.slice(1).map(g => { const c0 = D.groups.indexOf(g); return { type: "line", x0: c0 - 0.5, x1: c0 - 0.5, y0: 0, y1: 1, yref: "paper", line: { color: C("--surface-1"), width: 3 } }; });
      $("pwdrheat").style.height = Math.max(300, 140 + 16 * H.length) + "px";
      Plotly.react("pwdrheat", [{ type: "heatmap", z: Z, x: xs.map((s, c) => c), y: H.map((j, r) => r), text: H.map(j => xs.map((s, c) => nm(j) + "<br>" + s + " (" + D.groups[c] + ")")),
        hovertemplate: "%{text}<br>" + (z ? "z %{z:.2f}" : "log2 %{z:.2f}") + "<extra></extra>", colorscale: z ? [[0, DOWN], [0.5, "#f7f7f7"], [1, UP]] : "Viridis",
        zmin: z ? -zmax : undefined, zmax: z ? zmax : undefined, colorbar: { title: { text: z ? "z-score" : "log2" } }, xgap: 0, ygap: 1 }],
        Object.assign(bp(), { margin: { l: 10, r: 20, t: 34, b: 110 }, annotations: ann, shapes: shapes,
          xaxis: { tickvals: xs.map((s, c) => c), ticktext: xs, tickangle: -45, tickfont: { size: 10 } },
          yaxis: { tickvals: H.map((j, r) => r), ticktext: ys.map((y, r) => (P.call[H[r]] === "Up" ? "▲ " : P.call[H[r]] === "Down" ? "▼ " : "") + y), autorange: "reversed", automargin: true, tickfont: { size: H.length > 60 ? 8 : 10 } } }), cfg);
    }
    // table
    const T = ix.slice().sort((a, b) => (P.p[a] ?? 2) - (P.p[b] ?? 2));
    $("pwdrtt").textContent = "Proteins of the term (" + T.length + ")";
    $("pwdrtab").innerHTML = T.length ? "<table><thead><tr><th>Protein</th>" + (site ? "<th>Site</th>" : "") + "<th>Accession</th><th class=num>FC</th><th class=num>p</th><th>Call</th></tr></thead><tbody>" +
      T.slice(0, 500).map(j => "<tr><td>" + esc(P.name[j] || "") + "</td>" + (site ? "<td>" + esc(P.site[j] || "") + "</td>" : "") + "<td>" + (W.has_explore ? "<a href=\\"10e_explore.html#" + encodeURIComponent(P.acc[j]) + "\\" target=_blank>" + esc(P.acc[j]) + "</a>" : esc(P.acc[j])) +
        "</td><td class=num>" + (P.fc[j] != null ? (+P.fc[j]).toFixed(2) : "") + "</td><td class=num>" + (P.p[j] != null ? (+P.p[j]).toPrecision(2) : "") + "</td><td" + (P.call[j] === "Up" ? " class=xup" : P.call[j] === "Down" ? " class=xdown" : "") + ">" + esc(P.call[j] || "") + "</td></tr>").join("") +
      "</tbody></table>" + (T.length > 500 ? "<p class=note>First 500 of " + T.length + " (smallest p).</p>" : "") : "<p class=note>None of the tested proteins.</p>";
  }
  ["pwdrlab", "pwdrscale", "pwdrrows"].forEach(id => $(id).addEventListener("input", drill));
  // --- GO slim profile ---
  function slimDraw() {
    const S = W.comps[ci].slim, o = $("pwslont").value, d = S && S[o];
    if (!S) { $("pwslnote").textContent = "The GO slim profile was not made (switch it on in the options and run again)."; Plotly.purge("pwslchart"); $("pwslterms").innerHTML = ""; return; }
    $("pwslnote").textContent = W.comps[ci].name + ": " + d.annotated + " tested proteins with a GO " + { BP: "process", MF: "function", CC: "component" }[o] + " annotation (" + d.up + " Up, " + d.down + " Down).";
    let T = (d.terms || []).slice(); if (!T.length) { Plotly.purge("pwslchart"); $("pwslterms").innerHTML = "<p class=note>No category.</p>"; return; }
    const best = t => Math.min(t.fdr_up ?? 1, t.fdr_down ?? 1);
    if ($("pwslorder").value === "sig") T.sort((a, b) => best(a) - best(b) || b.tested - a.tested);
    T = T.slice(0, +$("pwsltopn").value);
    const star = (t, w) => (t["fdr_" + w] != null && t["fdr_" + w] <= 0.05) ? " *" : "";
    const ys = T.map(t => short(t.description, 50)).reverse(), R = T.slice().reverse();
    Plotly.react("pwslchart", [["tested", "pct_tested", C("--text-muted"), "tested"], ["Up", "pct_up", UP, "up"], ["Down", "pct_down", DOWN, "down"]].map(([nmv, f, col, w]) => ({
      type: "bar", orientation: "h", name: nmv, y: ys, x: R.map(t => t[f]), marker: { color: col }, text: R.map(t => w === "tested" ? "" : star(t, w)), textposition: "outside", cliponaxis: false, constraintext: "none", textfont: { size: 16, color: col },
      customdata: R.map(t => [t.term, t[w], w === "tested" ? null : t["fdr_" + w], t.description]), meta: { cd: ["GO id", "proteins", "FDR", "description"] },
      hovertemplate: "%{y}<br>" + nmv + ": %{customdata[1]} proteins (%{x:.1f}%)" + (w === "tested" ? "" : ", FDR %{customdata[2]:.2g}") + "<extra>%{customdata[0]}</extra>" })),
      Object.assign(bp(), { barmode: "group", bargap: 0.25, margin: { l: 330, r: 30, t: 10, b: 50 }, legend: { orientation: "h", y: 1.04 },
        xaxis: { title: "% of the annotated proteins of each set", gridcolor: C("--grid") }, yaxis: { categoryorder: "array", categoryarray: ys, automargin: true, tickfont: { size: 11 } } }), cfg);
    $("pwslchart").style.height = Math.max(360, 90 + 26 * T.length) + "px"; Plotly.Plots.resize("pwslchart");
    const sc = $("pwslchart"); if (!sc.__click) { sc.__click = true; sc.on("plotly_click", ev => { const pt = ev.points && ev.points[0]; if (!pt) return; const t = (W.comps[ci].slim[$("pwslont").value].terms || []).find(x => x.term === pt.customdata[0]); if (t) slimPick(t); }); }
    const fdr = v => v == null ? "" : (+v).toPrecision(2);
    $("pwslterms").innerHTML = "<table><thead><tr><th>GO slim category</th><th class=num>Tested</th><th class=num>%</th><th class=num>Up</th><th class=num>%</th><th class=num>FDR Up</th><th class=num>Down</th><th class=num>%</th><th class=num>FDR Down</th></tr></thead><tbody>" +
      T.map(t => "<tr data-term=\\"" + t.term + "\\" style=\\"cursor:pointer\\"><td>" + esc(t.description) + " <small>" + t.term + "</small></td><td class=num>" + t.tested + "</td><td class=num>" + (+t.pct_tested).toFixed(1) +
        "</td><td class=num>" + t.up + "</td><td class=num>" + (+t.pct_up).toFixed(1) + "</td><td class=num" + (star(t, "up") ? " style=\\"color:" + UP + ";font-weight:600\\"" : "") + ">" + fdr(t.fdr_up) +
        "</td><td class=num>" + t.down + "</td><td class=num>" + (+t.pct_down).toFixed(1) + "</td><td class=num" + (star(t, "down") ? " style=\\"color:" + DOWN + ";font-weight:600\\"" : "") + ">" + fdr(t.fdr_down) + "</td></tr>").join("") + "</tbody></table>";
  }
  function slimPick(t) { window.__pwDrill({ term: t.term, label: t.description + " (GO slim)", src: "GO slim", go: true }); $("pwdrill").scrollIntoView({ behavior: "smooth", block: "start" }); }
  $("pwslterms").addEventListener("click", e => { const tr = e.target.closest("tbody tr"); if (!tr) return; const t = (W.comps[ci].slim[$("pwslont").value].terms || []).find(x => x.term === tr.dataset.term); if (t) slimPick(t); });
  ["pwslont", "pwslorder", "pwsltopn"].forEach(id => $(id).addEventListener("input", slimDraw));
  ["pwchart", "pwtopn"].forEach(id => $(id).addEventListener("input", oraChart));
  $("pwsrc").addEventListener("input", oraChart); $("pwq").addEventListener("input", () => {});
  $("pwcomp").addEventListener("change", () => { gsel = null; window.__pwGoSel = null; open_.clear(); gseaSources(); gseaDraw(); oraChart(); goDraw(); slimDraw(); window.__pwDrill(null); });
  document.querySelectorAll(".pwdir:not(#pwgodir) button").forEach(b => b.addEventListener("click", oraChart));
  // --- downloads: every table (full content, not the shortened lists on screen) and the charts whose points need more columns ---
  const cname = () => W.comps[ci].name, D_ = () => dir === "up" ? "Up" : "Down", names = ids => (ids || []).map(nameOf);
  const accOf = id => { const h = (W.comps[ci].hits || {})[id]; return h ? h.acc : ""; };
  if (window.dpdpCsvButton) {
    dpdpCsvButton("pwterms", () => "Enriched terms " + cname() + " " + D_(), () => { const src = $("pwsrc").value, q = $("pwq").value.toLowerCase();
      return { cols: ["comparison", "direction", "source", "term", "description", "hits", "term size (tested)", "fold enrichment", "p", "FDR", "proteins", "accessions"],
        rows: terms().filter(t => (!src || t.source === src) && (!q || (t.description + " " + t.term).toLowerCase().includes(q)))
          .map(t => [cname(), D_(), t.source, t.term, t.description, t.hits, t.term_size, t.fold, t.p, t.fdr, names(t.ids), (t.ids || []).map(accOf)]) }; });
    dpdpCsvButton("pwgterms", () => "GSEA " + cname(), () => ({ cols: ["comparison", "source", "term", "description", "size", "NES", "p", "FDR", "leading edge"],
      rows: gseaRows().map(t => [cname(), t.source, t.term, t.description, t.size, t.NES, t.p, t.fdr, names(t.leading)]) }));
    dpdpCsvButton("pwgoterms", () => "GO topGO " + cname() + " " + $("pwgoont").value + " " + (godir === "up" ? "Up" : "Down"), () => {
      const d = goData(); return { cols: ["comparison", "ontology", "direction", "group", "representative", "term", "description", "annotated", "hits", "expected", "fold", "weight01 p", "classic p", "proteins"],
        rows: ((d && d.terms) || []).map(t => [cname(), $("pwgoont").value, godir === "up" ? "Up" : "Down", t.group, t.representative ? "yes" : "", t.term, t.description, t.annotated, t.hits, t.expected, t.fold, t.p, t.p_classic, names(t.ids)]) }; });
    dpdpCsvButton("pwslterms", () => "GO slim " + cname() + " " + $("pwslont").value, () => { const S = W.comps[ci].slim, d = S && S[$("pwslont").value];
      return { cols: ["comparison", "ontology", "term", "category", "tested", "% tested", "Up", "% Up", "FDR Up", "p Up", "Down", "% Down", "FDR Down", "p Down"],
        rows: ((d && d.terms) || []).map(t => [cname(), $("pwslont").value, t.term, t.description, t.tested, t.pct_tested, t.up, t.pct_up, t.fdr_up, t.p_up, t.down, t.pct_down, t.fdr_down, t.p_down]) }; });
    dpdpCsvButton("pwclusters", () => "Network clusters " + cname(), () => { const net = W.comps[ci].network || {};
      return { cols: ["comparison", "cluster", "proteins", "enriched for (top 3)", "members"], rows: (net.clusters || []).map(c => [cname(), c.cluster, c.size, c.terms || [], (net.nodes || []).filter(n => n.cluster === c.cluster).map(n => n.name)]) }; });
    dpdpCsvButton("pwmembers", () => $("pwmemt").textContent + " " + cname(), () => { const H = highlighted(), net = W.comps[ci].network;
      const byId = Object.assign({}, W.comps[ci].hits || {}); ((net && net.nodes) || []).forEach(n => { byId[n.id] = Object.assign({}, byId[n.id], n); });
      return { cols: ["protein", "accession", "STRING id", "FC", "p", "call"], rows: H ? [...H].map(id => { const n = byId[id] || {}; return [n.name || nameOf(id), n.acc || "", id, n.fc, n.p, n.call || ""]; }) : [] }; });
    dpdpCsvButton("pwdrtab", () => "Drill-down " + (drSel ? drSel.label : "") + " " + cname(), () => drRows());
  }
  // drill-down rows: every tested protein of the term (also used by the volcano download)
  function drRows(all) {
    const D = DR[ci]; if (!D || !drSel) return { cols: [], rows: [] };
    const P = D.prot, TT = drSel.go ? (D.go_terms || {}) : (D.terms || {}); let ix = drSel.term && TT[drSel.term];
    if (!ix) ix = (drSel.ids || []).map(x => D.sidIx[x]).filter(j => j != null);
    const inT = new Set(ix), J = all ? P.acc.map((_, j) => j) : [...inT].sort((a, b) => (P.p[a] ?? 2) - (P.p[b] ?? 2)), site = D.level === "site";
    return { cols: ["comparison", "term", "in term", "protein"].concat(site ? ["site"] : []).concat(["accession", "FC", "log2 FC", "p", "-log10 p", "call"]),
      rows: J.map(j => [cname(), drSel.label, inT.has(j) ? "yes" : "", P.name[j] || ""].concat(site ? [P.site[j] || ""] : []).concat([P.acc[j], P.fc[j],
        P.fc[j] == null ? "" : +(Math.sign(P.fc[j]) * Math.log2(Math.max(Math.abs(P.fc[j]), 1))).toFixed(4), P.p[j], P.p[j] == null ? "" : +(-Math.log10(Math.max(P.p[j], 1e-300))).toFixed(4), P.call[j] || ""])) };
  }
  $("pwdrvol").dpdpData = () => drRows(true);
  $("pwnet").dpdpData = () => { const net = W.comps[ci].network || {}, H = highlighted();
    return { cols: ["comparison", "protein", "accession", "gene", "STRING id", "FC", "p", "call", "cluster", "links", "highlighted", "x", "y"],
      rows: (net.nodes || []).map(n => [cname(), n.name, n.acc, n.gene, n.id, n.fc, n.p, n.call, n.cluster, n.degree, H && H.has(n.id) ? "yes" : "", n.x, n.y]) }; };
  fillSources(); drawTerms(); drawNet(); members(); oraChart(); gseaSources(); gseaDraw(); goDraw(); slimDraw(); ovSources(); overview();
})();
'
