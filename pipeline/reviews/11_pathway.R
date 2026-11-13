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
    '<p class="note" id="pwnote"></p><div id="pwterms" style="max-height:460px;overflow:auto"></div></section>',
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
  fillSources(); drawTerms(); drawNet(); members();
})();
'
