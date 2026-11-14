# Review page for Step 10e - Explore: search any protein (or peptide), plot its value per sample and the peptides of its protein
# (mean +- SD per group, log2 intensity or z-score), for a comparison's samples (with the pool) or all samples. The data
# (review/10e_explore_data.json.gz, written by steps/10e_explore.R) is loaded when the page opens; everything else runs in the browser.
# Link to a protein: 10e_explore.html#<accession>.

review_explore <- function(db_path) {
  run_dir <- dirname(db_path)
  p <- get_params(db_path)
  f <- file.path(run_dir, "steps", "10_stats", "explore.json")
  info <- if (file.exists(f)) fromJSON(f, simplifyVector = FALSE) else list()
  pep <- identical(info$level, "peptide")
  U <- if (pep) "peptide" else "protein"
  # The Stats Excel workbooks (step 10c) and each comparison's sheet in them
  outj <- file.path(run_dir, "steps", "10_stats", "output.json")
  out <- if (file.exists(outj)) fromJSON(outj, simplifyVector = FALSE) else list()
  xl <- unlist(out$excel$files)
  if (!length(xl)) xl <- list.files(file.path(run_dir, info$norm %||% p$stat_norm %||% ""), pattern = "\\.xlsx$", full.names = TRUE)
  xl <- xl[file.exists(xl)]
  sheets <- unlist(out$excel$sheets); sheet_of <- list()
  for (cn in unlist(info$comparisons)) { m <- sheets[sub("^Table[0-9]+ ", "", sheets) == cn]; if (length(m)) sheet_of[[cn]] <- m[1] }
  xl_html <- if (length(xl)) paste0('<p class="note">Excel: ', paste(vapply(xl, function(f) { u <- run_file_link(run_dir, f)
      if (nzchar(u)) sprintf('<a href="%s" download>download %s</a> (%s)', esc(u), esc(basename(f)), { b <- file.size(f); if (b >= 1e6) sprintf("%.1f MB", b / 1e6) else sprintf("%d KB", max(1L, as.integer(round(b / 1e3)))) }) else sprintf("<code>%s</code>", esc(f)) }, ""), collapse = " &middot; "),
      '<span id="xtsheet"></span></p>') else ""
  body <- paste0(
    '<section><h2>Find a ', U, '</h2>',
    sprintf('<p class="note">Search by gene, accession or description%s. %s of the tested data (%s normalization, %s %ss); FC and p are from the chosen comparison. Charts download from their toolbar (PNG / SVG).</p>',
            if (pep) ", sequence or site" else "", if (pep) "Every peptide" else "Every protein", esc(info$norm %||% p$stat_norm),
            fmt_n(info$rows %||% NA), U),
    '<div class="control"><input type="search" id="xq" placeholder="e.g. ACTB, P60709, actin" autocomplete="off" style="font:inherit;font-size:15px;padding:4px 10px;width:26em"/>',
    '<label for="xcomp">Samples</label><select id="xcomp"></select>',
    '<label><input type="checkbox" id="xpool" checked/> with the pool</label></div>',
    '<div id="xhits"></div><div id="xrecent" class="xchips"></div><p class="dmsg" id="xmsg">Loading the data...</p></section>',
    '<section id="xsel" hidden><h2 id="xtitle"></h2><p class="note" id="xinfo"></p>',
    '<div><h3 id="xpt">', if (pep) "Peptide" else "Protein", ' value per sample</h3>',
    '<div class="control"><label for="xpscale">Scale</label><select id="xpscale"><option value="int">intensity</option><option value="log2">log2 intensity</option></select></div>',
    '<div id="xprot" style="height:380px" role="img" aria-label="Value per sample"></div></div>',
    '<div><h3 id="xpept">Peptides of the protein</h3>',
    '<div class="control"><label for="xscale">Scale</label><select id="xscale"><option value="log2">log2 intensity (mean &plusmn; SD)</option><option value="z">z-score (mean &plusmn; SD)</option></select>',
    '<label for="xorder">Order</label><select id="xorder"><option value="pos">position in the protein</option><option value="int">intensity</option></select></div>',
    '<div id="xpep" style="height:480px" role="img" aria-label="Peptides of the protein, mean per group"></div><p class="note" id="xpepnote"></p></div>',
    '<h3>In each comparison</h3><div id="xstats"></div></section>',
    '<section id="xtsec"><h2 id="xttitle">Comparison table</h2><p class="note">Every row and column of the comparison&apos;s Excel sheet (the app&apos;s final stats table, without the _FC2 columns). Click a row to plot it; click a column name to sort.</p>', xl_html,
    '<div class="control"><input type="search" id="xtq" placeholder="filter rows (any column)" style="font:inherit;font-size:14px;padding:3px 8px;width:18em"/>',
    '<label><input type="checkbox" id="xtsig"/> Up/Down only</label><button type="button" id="xtcsv">Download CSV (these rows)</button><span class="note" id="xtcount"></span></div>',
    '<div id="xtwrap" style="overflow:auto;max-height:640px;border:1px solid var(--grid)"><table id="xtab" class="xtab"></table></div>',
    '<div class="control"><button type="button" id="xtprev">&lsaquo; previous</button><span id="xtpage"></span><button type="button" id="xtnext">next &rsaquo;</button></div></section>',
    '<style>#xhits button { display:block; width:100%; text-align:left; font:inherit; font-size:14px; padding:4px 8px; border:0; border-bottom:1px solid var(--grid); background:transparent; color:var(--text-primary); cursor:pointer; }
#xhits button:hover, #xhits button.on { background:var(--surface-0); } #xhits small { color:var(--text-secondary); } .xup { color:#c0392b; font-weight:600; } .xdown { color:#2e6fb7; font-weight:600; }
.xtab { border-collapse:collapse; font-size:12px; white-space:nowrap; } .xtab th { position:sticky; top:0; background:var(--surface-1); cursor:pointer; padding:4px 8px; border-bottom:2px solid var(--rule); text-align:left; }
.xtab td { padding:2px 8px; border-bottom:1px solid var(--grid); } .xtab td.n { text-align:right; font-variant-numeric:tabular-nums; } .xtab tbody tr { cursor:pointer; } .xtab tbody tr:hover { background:var(--surface-0); }
.xtab tr.sel td { background:var(--warn-bg); font-weight:600; } .xtab td.up { color:#c0392b; font-weight:600; } .xtab td.down { color:#2e6fb7; font-weight:600; }
.xchips button { font:inherit; font-size:12px; margin:6px 4px 0 0; padding:1px 8px; border:1px solid var(--rule); border-radius:10px; background:transparent; color:var(--text-primary); cursor:pointer; }</style>',
    '<script>', plotly_js(), '</script><script>window.XSHEETS = ', toJSON(if (length(sheet_of)) sheet_of else setNames(list(), character()), auto_unbox = TRUE), ';', EXPLORE_JS, '</script>')
  write_review_page(run_dir, "10e_explore", "Explore - protein and peptide plots", "Search any protein", "done", body)
}

EXPLORE_JS <- '
(function () {
  const $ = id => document.getElementById(id), msg = $("xmsg");
  const css = getComputedStyle(document.querySelector(".viz-root")), C = n => css.getPropertyValue(n).trim();
  const esc = s => String(s == null ? "" : s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\\"": "&quot;" }[c]));
  if (!location.protocol.startsWith("http")) { msg.textContent = "Open this page through the review server or the DPDP web app: it loads its data from the server."; return; }
  fetch("10e_explore_data.json.gz").then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }).then(init)
    .catch(e => { msg.textContent = "Could not load the data (" + e.message + "). Ask Claude to run Step 10e (Explore)."; });
  function init(D) {
    const R = D.rows, n = R.acc.length, pep = D.level === "peptide", P = D.pep;
    const lower = s => String(s || "").toLowerCase();
    const text = R.acc.map((a, i) => lower([R.gene[i], a, R.desc[i], pep ? R.seq[i] : "", pep ? R.site[i] : ""].join(" ")));
    const gnL = R.gene.map(lower), acL = R.acc.map(lower);
    // peptides per protein row (Protein output) / rows per accession (Peptide output)
    const byProt = {}; if (P) P.prot.forEach((k, j) => { (byProt[k] = byProt[k] || []).push(j); });
    const byAcc = {}; if (pep) R.acc.forEach((a, i) => { (byAcc[a] = byAcc[a] || []).push(i); });
    const comp = $("xcomp");
    D.comps.forEach((c, k) => comp.appendChild(Object.assign(document.createElement("option"), { value: k, textContent: c.name })));
    comp.appendChild(Object.assign(document.createElement("option"), { value: "all", textContent: "all samples" }));
    const cur = () => comp.value === "all" ? null : D.comps[+comp.value];
    const sampleSet = () => { const c = cur(); if (!c) return D.samples.map((_, i) => i); return c.N.concat(c.D, $("xpool").checked ? c.P : []); };
    const callTag = (c, i) => { if (!c || c.call[i] == null) return ""; const t = c.call[i]; return t === "Up" ? " <span class=xup>Up</span>" : t === "Down" ? " <span class=xdown>Down</span>" : ""; };
    const fmtP = v => v == null ? "-" : (+v).toPrecision(2), fmtF = v => v == null ? "-" : (+v).toFixed(2);
    const rowName = i => (R.gene[i] || R.acc[i]) + (pep ? " " + (R.site[i] || R.seq[i]) : "");
    let sel = -1, recent = [];
    // --- search: gene exact, gene starts, accession starts, anything contains; the hits of the chosen comparison first on ties
    function search() {
      const q = lower($("xq").value.trim()), box = $("xhits"); box.innerHTML = "";
      if (q.length < 2) return;
      const c = cur(), hits = [];
      for (let i = 0; i < n && hits.length < 4000; i++) {
        let s = gnL[i] === q ? 0 : gnL[i].startsWith(q) ? 1 : acL[i].startsWith(q) ? 2 : text[i].includes(q) ? 3 : -1;
        if (s >= 0) hits.push([s, c && c.p[i] != null ? c.p[i] : 2, i]);
      }
      hits.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
      hits.slice(0, 15).forEach(([s, pv, i]) => { const b = document.createElement("button");
        b.innerHTML = "<b>" + esc(rowName(i)) + "</b> " + esc(R.acc[i]) + " <small>" + esc(R.desc[i]) + (c ? " &middot; FC " + fmtF(c.fc[i]) + ", p " + fmtP(c.p[i]) : "") + "</small>" + callTag(c, i);
        b.onclick = () => choose(i); box.appendChild(b); });
      if (hits.length > 15) box.appendChild(Object.assign(document.createElement("p"), { className: "note", textContent: (hits.length >= 4000 ? "4,000+" : hits.length) + " matches; type more to narrow it down." }));
      if (!hits.length) box.innerHTML = "<p class=note>No match.</p>";
    }
    let t0; $("xq").addEventListener("input", () => { clearTimeout(t0); t0 = setTimeout(search, 120); });
    $("xq").addEventListener("keydown", e => { if (e.key === "Enter") { const b = $("xhits").querySelector("button"); if (b) b.click(); } });
    function choose(i) { sel = i; recent = [i].concat(recent.filter(x => x !== i)).slice(0, 10); $("xhits").innerHTML = ""; history.replaceState(null, "", "#" + encodeURIComponent(R.acc[i] + (pep ? "|" + R.seq[i] : ""))); draw(); }
    const base = () => ({ paper_bgcolor: C("--surface-1"), plot_bgcolor: C("--surface-1"), font: { color: C("--text-secondary"), size: 12 }, legend: { orientation: "h", y: 1.12 } });
    const ax = t => ({ title: t, gridcolor: C("--grid"), zeroline: false });
    const cfg = { displaylogo: false, responsive: true };
    const GC = g => DPDP_GCOL(g, D.order);
    function draw() {
      if (sel < 0) return;
      const i = sel, c = cur(), S = sampleSet(), smp = S.map(k => D.samples[k]), grs = S.map(k => D.groups[k]);
      $("xsel").hidden = false;
      $("xtitle").textContent = rowName(i) + "  " + R.acc[i];
      $("xinfo").innerHTML = esc(R.desc[i]) + (pep && R.seq[i] ? " &middot; " + esc(R.seq[i]) : "") + (c ? " &middot; <b>" + esc(c.name) + "</b>: FC " + fmtF(c.fc[i]) + ", p " + fmtP(c.p[i]) + ", adj p " + fmtP(c.adj[i]) + callTag(c, i) + (c.p[i] == null ? " (not tested: missing factor)" : "") : "");
      // value per sample
      const lin = $("xpscale").value === "int", vals = S.map(k => R.v[i][k] == null ? null : (lin ? Math.pow(2, R.v[i][k] / 100) : R.v[i][k] / 100));
      const groupsInOrder = D.order.filter(g => grs.includes(g));
      Plotly.react("xprot", groupsInOrder.map(g => { const ix = grs.map((x, j) => x === g ? j : -1).filter(j => j >= 0);
        return { type: "bar", name: g, x: ix.map(j => smp[j]), y: ix.map(j => vals[j]), marker: { color: GC(g) }, hovertemplate: "%{x}<br>%{y:.4s}<extra>" + g + "</extra>" }; }),
        Object.assign(base(), { margin: { l: 64, r: 10, t: 10, b: 90 }, bargap: 0.25, xaxis: { categoryorder: "array", categoryarray: smp, tickangle: -45 },
          yaxis: Object.assign(ax(lin ? "intensity" : "log2 intensity"), lin ? { rangemode: "tozero" } : {}) }), cfg);
      // peptides of the protein
      const rows = pep ? (byAcc[R.acc[i]] || []) : (byProt[i] || []), V = pep ? R.v : (P ? P.v : []);
      const seqOf = j => pep ? R.seq[j] : P.seq[j], posOf = j => pep ? R.pos[j] : P.pos[j];
      $("xpept").textContent = pep ? "Peptides of " + (R.gene[i] || R.acc[i]) : "Peptides of the protein";
      if (!rows.length) { Plotly.purge("xpep"); $("xpepnote").textContent = "No peptide data for this protein."; return; }
      const z = $("xscale").value === "z";
      const lv = j => S.map(k => V[j][k] == null ? null : V[j][k] / 100);
      const zrow = a => { const v = a.filter(x => x !== null), m = v.reduce((s, x) => s + x, 0) / v.length, sd = Math.sqrt(v.reduce((s, x) => s + (x - m) ** 2, 0) / Math.max(1, v.length - 1)); return a.map(x => x === null || !sd ? null : (x - m) / sd); };
      const meanI = j => { const v = lv(j).filter(x => x !== null); return v.length ? v.reduce((s, x) => s + x, 0) / v.length : -99; };
      let ord = rows.slice();
      if ($("xorder").value === "int") ord.sort((a, b) => meanI(b) - meanI(a)); else ord.sort((a, b) => (posOf(a) ?? 1e9) - (posOf(b) ?? 1e9));
      const MAX = 60, shown = ord.length > MAX ? (pep ? ord.slice().sort((a, b) => Math.abs((posOf(a) ?? 0) - (posOf(i) ?? 0)) - Math.abs((posOf(b) ?? 0) - (posOf(i) ?? 0))).slice(0, MAX).sort((a, b) => ord.indexOf(a) - ord.indexOf(b)) : ord.slice(0, MAX)) : ord;
      const lab = j => (pep && j === i ? "<b>&#9654; " : "") + (pep && R.site[j] ? R.site[j] + " " : "") + seqOf(j).replace(/\\[[^\\]]*\\]/g, "*").slice(0, 26) + (posOf(j) != null ? " (" + posOf(j) + ")" : "") + (pep && j === i ? "</b>" : "");
      const xs = shown.map(lab), data = shown.map(j => z ? zrow(lv(j)) : lv(j));
      const stat = (row, ix) => { const v = ix.map(k => row[k]).filter(x => x !== null); if (!v.length) return [null, null]; const m = v.reduce((s, x) => s + x, 0) / v.length; return [m, Math.sqrt(v.reduce((s, x) => s + (x - m) ** 2, 0) / Math.max(1, v.length - 1))]; };
      const shapes = z ? [0, 1, -1].map(yv => ({ type: "line", xref: "paper", x0: 0, x1: 1, y0: yv, y1: yv, line: { dash: "dot", width: 1, color: C("--text-secondary") } })) : [];
      Plotly.react("xpep", groupsInOrder.map(g => { const ix = grs.map((x, j) => x === g ? j : -1).filter(j => j >= 0);
        return { type: "bar", name: g, x: xs, y: data.map(r => stat(r, ix)[0]), error_y: { type: "data", array: data.map(r => stat(r, ix)[1]), color: C("--text-secondary"), thickness: 1 },
          marker: { color: GC(g) }, hovertemplate: "%{x}<br>mean %{y:.2f}<extra>" + g + "</extra>" }; }),
        Object.assign(base(), { barmode: "group", shapes, margin: { l: 64, r: 10, t: 10, b: 140 }, xaxis: { tickangle: -60, tickfont: { size: 9 } },
          yaxis: ax(z ? "z-score (mean ± SD)" : "log2 intensity (mean ± SD)") }), cfg);
      $("xpepnote").textContent = rows.length + " peptide" + (rows.length === 1 ? "" : "s") + (rows.length > MAX ? "; showing " + MAX + (pep ? " nearest to the selected one" : ", by the order chosen") : "") +
        (z ? ". z-score: each peptide centred on its mean over the shown samples, in SDs (as the app; on log2)." : ".");
      // every comparison
      $("xstats").innerHTML = "<table><thead><tr><th>Comparison</th><th class=num>FC</th><th class=num>p</th><th class=num>adj p</th><th>Call</th></tr></thead><tbody>" +
        D.comps.map(cc => "<tr><td>" + esc(cc.name) + "</td><td class=num>" + fmtF(cc.fc[i]) + "</td><td class=num>" + fmtP(cc.p[i]) + "</td><td class=num>" + fmtP(cc.adj[i]) + "</td><td>" + (cc.p[i] == null ? "not tested" : (cc.call[i] || "-")) + "</td></tr>").join("") + "</tbody></table>";
      if (window.dpdpExploreDrawn) window.dpdpExploreDrawn();
      $("xrecent").innerHTML = ""; if (recent.length > 1) recent.forEach(k => { const b = document.createElement("button"); b.textContent = rowName(k); b.onclick = () => choose(k); $("xrecent").appendChild(b); });
    }
    // --- the comparison table (its Excel sheet), loaded per comparison ---
    const T = { cache: {}, cur: null, rows: [], sortCol: -1, sortDir: 1, page: 0, per: 50 };
    const fmtCell = (name, v) => { if (v === null || v === undefined) return "";
      if (typeof v !== "number") { const x = String(v); return x.length > 60 && (name === "Description" || name === "Name") ? "<span title=" + JSON.stringify(esc(x)) + ">" + esc(x.slice(0, 57)) + "...</span>" : esc(x); }
      if (/_pval$|_adjpval$|_limma$/.test(name)) return v.toExponential(2);
      if (/_FC$|_CV$|_mf$|_cohensd$/.test(name)) return v.toFixed(2);
      return Math.abs(v) >= 1000 ? Math.round(v).toLocaleString("en-US") : String(v); };
    function loadTable() {
      const c = cur(); if (!c || !c.table) { $("xtsec").hidden = true; return; }
      $("xtsec").hidden = false; $("xttitle").textContent = "Comparison table: " + c.name;
      if ($("xtsheet")) $("xtsheet").textContent = (window.XSHEETS || {})[c.name] ? ". This comparison is the sheet \\u201c" + XSHEETS[c.name] + "\\u201d." : "";
      if (T.cache[c.table]) { T.cur = T.cache[c.table]; filterTable(true); return; }
      $("xtcount").textContent = "loading...";
      fetch(c.table).then(r => r.json()).then(t => { T.cache[c.table] = t; T.cur = t; T.sortCol = -1; filterTable(true); })
        .catch(() => { $("xtcount").textContent = "Could not load the table (re-run Step 10e)."; });
    }
    function filterTable(jump) {
      const t = T.cur; if (!t) return; const q = lower($("xtq").value.trim()), sig = $("xtsig").checked, nrow = t.row.length;
      const si = t.columns.indexOf("Stats"); let rows = [];
      for (let r = 0; r < nrow; r++) {
        if (sig && si >= 0 && !t.data[si][r]) continue;
        if (q && !t.data.some((col, j) => !t.numeric[j] && col[r] != null && lower(col[r]).includes(q))) continue;
        rows.push(r); }
      if (T.sortCol >= 0) { const col = t.data[T.sortCol], num = t.numeric[T.sortCol];
        rows.sort((a, b) => { const x = col[a], y = col[b]; if (x == null) return 1; if (y == null) return -1; return (num ? x - y : String(x).localeCompare(String(y))) * T.sortDir; }); }
      T.rows = rows; T.page = 0; renderTable(jump === true);
    }
    function renderTable(jump) {
      const t = T.cur; if (!t) return;
      if (jump && sel >= 0) { const at = T.rows.findIndex(r => t.row[r] === sel); if (at >= 0) T.page = Math.floor(at / T.per); }
      const pages = Math.max(1, Math.ceil(T.rows.length / T.per)); T.page = Math.min(T.page, pages - 1);
      const part = T.rows.slice(T.page * T.per, (T.page + 1) * T.per);
      $("xtab").innerHTML = "<thead><tr>" + t.columns.map((c, j) => "<th data-j=" + j + ">" + esc(c) + (T.sortCol === j ? (T.sortDir > 0 ? " &#9650;" : " &#9660;") : "") + "</th>").join("") + "</tr></thead><tbody>" +
        part.map(r => "<tr data-r=" + r + (t.row[r] === sel ? " class=sel" : "") + ">" + t.columns.map((c, j) => { const v = t.data[j][r];
          return "<td class=" + JSON.stringify((t.numeric[j] ? "n" : "") + (c === "Stats" && v ? " " + lower(v) : "")) + ">" + fmtCell(c, v) + "</td>"; }).join("") + "</tr>").join("") + "</tbody>";
      $("xtcount").textContent = T.rows.length.toLocaleString("en-US") + " of " + t.row.length.toLocaleString("en-US") + " rows";
      $("xtpage").textContent = " page " + (T.page + 1) + " of " + pages + " ";
      $("xtprev").disabled = T.page === 0; $("xtnext").disabled = T.page >= pages - 1;
    }
    $("xtab").addEventListener("click", e => { const th = e.target.closest("th"), tr = e.target.closest("tbody tr");
      if (th) { const j = +th.dataset.j; T.sortDir = T.sortCol === j ? -T.sortDir : (T.cur.numeric[j] ? 1 : 1); T.sortCol = j; filterTable(); return; }
      if (tr) { const k = T.cur.row[+tr.dataset.r]; if (k >= 0) { choose(k); document.getElementById("xsel").scrollIntoView({ behavior: "smooth" }); } } });
    let t1; $("xtq").addEventListener("input", () => { clearTimeout(t1); t1 = setTimeout(() => filterTable(), 200); });
    $("xtsig").addEventListener("change", () => filterTable());
    $("xtprev").onclick = () => { T.page--; renderTable(false); }; $("xtnext").onclick = () => { T.page++; renderTable(false); };
    $("xtcsv").onclick = () => { const t = T.cur; if (!t) return; const DQ = String.fromCharCode(34), NL = String.fromCharCode(10);
      const q = v => { if (v == null) return ""; const x = String(v); return (x.includes(DQ) || x.includes(",") || x.includes(NL)) ? DQ + x.split(DQ).join(DQ + DQ) + DQ : x; };
      const csv = [t.columns.map(q).join(",")].concat(T.rows.map(r => t.columns.map((_, k) => q(t.data[k][r])).join(","))).join(NL);
      const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" })); a.download = D.run + "_" + t.name + "_table.csv"; document.body.appendChild(a); a.click(); a.remove(); };
    window.dpdpExploreDrawn = () => { if (T.cur) renderTable(true); };
    ["xcomp", "xpool", "xscale", "xpscale", "xorder"].forEach(id => $(id).addEventListener("change", () => { if (id === "xcomp") { $("xpool").disabled = comp.value === "all"; loadTable(); } draw(); search(); }));
    loadTable();
    msg.textContent = "";
    // start: the protein in the URL (#accession), else the top hit of the first comparison
    const h = decodeURIComponent(location.hash.slice(1));
    let start = h ? R.acc.findIndex((a, k) => a === h.split("|")[0] && (!pep || !h.includes("|") || R.seq[k] === h.split("|")[1])) : -1;
    if (start < 0 && D.comps.length) { const c = D.comps[0]; let best = 2; c.p.forEach((pv, k) => { if (pv != null && c.call[k] && pv < best) { best = pv; start = k; } }); }
    choose(start >= 0 ? start : 0);
  }
})();
'
