# Review pages: one self-contained HTML file per step (graphs embedded, interactive parts inline),
# plus an index of all steps and recorded decisions. Written to <run_dir>/review/.
# The pages are read-only: the analyst decides in the Claude session, and the decision is logged to decisions.jsonl.

suppressPackageStartupMessages({ library(jsonlite); library(htmltools) })

review_dir <- function(run_dir) { d <- file.path(run_dir, "review"); dir.create(d, showWarnings = FALSE); d }

esc <- function(x) htmltools::htmlEscape(as.character(x))
fmt_n <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
fmt_pct <- function(x, d = 2) paste0(formatC(as.numeric(x), format = "f", digits = d), "%")

embed_png <- function(path, alt, caption = NULL) {
  path <- gsub("/+", "/", path)
  if (!file.exists(path)) return(sprintf('<p class="missing">Plot not found: %s</p>', esc(path)))
  sprintf('<figure><img alt="%s" src="data:image/png;base64,%s"/><figcaption>%s<span class="file">%s</span></figcaption></figure>',
          esc(alt), base64enc::base64encode(path), if (is.null(caption)) "" else paste0(esc(caption), "<br/>"), esc(path))
}

html_table <- function(df, numeric_cols = integer(0), class = "") {
  head <- paste0("<tr>", paste0("<th", ifelse(seq_along(df) %in% numeric_cols, ' class="num"', ""), ">", esc(names(df)), "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(df)), function(i)
    paste0("<tr>", paste0("<td", ifelse(seq_along(df) %in% numeric_cols, ' class="num"', ""), ">", vapply(df[i, ], function(v) as.character(v), ""), "</td>", collapse = ""), "</tr>"), "")
  sprintf('<table class="%s"><thead>%s</thead><tbody>%s</tbody></table>', class, head, paste(rows, collapse = ""))
}

read_decisions <- function(run_dir, step = NULL) {
  f <- file.path(run_dir, "decisions.jsonl")
  if (!file.exists(f)) return(list())
  recs <- lapply(readLines(f, warn = FALSE), function(l) tryCatch(fromJSON(l, simplifyVector = FALSE), error = function(e) NULL))
  recs <- Filter(Negate(is.null), recs)
  if (!is.null(step)) recs <- Filter(function(r) identical(r$step, step), recs)
  recs
}

plotly_js <- function() {
  f <- list.files(system.file("htmlwidgets/lib/plotlyjs", package = "plotly"), pattern = "\\.min\\.js$", full.names = TRUE)[1]
  if (is.na(f)) stop("plotly.js not found in the plotly R package")
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

REVIEW_CSS <- "
.viz-root { color-scheme: light; --surface-0:#f4f4f2; --surface-1:#fcfcfb; --text-primary:#0b0b0b; --text-secondary:#52514e;
  --text-muted:#6f6e69; --rule:#d9d8d3; --grid:#ebeae6; --series-1:#2a78d6; --ref-dynamic:#eb6834; --ref-fixed:#008300;
  --accent:#2a78d6; --warn-bg:#fff4d6; --warn-ink:#6b4a00; --ok-bg:#e3f4ea; --ok-ink:#0f5132; }
@media (prefers-color-scheme: dark) { :root:where(:not([data-theme=\"light\"])) .viz-root { color-scheme: dark; --surface-0:#111110;
  --surface-1:#1a1a19; --text-primary:#f0efec; --text-secondary:#c3c2b7; --text-muted:#9b9a92; --rule:#3a3a37; --grid:#2a2a28;
  --series-1:#3987e5; --ref-dynamic:#d95926; --ref-fixed:#2f9e2f; --accent:#3987e5; --warn-bg:#3a2e0c; --warn-ink:#f5d58a;
  --ok-bg:#12301f; --ok-ink:#a8e0bf; } }
body { margin:0; font: 15px/1.5 system-ui, -apple-system, 'Segoe UI', sans-serif; }
.viz-root { background:var(--surface-0); color:var(--text-primary); min-height:100vh; }
header { background:var(--surface-1); border-bottom:1px solid var(--rule); padding:14px 28px; display:flex; gap:24px; align-items:baseline; flex-wrap:wrap; }
header h1 { font-size:20px; margin:0; } header .meta { color:var(--text-secondary); font-size:13px; }
header a { color:var(--accent); }
main { max-width:1180px; margin:0 auto; padding:20px 28px 60px; }
section { background:var(--surface-1); border:1px solid var(--rule); border-radius:8px; padding:16px 20px; margin:16px 0; }
section h2 { font-size:16px; margin:0 0 10px; }
.status { display:inline-block; padding:2px 10px; border-radius:999px; font-size:13px; font-weight:600; }
.status.waiting { background:var(--warn-bg); color:var(--warn-ink); } .status.done { background:var(--ok-bg); color:var(--ok-ink); }
.tiles { display:flex; gap:12px; flex-wrap:wrap; }
.tile { border:1px solid var(--rule); border-radius:8px; padding:10px 14px; min-width:150px; }
.tile .v { font-size:22px; font-weight:600; font-variant-numeric: tabular-nums; } .tile .l { color:var(--text-secondary); font-size:13px; }
table { border-collapse:collapse; font-size:14px; margin:6px 0; max-width:100%; } td { overflow-wrap:anywhere; } td:nth-child(-n+3) { white-space:nowrap; overflow-wrap:normal; } th, td { padding:4px 10px; border-bottom:1px solid var(--grid); text-align:left; }
th { color:var(--text-secondary); font-weight:600; } td.num, th.num { text-align:right; font-variant-numeric: tabular-nums; }
figure { margin:8px 0 16px; } figure img { max-width:min(100%, 760px); border:1px solid var(--rule); border-radius:4px; background:#fff; }
figcaption { color:var(--text-secondary); font-size:13px; } .file { color:var(--text-muted); font-size:12px; }
.grid2 { display:grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap:16px; }
.options li { margin:6px 0; } code { background:var(--grid); padding:1px 5px; border-radius:4px; font-size:13px; }
pre { background:var(--grid); padding:10px; border-radius:6px; overflow-x:auto; font-size:13px; }
.note { color:var(--text-secondary); font-size:13px; } .missing { color:var(--warn-ink); }
.control { display:flex; gap:14px; align-items:center; flex-wrap:wrap; margin:6px 0 10px; }
.control input[type=range] { width:420px; } .control .val { font-weight:600; font-variant-numeric: tabular-nums; min-width:200px; }
.control button { font:inherit; font-size:13px; padding:3px 10px; border:1px solid var(--rule); border-radius:6px; background:var(--surface-1); color:var(--text-primary); cursor:pointer; }
"

write_review_page <- function(run_dir, step, title, status, status_class, body, extra_head = "") {
  p_prefix <- basename(run_dir)
  html <- sprintf('<!doctype html><html lang="en"><head><meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>%s - %s</title><style>%s</style>%s</head><body><div class="viz-root">
<header><h1>%s</h1><span class="status %s">%s</span><span class="meta">Run <b>%s</b> &middot; page made %s &middot; <a href="index.html">all steps</a></span></header>
<main>%s</main></div></body></html>',
    esc(title), esc(p_prefix), REVIEW_CSS, extra_head, esc(title), status_class, esc(status), esc(p_prefix),
    format(Sys.time(), "%Y-%m-%d %H:%M"), body)
  f <- file.path(review_dir(run_dir), paste0(step, ".html"))
  writeLines(html, f, useBytes = TRUE)
  write_review_index(run_dir)
  f
}

# index.html: every step page, plus every decision recorded in decisions.jsonl.
write_review_index <- function(run_dir) {
  d <- review_dir(run_dir)
  pages <- sort(setdiff(list.files(d, pattern = "\\.html$"), "index.html"))
  dec <- read_decisions(run_dir)
  dec_rows <- if (length(dec)) data.frame(
    Time = vapply(dec, function(r) r$time %||% "", ""), Step = vapply(dec, function(r) r$step %||% "", ""),
    Analyst = vapply(dec, function(r) r$analyst %||% "", ""),
    Settings = vapply(dec, function(r) esc(paste(names(r$settings), vapply(r$settings, function(v) paste(unlist(v), collapse = ","), ""), sep = "=", collapse = "; ")), ""),
    Note = vapply(dec, function(r) esc(r$note %||% ""), ""), stringsAsFactors = FALSE) else NULL
  links <- if (length(pages)) paste0("<ul>", paste0('<li><a href="', pages, '">', esc(sub("\\.html$", "", pages)), "</a></li>", collapse = ""), "</ul>") else "<p>No step pages yet.</p>"
  body <- paste0('<section><h2>Step review pages</h2>', links, '</section><section><h2>Decisions recorded</h2>',
                 if (is.null(dec_rows)) "<p>None yet.</p>" else html_table(dec_rows), '</section>')
  html <- sprintf('<!doctype html><html lang="en"><head><meta charset="utf-8"/><title>%s - review</title><style>%s</style></head><body><div class="viz-root">
<header><h1>Run %s</h1><span class="meta">updated %s</span></header><main>%s</main></div></body></html>',
    esc(basename(run_dir)), REVIEW_CSS, esc(basename(run_dir)), format(Sys.time(), "%Y-%m-%d %H:%M"), body)
  writeLines(html, file.path(d, "index.html"), useBytes = TRUE)
}
