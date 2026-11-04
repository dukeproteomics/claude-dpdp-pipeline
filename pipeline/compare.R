# Compare tables (and optionally params) between a Claude run and a Shiny reference run.
#
# Usage:
#   Rscript pipeline/compare.R --run=<db> --ref=<db> --tables=design,precursor_raw [--params=a,b,c] [--tol=1e-6]
#
# For each table: row/column counts, column names, and value agreement
# (numeric columns within tolerance, other columns exact). Rows are compared in order.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "common.R"))
args <- parse_args()
need(args, c("run", "ref", "tables"))
tol <- as.numeric(args$tol %||% 1e-6)

rd <- function(db, t) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db); on.exit(RSQLite::dbDisconnect(con))
  if (!RSQLite::dbExistsTable(con, t)) return(NULL)
  RSQLite::dbReadTable(con, t)
}

all_ok <- TRUE
for (t in strsplit(args$tables, ",")[[1]]) {
  a <- rd(args$run, t); b <- rd(args$ref, t)
  cat(sprintf("\n== %s\n", t))
  if (is.null(a) || is.null(b)) { cat("  missing in", if (is.null(a)) "run" else "ref", "\n"); all_ok <- FALSE; next }
  cat(sprintf("  run: %d x %d   ref: %d x %d\n", nrow(a), ncol(a), nrow(b), ncol(b)))
  if (!identical(names(a), names(b))) {
    cat("  column names differ\n    only in run:", paste(setdiff(names(a), names(b)), collapse = ", "),
        "\n    only in ref:", paste(setdiff(names(b), names(a)), collapse = ", "), "\n")
    all_ok <- FALSE
  }
  if (nrow(a) != nrow(b)) { cat("  row counts differ, value comparison skipped\n"); all_ok <- FALSE; next }
  bad <- character(0)
  for (cn in intersect(names(a), names(b))) {
    x <- a[[cn]]; y <- b[[cn]]
    if (is.numeric(x) && is.numeric(y)) {
      d <- abs(x - y) > tol * pmax(1, abs(y)); d[is.na(d)] <- xor(is.na(x), is.na(y))[is.na(d)]
    } else d <- !(as.character(x) == as.character(y) | (is.na(x) & is.na(y)))
    d[is.na(d)] <- TRUE
    if (any(d)) bad <- c(bad, sprintf("%s (%d rows differ)", cn, sum(d)))
  }
  if (length(bad)) { cat("  value differences:\n", paste0("    ", bad, collapse = "\n"), "\n"); all_ok <- FALSE } else cat("  all shared columns identical\n")
}

if (!is.null(args$params)) {
  pa <- rd(args$run, "params"); pb <- rd(args$ref, "params")
  cat("\n== params\n")
  for (p in strsplit(args$params, ",")[[1]]) {
    va <- as.character(pa[[p]] %||% NA); vb <- as.character(pb[[p]] %||% NA)
    flag <- if (identical(va, vb)) "  " else "!!"
    if (flag == "!!") all_ok <- FALSE
    cat(sprintf("  %s %-28s run=%-20s ref=%s\n", flag, p, va, vb))
  }
}
cat(if (all_ok) "\nRESULT: MATCH\n" else "\nRESULT: DIFFERENCES FOUND\n")
