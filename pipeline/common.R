# Shared helpers for the Claude DPDP pipeline.
# Each step script sources this file, then calls the existing shiny_DPDP functions.
# Run with plain `Rscript` (not --vanilla) so the shared library /home/R_Packages is used.

APP_DIR <- Sys.getenv("DPDP_APP", "/home/user/Claude_Shiny/shiny_DPDP")
START_WD <- getwd()   # use_app() changes the working directory; resolve user paths against this

abs_path <- function(x) normalizePath(if (grepl("^/", x)) x else file.path(START_WD, x), mustWork = TRUE)

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(jsonlite)
  library(RSQLite)   # app helpers call dbConnect() etc. unqualified
})

# App code sources sibling files by relative path, so it must run from the app directory.
use_app <- function(files) {
  setwd(APP_DIR)
  for (f in files) source(file.path(APP_DIR, f), local = globalenv())
  invisible(TRUE)
}

# Parse --key=value or --flag arguments into a named list.
parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (a in args) {
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a)
    kv <- sub("^--", "", a)
    if (grepl("=", kv)) out[[sub("=.*", "", kv)]] <- sub("^[^=]*=", "", kv) else out[[kv]] <- TRUE
  }
  out
}

need <- function(args, keys) {
  miss <- setdiff(keys, names(args))
  if (length(miss)) stop("Missing required argument(s): ", paste0("--", miss, collapse = ", "))
}

# Each step writes its outputs to <run_dir>/steps/<step>/.
step_dir <- function(run_dir, step) {
  d <- file.path(run_dir, "steps", step)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# summary.json is the main thing Claude reads after a step.
write_summary <- function(summary, dir) {
  summary$generated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  summary$app_commit <- tryCatch(
    system2("git", c("-C", APP_DIR, "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_)
  path <- file.path(dir, "summary.json")
  write_json(summary, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
  message("Summary written: ", path)
  invisible(path)
}

# Append every option choice to <run_dir>/decisions.jsonl.
log_decision <- function(run_dir, step, settings, note = "", analyst = NULL) {
  # A decision sent from the review page (review/pending_decision.json, written by pipeline/review_server.py) for this step is
  # archived to review/decisions_done/ and noted in the log. Its analyst/note fill in when the command line didn't give them.
  pend <- file.path(run_dir, "review", "pending_decision.json")
  page <- NULL
  if (file.exists(pend)) {
    pd <- tryCatch(fromJSON(pend), error = function(e) NULL)
    if (!is.null(pd) && identical(pd$step, step)) {
      done <- file.path(run_dir, "review", "decisions_done"); dir.create(done, showWarnings = FALSE)
      file.rename(pend, file.path(done, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", step, ".json")))
      page <- list(choice = pd$choice, submitted = pd$submitted)
      if (is.null(analyst) && nzchar(pd$analyst %||% "")) analyst <- pd$analyst
      if (!nzchar(note %||% "") && nzchar(pd$note %||% "")) note <- pd$note
    }
  }
  rec <- list(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step = step,
              analyst = analyst %||% Sys.getenv("USER"), settings = settings, note = note,
              source = if (is.null(page)) "chat" else "review page", page_choice = page$choice)
  cat(toJSON(rec, auto_unbox = TRUE), "\n", file = file.path(run_dir, "decisions.jsonl"),
      append = TRUE, sep = "")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Run one of the app's *_bg functions in a fresh R process, as the app does with callr::r_bg.
# stderr goes to <error_path>/error_<name>.txt; on failure the end of that log is shown.
run_bg <- function(func, args, name, db_path) {
  err <- file.path(get_param("error_path", db_path), str_c("error_", name, ".txt"))
  message("-> ", name)
  res <- tryCatch(
    callr::r(func, args = args, stderr = err, wd = APP_DIR,
             env = c(callr::rcmd_safe_env(), R_LIBS_USER = Sys.getenv("R_LIBS_USER"))),
    error = function(e) {
      if (file.exists(err)) message(paste(tail(readLines(err), 15), collapse = "\n"))
      stop(name, " failed: ", conditionMessage(e), call. = FALSE)
    })
  invisible(res)
}

# Convert "TRUE"/"1"/"yes" style CLI values to logical.
as_flag <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")

# Build step settings: app defaults < values copied from another .db (--settings_from) < CLI args.
# `spec` is a named list of defaults; the type of each default decides how CLI strings are converted.
resolve_settings <- function(spec, args) {
  out <- spec
  if (!is.null(args$settings_from)) {
    src <- abs_path(args$settings_from)
    con <- RSQLite::dbConnect(RSQLite::SQLite(), src); on.exit(RSQLite::dbDisconnect(con))
    sp <- RSQLite::dbReadTable(con, "params")
    for (k in intersect(names(spec), names(sp))) out[[k]] <- sp[[k]][1]
  }
  for (k in intersect(names(spec), names(args))) out[[k]] <- args[[k]]
  for (k in names(spec)) {
    if (is.logical(spec[[k]])) out[[k]] <- as_flag(out[[k]])
    else if (is.numeric(spec[[k]])) out[[k]] <- as.numeric(out[[k]])
    else out[[k]] <- as.character(out[[k]])
  }
  out
}

# All plots a step produced, as clean absolute paths (only those that exist).
existing_plots <- function(qc_path, files) {
  p <- gsub("/+", "/", file.path(qc_path, files))
  p[file.exists(p)]
}

# The imputation that was APPLIED (Step 7 summary), not the param (a preview also writes impute_type).
applied_impute_type <- function(db_path) {
  f <- file.path(dirname(db_path), "steps", "07_impute", "summary.json")
  if (!file.exists(f)) return(NA_character_)
  s <- jsonlite::fromJSON(f, simplifyVector = FALSE)
  if (identical(s$action, "apply")) s$settings$impute_type %||% NA_character_ else NA_character_
}

# Output level: Protein tables (protein_<norm>, stats_comp$Final_Table_Name) or Peptide tables (peptide_impute_<norm>,
# Final_Table_Name_Peptide; e.g. phospho). feature_key() identifies a row: the protein accession, or the peptide's grouping columns
# (the ones the app's collapse_precursor_ptm_raw / collapse_precursor_raw2 group by).
is_peptide_output <- function(p) identical(p$data_output, "Peptide")
base_table <- function(p, norm = p$stat_norm) if (is_peptide_output(p)) paste0("peptide_impute_", norm) else paste0("protein_", norm)
final_table <- function(stats_comp, i, p) if (is_peptide_output(p)) stats_comp$Final_Table_Name_Peptide[i] else stats_comp$Final_Table_Name[i]
PEPTIDE_KEY <- c("Accession", "Description", "Name", "Genes", "Sequence", "PTM_Loc", "Protein_PTM_Loc", "PeptidePosition")
# PTM-enriched runs (Parameters ptm = TRUE): normalization factors and Stats use only the rows whose Sequence matches ptm_grep
# (user rule 2026-11-05: enriched samples, the non-PTM peptides are not supposed to be there). Rows stay in the data tables.
is_ptm_run <- function(p) isTRUE(as_flag(p$ptm %||% FALSE)) && nzchar(p$ptm_grep %||% "")
is_ptm_row <- function(d, p) grepl(p$ptm_grep, d$Sequence, ignore.case = TRUE)
ptm_report <- function(p) is_peptide_output(p) && isTRUE(as_flag(p$checkbox_report_ptm %||% FALSE)) && !p$peptide_report_grep %in% c("", "Enter value here")
# The word for one output row on the pages: "protein", or "peptide" for Peptide output (phospho / PTM runs roll precursors up to peptides).
feat <- function(p, plural = TRUE, cap = FALSE) { w <- paste0(if (is_peptide_output(p)) "peptide" else "protein", if (plural) "s" else "")
  if (cap) paste0(toupper(substring(w, 1, 1)), substring(w, 2)) else w }
feature_key <- function(d, p) if (is_peptide_output(p)) do.call(paste, c(lapply(d[intersect(PEPTIDE_KEY, names(d))], as.character), sep = "\r")) else as.character(d$Accession)
