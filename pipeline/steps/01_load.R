# Step 1 — Load
# Same as the app's Load tab: design file -> new project + .db + `design` table,
# then data file -> source/format detection + `precursor_raw` (or protein_raw) table.
#
# Usage:
#   Rscript pipeline/steps/01_load.R --design=<xlsx> --data=<tsv> --prefix=<name> \
#       [--outdir=$DPDP_DATA (default /home/user/Claude_Shiny/data)] [--overwrite]
#       [--sp_fill=raw|background] [--sp_qvalue=0.01]
#
# Spectronaut precursor exports (pipeline/spectronaut.R): raw and background-signal exports are accepted, imputed (global / run-wise)
# exports are refused. The type is read from the run-wise EG.Qvalue columns (background = values present with q > sp_qvalue);
# --sp_fill states it when the export has no q-values (then only raw is accepted). Filled values are recorded in table
# precursor_fill (same rows as precursor_raw, 1 = Spectronaut background signal) and param sp_fill_mode.
#
# Output: <outdir>/<prefix>/<prefix>.db, plus steps/01_load/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "01_load.R"))
source(file.path(here, "..", "spectronaut.R"))

args <- parse_args()
need(args, c("design", "data", "prefix"))
design_file <- normalizePath(args$design, mustWork = TRUE)
data_file <- normalizePath(args$data, mustWork = TRUE)
outdir <- normalizePath(args$outdir %||% Sys.getenv("DPDP_DATA", "/home/user/Claude_Shiny/data"), mustWork = TRUE)
prefix <- args$prefix
if (!grepl("^[A-Za-z0-9_.-]+$", prefix)) stop("Prefix may only contain letters, digits, _ . -")

# --- 0. Spectronaut export type (before anything is written) ---------------
sp_fill <- NULL
if (!is.null(args$sp_fill) && !args$sp_fill %in% c("raw", "background")) stop("--sp_fill must be raw or background")
hdr <- names(data.table::fread(data_file, sep = "\t", nrows = 0))
if ("EG.PrecursorId" %in% hdr && length(grep("EG.TotalQuantity", hdr)) > 1) {
  df0 <- as.data.frame(data.table::fread(data_file, sep = "\t", fill = TRUE, showProgress = FALSE))
  sp_fill <- sp_detect_fill(df0, q_cut = as.numeric(args$sp_qvalue %||% 0.01), declared = args$sp_fill)
  rm(df0)
  if (!is.null(sp_fill$error)) stop("Spectronaut export not accepted: ", sp_fill$error, call. = FALSE)
  cat(sprintf("Spectronaut export type: %s (%s values present, %s filled by Spectronaut, %s gaps)\n", sp_fill$mode,
              format(sp_fill$values_present, big.mark = ","), format(sp_fill$filled, big.mark = ","), format(sp_fill$gaps, big.mark = ",")))
}

run_dir <- file.path(outdir, prefix)
if (dir.exists(run_dir)) {
  # The app's create_dir() deletes existing folders, so only allow this on request.
  if (!isTRUE(args$overwrite)) stop("Project folder already exists: ", run_dir, "  (use --overwrite to replace it)")
  unlink(run_dir, recursive = TRUE)
}

use_app(c("Shiny_File.R", "Shiny_Startup.R", "Shiny_Data.R"))

# --- 1a. Design file (app: load_design_file + file_set) ----------------------
params <- create_default_params(volumes = character(0), python_path = "")
params$file_prefix <- prefix
params$design_path <- str_c(dirname(design_file), "/")
params$design_file <- design_file
params$data_path <- str_c(run_dir, "/")
params$database_dir <- run_dir
db_path <- file.path(run_dir, str_c(prefix, ".db"))
params$database_path <- db_path
params$error_path <- str_c(params$data_path, "Error")

invisible(create_dir(params$data_path))
invisible(create_dir(params$error_path))
invisible(write_params(params, db_path))
excel_to_db(design_file, "design", db_path)
file_set(db_path)                  # Backup/, Extra/, QC/, String/, Phos/; archives app R files
save_data(design_file, db_path)    # backup copy of the design file

# --- 1b. Data file (app: load_data_file -> load_unknown_data_bg) -------------
params <- get_params(db_path)
params$data_source <- "unkown"     # same placeholder value the app uses
params$data_file <- basename(data_file)
write_table_try("params", params, db_path)
load_unknown_data_bg(list(datapath = data_file), db_path)
if (!is.null(sp_fill)) {
  pr <- read_table("precursor_raw", db_path)
  if (nrow(pr) != nrow(sp_fill$mask)) stop("precursor_raw has ", nrow(pr), " rows but the export has ", nrow(sp_fill$mask))
  fill_tbl <- data.frame(PrecursorId = pr$EG.PrecursorId, 1L * sp_fill$mask, check.names = FALSE)
  names(fill_tbl)[-1] <- grep("EG.TotalQuantity", names(pr), value = TRUE)
  write_table_try("precursor_fill", fill_tbl, db_path)
  params <- get_params(db_path)
  params$sp_fill_mode <- sp_fill$mode; params$sp_fill_qvalue <- sp_fill$q_cut
  write_table_try("params", params, db_path)
}

# --- Summary for Claude ------------------------------------------------------
params <- get_params(db_path)
design <- read_table("design", db_path)
raw_table <- if (params$raw_data_format == "protein") "protein_raw" else "precursor_raw"
if (!table_exists(raw_table, db_path)) stop("Data file was not recognised (no ", raw_table, " table written).")
df <- read_table(raw_table, db_path)

qcols <- grep(if (raw_table == "precursor_raw") "EG.TotalQuantity" else "PG.Quantity", names(df), value = TRUE)
d <- design[order(design$Raw_Order), ]

# Sample ID check. The app's version drops "_N" suffixes; this one also checks the full ID.
short_id <- gsub("_.+", "", d$ID)
full_ok <- mapply(function(id, col) {
  parts <- strsplit(id, "_")[[1]]
  pat <- if (length(parts) > 1) str_c(parts[1], "_0*", parts[2]) else parts[1]
  # also accept the replicate number at the end of the run name (e.g. ..._E5H50Y45_180K_2Th3p5ms_01.htrms)
  grepl(pat, col) || (length(parts) > 1 && grepl(parts[1], col, fixed = TRUE) && grepl(str_c("_0*", parts[2], "\\.(htrms|raw)"), col))
}, d$ID, qcols[seq_len(nrow(d))])
id_table <- data.frame(Raw_Order = d$Raw_Order, Label = d$Label, Group = d$Group, ID = d$ID,
                       data_column = sub("\\.htrms.*|\\.raw.*", "", qcols[seq_len(nrow(d))]),
                       app_check = mapply(grepl, short_id, qcols[seq_len(nrow(d))]),
                       full_id_check = unname(full_ok))

x <- suppressWarnings(sapply(df[qcols], as.numeric))
missing <- setNames(as.integer(colSums(is.na(x))), d$Label)

org <- if ("PG.Organisms" %in% names(df)) {
  u <- unique(df[c("PG.ProteinAccessions", "PG.Organisms")])
  sort(table(u$PG.Organisms), decreasing = TRUE)
} else NULL

warnings_out <- character(0)
one <- names(which(table(design$Group) == 1))
if (length(one)) warnings_out <- c(warnings_out, sprintf("STOP: group(s) with one sample: %s. The app's filter fails on a one-sample group (Shiny_Filter.R:274 drops to a vector). Merge it into another group or give it more samples in the sample sheet.", paste(one, collapse = ", ")))
if (!is.null(sp_fill)) warnings_out <- c(warnings_out, sp_fill$notes)
if (length(qcols) != nrow(design)) warnings_out <- c(warnings_out, sprintf("Design has %d samples but data has %d sample columns.", nrow(design), length(qcols)))
if (!all(id_table$app_check)) warnings_out <- c(warnings_out, "Sample IDs do not match the data columns in Raw_Order.")
if (!all(id_table$full_id_check)) warnings_out <- c(warnings_out, "Full-ID check failed for some samples (replicate suffix differs; OK for pooled QC).")
mf <- missing / nrow(df)
hi <- names(mf)[mf > median(mf) * 1.5]
if (length(hi)) warnings_out <- c(warnings_out, sprintf("Higher missingness than the median sample (>1.5x): %s", paste(hi, collapse = ", ")))

sd <- step_dir(run_dir, "01_load")
summary <- list(
  step = "01_load",
  run_dir = run_dir, db_path = db_path,
  inputs = list(design = design_file, data = data_file),
  detected = list(data_source = params$data_source, raw_data_format = params$raw_data_format,
                  data_table_format = params$data_table_format, ptm = as.logical(params$ptm),
                  data_output = params$data_output),
  counts = list(rows = nrow(df), sample_columns = length(qcols), design_samples = nrow(design),
                groups = as.list(table(design$Group)),
                proteins = if ("PG.ProteinAccessions" %in% names(df)) length(unique(df$PG.ProteinAccessions)) else NA,
                peptides = if ("PEP.StrippedSequence" %in% names(df)) length(unique(df$PEP.StrippedSequence)) else NA,
                precursors = if ("EG.PrecursorId" %in% names(df)) length(unique(df$EG.PrecursorId)) else NA,
                duplicate_precursors = if ("EG.PrecursorId" %in% names(df)) sum(duplicated(df$EG.PrecursorId)) else NA),
  spectronaut_fill = if (!is.null(sp_fill)) list(mode = sp_fill$mode, q_cut = sp_fill$q_cut, has_qvalue = sp_fill$has_qvalue,
    has_isimputed = sp_fill$has_isimputed, declared = sp_fill$declared, values_present = sp_fill$values_present, gaps = sp_fill$gaps,
    filled = sp_fill$filled, measured = sp_fill$measured %||% sp_fill$values_present, identified_but_empty = sp_fill$identified_but_empty,
    filled_per_sample = if (!is.null(sp_fill$filled_per_sample)) as.list(setNames(as.integer(sp_fill$filled_per_sample), d$Label)) else NULL) else NULL,
  missing_per_sample = as.list(missing),
  missing_fraction_per_sample = as.list(round(mf, 4)),
  organisms_top = if (!is.null(org)) as.list(head(org, 8)) else NULL,
  sample_id_check = id_table,
  warnings = warnings_out,
  tables = local({ con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path); on.exit(RSQLite::dbDisconnect(con)); RSQLite::dbListTables(con) })
)
write_summary(summary, sd)
log_decision(run_dir, "01_load", c(list(prefix = prefix, design = design_file, data = data_file),
  if (!is.null(sp_fill)) list(sp_fill_mode = sp_fill$mode, sp_fill_declared = sp_fill$declared, sp_qvalue = sp_fill$q_cut)))

cat("\n==== Load complete ====\n")
cat("DB:", db_path, "\n")
cat(sprintf("Detected: %s / %s / %s format / PTM=%s -> %s output\n", params$data_source, params$raw_data_format,
            params$data_table_format, as.logical(params$ptm), params$data_output))
cat(sprintf("%d rows, %d sample columns, %d design samples\n", nrow(df), length(qcols), nrow(design)))
if (length(warnings_out)) cat("Warnings:\n", paste0("  - ", warnings_out, collapse = "\n"), "\n") else cat("No warnings.\n")
cat("Review page:", review_load(db_path), "\n>>> Waiting for the analyst's review (accept with pipeline/accept.R).\n")
