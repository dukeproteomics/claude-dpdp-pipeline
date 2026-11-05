# Step 6 — Normalize
# Same as the app's Normalize tab:
#   --action=preview : "Set Norm Parameters" (norm_widget_save + norm_filter_bg + NormData bar plot).
#                      Builds precursor_normdata = the subset used to calculate normalization factors.
#   --action=apply   : "Apply Normalization" (norm_apply_widget_save + norm_apply_bg + bar plot).
#
# Usage:
#   Rscript pipeline/steps/06_normalize.R --db=<db> --action=preview [--settings_from=<db>]
#       [--norm_exclude=not_used|accession|description|name|genes|sequence] [--norm_exclude_grep=...]
#       [--norm_include=...] [--norm_include_grep=...]
#   Rscript pipeline/steps/06_normalize.R --db=<db> --action=apply --norm_type=<method> [--settings_from=<db>] \
#       --decided_by=user [--analyst=<name>] [--note="<why>"]
# Both actions (re)build the review page <run_dir>/review/06_normalize.html.
# Methods: sl, tmm, sltmm, directlfq, quantile, lr, loess, vsn, mi, ai
#   (sltmm: SL is applied here; the app applies TMM after imputation.)
#
# Output: precursor_normdata, precursor_norm_<method> (+ precursor_norm_impute on the first pass);
#         QC/Precursor_NormData_barplot.png, QC/Precursor_<method>_barplot.png; steps/06_normalize/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "background.R"))
source(file.path(here, "..", "reviews", "06_normalize.R"))

args <- parse_args()
need(args, c("db", "action"))
if (!args$action %in% c("preview", "apply")) stop("--action must be preview or apply")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_Norm.R", "Shiny_Plots.R"))
if (get_param("raw_data_format", db_path) != "precursor") stop("Normalize step is only wired up for precursor data so far.")
if (!table_exists("precursor_filter", db_path)) stop("Run Step 5 (Filter) first: precursor_filter is missing.")

fields <- c("not_used", "accession", "description", "name", "genes", "sequence")
methods <- c("sl", "tmm", "sltmm", "directlfq", "quantile", "lr", "loess", "vsn", "mi", "ai")
# Default exclude (user rule 2026-11-04): only added material -- digestion trypsin, human contaminant keratins,
# bovine alpha/beta casein spike. The app's default "trypsin|keratin|casein" also hits the sample's own proteins.
EXCLUDE_DEFAULT <- NORM_DEFAULT_EXCLUDE   # defined in pipeline/reviews/06_normalize.R
# PTM runs (user rule 2026-11-05): factors from the PTM peptides only (include sequence ~ ptm_grep), on top of the default exclude.
p0 <- get_params(db_path)
ptm <- is_ptm_run(p0)
spec <- list(norm_exclude = "description", norm_exclude_grep = EXCLUDE_DEFAULT,
             norm_include = if (ptm) "sequence" else "not_used", norm_include_grep = if (ptm) p0$ptm_grep else "Phospho",
             norm_ptm = FALSE, norm_ptm_grep = "Phospho")
settings <- resolve_settings(spec, args)
for (k in c("norm_exclude", "norm_include")) if (!settings[[k]] %in% fields)
  stop(sprintf("--%s must be one of: %s", k, paste(fields, collapse = ", ")))

if (settings$norm_exclude == "description" && tolower(settings$norm_exclude_grep) == "trypsin|keratin|casein")
  message("WARNING: norm_exclude_grep 'trypsin|keratin|casein' also excludes the sample's own proteins (casein kinases, ",
          "trypsin inhibitors, sample-species keratins). The user's rule is the default: --norm_exclude_grep='", EXCLUDE_DEFAULT, "'")
p <- get_params(db_path)
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))
qc_path <- p$qc_path
info_columns <- as.integer(p$info_col_precursor)
sd <- step_dir(run_dir, "06_normalize")

sample_cols <- function(df) df[(info_columns + 1):ncol(df)]
col_totals <- function(df) colSums(sample_cols(df), na.rm = TRUE)

# --- preview: build normalization subset (always rebuilt so it matches the current settings) ---
# Factors from measured values only: for background runs the subset is built from precursor_filter_measured (background blanked).
bg_run <- is_background_run(db_path)
run_bg(norm_filter_bg, list(if (bg_run) "precursor_filter_measured" else "precursor_filter", "precursor_normdata", p, db_path), "normfilter", db_path)
run_bg(bar_plot, list("precursor_normdata", "Precursor_NormData", qc_path, db_path), "normdatabarplot", db_path)

filt <- read_table("precursor_filter", db_path)
normdata <- read_table("precursor_normdata", db_path)
excluded <- filt[!filt$PrecursorId %in% normdata$PrecursorId, ]
top_excluded <- if (nrow(excluded)) {
  t <- sort(table(excluded$Description), decreasing = TRUE)
  head(data.frame(description = names(t), precursors = as.integer(t)), 15)
} else NULL
if (bg_run) filt <- read_table("precursor_filter_measured", db_path)   # same rows; exclusion list from the measured table
sl_target <- mean(col_totals(normdata))
sl_factors <- round(sl_target / col_totals(normdata), 4)

norm_subset <- list(precursors_in_filter = nrow(filt), precursors_used_for_factors = nrow(normdata),
                    precursors_excluded = nrow(excluded),
                    proteins_excluded = length(unique(excluded$Accession)),
                    top_excluded_descriptions = top_excluded,
                    sl_factors = as.list(sl_factors))

if (args$action == "preview") {
  summary <- list(step = "06_normalize", action = "preview", db_path = db_path, settings = settings,
                  normalization_subset = norm_subset,
                  plots = existing_plots(qc_path, "Precursor_NormData_barplot.png"),
                  decision_needed = "User reviews the excluded proteins and NormData bar plot, then approves --action=apply --norm_type=<method>.")
  write_summary(summary, sd)
  cat("\n==== Normalize preview ====\n")
  cat(sprintf("Norm subset: %d of %d precursors (%d excluded, %d proteins)\n", nrow(normdata), nrow(filt), nrow(excluded), norm_subset$proteins_excluded))
  page <- review_normalize(db_path)
  cat("Plot:", summary$plots, "\nReview page:", page, "\n>>> Waiting for the user's approval.\n")
  quit(save = "no")
}

# --- apply ---------------------------------------------------------------------
if (!identical(args$decided_by, "user")) stop("Apply needs the analyst's approval: run --action=preview, show the review page, then pass --decided_by=user.")
method <- args$norm_type
if (is.null(method) && !is.null(args$settings_from)) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), abs_path(args$settings_from))
  method <- trimws(strsplit(RSQLite::dbReadTable(con, "params")$norm_type, ",")[[1]][1]); RSQLite::dbDisconnect(con)
}
if (is.null(method)) stop("--norm_type is required for apply")
if (!method %in% methods) stop("--norm_type must be one of: ", paste(methods, collapse = ", "),
                               "  (the app's 'mti' option has no mti_normalize() function)")
if (method %in% c("ai", "mi")) stop("--norm_type=", method, ": the app's norm_apply_bg() calls ", method,
                                    "_normalize() with the wrong arguments, so it fails. Choose another method.")
if (bg_run && !method %in% c("sl", "tmm", "sltmm"))
  stop("--norm_type=", method, ": on a background export only factor methods (sl, tmm, sltmm) are allowed, because only they can be ",
       "calculated from measured values and then applied to the background values. ", method, " transforms the data as a whole.")
if (method == "directlfq") warning("directlfq needs python_path, which isn't set up in the pipeline yet.")

# app: norm_apply_widget_save -> new method is put in front of the existing list
p <- get_params(db_path)
norm_list <- unique(c(method, trimws(strsplit(p$norm_type, ",")[[1]])))
p$norm_type <- toString(norm_list)
invisible(write_params(p, db_path))
p <- get_params(db_path)

if (length(norm_list) == 2) {   # app: first pass also makes the impute-only table
  run_bg(norm_apply_bg, list("precursor_filter", "precursor_normdata", "precursor_norm_impute", info_columns, p, db_path, "impute"),
         "normapplyimpute", db_path)
}
new_table <- str_c("precursor_norm_", method)
run_bg(norm_apply_bg, list("precursor_filter", "precursor_normdata", new_table, info_columns, p, db_path, method), "normapply", db_path)
run_bg(bar_plot, list(new_table, str_c("Precursor_", method), qc_path, db_path), "normbarplot", db_path)

normed <- read_table(new_table, db_path)
before <- col_totals(read_table("precursor_filter", db_path)); after <- col_totals(normed)
summary <- list(
  step = "06_normalize", action = "apply", db_path = db_path, settings = settings,
  norm_type = method, norm_type_list = p$norm_type, table = new_table,
  note = if (method == "sltmm") "SL applied now; the app applies TMM after imputation (Shiny_Impute.R)." else NULL,
  normalization_subset = norm_subset,
  factors_from = if (bg_run) "measured values only; background values excluded, then normalized with the same factors" else "measured values",
  totals_before = as.list(round(before)), totals_after = as.list(round(after)),
  total_cv_pct_before = round(sd(before) / mean(before) * 100, 2),
  total_cv_pct_after = round(sd(after) / mean(after) * 100, 2),
  plots = existing_plots(qc_path, c("Precursor_NormData_barplot.png", str_c("Precursor_", method, "_barplot.png")))
)
write_summary(summary, sd)
log_decision(run_dir, "06_normalize", c(settings, list(norm_type = method)), analyst = args$analyst,
             note = paste(c(args$note, if (!is.null(args$settings_from)) paste("settings copied from", args$settings_from)), collapse = "; "))

cat("\n==== Normalization applied ====\n")
cat(sprintf("Method: %s -> %s   (norm_type list: %s)\n", method, new_table, p$norm_type))
cat(sprintf("Sample total CV: %.2f%% before -> %.2f%% after\n", summary$total_cv_pct_before, summary$total_cv_pct_after))
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
page <- review_normalize(db_path)
cat("Review page:", page, "\n")
