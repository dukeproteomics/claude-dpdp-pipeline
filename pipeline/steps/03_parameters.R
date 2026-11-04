# Step 3 — Parameters
# Same as the app's Parameters tab "Accept Parameters" button (server.R observeEvent(input$accept_parameters)):
#   parameter_widget_save -> set_sample_groups -> prepare_data -> order_rename_columns
#   -> meta_data("start") -> bar/box plots -> start histogram.
#
# Usage:
#   Rscript pipeline/steps/03_parameters.R --db=<db> [--primary_group=Full|Primary] [--data_output=Protein|Peptide]
#       [--ptm=FALSE] [--ptm_grep=Phospho] [--ptm_local=0.75] [--multi_tmt=FALSE]
#       [--peptide_select=Razor|Unique|Shared] [--use_isoform=FALSE]
# Options left out keep the value detected at Load (data_output, ptm) or the app default.
#
# Output: tables sample_groups, precursor_start; updated design; QC/Precursor_Start_{barplot,boxplot,Histogram}.png
#         steps/03_parameters/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "01_load.R"))
source(file.path(here, "..", "reviews", "03_parameters.R"))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
use_app(c("Shiny_File.R", "Shiny_Data.R", "Shiny_Setup.R", "Shiny_Plots.R"))

params <- get_params(db_path)
run_dir <- dirname(db_path)

# --- 3a. Widget values (app: parameter_widget_save) ---------------------------
choices <- list(primary_group = c("Full", "Primary"), data_output = c("Protein", "Peptide"),
                peptide_select = c("Razor", "Unique", "Shared"))
settings <- list(
  primary_group  = args$primary_group  %||% "Full",
  data_output    = args$data_output    %||% params$data_output,
  ptm            = if (!is.null(args$ptm)) as_flag(args$ptm) else as.logical(params$ptm),
  ptm_grep       = args$ptm_grep       %||% "Phospho",
  ptm_local      = as.numeric(args$ptm_local %||% 0.75),
  multi_tmt      = if (!is.null(args$multi_tmt)) as_flag(args$multi_tmt) else FALSE,
  peptide_select = args$peptide_select %||% "Razor",
  use_isoform    = if (!is.null(args$use_isoform)) as_flag(args$use_isoform) else FALSE
)
for (k in names(choices)) if (!settings[[k]] %in% choices[[k]])
  stop(sprintf("--%s must be one of: %s", k, paste(choices[[k]], collapse = ", ")))
if (settings$use_isoform) stop("use_isoform=TRUE is not supported yet in the pipeline (app's isoform_to_isoform() needs work).")
for (k in names(settings)) params[[k]] <- settings[[k]]
invisible(write_params(params, db_path))

# --- 3b. Sample groups (app: set_sample_groups) --------------------------------
run_bg(set_sample_groups_bg, list(db_path), "setsamplegroups", db_path)

# --- 3c. Prepare data (app: prepare_data). Only the Spectronaut precursor paths are wired up so far.
params <- get_params(db_path)
sample_error <- FALSE
if (params$raw_data_format == "precursor" && params$data_source == "SP" && !(params$data_output == "Peptide" && params$ptm)) {
  sample_error <- run_bg(precursor_to_precursor_bg, list(db_path), "preparedata", db_path)
} else if (params$raw_data_format == "precursor" && params$data_source == "SP" && params$data_output == "Peptide" && params$ptm) {
  sample_error <- run_bg(precursor_to_precursor_ptm_bg, list(db_path), "preparedata", db_path)
} else {
  stop(sprintf("Data type not supported yet in the pipeline: source=%s format=%s output=%s ptm=%s",
               params$data_source, params$raw_data_format, params$data_output, params$ptm))
}
params <- get_params(db_path)
params$current_data_format <- "precursor"
invisible(write_params(params, db_path))

# --- 3d. Order/rename sample columns + sample ID check (app: order_rename_columns)
res <- run_bg(order_rename_columns_bg, list("precursor_start", db_path), "orderrename", db_path)
param_update("info_col_precursor", res[[1]], db_path)
sample_id_alert <- isTRUE(res[[2]])

# --- 3e. Counts + plots (app: meta_data("start"), parameter_create_plots, filter_histogram_plot)
run_bg(meta_data_bg, list("precursor_start", "start", db_path), "start_meta", db_path)
qc_path <- get_param("qc_path", db_path)
run_bg(bar_plot, list("precursor_start", "Precursor_Start", qc_path, db_path), "barplot", db_path)
run_bg(box_plot, list("precursor_start", "Precursor_Start", qc_path, db_path), "boxplot", db_path)
run_bg(histogram_plot, list("precursor_start", "Precursor_Start_Histogram", db_path), "histogram", db_path)

# --- Summary for Claude ----------------------------------------------------------
params <- get_params(db_path)
groups <- read_table("sample_groups", db_path)
df <- read_table("precursor_start", db_path)
n_info <- as.integer(params$info_col_precursor)
x <- as.matrix(df[(n_info + 1):ncol(df)])
lx <- log2(x); lx[!is.finite(lx)] <- NA
per_sample <- data.frame(sample = colnames(x),
                         measured = colSums(!is.na(x) & x > 0, na.rm = TRUE),
                         missing = colSums(is.na(x) | x == 0, na.rm = TRUE),
                         total_log2 = round(log2(colSums(x, na.rm = TRUE)), 3),
                         median_log2 = round(apply(lx, 2, median, na.rm = TRUE), 3))
per_sample$total_vs_median <- round(2^(per_sample$total_log2 - median(per_sample$total_log2)), 3)

warnings_out <- character(0)
if (sample_id_alert) warnings_out <- c(warnings_out, "STOP: sample ID order does not match the sample list.")
if (isTRUE(sample_error)) warnings_out <- c(warnings_out, "STOP: number of sample columns extracted is not as expected.")
off <- per_sample$sample[abs(log2(per_sample$total_vs_median)) > log2(1.5)]
if (length(off)) warnings_out <- c(warnings_out, sprintf("Total intensity differs >1.5x from the median sample: %s", paste(off, collapse = ", ")))

sd <- step_dir(run_dir, "03_parameters")
plots <- file.path(qc_path, c("Precursor_Start_barplot.png", "Precursor_Start_boxplot.png", "Precursor_Start_Histogram.png"))
plots <- gsub("/+", "/", plots)
summary <- list(
  step = "03_parameters", db_path = db_path,
  settings = settings,
  info_columns = n_info, info_column_names = names(df)[seq_len(n_info)],
  counts_start = list(precursors = params$meta_precursor_start, peptides = params$meta_peptide_start,
                      proteins = params$meta_protein_start),
  sample_groups = groups[c("Group", "Count", "start", "end")],
  per_sample = per_sample,
  intensity = list(log2_mean = params$intensity_mean, log2_sd = params$intensity_sd,
                   intensity_cutoff = params$intensity_cutoff, intensity_cutoff_sd = params$intensity_cutoff_sd,
                   bottom_x = params$bottom_x, total_na = params$total_na, total_misaligned = params$total_misaligned),
  sample_id_alert = sample_id_alert, sample_error = isTRUE(sample_error),
  plots = plots[file.exists(plots)],
  warnings = warnings_out
)
write_summary(summary, sd)
log_decision(run_dir, "03_parameters", settings)

cat("\n==== Parameters complete ====\n")
cat(sprintf("Start: %s precursors, %s peptides, %s proteins | %d info columns\n",
            params$meta_precursor_start, params$meta_peptide_start, params$meta_protein_start, n_info))
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
if (length(warnings_out)) cat("Warnings:\n", paste0("  - ", warnings_out, collapse = "\n"), "\n") else cat("No warnings.\n")
cat("Review page:", review_parameters(db_path), "\n>>> Waiting for the analyst's review (accept with pipeline/accept.R).\n")
