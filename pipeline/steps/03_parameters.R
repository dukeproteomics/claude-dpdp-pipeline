# Step 3 — Parameters
# Same as the app's Parameters tab "Accept Parameters" button (server.R observeEvent(input$accept_parameters)):
#   parameter_widget_save -> set_sample_groups -> prepare_data -> order_rename_columns
#   -> meta_data("start") -> bar/box plots -> start histogram.
#
# Usage:
#   Rscript pipeline/steps/03_parameters.R --db=<db> [--from_page] [--force] [--primary_group=Full|Primary]
#       [--data_output=Protein|Peptide] [--ptm=FALSE] [--ptm_grep=Phospho] [--ptm_local=0.75]
#   (peptide_select, use_isoform: Proteome Discoverer only; multi_tmt: TMT only. Hidden and kept at the app default otherwise.)
# Every option, its default and its checks are in pipeline/options/03_parameters.R (also the page's Options panel).
# Options left out take their default: data_output and ptm as detected at Load, the rest as in the app.
# --from_page: take the values the analyst sent with "Run with these settings" from the page (review/pending_decision.json).
# Accept (pipeline/accept.R --step=03_parameters) only accepts the settings of the last run of this step.
# Later steps already run: the step stops (their results would no longer match); --force re-runs anyway.
#
# Output: tables sample_groups, precursor_start; updated design; QC/Precursor_Start_{barplot,boxplot,Histogram}.png
#         steps/03_parameters/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "background.R"))
source(file.path(here, "..", "reviews", "01_load.R"))
source(file.path(here, "..", "reviews", "03_parameters.R"))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
use_app(c("Shiny_File.R", "Shiny_Data.R", "Shiny_Setup.R", "Shiny_Plots.R"))

params <- get_params(db_path)
run_dir <- dirname(db_path)

# --- 3a. Widget values (app: parameter_widget_save) ---------------------------
source(file.path(here, "..", "options.R")); source(file.path(here, "..", "options", "03_parameters.R"))
load_summ <- file.path(run_dir, "steps", "01_load", "summary.json")
detected <- if (file.exists(load_summ)) fromJSON(load_summ, simplifyVector = FALSE)$detected else NULL
opts <- parameters_options(params, read_table("design", db_path), detected)
spec <- opts_spec(opts)
later <- list.files(file.path(run_dir, "steps"), pattern = "^(0[4-9]|1[0-9])")
later <- later[file.exists(file.path(run_dir, "steps", later, "summary.json"))]
if (length(later) && !as_flag(args$force %||% FALSE))
  stop("Later steps were already run on this run: ", paste(later, collapse = ", "), ". Re-running Parameters changes their input, so ",
       "they would all need re-running. Start a new run, or re-run with --force and then redo those steps.", call. = FALSE)
args <- page_args(args, run_dir, "03_parameters", opts, "preview")
settings <- resolve_settings(spec, args)
all_opts <- opts_byname(opts)
for (k in names(all_opts)) if (!is.null(all_opts[[k]]$hide) || !is.null(all_opts[[k]]$disabled)) settings[[k]] <- spec[[k]]
opts_validate(opts, settings)
if (settings$data_output == "Protein" && settings$ptm) settings$ptm <- FALSE   # PTM analysis needs Peptide output (the app hides it)
if (settings$use_isoform) stop("use_isoform=TRUE is not supported yet in the pipeline (app's isoform_to_isoform() needs work).")
for (k in names(settings)) params[[k]] <- settings[[k]]
invisible(write_params(params, db_path))

# --- 3b. Sample groups (app: set_sample_groups) --------------------------------
# The app gives each group one block of columns (sample_groups start..end), so the design must list every group's samples
# together. Its check (check_design_sort) only asks that each sample has a neighbour in the same group, so a sheet like
# A_1 A_1 B_1 B_1 A_2 A_2 passes with primary_group=Primary although group A is split. Then Filter and Impute treat the wrong
# samples as group A. Fix (pipeline): when the groups used (Primary or Full) or the full groups are not each in one block, sort
# the design by group in sheet order (primary group, then full group, then sheet order; SPQC last). The columns follow
# (order_rename_columns orders them by design$Raw_Order); the sample-to-ID match is unchanged.
design <- read_table("design", db_path)
gkey <- if (settings$primary_group == "Primary") sub("_.*$", "", design$Group) else design$Group
one_block <- function(k) !anyDuplicated(rle(k)$values)
design_sorted <- FALSE
if (!one_block(gkey) || !one_block(design$Group)) {
  ord <- order(gkey == "SPQC", match(gkey, unique(gkey)), match(design$Group, unique(design$Group)), seq_len(nrow(design)))
  design <- design[ord, ]; rownames(design) <- NULL
  write_table("design", design, db_path); design_sorted <- TRUE
  cat("NOTE: the sample list's groups were not each in one block; samples sorted by group (sheet order kept within a group).\n")
}
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
# Background export: precursor_start keeps measured values only; background values wait in precursor_start_bg until Step 5.
bg_info <- if (is_background_run(db_path)) bg_split(db_path) else NULL
if (!is.null(bg_info)) cat(sprintf("Background export: %d background values set aside (measured only from here to the filter)\n", bg_info$background_values))

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
if (design_sorted) warnings_out <- c(warnings_out, paste("The sample list's groups were not each in one block, so the samples were sorted by group",
  "(sheet order kept within a group). Sample IDs still match; plots and tables list samples in the new order."))
off <- per_sample$sample[abs(log2(per_sample$total_vs_median)) > log2(1.5)]
if (length(off)) warnings_out <- c(warnings_out, sprintf("Total intensity differs >1.5x from the median sample: %s", paste(off, collapse = ", ")))

sd <- step_dir(run_dir, "03_parameters")
plots <- file.path(qc_path, c("Precursor_Start_barplot.png", "Precursor_Start_boxplot.png", "Precursor_Start_Histogram.png"))
plots <- gsub("/+", "/", plots)
summary <- list(
  step = "03_parameters", action = "preview", db_path = db_path,   # "preview": accept.R accepts exactly these settings
  settings = settings, detected_at_load = detected,
  info_columns = n_info, info_column_names = names(df)[seq_len(n_info)],
  counts_start = list(precursors = params$meta_precursor_start, peptides = params$meta_peptide_start,
                      proteins = params$meta_protein_start),
  sample_groups = groups[c("Group", "Count", "start", "end")],
  per_sample = per_sample,
  intensity = list(log2_mean = params$intensity_mean, log2_sd = params$intensity_sd,
                   intensity_cutoff = params$intensity_cutoff, intensity_cutoff_sd = params$intensity_cutoff_sd,
                   bottom_x = params$bottom_x, total_na = params$total_na, total_misaligned = params$total_misaligned),
  background = bg_info, design_sorted = design_sorted,
  sample_id_alert = sample_id_alert, sample_error = isTRUE(sample_error),
  plots = plots[file.exists(plots)],
  warnings = warnings_out
)
write_summary(summary, sd)
log_decision(run_dir, "03_parameters", settings, note = args$note %||% "", analyst = args$analyst)

cat("\n==== Parameters complete ====\n")
cat(sprintf("Start: %s precursors, %s peptides, %s proteins | %d info columns\n",
            params$meta_precursor_start, params$meta_peptide_start, params$meta_protein_start, n_info))
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
if (length(warnings_out)) cat("Warnings:\n", paste0("  - ", warnings_out, collapse = "\n"), "\n") else cat("No warnings.\n")
cat("Review page:", review_parameters(db_path), "\n>>> Waiting for the analyst's review (accept with pipeline/accept.R; changes: re-run this step).\n")
