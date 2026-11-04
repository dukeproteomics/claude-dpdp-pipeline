# Step 5 — Filter
# Same as the app's Filter tab:
#   --action=preview : "Filter Cutoff" button. Saves settings, redraws the noise-filtered histogram
#                      (misaligned intensity cutoff = mean + intensity_cutoff_sd * SD).
#   --action=apply   : "Apply Filters" button. Saves settings, filter_data_bg, bar/box plots.
#                      If noise_keep is on, measured values removed as noise are put back for the precursors that
#                      pass the filter (pipeline/noise_keep.R), before the plots and missing-value tables.
#
# Usage:
#   Rscript pipeline/steps/05_filter.R --db=<db> --action=preview|apply [--settings_from=<shiny.db>] [--<setting>=<value> ...]
#       apply also needs --decided_by=user (after the analyst approves); optional --analyst=<name> --note="<why>"
# Both actions (re)build the review page <run_dir>/review/05_filter.html.
# Settings: app defaults, then values copied from --settings_from, then explicit --<setting> args.
#
# Output: precursor_filter, missing_values, missing_values_plots;
#         QC/Precursor_NoiseFiltered_Histogram.png, QC/Precursor_Filter_{barplot,boxplot}.png; steps/05_filter/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "noise_keep.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "05_filter.R"))

args <- parse_args()
need(args, c("db", "action"))
if (!args$action %in% c("preview", "apply")) stop("--action must be preview or apply")
db_path <- normalizePath(args$db, mustWork = TRUE)
stop_if_background(db_path, "05_filter")
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_Filter.R", "Shiny_Plots.R"))
if (get_param("raw_data_format", db_path) != "precursor") stop("Filter step is only wired up for precursor data so far.")
if (!table_exists("precursor_noise", db_path)) stop("Run Step 4 (Noise) first: precursor_noise is missing.")

# Widgets saved by the app's filter_widget_save(), with the app defaults.
spec <- list(
  filter_min_measured_all = 2, filter_x_percent = TRUE, filter_x_percent_value = 50,
  filter_cv = FALSE, filter_cv_group = "SPQC", filter_cv_value = 99, filter_ptm = FALSE,
  checkbox_misaligned = FALSE, misaligned_cutoff = 50, misaligned_target = "dataset", intensity_cutoff_sd = 0.5,
  precursor_quality = FALSE, precursor_quality_sd = 50, precursor_quality_intensity = 500, precursor_quality_min = 300,
  precursor_spqc_ratio = FALSE, precursor_spqc_intensity = 500, precursor_spqc_accuracy = 50
)
settings <- resolve_settings(spec, args)
if (!settings$misaligned_target %in% c("group", "dataset")) stop("--misaligned_target must be group or dataset")

p <- get_params(db_path)
groups <- read_table("sample_groups", db_path)
if (settings$filter_cv && !settings$filter_cv_group %in% groups$Group)
  stop("--filter_cv_group must be one of: ", paste(groups$Group, collapse = ", "))
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))

# Histogram (also recalculates intensity_cutoff, total_na, total_misaligned) — app does this on "Filter Cutoff".
run_bg(histogram_plot, list("precursor_noise", "Precursor_NoiseFiltered_Histogram", db_path), "histogram", db_path)
p <- get_params(db_path)
qc_path <- p$qc_path
sd <- step_dir(run_dir, "05_filter")
cutoff_info <- list(intensity_cutoff = p$intensity_cutoff, intensity_cutoff_log2 = round(log2(as.numeric(p$intensity_cutoff)), 2),
                    intensity_cutoff_sd = p$intensity_cutoff_sd, log2_mean = p$intensity_mean, log2_sd = p$intensity_sd,
                    total_na_by_group = setNames(as.list(as.numeric(strsplit(p$total_na, ",\\s*")[[1]])), groups$Group),
                    total_misaligned_by_group = setNames(as.list(as.numeric(strsplit(p$total_misaligned, ",\\s*")[[1]])), groups$Group))

if (args$action == "preview") {
  summary <- list(step = "05_filter", action = "preview", db_path = db_path, settings = settings,
                  cutoff = cutoff_info, plots = existing_plots(qc_path, "Precursor_NoiseFiltered_Histogram.png"),
                  decision_needed = "User reviews the histogram and settings, then approves --action=apply (or changes settings).")
  write_summary(summary, sd)
  cat("\n==== Filter preview ====\n")
  cat(sprintf("Misaligned intensity cutoff: %s (log2 %.2f, mean + %s SD)\n", p$intensity_cutoff, cutoff_info$intensity_cutoff_log2, p$intensity_cutoff_sd))
  page <- review_filter(db_path)
  cat("Plot:", summary$plots, "\nReview page:", page, "\n>>> Waiting for the user's approval.\n")
  quit(save = "no")
}

# --- apply -------------------------------------------------------------------
if (!identical(args$decided_by, "user")) stop("Apply needs the analyst's approval: run --action=preview, show the review page, then pass --decided_by=user.")
if (settings$checkbox_misaligned && settings$misaligned_target == "dataset" && sum(unlist(cutoff_info$total_misaligned_by_group)) == 0)
  stop("No misaligned rows expected; the app's filter would drop ALL rows in this case (df[-integer(0),]). Turn the misaligned filter off or use target=group.")

err_alert <- run_bg(filter_data_bg, list("precursor_noise", "precursor_filter", db_path), "filter", db_path)
if (nrow(read_table("precursor_filter", db_path)) == 0)
  stop("The filter removed every precursor. Known app behavior: with no rows failing the SPQC ratio (or misaligned, target=dataset) ",
       "filter, df[-integer(0),] drops all rows. Turn that filter off or loosen it and re-run.")
noise_keep_info <- NULL
if (as_flag(get_param("noise_keep", db_path)) && get_param("noise_type", db_path) != "none") {
  noise_keep_info <- noise_keep_restore(db_path)
  cat(sprintf("noise_keep: put back %d measured values (in %d precursors) that the noise step removed\n",
              noise_keep_info$values_restored, noise_keep_info$precursors_with_restored_values))
}
run_bg(bar_plot, list("precursor_filter", "Precursor_Filter", qc_path, db_path), "filterbarplot", db_path)
run_bg(box_plot, list("precursor_filter", "Precursor_Filter", qc_path, db_path), "filterboxplot", db_path)

# Filter cascade from the app's log ("step N - proteins peptides precursors" removed at each stage).
log_lines <- readLines(file.path(p$error_path, "error_filter.txt"))
removed <- regmatches(log_lines, regexpr("^step [0-9]+ - -?[0-9]+ -?[0-9]+ -?[0-9]+", log_lines))
stage_names <- c("ptm_isolate", "min_measured_all", "min_group_percent", "misaligned", "reapply_minimum",
                 "cv", "precursor_quality", "spqc_ratio", "remove_duplicates")
cascade <- lapply(seq_along(removed), function(i) {
  v <- as.integer(strsplit(sub("^step [0-9]+ - ", "", removed[i]), " ")[[1]])
  list(stage = stage_names[i] %||% paste0("stage_", i), proteins_removed = v[1], peptides_removed = v[2], precursors_removed = v[3])
})
mis_line <- grep("Misaligned rows -->|Misaligned dataset rows -->", log_lines, value = TRUE)

p <- get_params(db_path)
df <- read_table("precursor_filter", db_path)
start <- read_table("precursor_noise", db_path)
n_s <- as.integer(p$sample_number); cols <- (ncol(df) - n_s + 1):ncol(df)
final <- list(proteins = length(unique(df$Accession)), peptides = length(unique(df$Sequence)), precursors = nrow(df))
mv <- read_table("missing_values", db_path)
miss <- setNames(mv$num.missing, mv$key)[names(df)[cols]]
miss[is.na(miss)] <- 0

warnings_out <- character(0)
if (!identical(err_alert, "") && !is.null(err_alert)) warnings_out <- c(warnings_out, err_alert)
if (as.numeric(p$meta_precursor_filter) != final$precursors)
  warnings_out <- c(warnings_out, sprintf("App's stored filtered counts (meta_*_filter: %s precursors) are taken before the SPQC ratio + duplicate steps; the true final count is %d.", p$meta_precursor_filter, final$precursors))

summary <- list(
  step = "05_filter", action = "apply", db_path = db_path, settings = settings, cutoff = cutoff_info,
  counts_before = list(proteins = length(unique(start$Accession)), peptides = length(unique(start$Sequence)), precursors = nrow(start)),
  cascade = cascade, misaligned_log = mis_line,
  counts_final = final,
  app_reported_counts = list(proteins = p$meta_protein_filter, peptides = p$meta_peptide_filter, precursors = p$meta_precursor_filter),
  noise_keep = noise_keep_info,
  missing_after_filter_per_sample = as.list(miss),
  missing_after_filter_pct = as.list(round(miss / nrow(df) * 100, 2)),
  plots = existing_plots(qc_path, c("Precursor_NoiseFiltered_Histogram.png", "Precursor_Filter_barplot.png", "Precursor_Filter_boxplot.png")),
  warnings = warnings_out
)
write_summary(summary, sd)
log_decision(run_dir, "05_filter", settings, analyst = args$analyst,
             note = paste(c(args$note, if (!is.null(args$settings_from)) paste("settings copied from", args$settings_from)), collapse = "; "))

cat("\n==== Filter applied ====\n")
cat(sprintf("Before: %d proteins, %d peptides, %d precursors\n", summary$counts_before$proteins, summary$counts_before$peptides, summary$counts_before$precursors))
for (s in cascade) cat(sprintf("  %-18s -%d proteins  -%d peptides  -%d precursors\n", s$stage, s$proteins_removed, s$peptides_removed, s$precursors_removed))
cat(sprintf("After:  %d proteins, %d peptides, %d precursors\n", final$proteins, final$peptides, final$precursors))
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
if (length(warnings_out)) cat("Warnings:\n", paste0("  - ", warnings_out, collapse = "\n"), "\n")
page <- review_filter(db_path)
cat("Review page:", page, "\n")
