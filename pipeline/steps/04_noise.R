# Step 4 — Noise   (HUMAN DECISION STEP)
# Same as the app's Noise tab. Two actions, run one after the other:
#
#   1) --action=inflection : "Calculate Noise Inflection Point" (noise_inflection_bg)
#        Rscript pipeline/steps/04_noise.R --db=<db> --action=inflection [--noise_baseline_value=1]
#      -> QC/Inflection_Point.png + noise data point counts. Show these to the user; the user decides.
#
#   2) --action=apply : "Apply Noise" (noise_widget_save + noise_remove_bg + filtered histogram)
#        Rscript pipeline/steps/04_noise.R --db=<db> --action=apply --noise_type=none|fixed|dynamic
#            [--noise_baseline_value=1] [--noise_keep=FALSE] --decided_by=user
#      --noise_keep=TRUE: noise values are removed here and put back in Step 5 for precursors that pass the filter.
#      (noise_keep_value is saved for the app but not used.)
#      -> table precursor_noise + QC/Precursor_NoiseFiltered_Histogram.png
#
# --decided_by=user is required for apply: this step must not run on Claude's choice alone.
# --analyst=<name> and --note="<why>" are saved with the decision.
# Both actions (re)build the review page <run_dir>/review/04_noise.html.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "04_noise.R"))

args <- parse_args()
need(args, c("db", "action"))
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_Noise.R", "Shiny_Plots.R"))
if (get_param("raw_data_format", db_path) != "precursor") stop("Noise step only applies to precursor data (as in the app).")
sd <- step_dir(run_dir, "04_noise")

pct <- function(a, b) round(as.numeric(a) / as.numeric(b) * 100, 2)

if (args$action == "inflection") {
  baseline <- as.numeric(args$noise_baseline_value %||% 1)
  param_update("noise_baseline_value", baseline, db_path)
  run_bg(noise_inflection_bg, list("precursor_start", db_path), "inflection", db_path)

  p <- get_params(db_path)
  plot <- gsub("/+", "/", file.path(p$qc_path, "Inflection_Point.png"))
  summary <- list(
    step = "04_noise", action = "inflection", db_path = db_path,
    noise_baseline_value = baseline,
    noise_inflection = p$noise_inflection,
    noise_inflection_log2 = round(log2(as.numeric(p$noise_inflection)), 2),
    total_data_points = p$noise_total,
    dynamic = list(count = p$dynamic_noise_count, percent = pct(p$dynamic_noise_count, p$noise_total)),
    fixed = list(count = p$custom_noise_count, percent = pct(p$custom_noise_count, p$noise_total)),
    plots = plot,
    decision_needed = "User chooses noise_type (none / fixed / dynamic), and optionally noise_keep + noise_keep_value."
  )
  write_summary(summary, sd)
  cat("\n==== Noise inflection ====\n")
  cat(sprintf("Inflection point: %s (log2 %.2f)\n", p$noise_inflection, summary$noise_inflection_log2))
  cat(sprintf("Total data points: %s\nDynamic noise points: %s (%.2f%%)\nFixed (baseline %s) noise points: %s (%.2f%%)\n",
              p$noise_total, p$dynamic_noise_count, summary$dynamic$percent, baseline, p$custom_noise_count, summary$fixed$percent))
  page <- review_noise(db_path)
  cat("Plot:", plot, "\nReview page:", page, "\n>>> Waiting for the user's decision.\n")

} else if (args$action == "apply") {
  if (!identical(args$decided_by, "user")) stop("Noise is a human decision step: pass --decided_by=user after the user has chosen.")
  need(args, "noise_type")
  if (!args$noise_type %in% c("none", "fixed", "dynamic")) stop("--noise_type must be none, fixed or dynamic")
  if (args$noise_type == "dynamic" && as.numeric(get_param("noise_inflection", db_path)) == 0)
    stop("Run --action=inflection first (dynamic uses the inflection point).")

  settings <- list(noise_type = args$noise_type,
                   noise_baseline_value = as.numeric(args$noise_baseline_value %||% get_param("noise_baseline_value", db_path)),
                   noise_keep = if (!is.null(args$noise_keep)) as_flag(args$noise_keep) else FALSE,
                   noise_keep_value = as.numeric(args$noise_keep_value %||% 1))
  # noise_keep (user's rule, 2026-11-04): remove ALL noise values here so the filter treats them as missing; step 5 then
  # puts the measured values back for every precursor that survives the filter (pipeline/noise_keep.R).
  # The app's own noise_keep branch in noise_remove_bg() is not used (it cleans the wrong rows), so it runs with noise_keep off.
  p <- get_params(db_path)
  for (k in names(settings)) p[[k]] <- settings[[k]]
  p$noise_keep <- FALSE
  invisible(write_params(p, db_path))

  run_bg(noise_remove_bg, list("precursor_start", "precursor_noise", db_path), "noise", db_path)
  p <- get_params(db_path); p$noise_keep <- settings$noise_keep; invisible(write_params(p, db_path))
  run_bg(histogram_plot, list("precursor_noise", "Precursor_NoiseFiltered_Histogram", db_path), "histogram", db_path)

  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  before <- read_table("precursor_start", db_path); after <- read_table("precursor_noise", db_path)
  cols <- (ncol(before) - n_s + 1):ncol(before)
  na_before <- colSums(is.na(before[cols])); na_after <- colSums(is.na(after[cols]))
  removed <- setNames(as.integer(na_after - na_before), names(before)[cols])
  plot <- gsub("/+", "/", file.path(p$qc_path, "Precursor_NoiseFiltered_Histogram.png"))

  summary <- list(
    step = "04_noise", action = "apply", db_path = db_path, settings = settings,
    baseline_used = switch(settings$noise_type, dynamic = p$noise_inflection, fixed = settings$noise_baseline_value, none = 0),
    values_set_to_na = sum(removed), values_set_to_na_per_sample = as.list(removed),
    rows = nrow(after),
    intensity = list(log2_mean = p$intensity_mean, log2_sd = p$intensity_sd, intensity_cutoff = p$intensity_cutoff,
                     total_na = p$total_na, total_misaligned = p$total_misaligned),
    noise_keep = if (settings$noise_keep) "on: removed values are put back in Step 5 for precursors that pass the filter" else "off",
    plots = plot
  )
  write_summary(summary, sd)
  log_decision(run_dir, "04_noise", settings, note = args$note %||% "decided by user", analyst = args$analyst)
  cat("\n==== Noise applied ====\n")
  cat(sprintf("noise_type=%s  baseline used=%s  values set to NA: %d\n", settings$noise_type, summary$baseline_used, sum(removed)))
  page <- review_noise(db_path)
  cat("Plot:", plot, "\nReview page:", page, "\n")
} else stop("--action must be inflection or apply")
