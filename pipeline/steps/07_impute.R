# Step 7 — Impute
# Same as the app's Impute tab:
#   --action=preview : "Set Impute Parameters" (impute_apply_widget_save + impute_meta_data_bg + missing plots).
#                      Builds the impute_bin_<group> tables and the seeded `random` table (set.seed(123)).
#   --action=apply   : "Apply Imputation" (impute_apply_bg for every method in norm_type + create_imputed_df).
#                      For sltmm the app applies TMM here, after imputation.
#
# Usage:
#   Rscript pipeline/steps/07_impute.R --db=<db> --action=preview|apply [--settings_from=<db>]
#       [--impute_type=duke|bottomx] [--bottom_x=2] [--missing_cutoff=50]
#   apply also needs --decided_by=user (after the analyst approves); optional --analyst=<name> --note="<why>"
# Both actions (re)build the review page <run_dir>/review/07_impute.html.
#
# Output: impute_bin_<group>, random, precursor_impute_<norm> (each norm in norm_type), precursor_missing, protein_missing;
#         QC/missing_bar_plot.png, QC/missing_percent_plot.png, QC/Precursor_impute_<norm>_barplot.png (pipeline-added)

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "07_impute.R"))

args <- parse_args()
need(args, c("db", "action"))
if (!args$action %in% c("preview", "apply")) stop("--action must be preview or apply")
db_path <- normalizePath(args$db, mustWork = TRUE)
stop_if_background(db_path, "07_impute")
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_Data.R", "Shiny_Impute.R", "Shiny_MVA_Functions.R", "Shiny_Plots.R"))
if (get_param("raw_data_format", db_path) != "precursor") stop("Impute step is only wired up for precursor data so far.")

spec <- list(impute_type = "duke", bottom_x = 2, missing_cutoff = 50)
settings <- resolve_settings(spec, args)
# The app's dropdown lists more methods, but impute_apply_bg2() only implements these two.
if (!settings$impute_type %in% c("duke", "bottomx"))
  stop("--impute_type must be duke or bottomx (other app options are not implemented in impute_apply_bg2()).")

p <- get_params(db_path)
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))
qc_path <- p$qc_path
n_s <- as.integer(p$sample_number)
norms <- trimws(strsplit(p$norm_type, ",")[[1]])
for (n in norms) if (!table_exists(str_c("precursor_norm_", n), db_path)) stop("Run Step 6 (Normalize) first: precursor_norm_", n, " is missing.")
sd <- step_dir(run_dir, "07_impute")
groups <- read_table("sample_groups", db_path)

# Preview part always runs, so bins/random/plots match the current settings and filtered data.
run_bg(impute_meta_data_bg, list("precursor_filter", db_path), "impute_meta_data", db_path)
run_bg(missing_bar_plot, list(db_path), "missing_bar", db_path)
run_bg(missing_percent_plot, list(db_path), "missing_percent", db_path)

# Which method fills each missing value (same rule as impute_duke): per group, rows with
# 0 < missing <= Count * (100 - missing_cutoff)/100 -> duke (around the group mean); other missing -> bottom-X.
filt <- read_table("precursor_filter", db_path)
x <- filt[(ncol(filt) - n_s + 1):ncol(filt)]
fill_plan <- lapply(seq_len(nrow(groups)), function(i) {
  g <- x[groups$start[i]:groups$end[i]]
  miss <- rowSums(is.na(g)); max_miss <- groups$Count[i] * (100 - settings$missing_cutoff) / 100
  duke_rows <- miss > 0 & miss <= max_miss
  list(group = groups$Group[i], samples = groups$Count[i], missing_total = sum(miss),
       duke_values = if (settings$impute_type == "duke") sum(miss[duke_rows]) else 0,
       bottomx_values = if (settings$impute_type == "duke") sum(miss[!duke_rows]) else sum(miss),
       rows_missing_by_count = as.list(table(factor(miss[miss > 0], levels = seq_len(groups$Count[i])))))
})
bins <- setNames(lapply(groups$Group, function(g) read_table(str_c("impute_bin_", g), db_path)), groups$Group)
mv <- read_table("missing_values", db_path)
p <- get_params(db_path)

base_summary <- list(step = "07_impute", db_path = db_path, settings = settings, norm_type = norms,
                     missing_by_group = setNames(as.list(as.numeric(strsplit(p$meta_impute_na, ",\\s*")[[1]])), groups$Group),
                     fill_plan = fill_plan,
                     missing_per_sample = setNames(as.list(mv$num.missing), mv$key),
                     impute_bins = lapply(bins, function(b) b[c("bin", "min2", "max", "sd")]))

if (args$action == "preview") {
  summary <- c(base_summary, list(action = "preview",
               plots = existing_plots(qc_path, c("missing_bar_plot.png", "missing_percent_plot.png")),
               decision_needed = "User reviews the missing value plots and fill plan, then approves --action=apply."))
  write_summary(summary, sd)
  cat("\n==== Impute preview ====\n")
  for (f in fill_plan) cat(sprintf("  %-6s missing %6d  -> duke %6d, bottom-X %6d\n", f$group, f$missing_total, f$duke_values, f$bottomx_values))
  page <- review_impute(db_path)
  cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\nReview page:", page, "\n>>> Waiting for the user's approval.\n")
  quit(save = "no")
}

# --- apply -----------------------------------------------------------------------
if (!identical(args$decided_by, "user")) stop("Apply needs the analyst's approval: run --action=preview, show the review page, then pass --decided_by=user.")
run_bg(impute_apply_bg, list(as.list(norms), p, db_path), "impute_bg", db_path)
run_bg(create_imputed_df, list(p, db_path), "impute_df", db_path)
for (n in norms) run_bg(bar_plot, list(str_c("precursor_impute_", n), str_c("Precursor_impute_", n), qc_path, db_path),
                        str_c("imputebar_", n), db_path)

check <- lapply(norms, function(n) {
  d <- read_table(str_c("precursor_impute_", n), db_path); s <- d[(ncol(d) - n_s + 1):ncol(d)]
  tot <- colSums(s, na.rm = TRUE)
  list(table = str_c("precursor_impute_", n), rows = nrow(d), remaining_na = sum(is.na(s)),
       total_cv_pct = round(sd(tot) / mean(tot) * 100, 2))
})
summary <- c(base_summary, list(action = "apply", imputed_tables = check,
             note = if ("sltmm" %in% norms) "precursor_impute_sltmm has TMM applied after imputation (app behavior)." else NULL,
             plots = existing_plots(qc_path, c("missing_bar_plot.png", "missing_percent_plot.png",
                                               str_c("Precursor_impute_", norms, "_barplot.png"))),
             pipeline_added_plots = existing_plots(qc_path, str_c("Precursor_impute_", norms, "_barplot.png"))))
write_summary(summary, sd)
log_decision(run_dir, "07_impute", settings, analyst = args$analyst,
             note = paste(c(args$note, if (!is.null(args$settings_from)) paste("settings copied from", args$settings_from)), collapse = "; "))

cat("\n==== Imputation applied ====\n")
for (f in fill_plan) cat(sprintf("  %-6s missing %6d  -> duke %6d, bottom-X %6d\n", f$group, f$missing_total, f$duke_values, f$bottomx_values))
for (c in check) cat(sprintf("  %s: %d rows, %d NA left, sample total CV %.2f%%\n", c$table, c$rows, c$remaining_na, c$total_cv_pct))
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
page <- review_impute(db_path)
cat("Review page:", page, "\n")
