# Step 8 — Rollup
# Same as the app's Rollup tab "Apply Rollup" (server.R observeEvent(input$rollup_apply)), for Protein output:
#   rollup_widget_save -> rollup_apply_bg (precursor -> protein with the chosen method, precursor -> peptide by sum)
#   -> qc_stats_bg (group CVs, *_cv and *_final tables, summary_cv) -> cv_grouped_plot_bg -> norm_comp_plot_bg
#
# Usage:
#   Rscript pipeline/steps/08_rollup.R --db=<db> --action=preview [--settings_from=<db>] [--rollup_method=sum|iq_maxlfq] [--maxlfq_scale=TRUE]
#   Rscript pipeline/steps/08_rollup.R --db=<db> --action=apply   [...same...] --decided_by=user [--analyst=<name>] [--note="<why>"]
# preview: saves the settings and builds the review page (sum vs iq_maxlfq comparison); apply: runs the rollup.
# Only sum and iq_maxlfq work in the app: median, mean, median_polish and topn fail in Shiny_Rollup_Functions.R.
# limpa (pipeline/limpa.R): only after --impute_type=limpa. DPC-Quant on the gapped precursor tables; preview runs it and saves the result.
#
# Output: protein_<norm>, peptide_<norm>, protein_<norm>_cv, protein_<norm>_final, summary_cv;
#         QC/CV_barplot.png, QC/protein_<norm>_barplot.png; steps/08_rollup/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "background.R"))
source(file.path(here, "..", "limpa.R"))
source(file.path(here, "..", "reviews", "08_rollup.R"))

args <- parse_args()
need(args, c("db", "action"))
if (!args$action %in% c("preview", "apply")) stop("--action must be preview or apply")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_Rollup.R", "Shiny_QC.R", "Shiny_Plots.R"))

spec <- list(rollup_method = "sum", rollup_topn = 3, maxlfq_scale = FALSE, limpa_background = "missing")
settings <- resolve_settings(spec, args)
methods <- c("sum", "median", "median_polish", "mean", "iq_maxlfq", "topn", "limpa")
if (!settings$rollup_method %in% methods) stop("--rollup_method must be one of: ", paste(methods, collapse = ", "))
is_limpa <- identical(get_param("impute_type", db_path), "limpa")
if (is_limpa && settings$rollup_method != "limpa") stop("Impute was limpa (no imputation), so the rollup must be --rollup_method=limpa.")
if (!is_limpa && settings$rollup_method == "limpa") stop("--rollup_method=limpa needs --impute_type=limpa in Step 7 (the gaps must be kept).")
if (!settings$rollup_method %in% c(ROLLUP_WORKING, "limpa"))
  stop("--rollup_method=", settings$rollup_method, ": the app's function for this method fails (tested on 11488). Use sum or iq_maxlfq.")
if (!settings$rollup_topn %in% 1:5) stop("--rollup_topn must be 1-5")
if (!settings$limpa_background %in% c("missing", "use")) stop("--limpa_background must be missing or use")

p <- get_params(db_path)
if (p$data_output != "Protein") stop("Rollup tab only applies to Protein output (the app rolls up Peptide output during Impute).")
norms <- trimws(strsplit(p$norm_type, ",")[[1]])
for (n in norms) if (!table_exists(str_c("precursor_impute_", n), db_path)) stop("Run Step 7 (Impute) first: precursor_impute_", n, " is missing.")
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))
qc_path <- p$qc_path
n_s <- as.integer(p$sample_number)

limpa_info <- NULL
if (is_limpa) {   # DPC-Quant at preview (the page shows it) and again at apply, so apply always uses the current tables
  limpa_info <- lapply(setNames(norms, norms), function(n) {
    t1 <- Sys.time(); obj <- limpa_quant(db_path, n, background = settings$limpa_background)
    cat(sprintf("limpa DPC-Quant (%s): %d proteins, %.0f s; DPC beta0 %.3f beta1 %.3f\n", n, nrow(obj$Q$E),
                as.numeric(difftime(Sys.time(), t1, units = "secs")), obj$dpc[1], obj$dpc[2]))
    list(norm = n, proteins = nrow(obj$Q$E), background = obj$background, dpc = as.list(round(obj$dpc, 4)),
         se_median = round(median(obj$Q$other$standard.error), 3), se_q90 = round(unname(quantile(obj$Q$other$standard.error, 0.9)), 3))
  })
}
if (args$action == "preview") {
  sd <- step_dir(run_dir, "08_rollup")
  write_summary(list(step = "08_rollup", action = "preview", db_path = db_path, settings = settings, norm_type = norms, limpa = limpa_info,
                     decision_needed = "Analyst picks rollup_method (sum or iq_maxlfq) and maxlfq_scale from the review page."), sd)
  page <- review_rollup(db_path)
  cat("Review page:", page, "\n>>> Waiting for the user's decision.\n")
  quit(save = "no")
}
if (!identical(args$decided_by, "user")) stop("Apply needs the analyst's approval: run --action=preview, show the review page, then pass --decided_by=user.")
if (is_limpa) { for (n in norms) limpa_write_tables(db_path, n, readRDS(limpa_rds(run_dir, n)))
} else run_bg(rollup_apply_bg, list(db_path), "rollup_bg", db_path)
run_bg(qc_stats_bg, list(db_path), "qc", db_path)
run_bg(cv_grouped_plot_bg, list(db_path), "cv_plot", db_path)
run_bg(norm_comp_plot_bg, list(db_path), "norm_comp_plot", db_path)

# --- Summary ---
groups <- read_table("sample_groups", db_path)
cv <- read_table("summary_cv", db_path)
tables <- lapply(norms, function(n) {
  prot <- read_table(str_c("protein_", n), db_path); pep <- read_table(str_c("peptide_", n), db_path)
  pcv <- read_table(str_c("protein_", n, "_cv"), db_path)
  tot <- colSums(prot[(ncol(prot) - n_s + 1):ncol(prot)], na.rm = TRUE)
  cv_cols <- grep("_CV$", names(pcv), value = TRUE)
  list(norm = n, proteins = nrow(prot), peptides = nrow(pep),
       protein_total_cv_pct = round(sd(tot) / mean(tot) * 100, 2),
       median_cv_by_group = as.list(round(sapply(pcv[cv_cols], median, na.rm = TRUE), 1)),
       proteins_cv_under_20_pct = as.list(round(sapply(pcv[cv_cols], function(v) mean(v < 20, na.rm = TRUE) * 100), 1)))
})
summary <- list(
  step = "08_rollup", action = "apply", db_path = db_path, settings = settings, norm_type = norms,
  summary_cv = cv, tables = tables, limpa = limpa_info,
  plots = existing_plots(qc_path, c("CV_barplot.png", str_c("protein_", norms, "_barplot.png")))
)
sd <- step_dir(run_dir, "08_rollup")
write_summary(summary, sd)
log_decision(run_dir, "08_rollup", settings, analyst = args$analyst,
             note = paste(c(args$note, if (!is.null(args$settings_from)) paste("settings copied from", args$settings_from)), collapse = "; "))

cat("\n==== Rollup applied ====\n")
cat(sprintf("Method: %s (topn %s, maxlfq_scale %s)\n", settings$rollup_method, settings$rollup_topn, settings$maxlfq_scale))
for (t in tables) cat(sprintf("  %-7s %d proteins, %d peptides, protein total CV %.2f%%\n", t$norm, t$proteins, t$peptides, t$protein_total_cv_pct))
cat("Average CV% by group (summary_cv):\n"); print(cv, row.names = FALSE)
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
page <- review_rollup(db_path)
cat("Review page:", page, "\n")
