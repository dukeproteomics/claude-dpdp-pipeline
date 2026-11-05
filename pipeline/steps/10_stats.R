# Step 10 — Stats (Setup + Comparisons)
# Same as the app's Stats tabs:
#   Setup        "Save Stat Options" (save_stat_options)
#   Comparisons  "Check Comparisons" (stat_widget_save + check_comp_names / check_comp_names_bg) then
#                "Start Stats" (stat_calc_bg for each comparison)
#
# Usage:
#   Rscript pipeline/steps/10_stats.R --db=<db> --action=preview|apply|accept [--settings_from=<db>]
#       apply/accept need --decided_by=user; optional --analyst=<name> --note="<why>"
#       [--stats_norm_type=sltmm] [--comp_spqc=SPQC]
#       [--comp1_N=B --comp1_D=A [--comp2_N=... --comp2_D=...]]     (groups; several groups: comma separated)
#       [--pvalue_cutoff=0.05 --foldchange_cutoff=1.5 --checkbox_adjpval=TRUE --padjust_options=fdr ...]
# With --settings_from and no --compN_* args, the comparisons are copied from that run's stats_comp table.
#
# Output: stats_comp, protein_<norm>_<N>_v_<D>, protein_<norm>_<N>_v_<D>_final, peptide_impute_<norm>_<N>_v_<D>_final;
#         steps/10_stats/summary.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "contaminants.R"))
source(file.path(here, "..", "background.R"))
source(file.path(here, "..", "limpa.R"))
source(file.path(here, "..", "refilter_flags.R"))
source(file.path(here, "..", "reviews", "10_stats.R"))

args <- parse_args()
need(args, c("db", "action"))
if (!args$action %in% c("preview", "apply", "accept")) stop("--action must be preview, apply or accept")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_MVA.R", "Shiny_MVA_Functions.R"))
if (args$action == "accept") {
  if (!identical(args$decided_by, "user")) stop("accept needs --decided_by=user (after the analyst has reviewed the stats page).")
  log_decision(run_dir, "10_stats", list(accepted = TRUE), note = args$note %||% "", analyst = args$analyst)
  cat("Stats results accepted\nReview page:", review_stats(db_path), "\n"); quit(save = "no")
}
if (args$action == "apply" && !identical(args$decided_by, "user"))
  stop("Apply needs the analyst's approval: run --action=preview, show the review page, then pass --decided_by=user.")

# Stats Setup widgets (values = 11488 reference / app UI defaults)
spec <- list(pvalue_cutoff = 0.05, pair_comp = FALSE, checkbox_adjpval = TRUE, padjust_options = "fdr",
             foldchange_cutoff = 1.5, missing_factor = 0.6,
             peptide_refilter = TRUE, peptide_missing_filter = FALSE, peptide_missing_factor = 0.1,
             peptide_cv_filter = FALSE, peptide_cv_factor = 0.5,
             stats_spqc_cv_filter = FALSE, stats_spqc_cv_filter_factor = 50,
             stats_comp_cv_filter = FALSE, stats_comp_cv_filter_factor = 50,
             stats_peptide_minimum = FALSE, stats_peptide_minimum_factor = 1,
             checkbox_filter_adjpval = FALSE, checkbox_cohensd = FALSE, checkbox_cohensd_hedges = FALSE,
             checkbox_limmapvalue = FALSE,
             checkbox_report_ptm = FALSE, peptide_report_grep = "Enter value here",
             checkbox_report_accession = FALSE, report_accession = "Enter value")
is_limpa <- identical(applied_impute_type(db_path), "limpa") && identical(get_param("rollup_method", db_path), "limpa")
# Peptide refilter is replaced by flags (user decision 2026-11-04): Stats never re-rolls the comparison's samples; the proteins the
# refilter would have removed are flagged (pipeline/refilter_flags.R). peptide_missing_filter / peptide_cv_filter now set flag rules.
if (!is.null(args$peptide_refilter) && as_flag(args$peptide_refilter))
  stop("--peptide_refilter=TRUE is no longer used: Stats doesn't re-roll per comparison; proteins the refilter would remove are flagged.")
args$peptide_refilter <- "FALSE"
settings <- resolve_settings(spec, args)
settings$peptide_refilter <- FALSE
if (is_limpa && as_flag(settings$pair_comp)) stop("limpa run: paired comparisons (--pair_comp) are not built for dpcDE yet.")
padj <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
if (!settings$padjust_options %in% padj) stop("--padjust_options must be one of: ", paste(padj, collapse = ", "))

p <- get_params(db_path)
groups <- read_table("sample_groups", db_path)
norms <- trimws(strsplit(p$norm_type, ",")[[1]])

# Comparison setup: stats_norm_type, comp_spqc, comparisons
ref_comp <- NULL
if (!is.null(args$settings_from)) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), abs_path(args$settings_from))
  rp <- RSQLite::dbReadTable(con, "params")
  if (RSQLite::dbExistsTable(con, "stats_comp")) ref_comp <- RSQLite::dbReadTable(con, "stats_comp")
  RSQLite::dbDisconnect(con)
}
stat_norm <- args$stats_norm_type %||% (if (!is.null(args$settings_from)) rp$stat_norm else norms[1])
# Pooled QC group: SPQC when the design has one; otherwise none ("0"). The pipeline works without a pool.
comp_spqc <- args$comp_spqc %||% (if (!is.null(args$settings_from)) rp$comp_spqc else if ("SPQC" %in% groups$Group) "SPQC" else "0")
comp_args <- grep("^comp[0-9]+_[ND]$", names(args), value = TRUE)
if (length(comp_args)) {
  ids <- sort(unique(as.integer(sub("^comp([0-9]+)_.*", "\\1", comp_args))))
  comps <- lapply(ids, function(i) list(N = trimws(strsplit(args[[str_c("comp", i, "_N")]], ",")[[1]]),
                                        D = trimws(strsplit(args[[str_c("comp", i, "_D")]], ",")[[1]])))
} else if (!is.null(ref_comp)) {
  comps <- lapply(seq_len(nrow(ref_comp)), function(i) list(N = trimws(strsplit(ref_comp$FactorsN[i], ",")[[1]]),
                                                            D = trimws(strsplit(ref_comp$FactorsD[i], ",")[[1]])))
} else stop("Give comparisons with --comp1_N=<group> --comp1_D=<group> (or --settings_from=<db>).")

if (!stat_norm %in% norms) stop("--stats_norm_type must be one of: ", paste(norms, collapse = ", "))
if (comp_spqc != "0" && !comp_spqc %in% groups$Group) stop("--comp_spqc must be one of: ", paste(groups$Group, collapse = ", "), " (or 0)")
for (cp in comps) { bad <- setdiff(c(cp$N, cp$D), groups$Group); if (length(bad)) stop("Unknown group(s): ", paste(bad, collapse = ", ")) }

# --- Setup + widget save (save_stat_options + stat_widget_save) ---
for (k in names(settings)) p[[k]] <- settings[[k]]
# --- check_comp_names ---
table_name <- str_c("protein_", stat_norm)
table_name_peptide <- str_c("peptide_impute_", stat_norm)
p$stat_norm <- stat_norm
p$comp_spqc <- comp_spqc
p$comp_number <- length(comps)
invisible(write_params(p, db_path))
write_table_try("stats_comp", data.frame(matrix(ncol = 13, nrow = 0)), db_path)
p <- get_params(db_path)
for (i in seq_along(comps))
  run_bg(check_comp_names_bg, list(p, db_path, table_name, table_name_peptide, i, comps[[i]]$N, comps[[i]]$D),
         str_c("stat_groups_", i), db_path)
stats_comp <- read_table("stats_comp", db_path)
long <- stats_comp$Name[nchar(stats_comp$Name) > 31]
if (length(long)) stop("Comparison name(s) longer than 31 characters (Excel sheet limit): ", paste(long, collapse = ", "))
if (args$action == "preview") {
  sd <- step_dir(run_dir, "10_stats")
  write_summary(list(step = "10_stats", action = "preview", db_path = db_path, settings = settings, stat_norm = stat_norm,
                     comp_spqc = comp_spqc, stats_comp = stats_comp,
                     decision_needed = "Analyst approves the stats settings and comparisons, then --action=apply."), sd)
  cat("Comparisons:", paste(stats_comp$Name, collapse = ", "), "\nReview page:", review_stats(db_path), "\n>>> Waiting for the user's approval.\n")
  quit(save = "no")
}

# --- stat_calc ---
# Without refilter the pipeline runs the stats itself (pipeline/refilter_flags.R::stat_calc_norefilter_bg, the app's stat_add core),
# because the app's no-refilter branch reads the wrong columns and fails when a comparison doesn't use every sample.
for (i in seq_len(nrow(stats_comp)))
  run_bg(stat_calc_norefilter_bg, list(db_path, i, stats_comp), str_c("stat_calc_", i), db_path)
# limpa: the app built the tables (layout, CVs, missing factor); FC / pval / adjpval / Stats now come from dpcDE
limpa_test <- if (is_limpa) lapply(seq_len(nrow(stats_comp)), function(i) limpa_stats(db_path, stats_comp, i, get_params(db_path))) else NULL
refilter <- lapply(seq_len(nrow(stats_comp)), function(i) refilter_flags(db_path, stats_comp, i, get_params(db_path)))

# --- Summary ---
p <- get_params(db_path)
comp_results <- lapply(seq_len(nrow(stats_comp)), function(i) {
  nm <- stats_comp$Name[i]
  df <- read_table(stats_comp$Final_Table_Name[i], db_path)
  fc_col <- str_c(nm, "_FC"); p_col <- str_c(nm, "_pval"); adj_col <- str_c(nm, "_adjpval"); mf_col <- str_c(nm, "_mf")
  fc <- df[[fc_col]]; pv <- df[[p_col]]; adj <- if (adj_col %in% names(df)) df[[adj_col]] else rep(NA_real_, nrow(df))
  fcc <- as.numeric(p$foldchange_cutoff); pc <- as.numeric(p$pvalue_cutoff)
  # The app's own call (Stats column): raw p unless checkbox_filter_adjpval is on.
  up <- df$Stats == "Up"; down <- df$Stats == "Down"; sig <- up | down
  adj_sig <- !is.na(adj) & adj <= pc & abs(fc) >= fcc
  top <- df[sig, intersect(c("Accession", "Name", "Genes", fc_col, p_col, adj_col, mf_col, str_c(nm, "_refilter_flag"), "Stats"), names(df))]
  top <- head(top[order(top[[p_col]]), ], 15)
  pv_dist <- list(raw_p_lt_0.05 = sum(pv < 0.05, na.rm = TRUE), raw_p_lt_0.01 = sum(pv < 0.01, na.rm = TRUE),
                  raw_p_lt_0.001 = sum(pv < 0.001, na.rm = TRUE),
                  min_adjpval = if (all(is.na(adj))) NA else signif(min(adj, na.rm = TRUE), 3),
                  adjp_lt_0.05 = sum(adj < 0.05, na.rm = TRUE), adjp_lt_0.1 = sum(adj < 0.1, na.rm = TRUE),
                  adjp_lt_0.2 = sum(adj < 0.2, na.rm = TRUE))
  removed_mf <- nrow(read_table(stats_comp$Table_Name[i], db_path)) - nrow(df)
  sp <- qc_spike_settings(run_dir)
  flag <- flag_contaminants(df$Accession, df$Description, sp$spike, sp$constant)
  flagged_sig <- df[sig & flag != "", intersect(c("Accession", "Genes", fc_col, p_col, "Stats"), names(df))]
  flagged_sig$Flag <- flag[sig & flag != ""]
  list(comparison = nm, N = stats_comp$FactorsN[i], D = stats_comp$FactorsD[i],
       samples = list(N = as.integer(stats_comp$N[i]), D = as.integer(stats_comp$D[i]), SPQC = as.integer(stats_comp$SPQC[i])),
       table = stats_comp$Final_Table_Name[i], proteins = nrow(df),
       columns = names(df)[!names(df) %in% names(read_table(str_c("protein_", stat_norm), db_path))],
       removed_by_missing_factor = removed_mf,
       significant = sum(sig), up = sum(up), down = sum(down),
       pvalue_used_for_calls = if (as.logical(p$checkbox_filter_adjpval)) str_c("adjusted (", p$padjust_options, ")") else "raw (checkbox_filter_adjpval is off)",
       significant_if_adjusted_p = sum(adj_sig), pvalue_distribution = pv_dist,
       flagged = list(note = "Flagged, not removed (lab rule): digestion trypsin, human keratins (non-human samples), QC spike-ins.",
                      tested = as.list(table(flag[flag != ""])), significant = flagged_sig),
       refilter_flags = refilter[[i]], top_by_pvalue = top)
})
summary <- list(step = "10_stats", action = "apply", db_path = db_path, settings = settings,
                test = if (is_limpa) "limpa dpcDE: limma with DPC-Quant precision weights; FC from the model (log2 means)" else "app: t-test per protein; FC from the ratio of means",
                limpa = limpa_test, refilter_flags = refilter,
                stat_norm = stat_norm, comp_spqc = comp_spqc, stats_comp = stats_comp, comparisons = comp_results)
sd <- step_dir(run_dir, "10_stats")
write_summary(summary, sd)
log_decision(run_dir, "10_stats", analyst = args$analyst, settings = c(settings, list(stat_norm = stat_norm, comp_spqc = comp_spqc,
             comparisons = vapply(comps, function(cp) str_c(paste(cp$N, collapse = "_"), "_v_", paste(cp$D, collapse = "_")), ""))),
             note = paste(c(args$note, if (!is.null(args$settings_from)) paste("settings copied from", args$settings_from)), collapse = "; "))

cat("\n==== Stats calculated ====\n")
cat("Test:", summary$test, "\n")
cat(sprintf("Normalization: %s | SPQC group: %s | calls use %s p <= %s | FC >= %s | missing factor >= %s\n",
            stat_norm, comp_spqc, comp_results[[1]]$pvalue_used_for_calls, p$pvalue_cutoff, p$foldchange_cutoff, p$missing_factor))
for (r in comp_results) {
  cat(sprintf("  %s: %d proteins (%d removed by missing factor), %d significant (%d up, %d down) -> %s\n",
              r$comparison, r$proteins, r$removed_by_missing_factor, r$significant, r$up, r$down, r$table))
  cat(sprintf("     with adjusted p instead: %d significant (min adj p %s)\n", r$significant_if_adjusted_p, r$pvalue_distribution$min_adjpval))
  rf <- refilter[[which(stats_comp$Name == r$comparison)]]
  cat(sprintf("     refilter flags (not removed): all precursors %d (%d significant), half or more %d (%d significant)\n",
              rf$flagged_all, rf$flagged_all_significant, rf$flagged_most, rf$flagged_most_significant))
  if (nrow(r$flagged$significant)) cat(sprintf("     flagged among significant: %s\n", paste(sprintf("%s (%s)", r$flagged$significant$Genes, r$flagged$significant$Flag), collapse = ", ")))
}
cat("Review page:", review_stats(db_path), "\n")
