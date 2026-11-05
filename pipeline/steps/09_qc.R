# Step 9 — QC
# Same as the app's QC tab plot buttons:
#   Spike plot   (qc_spike_plots_bg):   summed intensity of the spike accessions by design spike level (QC Spike Level).
#   Protein plot (qc_protein_plots_bg): one accession across all samples (one plot per accession).
# Run for each normalization in norm_type (e.g. sltmm and impute-only for comparison).
#
# Usage:
#   Rscript pipeline/steps/09_qc.R --db=<db> [--spike_accession=P02662,P02663] [--constant_accession=P00330,P00331]
#       [--protein_accession=<acc>,<acc>] [--norms=sltmm,impute]
#   spike_accession    = spiked at levels given in the design (lab: bovine casein)
#   constant_accession = spiked at the same level in every sample (lab: yeast ADH); plotted per sample, should be flat
#   protein_accession  = any other proteins to plot per sample
#
#   Rscript pipeline/steps/09_qc.R --db=<db> --action=accept --decided_by=user [--analyst=<name>] [--flag_samples=A_2,B_1] [--note="<why>"]
#     records the analyst's review (QC changes no data) and rebuilds the page.
#
# Output: QC/QC_Spike_<norm>_barplot.png, QC/<acc>_<norm>_barplot.png; steps/09_qc/summary.json; review/09_qc.html

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "09_qc.R"))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_QC.R", "Shiny_Plots.R"))

if (identical(args$action, "accept")) {
  if (!identical(args$decided_by, "user")) stop("accept needs --decided_by=user (after the analyst has reviewed the QC page).")
  flags <- if (is.null(args$flag_samples)) character(0) else trimws(strsplit(args$flag_samples, ",")[[1]])
  log_decision(run_dir, "09_qc", list(accepted = TRUE, flagged_samples = flags), note = args$note %||% "", analyst = args$analyst)
  page <- review_qc(db_path)
  cat("QC accepted", if (length(flags)) paste("- flagged:", paste(flags, collapse = ", ")) else "", "\nReview page:", page, "\n")
  quit(save = "no")
}

split_acc <- function(x) if (is.null(x) || identical(x, "")) character(0) else trimws(strsplit(x, ",")[[1]])
spike_acc <- split_acc(args$spike_accession %||% "P02662,P02663")
const_acc <- split_acc(args$constant_accession %||% "P00330,P00331")
prot_acc <- split_acc(args$protein_accession)

p <- get_params(db_path)
is_peptide <- p$data_output == "Peptide"
# Peptide output (e.g. phospho): the spike / constant / protein checks are protein-level, so they run on protein_qc_<norm>, the sum of each
# protein's peptides in peptide_impute_<norm> (QC only; nothing downstream uses it).
norms <- split_acc(args$norms %||% p$norm_type)
qc_path <- p$qc_path
n_s <- as.integer(p$sample_number)
design <- read_table("design", db_path); design <- design[order(design$Raw_Order), ]
level_col <- grep("Level", names(design), value = TRUE)[1]
levels <- as.numeric(design[[level_col]])
plots <- character(0)

results <- lapply(norms, function(n) {
  tbl <- str_c("protein_", n)
  if (is_peptide) {
    pt <- str_c("peptide_impute_", n)
    if (!table_exists(pt, db_path)) stop("Run Step 8 (Rollup) first: ", pt, " is missing.")
    pd <- read_table(pt, db_path); X <- as.matrix(pd[(ncol(pd) - n_s + 1):ncol(pd)])
    S <- rowsum(X, pd$Accession, reorder = FALSE, na.rm = TRUE); first <- match(rownames(S), pd$Accession)
    tbl <- str_c("protein_qc_", n)
    write_table_try(tbl, data.frame(pd[first, c("Accession", "Description", "Name", "Genes")], Precursors = as.integer(rowsum(pd$Precursors, pd$Accession, reorder = FALSE)),
                                    S, check.names = FALSE, row.names = NULL), db_path)
  }
  if (!table_exists(tbl, db_path)) stop("Run Step 8 (Rollup) first: ", tbl, " is missing.")
  df <- read_table(tbl, db_path)
  sc <- (ncol(df) - n_s + 1):ncol(df)
  out <- list(norm = n)

  # --- spike (levels from design) ---
  found <- intersect(spike_acc, df$Accession)
  if (length(found)) {
    title <- str_c("QC_Spike_", n)
    run_bg(qc_spike_plots_bg, list(tbl, title, paste(found, collapse = ","), db_path), str_c("qcspike_", n), db_path)
    plots <<- c(plots, existing_plots(qc_path, str_c(title, "_barplot.png")))
    inten <- colSums(df[df$Accession %in% found, sc])
    by_level <- aggregate(inten, by = list(level = levels), FUN = mean)
    names(by_level)[2] <- "mean_intensity"
    by_level$measured_ratio_to_lowest <- round(by_level$mean_intensity / by_level$mean_intensity[1], 2)
    by_level$expected_ratio_to_lowest <- round(by_level$level / by_level$level[1], 2)
    ref <- mean(inten[levels == min(levels)])
    per_sample <- data.frame(label = design$Label, level = levels, intensity = round(inten),
                             measured_ratio = round(inten / ref, 2), expected_ratio = round(levels / min(levels), 2))
    per_sample$error_pct <- round((per_sample$measured_ratio / per_sample$expected_ratio - 1) * 100, 1)
    out$spike <- list(accessions = found, by_level = by_level, per_sample = per_sample)
  } else out$spike <- list(accessions = character(0), note = "spike accessions not found")

  # --- constant spike + other proteins (per-sample plots) ---
  per_protein <- function(acc_list) lapply(intersect(acc_list, df$Accession), function(a) {
    title <- str_c(a, "_", n)
    run_bg(qc_protein_plots_bg, list(tbl, title, a, qc_path, db_path), str_c("qcprotein_", a, "_", n), db_path)
    plots <<- c(plots, existing_plots(qc_path, str_c(title, "_barplot.png")))
    v <- unlist(df[df$Accession == a, sc])
    list(accession = a, name = df$Name[df$Accession == a][1],
         values = setNames(as.list(round(v)), design$Label),
         cv_pct_all_samples = round(sd(v) / mean(v) * 100, 1),
         max_over_min = round(max(v) / min(v), 2))
  })
  out$constant <- per_protein(const_acc)
  out$proteins <- per_protein(prot_acc)
  out
})

summary <- list(step = "09_qc", db_path = db_path,
                settings = list(spike_accession = spike_acc, constant_accession = const_acc, protein_accession = prot_acc, norms = norms),
                design_spike_levels = setNames(as.list(levels), design$Label),
                results = results, plots = unique(plots))
sd <- step_dir(run_dir, "09_qc")
write_summary(summary, sd)
page <- review_qc(db_path)

cat("\n==== QC plots ====\n")
for (r in results) {
  cat(sprintf("[%s]\n", r$norm))
  if (length(r$spike$accessions)) { cat("  Spike by level:\n"); print(r$spike$by_level, row.names = FALSE) }
  for (c in r$constant) cat(sprintf("  constant %s (%s): CV %.1f%%, max/min %.2f\n", c$accession, c$name, c$cv_pct_all_samples, c$max_over_min))
  for (c in r$proteins) cat(sprintf("  protein  %s (%s): CV %.1f%%, max/min %.2f\n", c$accession, c$name, c$cv_pct_all_samples, c$max_over_min))
}
cat("Plots:\n", paste0("  ", summary$plots, collapse = "\n"), "\n")
cat("Review page:", page, "\n>>> Waiting for the analyst's review.\n")
