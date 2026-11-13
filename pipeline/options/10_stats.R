# Step 10 - Stats: the app's Stats Setup and Comparisons tabs (ui.R stats_setup, stats_compare). Used by steps/10_stats.R (defaults,
# checks) and reviews/10_stats.R (the Options panel, with the comparison builder).
# p: run params; design: design table (comparisons use the full group names, as the app); norms: the normalizations available
# (norm_type); limpa: TRUE after limpa in Impute and Rollup.
#
# Comparisons are one option, "comparisons": "<N groups>|<D groups>;..." with groups comma separated, e.g. "KO_30|WT_30;KO_37|WT_37"
# or "WT_30,WT_37|KO_30,KO_37" (several groups on one side are pooled, as the app's multi-select).
stats_groups <- function(design) unique(design$Group)
comps_encode <- function(comps) paste(vapply(comps, function(c) paste(paste(c$N, collapse = ","), paste(c$D, collapse = ","), sep = "|"), ""), collapse = ";")
comps_decode <- function(x) {
  x <- trimws(as.character(x %||% "")); if (!nzchar(x)) return(list())
  lapply(strsplit(x, ";", fixed = TRUE)[[1]], function(c) { s <- strsplit(c, "|", fixed = TRUE)[[1]]
    list(N = trimws(strsplit(s[1] %||% "", ",", fixed = TRUE)[[1]]), D = trimws(strsplit(s[2] %||% "", ",", fixed = TRUE)[[1]])) })
}
# The app's comparison name: the parts of the group names (split at "_"), each once, joined by "_"; then N_v_D.
comp_name <- function(N, D) { f <- function(g) paste(unique(unlist(strsplit(g, "_", fixed = TRUE))), collapse = "_"); paste0(f(N), "_v_", f(D)) }
# Default: the second group against the first (sheet order), pool (SPQC) left out.
default_comparisons <- function(groups, spqc = "SPQC") {
  g <- setdiff(groups, spqc)
  if (length(g) >= 2) comps_encode(list(list(N = g[2], D = g[1]))) else ""
}

stats_options <- function(p = NULL, design = NULL, norms = NULL, limpa = FALSE) {
  grp <- if (!is.null(design)) stats_groups(design) else NULL
  has_spqc <- !is.null(grp) && "SPQC" %in% grp
  counts <- if (!is.null(design)) as.list(table(factor(design$Group, levels = grp))) else NULL
  ptm_pep <- !is.null(p) && is_ptm_run(p) && is_peptide_output(p)
  list(
    opt_section("Comparisons", list(
      opt("comparisons", "Comparisons", "comparisons", default_comparisons(grp),
          paste("Each comparison: numerator group(s) <b>v</b> denominator group(s); fold change = numerator / denominator. Pick several",
                "groups on one side to pool them. <b>Swap</b> turns a comparison around. Long names are shortened in the Excel sheet names only (as the app does)."),
          choices = grp, data = list(counts = counts)),
      opt("comp_spqc", "SPQC Group?", "select", if (has_spqc) "SPQC" else "0",
          "Pooled QC carried along in the tables and plots; not part of the test. 0 = none.", choices = c("0", grp)),
      opt("stats_norm_type", "Normalization", "select", norms[1] %||% "",
          "Which normalized data the tests use.", choices = norms))),
    opt_section("Significance", list(
      opt("pvalue_cutoff", "Pvalue Cutoff", "number", 0.05, "p at or below this (with the fold change) is Up/Down.", min = 0, max = 1, step = "any"),
      opt("foldchange_cutoff", "Fold Change cutoff", "number", 1.5, "Absolute fold change needed for Up/Down.", min = 1, step = "any"),
      opt("missing_factor", "Measured % (decimal)", "number", 0.6,
          "Requires this share of measured (not imputed) values in one group; rows below it are removed from the comparison.", min = 0, max = 1, step = "any"),
      opt("checkbox_adjpval", "Include adjusted pvalue?", "checkbox", TRUE, "Adds the adjusted p column to the tables."),
      opt("padjust_options", "p.adjust method", "select", "fdr", "", choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
          depends = "checkbox_adjpval"),
      opt("checkbox_filter_adjpval", "Filter with adjusted pvalue?", "checkbox", FALSE,
          "Up/Down calls use the adjusted p (controls false discoveries) instead of the raw p (more calls, no FDR control).",
          depends = "checkbox_adjpval"))),
    opt_section("Peptide refilter (flags)", note = "Stats doesn't re-roll per comparison (user decision 2026-11-04); these set the flag rules instead.", list(
      opt("peptide_refilter", "Peptide refilter", "checkbox", FALSE, "", hide = "retired: replaced by flags"),
      opt("peptide_missing_filter", "Refilter precursors/peptides by requiring X% measured values in one group?", "checkbox", FALSE,
          "Flags rows whose precursors are measured in fewer than X of one group's samples."),
      opt("peptide_missing_factor", "Peptide X% measured cutoff (decimal)", "number", 0.1, "", depends = "peptide_missing_filter", min = 0, max = 1, step = "any"),
      opt("peptide_cv_filter", "Refilter precursors/peptides by requiring X %CV one group?", "checkbox", FALSE, "Flag rule on the precursors' CV."),
      opt("peptide_cv_factor", "Peptide X CV% cutoff (decimal)", "number", 0.5, "", depends = "peptide_cv_filter", min = 0, step = "any"))),
    opt_section("Filters", list(
      opt("stats_spqc_cv_filter", "Filter by SPQC CV", "checkbox", FALSE, "Removes rows whose pool CV is above the cutoff.",
          disabled = if (!is.null(grp) && !has_spqc) "No SPQC group in this run."),
      opt("stats_spqc_cv_filter_factor", "SPQC %CV Cutoff", "number", 50, "", depends = "stats_spqc_cv_filter", min = 0),
      opt("stats_comp_cv_filter", "Require one group CV below cuttoff", "checkbox", FALSE, "Removes rows where no group of the comparison has a CV below the cutoff."),
      opt("stats_comp_cv_filter_factor", "Comp %CV Cutoff", "number", 50, "", depends = "stats_comp_cv_filter", min = 0),
      opt("stats_peptide_minimum", "Require minimum # of peptides per protein", "checkbox", FALSE, "Removes proteins with fewer peptides."),
      opt("stats_peptide_minimum_factor", "Peptide Minimum", "number", 1, "", depends = "stats_peptide_minimum", min = 1, step = 1))),
    opt_section("Extra columns", list(
      opt("pair_comp", "Pairwise Comparisons", "checkbox", FALSE,
          paste("Paired t-test: the samples of the two sides are paired by their order in the sample sheet (1st with 1st, ...), so both",
                "sides need the same number of samples; fold change = mean of the per-pair ratios. Cohen&apos;s d and the limma p stay unpaired (as in the app)."),
          disabled = if (isTRUE(limpa)) "Not built for limpa (dpcDE) yet."),
      opt("checkbox_cohensd", "Include Cohen's D?", "checkbox", FALSE, "Effect size (R package effsize)."),
      opt("checkbox_cohensd_hedges", "Use Hedge's Correction (low N)?", "checkbox", FALSE, "Correction for small sample size.", depends = "checkbox_cohensd"),
      opt("checkbox_limmapvalue", "Include Limma Pvalue?", "checkbox", FALSE, "Adds limma's p (R package limma)."))),
    opt_section("Report", list(
      opt("checkbox_report_ptm", "Report Only PTM?", "checkbox", ptm_pep,
          "Tests and reports only the peptides whose Sequence matches the grep; adjusted p covers those only (the app adjusts over all, then hides the rest).",
          disabled = if (!is.null(p) && !is_peptide_output(p)) "Protein output."),
      opt("peptide_report_grep", "Report PTM grep", "text", if (ptm_pep) p$ptm_grep else "Enter value here", "", depends = "checkbox_report_ptm"),
      opt("checkbox_report_accession", "Report Specific Accession(s) Only", "checkbox", FALSE, ""),
      opt("report_accession", "Protein Accessions for Final Report", "text", "Enter value", "", depends = "checkbox_report_accession")))
  )
}

# Check comparisons against the groups: each side at least one group, no group on both sides, no duplicates. (Long names are fine:
# the app's Excel writer cuts sheet names to 31 characters, and each starts with "Table<n> ", so they stay unique.)
# pair_comp: each comparison needs the same number of samples on both sides (the app pairs them by order and fails otherwise).
comps_pair_check <- function(comps, design) {
  n <- function(g) sum(design$Group %in% g)
  bad <- Filter(function(c) n(c$N) != n(c$D), comps)
  if (length(bad)) stop("Pairwise comparisons need the same number of samples on both sides: ",
    paste(vapply(bad, function(c) sprintf("%s (%d v %d)", comp_name(c$N, c$D), n(c$N), n(c$D)), ""), collapse = ", "),
    ". Turn Pairwise Comparisons off, or change these comparisons.", call. = FALSE)
  invisible(TRUE)
}
# The sample pairs of a paired comparison, in the app's order (design order within each side).
comps_pairs <- function(c, design) { a <- design$ID[design$Group %in% c$N]; b <- design$ID[design$Group %in% c$D]
  if (length(a) != length(b)) return("") ; paste(a, b, sep = " ~ ", collapse = ", ") }

comps_check <- function(comps, groups) {
  if (!length(comps)) stop("No comparison: add at least one (numerator v denominator).", call. = FALSE)
  for (c in comps) {
    if (!length(c$N) || !length(c$D) || any(!nzchar(c(c$N, c$D)))) stop("A comparison needs at least one group on each side.", call. = FALSE)
    bad <- setdiff(c(c$N, c$D), groups); if (length(bad)) stop("Unknown group(s): ", paste(bad, collapse = ", "), call. = FALSE)
    both <- intersect(c$N, c$D); if (length(both)) stop("Group on both sides of a comparison: ", paste(both, collapse = ", "), call. = FALSE)
  }
  nms <- vapply(comps, function(c) comp_name(c$N, c$D), ""); if (anyDuplicated(nms)) stop("The same comparison twice: ", nms[duplicated(nms)][1], call. = FALSE)
  invisible(TRUE)
}
