# Step 5 - Filter: every option of the app's Filter tab (ui.R, filter_widget_save). Precursor quality uses the pipeline's profile method.
# Used by steps/05_filter.R (defaults, checks) and reviews/05_filter.R (the Options panel). Order and labels follow the app.
# p: run params (NULL = no run context), groups: sample_groups table.
filter_options <- function(p = NULL, groups = NULL) {
  pep <- !is.null(p) && is_peptide_output(p)
  ptm_run <- !is.null(p) && as_flag(p$ptm %||% FALSE)
  grp <- if (!is.null(groups)) groups$Group else NULL
  n_s <- if (!is.null(p)) as.integer(p$sample_number) else NULL
  has_spqc <- is.null(grp) || any(grepl("SPQC", grp, ignore.case = TRUE))
  cut_info <- if (!is.null(p) && !is.null(p$intensity_cutoff) && !is.na(suppressWarnings(as.numeric(p$intensity_cutoff))))
    sprintf("At the last preview: mean + %s SD = %s.", p$intensity_cutoff_sd %||% "?", format(round(as.numeric(p$intensity_cutoff)), big.mark = ","))
  list(
    opt_section("Minimum measured", list(
      opt("filter_min_measured_all", "Minimum # measured values (all samples)", "number", 2,
          "Removes precursors measured in fewer samples than this, across all samples.", min = 0, max = n_s, step = 1))),
    opt_section("Measured in X% of a group", list(
      opt("filter_x_percent", "Require X% measured values in at least one group", "checkbox", TRUE,
          "Keeps a precursor only if at least one group has it in X% of its samples."),
      opt("filter_x_percent_value", "X% measured values", "number", 50, "Percent of a group's samples.", depends = "filter_x_percent",
          min = 0, max = 100, step = 1))),
    opt_section("CV of one group", list(
      opt("filter_cv", "Filter on group CV", "checkbox", FALSE, "Removes precursors that vary too much within the chosen group (usually SPQC)."),
      opt("filter_cv_group", "Group", "select", if (!is.null(grp) && !"SPQC" %in% grp) grp[1] else "SPQC", "The group whose CV is used.",
          choices = grp, depends = "filter_cv"),
      opt("filter_cv_value", "Cutoff %", "number", 99, "Precursors with a CV at or above this are removed.", depends = "filter_cv", min = 0))),
    opt_section("PTM only", list(
      opt("filter_ptm", "Keep only PTM", "checkbox", FALSE,
          paste0("Keeps only sequences matching the PTM pattern", if (!is.null(p$ptm_grep)) paste0(" (<b>", esc(p$ptm_grep), "</b>)"), "."),
          disabled = if (!is.null(p) && !ptm_run) "Not a PTM run (Parameters)."))),
    opt_section("Misaligned", note = "On by default with target dataset (triple-spike test, 2026-11-06).", list(
      opt("checkbox_misaligned", "Misaligned filter", "checkbox", TRUE,
          "A precursor that is high where it is measured but missing in much of a group is likely misaligned."),
      opt("intensity_cutoff_sd", "Misaligned intensity cutoff = mean + (x * SD)", "number", 0.5,
          "x: how high (in SDs of all log2 values) a group's mean must be to count.", depends = "checkbox_misaligned", step = 0.05, info = cut_info),
      opt("misaligned_cutoff", "X% missing values to be considered misaligned if the average > intensity cutoff", "number", 50,
          "A group with more missing than this (and not all missing) is tested.", depends = "checkbox_misaligned", min = 0, max = 100, step = 1),
      opt("misaligned_target", "Remove misalignment from group or dataset", "select", "dataset",
          "dataset removes the precursor; group blanks that group, and limpa (Step 7) is then not available (blanked groups read as low).",
          choices = c("dataset", "group"), depends = "checkbox_misaligned"))),
    opt_section("Precursor quality", list(
      opt("precursor_quality", "Precursor quality filter", "checkbox", FALSE,
          paste(if (pep) "Within each peptide (3+ charge states), removes charge states that don't follow the peptide's median-polish profile. Phosphopeptides of one protein are never compared with each other."
                else "Within each protein (3+ precursors), removes precursors that don't follow the protein's median-polish profile.",
                "Pipeline method (profile); the app's method was retired (2026-11-09): it removed good precursors of proteins that change.")),
      opt("precursor_quality_sd", "Percent SD cutoff", "number", 50, "Precursors whose residual CV % is above this are removed (worst first, at least 2 kept).",
          depends = "precursor_quality", min = 0),
      opt("precursor_quality_intensity", "Minimum average intensity", "number", 500, "Only precursors with this mean intensity or more are tested.",
          depends = "precursor_quality", min = 0),
      opt("precursor_quality_min", "Minimum intensity for consideration", "number", 300, "Values below this are ignored in the test.",
          depends = "precursor_quality", min = 0))),
    opt_section("SPQC ratio", list(
      opt("precursor_spqc_ratio", "Precursor SPQC accuracy filter", "checkbox", FALSE,
          "The pool should sit near the sample average; precursors where it doesn't are removed.",
          disabled = if (!has_spqc) "No SPQC group in this run."),
      opt("precursor_spqc_accuracy", "+- percent SPQC accuracy", "number", 50, "Allowed distance of SPQC mean / sample mean from 100%.",
          depends = "precursor_spqc_ratio", min = 0),
      opt("precursor_spqc_intensity", "Minimum intensity of SPQC samples", "number", 500, "Only precursors whose SPQC mean is above this are tested.",
          depends = "precursor_spqc_ratio", min = 0)))
  )
}
