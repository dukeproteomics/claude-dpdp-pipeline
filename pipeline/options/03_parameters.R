# Step 3 - Parameters: the app's Parameters tab (ui.R param_box, parameter_widget_save). Used by steps/03_parameters.R
# (defaults, checks) and reviews/03_parameters.R (the Options panel).
# p: run params; design: the design table (Group column); detected: what Load found (steps/01_load/summary.json $detected).
# Options that never apply to the data source are hidden, as the app hides them (Shiny_Hide.R): peptide_select and use_isoform
# are Proteome Discoverer only; multi_tmt is TMT only (TMT is not supported by the pipeline).
parameters_options <- function(p = NULL, design = NULL, detected = NULL) {
  src <- detected$data_source %||% p$data_source %||% ""
  ptm_found <- isTRUE(as.logical(detected$ptm %||% p$ptm %||% FALSE))
  grp <- if (!is.null(design)) unique(design$Group) else NULL
  prim <- if (!is.null(grp)) unique(sub("_.*$", "", grp)) else NULL
  same <- !is.null(grp) && length(prim) == length(grp)
  grp_info <- if (!is.null(grp)) sprintf("Full: %d groups (%s). Primary: %d groups (%s).", length(grp), esc(paste(grp, collapse = ", ")),
                                         length(prim), esc(paste(prim, collapse = ", ")))
  list(
    opt_section("Groups", list(
      opt("primary_group", "Full or primary group for filter and impute", "select", "Full",
          paste("Full: every step uses the groups as named in the sample sheet. Primary: the Filter (X% of a group, misaligned, CV) and",
                "Impute use the part of the group name before the first <b>_</b>, so related groups (e.g. WT_30 and WT_37) are pooled;",
                "Stats comparisons are still built from the full group names."),
          choices = c("Full", "Primary"), info = grp_info,
          disabled = if (same) "No group name has a <b>_</b>, so Primary would be the same as Full."))),
    opt_section("Output", list(
      opt("data_output", "Select Output Data Type", "select", detected$data_output %||% "Protein",
          paste("Protein: precursors roll up to proteins (Step 8) and Stats tests proteins. Peptide: Stats tests peptides (charge states",
                "summed), for phospho/PTM sets or when peptide-level results are wanted in a set that is not enriched.",
                "The default is what Load detected."),
          choices = c("Protein", "Peptide")),
      opt("ptm", "Phos Analysis?", "checkbox", ptm_found,
          "PTM (phospho) analysis: keeps the PTM localization per site and reports PTM peptides. Needs Peptide output.",
          depends = "data_output=Peptide",
          disabled = if (!is.null(detected) && !ptm_found) "The report has no PTM localization column (EG.PTMProbabilities)."),
      opt("ptm_grep", "PTM grep", "text", "Phospho", "The modification name to look for in the sequences (e.g. Phospho).", depends = "ptm"),
      opt("ptm_local", "PTM localization min", "number", 0.75, "Minimum site localization probability for a site to count as localized.",
          depends = "ptm", min = 0, max = 1, step = 0.05))),
    opt_section("Other data sources", list(
      opt("multi_tmt", "SPQC Normalized TMT sets", "checkbox", FALSE, "Several TMT sets normalized to their SPQC.",
          hide = "TMT only"),
      opt("peptide_select", "Peptides to Use", "select", "Razor", "Which peptides count for a protein: Razor, Unique or Shared.",
          choices = c("Razor", "Unique", "Shared"), hide = if (src != "PD") "Proteome Discoverer only"),
      opt("use_isoform", "Use Peptide Isoform?", "checkbox", FALSE, "Keep protein isoforms apart.",
          hide = if (src != "PD") "Proteome Discoverer only", disabled = "Not supported in the pipeline yet.")))
  )
}
