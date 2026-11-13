# Step 4 - Noise: the app's Noise tab (ui.R noise_box, noise_widget_save). Used by steps/04_noise.R (defaults, checks) and
# reviews/04_noise.R (the Options panel). summ: the last preview's summary (counts at the previewed baseline), or NULL.
noise_options <- function(p = NULL, summ = NULL) {
  infl <- suppressWarnings(as.numeric(p$noise_inflection %||% NA))
  have_infl <- !is.na(infl) && infl > 0
  cnt <- function(x) if (is.null(x$count)) "?" else sprintf("%s values (%s%%)", format(as.numeric(x$count), big.mark = ","), x$percent)
  type_info <- if (!is.null(summ)) sprintf("At the last preview: dynamic removes %s below the inflection %s (log2 %s); fixed at %s removes %s.",
    cnt(summ$dynamic), format(summ$noise_inflection %||% "?"), summ$noise_inflection_log2 %||% "?", summ$noise_baseline_value %||% "?", cnt(summ$fixed))
  sections <- list(
    opt_section("Noise reduction", list(
      opt("noise_type", "Select noise reduction strategy", "select", "none",
          paste("none: keep every value. fixed: remove values below the baseline you enter. dynamic: remove values below the inflection",
                "point of the sorted intensity curve (the app's Calculate Noise Inflection Point; recalculated at every preview)."),
          choices = c("none", "fixed", "dynamic"), info = type_info),
      opt("noise_baseline_value", "Enter fixed baseline", "number", 1,
          "Raw intensity (not log2). Values below it are removed. Tip: the reproducibility panel below can put a cutoff here.",
          depends = "noise_type=fixed", min = 0, step = "any"),
      opt("noise_keep", "Keep noise if x value above noise", "checkbox", FALSE,
          paste("Noise values are removed so the filter treats them as missing (precursors that are mostly noise then fail the filter);",
                "for every precursor that passes the filter, its measured values are put back, so they are not imputed (Step 5)."),
          depends = "noise_type!=none"),
      opt("noise_keep_value", "Number of values above noise", "number", 1, "",
          depends = "noise_keep", min = 0, step = 1,
          disabled = "Not used: the pipeline puts back every measured value of a precursor that passes the filter (rule 2026-11-04; the app's version cleans the wrong rows).")))
  )
  attr(sections, "no_preview") <- list(noise_type = "none")   # none removes nothing: apply without a preview (user, 2026-11-13)
  sections
}
