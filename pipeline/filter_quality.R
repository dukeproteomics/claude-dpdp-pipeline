# Precursor quality filter, pipeline version ("profile"): replaces the app's precursor_quality() when
# --precursor_quality_method=profile (Step 5).
#
# The app (Shiny_Filter.R precursor_quality) takes each precursor's share of (itself + the weaker precursors of its protein) and its CV
# across all samples. Missing values count as 0 in that sum, so when a protein drops to a low level in some groups and its weaker
# precursors go missing, the share jumps there: the app removes good precursors of proteins that change (triple-spike test, PLAN.md).
#
# Profile version: per protein, fit log2 values with median polish (precursor effect + sample effect, missing values skipped). A protein
# that changes between groups is in the sample effects, so only a precursor that does not follow its protein has large residuals.
# Score = SD of a precursor's residuals across its measured samples, as a CV % (sqrt(exp((SD * ln 2)^2) - 1) * 100), so the app's
# precursor_quality_sd cutoff keeps its meaning. Same limits as the app: proteins with 3+ precursors, values below
# precursor_quality_min ignored, precursors with mean intensity >= precursor_quality_intensity and 3+ values are tested. Worst first,
# and a protein always keeps at least 2 precursors.

pq_profile <- function(f, n_s, cv_cut = 50, int_min = 500, val_min = 300) {
  X <- as.matrix(f[(ncol(f) - n_s + 1):ncol(f)]); X[!is.finite(X) | X < val_min] <- NA
  L <- log2(X); avg <- rowMeans(X, na.rm = TRUE); nval <- rowSums(!is.na(X))
  cv <- rep(NA_real_, nrow(f)); bad <- rep(FALSE, nrow(f))
  for (i in split(seq_len(nrow(f)), f$Accession)) {
    i <- i[nval[i] > 0]; if (length(i) < 3) next
    M <- L[i, , drop = FALSE]
    mp <- suppressWarnings(stats::medpolish(M, na.rm = TRUE, trace.iter = FALSE, maxiter = 20))   # non-convergence warnings: residuals still usable
    s <- apply(mp$residuals, 1, stats::sd, na.rm = TRUE)   # stats:: because 05_filter.R has a variable named sd
    c_i <- sqrt(exp((s * log(2))^2) - 1) * 100
    test <- nval[i] >= 3 & avg[i] >= int_min & is.finite(c_i)
    cv[i[test]] <- c_i[test]
    cand <- which(test & c_i > cv_cut); if (!length(cand)) next
    cand <- cand[order(-c_i[cand])][seq_len(min(length(cand), length(i) - 2))]
    bad[i[cand]] <- TRUE
  }
  list(bad = bad, cv = cv, tested = sum(!is.na(cv)))
}

# Apply to precursor_filter (measured values, before noise_keep / background restore) and rebuild the missing-value tables.
pq_profile_apply <- function(db_path, p) {
  n_s <- as.integer(p$sample_number)
  f <- read_table("precursor_filter", db_path); info <- ncol(f) - n_s
  r <- pq_profile(f, n_s, as.numeric(p$precursor_quality_sd), as.numeric(p$precursor_quality_intensity), as.numeric(p$precursor_quality_min))
  kept <- f[!r$bad, , drop = FALSE]
  out <- list(stage = "precursor_quality", method = "profile",
              proteins_removed = length(unique(f$Accession)) - length(unique(kept$Accession)),
              peptides_removed = length(unique(f$Sequence)) - length(unique(kept$Sequence)),
              precursors_removed = sum(r$bad), precursors_tested = r$tested)
  write_table_try("precursor_filter", kept, db_path)
  rebuild_missing_tables(kept, info, db_path)
  out
}
