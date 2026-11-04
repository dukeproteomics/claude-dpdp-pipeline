# Imputation methods for the DPDP pipeline (raw Spectronaut input), all on a log2 precursor matrix (rows = precursors, cols = samples).
# Each function returns a complete log2 matrix of the same shape. `groups` is the group label of each column.
# MAR vs MNAR rule used by the mixed methods (same idea as the app's duke): a gap is treated as missing at random (MAR) when its precursor
# is measured in at least `mar_min` (default 50%) of its group's samples; otherwise as missing because it is low (MNAR).

mar_mask <- function(x, groups, mar_min = 0.5) {
  m <- matrix(FALSE, nrow(x), ncol(x))
  for (g in unique(groups)) { cc <- which(groups == g); frac <- rowMeans(!is.na(x[, cc, drop = FALSE]))
    m[, cc] <- is.na(x[, cc, drop = FALSE]) & (frac >= mar_min) }
  m
}

# Perseus-style down-shifted normal, per sample: N(mean - shift*sd, (width*sd)^2) of the sample's measured values.
imp_perseus <- function(x, groups = NULL, shift = 1.8, width = 0.3, seed = 123) {
  set.seed(seed)
  for (j in seq_len(ncol(x))) { v <- x[, j]; o <- v[!is.na(v)]; k <- is.na(v)
    x[k, j] <- rnorm(sum(k), mean(o) - shift * sd(o), width * sd(o)) }
  x
}

imp_minprob <- function(x, groups = NULL, q = 0.01, seed = 123) { set.seed(seed); invisible(utils::capture.output(r <- imputeLCMD::impute.MinProb(x, q = q, tune.sigma = 1))); r }

imp_qrilc <- function(x, groups = NULL, seed = 123) { set.seed(seed); imputeLCMD::impute.QRILC(x, tune.sigma = 1)[[1]] }

# kNN on row-centred profiles: neighbours are precursors with the same shape across samples; the precursor's own level is added back.
# (Plain impute.knn on raw log2 picks neighbours by absolute level and did much worse in the 11488 masking test: RMSE 1.94 vs 0.50.)
imp_knn <- function(x, groups = NULL, k = 10) {
  mu <- rowMeans(x, na.rm = TRUE); mu[!is.finite(mu)] <- 0
  invisible(utils::capture.output(r <- suppressWarnings(impute::impute.knn(x - mu, k = k, rowmax = 0.9, colmax = 0.9, rng.seed = 123))))
  r$data + mu
}

# Mixed: kNN for MAR gaps, MinProb for MNAR gaps.
imp_mixed <- function(x, groups, mar_min = 0.5, k = 10, seed = 123) {
  mar <- mar_mask(x, groups, mar_min)
  a <- imp_knn(x, k = k); b <- imp_minprob(x, seed = seed)
  out <- x; out[mar] <- a[mar]; mnar <- is.na(x) & !mar; out[mnar] <- b[mnar]
  out
}

# The app's duke + bottom-X (impute_duke / impute_bottomx / impute_meta_data_bg), re-implemented on log2:
#   duke: per group, rows missing in <= Count*(100-missing_cutoff)/100 samples get group mean + U(-1,1) * SD of the group's
#         intensity bin (20 bins of the row average; SD = mean of the bin's row SDs).
#   bottom-X: every remaining gap gets U(min, X-th percentile) of all measured values.
imp_duke_bottomx <- function(x, groups, missing_cutoff = 50, bottom_x = 2, seed = 123) {
  set.seed(seed)
  out <- x
  for (g in unique(groups)) { cc <- which(groups == g); sub <- x[, cc, drop = FALSE]; n <- length(cc)
    if (n < 2) next
    avg <- rowMeans(sub, na.rm = TRUE); sdv <- apply(sub, 1, sd, na.rm = TRUE)
    ok <- !is.na(sdv); bin <- rep(NA_integer_, length(avg)); bin[ok] <- dplyr::ntile(avg[ok], 20)
    bsd <- tapply(sdv[ok], bin[ok], mean); brk <- tapply(avg[ok], bin[ok], min)
    miss <- rowSums(is.na(sub)); rows <- which(miss > 0 & miss <= n * (100 - missing_cutoff) / 100)
    if (length(rows)) { b <- findInterval(avg[rows], brk); b[b < 1] <- 1; s <- bsd[b]
      for (r in seq_along(rows)) { i <- rows[r]; k <- which(is.na(sub[i, ])); out[i, cc[k]] <- avg[i] + runif(length(k), -1, 1) * s[r] } }
  }
  o <- x[!is.na(x) & x > 0]; lo <- min(o); hi <- sort(o)[ceiling(bottom_x * length(o) / 100)]
  k <- is.na(out); out[k] <- lo + runif(sum(k)) * (hi - lo)
  out
}

IMPUTE_METHODS <- list(
  duke_bottomx = list(fun = imp_duke_bottomx, label = "app: duke + bottom-X 2%"),
  mixed        = list(fun = imp_mixed,        label = "mixed: kNN (MAR) + MinProb (MNAR)"),
  perseus      = list(fun = imp_perseus,      label = "Perseus down-shift (1.8 SD, width 0.3)"),
  minprob      = list(fun = imp_minprob,      label = "MinProb"),
  qrilc        = list(fun = imp_qrilc,        label = "QRILC"),
  knn          = list(fun = imp_knn,          label = "kNN (MAR only)"))
