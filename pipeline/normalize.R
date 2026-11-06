# Normalization methods for the DPDP pipeline (Step 6 and the 06b comparison).
#
# Every method has two inputs:
#   R  the normalization subset (precursor_normdata: exclude/include rule applied; background runs: measured values only), log2
#   L  the table to normalize (precursor_filter, all values incl. Spectronaut background), log2, same columns as R
# The method is fitted on R and applied to L, so "what is assumed constant" (the subset) is a separate choice from "what bias is
# corrected" (the method). Missing values stay missing. Returns the normalized log2 L.
#
# The app's sl and sltmm are still run through the app (Shiny_Norm.R) in Step 6 so the tables match a Shiny run; norm_sl() here is the
# same formula, used by the comparison. App methods that were broken or missing (2026-11-05 review):
#   tmm   app divides by the TMM factor only; edgeR factors are relative to each sample's total, so loading was never corrected
#         -> fixed: effective library size (total x factor).
#   mi    called with the wrong arguments -> "median".          lr        used old globals (dpmsr_set) -> "rlr".
#   loess ignored the subset (whole table)  -> fit on subset.  vsn       ignored the subset; ad hoc "/10" rescale -> fit on subset.
#   quantile lost the sample column names and ignored the subset -> target distribution from the subset, columns kept.
#   ai, ti  removed: ti is SL with another target (identical result); ai is SL from column means (differs only through missingness).
#   mti     removed: the app never had the function.          directlfq is a protein rollup method (belongs in Step 8).
# New (current practice): medianratio (median of log ratios to a reference, robust to different missingness), rt (retention-time
# dependent, as DIA-NN / Spectronaut local normalization), none (for designs where total signal is the biology).

NORM_METHODS <- list(
  sltmm = list(label = "SL + TMM", family = "scaling", desc = "SL now, TMM after imputation (lab default)",
               when = "Most proteins unchanged; a minority may change strongly or one-sidedly."),
  sl = list(label = "SL (total)", family = "scaling", desc = "same total intensity of the subset in every sample",
            when = "Equal loading expected and most of the subset unchanged."),
  tmm = list(label = "TMM", family = "scaling", desc = "edgeR trimmed mean of log ratios, with the loading correction",
             when = "Part of the proteome changes, possibly in one direction (composition changes)."),
  median = list(label = "Median", family = "scaling", desc = "same median log2 intensity of the subset in every sample",
                when = "Simple standard (MSstats, Perseus); best when samples have similar numbers of missing values."),
  medianratio = list(label = "Median ratio", family = "scaling", desc = "median log ratio to a reference profile, shared values only",
                     when = "Like median, but robust to samples with different missingness (DESeq-style)."),
  rlr = list(label = "Robust linear regression", family = "intensity", desc = "per-sample robust line against the reference profile",
             when = "Intensity-dependent bias with a small subset (e.g. one spike-in protein) where loess would overfit."),
  loess = list(label = "Loess (intensity)", family = "intensity", desc = "per-sample loess of log ratio vs intensity, fitted on the subset",
               when = "Intensity-dependent bias (saturation, different columns or loads)."),
  rt = list(label = "RT-dependent", family = "rt", desc = "per-sample loess of log ratio vs retention time, fitted on the subset",
            when = "Spray or loading drift along the gradient; long sequences. Needs EG.ApexRT in the export."),
  rt_loess = list(label = "RT + loess", family = "rt", desc = "RT-dependent, then intensity loess",
                  when = "Both drift along the gradient and an intensity-dependent bias. Needs EG.ApexRT in the export."),
  vsn = list(label = "VSN", family = "distribution", desc = "variance stabilizing (glog) fit on the subset",
             when = "Cohorts where variance changes with intensity. Shrinks fold changes at low intensity."),
  quantile = list(label = "Quantile", family = "distribution", desc = "every sample gets the subset's average distribution",
                  when = "Large cohorts of similar samples (plasma, biobanks); not for designs with real global changes."),
  none = list(label = "None", family = "none", desc = "no normalization",
              when = "Totals are the biology: IP / AP-MS, proximity labeling, secretome, or a spike-in already used as the subset."))
norm_method_names <- function() names(NORM_METHODS)
APP_NORM_METHODS <- c("sl", "sltmm")      # run through the app in Step 6

center <- function(v) v - mean(v)
min_present <- function(n) max(2L, ceiling(n / 2))
# reference profile: row mean of rows measured in at least half of the samples (NA elsewhere)
ref_profile <- function(R) { k <- rowSums(!is.na(R)); r <- rowMeans(R, na.rm = TRUE); r[k < min_present(ncol(R))] <- NA; r }
shift_cols <- function(L, s) sweep(L, 2, s, "-")   # log2: subtract a per-sample offset

# --- scaling methods: return per-sample log2 offsets (subtract) ---
off_sl <- function(R) { t <- log2(colSums(2^R, na.rm = TRUE)); center(t) }
off_median <- function(R) center(matrixStats::colMedians(R, na.rm = TRUE))
off_tmm <- function(R) { Y <- 2^R; Y[is.na(Y)] <- 0
  f <- suppressMessages(edgeR::calcNormFactors(Y, method = "TMM", sumTrim = 0.1)); center(log2(colSums(Y) * f)) }
off_medianratio <- function(R, iter = 3) { s <- rep(0, ncol(R))
  for (i in seq_len(iter)) { Rs <- shift_cols(R, s); r <- ref_profile(Rs)
    s <- s + center(apply(Rs - r, 2, median, na.rm = TRUE)) }
  center(s) }
norm_offsets <- function(method, R) switch(method, sl = off_sl(R), sltmm = off_sl(R), median = off_median(R), tmm = off_tmm(R),
                                           medianratio = off_medianratio(R), none = rep(0, ncol(R)), NULL)

# --- curve methods: per-sample lowess of the log ratio to the reference against x (intensity or RT), fitted on the subset ---
lowess_fun <- function(x, y, f) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 50) return(function(z) rep(median(y[ok]), length(z)))
  lo <- stats::lowess(x[ok], y[ok], f = f, iter = 3, delta = 0.01 * diff(range(x[ok])))
  u <- !duplicated(lo$x); xs <- lo$x[u]; ys <- lo$y[u]
  function(z) { o <- stats::approx(xs, ys, xout = z, rule = 2)$y; o[is.na(z)] <- 0; o } }
# x_fit / x_app: NULL = intensity (row reference level, recomputed each pass); else a fixed covariate (RT) per row.
norm_curve <- function(R, L, x_fit = NULL, x_app = NULL, f = 0.4, iter = 3) {
  s <- norm_offsets("medianratio", R); R <- shift_cols(R, s); L <- shift_cols(L, s)
  for (i in seq_len(iter)) {
    r <- ref_profile(R); ra <- rowMeans(L, na.rm = TRUE)
    xf <- if (is.null(x_fit)) r else x_fit; xa <- if (is.null(x_app)) ra else x_app
    for (j in seq_len(ncol(R))) { g <- lowess_fun(xf, R[, j] - r, f)
      R[, j] <- R[, j] - g(xf); L[, j] <- L[, j] - g(xa) }
  }
  L }

norm_rlr <- function(R, L) {
  s <- norm_offsets("medianratio", R); R <- shift_cols(R, s); L <- shift_cols(L, s); r <- ref_profile(R)
  for (j in seq_len(ncol(R))) { ok <- is.finite(r) & is.finite(R[, j])
    cf <- stats::coef(MASS::rlm(R[ok, j] ~ r[ok], maxit = 50)); L[, j] <- (L[, j] - cf[1]) / cf[2] }
  L }

norm_vsn <- function(R, L) {
  fit <- suppressMessages(vsn::vsn2(2^R, verbose = FALSE))
  out <- vsn::predict(fit, newdata = 2^L, log2scale = TRUE); out[is.na(L)] <- NA
  out - median(out, na.rm = TRUE) + median(L, na.rm = TRUE)   # keep the data's overall level
}

norm_quantile <- function(R, L) {
  Q <- preprocessCore::normalize.quantiles(R)                 # target = the subset's mean distribution
  for (j in seq_len(ncol(L))) { ok <- !is.na(R[, j]); x <- sort(R[ok, j]); y <- sort(Q[ok, j])
    u <- !duplicated(x); x <- x[u]; y <- y[u]; z <- L[, j]
    o <- stats::approx(x, y, xout = z, rule = 2)$y
    lo <- !is.na(z) & z < x[1]; hi <- !is.na(z) & z > x[length(x)]      # beyond the subset's range: keep the end offset
    o[lo] <- z[lo] + (y[1] - x[1]); o[hi] <- z[hi] + (y[length(y)] - x[length(x)]); L[, j] <- o }
  L }

# R, L: log2 matrices (same columns). ref_rows: which rows of L are the rows of R (for rt). rt: per-row RT of L (NA allowed).
norm_apply_method <- function(method, R, L, rt = NULL, ref_rows = NULL) {
  if (!method %in% names(NORM_METHODS)) stop("unknown normalization method: ", method)
  off <- norm_offsets(method, R)
  if (!is.null(off)) return(shift_cols(L, off))
  switch(method,
    loess = norm_curve(R, L, f = 0.4),
    rlr = norm_rlr(R, L),
    vsn = norm_vsn(R, L),
    quantile = norm_quantile(R, L),
    rt = , rt_loess = { if (is.null(rt) || all(is.na(rt))) stop("rt normalization needs retention times (EG.ApexRT in the Spectronaut export).")
           if (is.null(ref_rows)) stop("rt normalization needs ref_rows")
           Lr <- norm_curve(R, L, x_fit = rt[ref_rows], x_app = rt, f = 0.1)
           noRT <- is.na(rt); if (any(noRT)) Lr[noRT, ] <- shift_cols(L[noRT, , drop = FALSE], norm_offsets("medianratio", R))
           if (method == "rt_loess") { Rr <- norm_curve(R, R, x_fit = rt[ref_rows], x_app = rt[ref_rows], f = 0.1); Lr <- norm_curve(Rr, Lr, f = 0.4) }
           Lr })
}

# Per-precursor retention time: median EG.ApexRT over runs, from precursor_raw. NULL if the export has none.
precursor_rt <- function(db_path, precursor_ids) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path); on.exit(RSQLite::dbDisconnect(con))
  if (!RSQLite::dbExistsTable(con, "precursor_raw")) return(NULL)
  f <- RSQLite::dbListFields(con, "precursor_raw"); rc <- grep("EG\\.ApexRT$", f, value = TRUE)
  if (!length(rc) || !"EG.PrecursorId" %in% f) return(NULL)
  q <- sprintf("SELECT %s FROM precursor_raw", paste0('"', c("EG.PrecursorId", rc), '"', collapse = ", "))
  d <- RSQLite::dbGetQuery(con, q)
  m <- sapply(d[-1], function(x) suppressWarnings(as.numeric(as.character(x))))
  rt <- matrixStats::rowMedians(m, na.rm = TRUE)
  rt <- tapply(rt, d$EG.PrecursorId, median, na.rm = TRUE)
  out <- as.numeric(rt[as.character(precursor_ids)]); out[!is.finite(out)] <- NA
  if (all(is.na(out))) NULL else out
}

# Which methods can run on this run (rt needs retention times).
norm_available <- function(db_path, precursor_ids = NULL) {
  m <- names(NORM_METHODS)
  has_rt <- !is.null(precursor_ids) && !is.null(precursor_rt(db_path, precursor_ids[seq_len(min(length(precursor_ids), 2000))]))
  if (!has_rt) m <- setdiff(m, c("rt", "rt_loess"))
  m
}

log2_samples <- function(df, n_s) { x <- as.matrix(df[(ncol(df) - n_s + 1):ncol(df)]); storage.mode(x) <- "numeric"; x[!is.na(x) & x <= 0] <- NA; log2(x) }
