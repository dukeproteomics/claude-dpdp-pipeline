# Step 6b - Normalization comparison (no decision; informs the Normalize decision).
# Runs every available normalization method (pipeline/normalize.R) on this run's filtered precursors, fitted on the current
# normalization subset (precursor_normdata, so run Step 6 --action=preview first), and scores each one on measured values only.
#   Scores (no single winner; datasets differ):
#     - median precursor CV within groups, and in the pool (SPQC) group
#     - CV of the constant yeast ADH spike (P00330) across all samples, when present (information: a constant spike is expected to
#       vary after a correct loading normalization when protein loads differ)
#     - intensity trend: in each intensity decile, the per-sample median log ratio to the reference profile, centred across samples;
#       score = mean SD between samples over the deciles (0 = samples agree at every intensity)
#     - RT trend: the same across retention-time deciles (exports with EG.ApexRT)
#     - fold-change spread vs SL: IQR of the between-group log2 FCs relative to SL (1 = same; < 1 = the method shrinks differences,
#       which also lowers CVs, so read it with the CV columns)
#     - median log2 FC of each comparison: what the method assumes is unchanged
#     - known truth, when the group names give species mix ratios (E<n>H<n>Y<n>, three-species spike-in): mean |log2 FC error|
#
# Usage: Rscript pipeline/steps/06b_normalize_compare.R --db=<db> [--methods=sl,tmm,...]
# Output: steps/06_normalize/compare.json, review/06b_normalize_compare.html

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "background.R"))
source(file.path(here, "..", "normalize.R"))
source(file.path(here, "..", "reviews", "06b_normalize_compare.R"))
suppressPackageStartupMessages(library(matrixStats))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE); run_dir <- dirname(db_path)
use_app("Shiny_File.R")
if (!table_exists("precursor_normdata", db_path)) stop("Run Step 6 --action=preview first: precursor_normdata is missing.")
p <- get_params(db_path); n_s <- as.integer(p$sample_number)
design <- read_table("design", db_path); design <- design[order(design$Raw_Order), ]
groups_tbl <- read_table("sample_groups", db_path)
grp <- rep(groups_tbl$Group, as.integer(groups_tbl$Count)); labels <- design$Label
pool <- if (!is.null(p$comp_spqc) && p$comp_spqc %in% grp) p$comp_spqc else if ("SPQC" %in% grp) "SPQC" else NA
t0 <- Sys.time()

filt <- read_table("precursor_filter", db_path); nd <- read_table("precursor_normdata", db_path)
L <- log2_samples(filt, n_s); R <- log2_samples(nd, n_s); colnames(L) <- colnames(R) <- labels
ref_rows <- match(nd$PrecursorId, filt$PrecursorId)
meas <- !is.na(L)
if (is_background_run(db_path)) meas <- meas & samples_of(read_table("precursor_filter_bgmask", db_path), n_s) != 1
rt <- precursor_rt(db_path, filt$PrecursorId)
avail <- names(NORM_METHODS); if (is.null(rt)) avail <- setdiff(avail, c("rt", "rt_loess"))
methods <- if (!is.null(args$methods)) intersect(trimws(strsplit(args$methods, ",")[[1]]), avail) else avail
cat(sprintf("Comparing %d methods on %s precursors (%s in the normalization subset)%s\n", length(methods), fmt_n(nrow(L)), fmt_n(nrow(R)),
            if (is.null(rt)) "; no RT in this export, so rt methods are skipped" else ""))

# --- reference for the trend scores: the subset rows, measured in at least half of the samples ---
in_ref <- seq_len(nrow(L)) %in% ref_rows
# ADH1 is a constant spike only when yeast is a spike, not part of the sample (yeast proteome or species-mix designs)
adh <- grepl("P00330", filt$Accession) & mean(grepl("Saccharomyces", filt$Description)) < 0.02
species <- ifelse(grepl("Homo sapiens", filt$Description), "H", ifelse(grepl("Escherichia", filt$Description), "E",
             ifelse(grepl("Saccharomyces", filt$Description), "Y", NA)))
truth_groups <- grepl("^E[0-9]+H[0-9]+Y[0-9]+$", unique(grp))
use_truth <- sum(truth_groups) >= 2 && mean(!is.na(species)) > 0.9
test_groups <- setdiff(unique(grp), pool); test_groups <- test_groups[sapply(test_groups, function(g) sum(grp == g) >= 2)]
pairs <- if (length(test_groups) >= 2) utils::combn(test_groups, 2, simplify = FALSE) else list()
if (length(pairs) > 10) pairs <- pairs[seq_len(10)]

group_mean <- function(X, g, min_n = 2) { s <- X[, grp == g, drop = FALSE]; m <- rowMeans(s, na.rm = TRUE); m[rowSums(!is.na(s)) < min_n] <- NA; m }
med_cv <- function(X, g) { s <- 2^X[, grp == g, drop = FALSE]; ok <- rowSums(!is.na(s)) == ncol(s); if (sum(ok) < 20) return(NA)
  round(median(rowSds(s[ok, , drop = FALSE]) / rowMeans(s[ok, , drop = FALSE])) * 100, 2) }
trend <- function(X, x = NULL, nb = 10) {   # per-sample median residual by decile of x (subset rows); x NULL = intensity (reference level)
  r <- ref_profile(X); if (is.null(x)) x <- r; k <- in_ref & is.finite(r) & is.finite(x)
  br <- unique(quantile(x[k], seq(0, 1, length.out = nb + 1))); b <- cut(x[k], br, include.lowest = TRUE)
  M <- sapply(seq_len(ncol(X)), function(j) tapply(X[k, j] - r[k], b, median, na.rm = TRUE))
  M <- M - rowMeans(M, na.rm = TRUE)   # a shape shared by all samples (e.g. from missing values at low intensity) is not a bias between samples
  list(mid = as.numeric(tapply(x[k], b, median)), curves = unname(split(t(round(M, 4)), seq_len(ncol(M)))),
       score = round(mean(apply(M, 1, sd, na.rm = TRUE)), 3)) }

score_method <- function(m) {
  Ln <- if (m == "sltmm") { o1 <- off_sl(R); Rs <- shift_cols(R, o1); shift_cols(L, o1 + off_tmm(Rs)) }
        else norm_apply_method(m, R, L, rt = rt, ref_rows = ref_rows)
  X <- Ln; X[!meas] <- NA
  out <- list(label = NORM_METHODS[[m]]$label, family = NORM_METHODS[[m]]$family)
  out$cv_within <- if (length(test_groups)) round(median(sapply(test_groups, function(g) med_cv(X, g)), na.rm = TRUE), 2) else NA
  out$cv_pool <- if (!is.na(pool)) med_cv(X, pool) else NA
  if (sum(adh) >= 3) { A <- 2^X[adh, , drop = FALSE]; A <- A[rowSums(is.na(A)) == 0, , drop = FALSE]
    if (nrow(A) >= 3) { tot <- colSums(A); out$adh_cv <- round(sd(tot) / mean(tot) * 100, 2); out$adh_precursors <- nrow(A) } }
  ti <- trend(X); out$int_trend <- ti$score; out$int_curves <- ti$curves; out$int_mid <- round(ti$mid, 3)
  if (!is.null(rt)) { tr <- trend(X, rt); out$rt_trend <- tr$score; out$rt_curves <- tr$curves; out$rt_mid <- round(tr$mid, 2) }
  fc <- lapply(pairs, function(pr) group_mean(X, pr[2]) - group_mean(X, pr[1]))
  out$fc_iqr <- sapply(fc, function(v) IQR(v, na.rm = TRUE)); out$fc_median <- round(sapply(fc, median, na.rm = TRUE), 3)
  if (use_truth) {
    lev <- function(g, s) as.numeric(sub(paste0(".*", s, "([0-9]+).*"), "\\1", g))
    err <- unlist(lapply(seq_along(pairs), function(i) { pr <- pairs[[i]]
      if (!all(grepl("^E[0-9]+H[0-9]+Y[0-9]+$", pr))) return(NULL)
      sapply(c("H", "E", "Y"), function(s) median(fc[[i]][species %in% s], na.rm = TRUE) - log2(lev(pr[2], s) / lev(pr[1], s))) }))
    out$truth_error <- round(mean(abs(err)), 3); out$truth_bias <- round(err, 3)
  }
  rm(Ln, X); gc(verbose = FALSE)
  out
}

res <- list()
for (m in methods) { t1 <- Sys.time()
  r <- tryCatch(score_method(m), error = function(e) { message(m, " failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  r$seconds <- round(as.numeric(Sys.time() - t1, units = "secs"), 1); res[[m]] <- r
  cat(sprintf("  %-11s CV within %5s  pool %5s  ADH %5s  intensity trend %5.3f  RT trend %5s  %4.1fs\n", m, r$cv_within %||% "-", r$cv_pool %||% "-",
              r$adh_cv %||% "-", r$int_trend, r$rt_trend %||% "-", r$seconds)) }
base_iqr <- res$sl$fc_iqr %||% res[[1]]$fc_iqr
for (m in names(res)) { res[[m]]$fc_spread <- if (length(pairs)) round(mean(res[[m]]$fc_iqr / base_iqr), 3) else NA; res[[m]]$fc_iqr <- NULL }

# ADH CV is circular when the subset is (mostly) ADH itself (--norm_preset=adh): the factors were made to equalize it
adh_circular <- any(adh) && mean(in_ref[adh]) > 0.5 && sum(in_ref) < 0.2 * nrow(L)
out <- list(step = "06b_normalize_compare", db_path = db_path, samples = labels, groups = grp, pool = pool, adh_circular = adh_circular,
            subset = list(precursors = nrow(L), used_for_factors = nrow(R), rule_exclude = p$norm_exclude_grep, rule_include = paste(p$norm_include, p$norm_include_grep)),
            has_rt = !is.null(rt), truth = use_truth, pairs = sapply(pairs, function(pr) paste0(pr[2], " / ", pr[1])),
            seconds = round(as.numeric(Sys.time() - t0, units = "secs")), methods = res)
sd <- step_dir(run_dir, "06_normalize")
write_json(out, file.path(sd, "compare.json"), auto_unbox = TRUE, pretty = FALSE, digits = NA, na = "null")
page <- review_normalize_compare(db_path)
cat("Comparison:", file.path(sd, "compare.json"), "\nReview page:", page, "\n")
