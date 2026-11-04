# Step 7b - Imputation comparison (no decision; informs the Impute decision).
# Runs several imputation strategies on the same filtered, normalized precursors and scores them, so the analyst can choose with evidence.
#   Start: precursor_filter, scaled by SL factors calculated from measured values (precursor_normdata), log2. Same for every method.
#   Methods: pipeline/impute_methods.R (app duke+bottom-X, mixed kNN/MinProb, Perseus, MinProb, QRILC, kNN), Spectronaut background
#            (--background=<tsv>; remaining gaps by Perseus), no imputation (MaxLFQ with gaps), limpa (DPC-Quant, no imputation).
#   Each: MaxLFQ rollup (iq; limpa uses its own DPC-Quant), limma for the comparisons (limpa: dpcDE).
#   Scores: masking test (hide measured values the way the data loses them, refill, compare), false positives in no-difference splits,
#           SPQC protein CV, ADH1 CV, B vs A calls and how many are measured in one group only.
#
# Usage: Rscript pipeline/steps/07b_impute_compare.R --db=<db> [--background=<spectronaut background tsv>] [--comp_N=B --comp_D=A]
#        [--methods=duke_bottomx,mixed,perseus,minprob,qrilc,knn,background,none,limpa] [--mask_n=10000]
# Output: steps/07_impute/compare.json, review/07b_impute_compare.html

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "impute_methods.R"))
source(file.path(here, "..", "reviews", "07b_impute_compare.R"))
suppressPackageStartupMessages({ library(limma); library(data.table) })

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE); run_dir <- dirname(db_path)
use_app("Shiny_File.R")
p <- get_params(db_path); n_s <- as.integer(p$sample_number)
design <- read_table("design", db_path); design <- design[order(design$Raw_Order), ]
groups_tbl <- read_table("sample_groups", db_path)
grp <- rep(groups_tbl$Group, as.integer(groups_tbl$Count)); labels <- design$Label
pool <- if (!is.null(p$comp_spqc) && p$comp_spqc %in% grp) p$comp_spqc else if ("SPQC" %in% grp) "SPQC" else NA
compN <- args$comp_N %||% "B"; compD <- args$comp_D %||% "A"
methods <- trimws(strsplit(args$methods %||% "duke_bottomx,mixed,perseus,minprob,qrilc,knn,background,none,limpa", ",")[[1]])
mask_n <- as.integer(args$mask_n %||% 10000)
t0 <- Sys.time()

# --- starting matrix: filtered precursors, SL factors from measured values, log2 ---
filt <- read_table("precursor_filter", db_path); info <- ncol(filt) - n_s
X <- as.matrix(filt[, (info + 1):ncol(filt)]); X[!is.na(X) & X <= 0] <- NA; colnames(X) <- labels
nd <- read_table("precursor_normdata", db_path); ND <- as.matrix(nd[, (ncol(nd) - n_s + 1):ncol(nd)])
sl <- mean(colSums(ND, na.rm = TRUE)) / colSums(ND, na.rm = TRUE)
L <- log2(sweep(X, 2, sl, "*"))
prot <- filt$Accession
cat(sprintf("Start: %d precursors, %d proteins, %d gaps (%.1f%%); SL factors from measured values\n", nrow(L), length(unique(prot)), sum(is.na(L)), 100 * mean(is.na(L))))

# --- helpers ---
lfq <- function(M) {
  long <- data.frame(protein_list = rep(prot, ncol(M)), sample_list = rep(colnames(M), each = nrow(M)), quant = as.vector(M), id = rep(seq_len(nrow(M)), ncol(M)))
  long <- long[!is.na(long$quant), ]
  invisible(utils::capture.output(r <- iq::fast_MaxLFQ(long)))
  r$estimate[, colnames(M), drop = FALSE]
}
lim <- function(P, a, b) {
  y <- P[, c(a, b), drop = FALSE]; f <- factor(c(rep("a", length(a)), rep("b", length(b))))
  fit <- eBayes(lmFit(y, model.matrix(~f))); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, p = tt$P.Value, adj = tt$adj.P.Val, row.names = rownames(tt))
}
lim_limpa <- function(Q, a, b) {
  f <- factor(c(rep("a", length(a)), rep("b", length(b))))
  fit <- limpa::dpcDE(Q[, c(a, b)], model.matrix(~f), plot = FALSE)
  if (is.null(fit$t)) fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, p = tt$P.Value, adj = tt$adj.P.Val, row.names = rownames(tt))
}
cvp <- function(P, cols) { v <- 2^P[, cols, drop = FALSE]; ok <- rowSums(is.na(v)) == 0; apply(v[ok, , drop = FALSE], 1, function(r) sd(r) / mean(r) * 100) }
iA <- which(grp == compD); iB <- which(grp == compN)
splits <- if (!is.na(pool)) { iS <- which(grp == pool); lapply(list(c(1, 2), c(1, 3), c(2, 3)), function(s) list(a = c(iA[s[1]], iB[s[1]], iS[s[1]]), b = c(iA[s[2]], iB[s[2]], iS[s[2]]))) } else list()
# measured mask at protein level (any measured precursor in that sample)
meas_prot <- as.matrix(data.table(prot = prot, !is.na(L))[, lapply(.SD, any), by = prot][, -1]); rownames(meas_prot) <- unique(prot)
onegroup <- function(ids) { m <- meas_prot[ids, , drop = FALSE]; (rowSums(m[, iA, drop = FALSE]) == 0) != (rowSums(m[, iB, drop = FALSE]) == 0) }

fake_groups <- function(s) { g <- rep("Z", n_s); g[s$a] <- "X"; g[s$b] <- "Y"; g }
score_proteins <- function(P, stats_fun, make_P = NULL) {
  Pm <- if (is.matrix(P)) P else P$E
  bva <- stats_fun(P, iA, iB)
  fcc <- log2(as.numeric(p$foldchange_cutoff %||% 1.5)); pc <- as.numeric(p$pvalue_cutoff %||% 0.05)
  sig <- !is.na(bva$p) & bva$p <= pc & abs(bva$logFC) >= fcc
  og <- onegroup(rownames(bva))
  nulls <- lapply(splits, function(s) { Ps <- if (is.null(make_P)) P else make_P(fake_groups(s)); stats_fun(Ps, s$a, s$b) })
  null_p01 <- if (length(nulls)) mean(sapply(nulls, function(r) mean(r$p < 0.01, na.rm = TRUE) * 100)) else NA
  null_calls <- if (length(nulls)) mean(sapply(nulls, function(r) sum(!is.na(r$p) & r$p <= pc & abs(r$logFC) >= fcc))) else NA
  spqc <- if (!is.na(pool)) cvp(Pm, which(grp == pool)) else NA
  adh <- if ("P00330" %in% rownames(Pm)) { v <- 2^Pm["P00330", ]; round(sd(v) / mean(v) * 100, 1) } else NA
  list(proteins = nrow(bva), calls = sum(sig), calls_up = sum(sig & bva$logFC > 0), calls_down = sum(sig & bva$logFC < 0),
       calls_adj = sum(!is.na(bva$adj) & bva$adj <= 0.05 & abs(bva$logFC) >= fcc),
       calls_onegroup = sum(sig & og), onegroup_median_abs_log2fc = if (any(sig & og)) round(median(abs(bva$logFC[sig & og])), 2) else NA,
       pi0 = round(propTrueNull(bva$p[!is.na(bva$p)]), 3), null_pct_p01 = round(null_p01, 2), null_calls = round(null_calls, 1),
       spqc_median_cv = round(median(spqc, na.rm = TRUE), 1), spqc_pct_cv_lt20 = round(100 * mean(spqc < 20, na.rm = TRUE), 1), adh1_cv = adh,
       phist = as.integer(tabulate(pmin(20, floor(bva$p[!is.na(bva$p)] * 20) + 1), 20)))
}

# --- background values (Spectronaut background export), same SL factors ---
bg_matrix <- NULL
if ("background" %in% methods) {
  if (is.null(args$background)) { methods <- setdiff(methods, "background"); message("No --background export given: skipping the background method.") }
  else {
    d <- fread(abs_path(args$background), sep = "\t", na.strings = c("", "NA", "NaN", "Filtered"), showProgress = FALSE)
    qc <- grep("EG.TotalQuantity", names(d), value = TRUE)
    B <- as.matrix(d[, ..qc][, lapply(.SD, function(x) suppressWarnings(as.numeric(x)))]); B[!is.na(B) & B <= 0] <- NA
    idx <- match(filt$PrecursorId, d$EG.PrecursorId)
    cat(sprintf("Background export: %d of %d filtered precursors matched by PrecursorId\n", sum(!is.na(idx)), length(idx)))
    bg_matrix <- log2(sweep(B[idx, , drop = FALSE], 2, sl, "*")); colnames(bg_matrix) <- labels
  }
}
fill_background <- function(M) { out <- M; k <- is.na(M) & !is.na(bg_matrix); out[k] <- bg_matrix[k]; imp_perseus(out) }

# --- masking test: hide measured values in complete precursors with the data's own intensity-dependent missingness ---
complete <- which(rowSums(is.na(L)) == 0)
set.seed(1)
sm <- sample(length(L), min(300000, length(L))); rm_ <- rowMeans(L, na.rm = TRUE)
dfm <- data.frame(miss = as.integer(is.na(L))[sm], mu = rm_[((sm - 1) %% nrow(L)) + 1])
gm <- glm(miss ~ mu, data = dfm, family = binomial())
pr <- predict(gm, newdata = data.frame(mu = rm_[complete]), type = "response")
scale <- mask_n / (sum(pr) * n_s)
mask <- matrix(FALSE, nrow(L), ncol(L))
mask[complete, ] <- matrix(runif(length(complete) * n_s), ncol = n_s) < pmin(0.9, pr * scale)
Lm <- L; Lm[mask] <- NA
truth <- L[mask]; mar_cells <- mar_mask(Lm, grp)[mask]
cat(sprintf("Masking test: %d values hidden in %d complete precursors (missingness model: logit(miss) = %.2f %+.2f x mean log2)\n",
            sum(mask), length(complete), coef(gm)[1], coef(gm)[2]))

# --- run the methods ---
res <- list(); masked <- list()
for (m in methods) {
  t1 <- Sys.time(); cat("Method:", m, "... ")
  if (m %in% names(IMPUTE_METHODS)) {
    f <- IMPUTE_METHODS[[m]]$fun
    I <- f(L, grp); P <- lfq(I); res[[m]] <- score_proteins(P, lim, make_P = function(g) lfq(f(L, g)))
    Im <- f(Lm, grp); e <- Im[mask] - truth
    masked[[m]] <- list(rmse = round(sqrt(mean(e^2)), 3), bias = round(mean(e), 3), rmse_mar = round(sqrt(mean(e[mar_cells]^2)), 3),
                        rmse_mnar = round(sqrt(mean(e[!mar_cells]^2)), 3), cor = round(cor(Im[mask], truth), 3))
  } else if (m == "background") {
    I <- fill_background(L); P <- lfq(I); res[[m]] <- score_proteins(P, lim)
  } else if (m == "none") {
    P <- lfq(L); res[[m]] <- score_proteins(P, lim)
  } else if (m == "limpa") {
    dp <- limpa::dpc(L)
    Q <- limpa::dpcQuant(L, protein.id = prot, dpc = dp, verbose = FALSE)
    res[[m]] <- score_proteins(Q, lim_limpa)
    res[[m]]$spqc_median_cv <- { v <- cvp(Q$E, which(grp == pool)); round(median(v), 1) }
    res[[m]]$spqc_pct_cv_lt20 <- { v <- cvp(Q$E, which(grp == pool)); round(100 * mean(v < 20), 1) }
    res[[m]]$adh1_cv <- if ("P00330" %in% rownames(Q$E)) { v <- 2^Q$E["P00330", ]; round(sd(v) / mean(v) * 100, 1) } else NA
  }
  cat(sprintf("%.0f s\n", as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}

labels_m <- c(sapply(IMPUTE_METHODS, `[[`, "label"), background = "Spectronaut background (+ Perseus for empty gaps)",
              none = "no imputation (MaxLFQ with gaps)", limpa = "limpa DPC-Quant (no imputation)")
out <- list(step = "07b_impute_compare", db_path = db_path, comparison = paste0(compN, "_v_", compD), pool = pool,
            start = list(precursors = nrow(L), proteins = length(unique(prot)), gaps = sum(is.na(L))),
            masking = list(hidden = sum(mask), complete_precursors = length(complete), model = as.list(round(coef(gm), 3))),
            methods = lapply(setNames(methods, methods), function(m) c(list(label = unname(labels_m[m])), res[[m]], masked = list(masked[[m]]))),
            seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
sd_ <- step_dir(run_dir, "07_impute")
write_json(out, file.path(sd_, "compare.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
page <- review_impute_compare(db_path)
cat("\n==== Imputation comparison ====\n")
tab <- do.call(rbind, lapply(methods, function(m) { r <- res[[m]]; mk <- masked[[m]]
  data.frame(method = m, mask_rmse = mk$rmse %||% NA, null_p01 = r$null_pct_p01, spqc_cv = r$spqc_median_cv, adh_cv = r$adh1_cv,
             calls = r$calls, onegroup = r$calls_onegroup, adj = r$calls_adj) }))
print(tab, row.names = FALSE)
cat("Review page:", page, "\n")
