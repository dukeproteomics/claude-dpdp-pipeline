# Known-ratio imputation benchmark: three-species spike-in (human constant, E. coli and yeast at known levels; design after
# Guzman et al. 2024, Nat Biotechnol, the "Olsen" mixed proteome). Scores the pipeline's imputation options against the true ratios.
#
# Input: a raw Spectronaut precursor export (empty/"Filtered" = missing) and, optionally, the background export of the same analysis
# (filled = empty in raw, present in background), plus the sample list (ID, Label, Group). Group names must look like E5H50Y45.
# Columns whose ID is not in the sample list's real groups (for example a synthetic SPQC) are dropped.
#
# Same processing for every method:
#   filter: precursor measured in >= 2 of 3 replicates in at least one group; single-species protein groups only.
#   normalize: one factor per sample so the human precursors measured in every sample have equal totals (human is constant).
#   rollup: MaxLFQ (iq); limpa uses DPC-Quant. Test: limma (limpa: dpcDE). Calls: adj p <= 0.05 and |FC| >= 1.5.
# Scores per comparison (all 6 group pairs):
#   true calls: E. coli / yeast proteins with an expected |log2 FC| >= log2(1.5), called in the right direction;
#   false calls: human proteins called, or a spike protein called in the wrong direction; observed FDR = false / all calls;
#   log2 FC accuracy per species: median error vs expected, and spread (IQR).
# Plus the masking test (as in pipeline step 07b).
#
# Usage: Rscript tools/impute_benchmark/triple_spike.R --raw=<tsv> [--background=<tsv>] --samples=<xlsx> [--out=data/<name>] [--methods=...]

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
pl <- file.path(here, "..", "..", "pipeline")
source(file.path(pl, "common.R")); source(file.path(pl, "review.R")); source(file.path(pl, "impute_methods.R"))
suppressPackageStartupMessages({ library(limma); library(data.table) })

args <- parse_args()
need(args, "raw"); need(args, "samples")
out_dir <- args$out %||% "data/triple_spike_bench"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); out_dir <- abs_path(out_dir)
methods <- trimws(strsplit(args$methods %||% "duke_bottomx,mixed,mixed_blind,perseus,minprob,qrilc,knn,background,none,limpa", ",")[[1]])
if (is.null(args$background)) methods <- setdiff(methods, "background")
FC <- log2(1.5); t0 <- Sys.time()

# --- read ---
sl_tab <- as.data.frame(readxl::read_excel(abs_path(args$samples)))
sl_tab <- sl_tab[grepl("^E\\d+H\\d+Y\\d+$", sl_tab$Group), ]
read_q <- function(f) {
  d <- fread(abs_path(f), sep = "\t", na.strings = c("", "NA", "NaN", "Filtered"), showProgress = FALSE)
  qc <- grep("EG.TotalQuantity", names(d), value = TRUE)
  id <- sub("^\\[\\d+\\] ID(\\d+)_\\d+\\..*$", "\\1", qc)
  keep <- match(as.character(sl_tab$ID), id)
  if (anyNA(keep)) stop("Sample list IDs not found in the export: ", paste(sl_tab$ID[is.na(keep)], collapse = ", "))
  M <- as.matrix(d[, qc[keep], with = FALSE][, lapply(.SD, function(x) suppressWarnings(as.numeric(x)))])
  M[!is.na(M) & M <= 0] <- NA; colnames(M) <- sl_tab$Label
  list(info = d[, .(PrecursorId = EG.PrecursorId, Accession = PG.ProteinAccessions, Organisms = PG.Organisms)], M = M,
       dropped = setdiff(qc, qc[keep]))
}
R <- read_q(args$raw)
cat(sprintf("Raw export: %d precursors, %d samples used; dropped columns: %s\n", nrow(R$M), ncol(R$M), paste(R$dropped, collapse = ", ")))
grp <- sl_tab$Group; labels <- sl_tab$Label; n_s <- length(grp); glev <- unique(grp)

species <- function(o) { parts <- strsplit(o, ";", fixed = TRUE)
  vapply(parts, function(p) { s <- unique(ifelse(grepl("Homo sapiens", p), "human", ifelse(grepl("Saccharomyces", p), "yeast",
                                ifelse(grepl("Escherichia coli", p), "ecoli", "other")))); if (length(s) == 1) s else "mixed" }, "") }
sp_prec <- species(R$info$Organisms)

# --- filter + normalize ---
meas23 <- sapply(glev, function(g) rowSums(!is.na(R$M[, grp == g, drop = FALSE])) >= 2)
keep <- rowSums(meas23) > 0 & sp_prec %in% c("human", "yeast", "ecoli")
X <- R$M[keep, ]; info <- R$info[keep]; sp_p <- sp_prec[keep]; prot <- info$Accession
hum_full <- sp_p == "human" & rowSums(is.na(X)) == 0
nf <- mean(colSums(X[hum_full, ])) / colSums(X[hum_full, ])
L <- log2(sweep(X, 2, nf, "*"))
sp_prot <- tapply(sp_p, prot, `[`, 1)
cat(sprintf("Start: %d precursors, %d proteins (%s), %d gaps (%.1f%%); factors from %d complete human precursors: %s\n",
            nrow(L), length(unique(prot)), paste(names(table(sp_prot)), table(sp_prot), sep = " ", collapse = ", "),
            sum(is.na(L)), 100 * mean(is.na(L)), sum(hum_full), paste(round(nf, 3), collapse = " ")))
gaps_per_sample <- setNames(colSums(is.na(L)), labels)

bg <- NULL
if ("background" %in% methods) {
  Bx <- read_q(args$background)
  idx <- match(info$PrecursorId, Bx$info$PrecursorId)
  bg <- log2(sweep(Bx$M[idx, , drop = FALSE], 2, nf, "*"))
  meas_same <- !is.na(L) & !is.na(bg)
  cat(sprintf("Background export: %d of %d gaps filled; measured values identical: %s (max |diff| %.2g log2)\n",
              sum(is.na(L) & !is.na(bg)), sum(is.na(L)), all(abs(L[meas_same] - bg[meas_same]) < 1e-6), max(abs(L[meas_same] - bg[meas_same]))))
}

# --- expected ratios ---
lev <- t(sapply(glev, function(g) as.numeric(regmatches(g, regexec("^E(\\d+)H(\\d+)Y(\\d+)$", g))[[1]][2:4])))
colnames(lev) <- c("ecoli", "human", "yeast")
pairs <- combn(seq_along(glev), 2, simplify = FALSE)
comps <- lapply(pairs, function(pp) { d <- glev[pp[1]]; n <- glev[pp[2]]
  list(name = paste0(n, "_v_", d), N = n, D = d, expected = as.list(round(log2(lev[n, ] / lev[d, ]), 3))) })

# --- helpers (as in step 07b) ---
lfq <- function(M) {
  long <- data.frame(protein_list = rep(prot, ncol(M)), sample_list = rep(colnames(M), each = nrow(M)), quant = as.vector(M), id = rep(seq_len(nrow(M)), ncol(M)))
  long <- long[!is.na(long$quant), ]
  invisible(utils::capture.output(r <- iq::fast_MaxLFQ(long)))
  r$estimate[, colnames(M), drop = FALSE]
}
lim <- function(P, a, b) {
  y <- P[, c(a, b), drop = FALSE]; f <- factor(c(rep("a", length(a)), rep("b", length(b))))
  ok <- rowSums(!is.na(y[, seq_along(a), drop = FALSE])) >= 1 & rowSums(!is.na(y[, length(a) + seq_along(b), drop = FALSE])) >= 1
  fit <- suppressWarnings(eBayes(lmFit(y[ok, , drop = FALSE], model.matrix(~f)))); tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, p = tt$P.Value, adj = tt$adj.P.Val, row.names = rownames(tt))
}
lim_limpa <- function(Q, a, b) {
  f <- factor(c(rep("a", length(a)), rep("b", length(b))))
  fit <- limpa::dpcDE(Q[, c(a, b)], model.matrix(~f), plot = FALSE)
  if (is.null(fit$t)) fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  data.frame(logFC = tt$logFC, p = tt$P.Value, adj = tt$adj.P.Val, row.names = rownames(tt))
}
q5 <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return(NULL); round(unname(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95))), 3) }

score <- function(P, stats_fun) {
  Pm <- if (is.matrix(P)) P else P$E
  per <- lapply(comps, function(cp) {
    r <- stats_fun(P, which(grp == cp$D), which(grp == cp$N)); s <- sp_prot[rownames(r)]
    ex <- unlist(cp$expected)[s]; called <- !is.na(r$adj) & r$adj <= 0.05 & abs(r$logFC) >= FC
    should <- s != "human" & abs(ex) >= FC
    tp <- sum(called & should & sign(r$logFC) == sign(ex))
    fp_h <- sum(called & s == "human"); fp_dir <- sum(called & s != "human" & sign(r$logFC) != sign(ex))
    acc <- lapply(c(human = "human", ecoli = "ecoli", yeast = "yeast"), function(k) { e <- r$logFC[s == k] - cp$expected[[k]]
      list(n = sum(s == k & !is.na(r$logFC)), median_error = round(median(e, na.rm = TRUE), 3), iqr = round(IQR(e, na.rm = TRUE), 3), q = q5(r$logFC[s == k])) })
    list(name = cp$name, tp = tp, possible = sum(should), sensitivity = round(100 * tp / sum(should), 1), fp_human = fp_h, fp_direction = fp_dir,
         fdr = round(100 * (fp_h + fp_dir) / max(1, sum(called)), 2), human_p01 = round(100 * mean(r$p[s == "human"] < 0.01, na.rm = TRUE), 2), species = acc)
  })
  names(per) <- vapply(comps, `[[`, "", "name")
  big <- vapply(comps, function(cp) max(abs(unlist(cp$expected[c("ecoli", "yeast")]))) >= FC, TRUE)
  list(proteins = nrow(Pm), per = per,
       tp_total = sum(sapply(per, `[[`, "tp")), fp_total = sum(sapply(per, function(x) x$fp_human + x$fp_direction)),
       fdr_mean = round(mean(sapply(per, `[[`, "fdr")), 2), fdr_max = max(sapply(per, `[[`, "fdr")),
       sensitivity_mean = round(mean(sapply(per[big], `[[`, "sensitivity")), 1),
       human_p01_mean = round(mean(sapply(per, `[[`, "human_p01")), 2),
       abs_error_mean = round(mean(unlist(lapply(per, function(x) abs(c(x$species$ecoli$median_error, x$species$yeast$median_error))))), 3),
       iqr_mean = round(mean(unlist(lapply(per, function(x) c(x$species$ecoli$iqr, x$species$yeast$iqr)))), 3),
       human_iqr_mean = round(mean(sapply(per, function(x) x$species$human$iqr)), 3))
}

# --- masking test ---
complete <- which(rowSums(is.na(L)) == 0)
set.seed(1); sm <- sample(length(L), min(300000, length(L))); rm_ <- rowMeans(L, na.rm = TRUE)
gm <- glm(miss ~ mu, data = data.frame(miss = as.integer(is.na(L))[sm], mu = rm_[((sm - 1) %% nrow(L)) + 1]), family = binomial())
pr <- predict(gm, newdata = data.frame(mu = rm_[complete]), type = "response")
mask_n <- as.integer(args$mask_n %||% 20000); sc <- mask_n / (sum(pr) * n_s)
mask <- matrix(FALSE, nrow(L), ncol(L)); mask[complete, ] <- matrix(runif(length(complete) * n_s), ncol = n_s) < pmin(0.9, pr * sc)
Lm <- L; Lm[mask] <- NA; truth <- L[mask]
cat(sprintf("Masking test: %d values hidden in %d complete precursors\n", sum(mask), length(complete)))

# --- run ---
res <- list(); masked <- list()
for (m in methods) {
  t1 <- Sys.time(); cat("Method:", m, "... ")
  if (m %in% names(IMPUTE_METHODS)) {
    f <- IMPUTE_METHODS[[m]]$fun
    res[[m]] <- score(lfq(f(L, grp)), lim)
    e <- f(Lm, grp)[mask] - truth
    masked[[m]] <- list(rmse = round(sqrt(mean(e^2)), 3), bias = round(mean(e), 3))
  } else if (m == "background") {
    I <- L; k <- is.na(L) & !is.na(bg); I[k] <- bg[k]; res[[m]] <- score(lfq(imp_perseus(I)), lim)
  } else if (m == "none") {
    res[[m]] <- score(lfq(L), lim)
  } else if (m == "limpa") {
    Q <- limpa::dpcQuant(L, protein.id = prot, dpc = limpa::dpc(L), verbose = FALSE); res[[m]] <- score(Q, lim_limpa)
  }
  cat(sprintf("%.0f s\n", as.numeric(difftime(Sys.time(), t1, units = "secs"))))
}

labels_m <- c(sapply(IMPUTE_METHODS, `[[`, "label"), background = "Spectronaut background (+ Perseus for empty gaps)",
              none = "no imputation (MaxLFQ with gaps)", limpa = "limpa DPC-Quant (no imputation)")
out <- list(raw = abs_path(args$raw), background = if (!is.null(args$background)) abs_path(args$background) else NULL,
            groups = glev, comparisons = comps, gaps_per_sample = as.list(gaps_per_sample), norm_factors = as.list(setNames(round(nf, 4), labels)),
            start = list(precursors = nrow(L), proteins = length(unique(prot)), species = as.list(table(sp_prot)), gaps = sum(is.na(L))),
            masking = list(hidden = sum(mask), complete_precursors = length(complete)),
            methods = lapply(setNames(methods, methods), function(m) c(list(label = unname(labels_m[m])), res[[m]], masked = list(masked[[m]]))),
            seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
write_json(out, file.path(out_dir, "benchmark.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")

cat("\n==== Triple-spike benchmark (adj p <= 0.05, |FC| >= 1.5; 6 comparisons) ====\n")
tab <- do.call(rbind, lapply(methods, function(m) { r <- res[[m]]
  data.frame(method = m, mask_rmse = masked[[m]]$rmse %||% NA, proteins = r$proteins, true_calls = r$tp_total, false_calls = r$fp_total,
             fdr_mean = r$fdr_mean, fdr_max = r$fdr_max, sens_mean = r$sensitivity_mean, human_p01 = r$human_p01_mean,
             fc_abs_err = r$abs_error_mean, fc_iqr = r$iqr_mean, human_iqr = r$human_iqr_mean) }))
print(tab, row.names = FALSE)
source(file.path(here, "review_triple_spike.R"))
cat("Review page:", review_triple_spike(out_dir), "\n")
