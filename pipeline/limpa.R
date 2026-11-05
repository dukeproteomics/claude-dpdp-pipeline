# limpa (Bioconductor; Li & Smyth): no imputation. A detection probability curve (DPC) says how likely a precursor is to be detected
# at a given intensity; DPC-Quant uses it to estimate every protein in every sample with a standard error; dpcDE is limma with
# precision weights from those standard errors. In the triple-spike benchmark the weights carry most of the benefit
# (false call rate 3.1 / 5.1% mean / worst with dpcDE, 4.2 / 8.3% with plain limma on the same values).
#
# Pipeline wiring:
#   Step 7 --impute_type=limpa : precursor_impute_<norm> = the normalized table with its gaps (sltmm: TMM from measured values).
#   Step 8 --rollup_method=limpa : DPC-Quant per normalization -> protein_<norm> (2^estimate, app table layout), peptide_<norm> (sum of
#          measured precursors), steps/08_rollup/limpa_<norm>.rds (EList with standard errors, DPC).
#   Step 10 : the app's stat_calc_bg builds the tables (CVs, missing factor filter, layout); then FC, FC2, pval, adjpval and Stats are
#          replaced with dpcDE results for the comparison's samples. peptide_refilter must be off (it re-rolls precursors).
# The DPC is estimated from measured values only (background values blanked). In DPC-Quant background values are left missing by default
# (limpa_background=missing; see limpa_quant).

limpa_rds <- function(run_dir, norm) file.path(run_dir, "steps", "08_rollup", paste0("limpa_", norm, ".rds"))

limpa_dpc <- function(L, bgmask = NULL) { Lm <- L; if (!is.null(bgmask)) Lm[bgmask] <- NA; limpa::dpc(Lm) }

# DPC-Quant for one normalization; saves and returns list(Q = EList, dpc, n_measured = measured precursors per protein and sample)
# background = "missing" (default): background values are left to the detection curve; "use": they count as observations.
# Triple-spike background export, E45/E5 and E10/E5 false call rate: use 6.9 / 9.2%, missing 4.9 / 5.6%, raw export 4.0 / 5.6%;
# E. coli log2 FC (true 3.17): use +2.37, missing +2.82, raw +3.03. Background signal is too noisy to count as a measurement here.
limpa_quant <- function(db_path, norm, background = "missing") {
  n_s <- as.integer(get_param("sample_number", db_path))
  d <- read_table(paste0("precursor_impute_", norm), db_path); info <- ncol(d) - n_s
  X <- as.matrix(d[(info + 1):ncol(d)]); X[!is.na(X) & X <= 0] <- NA; L <- log2(X)
  bgm <- if (is_background_run(db_path)) as.matrix(read_table("precursor_filter_bgmask", db_path)[(info + 1):ncol(d)]) == 1 else NULL
  if (!anyNA(L)) stop("limpa_quant: precursor_impute_", norm, " has no gaps, so it was imputed. Apply Step 7 with --impute_type=limpa first.")
  dp <- limpa_dpc(L, bgm)
  if (!is.null(bgm) && background == "missing") L[bgm] <- NA   # background values left to the detection curve
  p <- get_params(db_path); id <- feature_key(d, p)   # protein accession, or the peptide's grouping columns for Peptide output
  Q <- limpa::dpcQuant(L, protein.id = id, dpc = dp, verbose = FALSE)
  meas <- !is.na(L); if (!is.null(bgm)) meas <- meas & !bgm
  nm <- rowsum(1L * meas, id, reorder = TRUE)[rownames(Q$E), , drop = FALSE]
  obj <- list(Q = Q, dpc = dp$dpc, background = if (is.null(bgm)) "none" else background, n_measured = nm, norm = norm, made = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  dir.create(dirname(limpa_rds(dirname(db_path), norm)), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, limpa_rds(dirname(db_path), norm))
  obj
}

# protein_<norm> and peptide_<norm> in the app's layout from the DPC-Quant result
limpa_write_tables <- function(db_path, norm, obj) {
  n_s <- as.integer(get_param("sample_number", db_path)); p <- get_params(db_path)
  if (is_peptide_output(p)) {   # the app's peptide layout (made by rollup_apply_bg from the gapped table), sample values from DPC-Quant
    t <- paste0("peptide_impute_", norm); pd <- read_table(t, db_path); k <- match(feature_key(pd, p), rownames(obj$Q$E))
    if (anyNA(k)) stop("limpa_write_tables: ", sum(is.na(k)), " peptides of ", t, " are not in the DPC-Quant result.")
    sc <- (ncol(pd) - n_s + 1):ncol(pd)
    pd[sc] <- as.data.frame(2^obj$Q$E[k, , drop = FALSE])
    pd$Detected_Imputed <- apply(obj$n_measured[k, , drop = FALSE], 1, function(r) paste(ifelse(r == 0, "-", as.character(r)), collapse = "."))
    write_table_try(t, pd, db_path)
    return(list(peptides = nrow(pd)))
  }
  d <- read_table(paste0("precursor_impute_", norm), db_path); info <- ncol(d) - n_s; sc <- names(d)[(info + 1):ncol(d)]
  acc <- rownames(obj$Q$E); first <- match(acc, d$Accession)
  # Detected_Imputed: measured precursors per sample, "n.n.n" with 0 as "-" (the app's reduce_imputed_df format); the app's Stats reads it
  # for the missing factor when peptide_refilter is off
  nm <- obj$n_measured[acc, , drop = FALSE]
  di <- apply(nm, 1, function(r) paste(ifelse(r == 0, "-", as.character(r)), collapse = "."))
  prot <- data.frame(d[first, c("Accession", "Description", "Name", "Genes")], Precursors = as.integer(table(d$Accession)[acc]),
                     Detected_Imputed = di, 2^obj$Q$E, check.names = FALSE, row.names = NULL)
  names(prot)[7:ncol(prot)] <- sc
  write_table_try(paste0("protein_", norm), prot, db_path)
  key <- c("Accession", "Description", "Name", "Genes", "Sequence", "PeptidePosition")
  X <- as.matrix(d[sc]); g <- do.call(paste, c(d[key], sep = "\r"))
  s <- rowsum(X, g, reorder = FALSE, na.rm = TRUE); cnt <- rowsum(1L * !is.na(X), g, reorder = FALSE); s[cnt == 0] <- NA
  firstp <- match(rownames(s), g)
  pep <- data.frame(d[firstp, key], Precursors = as.integer(table(g)[rownames(s)]), s, check.names = FALSE, row.names = NULL)
  names(pep)[(length(key) + 2):ncol(pep)] <- sc
  write_table_try(paste0("peptide_", norm), pep, db_path)
  list(proteins = nrow(prot), peptides = nrow(pep))
}

# Replace the app's FC / FC2 / pval / adjpval / Stats in comparison i's final table with dpcDE results.
limpa_stats <- function(db_path, stats_comp, i, p) {
  obj <- readRDS(limpa_rds(dirname(db_path), p$stat_norm)); Q <- obj$Q
  loc <- function(x) as.integer(trimws(strsplit(as.character(x), ",")[[1]]))
  nl <- loc(stats_comp$N_loc[i]); dl <- loc(stats_comp$D_loc[i]); nm <- stats_comp$Name[i]
  f <- factor(c(rep("D", length(dl)), rep("N", length(nl))), levels = c("D", "N"))
  fit <- limpa::dpcDE(Q[, c(dl, nl)], stats::model.matrix(~f), plot = FALSE)
  if (is.null(fit$t)) fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  adj_all <- stats::p.adjust(tt$P.Value, method = p$padjust_options)
  tab <- final_table(stats_comp, i, p); df <- read_table(tab, db_path); m <- match(feature_key(df, p), rownames(tt))
  if (anyNA(m)) stop("limpa_stats: ", sum(is.na(m)), " proteins in ", tab, " are not in the DPC-Quant result.")
  ratio <- 2^tt$logFC[m]
  df[[paste0(nm, "_FC")]] <- signif(ifelse(ratio >= 1, ratio, -1 / ratio), 7)
  df[[paste0(nm, "_FC2")]] <- signif(ratio, 7)
  df[[paste0(nm, "_pval")]] <- tt$P.Value[m]
  if (as_flag(p$checkbox_adjpval)) df[[paste0(nm, "_adjpval")]] <- adj_all[m]
  fcc <- as.numeric(p$foldchange_cutoff); pc <- as.numeric(p$pvalue_cutoff)
  pv <- if (as_flag(p$checkbox_filter_adjpval)) df[[paste0(nm, "_adjpval")]] else df[[paste0(nm, "_pval")]]
  df$Stats <- ""; df$Stats[df[[paste0(nm, "_FC")]] >= fcc & pv <= pc] <- "Up"; df$Stats[df[[paste0(nm, "_FC")]] <= -fcc & pv <= pc] <- "Down"
  df <- df[order(df$Stats, -df[[paste0(nm, "_pval")]], decreasing = TRUE), ]
  write_table_try(tab, df, db_path)
  list(test = "limpa dpcDE (limma with DPC-Quant precision weights)", samples = c(N = length(nl), D = length(dl)),
       proteins_tested = nrow(tt), df_prior = round(fit$df.prior[1], 2))
}
