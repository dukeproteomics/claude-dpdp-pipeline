# Peptide refilter -> flags (user decision 2026-11-04).
# The app's peptide refilter (Stats) removed, per comparison, the precursors with no measured value in the comparison's samples (plus,
# when switched on, precursors failing the peptide missing-factor or CV filters), then rolled the proteins up again from those samples
# only. Re-rolling a subset gives different MaxLFQ values than the whole-dataset rollup, which confused results. The pipeline now always
# uses the whole-dataset protein values (no re-rollup) and FLAGS what the refilter would have removed, with the same rules.
#
#   add_detected_column(): protein_<norm> (and _final) get Detected_Imputed (measured precursors per sample, "n.n.n", "-" = 0), which
#     the app's Stats reads for the missing factor when it doesn't refilter. Background values don't count as measured.
#   refilter_flags(): for comparison i, per protein: precursors the refilter would remove, of all its precursors, and a flag
#     ("all precursors" = the refilter would have removed the protein; "half or more precursors" = at least half removed). On 11488 B vs A,
#     45 proteins lose at least one precursor; the 5 whose FC moved > 0.1 log2 all lose half or more (Q9JJ61: FC 4.8 whole-set vs 29 refiltered).

detected_string <- function(counts) apply(counts, 1, function(r) paste(ifelse(r == 0, "-", as.character(r)), collapse = "."))

add_detected_column <- function(db_path, norm) {
  n_s <- as.integer(get_param("sample_number", db_path))
  pf <- read_table("precursor_filter", db_path); pm <- as.matrix(read_table("precursor_missing", db_path))
  if (nrow(pm) != nrow(pf)) stop("add_detected_column: precursor_missing and precursor_filter have different rows")
  cnt <- rowsum(1L * (pm > 0), pf$Accession, reorder = FALSE)
  for (t in paste0("protein_", norm, c("", "_final"))) {
    if (!table_exists(t, db_path)) next
    d <- read_table(t, db_path)
    if ("Detected_Imputed" %in% names(d)) next
    di <- detected_string(cnt[match(d$Accession, rownames(cnt)), , drop = FALSE])
    pos <- which(names(d) == "Precursors"); if (!length(pos)) pos <- ncol(d) - n_s
    d <- tibble::add_column(d, Detected_Imputed = di, .after = pos)
    write_table_try(t, d, db_path)
  }
}

refilter_flags <- function(db_path, stats_comp, i, p) {
  loc <- function(x) as.integer(trimws(strsplit(as.character(x), ",")[[1]]))
  nl <- loc(stats_comp$N_loc[i]); dl <- loc(stats_comp$D_loc[i]); nm <- stats_comp$Name[i]; cols <- c(nl, dl)
  pf <- read_table("precursor_filter", db_path); pm <- as.matrix(read_table("precursor_missing", db_path)) > 0
  n_s <- as.integer(p$sample_number)
  rm_all <- rowSums(pm[, cols, drop = FALSE]) == 0                          # no measured value in the comparison's samples
  rm_mf <- if (as_flag(p$peptide_missing_filter)) {
    pmax(rowMeans(pm[, nl, drop = FALSE]), rowMeans(pm[, dl, drop = FALSE])) < as.numeric(p$peptide_missing_factor) } else FALSE
  rm_cv <- if (as_flag(p$peptide_cv_filter)) {
    x <- as.matrix(read_table(paste0("precursor_impute_", p$stat_norm), db_path)); x <- x[, (ncol(x) - n_s + 1):ncol(x)]; storage.mode(x) <- "numeric"
    cv <- function(cc) 100 * apply(x[, cc, drop = FALSE], 1, stats::sd, na.rm = TRUE) / rowMeans(x[, cc, drop = FALSE], na.rm = TRUE)
    pmin(cv(nl), cv(dl)) > as.numeric(p$peptide_cv_factor) * 100 } else FALSE
  rm <- rm_all | rm_mf | rm_cv; rm[is.na(rm)] <- FALSE
  pk <- feature_key(pf, p); n_rm <- tapply(rm, pk, sum); n_all <- tapply(rm, pk, length)
  tab <- final_table(stats_comp, i, p); df <- read_table(tab, db_path); a <- feature_key(df, p)
  r <- as.integer(n_rm[a]); n <- as.integer(n_all[a])
  flag <- ifelse(r == n, "refilter: all precursors", ifelse(r >= n / 2, "refilter: half or more precursors", ""))
  df[[comp_col(nm, "_refilter_precursors")]] <- sprintf("%d of %d", r, n)
  df[[comp_col(nm, "_refilter_flag")]] <- flag
  # keep Stats last, as in the app's layout
  df <- df[, c(setdiff(names(df), "Stats"), "Stats")]
  write_table_try(tab, df, db_path)
  sig <- df$Stats != ""
  list(rules = c("no measured value in the comparison's samples",
                 if (as_flag(p$peptide_missing_filter)) sprintf("peptide missing factor < %s", p$peptide_missing_factor),
                 if (as_flag(p$peptide_cv_filter)) sprintf("peptide CV > %s%% in both groups", as.numeric(p$peptide_cv_factor) * 100)),
       flagged_all = sum(flag == "refilter: all precursors"), flagged_most = sum(flag == "refilter: half or more precursors"), with_any_removed = sum(r > 0),
       flagged_all_significant = sum(sig & flag == "refilter: all precursors"), flagged_most_significant = sum(sig & flag == "refilter: half or more precursors"))
}

# Stats without refilter, run by the pipeline instead of the app's stat_calc_bg (whose no-refilter branch is broken twice: it takes the
# last sample_number columns of protein_<norm>_final, which end with CV columns, and its detected counts assume every sample is in the
# comparison, so it fails when other groups exist). Same core as the app: stat_create_comp_df + the app's stat_add (CVs, FC, p, adj p,
# missing factor filter, Up/Down), on the whole-dataset protein_<norm> table. Runs in callr with the app directory as working directory.
# Effect size and limma p for a comparison (pipeline, user 2026-11-13: "statistically and programmatically correct"). Replaces
# the app's versions, which took Cohen's d on raw intensities with the whole-sample SD (effsize pooled = FALSE) and ignored pairing.
#   Cohen's d on log2 (the scale of the t-test), numerator minus denominator, over the pooled SD of the two groups;
#   Hedges' g = d x J, J = 1 - 3 / (4 (nN + nD) - 9). Paired: d_z = mean of the per-pair log2 differences / their SD,
#   J = 1 - 3 / (4 (n - 1) - 1). NA where an SD is 0 or there are fewer than 2 values (pairs) per side.
#   limma p: limma on log2 (lmFit + eBayes over all rows), group N v D; paired: the pair as a blocking factor (~ pair + group).
# xN, xD: numeric matrices (rows = proteins/peptides, columns = samples in design order; paired: column k of xN with column k of xD).
stat_calc_norefilter_bg <- function(db_path, comp_number, stats_comp) {
  # run_bg runs this in a fresh R process: the pipeline helpers it needs are defined here.
  comp_col <- function(nm, suffix) make.names(paste0(nm, suffix))
  stat_effect_size <- function(xN, xD, paired = FALSE, hedges = FALSE) {
    a <- log2(as.matrix(xN)); b <- log2(as.matrix(xD)); a[!is.finite(a)] <- NA; b[!is.finite(b)] <- NA
    rv <- function(m) { k <- rowSums(!is.na(m)); mu <- rowMeans(m, na.rm = TRUE); v <- rowSums((m - mu)^2, na.rm = TRUE) / (k - 1); list(n = k, mu = mu, v = v) }
    if (paired) {
      dd <- rv(a - b); d <- dd$mu / sqrt(dd$v); J <- 1 - 3 / (4 * (dd$n - 1) - 1); ok <- dd$n >= 2
    } else {
      A <- rv(a); B <- rv(b); sp <- sqrt(((A$n - 1) * A$v + (B$n - 1) * B$v) / (A$n + B$n - 2))
      d <- (A$mu - B$mu) / sp; J <- 1 - 3 / (4 * (A$n + B$n) - 9); ok <- A$n >= 2 & B$n >= 2
    }
    if (isTRUE(hedges)) d <- d * J
    d[!ok | !is.finite(d)] <- NA
    signif(d, 4)
  }
  stat_limma_p <- function(xN, xD, paired = FALSE) {
    y <- log2(cbind(as.matrix(xN), as.matrix(xD))); y[!is.finite(y)] <- NA
    g <- factor(c(rep("N", ncol(xN)), rep("D", ncol(xD))), levels = c("D", "N"))
    des <- if (paired) model.matrix(~ factor(c(seq_len(ncol(xN)), seq_len(ncol(xD)))) + g) else model.matrix(~ g)
    fit <- limma::eBayes(limma::lmFit(y, des))
    signif(as.numeric(fit$p.value[, "gN"]), 7)
  }

  source("Shiny_MVA_Functions.R"); source("Shiny_Misc_Functions.R"); source("Shiny_File.R")
  params <- get_params(db_path); df_design <- read_table_try("design", db_path)
  base <- if (params$data_output == "Peptide") stringr::str_c("peptide_impute_", trimws(params$stat_norm)) else stringr::str_c("protein_", trimws(params$stat_norm))
  df <- read_table_try(base, db_path)
  if (!"Detected_Imputed" %in% names(df)) stop(base, " has no Detected_Imputed column (re-run Step 8).")
  cnt <- do.call(rbind, lapply(strsplit(as.character(df$Detected_Imputed), ".", fixed = TRUE), function(x) as.numeric(gsub("-", "0", x))))
  if (ncol(cnt) != params$sample_number) stop("Detected_Imputed has ", ncol(cnt), " counts per protein, expected ", params$sample_number)
  # PTM report (pipeline): drop the non-PTM peptides before testing, so adj p covers the PTM peptides only
  if (params$data_output == "Peptide" && as.logical(params$checkbox_report_ptm) && !params$peptide_report_grep %in% c("", "Enter value here")) {
    k <- grepl(params$peptide_report_grep, df$Sequence); df <- df[k, ]; cnt <- cnt[k, , drop = FALSE]
  }
  df_list <- stat_create_comp_df(df, stats_comp, comp_number, params, df_design)
  spqc <- as.numeric(stats_comp$SPQC[comp_number]) != 0
  df <- if (spqc) cbind(df_list[[1]], df_list[[2]], df_list[[3]], df_list[[4]]) else cbind(df_list[[1]], df_list[[2]], df_list[[3]])
  nl <- str_to_numlist(stats_comp$N_loc[comp_number]); dl <- str_to_numlist(stats_comp$D_loc[comp_number])
  df_missing <- as.data.frame(cnt[, c(nl, dl), drop = FALSE])
  # Cohen's d and limma p: the pipeline's (stat_effect_size / stat_limma_p above), on the same rows as the app's (before its filters)
  want_d <- as.logical(params$checkbox_cohensd); want_l <- as.logical(params$checkbox_limmapvalue)
  paired <- as.logical(params$pair_comp)
  extra <- list()
  if (want_d) extra$cohensd <- stat_effect_size(df_list[[2]], df_list[[3]], paired, as.logical(params$checkbox_cohensd_hedges))
  if (want_l) extra$limma <- stat_limma_p(df_list[[2]], df_list[[3]], paired)
  params$checkbox_cohensd <- FALSE; params$checkbox_limmapvalue <- FALSE   # the app's versions are not used
  df <- cbind(.row = seq_len(nrow(df)), df)   # first, not last: stat_add finds the sample columns counting from the end
  df <- stat_add(df, df_missing, params, comp_number, stats_comp, df_design)
  nm <- stats_comp$Name[comp_number]
  if (length(extra)) {   # put them where the app puts them: after _pval / _adjpval
    at <- max(which(names(df) %in% c(comp_col(nm, "_pval"), comp_col(nm, "_adjpval"))))
    add <- as.data.frame(lapply(extra, function(v) v[df$.row]))
    names(add) <- vapply(names(extra), function(k) comp_col(nm, paste0("_", k)), "")
    df <- cbind(df[, seq_len(at), drop = FALSE], add, df[, -seq_len(at), drop = FALSE])
  }
  df$.row <- NULL
  if (params$stat_norm == "impute") colnames(df) <- gsub("Normalized", "Imputed", colnames(df))
  write_table_try(if (params$data_output == "Peptide") stats_comp$Final_Table_Name_Peptide[comp_number] else stats_comp$Final_Table_Name[comp_number], df, db_path)
  invisible(TRUE)
}
