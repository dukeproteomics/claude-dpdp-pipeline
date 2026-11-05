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
  n_rm <- tapply(rm, pf$Accession, sum); n_all <- tapply(rm, pf$Accession, length)
  tab <- stats_comp$Final_Table_Name[i]; df <- read_table(tab, db_path); a <- df$Accession
  r <- as.integer(n_rm[a]); n <- as.integer(n_all[a])
  flag <- ifelse(r == n, "refilter: all precursors", ifelse(r >= n / 2, "refilter: half or more precursors", ""))
  df[[paste0(nm, "_refilter_precursors")]] <- sprintf("%d of %d", r, n)
  df[[paste0(nm, "_refilter_flag")]] <- flag
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
stat_calc_norefilter_bg <- function(db_path, comp_number, stats_comp) {
  source("Shiny_MVA_Functions.R"); source("Shiny_Misc_Functions.R"); source("Shiny_File.R")
  params <- get_params(db_path); df_design <- read_table_try("design", db_path)
  df <- read_table_try(stringr::str_c("protein_", trimws(params$stat_norm)), db_path)
  if (!"Detected_Imputed" %in% names(df)) stop("protein_", params$stat_norm, " has no Detected_Imputed column (re-run Step 8).")
  cnt <- do.call(rbind, lapply(strsplit(as.character(df$Detected_Imputed), ".", fixed = TRUE), function(x) as.numeric(gsub("-", "0", x))))
  if (ncol(cnt) != params$sample_number) stop("Detected_Imputed has ", ncol(cnt), " counts per protein, expected ", params$sample_number)
  df_list <- stat_create_comp_df(df, stats_comp, comp_number, params, df_design)
  spqc <- as.numeric(stats_comp$SPQC[comp_number]) != 0
  df <- if (spqc) cbind(df_list[[1]], df_list[[2]], df_list[[3]], df_list[[4]]) else cbind(df_list[[1]], df_list[[2]], df_list[[3]])
  nl <- str_to_numlist(stats_comp$N_loc[comp_number]); dl <- str_to_numlist(stats_comp$D_loc[comp_number])
  df_missing <- as.data.frame(cnt[, c(nl, dl), drop = FALSE])
  df <- stat_add(df, df_missing, params, comp_number, stats_comp, df_design)
  if (params$stat_norm == "impute") colnames(df) <- gsub("Normalized", "Imputed", colnames(df))
  write_table_try(stats_comp$Final_Table_Name[comp_number], df, db_path)
  invisible(TRUE)
}
