# Spectronaut background-signal exports in Steps 3-7 (user rule, 2026-11-04): background values are signal, not identifications.
#   - Quantification: kept. They go back into the data after the filter and stay through Impute, rollup and stats.
#   - Normalization factors: measured values only (background excluded); the factors are then applied to everything.
#   - Noise and filter decisions: background counts as not detected, like an empty value.
#   - Detected/imputed tracking (precursor_missing, protein_missing, so the stats missing factor): background counts as not detected.
#
# Tables (background runs only):
#   precursor_start_bg       precursor_start with the background values (Step 3); precursor_start itself is measured only.
#   precursor_start_bgmask   1 = background value, same rows/columns as precursor_start.
#   precursor_filter_measured  precursor_filter with background blanked (Step 5), used for normalization factors.
#   precursor_filter_bgmask  1 = background value, same rows/columns as precursor_filter.

is_background_run <- function(db_path) {   # runs loaded before 2026-11-04 have no sp_fill_mode param
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path); on.exit(RSQLite::dbDisconnect(con))
  "sp_fill_mode" %in% RSQLite::dbListFields(con, "params") &&
    identical(RSQLite::dbGetQuery(con, "SELECT sp_fill_mode FROM params")$sp_fill_mode[1], "background")
}

samples_of <- function(df, n_s) as.matrix(df[(ncol(df) - n_s + 1):ncol(df)])
row_key <- function(d, info) do.call(paste, c(lapply(d[seq_len(info)], as.character), sep = "\r"))

# Step 3, right after precursor_start is made: move the background values out of precursor_start.
# precursor_fill has the export's row order and run columns; precursor_start may be reordered (rows by PrecursorId, columns by the
# design), so columns are matched by their values.
bg_split <- function(db_path) {
  n_s <- as.integer(get_param("sample_number", db_path))
  st <- read_table("precursor_start", db_path); fill <- read_table("precursor_fill", db_path); raw <- read_table("precursor_raw", db_path)
  if (anyDuplicated(st$PrecursorId) || anyDuplicated(fill$PrecursorId)) stop("bg_split: PrecursorId is not unique; can't map background values.")
  ri <- match(st$PrecursorId, fill$PrecursorId)
  if (anyNA(ri)) stop("bg_split: ", sum(is.na(ri)), " precursors in precursor_start are not in precursor_fill.")
  S <- samples_of(st, n_s)
  qc <- grep("EG.TotalQuantity", names(raw), value = TRUE)
  R <- sapply(raw[ri, qc], function(x) suppressWarnings(as.numeric(as.character(x)))); R[!is.na(R) & R <= 0] <- NA
  Fm <- as.matrix(fill[ri, -1]) == 1
  # column j of precursor_start = the raw column with the same values
  cj <- vapply(seq_len(n_s), function(j) { s <- S[, j]
    hits <- which(vapply(seq_len(ncol(R)), function(k) { r <- R[, k]; ok <- !is.na(s) | !is.na(r)
      all(is.na(s[ok]) == is.na(r[ok])) && isTRUE(all.equal(s[ok & !is.na(s)], r[ok & !is.na(s)], tolerance = 1e-6)) }, TRUE))
    if (length(hits) != 1) stop("bg_split: sample column ", colnames(S)[j], " matches ", length(hits), " export columns.")
    hits }, 1L)
  M <- Fm[, cj, drop = FALSE]
  write_table_try("precursor_start_bg", st, db_path)
  S[M] <- NA
  st[(ncol(st) - n_s + 1):ncol(st)] <- as.data.frame(S)
  write_table_try("precursor_start", st, db_path)
  mk <- st; mk[(ncol(mk) - n_s + 1):ncol(mk)] <- as.data.frame(1L * M)
  write_table_try("precursor_start_bgmask", mk, db_path)
  list(background_values = sum(M), per_sample = as.list(setNames(as.integer(colSums(M)), colnames(S))),
       measured_values = sum(!is.na(S)))
}

# Step 5, after the filter (and noise_keep): put background values back for the precursors that passed.
# Not put back: cells the misaligned filter blanked (its own decision), and background values below the noise baseline when noise
# removal is on and noise_keep is off (the same rule measured values get).
bg_restore <- function(db_path) {
  p <- get_params(db_path); n_s <- as.integer(p$sample_number)
  f <- read_table("precursor_filter", db_path); bgs <- read_table("precursor_start_bg", db_path)
  mk <- read_table("precursor_start_bgmask", db_path); nz <- read_table("precursor_noise", db_path)
  info <- ncol(f) - n_s
  idx <- match(row_key(f, info), row_key(bgs, info))
  if (anyNA(idx)) stop("bg_restore: ", sum(is.na(idx)), " filtered rows not found in precursor_start_bg.")
  B <- samples_of(bgs, n_s)[idx, , drop = FALSE]; M <- samples_of(mk, n_s)[idx, , drop = FALSE] == 1
  N <- samples_of(nz, n_s)[idx, , drop = FALSE]; Fm <- samples_of(f, n_s)
  groups <- read_table("sample_groups", db_path)
  blanked <- matrix(FALSE, nrow(Fm), ncol(Fm))
  for (g in seq_len(nrow(groups))) { gc <- as.integer(groups$start[g]):as.integer(groups$end[g])
    rows <- rowSums(!is.na(Fm[, gc, drop = FALSE])) == 0 & rowSums(!is.na(N[, gc, drop = FALSE])) > 0
    blanked[rows, gc] <- TRUE }
  base <- switch(p$noise_type, dynamic = as.numeric(p$noise_inflection), fixed = as.numeric(p$noise_baseline_value), 0)
  below_noise <- if (base > 0 && !as_flag(p$noise_keep)) M & !is.na(B) & B < base else matrix(FALSE, nrow(B), ncol(B))
  put <- M & is.na(Fm) & !blanked & !below_noise
  meas <- f
  Fm[put] <- B[put]; f[(info + 1):ncol(f)] <- as.data.frame(Fm)
  write_table_try("precursor_filter", f, db_path)
  write_table_try("precursor_filter_measured", meas, db_path)
  fm <- f; fm[(info + 1):ncol(fm)] <- as.data.frame(1L * put); write_table_try("precursor_filter_bgmask", fm, db_path)
  rebuild_missing_tables(f, info, db_path)
  list(background_values_restored = sum(put), per_sample = as.list(setNames(as.integer(colSums(put)), names(f)[(info + 1):ncol(f)])),
       background_in_kept_rows = sum(M), left_out_misaligned = sum(M & blanked), left_out_below_noise = sum(below_noise & !blanked),
       empty_after = sum(is.na(Fm)))
}

# The Impute page's missing-value tables, as filter_data_bg() builds them.
rebuild_missing_tables <- function(f, info, db_path) {
  df_samples <- f[(info + 1):ncol(f)]
  mv <- df_samples |> tidyr::gather(key = "key", value = "val") |> dplyr::mutate(is.missing = is.na(val)) |>
    dplyr::group_by(key, is.missing) |> dplyr::summarise(num.missing = dplyr::n(), .groups = "drop_last") |>
    dplyr::filter(is.missing == TRUE) |> dplyr::select(-is.missing) |> dplyr::arrange(desc(num.missing))
  write_table_try("missing_values", mv, db_path)
  mv2 <- df_samples |> tidyr::gather(key = "key", value = "val") |> dplyr::mutate(isna = is.na(val)) |>
    dplyr::group_by(key) |> dplyr::mutate(total = dplyr::n()) |> dplyr::group_by(key, total, isna) |>
    dplyr::summarise(num.isna = dplyr::n(), .groups = "drop_last") |> dplyr::mutate(pct = num.isna / total * 100)
  write_table_try("missing_values_plots", mv2, db_path)
}

# Step 7: precursor_missing / protein_missing as create_imputed_df() makes them, from measured values only (background = 0).
bg_detected_tables <- function(db_path) {
  n_s <- as.integer(get_param("sample_number", db_path))
  df <- read_table("precursor_filter_measured", db_path); info <- ncol(df) - n_s
  x <- df[(info + 1):ncol(df)]; x[!is.na(x) & x > 0] <- 1; x[is.na(x)] <- 0
  prot <- df[c("Accession", "Description", "Name", "Genes")]
  dp <- cbind(prot, x) |> dplyr::group_by(Accession, Description, Name, Genes) |> dplyr::summarise_all(list(sum)) |> dplyr::ungroup() |> data.frame()
  write_table_try("precursor_missing", x, db_path)
  write_table_try("protein_missing", dp[5:ncol(dp)], db_path)
  invisible(TRUE)
}

# sltmm: the app computes the TMM factors after imputation, on the imputed table (so imputed values count). The pipeline computes
# them from the measured values of the SL-normalized table instead (norm exclude/include applied, background blanked), and applies
# them to the imputed table. Returns the factors.
tmm_measured_factors <- function(db_path, p) {
  n_s <- as.integer(p$sample_number)
  sl <- read_table("precursor_norm_sltmm", db_path); info <- ncol(sl) - n_s
  if (is_background_run(db_path)) {
    mk <- samples_of(read_table("precursor_filter_bgmask", db_path), n_s) == 1
    X <- samples_of(sl, n_s); X[mk] <- NA; sl[(info + 1):ncol(sl)] <- as.data.frame(X)
  }
  nd <- norm_filter_exclude_include(sl, p)
  X <- samples_of(nd, n_s); X[is.na(X)] <- 0
  edgeR::calcNormFactors(X, method = "TMM", sumTrim = 0.1)
}
