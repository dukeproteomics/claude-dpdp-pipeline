# noise_keep, as the user intends it (2026-11-04):
#   1. Noise step: every value below the noise baseline becomes NA (plain removal), so the filter step treats noise
#      as missing. A precursor that is noise in (almost) every sample then fails the filter and is dropped.
#   2. After the filter: for each precursor that survived, put the measured values that the noise step removed back
#      into the data, so they are not imputed later. (Removing a measured value and imputing noise in its place makes no sense.)
# The app's noise_remove_bg() noise_keep branch does something else (it cleans rows with >= noise_keep_value values
# above noise and leaves the others), so the pipeline doesn't use that branch.
#
# Values blanked by the misaligned filter (target = group/sample) stay NA: that is the filter's own decision.

noise_keep_restore <- function(db_path, start_table = "precursor_start", noise_table = "precursor_noise",
                               filter_table = "precursor_filter") {
  p <- get_params(db_path)
  n_s <- as.integer(p$sample_number)
  s <- read_table(start_table, db_path); n <- read_table(noise_table, db_path); f <- read_table(filter_table, db_path)
  info <- ncol(f) - n_s
  cols <- (info + 1):ncol(f)
  key <- function(d) do.call(paste, c(lapply(d[seq_len(info)], as.character), sep = "\r"))
  ks <- key(s)
  if (anyDuplicated(ks)) stop("noise_keep_restore: precursor rows in ", start_table, " are not unique by their info columns.")
  idx <- match(key(f), ks)
  if (anyNA(idx)) stop("noise_keep_restore: ", sum(is.na(idx)), " filtered rows not found in ", start_table)
  if (!identical(key(n), ks)) stop("noise_keep_restore: ", noise_table, " rows are not in the same order as ", start_table)

  S <- as.matrix(s[idx, (ncol(s) - n_s + 1):ncol(s)]); N <- as.matrix(n[idx, (ncol(n) - n_s + 1):ncol(n)])
  Fm <- as.matrix(f[, cols])
  noise_removed <- is.na(N) & !is.na(S) & S > 0
  # misaligned filter: a whole group blanked in a row where the noise table still had values
  blanked <- matrix(FALSE, nrow(Fm), ncol(Fm))
  groups <- read_table("sample_groups", db_path)
  for (g in seq_len(nrow(groups))) {
    gc <- as.integer(groups$start[g]):as.integer(groups$end[g])
    rows <- rowSums(!is.na(Fm[, gc, drop = FALSE])) == 0 & rowSums(!is.na(N[, gc, drop = FALSE]) & N[, gc, drop = FALSE] > 0) > 0
    blanked[rows, gc] <- TRUE
  }
  restore <- noise_removed & is.na(Fm) & !blanked
  Fm[restore] <- S[restore]
  f[, cols] <- Fm
  write_table_try(filter_table, f, db_path)

  # Rebuild the Impute-page missing-value tables the same way filter_data_bg() does.
  df_samples <- f[cols]
  mv <- df_samples |> tidyr::gather(key = "key", value = "val") |> dplyr::mutate(is.missing = is.na(val)) |>
    dplyr::group_by(key, is.missing) |> dplyr::summarise(num.missing = dplyr::n(), .groups = "drop_last") |>
    dplyr::filter(is.missing == TRUE) |> dplyr::select(-is.missing) |> dplyr::arrange(desc(num.missing))
  write_table_try("missing_values", mv, db_path)
  mv2 <- df_samples |> tidyr::gather(key = "key", value = "val") |> dplyr::mutate(isna = is.na(val)) |>
    dplyr::group_by(key) |> dplyr::mutate(total = dplyr::n()) |> dplyr::group_by(key, total, isna) |>
    dplyr::summarise(num.isna = dplyr::n(), .groups = "drop_last") |> dplyr::mutate(pct = num.isna / total * 100)
  write_table_try("missing_values_plots", mv2, db_path)

  list(values_restored = sum(restore),
       values_restored_per_sample = as.list(setNames(as.integer(colSums(restore)), names(f)[cols])),
       noise_values_in_kept_rows = sum(noise_removed),
       left_na_by_misaligned_filter = sum(noise_removed & blanked),
       precursors_with_restored_values = sum(rowSums(restore) > 0))
}
