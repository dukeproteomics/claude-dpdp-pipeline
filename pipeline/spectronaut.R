# Spectronaut precursor exports: which values are measured and which were filled by Spectronaut.
# Accepted (user decision 2026-11-04): raw exports (no background, no imputation) and background-signal exports, one TSV each.
#   raw:        empty / "Filtered" = missing; every value present is measured.
#   background: a value with run-wise EG.Qvalue > q_cut (or no q-value) is Spectronaut background signal (filled); q <= q_cut is measured.
#               Gaps with no value at all stay missing for the pipeline to fill.
#   global / run-wise imputation (EG.IsImputed = True anywhere): refused, the flag also marks measured values (11488: 79,547).
# Without run-wise EG.Qvalue columns the two accepted types can't be told apart, so the analyst must say which it is, and a
# background export is refused (its filled values could not be tracked).

sp_run_prefix <- function(cols, what) sub(paste0("\\.EG\\.", what, ".*$"), "", cols)

# df: the precursor table as exported (precursor_raw). Returns a list with mode, mask (TRUE = filled by Spectronaut), counts, notes.
sp_detect_fill <- function(df, q_cut = 0.01, declared = NULL) {
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  qcols <- grep("EG\\.TotalQuantity", names(df), value = TRUE)
  vcols <- grep("\\.EG\\.Qvalue$", names(df), value = TRUE)
  icols <- grep("\\.EG\\.IsImputed$", names(df), value = TRUE)
  Q <- sapply(df[qcols], num); Q[!is.na(Q) & Q <= 0] <- NA
  present <- !is.na(Q)
  out <- list(q_cut = q_cut, has_qvalue = length(vcols) > 0, has_isimputed = length(icols) > 0, declared = declared %||% "",
              values_present = sum(present), gaps = sum(!present), notes = character(0))
  if (length(icols)) {
    imp <- sapply(df[icols], function(x) toupper(as.character(x)) == "TRUE")
    out$isimputed_true <- sum(imp, na.rm = TRUE)
    if (out$isimputed_true > 0) {
      out$mode <- "imputed"
      out$error <- sprintf("This export has Spectronaut imputation (EG.IsImputed = True for %s values). Only raw or background-signal exports are accepted: re-export with imputation off (raw) or with background signal.",
                           format(out$isimputed_true, big.mark = ","))
      return(out)
    }
  }
  if (!length(vcols)) {
    if (is.null(declared)) {
      out$mode <- "unknown"
      out$error <- "This export has no run-wise EG.Qvalue columns, so raw and background exports can't be told apart. Re-export with EG.Qvalue, or, if it is a raw export (no background, no imputation), load it with --sp_fill=raw."
      return(out)
    }
    if (declared == "background") {
      out$mode <- "background"
      out$error <- "A background export needs the run-wise EG.Qvalue columns to tell measured from filled values. Re-export with EG.Qvalue."
      return(out)
    }
    out$mode <- "raw"; out$mask <- matrix(FALSE, nrow(Q), ncol(Q)); out$filled <- 0L
    out$notes <- c(out$notes, "No run-wise q-values in the export: declared raw by the analyst, so every value present is treated as measured.")
    return(out)
  }
  key_q <- sp_run_prefix(qcols, "TotalQuantity"); key_v <- sp_run_prefix(vcols, "Qvalue")
  if (!setequal(key_q, key_v)) stop("EG.Qvalue and EG.TotalQuantity columns name different runs.")
  V <- sapply(df[vcols[match(key_q, key_v)]], num)
  mask <- present & (is.na(V) | V > q_cut)
  out$filled <- sum(mask); out$measured <- sum(present & !mask)
  out$identified_but_empty <- sum(!present & !is.na(V) & V <= q_cut)
  out$mode <- if (out$filled > 0) "background" else "raw"
  if (!is.null(declared) && declared != out$mode) {
    out$error <- sprintf("You said this is a %s export, but the q-values say %s (%s values present with q > %s).", declared, out$mode,
                         format(out$filled, big.mark = ","), q_cut)
    return(out)
  }
  out$mask <- mask
  out$filled_per_sample <- colSums(mask)
  out
}
