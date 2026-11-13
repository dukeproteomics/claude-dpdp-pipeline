# Going back to an earlier step and redoing the later ones (user, 2026-11-13). Used by pipeline/revisit.R and pipeline/replay.R.
#
# How a step is redone with the analyst's earlier settings: each step reads them from the backup .db (--settings_from, the
# params as they were when the analyst went back), so nothing has to be retyped and nothing is chosen by Claude.
# A step the analyst had applied before is previewed and applied; the step they were on (not applied yet) is only previewed.
# Stats results are never accepted automatically: they changed, so the analyst looks at them again.

# The commands that redo one step. b_db: the backup .db; b_dec: decisions from the backup; b_dir: the backup folder.
# Each command: list(script, args, apply = TRUE if it applies a decision).
redo_cmds <- function(step, db, b_db, b_dec, b_dir, applied = TRUE) {
  sf <- paste0("--settings_from=", b_db)
  st <- function(s) file.path("pipeline", "steps", s)
  pre <- function(s, ...) list(script = st(s), args = c(paste0("--db=", db), "--action=preview", sf, ...), apply = FALSE)
  app <- function(s, ...) list(script = st(s), args = c(paste0("--db=", db), "--action=apply", sf, ...), apply = TRUE)
  cmds <- switch(step,
    "03_parameters" = list(list(script = st("03_parameters.R"), args = c(paste0("--db=", db), "--force", sf), apply = FALSE),
                           list(script = "pipeline/accept.R", args = c(paste0("--db=", db), "--step=03_parameters"), apply = TRUE)),
    "04_noise" = list(pre("04_noise.R"), app("04_noise.R")),
    "05_filter" = list(pre("05_filter.R"), app("05_filter.R")),
    "06_normalize" = list(pre("06_normalize.R"), app("06_normalize.R")),
    "07_impute" = list(pre("07_impute.R"), app("07_impute.R")),
    "08_rollup" = list(pre("08_rollup.R"), app("08_rollup.R")),
    "09_qc" = {
      qs <- tryCatch(fromJSON(file.path(b_dir, "steps", "09_qc", "summary.json"), simplifyVector = FALSE)$settings, error = function(e) NULL)
      qa <- function(k) { v <- unlist(qs[[k]]); if (length(v) && !identical(v, "0")) paste0("--", k, "=", paste(v, collapse = ",")) }
      acc <- Filter(function(r) identical(r$step, "09_qc") && isTRUE(r$settings$accepted), b_dec)
      fl <- if (length(acc)) unlist(acc[[length(acc)]]$settings$flagged_samples) else NULL
      list(list(script = st("09_qc.R"), args = c(paste0("--db=", db), qa("spike_accession"), qa("constant_accession"), qa("protein_accession")), apply = FALSE),
           list(script = st("09_qc.R"), args = c(paste0("--db=", db), "--action=accept", if (length(fl)) paste0("--flag_samples=", paste(fl, collapse = ","))), apply = TRUE))
    },
    "10_stats" = c(list(pre("10_stats.R"), app("10_stats.R")),
                   if (file.exists(file.path(b_dir, "steps", "10_stats", "graphs.json"))) list(list(script = st("10b_stats_graphs.R"), args = paste0("--db=", db), apply = FALSE)),
                   if (file.exists(file.path(b_dir, "steps", "10_stats", "output.json"))) list(list(script = st("10c_stats_output.R"), args = paste0("--db=", db), apply = FALSE)),
                   if (file.exists(file.path(b_dir, "steps", "10_stats", "explore.json"))) list(list(script = st("10e_explore.R"), args = paste0("--db=", db), apply = FALSE))),
    "11_pathway" = list(pre("11_pathway.R"), app("11_pathway.R")),
    stop("No redo commands for ", step))
  if (!applied) cmds <- cmds[1]   # the step the analyst was on: preview only
  cmds
}

# Did the analyst apply (or accept) this step? decisions: from read_decisions.
step_was_applied <- function(decisions, step) {
  if (step == "11_pathway") return(length(Filter(function(r) identical(r$step, step) && isTRUE(r$settings$accepted), decisions)) > 0)
  if (step %in% c("03_parameters", "09_qc", "10_stats") &&
      length(Filter(function(r) identical(r$step, step) && isTRUE(r$settings$accepted), decisions))) return(TRUE)
  step != "03_parameters" && step != "09_qc" && !is.null(applied_settings(decisions, step))
}

# Run one command through run_step.sh (stderr kept in <run>/Error/). Extra args: analyst, note, decided_by for applies.
run_redo_cmd <- function(cmd, root, extra = character(0)) {
  a <- c(cmd$args, if (cmd$apply) extra)
  cat("\n>>", cmd$script, paste(a, collapse = " "), "\n")
  status <- system2(file.path(root, "pipeline", "run_step.sh"), c(file.path(root, cmd$script), shQuote(a)))
  invisible(status)
}

# Write <run>/revisit.json.
write_revisit <- function(run_dir, rv) write_json(rv, file.path(run_dir, "revisit.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")

# A record in decisions.jsonl that doesn't count as a step decision (step "revisit"), so it doesn't mark steps out of date.
log_revisit_event <- function(run_dir, settings, note = "", analyst = NULL) {
  rec <- list(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step = "revisit", analyst = analyst %||% Sys.getenv("USER"),
              settings = settings, note = note, source = "pipeline")
  cat(toJSON(rec, auto_unbox = TRUE, null = "null"), "\n", file = file.path(run_dir, "decisions.jsonl"), append = TRUE, sep = "")
}

# Take the page's pending decision for a step (archive it); returns it or NULL.
take_page_decision <- function(run_dir, step, choice) {
  pend <- file.path(run_dir, "review", "pending_decision.json")
  if (!file.exists(pend)) stop("--from_page: no decision from the review page (", pend, ")", call. = FALSE)
  pd <- fromJSON(pend, simplifyVector = FALSE)
  if (!identical(pd$step, step) || !identical(pd$choice, choice))
    stop(sprintf("--from_page: the page asked for '%s' on %s, not '%s' on %s", pd$choice, pd$step, choice, step), call. = FALSE)
  done <- file.path(run_dir, "review", "decisions_done"); dir.create(done, showWarnings = FALSE)
  file.rename(pend, file.path(done, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", step, "_", choice, ".json")))
  pd
}
