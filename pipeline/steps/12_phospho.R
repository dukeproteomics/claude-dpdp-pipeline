# Step 12 - Phospho motifs and kinases (pipeline; user 2026-11-13: better than the app's rmotifx / MoMo).
# For each comparison, the Up and the Down phosphosites against the tested sites:
#   Kinase Library kinase enrichment (pipeline/kinase_library_run.py, the kinase-library Python package, local),
#   enrichment logos (pLogo-style) and motif-x motifs (pipeline/phospho_motifs.R).
# Sequences: UniProt by accession (rest.uniprot.org; only accessions are sent), cached in steps/12_phospho/uniprot.fasta.
#
#   Rscript pipeline/steps/12_phospho.R --db=<db> --action=preview [--from_page] [--phospho_<option>=...]
#   Rscript pipeline/steps/12_phospho.R --db=<db> --action=apply --decided_by=user [--from_page]     (accept the last run)
# Output: steps/12_phospho/summary.json, results.json, sites.tsv; Phospho/Phospho_motifs_<norm>.xlsx; review/12_phospho.html

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
root <- normalizePath(file.path(here, "..", ".."))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "phospho_motifs.R"))
source(file.path(here, "..", "options.R")); source(file.path(here, "..", "options", "12_phospho.R"))
source(file.path(here, "..", "reviews", "12_phospho.R"))
args <- parse_args()
need(args, c("db", "action"))
if (args$action == "run") args$action <- "preview"
if (!args$action %in% c("preview", "apply")) stop("--action must be preview (run) or apply (accept)")
if (args$action == "apply" && !identical(args$decided_by, "user")) stop("Accepting needs --decided_by=user (after the analyst has reviewed the page).")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app("Shiny_File.R")
p <- get_params(db_path)
if (!(is_ptm_run(p) && is_peptide_output(p))) stop("Phospho motifs need a phospho / PTM run with Peptide output (Parameters).")
if (!table_exists("stats_comp", db_path)) stop("Run Stats (Step 10) first.")
sc <- read_table("stats_comp", db_path)
sd <- step_dir(run_dir, "12_phospho")
opts <- phospho_options(p); spec <- opts_spec(opts)
args <- page_args(args, run_dir, "12_phospho", opts, args$action)
settings <- resolve_settings(spec, args)
if (args$action == "apply") settings <- apply_settings(run_dir, "12_phospho", spec, args, settings, opts)
opts_validate(opts, settings)
if (args$action == "apply") {
  log_decision(run_dir, "12_phospho", c(list(accepted = TRUE), settings), note = args$note %||% "", analyst = args$analyst)
  cat("Phospho results accepted\nReview page:", review_phospho(db_path), "\n"); quit(save = "no")
}
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))
py <- kl_python(root)
if (is.na(py)) stop("The kinase-library Python package was not found (run sudo bash web/install.sh, or set DPDP_KL_PYTHON).")

# sites per comparison
sites <- lapply(seq_len(nrow(sc)), function(i) {
  st <- read_table(final_table(sc, i, p), db_path)
  d <- phospho_sites(st, sc$Name[i], settings$phospho_hits, settings$phospho_p_type, settings$phospho_p_cutoff, settings$phospho_fc_cutoff,
                     as_flag(settings$phospho_localized))
  d$comp <- sc$Name[i]; d })
all_acc <- unique(unlist(lapply(sites, `[[`, "acc")))
cat(sprintf("Sites: %s; proteins: %d; fetching sequences from UniProt...\n", paste(vapply(sites, nrow, 1L), collapse = " / "), length(all_acc)))
seqs <- uniprot_sequences(all_acc, file.path(sd, "uniprot.fasta"))
sites <- lapply(sites, function(d) { d$window <- site_window(unname(seqs[d$acc]), d$pos, d$residue); d })
miss <- sum(vapply(sites, function(d) sum(is.na(d$window)), 1L)); tot <- sum(vapply(sites, nrow, 1L))
cat(sprintf("Windows: %d of %d sites (%d without a matching UniProt sequence)\n", tot - miss, tot, miss))

# Kinase Library (Python)
tab <- do.call(rbind, lapply(sites, function(d) { d <- d[!is.na(d$window), ]
  data.frame(comp = d$comp, site = d$key, window = d$window, label = ifelse(d$call == "Up", 1L, ifelse(d$call == "Down", -1L, 0L))) }))
f_in <- file.path(sd, "sites.tsv"); f_out <- file.path(sd, "kinase_library.json")
write.table(tab, f_in, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Kinase Library...\n")
st_ <- system2(py, c(file.path(here, "..", "kinase_library_run.py"), "--sites", f_in, "--out", f_out, "--rank", settings$phospho_kl_rank), stdout = TRUE, stderr = TRUE)
if (!file.exists(f_out)) stop("Kinase Library failed:\n", paste(tail(st_, 15), collapse = "\n"))
KL <- fromJSON(f_out, simplifyVector = FALSE)
cat(tail(st_, 1), "\n")

# logos and motif-x per comparison and direction
results <- lapply(sites, function(d) {
  nm <- d$comp[1]; w <- d[!is.na(d$window), ]
  bg <- w$window
  per_dir <- lapply(c(up = "Up", down = "Down"), function(dirn) {
    fg <- w$window[w$call == dirn]
    st_fg <- fg[toupper(substr(fg, 8, 8)) %in% c("S", "T")]; st_bg <- bg[toupper(substr(bg, 8, 8)) %in% c("S", "T")]
    motifs <- do.call(rbind, lapply(c("S", "T", "Y"), function(cr) { m <- motif_x(fg, bg, cr, settings$phospho_motif_min, settings$phospho_motif_sig); if (NROW(m)) cbind(central = cr, m) else NULL }))
    list(n = sum(w$call == dirn), logo = logo_stats(st_fg, st_bg), motifs = if (is.null(motifs)) list() else motifs)
  })
  kl <- KL$comparisons[[nm]]
  cat(sprintf("== %s: %d sites (%d Up, %d Down); motifs Up %d, Down %d\n", nm, nrow(w), per_dir$up$n, per_dir$down$n, NROW(per_dir$up$motifs), NROW(per_dir$down$motifs)))
  top <- KL$top_kinases
  w$top_kinases <- vapply(toupper(w$window), function(x) { t <- top[[x]]; if (is.null(t)) "" else paste(vapply(t, function(k) sprintf("%s (%s)", k[[1]], k[[2]]), ""), collapse = ", ") }, "")
  list(name = nm, sites = nrow(d), with_window = nrow(w), up = per_dir$up, down = per_dir$down, kinases = kl,
       site_table = w[order(w$call == "", w$p), c("acc", "gene", "site", "window", "fc", "p", "call", "peptides", "top_kinases")])
})

# Excel
xdir <- file.path(run_dir, "Phospho"); dir.create(xdir, showWarnings = FALSE)
xf <- file.path(xdir, sprintf("Phospho_motifs_%s.xlsx", trimws(p$stat_norm)))
wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "About")
openxlsx::writeData(wb, "About", data.frame(Item = c("Run", "Kinase Library", "Method", "Sequences", "Background", "Hits", "motif-x", "Made"),
  Value = c(basename(run_dir), paste("kinase-library", KL$version), sprintf("percentile rank, top %s kinases per site; Fisher's exact (BH), Up or Down v unchanged sites", KL$rank),
            "UniProt (first accession), 15-residue windows", "tested sites of the comparison (logos, motif-x); unchanged sites (Kinase Library)",
            if (identical(settings$phospho_hits, "custom")) sprintf("custom: %s p <= %s, |FC| >= %s", settings$phospho_p_type, settings$phospho_p_cutoff, settings$phospho_fc_cutoff) else "Stats calls",
            sprintf("min %s occurrences, p < %s", settings$phospho_motif_min, settings$phospho_motif_sig), format(Sys.time(), "%Y-%m-%d %H:%M"))))
for (i in seq_along(results)) { r <- results[[i]]
  for (kt in c("ser_thr", "tyrosine")) for (dirn in c("up", "down")) { t <- r$kinases[[kt]][[dirn]]
    if (length(t)) { sh <- sprintf("C%d %s kin %s", i, if (dirn == "up") "Up" else "Down", if (kt == "ser_thr") "ST" else "Y"); openxlsx::addWorksheet(wb, sh)
      openxlsx::writeDataTable(wb, sh, do.call(rbind, lapply(t, as.data.frame)), tableStyle = "TableStyleLight9") } }
  for (dirn in c("up", "down")) if (NROW(r[[dirn]]$motifs)) { sh <- sprintf("C%d %s motifs", i, if (dirn == "up") "Up" else "Down"); openxlsx::addWorksheet(wb, sh); openxlsx::writeDataTable(wb, sh, r[[dirn]]$motifs, tableStyle = "TableStyleLight9") }
  sh <- sprintf("C%d sites", i); openxlsx::addWorksheet(wb, sh); openxlsx::writeDataTable(wb, sh, r$site_table, tableStyle = "TableStyleLight9")
}
openxlsx::saveWorkbook(wb, xf, overwrite = TRUE)

summary <- list(step = "12_phospho", action = "preview", db_path = db_path, settings = settings, kinase_library = KL$version, excel = xf,
                windows = list(sites = tot, with_window = tot - miss),
                comparisons = lapply(results, function(r) { top <- function(t) if (length(t)) head(vapply(Filter(function(x) x$fisher_adj_pval <= settings$phospho_kl_fdr, t), function(x) sprintf("%s (%s)", x$kinase, signif(x$fisher_adj_pval, 2)), ""), 8) else character(0)
                  list(name = r$name, sites = r$sites, with_window = r$with_window, up = r$up$n, down = r$down$n,
                       up_kinases = top(r$kinases$ser_thr$up), down_kinases = top(r$kinases$ser_thr$down),
                       up_motifs = if (NROW(r$up$motifs)) head(r$up$motifs$motif, 6) else list(), down_motifs = if (NROW(r$down$motifs)) head(r$down$motifs$motif, 6) else list()) }))
write_summary(summary, sd)
# the page: everything but the site tables, which go in one gzip file per comparison (loaded when shown)
for (i in seq_along(results)) {
  con <- gzfile(file.path(review_dir(run_dir), sprintf("12_phospho_sites_%d.json.gz", i)), "wb", compression = 6)
  writeLines(toJSON(results[[i]]$site_table, digits = NA, na = "null"), con); close(con)
}
write_json(lapply(seq_along(results), function(i) { r <- results[[i]]; r$site_table <- NULL; r$site_file <- sprintf("12_phospho_sites_%d.json.gz", i); r }),
           file.path(sd, "results.json"), auto_unbox = TRUE, digits = NA, na = "null", null = "null")
cat("\nExcel:", xf, "\nReview page:", review_phospho(db_path), "\n>>> Waiting for the analyst's review.\n")
