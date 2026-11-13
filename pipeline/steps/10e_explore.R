# Step 10e - Explore: search any protein (or peptide) and plot it, with the peptides of its protein (user, 2026-11-13).
# The app's Stats > Protein Plot and Peptide Plot tabs (one protein's value per sample; its peptides, mean per group, as intensity or
# z-score), for any row of the tested data, not only the significant ones. The page searches and draws in the browser.
#
# Usage: Rscript pipeline/steps/10e_explore.R --db=<db>
# Run after Stats (10_stats apply); uses the tested normalization (stat_norm) and the comparisons (stats_comp) for the sample sets and
# the FC / p shown next to each hit. Protein output: the protein table and the peptide table of that normalization. Peptide output
# (phospho): the tested peptide table; a protein's peptides are its rows.
#
# Output: review/10e_explore_data.json.gz (the data, loaded by the page; log2 x 100 as integers), steps/10_stats/explore.json,
#         review/10e_explore.html

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "reviews", "10e_explore.R"))
args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app("Shiny_File.R")

p <- get_params(db_path)
if (!table_exists("stats_comp", db_path)) stop("Run Stats (10_stats.R --action=apply) first.")
sc <- read_table("stats_comp", db_path)
design <- read_table("design", db_path)
n_s <- as.integer(p$sample_number)
pep_out <- is_peptide_output(p)
log2i <- function(m) { m <- as.matrix(m); storage.mode(m) <- "numeric"; m[!is.na(m) & m <= 0] <- NA; v <- round(log2(m) * 100); storage.mode(v) <- "integer"; v }
samples_of <- function(df) df[, (ncol(df) - n_s + 1):ncol(df), drop = FALSE]
clean_seq <- function(x) gsub("_", "", x %||% "")
locs0 <- function(x) if (is.na(x) || !nzchar(x)) integer(0) else as.integer(trimws(strsplit(as.character(x), ",")[[1]])) - 1L

tab <- read_table(base_table(p), db_path)
main <- list(acc = tab$Accession, gene = ifelse(is.na(tab$Genes), "", tab$Genes),
             desc = substr(sub(" OS=.*$", "", ifelse(is.na(tab$Description), "", tab$Description)), 1, 90), v = log2i(samples_of(tab)))
if (pep_out) { main$seq <- clean_seq(tab$Sequence); main$site <- if (!is.null(tab$Protein_PTM_Loc)) ifelse(is.na(tab$Protein_PTM_Loc), "", tab$Protein_PTM_Loc) else rep("", nrow(tab))
  main$pos <- suppressWarnings(as.integer(sub("[^0-9].*", "", tab$PeptidePosition %||% rep(NA, nrow(tab))))) }
key_main <- if (pep_out) paste(tab$Accession, tab$Sequence, tab$Protein_PTM_Loc %||% "") else tab$Accession

# Protein output: the peptides of each protein (whole-dataset peptide table of the tested normalization, as 10c)
pep <- NULL
if (!pep_out) {
  pt <- paste0("peptide_", trimws(p$stat_norm))
  if (table_exists(pt, db_path)) {
    pd <- read_table(pt, db_path)
    k <- match(pd$Accession, tab$Accession); keep <- !is.na(k)
    pd <- pd[keep, ]; k <- k[keep]
    pep <- list(prot = k - 1L, seq = clean_seq(pd$Sequence),
                pos = suppressWarnings(as.integer(sub("[^0-9].*", "", pd$PeptidePosition %||% rep(NA, nrow(pd))))), v = log2i(samples_of(pd)))
  }
}

# FC / p / adj p / call per comparison, aligned to the main rows (null where the row wasn't tested, e.g. the missing factor)
comps <- lapply(seq_len(nrow(sc)), function(i) {
  nm <- sc$Name[i]; ft <- final_table(sc, i, p)
  st <- if (table_exists(ft, db_path)) read_table(ft, db_path) else NULL
  m <- if (!is.null(st)) match(key_main, if (pep_out) paste(st$Accession, st$Sequence, st$Protein_PTM_Loc %||% "") else st$Accession) else rep(NA_integer_, nrow(tab))
  g <- function(s) { v <- if (!is.null(st)) comp_get(st, nm, s) else NULL; if (is.null(v)) rep(NA, nrow(tab)) else v[m] }
  list(name = nm, N = locs0(sc$N_loc[i]), D = locs0(sc$D_loc[i]), P = locs0(sc$SPQC_loc[i]), Nname = sc$FactorsN[i], Dname = sc$FactorsD[i],
       fc = signif(g("_FC"), 4), p = signif(g("_pval"), 3), adj = signif(g("_adjpval"), 3), call = if (!is.null(st)) st$Stats[m] else rep(NA, nrow(tab)))
})

# Each comparison's full table, as its Excel sheet (the app's stats_Final_Excel: the _final table without the _FC2 columns;
# "Normalized" -> "Imputed" for the impute normalization), one file per comparison, loaded when the page shows that comparison.
# row = the main row (0-based) each table row belongs to, so a click plots it.
tables <- vapply(seq_len(nrow(sc)), function(i) {
  ft <- final_table(sc, i, p); if (!table_exists(ft, db_path)) return("")
  st <- read_table(ft, db_path); st <- st[, !grepl("_FC2", names(st)), drop = FALSE]
  if (trimws(p$stat_norm) == "impute") names(st) <- gsub("Normalized", "Imputed", names(st))
  st_names <- display_names(names(st), design, sc)   # as written in the design (no "X55Ex.152025...")
  k <- match(if (pep_out) paste(st$Accession, st$Sequence, st$Protein_PTM_Loc %||% "") else st$Accession, key_main) - 1L
  cols <- lapply(st, function(v) if (is.numeric(v)) signif(v, 6) else as.character(v))
  fn <- sprintf("10e_explore_table_%d.json.gz", i)
  con <- gzfile(file.path(review_dir(run_dir), fn), "wb", compression = 6)
  writeLines(toJSON(list(name = sc$Name[i], columns = st_names, numeric = unname(vapply(st, is.numeric, TRUE)), row = k, data = unname(cols)),
                    auto_unbox = TRUE, digits = NA, na = "null", null = "null"), con); close(con)
  fn }, "")
for (i in seq_along(comps)) comps[[i]]$table <- tables[i]

data <- list(run = basename(run_dir), level = if (pep_out) "peptide" else "protein", norm = trimws(p$stat_norm), spqc = p$comp_spqc,
             samples = design$Label, ids = design$ID, groups = design$Group, order = unique(design$Group), comps = comps, rows = main, pep = pep)
f <- file.path(review_dir(run_dir), "10e_explore_data.json.gz")
con <- gzfile(f, "wb", compression = 6); writeLines(toJSON(data, auto_unbox = TRUE, digits = NA, na = "null", null = "null"), con); close(con)
info <- list(step = "10e_explore", rows = nrow(tab), peptides = if (!is.null(pep)) length(pep$seq) else nrow(tab), level = data$level,
             norm = data$norm, comparisons = sc$Name, data_file = f, data_mb = round(file.size(f) / 1e6, 1))
write_json(c(info, generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file.path(step_dir(run_dir, "10_stats"), "explore.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Explore data: %d %ss, %d peptides, %s MB\n", info$rows, data$level, info$peptides, info$data_mb))
cat("Review page:", review_explore(db_path), "\n")
