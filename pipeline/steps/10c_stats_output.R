# Step 10c — Stats output: Excel files + Protein Plots
#   Excel: the app's stats_Final_Excel_bg -> <data_path>/<norm>/Final_<norm>_stats.xlsx and <prefix>_Parameters.xlsx
#          (Table1 Sample Info, Table2 Raw Precursor Data, Table3 Raw Peptide Data, Table4 Normalized Data, Table5.. comparisons)
#   Protein plots: like the app's Stats > Protein Plots (per-sample protein bar plot + per-peptide grouped bar plot).
#          The app's create_stats_oneprotein_plots_bg finds sample columns with grep(group name, column names), which
#          breaks for single-letter groups (A matches "Accession"), so this uses the design positions from stats_comp instead.
#
# Usage:
#   Rscript pipeline/steps/10c_stats_output.R --db=<db> [--excel=TRUE] [--final_stats_name=Final_<norm>_stats.xlsx]
#       [--comp=B_v_A] [--proteins=<acc>,<acc>] [--top=6] [--plot_spqc=TRUE]
#   --proteins  accessions to plot; default = the top N significant proteins by p-value (--top)
#
# Output: Excel files; QC/Protein_<comp>_<acc>_{bar,peptides}.png; steps/10_stats/output.json

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "contaminants.R"))
source(file.path(here, "..", "reviews", "10_stats.R"))
suppressPackageStartupMessages(library(ggplot2))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R", "Shiny_MVA.R", "Shiny_MVA_Functions.R"))

p <- get_params(db_path)
stats_comp <- read_table("stats_comp", db_path)
design <- read_table("design", db_path)
out <- list(step = "10c_stats_output", db_path = db_path)

# --- Excel (app: stats_Final_Excel) ---
if (is.null(args$excel) || as_flag(args$excel)) {
  norm <- trimws(p$stat_norm)
  final_name <- args$final_stats_name %||% str_c("Final_", norm, "_stats.xlsx")
  file_dir <- str_c(p$data_path, norm)
  filename <- str_c(p$data_path, norm, "//", final_name)
  filename_params <- str_c(p$data_path, norm, "//", p$file_prefix, "_Parameters.xlsx")
  if (!dir.exists(file_dir)) dir.create(file_dir, recursive = TRUE)
  run_bg(stats_Final_Excel_bg, list(file_dir, filename, filename_params, p), "finalexcel", db_path)
  fn <- gsub("/+", "/", c(filename, filename_params))
  # Pipeline addition: a "Flagged proteins" sheet (contaminants / spike-ins, flagged not removed). The app's sheets are unchanged.
  sp <- qc_spike_settings(run_dir)
  flag_rows <- do.call(rbind, lapply(seq_len(nrow(stats_comp)), function(i) {
    nm <- stats_comp$Name[i]; d <- read_table(final_table(stats_comp, i, p), db_path)
    fl <- flag_contaminants(d$Accession, d$Description, sp$spike, sp$constant); k <- fl != ""
    if (!any(k)) return(NULL)
    data.frame(Comparison = nm, Accession = d$Accession[k], Description = d$Description[k], Genes = d$Genes[k], Flag = fl[k],
               FC = d[[paste0(nm, "_FC")]][k], pval = d[[paste0(nm, "_pval")]][k], adjpval = d[[paste0(nm, "_adjpval")]][k] %||% NA,
               Stats = d$Stats[k], stringsAsFactors = FALSE)
  }))
  if (!is.null(flag_rows) && nrow(flag_rows)) {
    wb <- openxlsx::loadWorkbook(fn[1])
    if ("Flagged proteins" %in% names(wb)) openxlsx::removeWorksheet(wb, "Flagged proteins")
    openxlsx::addWorksheet(wb, "Flagged proteins")
    openxlsx::writeData(wb, "Flagged proteins", data.frame(Note = "Flagged, not removed: digestion trypsin, human keratins (non-human samples), QC spike-ins. Their rows stay in every sheet."), startRow = 1)
    openxlsx::writeDataTable(wb, "Flagged proteins", flag_rows[order(flag_rows$Stats == "", flag_rows$pval), ], startRow = 3, tableStyle = "TableStyleLight9")
    openxlsx::setColWidths(wb, "Flagged proteins", cols = 1:9, widths = c(12, 12, 50, 12, 24, 10, 12, 12, 8))
    openxlsx::saveWorkbook(wb, fn[1], overwrite = TRUE)
  }
  sheets <- openxlsx::getSheetNames(fn[1])
  sheet_rows <- setNames(lapply(sheets, function(s) { d <- openxlsx::read.xlsx(fn[1], sheet = s, rows = 1:2); if (is.null(d)) 0 else ncol(d) }), sheets)
  out$excel <- list(files = fn, sheets = sheets, sheet_columns = sheet_rows, flagged_rows = if (is.null(flag_rows)) 0 else nrow(flag_rows),
                    note = "Table1 Sample Info is created empty by the app (no data is written to it). 'Flagged proteins' is added by the pipeline.")
  cat("Excel:\n", paste0("  ", fn, collapse = "\n"), "\n  sheets:", paste(sheets, collapse = " | "), "\n")
}

# --- Protein plots ---
comp <- args$comp %||% stats_comp$Name[1]
ci <- which(stats_comp$Name == comp); if (!length(ci)) stop("Unknown --comp: ", comp)
plot_spqc <- if (!is.null(args$plot_spqc)) as_flag(args$plot_spqc) else TRUE
locs <- function(x) as.integer(trimws(strsplit(x, ",")[[1]]))
pos <- c(locs(stats_comp$N_loc[ci]), locs(stats_comp$D_loc[ci]), if (plot_spqc && p$comp_spqc != "0") locs(stats_comp$SPQC_loc[ci]))
grp_of <- function(i) if (i %in% locs(stats_comp$N_loc[ci])) stats_comp$FactorsN[ci] else if (i %in% locs(stats_comp$D_loc[ci])) stats_comp$FactorsD[ci] else p$comp_spqc
prot <- read_table(final_table(stats_comp, ci, p), db_path)
# per-comparison peptide table only exists with the old peptide refilter; otherwise the whole-dataset peptide table (Step 8)
pep_tab <- stats_comp$Final_Table_Name_Peptide[ci]
pep <- read_table(if (table_exists(pep_tab, db_path)) pep_tab else str_c("peptide_", p$stat_norm), db_path)
id_cols <- function(df) sapply(design$ID[pos], function(id) grep(str_c("\\.", id, "(\\.|$)"), names(df))[1])
pc <- id_cols(prot); qc <- id_cols(pep)
if (anyNA(pc) || anyNA(qc)) stop("Could not match sample columns by design ID.")

# Rows to plot: the top significant rows (proteins, or peptides for Peptide output), or every row of the given accessions.
pep_out <- is_peptide_output(p)
rows <- if (!is.null(args$proteins)) which(prot$Accession %in% trimws(strsplit(args$proteins, ",")[[1]])) else {
  s <- which(prot$Stats != ""); head(s[order(prot[[str_c(comp, "_pval")]][s])], as.integer(args$top %||% 6)) }
if (!is.null(args$proteins) && !pep_out) rows <- rows[!duplicated(prot$Accession[rows])]
labels <- design$Label[pos]; groups <- vapply(pos, grp_of, ""); cols <- design$colorlist[pos]
group_cols <- setNames(cols[!duplicated(groups)], groups[!duplicated(groups)])
plots <- character(0); prot_info <- list()
for (ri in rows) {
  r <- prot[ri, ]; a <- r$Accession
  tag <- if (pep_out) str_c(a, "_", ri) else a
  seq_lab <- if (pep_out) str_c(" ", gsub("_", "", r$Sequence), if (!is.null(r$Protein_PTM_Loc) && nzchar(r$Protein_PTM_Loc)) str_c(" ", r$Protein_PTM_Loc) else "") else ""
  d <- data.frame(Sample = factor(labels, levels = labels), Group = groups, Intensity = unlist(r[1, pc]))
  ttl <- sprintf("%s %s%s\nFC %.2f  p %.2g  adj p %.2g  %s", a, r$Genes, seq_lab, r[[str_c(comp, "_FC")]], r[[str_c(comp, "_pval")]],
                 r[[str_c(comp, "_adjpval")]] %||% NA, r$Stats)
  f1 <- gsub("/+", "/", file.path(p$qc_path, str_c(if (pep_out) "Peptide_" else "Protein_", comp, "_", tag, "_bar.png")))
  ggsave(f1, ggplot(d, aes(Sample, Intensity, fill = Group)) + geom_col() + scale_fill_manual(values = group_cols) +
           theme_classic() + ggtitle(ttl) + xlab(NULL) +
           theme(plot.title = element_text(hjust = 0.5, size = 10), axis.text.x = element_text(angle = 90)), width = 7, height = 5, dpi = 150)
  pr <- pep[pep$Accession == a, ]
  f2 <- NULL
  if (nrow(pr)) {
    lg <- log2(as.matrix(pr[, qc])); colnames(lg) <- labels
    long <- data.frame(Sequence = rep(pr$Sequence, times = length(pos)), Group = rep(groups, each = nrow(pr)), y = as.vector(lg))
    agg <- aggregate(y ~ Sequence + Group, long, function(v) c(mean = mean(v), sd = sd(v)))
    agg <- data.frame(Sequence = agg$Sequence, Group = factor(agg$Group, levels = unique(groups)), y_mean = agg$y[, "mean"], sd = agg$y[, "sd"])
    seq_order <- if ("PeptidePosition" %in% names(pr)) pr$Sequence[order(suppressWarnings(as.numeric(sub("[^0-9].*", "", pr$PeptidePosition))))] else pr$Sequence
    agg$Sequence <- factor(agg$Sequence, levels = unique(seq_order))
    f2 <- gsub("/+", "/", file.path(p$qc_path, str_c(if (pep_out) "Peptide_" else "Protein_", comp, "_", tag, "_peptides.png")))
    ggsave(f2, ggplot(agg, aes(Sequence, y_mean, fill = Group)) + geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
             geom_errorbar(aes(ymin = y_mean - sd, ymax = y_mean + sd), position = position_dodge(0.8), width = 0.25) +
             scale_fill_manual(values = group_cols) + theme_classic() + ylab("log2 intensity (mean ± SD)") + xlab(NULL) +
             ggtitle(str_c(a, " ", r$Genes, " - peptides (", nrow(pr), ")")) +
             theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, size = 6)),
           width = max(7, 0.25 * nrow(pr) + 3), height = 5, dpi = 150, limitsize = FALSE)
  }
  plots <- c(plots, f1, f2)
  prot_info[[tag]] <- list(accession = a, sequence = if (pep_out) r$Sequence else NULL, genes = r$Genes, name = r$Name, stats = r$Stats, FC = r[[str_c(comp, "_FC")]],
                         pval = r[[str_c(comp, "_pval")]], adjpval = r[[str_c(comp, "_adjpval")]], mf = r[[str_c(comp, "_mf")]],
                         detected_imputed = r$Detected_Imputed, precursors = r$Precursors, peptides = nrow(pr),
                         values = setNames(as.list(round(unlist(r[1, pc]))), labels))
}
out$protein_plots <- list(comp = comp, proteins = prot_info, plots = plots,
                          note = "Pipeline version of Stats > Protein Plots (the app's version fails for group names that match annotation columns).")
out$plots <- plots
sd <- step_dir(run_dir, "10_stats")
write_json(out, file.path(sd, "output.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Protein plots:\n", paste0("  ", plots, collapse = "\n"), "\n")
cat("Review page:", review_stats(db_path), "\n")
