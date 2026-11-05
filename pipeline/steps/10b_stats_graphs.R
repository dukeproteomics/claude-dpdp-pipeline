# Step 10b — Stats Graphs
# Data for the interactive versions of the app's Stats > Graphs tab, drawn on the Stats review page: total intensity and
# log2 distribution per sample, PCA, sample distance, and a heatmap. Volcanoes are drawn on the page from the stats tables.
# Same data as the app's create_plot_bg (the whole-dataset table of the tested normalization, each comparison's samples), with two
# changes: PTM runs use the PTM peptides only (as Stats), and the heatmap uses the 1,000 most variable rows (the app clusters every
# row, which needs an n x n distance matrix: ~118k peptides ran out of memory).
#
# Usage:
#   Rscript pipeline/steps/10b_stats_graphs.R --db=<db> [--comps=B_v_A] [--include_spqc=TRUE] [--heatmap_rows=1000] [--png=FALSE]
#   --comps         comparisons (default: all in stats_comp). Each comparison's samples is one set, plus an "all samples" set
#                   (comparisons + SPQC) when --include_spqc=TRUE.
#   --png=TRUE      also write the app-style static PNGs (bar, box, PCA, sample distance, heatmap) to QC/Stats_<set>_<plot>.png
#
# Output: steps/10_stats/graphs.json (chart data); QC/Stats_*.png only with --png=TRUE

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "contaminants.R"))
source(file.path(here, "..", "reviews", "10_stats.R"))

args <- parse_args()
need(args, "db")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app(c("Shiny_File.R"))

p <- get_params(db_path)
stats_comp <- read_table("stats_comp", db_path)
design <- read_table("design", db_path)
group_order <- read_table("sample_groups", db_path)$Group
comps <- if (!is.null(args$comps)) trimws(strsplit(args$comps, ",")[[1]]) else stats_comp$Name
include_spqc <- if (!is.null(args$include_spqc)) as_flag(args$include_spqc) else (p$comp_spqc != "0")
n_heat <- as.integer(args$heatmap_rows %||% 1000)
make_png <- as_flag(args$png %||% FALSE)
n_s <- as.integer(p$sample_number)
locs <- function(x) if (is.na(x) || !nzchar(x)) integer(0) else as.integer(trimws(strsplit(as.character(x), ",")[[1]]))

# whole-dataset table of the tested normalization (as create_plot_bg); PTM runs: PTM rows only
tab <- read_table(base_table(p), db_path)
ptm_only <- ptm_report(p)
if (ptm_only) tab <- tab[grepl(p$peptide_report_grep, tab$Sequence), ]
X <- as.matrix(tab[, (ncol(tab) - n_s + 1):ncol(tab)]); storage.mode(X) <- "numeric"
X[!is.na(X) & X <= 0] <- NA
row_label <- if (is_peptide_output(p)) {
  site <- if (!is.null(tab$Protein_PTM_Loc)) tab$Protein_PTM_Loc else ""
  trimws(paste(ifelse(is.na(tab$Genes) | tab$Genes == "", tab$Accession, tab$Genes), ifelse(nzchar(site), site, gsub("_", "", tab$Sequence))))
} else ifelse(is.na(tab$Genes) | tab$Genes == "", tab$Accession, paste(tab$Genes, tab$Accession))

sets <- lapply(comps, function(cp) { i <- which(stats_comp$Name == cp); c(locs(stats_comp$N_loc[i]), locs(stats_comp$D_loc[i])) })
names(sets) <- comps
if (include_spqc) {
  sp <- unique(unlist(lapply(comps, function(cp) locs(stats_comp$SPQC_loc[stats_comp$Name == cp]))))
  sets[["all samples"]] <- sort(unique(c(unlist(sets), sp)))
}
grp_of <- rep(read_table("sample_groups", db_path)$Group, as.integer(read_table("sample_groups", db_path)$Count))

info <- list(); made <- character(0)
for (set_name in names(sets)) {
  pos <- sets[[set_name]]
  m <- X[, pos, drop = FALSE]; m <- m[rowSums(is.na(m)) == 0, , drop = FALSE]   # complete rows, as the app's plots (imputed data)
  lab <- design$Label[pos]; grp <- grp_of[pos]; lm <- log2(m)
  q <- t(apply(lm, 2, function(v) as.numeric(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95)))))
  sd_r <- apply(lm, 1, sd); lv <- lm[sd_r > 0, , drop = FALSE]   # constant rows can't be scaled (PCA, distance)
  pc <- prcomp(t(lv), scale. = TRUE); ve <- round(pc$sdev^2 / sum(pc$sdev^2) * 100, 1)
  dm <- as.matrix(dist(scale(t(lv))))   # app: factoextra::get_dist(scale(t(data)), "euclidean")
  keep <- order(sd_r, decreasing = TRUE)[seq_len(min(n_heat, sum(sd_r > 0)))]
  h <- lm[keep, , drop = FALSE]; z <- t(scale(t(h)))
  ro <- hclust(as.dist(1 - cor(t(h), method = "pearson")), method = "complete")$order
  co <- hclust(as.dist(1 - cor(h, method = "spearman")), method = "complete")$order
  info[[set_name]] <- list(
    samples = lab, groups = grp, rows_used = nrow(m), ptm_only = ptm_only,
    total = round(colSums(m), 0), q = unname(split(round(q, 3), row(q))),
    pca = list(variance_pct = as.list(setNames(head(ve, 5), paste0("PC", seq_len(min(5, length(ve)))))),
               scores = data.frame(Sample = lab, group = grp, PC1 = round(pc$x[, 1], 2), PC2 = round(pc$x[, 2], 2),
                                   PC3 = if (ncol(pc$x) >= 3) round(pc$x[, 3], 2) else 0)),
    distance = list(samples = lab, d = unname(round(dm, 2))),
    nearest_sample = setNames(as.list(vapply(seq_len(nrow(dm)), function(k) { v <- dm[k, ]; v[k] <- Inf; lab[which.min(v)] }, "")), lab),
    heatmap = list(rows = nrow(h), of = nrow(m), labels = row_label[rowSums(is.na(X[, pos, drop = FALSE])) == 0][keep][ro],
                   samples = lab[co], groups = grp[co], z = unname(round(z[ro, co, drop = FALSE], 2))))

  if (make_png) {
    suppressPackageStartupMessages({ library(ggplot2) })
    cols <- design$colorlist[pos]
    sv <- function(nm, g, w = 8, h = 6) { f <- gsub("/+", "/", file.path(p$qc_path, str_c("Stats_", gsub(" ", "_", set_name), "_", nm, ".png"))); ggsave(f, g, width = w, height = h, dpi = 150); made <<- c(made, f) }
    sv("bar", ggplot(data.frame(Sample = factor(lab, lab), T = colSums(m)), aes(Sample, T)) + geom_col(fill = cols) + theme_classic() + xlab(NULL) + ylab("total intensity") + theme(axis.text.x = element_text(angle = 90)))
    d3 <- data.frame(Sample = factor(rep(lab, each = nrow(lm)), rev(lab)), I = as.vector(lm))
    sv("box", ggplot(d3, aes(Sample, I)) + geom_boxplot(fill = rev(cols)[match(levels(d3$Sample), rev(lab))]) + coord_flip() + theme_classic() + ylab("log2 intensity"))
    sv("pca2d", ggplot(info[[set_name]]$pca$scores, aes(PC1, PC2, color = group)) + geom_point(size = 4) + geom_text(aes(label = Sample), vjust = -1, size = 3) + theme_minimal())
  }
}

out <- list(step = "10b_stats_graphs", db_path = db_path, stat_norm = p$stat_norm, group_order = group_order,
            sets = lapply(sets, function(s) design$Label[s]), info = info, plots = made,
            note = sprintf("Chart data for the Stats page (interactive).%s Heatmap: the %d most variable rows, z-scored per row.",
                           if (ptm_only) sprintf(" PTM peptides only (%s).", p$peptide_report_grep) else "", n_heat))
sd <- step_dir(run_dir, "10_stats")
write_json(out, file.path(sd, "graphs.json"), auto_unbox = TRUE, pretty = FALSE, digits = NA)
cat("\n==== Stats graphs (data) ====\n")
for (s in names(info)) cat(sprintf("  %s: %d samples, %d rows; PCA PC1 %s%%, PC2 %s%%; heatmap %d rows\n", s, length(info[[s]]$samples),
                                   info[[s]]$rows_used, info[[s]]$pca$variance_pct$PC1, info[[s]]$pca$variance_pct$PC2, info[[s]]$heatmap$rows))
if (length(made)) cat("PNGs:\n", paste0("  ", made, collapse = "\n"), "\n")
cat("Review page:", review_stats(db_path), "\n")
