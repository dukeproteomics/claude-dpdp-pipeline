# Step 10b — Stats Graphs
# Static PNG versions of the app's Stats > Graphs tab. Data prep is the app's own create_plot_bg / create_volcano_bg;
# the plot code follows Shiny_Interactive.R (interactive_barplot, _boxplot, _pca2d, _cluster, _heatmap, _stats_volcano)
# with the widget values fixed to sensible defaults. PCA_3D (interactive rgl) is not made.
#
# Usage:
#   Rscript pipeline/steps/10b_stats_graphs.R --db=<db> [--comps=B_v_A] [--include_spqc=TRUE]
#       [--plots=bar,box,pca2d,cluster,heatmap,volcano] [--volcano_highlight=<description text,...>]
#   --comps         comparisons to plot (default: all in stats_comp). Bar/box/PCA/cluster/heatmap use each comparison's
#                   samples, plus an "all samples" set (comparisons + SPQC) when --include_spqc=TRUE.
#
# Output: QC/Stats_<set>_<plot>.png; steps/10_stats/graphs.json

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
use_app(c("Shiny_File.R", "Shiny_Plots.R"))
suppressPackageStartupMessages({ library(ggplot2); library(gplots) })

p <- get_params(db_path)
stats_comp <- read_table("stats_comp", db_path)
comps <- if (!is.null(args$comps)) trimws(strsplit(args$comps, ",")[[1]]) else stats_comp$Name
plot_types <- trimws(strsplit(args$plots %||% "bar,box,pca2d,cluster,heatmap,volcano", ",")[[1]])
include_spqc <- if (!is.null(args$include_spqc)) as_flag(args$include_spqc) else (p$comp_spqc != "0")
out_dir <- p$qc_path
made <- character(0); info <- list()
save_png <- function(name, plot, w = 8, h = 6) {
  f <- gsub("/+", "/", file.path(out_dir, str_c("Stats_", name, ".png")))
  ggsave(f, plot = plot, width = w, height = h, dpi = 150); made <<- c(made, f); f
}

sets <- setNames(lapply(comps, function(x) x), comps)
if (include_spqc) sets[[str_c(paste(comps, collapse = "_"), "_SPQC")]] <- c(comps, p$comp_spqc)

for (set_name in names(sets)) {
  dat <- run_bg(create_plot_bg, list(sets[[set_name]], p$stat_norm, p, db_path, 1), str_c("create_plot_", set_name), db_path)
  df <- dat[[1]]; namex <- dat[[2]]; color_list <- dat[[3]]; groupx <- dat[[4]]
  set_info <- list(samples = namex, proteins_used = nrow(df))

  if ("bar" %in% plot_types) {
    d2 <- data.frame(Sample = factor(namex, levels = namex), Total_Intensity = colSums(df, na.rm = TRUE))
    save_png(str_c(set_name, "_bar"), ggplot(d2, aes(Sample, Total_Intensity)) + geom_bar(stat = "identity", fill = color_list) +
               theme_classic() + ggtitle(str_c("Total intensity - ", set_name)) + xlab(NULL) +
               theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, color = "black")))
  }
  if ("box" %in% plot_types) {
    d <- df; colnames(d) <- namex
    d3 <- tidyr::gather(log2(d), Sample, Intensity, colnames(d)); d3$Sample <- factor(d3$Sample, levels = rev(namex))
    save_png(str_c(set_name, "_box"), ggplot(d3, aes(Sample, Intensity)) +
               geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 1, outlier.size = 1, fill = rev(color_list)) +
               theme_classic() + coord_flip() + ggtitle(str_c("log2 intensity - ", set_name)) + theme(plot.title = element_text(hjust = 0.5)))
  }
  if ("pca2d" %in% plot_types) {
    x_pca <- prcomp(t(df), scale = TRUE)
    ve <- round(x_pca$sdev^2 / sum(x_pca$sdev^2) * 100, 1)
    d_out <- data.frame(x_pca$x[, 1:2], group = factor(groupx), Sample = namex)
    group_cols <- setNames(color_list[!duplicated(groupx)], groupx[!duplicated(groupx)])
    save_png(str_c(set_name, "_pca2d"), ggplot(d_out, aes(PC1, PC2, color = group)) + geom_point(alpha = 0.8, size = 4) +
               geom_text(aes(label = Sample), vjust = -1, size = 3, show.legend = FALSE) +
               scale_color_manual(values = group_cols) + theme_minimal() +
               xlab(str_c("PC1 (", ve[1], "%)")) + ylab(str_c("PC2 (", ve[2], "%)")) +
               ggtitle(str_c("PCA - ", set_name)) + theme(plot.title = element_text(hjust = 0.5)))
    set_info$pca <- list(variance_pct = as.list(setNames(head(ve, 5), str_c("PC", 1:min(5, length(ve))))),
                         scores = data.frame(Sample = namex, group = groupx, PC1 = round(x_pca$x[, 1], 2), PC2 = round(x_pca$x[, 2], 2)))
  }
  if ("cluster" %in% plot_types) {
    d <- df; colnames(d) <- namex; d <- scale(data.frame(t(d)))
    dist_m <- factoextra::get_dist(d, method = "euclidean")
    save_png(str_c(set_name, "_cluster"), factoextra::fviz_dist(dist_m, show_labels = TRUE, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) +
               ggtitle(str_c("Sample distance - ", set_name)) + theme(plot.title = element_text(hjust = 0.5)), w = 7, h = 6)
    dm <- as.matrix(dist_m)
    set_info$nearest_sample <- setNames(as.list(sapply(rownames(dm), function(s) { v <- dm[s, ]; v[s] <- Inf; names(which.min(v)) })), rownames(dm))
    set_info$cor_pearson_log2 <- round(cor(log2(df)), 3); dimnames(set_info$cor_pearson_log2) <- list(namex, namex)
  }
  if ("heatmap" %in% plot_types) {
    d <- df; d$stdev <- apply(d, 1, sd); d <- d[d$stdev > 0.1, seq_len(ncol(df))]
    colnames(d) <- namex; d <- data.matrix(log2(d))
    hr <- hclust(as.dist(1 - cor(t(d), method = "pearson")), method = "complete")
    hc <- hclust(as.dist(1 - cor(d, method = "spearman")), method = "complete")
    mycl <- cutree(hr, h = max(hr$height) / 1.5); rc <- rainbow(length(unique(mycl)), start = 0.1, end = 0.9)[as.vector(mycl)]
    f <- gsub("/+", "/", file.path(out_dir, str_c("Stats_", set_name, "_heatmap.png")))
    png(f, width = 1600, height = 1600, res = 150)
    heatmap.2(d, Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc), col = redgreen(75), labCol = groupx, scale = "row",
              density.info = "none", trace = "none", RowSideColors = rc, main = str_c("Heatmap - ", set_name), margins = c(10, 10))
    dev.off(); made <- c(made, f)
    set_info$column_dendrogram_order <- namex[hc$order]
  }
  info[[set_name]] <- set_info
}

if ("volcano" %in% plot_types) for (cp in comps) {
  hl <- args$volcano_highlight %||% ""
  v <- run_bg(create_volcano_bg, list(cp, p$stat_norm, as.logical(p$checkbox_filter_adjpval), FALSE, 0, 0, hl, TRUE, TRUE, p, db_path, 1),
              str_c("create_volcano_", cp), db_path)
  df <- v[[1]]; xmax <- v[[2]]; ymax <- v[[3]]
  xmax <- max(abs(df$log_fc), na.rm = TRUE)  # app uses max(log_fc); symmetric limits keep all points
  g <- ggplot(df, aes(log_fc, log_pvalue)) + theme_minimal() + geom_point(alpha = 0.4, size = 1, color = "grey40") +
    xlab("log2 fold change") + ylab(if (as.logical(p$checkbox_filter_adjpval)) "-log10 adjusted p" else "-log10 p") +
    ggtitle(str_c("Volcano - ", cp)) + xlim(-xmax, xmax) + ylim(0, ymax) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    geom_vline(xintercept = c(-1, 1) * log2(as.numeric(p$foldchange_cutoff)), linetype = "dotted") +
    geom_hline(yintercept = -log10(as.numeric(p$pvalue_cutoff)), linetype = "dotted") +
    geom_point(data = v[[6]], color = "red", size = 1.5, alpha = 0.8) + geom_point(data = v[[7]], color = "blue", size = 1.5, alpha = 0.8)
  if (nrow(v[[5]])) g <- g + geom_point(data = v[[5]], color = "orange", size = 2.5)
  save_png(str_c(cp, "_volcano"), g, w = 8, h = 7)
  info[[str_c(cp, "_volcano")]] <- list(points = nrow(df), up = nrow(v[[6]]), down = nrow(v[[7]]),
                                         pvalue_axis = if (as.logical(p$checkbox_filter_adjpval)) "adjusted" else "raw",
                                         highlighted = nrow(v[[5]]))
}

out <- list(step = "10b_stats_graphs", db_path = db_path, stat_norm = p$stat_norm, sets = lapply(sets, unlist),
            info = info, plots = made, note = "Static versions of the app's interactive Stats graphs (PCA_3D not made).")
sd <- step_dir(run_dir, "10_stats")
write_json(out, file.path(sd, "graphs.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("\n==== Stats graphs ====\n")
for (s in names(info)) if (!is.null(info[[s]]$pca)) cat(sprintf("  %s PCA: PC1 %s%%, PC2 %s%%\n", s, info[[s]]$pca$variance_pct$PC1, info[[s]]$pca$variance_pct$PC2))
cat("Plots:\n", paste0("  ", made, collapse = "\n"), "\n")
cat("Review page:", review_stats(db_path), "\n")
