# Step 11 - Pathway (pipeline version of the app's Pathway tabs; user, 2026-11-13: STRING and GO, best practice).
# For each comparison, the Up and the Down proteins separately:
#   over-representation through STRING (GO, KEGG, Reactome, WikiPathways, UniProt keywords, STRING clusters, ...) against the
#   proteins tested in that comparison (the measured proteome, not the genome), and the STRING network of the hits with clusters
#   (Louvain on the confidence scores) and what each cluster is enriched for.
# Only accessions / STRING ids go to string-db.org (the version's stable address); no intensities.
#
#   Rscript pipeline/steps/11_pathway.R --db=<db> --action=preview [--from_page] [--settings_from=<db>] [--pathway_<option>=...]
#   Rscript pipeline/steps/11_pathway.R --db=<db> --action=apply --decided_by=user [--from_page]   (accept the results of the last run)
# Options: pipeline/options/11_pathway.R (also the page's Options panel). Output: steps/11_pathway/summary.json, results.json,
# Pathway/Pathway_STRING_<norm>.xlsx, review/11_pathway.html.

args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "..", "common.R"))
source(file.path(here, "..", "review.R"))
source(file.path(here, "..", "contaminants.R"))
source(file.path(here, "..", "string_api.R"))
source(file.path(here, "..", "options.R")); source(file.path(here, "..", "options", "11_pathway.R"))
source(file.path(here, "..", "reviews", "11_pathway.R"))
args <- parse_args()
need(args, c("db", "action"))
if (args$action == "run") args$action <- "preview"
if (!args$action %in% c("preview", "apply")) stop("--action must be preview (run) or apply (accept)")
if (args$action == "apply" && !identical(args$decided_by, "user")) stop("Accepting needs --decided_by=user (after the analyst has reviewed the page).")
db_path <- normalizePath(args$db, mustWork = TRUE)
run_dir <- dirname(db_path)
use_app("Shiny_File.R")
p <- get_params(db_path)
if (!table_exists("stats_comp", db_path)) stop("Run Stats (Step 10) first.")
sc <- read_table("stats_comp", db_path)
pep_out <- is_peptide_output(p)
sd <- step_dir(run_dir, "11_pathway")

# organisms in the tested data (most first)
base_tab <- read_table(base_table(p), db_path)
orgs <- sort(table(organism_of(base_tab$Description)[nzchar(organism_of(base_tab$Description))]), decreasing = TRUE)
orgs <- setNames(as.integer(orgs), names(orgs))
opts <- pathway_options(p, orgs)
spec <- opts_spec(opts)
args <- page_args(args, run_dir, "11_pathway", opts, args$action)
settings <- resolve_settings(spec, args)
if (args$action == "apply") settings <- apply_settings(run_dir, "11_pathway", spec, args, settings, opts)
opts_validate(opts, settings)

if (args$action == "apply") {
  log_decision(run_dir, "11_pathway", c(list(accepted = TRUE), settings), note = args$note %||% "", analyst = args$analyst)
  cat("Pathway results accepted\nReview page:", review_pathway(db_path), "\n"); quit(save = "no")
}
for (k in names(settings)) p[[k]] <- settings[[k]]
invisible(write_params(p, db_path))

taxon <- as.integer(settings$pathway_species)
org_name <- if (length(orgs) && identical(unname(STRING_TAXA[names(orgs)[1]]), taxon)) names(orgs)[1] else names(STRING_TAXA)[match(taxon, STRING_TAXA)][1]
cats <- pathway_categories(settings)
sp <- qc_spike_settings(run_dir)
v <- string_version()
cat(sprintf("STRING %s (%s), species %d%s\n", v$version, v$base, taxon, if (!is.na(org_name)) paste0(" (", org_name, ")") else ""))

# Per comparison: the tested proteins (background) and the hits, one row per protein (Peptide output: collapsed by accession)
comp_sets <- lapply(seq_len(nrow(sc)), function(i) {
  nm <- sc$Name[i]; st <- read_table(final_table(sc, i, p), db_path)
  acc <- sub(";.*", "", st$Accession)
  keep <- rep(TRUE, nrow(st))
  if (as_flag(settings$pathway_species_only)) {   # by taxon id (two names can share one, e.g. Cryptococcus neoformans / deneoformans)
    o <- organism_of(st$Description); ot <- unname(STRING_TAXA[o]); keep <- keep & (!nzchar(o) | (!is.na(ot) & ot == taxon) | (is.na(ot) & o == names(orgs)[1] & taxon == (STRING_TAXA[names(orgs)[1]] %||% -1))) }
  if (as_flag(settings$pathway_exclude_flagged)) keep <- keep & flag_contaminants(st$Accession, st$Description, sp$spike, sp$constant) == ""
  fc <- comp_get(st, nm, "_FC"); pv <- comp_get(st, nm, if (identical(settings$pathway_p_type, "raw")) "_pval" else "_adjpval") %||% comp_get(st, nm, "_pval")
  if (identical(settings$pathway_hits, "custom")) {
    call <- ifelse(!is.na(pv) & pv <= settings$pathway_p_cutoff & fc >= settings$pathway_fc_cutoff, "Up",
                   ifelse(!is.na(pv) & pv <= settings$pathway_p_cutoff & fc <= -settings$pathway_fc_cutoff, "Down", ""))
  } else call <- ifelse(is.na(st$Stats), "", st$Stats)
  praw <- comp_get(st, nm, "_pval")
  d <- data.frame(acc = acc, gene = ifelse(is.na(st$Genes), "", st$Genes), fc = fc, p = praw, call = call, stringsAsFactors = FALSE)[keep, ]
  # one row per protein: the most significant row (sites in Peptide output); Up / Down if any row is
  d <- d[order(d$p), ]
  up <- unique(d$acc[d$call == "Up"]); down <- unique(d$acc[d$call == "Down"])
  d <- d[!duplicated(d$acc), ]
  list(name = nm, N = sc$FactorsN[i], D = sc$FactorsD[i], prot = d, up = up, down = down, left_out = sum(!keep))
})

# STRING ids for every tested protein (cached per species)
all_acc <- unique(unlist(lapply(comp_sets, function(x) x$prot$acc)))
cache <- file.path(sd, sprintf("string_ids_%d_%s.csv", taxon, gsub("[^0-9.]", "", v$version)))
m <- if (file.exists(cache)) read.csv(cache, stringsAsFactors = FALSE) else data.frame(query = character(0), string_id = character(0), name = character(0))
todo <- setdiff(all_acc, m$query)
if (length(todo)) { cat(sprintf("Mapping %d accessions to STRING...\n", length(todo))); m <- rbind(m, string_map(v$base, todo, taxon)); write.csv(m, cache, row.names = FALSE) }
sid <- setNames(m$string_id, m$query); sname <- setNames(m$name, m$string_id)
# names to show: STRING's preferred name, or the gene from the data when STRING only has a UniProt code (e.g. Q5KPY0_CRYNJ)
gene_of <- unlist(lapply(comp_sets, function(x) setNames(x$prot$gene, x$prot$acc))); gene_of <- gene_of[!duplicated(names(gene_of))]
code <- grepl("^[A-Z0-9]+_[A-Z0-9]+$", sname); g_by_id <- setNames(gene_of[m$query], m$string_id)
sname[code & !is.na(g_by_id[names(sname)]) & nzchar(g_by_id[names(sname)])] <- g_by_id[names(sname)][code & !is.na(g_by_id[names(sname)]) & nzchar(g_by_id[names(sname)])]
cat(sprintf("STRING ids: %d of %d proteins\n", sum(all_acc %in% m$query), length(all_acc)))

enrich_table <- function(ids, bg) {
  if (length(ids) < 3) return(list(terms = data.frame(), note = sprintf("%d mapped hit%s: too few to test (3 or more needed).", length(ids), if (length(ids) == 1) "" else "s")))
  e <- string_enrich(v$base, ids, bg, taxon)
  if (!NROW(e)) return(list(terms = data.frame(), note = "No term is enriched (FDR < 0.05 in STRING)."))
  e <- e[e$category %in% cats & e$fdr <= settings$pathway_fdr & e$number_of_genes >= settings$pathway_min_genes, ]
  if (!nrow(e)) return(list(terms = data.frame(), note = "No term passes the FDR and minimum-hits cutoffs."))
  n <- length(ids); N <- length(bg)
  t <- data.frame(category = e$category, source = unname(PATHWAY_CATEGORY_NAMES[e$category]), term = e$term, description = e$description,
                  hits = e$number_of_genes, term_size = e$number_of_genes_in_background,
                  fold = round((e$number_of_genes / n) / (e$number_of_genes_in_background / N), 2), p = signif(e$p_value, 3), fdr = signif(e$fdr, 3),
                  stringsAsFactors = FALSE)
  t$ids <- lapply(e$inputGenes, function(x) unlist(strsplit(paste(unlist(x), collapse = ","), ",")))
  t <- t[order(t$fdr, -t$fold), ]
  list(terms = t, note = sprintf("%d terms (FDR <= %s, >= %d hits).", nrow(t), settings$pathway_fdr, settings$pathway_min_genes))
}

results <- lapply(comp_sets, function(cs) {
  cat(sprintf("\n== %s: %d tested, %d Up, %d Down\n", cs$name, nrow(cs$prot), length(cs$up), length(cs$down)))
  bg <- unique(na.omit(sid[cs$prot$acc])); up <- unique(na.omit(sid[cs$up])); dn <- unique(na.omit(sid[cs$down]))
  eu <- enrich_table(up, bg); ed <- enrich_table(dn, bg)
  cat(sprintf("  Up: %s  Down: %s\n", eu$note, ed$note))
  # network of the hits (Up and Down; the smallest p first when there are more than network_max)
  hits <- cs$prot[cs$prot$acc %in% c(cs$up, cs$down) & cs$prot$acc %in% names(sid), ]
  hits <- head(hits[order(hits$p), ], settings$pathway_network_max)
  net <- NULL
  if (nrow(hits) >= 2) {
    ids <- unname(sid[hits$acc])
    ed_ <- string_network(v$base, ids, taxon, settings$pathway_score, settings$pathway_network_type)
    ed_ <- ed_[ed_$a %in% ids & ed_$b %in% ids & ed_$a != ed_$b, ]
    ed_ <- ed_[!duplicated(t(apply(ed_[, c("a", "b")], 1, sort))), , drop = FALSE]
    g <- igraph::graph_from_data_frame(ed_, directed = FALSE, vertices = data.frame(name = ids))
    set.seed(1)
    cl <- if (igraph::ecount(g)) igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$score)) else setNames(seq_along(ids), ids)
    sizes <- table(cl); cl_rank <- match(as.character(cl), names(sort(sizes, decreasing = TRUE)))   # cluster 1 = largest
    deg <- igraph::degree(g)
    # layout: the connected proteins (Fruchterman-Reingold); proteins without interactions in a grid below them
    xy <- matrix(0, length(ids), 2); con <- deg > 0
    if (sum(con) >= 2) { set.seed(1); g2 <- igraph::induced_subgraph(g, which(con)); xy[con, ] <- igraph::layout_with_fr(g2, weights = igraph::E(g2)$score) }
    if (any(!con)) { k <- which(!con); rng <- if (any(con)) range(xy[con, 1]) else c(0, 10); per <- max(1, ceiling(sqrt(length(k) * 3)))
      y0 <- if (any(con)) min(xy[con, 2]) - 0.08 * diff(range(xy[con, 2]) + c(0, 1)) - 1 else 0
      xy[k, 1] <- rng[1] + ((seq_along(k) - 1) %% per) * diff(rng + c(0, 1e-9)) / max(1, per - 1)
      xy[k, 2] <- y0 - ((seq_along(k) - 1) %/% per) * max(0.5, diff(rng) / per) }
    # what each cluster (5+ proteins) is enriched for: the top terms (at most 8 clusters)
    big <- sort(unique(cl_rank[cl_rank %in% which(as.integer(sort(sizes, decreasing = TRUE)) >= 5)]))
    clusters <- lapply(head(big, 8), function(k) { mem <- ids[cl_rank == k]
      e <- tryCatch(string_enrich(v$base, mem, bg, taxon), error = function(e) data.frame())
      top <- if (NROW(e)) { e <- e[e$category %in% cats & e$fdr <= settings$pathway_fdr, ]; head(e[order(e$fdr), ], 3) } else data.frame()
      list(cluster = k, size = length(mem), terms = if (NROW(top)) paste0(top$description, " (", unname(PATHWAY_CATEGORY_NAMES[top$category]), ", FDR ", signif(top$fdr, 2), ")") else character(0)) })
    nodes <- data.frame(id = ids, name = unname(sname[ids]), acc = hits$acc, gene = hits$gene, fc = signif(hits$fc, 4), p = signif(hits$p, 3), call = hits$call,
                        cluster = cl_rank, degree = as.integer(deg[ids]), x = round(xy[, 1], 3), y = round(xy[, 2], 3), stringsAsFactors = FALSE)
    nodes$cluster[nodes$degree == 0] <- 0L   # no interactions: no cluster
    net <- list(nodes = nodes, edges = ed_, clusters = clusters, link = string_link(v$base, ids, taxon, settings$pathway_score, settings$pathway_network_type),
                shown = nrow(hits), of = length(unique(c(cs$up, cs$down))))
    cat(sprintf("  Network: %d proteins, %d interactions (score >= %s), %d clusters of 5+\n", nrow(nodes), nrow(ed_), settings$pathway_score, length(big)))
  }
  list(name = cs$name, N = cs$N, D = cs$D, tested = nrow(cs$prot), mapped_background = length(bg), left_out = cs$left_out,
       up = list(n = length(cs$up), mapped = length(up), terms = eu$terms, note = eu$note),
       down = list(n = length(cs$down), mapped = length(dn), terms = ed$terms, note = ed$note),
       network = net, names = as.list(sname[unique(c(up, dn))]),
       hits = { h <- cs$prot[cs$prot$acc %in% c(cs$up, cs$down) & cs$prot$acc %in% names(sid), ]   # every mapped hit: for the term protein lists
                setNames(lapply(seq_len(nrow(h)), function(j) list(acc = h$acc[j], name = unname(sname[sid[h$acc[j]]]), fc = signif(h$fc[j], 4), p = signif(h$p[j], 3),
                  call = if (h$acc[j] %in% cs$up && h$acc[j] %in% cs$down) "Up and Down (sites)" else h$call[j])), unname(sid[h$acc])) })
})

# --- Ranked analysis (GSEA, fgsea): every tested protein ranked by sign(FC) x -log10(p), no cutoff; STRING's annotation of the
# species (same terms as the over-representation, except KEGG); gene sets restricted to the tested proteins, min/max size.
gsea_note <- NULL
if (as_flag(settings$pathway_gsea)) {
  gcats <- setdiff(cats, "KEGG")
  tt <- tryCatch(string_terms(taxon, file.path(dirname(run_dir), ".dpdp_cache"), gcats), error = function(e) { gsea_note <<- conditionMessage(e); NULL })
  if (!is.null(tt)) {
    tdesc <- tt[!duplicated(tt$term), c("term", "category", "description")]
    sets_all <- split(tt$id, tt$term)
    for (i in seq_along(results)) {
      cs <- comp_sets[[i]]; d <- cs$prot[cs$prot$acc %in% names(sid) & !is.na(cs$prot$p) & !is.na(cs$prot$fc), ]
      ids <- unname(sid[d$acc]); keep <- !duplicated(ids); d <- d[keep, ]; ids <- ids[keep]
      l2 <- sign(d$fc) * log2(pmax(abs(d$fc), 1))
      stat <- setNames(sign(d$fc) * -log10(pmax(d$p, 1e-300)) + l2 * 1e-6, ids)   # ties broken by the fold change
      sets <- lapply(sets_all, function(x) intersect(x, ids))
      sets <- sets[lengths(sets) >= settings$pathway_gsea_min & lengths(sets) <= settings$pathway_gsea_max]
      set.seed(1)
      g <- suppressWarnings(fgsea::fgseaMultilevel(sets, stat, minSize = settings$pathway_gsea_min, maxSize = settings$pathway_gsea_max, eps = 0, nproc = 1))
      g <- g[!is.na(g$padj), ]; g <- g[order(g$padj, -abs(g$NES)), ]
      m <- match(g$pathway, tdesc$term)
      tab <- data.frame(category = tdesc$category[m], source = unname(PATHWAY_CATEGORY_NAMES[tdesc$category[m]]), term = g$pathway, description = tdesc$description[m],
                        size = g$size, NES = round(g$NES, 3), p = signif(g$pval, 3), fdr = signif(g$padj, 3), stringsAsFactors = FALSE)
      tab$leading <- lapply(g$leadingEdge, function(x) unname(x))
      sig <- tab[tab$fdr <= settings$pathway_fdr, ]
      # the ranked list and the members' ranks of the top terms, for the running-score plot on the page
      ord <- order(stat, decreasing = TRUE); ranked <- names(stat)[ord]
      top_terms <- head(tab$term, 80)
      results[[i]]$gsea <- list(n = length(stat), sets_tested = nrow(tab), terms = sig, all_top = head(tab, 300),
                                ranked_stat = round(unname(stat[ord]), 3),
                                members = setNames(lapply(top_terms, function(t) match(sets[[t]], ranked) - 1L), top_terms),
                                note = sprintf("%d of %d gene sets with FDR <= %s (%d up, %d down).", nrow(sig), nrow(tab), settings$pathway_fdr, sum(sig$NES > 0), sum(sig$NES < 0)))
      for (x in unlist(tab$leading)) if (is.null(results[[i]]$names[[x]])) results[[i]]$names[[x]] <- unname(sname[x])
      cat(sprintf("  GSEA %s: %s\n", cs$name, results[[i]]$gsea$note))
    }
  }
  if (!is.null(gsea_note)) cat("GSEA not run:", gsea_note, "\n")
}

# --- GO with the hierarchy (topGO weight01, Alexa et al. 2006: a parent term's signal coming from a child term is discounted),
# Up and Down separately against the tested proteins; then terms that mean nearly the same thing are grouped (GOSemSim, Wang
# graph similarity, which needs no organism database): one representative per group (the most significant).
go_note <- NULL
if (as_flag(settings$pathway_topgo) && as_flag(settings$pathway_go)) {
  suppressPackageStartupMessages({ library(topGO); library(GOSemSim) })
  ga <- tryCatch(string_terms(taxon, file.path(dirname(run_dir), ".dpdp_cache"), c("Process", "Function", "Component")), error = function(e) { go_note <<- conditionMessage(e); NULL })
  if (!is.null(ga)) {
    onts <- c(BP = "Process", MF = "Function", CC = "Component")
    semd <- list()
    for (i in seq_along(results)) {
      cs <- comp_sets[[i]]; bg <- unique(na.omit(sid[cs$prot$acc]))
      hits <- list(up = unique(na.omit(sid[cs$up])), down = unique(na.omit(sid[cs$down])))
      res_go <- list()
      for (o in names(onts)) {
        an <- ga[ga$category == onts[[o]] & ga$id %in% bg, ]
        g2go <- split(an$term, an$id); uni <- names(g2go)
        res_go[[o]] <- lapply(hits, function(h) {
          h <- intersect(h, uni)
          if (length(h) < 3) return(list(terms = data.frame(), note = sprintf("%d annotated hits: too few to test.", length(h))))
          allg <- factor(as.integer(uni %in% h), levels = c(0, 1)); names(allg) <- uni
          gd <- suppressMessages(new("topGOdata", ontology = o, allGenes = allg, annot = annFUN.gene2GO, gene2GO = g2go, nodeSize = 5))
          rw <- suppressMessages(runTest(gd, algorithm = "weight01", statistic = "fisher")); rc <- suppressMessages(runTest(gd, algorithm = "classic", statistic = "fisher"))
          pw <- score(rw); pc <- score(rc)[names(pw)]
          st <- termStat(gd, names(pw))
          k <- names(pw)[pw <= settings$pathway_topgo_p & st[names(pw), "Significant"] >= settings$pathway_min_genes]
          if (!length(k)) return(list(terms = data.frame(), note = sprintf("No GO term with weight01 p <= %s.", settings$pathway_topgo_p)))
          t <- data.frame(term = k, description = unname(vapply(k, function(x) tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[x]]), error = function(e) x), "")),
                          annotated = st[k, "Annotated"], hits = st[k, "Significant"], expected = round(st[k, "Expected"], 2),
                          fold = round(st[k, "Significant"] / pmax(st[k, "Expected"], 1e-9), 2), p = signif(pw[k], 3), p_classic = signif(pc[k], 3), stringsAsFactors = FALSE)
          t <- t[order(t$p, -t$fold), ]
          mem <- genesInTerm(gd, t$term); t$ids <- lapply(t$term, function(x) intersect(mem[[x]], h))
          # group similar terms
          t$group <- seq_len(nrow(t)); t$representative <- TRUE
          if (nrow(t) >= 2) {
            if (is.null(semd[[o]])) semd[[o]] <<- suppressMessages(godata(ont = o, computeIC = FALSE))
            m <- suppressWarnings(mgoSim(t$term, t$term, semData = semd[[o]], measure = "Wang", combine = NULL))
            m[is.na(m)] <- 0; diag(m) <- 1
            cl <- cutree(hclust(as.dist(1 - m[t$term, t$term]), method = "average"), h = 1 - settings$pathway_topgo_sim)
            rep_of <- tapply(seq_len(nrow(t)), cl, min)          # rows are sorted by p: the first of a group is its most significant term
            t$group <- match(cl, unique(cl))                     # groups numbered by their most significant term
            t$representative <- seq_len(nrow(t)) %in% rep_of
          }
          list(terms = t, note = sprintf("%d GO terms with weight01 p <= %s in %d groups (%d annotated hits of %d tested proteins).",
                                         nrow(t), settings$pathway_topgo_p, length(unique(t$group)), length(h), length(uni)))
        })
      }
      results[[i]]$go <- res_go
      for (o in names(res_go)) for (d in c("up", "down")) for (x in unlist(res_go[[o]][[d]]$terms$ids)) if (is.null(results[[i]]$names[[x]])) results[[i]]$names[[x]] <- unname(sname[x])
      cat(sprintf("  GO %s: %s\n", cs$name, paste(vapply(names(res_go), function(o) sprintf("%s up %d / down %d", o, NROW(res_go[[o]]$up$terms), NROW(res_go[[o]]$down$terms)), ""), collapse = ", ")))
    }
  }
  if (!is.null(go_note)) cat("GO (topGO) not run:", go_note, "\n")
}

# Excel: the term tables (Up, Down) and the network proteins of each comparison
xdir <- file.path(run_dir, "Pathway"); dir.create(xdir, showWarnings = FALSE)
xf <- file.path(xdir, sprintf("Pathway_STRING_%s.xlsx", trimws(p$stat_norm)))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "About")
openxlsx::writeData(wb, "About", data.frame(Item = c("Run", "STRING version", "STRING address", "Species (taxon)", "Background", "Hits", "Categories", "FDR cutoff", "Minimum hits", "Network confidence", "Network type", "Made"),
  Value = c(basename(run_dir), v$version, v$base, sprintf("%d %s", taxon, org_name %||% ""), "proteins tested in each comparison",
            if (identical(settings$pathway_hits, "custom")) sprintf("custom: %s p <= %s, |FC| >= %s", settings$pathway_p_type, settings$pathway_p_cutoff, settings$pathway_fc_cutoff) else "Stats calls",
            paste(PATHWAY_CATEGORY_NAMES[cats], collapse = ", "), settings$pathway_fdr, settings$pathway_min_genes, settings$pathway_score, settings$pathway_network_type,
            format(Sys.time(), "%Y-%m-%d %H:%M"))))
for (i in seq_along(results)) {
  r <- results[[i]]
  for (dir in c("up", "down")) {
    t <- r[[dir]]$terms; sh <- sprintf("C%d %s", i, if (dir == "up") "Up" else "Down")
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, data.frame(Comparison = r$name, Direction = dir, Note = r[[dir]]$note), startRow = 1)
    if (NROW(t)) { t$proteins <- vapply(t$ids, function(x) paste(unname(sname[x]), collapse = ", "), ""); t$ids <- NULL
      openxlsx::writeDataTable(wb, sh, t, startRow = 3, tableStyle = "TableStyleLight9") }
  }
  if (!is.null(r$gsea) && NROW(r$gsea$all_top)) { sh <- sprintf("C%d GSEA", i); openxlsx::addWorksheet(wb, sh); t <- r$gsea$all_top
    t$leading_edge <- vapply(t$leading, function(x) paste(unname(sname[x]), collapse = ", "), ""); t$leading <- NULL
    openxlsx::writeData(wb, sh, data.frame(Comparison = r$name, Note = r$gsea$note), startRow = 1); openxlsx::writeDataTable(wb, sh, t, startRow = 3, tableStyle = "TableStyleLight9") }
  if (!is.null(r$go)) { gt <- do.call(rbind, lapply(names(r$go), function(o) do.call(rbind, lapply(c("up", "down"), function(d) { t <- r$go[[o]][[d]]$terms
      if (!NROW(t)) return(NULL); t$proteins <- vapply(t$ids, function(x) paste(unname(sname[x]), collapse = ", "), ""); t$ids <- NULL
      cbind(ontology = o, direction = d, t) }))))
    if (NROW(gt)) { sh <- sprintf("C%d GO topGO", i); openxlsx::addWorksheet(wb, sh); openxlsx::writeDataTable(wb, sh, gt, tableStyle = "TableStyleLight9") } }
  if (!is.null(r$network)) { sh <- sprintf("C%d network", i); openxlsx::addWorksheet(wb, sh); openxlsx::writeDataTable(wb, sh, r$network$nodes[, c("name", "acc", "gene", "fc", "p", "call", "cluster", "degree")], tableStyle = "TableStyleLight9") }
}
openxlsx::saveWorkbook(wb, xf, overwrite = TRUE)

summary <- list(step = "11_pathway", action = "preview", db_path = db_path, settings = settings, string = v, species = taxon, organism = org_name,
                organisms_in_data = as.list(head(orgs, 6)), categories = cats, excel = xf,
                comparisons = lapply(results, function(r) list(name = r$name, tested = r$tested, mapped_background = r$mapped_background,
                  up = r$up$n, up_terms = NROW(r$up$terms), down = r$down$n, down_terms = NROW(r$down$terms),
                  top_up = if (NROW(r$up$terms)) head(paste0(r$up$terms$description, " (", r$up$terms$source, ", FDR ", r$up$terms$fdr, ")"), 5) else list(),
                  top_down = if (NROW(r$down$terms)) head(paste0(r$down$terms$description, " (", r$down$terms$source, ", FDR ", r$down$terms$fdr, ")"), 5) else list(),
                  network_proteins = NROW(r$network$nodes), network_edges = NROW(r$network$edges),
                  gsea_terms = NROW(r$gsea$terms), top_gsea = if (NROW(r$gsea$terms)) head(paste0(r$gsea$terms$description, " (NES ", r$gsea$terms$NES, ", FDR ", r$gsea$terms$fdr, ")"), 5) else list())),
                gsea_note = gsea_note, go_note = go_note)
write_summary(summary, sd)
write_json(results, file.path(sd, "results.json"), auto_unbox = TRUE, digits = NA, na = "null", null = "null")
cat("\nExcel:", xf, "\nReview page:", review_pathway(db_path), "\n>>> Waiting for the analyst's review.\n")
