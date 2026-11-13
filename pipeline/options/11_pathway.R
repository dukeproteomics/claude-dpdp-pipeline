# Step 11 - Pathway (pipeline; the app's Pathway tabs, improved): over-representation through STRING for the Up and the Down
# proteins of each comparison separately, against the proteins tested in that comparison, and the STRING network of the hits.
# p: run params; orgs: named counts of the organisms in the tested data (organism_of on the descriptions).
pathway_options <- function(p = NULL, orgs = NULL) {
  main <- if (length(orgs)) names(orgs)[1] else ""
  taxon <- if (nzchar(main) && main %in% names(STRING_TAXA)) unname(STRING_TAXA[main]) else 9606
  mixed <- length(orgs) > 1 && orgs[1] / sum(orgs) < 0.98
  org_info <- if (length(orgs)) paste0("In the tested data: ", paste(sprintf("%s %s", head(names(orgs), 4), format(head(orgs, 4), big.mark = ",")), collapse = ", "), ".")
  list(
    opt_section("Proteins", list(
      opt("pathway_species", "Organism (NCBI taxon id)", "number", taxon,
          paste0("STRING species (", if (nzchar(main)) main else "not detected", "). Common: human 9606, mouse 10090, rat 10116, yeast 4932, E. coli 511145, Cryptococcus JEC21 214684."),
          min = 1, step = 1, app = FALSE, info = org_info),
      opt("pathway_species_only", "Only proteins of this organism", "checkbox", TRUE,
          "Leaves out proteins of other species (spike-ins, mixed-species designs, contaminants) from the hits and the background.", app = FALSE),
      opt("pathway_exclude_flagged", "Leave out flagged proteins", "checkbox", TRUE,
          "Digestion trypsin, human keratins (non-human samples) and the QC spike-ins (Stats flags).", app = FALSE))),
    opt_section("Hits", note = "Up and Down are tested separately (pooling them hides or invents signal). Phospho runs: a protein is a hit if any of its tested sites is.", list(
      opt("pathway_hits", "Which proteins are hits", "select", "stats calls",
          "stats calls: the Up / Down calls of the Stats step. custom: your own p and fold change cutoffs.", choices = c("stats calls", "custom"), app = FALSE),
      opt("pathway_p_type", "p value", "select", "adjusted", "", choices = c("adjusted", "raw"), depends = "pathway_hits=custom", app = FALSE),
      opt("pathway_p_cutoff", "p cutoff", "number", 0.05, "", depends = "pathway_hits=custom", min = 0, max = 1, step = "any", app = FALSE),
      opt("pathway_fc_cutoff", "Fold change cutoff", "number", 1.5, "", depends = "pathway_hits=custom", min = 1, step = "any", app = FALSE))),
    opt_section("Enrichment (STRING)", note = "Background: every protein tested in the comparison (what the experiment could see), not the genome.", list(
      opt("pathway_go", "GO (biological process, molecular function, cellular component)", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_kegg", "KEGG pathways", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_reactome", "Reactome pathways", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_wiki", "WikiPathways", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_keywords", "UniProt keywords", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_string_clusters", "STRING clusters (local network clusters)", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_domains", "Protein domains (Pfam, InterPro, SMART)", "checkbox", FALSE, "", app = FALSE),
      opt("pathway_other", "Compartments, tissues, diseases, publications", "checkbox", FALSE, "", app = FALSE),
      opt("pathway_fdr", "FDR cutoff", "number", 0.05, "Terms with a higher FDR are not shown.", min = 0, max = 1, step = "any", app = FALSE),
      opt("pathway_min_genes", "Minimum hits in a term", "number", 3, "Terms with fewer hits are not shown.", min = 1, step = 1, app = FALSE))),
    opt_section("Network (STRING)", list(
      opt("pathway_score", "Interaction confidence", "select", "400",
          "STRING combined score: 150 low, 400 medium, 700 high, 900 highest.", choices = c("150", "400", "700", "900"), app = FALSE),
      opt("pathway_network_type", "Network type", "select", "functional",
          "functional: all evidence (co-expression, databases, text mining, experiments); physical: proteins in the same complex only.",
          choices = c("functional", "physical"), app = FALSE),
      opt("pathway_network_max", "Most proteins in the network", "number", 300,
          "The hits with the smallest p are drawn when there are more (Up and Down together).", min = 10, max = 2000, step = 10, app = FALSE)))
  )
}

# STRING categories for the switches above
pathway_categories <- function(s) c(
  if (as_flag(s$pathway_go)) c("Process", "Function", "Component"), if (as_flag(s$pathway_kegg)) "KEGG",
  if (as_flag(s$pathway_reactome)) "RCTM", if (as_flag(s$pathway_wiki)) "WikiPathways", if (as_flag(s$pathway_keywords)) "Keyword",
  if (as_flag(s$pathway_string_clusters)) "NetworkNeighborAL", if (as_flag(s$pathway_domains)) c("Pfam", "InterPro", "SMART"),
  if (as_flag(s$pathway_other)) c("COMPARTMENTS", "TISSUES", "DISEASES", "HPO", "PMID"))
PATHWAY_CATEGORY_NAMES <- c(Process = "GO process", Function = "GO function", Component = "GO component", KEGG = "KEGG", RCTM = "Reactome",
  WikiPathways = "WikiPathways", Keyword = "UniProt keyword", NetworkNeighborAL = "STRING cluster", Pfam = "Pfam", InterPro = "InterPro",
  SMART = "SMART", COMPARTMENTS = "Compartments", TISSUES = "Tissues", DISEASES = "Diseases", HPO = "Human phenotype", PMID = "Publications")
