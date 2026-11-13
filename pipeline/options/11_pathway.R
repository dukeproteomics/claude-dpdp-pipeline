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
    opt_section("Ranked analysis (GSEA)", note = "Every tested protein ranked by sign(FC) x -log10(p): no cutoff, finds pathways where many proteins shift a little. fgsea (Korotkevich et al.), the databases chosen above from STRING's annotation of the species (KEGG is not in STRING's downloads).", list(
      opt("pathway_gsea", "Run GSEA", "checkbox", TRUE, "", app = FALSE),
      opt("pathway_gsea_min", "Smallest gene set", "number", 10, "Gene sets with fewer tested proteins are skipped.", depends = "pathway_gsea", min = 3, step = 1, app = FALSE),
      opt("pathway_gsea_max", "Largest gene set", "number", 500, "Very broad sets (e.g. 'Cytoplasm') are skipped.", depends = "pathway_gsea", min = 20, step = 10, app = FALSE))),
    opt_section("GO with the hierarchy (topGO)", note = "topGO weight01 (Alexa et al. 2006; the method of the app's ViSEAGO GO analysis): a parent term whose signal comes from a child term is discounted, so the list is shorter and more specific. Then similar terms are grouped (GOSemSim, Wang similarity). GO annotations: STRING's file of the species.", list(
      opt("pathway_topgo", "Run GO (topGO)", "checkbox", TRUE, "Needs the GO switch above.", app = FALSE),
      opt("pathway_topgo_p", "p cutoff (weight01)", "number", 0.01,
          "topGO's authors advise not to correct weight01 p for multiple testing (the tests are not independent); 0.01 is their usual cutoff.",
          depends = "pathway_topgo", min = 0, max = 0.1, step = "any", app = FALSE),
      opt("pathway_topgo_sim", "Group terms more similar than", "number", 0.7, "GOSemSim Wang similarity (0 to 1): 0.7 groups close relatives, 0.5 groups more broadly.",
          depends = "pathway_topgo", min = 0.3, max = 0.95, step = 0.05, app = FALSE))),
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
