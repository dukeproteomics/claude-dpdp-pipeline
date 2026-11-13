# Step 12 - Phospho (pipeline; replaces the app's MotifX / MoMo tabs): kinase enrichment with the Kinase Library, enrichment logos
# and motif-x motifs for the Up and the Down sites of each comparison, against the tested sites. Phospho / PTM runs only.
phospho_options <- function(p = NULL) list(
  opt_section("Sites", note = "Sites: every localized site of the tested phosphopeptides; windows of 15 residues from UniProt (first accession of the protein group).", list(
    opt("phospho_localized", "Localized sites only", "checkbox", TRUE,
        sprintf("Only peptides whose sites pass the localization cutoff (Parameters: %s).", p$ptm_local %||% 0.75), app = FALSE),
    opt("phospho_hits", "Which sites are hits", "select", "stats calls", "stats calls: the Up / Down calls of the Stats step. custom: your own cutoffs.",
        choices = c("stats calls", "custom"), app = FALSE),
    opt("phospho_p_type", "p value", "select", "adjusted", "", choices = c("adjusted", "raw"), depends = "phospho_hits=custom", app = FALSE),
    opt("phospho_p_cutoff", "p cutoff", "number", 0.05, "", depends = "phospho_hits=custom", min = 0, max = 1, step = "any", app = FALSE),
    opt("phospho_fc_cutoff", "Fold change cutoff", "number", 1.5, "", depends = "phospho_hits=custom", min = 1, step = "any", app = FALSE))),
  opt_section("Kinase Library", note = "Johnson et al., Nature 2023 (303 Ser/Thr kinases); Yaron-Barir et al., Nature 2024 (Tyr kinases). Human kinase specificities: for other organisms read them as motif classes (e.g. proline-directed CDK / MAPK, basophilic PKA / AKT, acidophilic CK2).", list(
    opt("phospho_kl_rank", "Top kinases per site (percentile rank)", "number", 15,
        "A site counts as a kinase's substrate when the kinase is among its top this many (the authors' recommended method and default).",
        min = 1, max = 50, step = 1, app = FALSE),
    opt("phospho_kl_fdr", "FDR cutoff for the kinase tables", "number", 0.1, "Kinases with a higher adjusted p are greyed on the volcano.", min = 0, max = 1, step = "any", app = FALSE))),
  opt_section("motif-x", note = "Schwartz & Gygi, Nat Biotechnol 2005 (re-implemented); foreground = Up or Down sites, background = the tested sites, one central residue at a time.", list(
    opt("phospho_motif_min", "Minimum occurrences", "number", 20, "A motif needs at least this many foreground sites (lower it for small sets).", min = 3, step = 1, app = FALSE),
    opt("phospho_motif_sig", "Significance", "number", 1e-6, "Binomial p needed to fix a residue (motif-x default 0.000001).", min = 0, max = 0.01, step = "any", app = FALSE)))
)
