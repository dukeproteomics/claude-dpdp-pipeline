---
name: dpdp-pathway
description: Step 11 of the DPDP proteomics pipeline (improved version of the Shiny app's Pathway tabs). Use after dpdp-stats. Over-representation through STRING (GO, KEGG, Reactome, WikiPathways, UniProt keywords, STRING clusters) for the Up and the Down proteins of each comparison separately, against the proteins tested in that comparison, plus the STRING network of the hits with clusters. The analyst sets the options on the page, runs it, and accepts the results.
---

# DPDP Step 11 — Pathway (STRING)

Best practice chosen with the user (2026-11-13; the STRING website is fine to use, only accessions are sent):
- **Up and Down tested separately**, never pooled.
- **Background = the proteins tested in that comparison** (the measured proteome), not the genome. Against the genome, anything
  abundant looks enriched.
- Other species (spike-ins, mixed designs), contaminants and QC spike-ins left out (switches). Phospho / Peptide output: a protein
  is a hit if any of its tested sites is; background = proteins with a tested site.
- Fewer than 3 mapped hits: not tested (the page says so). With few hits, say ORA has little power (a ranked analysis, GSEA, is
  planned next).
- STRING version's stable address (e.g. version-12-5) and every setting are recorded (summary, Excel "About" sheet).

## Options (`pipeline/options/11_pathway.R`, the page's Options panel)
Organism (NCBI taxon; detected from "org=" / "OS=" in the descriptions: human 9606, mouse 10090, rat 10116, yeast 4932, E. coli
511145, Cryptococcus JEC21 214684), only this organism, leave out flagged; hits = Stats calls or custom (p type, p, FC); databases
(GO, KEGG, Reactome, WikiPathways, keywords, STRING clusters on; domains, compartments/tissues/diseases/publications off); FDR 0.05,
minimum 3 hits per term; network confidence 400 (150 / 400 / 700 / 900), functional or physical, at most 300 proteins (smallest p).

## Run
```bash
Rscript pipeline/steps/11_pathway.R --db=<db> --action=preview [--from_page] [--pathway_<option>=...]   # asks string-db.org, ~1 min per comparison
Rscript pipeline/steps/11_pathway.R --db=<db> --action=apply --from_page --decided_by=user               # accept (same settings as the last run)
```
Page choices: `preview` ("Run with these settings") -> `--action=preview --from_page`, then wait again; `apply` ("Accept these
settings") -> `--action=apply --from_page --decided_by=user`; `other` -> discuss. Accept with different settings is refused.

## Output
- `review/11_pathway.html`: per comparison and direction the enriched terms (source, hits / term size, fold enrichment, FDR,
  proteins), filter by source or text; the network of the hits (Up red, Down blue, size = fold change; proteins without
  interactions hidden, a switch shows them; colour by cluster), the clusters (Louvain on the confidence scores) with what each is
  enriched for, and a link to the same network on string-db.org. Click a term or a cluster to highlight its proteins and list them
  (FC, p, call, Explore link); click a protein in the network to open it in the Explore tab.
- `Pathway/Pathway_STRING_<norm>.xlsx`: About (versions, settings), C<n> Up / C<n> Down (terms with their proteins), C<n> network.
- `steps/11_pathway/summary.json` (counts, top terms) and `results.json` (everything the page shows); STRING id map cached per
  species and STRING version (`string_ids_<taxon>_<version>.csv`).

## Read the result
Report per comparison: how many Up / Down were tested, the top terms per direction (fold enrichment and FDR, not only FDR), the
biggest clusters and what they are enriched for. Say when nothing is enriched, and what that means for the design (e.g. the
triple-spike's human proteins don't change, so its few human hits enrich nothing). Broad terms ("Phosphoprotein", "Cytoplasm") are
expected with large hit lists; point to the specific ones.

## Next
Report (not built yet). Planned for this step: ranked analysis (GSEA, no cutoff), grouping of redundant terms, term drill-down
(volcano, heatmap).
