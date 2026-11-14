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

## Ranked analysis (GSEA, added 2026-11-13)
Every tested protein ranked by sign(FC) x -log10(p) (no cutoff), fgsea (multilevel) on the chosen databases from STRING's
annotation file of the species (`<data>/.dpdp_cache/<taxon>.protein.enrichment.terms.v12.0.txt.gz`, downloaded once; KEGG is not
in STRING's downloads), gene sets of 10-500 tested proteins (options). Report NES (+ up, - down), FDR and the leading-edge
proteins. GSEA finds pathways where many proteins shift a little; say so when it finds terms ORA misses (few hits) or vice versa.

## GO with the hierarchy (topGO, added 2026-11-13)
topGO weight01 (the method of the app's ViSEAGO GO analysis): Up and Down separately against the tested proteins, for BP, MF and
CC, GO annotations from STRING's species file; terms with weight01 p <= 0.01 (not FDR-corrected, as topGO's authors advise; the
classic Fisher p is shown next to it) and >= the minimum hits. Similar terms are grouped (GOSemSim, Wang similarity >= 0.7, no
organism database needed); the most significant term represents its group. 11517 (KO_37_v_KO_30, Down BP): STRING ORA 34 GO
process terms, topGO 9 specific ones (nuclear migration, spindle localization, chromosome condensation, transcription by RNA pol II).
Report the group representatives; mention the grouped terms only when they add meaning.

## GO slim profile and term drill-down (added 2026-11-13)
- GO slim: the GO Consortium's generic slim (broad categories; a protein counts in every category above its terms). Per comparison
  and ontology, the share of the annotated tested, Up and Down proteins in each category; Fisher (Up or Down v tested), BH. A
  descriptive overview for reports ("what kinds of proteins changed"); the specific answer is in topGO / ORA.
- Drill-down (page only): clicking any term (ORA, GSEA, GO, slim) or network cluster shows its tested proteins on the volcano and
  in a heatmap of the comparison samples, plus a member table. Members come from STRING's annotation file of the same version as
  the API; for KEGG (not in the downloads) only the hits are known ("partial").
- Annotation files follow the API's STRING version (`protein.enrichment.terms.v<version>`, cached in `<data>/.dpdp_cache`).

## Charts for publication (all with PNG / SVG download)
- GO: group representatives as bars (-log10 weight01 p, "+n" = similar terms), or a treemap of groups and terms (REVIGO-style).
- Over-representation: top terms as a diverging bar chart (Down left, Up right, -log10 FDR) or a dot plot (hit ratio, size =
  hits, colour = FDR); follows the source filter.
- GSEA: NES bar chart (top up / down), the table, and the classic running enrichment-score plot of a clicked gene set.
- Overview across comparisons: dot plot of the top terms of every comparison (ORA: colour = direction x -log10 FDR, size = hits;
  GSEA: colour = NES, size = -log10 FDR; GO: the group representatives of biological process).

## Output
- `review/11_pathway.html`: per comparison and direction the enriched terms (source, hits / term size, fold enrichment, FDR,
  proteins), filter by source or text; the network of the hits (Up red, Down blue, size = fold change; proteins without
  interactions hidden, a switch shows them; colour by cluster), the clusters (Louvain on the confidence scores) with what each is
  enriched for, and a link to the same network on string-db.org. Click a term or a cluster to highlight its proteins and list them
  (FC, p, call, Explore link); click a protein in the network to open it in the Explore tab.
- `Pathway/Pathway_STRING_<norm>.xlsx`: About (versions, settings), C<n> Up / C<n> Down (terms with their proteins), C<n> GSEA,
  C<n> GO topGO, C<n> GO slim, C<n> network. `review/11_pathway_drill_<n>.json.gz`: the drill-down data per comparison.
- `steps/11_pathway/summary.json` (counts, top terms) and `results.json` (everything the page shows); STRING id map cached per
  species and STRING version (`string_ids_<taxon>_<version>.csv`).

## Read the result
Report per comparison: how many Up / Down were tested, the top terms per direction (fold enrichment and FDR, not only FDR), the
biggest clusters and what they are enriched for. Say when nothing is enriched, and what that means for the design (e.g. the
triple-spike's human proteins don't change, so its few human hits enrich nothing). Broad terms ("Phosphoprotein", "Cytoplasm") are
expected with large hit lists; point to the specific ones.

## Next
Report (not built yet).
