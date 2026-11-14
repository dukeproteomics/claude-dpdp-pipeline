---
name: dpdp-phospho
description: Step 12 of the DPDP proteomics pipeline (replaces the Shiny app's MotifX / MoMo tabs). Phospho / PTM runs with Peptide output, after dpdp-stats. Kinase enrichment with the Kinase Library (Cantley lab; 303 Ser/Thr and the Tyr kinases), enrichment logos and motif-x motifs for the Up and the Down phosphosites of each comparison against the tested sites. The analyst sets the options on the page, runs it, and accepts the results.
---

# DPDP Step 12 — Phospho motifs and kinases

Chosen with the user (2026-11-13): sequences from UniProt (only accessions sent), the `kinase-library` Python package installed
(web app venv via `web/install.sh`; dev `.venv-kl/`; or `DPDP_KL_PYTHON`).

## What it does
- **Sites:** every localized site (`Local2 = "Y"`) of the tested phosphopeptides, one row per site (`Protein_PTM_Loc`, e.g. S48);
  a site on several peptides takes the most significant one. Up / Down: the Stats calls, or custom cutoffs.
- **Windows:** 15 residues (site +-7) from the UniProt sequence of the first accession; the site must be that residue (11517:
  141,757 of 141,908 sites matched). Cached in `steps/12_phospho/uniprot.fasta`.
- **Kinase Library** (`pipeline/kinase_library_run.py`): every site scored against all kinases (percentile); a site is a kinase's
  substrate when the kinase is in its top 15 (percentile rank, the authors' method); Fisher's exact test (BH) of each kinase's
  substrate share in the Up and in the Down sites v the unchanged sites (the package's differential-phosphorylation design).
  S/T sites -> Ser/Thr kinases, Y sites -> Tyr kinases. Per site: its top 3 kinases.
- **Enrichment logo** (pLogo-style): binomial significance of each residue at each position, Up or Down S/T sites v every tested
  S/T site; letters where |log10 p| passes Bonferroni.
- **motif-x** (re-implemented, `pipeline/phospho_motifs.R`): Up or Down v the tested sites, one central residue at a time,
  minimum 20 occurrences, p < 1e-6 (options).

Human kinase specificities: for other organisms (yeast, Cryptococcus) report kinases as motif classes (CDK / MAPK =
proline-directed, PKA / AKT / RSK = basophilic, CK2 / GRK = acidophilic), not as the organism's own kinases.

## Run
```bash
Rscript pipeline/steps/12_phospho.R --db=<db> --action=preview [--from_page] [--phospho_<option>=...]    # ~1-2 min
Rscript pipeline/steps/12_phospho.R --db=<db> --action=apply --from_page --decided_by=user                # accept (same settings)
```
Page choices: `preview` -> `--action=preview --from_page` (then wait again); `apply` -> `--action=apply --from_page --decided_by=user`;
`other` -> discuss. Phospho / PTM runs with Peptide output only (the step stops otherwise).

## Output
`review/12_phospho.html` (kinase volcano + table per comparison, Up / Down, Ser/Thr / Tyr; logo heatmap; motifs; the site table
with windows and top kinases, loaded per comparison from `review/12_phospho_sites_<n>.json.gz`), `Phospho/Phospho_motifs_<norm>.xlsx`
(About, kinase tables, motifs, sites), `steps/12_phospho/summary.json` (top kinases and motifs), `results.json`, `kinase_library.json`.

## Read the result
Report per comparison the top kinases (adj p and the share of sites), whether the logo and motif-x agree with them (e.g. P at +1 /
`sP` with CDK / MAPK; acidic +1..+3 / `s..E`, `sD` with CK2), and the biggest differences between Up and Down. 11517 (Cryptococcus,
30 -> 37 C): Up = p38 / JNK / ERK (proline-directed stress MAPK, adj p ~1e-22, `sP`), Down = CK2 (`s..E`, `sD`).

## Downloads (user, 2026-11-13)
CSV for the summary, kinases (every kinase, not only the enriched ones shown), motifs and sites (every site of the filter, not the
first 300); the kinase volcano and logo charts have "Data (CSV)"; the About section links the Excel workbook.
