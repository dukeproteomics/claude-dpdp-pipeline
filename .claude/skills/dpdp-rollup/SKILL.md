---
name: dpdp-rollup
description: Step 8 of the DPDP proteomics pipeline (same as the Shiny app's Rollup tab "Apply Rollup", Protein output). Use after dpdp-impute. Rolls imputed precursors up to proteins (sum, median, median_polish, mean, iq_maxlfq, topn) and to peptides (sum), then calculates group CVs and makes the CV bar plot and the protein bar plot for each normalization. The user approves the settings and reviews every plot.
---

# DPDP Step 8 — Rollup

**Show the user every graph** and the CV summary; the user decides.

## Settings
| Setting | Choices | App default | 11488 reference |
|---|---|---|---|
| `rollup_method` | **sum, iq_maxlfq** (median, median_polish, mean, topn fail in the app; the script refuses them) | sum | **iq_maxlfq** |
| `rollup_topn` | 1–5 (only for topn, which doesn't work) | 3 | 3 |
| `maxlfq_scale` | TRUE/FALSE (scale MaxLFQ output to summed protein intensity) | FALSE | **TRUE** |

## Peptide output (phospho / PTM runs)
The unit is the **peptide** (modified sequence + site localization, `feature_key()`), not the protein: precursors roll up to peptides
and nothing is rolled up to proteins. A peptide's precursors are its charge states.
- Methods: `sum` (adds the charge states) or `limpa` after limpa in Impute (DPC-Quant per peptide). iq_maxlfq is refused: it combines
  peptides of a protein. The page form offers only these.
- The page follows the unit: peptide CV by group, precursors (charge states) per peptide, and a **peptide explorer** (charge states +
  sum / DPC-Quant; spike-in proteins show their 2 peptides with the most charge states). Talk about peptides, not proteins, with the analyst.
- Summary fields: `level`, `features` (output rows), `total_cv_pct`, `features_cv_under_20_pct`; `proteins` is the number of protein
  accessions behind them. The old names `protein_total_cv_pct` / `proteins_cv_under_20_pct` are kept with the same values.

## limpa (after `--impute_type=limpa`)
`--rollup_method=limpa` is then the only allowed rollup (`pipeline/limpa.R`). Preview and apply both run DPC-Quant for every normalization
(~20 s per table on 11488, ~60 s on triple-spike): every protein gets a value in every sample plus a standard error, saved in
`steps/08_rollup/limpa_<norm>.rds` for Stats. `protein_<norm>` gets `2^estimate` in the app layout plus a `Detected_Imputed` column
(measured precursors per sample, background not counted), which the app's Stats reads for the missing factor. `peptide_<norm>` is the
sum of measured precursors (it can have gaps). The detection curve comes from measured values only.
`--limpa_background=missing` (default) leaves Spectronaut background values to the detection curve; `use` counts them as measurements.
Triple-spike background export, false call rate E45/E5 and E10/E5: missing 4.9 / 5.6%, use 6.9 / 9.2%, raw export 4.0 / 5.6%
(E. coli log2 FC, true 3.17: +2.82, +2.37, +3.03). Recommend `missing`.
CVs of limpa values are higher than imputed data (11488 SPQC median 8.2% vs 5.4% iq_maxlfq): imputation shrinks spread artificially.

## Run
```bash
cd /home/user/Claude_Shiny
# 1) preview: saves settings, builds the review page (no rollup yet)
Rscript pipeline/steps/08_rollup.R --db=<db> --action=preview [--settings_from=<db>] [--rollup_method=iq_maxlfq --maxlfq_scale=TRUE]
# 2) after the analyst decides:
Rscript pipeline/steps/08_rollup.R --db=<db> --action=apply [--settings_from=<db>] [--rollup_method=...] --decided_by=user --analyst=<name> [--note="<reason>"]
```
Like the app, this also runs the QC statistics (group CVs, `protein_<norm>_cv`, `protein_<norm>_final`, `summary_cv`)
for every normalization in `norm_type`. The impute-only table gives the "no normalization" comparison.

## Review page and decision
Both actions (re)build **`data/<prefix>/review/08_rollup.html`** (rebuild: `Rscript pipeline/review_page.R --db=<db> --step=08_rollup`).
The page has:
- **sum vs iq_maxlfq comparison,** calculated on the page from the imputed precursor tables (MaxLFQ takes under 1 s):
  protein CV curves per group, plus a table of median CV, mean CV and % of proteins under 20% CV, with a list to pick the normalization.
  The page's iq_maxlfq numbers matched the app's run exactly for 11488. For 11488, sum gave slightly lower CVs (median 9.4 / 12.2 / 5.4
  vs 9.8 / 12.4 / 5.7).
- **Precursors per protein,** including how many proteins have one precursor (for those, all methods give the same value).
- **Protein explorer:** precursor profiles plus both rollups for spike-ins, the largest proteins and a random set. The relative view
  centers each line on its own mean.
- **After apply:** the CV table (the app's mean plus median and % under 20%), the CV bar plot, and the protein bar plots.

Report the comparison neutrally, then ask with **AskUserQuestion** (header "Rollup"), no option marked as recommended:
"Which rollup for <prefix>?" Options: "iq_maxlfq, scaled" / "iq_maxlfq, not scaled" / "sum".

## Show the user
From `steps/08_rollup/summary.json`, and **open every file in `plots`**:
1. **CV bar plot** (`CV_barplot.png`): average CV% per group for each normalization (same as the app's QC tab).
   Also give the **median** CV and % of proteins under 20% CV (`tables[].median_cv_by_group`, `proteins_cv_under_20_pct`).
   The app's average is pulled up by a high-CV tail, so the median is usually much lower.
   Pooled QC (SPQC) should have the lowest CV; it measures technical variation only.
2. **Protein bar plots** (`protein_<norm>_barplot.png`): protein totals per sample for each normalization.
3. Protein and peptide counts (should equal the filtered protein count from Step 5).
Then ask the user to accept, or change settings and re-run.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db \
  --tables=protein_<norm>,peptide_<norm>,protein_<norm>_cv,protein_<norm>_final,summary_cv --params=rollup_method,rollup_topn,maxlfq_scale
```

## Next
Step 9 — QC (protein plots and spike-in plots).
