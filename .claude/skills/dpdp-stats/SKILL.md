---
name: dpdp-stats
description: Step 10 of the DPDP proteomics pipeline (same as the Shiny app's Stats Setup, Comparisons and Graphs tabs). Use after dpdp-qc. Sets stats options and comparisons (e.g. B_v_A plus SPQC), runs the app's stat calculation (peptide refilter, rollup, fold change, p-value, missing factor, Up/Down calls), then makes static versions of the Stats graphs (bar, box, PCA, sample distance, heatmap, volcano). The user approves the settings and reviews every plot.
---

# DPDP Step 10 — Stats

**Show the user every graph and the significance counts**; the user decides.

## Settings (Stats Setup + Comparisons)
| Setting | 11488 reference | Notes |
|---|---|---|
| `stats_norm_type` | sltmm | which normalization to test |
| `comp_spqc` | SPQC if the design has one, else `0` | pooled QC group carried into the tables and plots (`0` = none; no SPQC is fine) |
| comparisons | `--comp1_N=B --comp1_D=A` → `B_v_A` | several groups per side: comma separated. Names must be ≤ 31 characters (Excel) |
| `pvalue_cutoff` / `foldchange_cutoff` | 0.05 / 1.5 | |
| `checkbox_adjpval` / `padjust_options` | TRUE / fdr | **only adds an `_adjpval` column** |
| `checkbox_filter_adjpval` | FALSE | **decides whether Up/Down calls use the adjusted p.** FALSE = raw p |
| `missing_factor` | 0.6 | proteins with `_mf` below this are **removed** from the comparison table |
| `peptide_refilter` | TRUE | re-filters precursors for the comparison's samples, then rolls up again |
| others (`peptide_*`, `stats_*_cv_filter*`, `checkbox_cohensd*`, `checkbox_limmapvalue`) | off | see `10_stats.R` spec |

`--settings_from=<db>` copies the settings **and** the comparisons from a previous run.

## limpa runs
After limpa in Impute and Rollup, Stats sets `peptide_refilter=FALSE` (refilter re-rolls precursors and would bypass DPC-Quant) and refuses
`--pair_comp`. The app's `stat_calc_bg` still builds the tables (layout, CVs, missing factor filter, flags); then `limpa_stats()` replaces
`<comp>_FC`, `_FC2`, `_pval`, `_adjpval` and `Stats` with limpa's dpcDE (limma with precision weights) on the comparison's samples.
FC comes from the model (difference of log2 means), not the ratio of means. Adjusted p is over all tested proteins, as in the app.
The summary's `test` field says which test was used: report it. 11488 B vs A: 229 calls by raw p, 1 by adjusted p (duke run: 221 / 72;
the duke adjusted calls are mostly false positives, step 07b).

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/10_stats.R --db=<db> --action=preview [--settings_from=<db>] [--comp1_N=B --comp1_D=A] [--name=value ...]
Rscript pipeline/steps/10_stats.R --db=<db> --action=apply   [...same...] --decided_by=user --analyst=<name> [--note="<reason>"]
Rscript pipeline/steps/10b_stats_graphs.R --db=<db> [--comps=B_v_A] [--include_spqc=TRUE] [--volcano_highlight=text,...]
```
Graphs are static PNGs made from the app's own data prep (`create_plot_bg`, `create_volcano_bg`), with the plot code
from `Shiny_Interactive.R`. PCA_3D is not made. Graph sets: each comparison's samples, plus "all" (comparisons + SPQC).

## Review page and decisions
`10_stats.R` (preview, apply, accept), `10b` and `10c` all rebuild **`data/<prefix>/review/10_stats.html`**
(rebuild: `Rscript pipeline/review_page.R --db=<db> --step=10_stats`). The page has:
- **Settings and comparisons** in plain words, including which p the Up/Down calls use.
- **Significance per comparison:** calls by raw p next to what would pass with adjusted p, the smallest adjusted p, and raw p counts.
  A warning appears when the calls use raw p and none pass FDR, with the number expected by chance.
- **Volcano explorer:** switch raw/adjusted p, the p cutoff and the FC cutoff, and leave out proteins measured in one group only.
  Counts and the list of significant proteins update. It uses the same rule as the app (it reproduced 221 = 133 / 88 for 11488).
  Proteins measured in only one group are ringed.
- **Raw p histogram** with a flat level estimated from p > 0.5 and a rough count of real differences (11488: about 387 above chance at p < 0.05).
- **PCA** (sample sets, nearest sample), the static graphs from 10b, the protein plots from 10c, and the Excel paths.

Decisions with **AskUserQuestion**, no option marked as recommended:
1. Before apply (header "Stats"): "Run stats for <comparisons> with these settings?" Options: "Run as shown" / "Change settings"
   (p type for calls, cutoffs, missing factor) / "Change comparisons".
2. After apply and the graphs (header "Results"): "Accept the stats results for <prefix>?" Options: "Accept" / "Change cutoffs and re-run" /
   "Change comparisons". Then record it:
   `Rscript pipeline/steps/10_stats.R --db=<db> --action=accept --decided_by=user --analyst=<name> [--note="<reason>"]`
**Contaminants are flagged, never removed** (lab rule, 2026-11-04; `pipeline/contaminants.R`). The flags are digestion trypsin (pig/cow),
human keratins (when the sample isn't human), and the QC spike-ins (casein, ADH). They appear as diamonds on the volcano, a Flag column in
the significant list, `flagged` in `summary.json`, and a **"Flagged proteins"** sheet that 10c adds to the Excel file (the app's sheets are
unchanged). Report how many significant calls are flagged, and whether the flagged proteins lean to one group: that suggests that group
has more contamination (11488: 7 human keratins, all Down, i.e. higher in A).

## Show the user
From `steps/10_stats/summary.json` and `steps/10_stats/graphs.json`, and **open every plot**:
1. **Significance:** the app's Up/Down counts and which p-value they use (`pvalue_used_for_calls`).
   Always also give `significant_if_adjusted_p` and `pvalue_distribution` (raw p < 0.05/0.01/0.001; smallest adjusted p).
   If the calls use raw p and nothing passes FDR, **say so clearly**.
2. **Missing factor:** how many proteins were removed (`removed_by_missing_factor`).
3. **Large fold changes:** check `Detected_Imputed` for the top hits. Very large FCs (≥ 10×) usually mean
   detected in one group and imputed (bottom-X) in the other. Report how many significant proteins are like that.
4. **Volcano** (`Stats_<comp>_volcano.png`): red = Up, blue = Down; dotted lines = FC and p cutoffs.
5. **PCA** (`Stats_<set>_pca2d.png`): do the groups separate on PC1 or PC2? SPQC should sit tightly in the middle (a pool).
   Point out samples far from their group (`graphs.json` scores, `nearest_sample`).
6. **Sample distance** (`_cluster.png`) and **heatmap** (`_heatmap.png`, column dendrogram order): which samples cluster together.
7. **Bar / box** plots: totals and distributions of the samples used.
Then ask the user to accept, or change settings or comparisons and re-run.

## Excel output and Protein Plots (step 10c)
```bash
Rscript pipeline/steps/10c_stats_output.R --db=<db> [--excel=TRUE] [--comp=B_v_A] [--proteins=<acc>,...] [--top=6] [--plot_spqc=TRUE]
```
- **Excel:** runs the app's `stats_Final_Excel_bg` → `data/<prefix>/<norm>/Final_<norm>_stats.xlsx`
  (Table1 Sample Info [created empty by the app], Table2 Raw Precursor, Table3 Raw Peptide, Table4 Normalized, Table5+ one per comparison)
  and `<prefix>_Parameters.xlsx`. Tell the user the file paths. Remind them the Up/Down column uses raw p if `checkbox_filter_adjpval` is off.
- **Protein Plots:** per-sample protein bar plot (title shows FC, p, adj p, call) and per-peptide grouped bar plot (log2 mean ± SD by group).
  Default = top 6 significant proteins by p. Ask the user which proteins they want. Good extras: constant spike (ADH P00330),
  and on/off hits, to show that A is at the bottom-X level.
  The app's own version (`create_stats_oneprotein_plots_bg`) fails when a group name matches an annotation column
  (`A` → `Accession`) and writes a debug file `z1b24` into the app folder, so the pipeline uses its own data prep.
- The Stats **Data** tab is a filtered table browser; the Excel comparison sheet has the same data.

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db \
  --tables=stats_comp,protein_<norm>_<comp>,protein_<norm>_<comp>_final,peptide_impute_<norm>_<comp>_final \
  --params=stat_norm,comp_spqc,pvalue_cutoff,checkbox_adjpval,padjust_options,foldchange_cutoff,missing_factor,peptide_refilter
```

## Next
Stats Data / Protein Plots / Excel output. Pathway (11) and Report (13) are skipped for now, so the first pass ends here.
