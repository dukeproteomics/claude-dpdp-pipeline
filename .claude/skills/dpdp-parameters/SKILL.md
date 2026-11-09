---
name: dpdp-parameters
description: Step 3 of the DPDP proteomics pipeline (same as the Shiny app's Parameters tab, "Accept Parameters"). Use after dpdp-load. Sets analysis options (Full/Primary groups, Protein/Peptide output, PTM, peptide selection), builds sample groups, prepares and renames sample columns, runs the sample ID check, and makes the start bar, box and intensity histogram plots.
---

# DPDP Step 3 — Parameters

(Step 2, TIC, is skipped for now.)

## Options (same as the app's Parameters tab)
Every option, its default, choices and meaning are in `pipeline/options/03_parameters.R` (also the page's Options panel).
| Option | Choices | Default | Notes |
|---|---|---|---|
| `primary_group` | Full, Primary | Full | Primary: Filter (X% of a group, misaligned, CV) and Impute use the part of `Group` before the first `_` (WT_30 + WT_37 -> WT). Stats comparisons still use the full group names (as in the app), so Primary gives pooled groups for filtering/imputing while keeping every comparison choice. Not offered when no group has a `_`. |
| `data_output` | Protein, Peptide | detected at Load | Peptide is detected for PTM data. Peptide output also works without PTM (a set that is not enriched but should be reported per peptide; tested on 11531, 2026-11-09). |
| `ptm` | TRUE/FALSE | detected at Load | Only with Peptide output (set to FALSE for Protein output). Not offered when the report has no `EG.PTMProbabilities` column. |
| `ptm_grep` | text | Phospho | only with `ptm` |
| `ptm_local` | 0-1 | 0.75 | PTM localization minimum; only with `ptm` |
| `multi_tmt` | TRUE/FALSE | FALSE | TMT only: hidden (TMT is not supported) |
| `peptide_select` | Razor, Unique, Shared | Razor | Proteome Discoverer only: hidden for Spectronaut |
| `use_isoform` | TRUE/FALSE | FALSE | Proteome Discoverer only; not supported in the pipeline yet |

Use the defaults unless the user says otherwise; the 11488 reference run used all defaults.
The analyst usually sets these on the page (Options panel); in the chat, tell the user which options you will use and confirm them.

**Group order (pipeline fix, 2026-11-09):** each group must be one block of samples (sample_groups start..end). The app's sort check
only asks that each sample has a same-group neighbour, so e.g. CreNeg_Biotin x3, CrePos_Biotin x6, CreNeg_NoBiotin x3 passes with
Primary although CreNeg is split, and Filter/Impute then use the wrong samples. When the groups (Primary or Full) are not each one
block, the step sorts the samples by group (sheet order kept within a group, SPQC last), warns, and sets `design_sorted` in the summary.

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/03_parameters.R --db=data/<prefix>/<prefix>.db [--from_page] [--force] [--primary_group=Full] [--data_output=Peptide] ...
```
- `--from_page`: the values the analyst sent with **Run with these settings** (`review/pending_decision.json`); archives that request.
- The step stops if later steps were already run (their results would no longer match). `--force` re-runs anyway; then redo those steps.
- Accept takes the settings of the last run: `Rscript pipeline/accept.R --db=<db> --step=03_parameters --decided_by=user [--from_page]`.
  Option values that differ from the last run (from the page or the command line) are refused: re-run the step with them first.
Supported so far: Spectronaut precursor data (Protein output, or Peptide output with or without PTM).
Other data types stop with a clear message.

## Read the result
Read `data/<prefix>/steps/03_parameters/summary.json`, then **look at the three plots** in `plots`:
1. **Stop conditions:** `sample_id_alert` or `sample_error` TRUE means the design order or the export is wrong. Stop and tell the user.
2. **Counts:** starting precursors, peptides and proteins (`counts_start`); 8 info columns is normal for SP precursor data.
3. **Bar plot (total intensity per sample):** look for samples far from the others. `per_sample$total_vs_median` > 1.5 or < 0.67 is a warning (low load or a bad injection). Small differences are normal; normalization corrects them.
4. **Box plot (log2 intensity per sample):** medians should line up. A shifted box suggests a loading difference; a narrow or skewed box suggests a bad run.
5. **Histogram (log2 intensity, all samples):** should be one roughly normal peak.
   - A separate low-intensity peak or shoulder suggests noise. Point this out when recommending Step 4 (Noise).
   - Lines: orange = bottom-X percentile (used by bottom-X imputation), black = mean ± SD, green = misaligned intensity cutoff (mean + `intensity_cutoff_sd` × SD, used by the filter step).
6. Compare `missing` per sample with total intensity. **High missingness with normal total intensity** suggests fewer identifications, not lower load. Mention it and carry it to the Filter and QC steps.

## Review page and decision
The step builds **`data/<prefix>/review/03_parameters.html`**: the **Parameters** panel (every option that applies to this data
source, editable; options that never apply are listed as not shown, options that depend on another are greyed), start counts, total
intensity per sample relative to the median sample, and the app's start bar/box plots and histogram.
Page choices: `preview` = **Run with these settings** (re-run the step with `--from_page`, then wait again); `apply` = **Accept these
settings** (`accept.R --step=03_parameters --from_page --decided_by=user`); `other` = free text, discuss it.
In the chat: ask with **AskUserQuestion** (header "Parameters"): "Keep these parameters for <prefix>?" Options: "Yes, continue" /
"Change settings" (re-run the step with `--name=value`, then ask again). Accept: `Rscript pipeline/accept.R --db=<db> --step=03_parameters --decided_by=user --analyst=<name>`

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db --tables=design,sample_groups,precursor_start \
  --params=primary_group,data_output,ptm,peptide_select,info_col_precursor,meta_precursor_start,meta_peptide_start,meta_protein_start,intensity_mean,intensity_sd,intensity_cutoff
```
Expected difference: `colorlist` in `design` and `sample_groups`. Group colors are random on every run (`randomcoloR`).

## Next
Step 4 — Noise.
