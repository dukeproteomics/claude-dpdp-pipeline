---
name: dpdp-parameters
description: Step 3 of the DPDP proteomics pipeline (same as the Shiny app's Parameters tab, "Accept Parameters"). Use after dpdp-load. Sets analysis options (Full/Primary groups, Protein/Peptide output, PTM, peptide selection), builds sample groups, prepares and renames sample columns, runs the sample ID check, and makes the start bar, box and intensity histogram plots.
---

# DPDP Step 3 — Parameters

(Step 2, TIC, is skipped for now.)

## Options (same as the app's widgets)
| Option | Choices | Default | Notes |
|---|---|---|---|
| `primary_group` | Full, Primary | Full | Primary = use the part of `Group` before the first `_` for filter/impute grouping |
| `data_output` | Protein, Peptide | detected at Load | Peptide is detected automatically for PTM data |
| `ptm` | TRUE/FALSE | detected at Load | detected from `EG.PTMProbabilities` |
| `ptm_grep` | text | Phospho | |
| `ptm_local` | number | 0.75 | PTM localization minimum |
| `multi_tmt` | TRUE/FALSE | FALSE | SPQC-normalized TMT sets |
| `peptide_select` | Razor, Unique, Shared | Razor | |
| `use_isoform` | TRUE/FALSE | FALSE | not supported in the pipeline yet |

Use the defaults unless the user says otherwise; the 11488 reference run used all defaults.
Before running, tell the user which options you will use and confirm them.

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/03_parameters.R --db=data/<prefix>/<prefix>.db [--primary_group=Full] [--peptide_select=Razor] ...
```
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
The step builds **`data/<prefix>/review/03_parameters.html`**: the settings in plain words, start counts, total intensity per sample
relative to the median sample, and the app's start bar/box plots and histogram.
Ask with **AskUserQuestion** (header "Parameters"): "Keep these parameters for <prefix>?" Options: "Yes, continue" / "Change settings"
(re-run with `--name=value`). Record it: `Rscript pipeline/accept.R --db=<db> --step=03_parameters --decided_by=user --analyst=<name>`

## Check against a Shiny reference run
```bash
Rscript pipeline/compare.R --run=<db> --ref=<shiny>.db --tables=design,sample_groups,precursor_start \
  --params=primary_group,data_output,ptm,peptide_select,info_col_precursor,meta_precursor_start,meta_peptide_start,meta_protein_start,intensity_mean,intensity_sd,intensity_cutoff
```
Expected difference: `colorlist` in `design` and `sample_groups`. Group colors are random on every run (`randomcoloR`).

## Next
Step 4 — Noise.
