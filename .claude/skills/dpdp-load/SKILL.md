---
name: dpdp-load
description: Step 1 of the DPDP proteomics pipeline (same as the Shiny app's Load tab). Use when starting a new proteomics project from a search-engine report (Spectronaut, Proteome Discoverer, FragPipe) and a sample description (design) Excel file. Creates the project folder and SQLite .db, detects data source and format, checks sample IDs, and reports per-sample missingness.
---

# DPDP Step 1 — Load

## Inputs
- **Design file** (.xlsx, sheet `SampleList`): columns `Set, Raw_Order, Replicate, QC Spike Level, ID, Label, Group`.
- **Data file**: Spectronaut report (.tsv), or PD PSMs export (expects a matching `PeptideGroups` file next to it).
- **Prefix**: project name. Convention: `<project#>_<MMDDYY>`, e.g. `11488_110326`. Ask the user if unsure.

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/01_load.R --design=<xlsx> --data=<report> --prefix=<prefix> [--outdir=data] [--overwrite] \
  [--sp_fill=raw|background] [--sp_qvalue=0.01]
```
- **Spectronaut export type** (user decision 2026-11-04; `pipeline/spectronaut.R`). Accepted: **raw** (no background, no imputation)
  and **background-signal** exports, one TSV. Load checks before creating anything:
  - `EG.IsImputed = True` anywhere: global / run-wise imputation export, **refused** (the flag also marks measured values).
  - Run-wise `EG.Qvalue` present: values with q > `sp_qvalue` (or no q) are background signal, so the export is `background`;
    none means `raw`. `--sp_fill` may state the type; a mismatch stops Load.
  - No run-wise `EG.Qvalue`: Load stops and asks. If the user says it is raw, rerun with `--sp_fill=raw`. A background export without
    q-values is refused (re-export with `EG.Qvalue`).
  - Output: param `sp_fill_mode` (raw/background) and table `precursor_fill` (same rows as `precursor_raw`, 1 = background signal).
  - **Steps 4-7 refuse background runs for now** (`stop_if_background()` in `common.R`) until they use the mask.
- Use plain `Rscript` (not `--vanilla`); packages are in `/home/R_Packages`.
- Output folder: `data/<prefix>/` with `<prefix>.db`, `Backup/` (design copy + archived app code), `QC/`, `Error/`, etc.
- **`--overwrite` deletes the existing project folder.** Only use it with the user's OK.

## Read the result
Read `data/<prefix>/steps/01_load/summary.json` and report:
1. **Detected:** `data_source` (SP/PD), `raw_data_format` (precursor/protein), `data_table_format` (short = one column per sample), `ptm`, `data_output` (Protein, or Peptide when PTM).
   PTM is detected from `EG.PTMProbabilities` columns. If that doesn't match what the user expects, stop and ask.
2. **Counts:** rows, proteins, peptides, precursors, duplicate precursors (should be 0).
3. **Sample ID check:** `app_check` must be TRUE for every sample, otherwise **stop**. The design order or the export is wrong.
   `full_id_check` FALSE for **SPQC** samples is fine (pooled QC, order within SPQC doesn't matter).
4. **Missingness per sample:** point out samples well above the others. Don't decide anything here; it feeds TIC (Step 2) and QC.
5. **Spectronaut export** (`spectronaut_fill`): type, measured, filled by Spectronaut, empty. Say which it is and what it means.
6. **Organisms:** the main species plus a few others is normal.

## Review page and decision
The step builds **`data/<prefix>/review/01_load.html`**: detected format, counts, the sample sheet with each sample's matched data
column (a warning appears if any did not match), missing % per sample, and proteins per organism (spike-ins and contaminants show here).
Ask with **AskUserQuestion** (header "Load"): "Does the design match the data for <prefix>?" Options: "Yes, continue" / "No, fix the sample sheet".
Record it: `Rscript pipeline/accept.R --db=<db> --step=01_load --decided_by=user --analyst=<name> [--note="<reason>"]`

## Lab rules (from the user)
- **Yeast proteins are a spike-in**, not contamination (design column `QC Spike Level`). Don't suggest removing them.
- Human, bovine, porcine, keratin and similar are usually contaminants. Normalization later excludes `trypsin|keratin|casein`.
- SPQC = pooled QC; replicate order within SPQC doesn't matter.

## Check against a Shiny reference run (if one exists)
```bash
Rscript pipeline/compare.R --run=data/<prefix>/<prefix>.db --ref=<shiny>.db --tables=design,precursor_raw \
  --params=data_source,raw_data_format,data_table_format,ptm,data_output
```
After Load only, the reference `design` has 6 extra columns (`PrimaryGroup, Count, colorlist, Header1-3`). They are added in Step 3 (Parameters), so that difference is expected.

## Next
Step 2 — TIC.
