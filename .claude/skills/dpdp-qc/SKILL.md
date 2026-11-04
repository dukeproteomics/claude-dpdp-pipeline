---
name: dpdp-qc
description: Step 9 of the DPDP proteomics pipeline (same as the Shiny app's QC tab plot buttons). Use after dpdp-rollup. Makes spike-in plots by design spike level (casein) and per-sample plots for constant spike-ins (yeast ADH) and other chosen proteins, for each normalization. Compares measured spike ratios with the design. The user reviews every plot.
---

# DPDP Step 9 — QC

**Show the user every graph and the spike numbers.** Describe what matches or doesn't match the design; the user decides what it means.

## Lab spike-ins (from the user)
| Spike | Accessions | How it's spiked | What to expect |
|---|---|---|---|
| Bovine casein | P02662 (alpha-S1), P02663 (alpha-S2) | **different levels**, from the design column `QC Spike Level` | intensity ratios follow the levels (e.g. 1 : 1.5 : 2) |
| Yeast ADH | P00330 (ADH1), P00331 (ADH2) | **same volume** added before injection, loaded on the EvoSep tip with the sample | flat across samples **without normalization** |

User rules (2026-11-04):
- **Casein:** assume the sample description sheet is correct. If the measured levels don't match, the problem is in the
  sample, not the sheet. Report which samples deviate.
- **ADH:** it is a volumetric injection control, not part of the sample's protein, so judge it on the **un-normalized**
  (impute) values. Sample normalization can make ADH less flat; that is expected, not a problem with normalization.

These are the script defaults. Use `--spike_accession` / `--constant_accession` for other projects.

## Run
```bash
cd /home/user/Claude_Shiny
Rscript pipeline/steps/09_qc.R --db=<db> [--spike_accession=P02662,P02663] [--constant_accession=P00330,P00331] \
  [--protein_accession=<acc>,...] [--norms=sltmm,impute]
```
Runs for every normalization in `norm_type`. The impute-only table shows the "before normalization" picture.

## Review page and decision
Running the step builds **`data/<prefix>/review/09_qc.html`** (rebuild: `Rscript pipeline/review_page.R --db=<db> --step=09_qc`). The page has:
- **Casein, measured vs expected per sample,** with a normalization list and a tolerance slider (default ±50%). Study samples are
  compared with their design level; **pool samples (SPQC) are compared with the mean of the study samples**, because that is what
  pooling gives. Also: measured spike precursors per sample, and flags ("check sample").
  For 11488: A_2 is high (8.3×), B_1 and B_3 are low (0.37×, 0.40×), and SPQC is 0.94 of the study mean.
- **ADH per sample,** with the un-normalized (impute) view by default, plus CV and max/min.
- **The app's QC plots.**

QC changes no data. Ask with **AskUserQuestion** (header "QC"): "Accept the QC for <prefix> and continue to Stats?" Options:
"Accept" / "Accept and flag samples" (list the samples the page flags; let them edit) / "Stop here". Then record it:
```bash
Rscript pipeline/steps/09_qc.R --db=<db> --action=accept --decided_by=user --analyst=<name> [--flag_samples=A_2,B_1,B_3] [--note="<reason>"]
```
Flagged samples are recorded for the report; they are not removed.

## Show the user
From `steps/09_qc/summary.json`, and **open every file in `plots`**:
1. **Spike plot** (`QC_Spike_<norm>_barplot.png`): bars per sample, grouped by design spike level.
   Report `by_level` (measured vs expected ratio to the lowest level) and `per_sample` (error %).
   - Look for single samples that don't follow their level.
   - SPQC is a pool, so its spike level reflects the pool, not a separate addition.
     If one sample stands out, check whether SPQC ≈ (that sample)/(number of pooled samples).
   - Very low spike intensities with many missing precursors are near the detection limit, and imputation affects them.
2. **Constant spike plots** (`P00330_<norm>_barplot.png` …): CV across samples and max/min, before (impute) and after normalization.
   Lead with the un-normalized (impute) CV: that is the injection check. A higher CV after normalization is expected,
   because ADH is added by volume, not in proportion to sample protein.
3. Any `--protein_accession` plots.

## Next
Step 10 — Stats.
