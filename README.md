# EWAS–CLPM Analysis of DNAm and TyG Index

This repository contains the R scripts used for a longitudinal epigenome-wide association study (EWAS) examining the relationship between DNA methylation (DNAm) and the triglyceride-glucose (TyG) index, along with cross-lagged panel modelling (CLPM) and downstream analyses.

---

## Repository Contents

### 1. `longi_lmer.R` — Longitudinal EWAS (Delta Model)
Runs a genome-wide delta-EWAS using linear mixed-effects models (`lme4`/`lmerTest`). For each CpG site, it models the change in methylation (ΔM = M_T1 − M_T0) as a function of the change in TyG index (ΔTyG), adjusting for baseline age and sex, with a random intercept for family clustering. Outputs full results, FDR-significant hits, and the genomic inflation factor (lambda).

### 2. `CLPM.R` — Cross-Lagged Panel Model (CLPM)
Fits cross-lagged panel models using `lavaan` for each CpG site. The model estimates bidirectional prospective associations between baseline DNAm and follow-up TyG (and vice versa), conditioning on the autoregressive stability of each variable. This captures a fundamentally different estimand from the delta-EWAS: the association between a baseline level and a future level after partialling out prior values.

### 3. `MethylGSA_enrichment.R` — Pathway Enrichment Analysis
Performs GO, KEGG, and Reactome enrichment analysis on EWAS-significant CpG sites using `methylGSA`. Uses `methylglm` (logistic regression, recommended) with `methylgometh` (hypergeometric test) as a fallback. Outputs per-database enrichment tables and a combined summary of significant pathways (adjusted p < 0.05).

### 4. `simulate.R` — Simulation Study (Reviewer Diagnostic)
A simulation study conducted in response to a reviewer comment on the consistent sign reversal observed between EWAS and CLPM coefficients across all reported loci. Using data-generating parameters informed by the real cohort (n = 100; autoregressive coefficients derived from the CLPM output), the simulation demonstrates that sign discordance between the two estimands can arise systematically from the structural difference between a difference-score model and a cross-lagged level model, rather than from any coding error. A sensitivity analysis across a grid of autoregressive and effect-size parameters is also included.

---

## Input Data Format

All scripts expect a CSV input with one row per subject and the following key columns:

| Column | Description |
|---|---|
| `unit` | CpG probe ID (e.g. `cg00000029`) |
| `P.Value` | Nominal p-value from the EWAS model |
| `adj.P.Val` | BH-adjusted p-value |
| `<CpGid>_T0` / `<CpGid>_T1` | Methylation values at baseline / follow-up |
| `T0_TyG`, `T1_TyG` | TyG index at baseline / follow-up |
| `T0_age`, `T0_gender` | Baseline covariates |
| `family_ID` | Family identifier for random effects |

---

## Usage

Set the `input_file` and `out_dir` variables at the top of each script to your local paths before running. Scripts are intended to be run in the following order:

```
longi_lmer.R → MethylGSA_enrichment.R
CLPM_1_3.R
simulate.R   (independent; reviewer diagnostic only)
```

---

## Dependencies

```r
install.packages(c("data.table", "lme4", "lmerTest", "lavaan", "MASS", "parallel"))
BiocManager::install(c("methylGSA",
                       "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       "org.Hs.eg.db",
                       "reactome.db"))
```
