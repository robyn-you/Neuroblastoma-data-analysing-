# Integrative Bioinformatic Dissection of ADRN-MES Transcriptional Plasticity in Neuroblastoma

Research project for the Master of Professional Engineering (Biomedical), University of Technology Sydney.
Supervisor: Daniel Carter.

## Overview

This project identifies transcription factors (TFs) that regulate Adrenergic (ADRN) and
Mesenchymal (MES) cell identity in neuroblastoma, by integrating:

- **Aim 1** — genome-scale CRISPR-Cas9 dependency screens (DepMap) to identify
  state-selective gene dependencies among neuroblastoma cell lines.
- **Aim 2** — clinical transcriptomic and survival data (KOCAK cohort) to evaluate
  the prognostic significance of candidate TFs in patients.

Results from both aims are combined into a single evidence-ranked TF list (see
`scripts/04_crossdataset_integration.R`).

## Repository structure

```
.
data_old/     # supervisor-provided raw data export (not committed — see Data Availability)
data_raw/     # supervisor-provided cleaned data export (not committed — see Data Availability)
scripts/      # Analysis scripts, run in numbered order
01_reconcile_data.R
# Dataset reconciliation and quality control
 
02a_aim1_differential_dependency.R
# Primary Aim 1 analysis: binary ADRN versus MES comparison using limma
 
02b_aim1_sensitivity_nonstringent.R
# Sensitivity analysis using the non-stringent ADRN/MES classification
 
02d_aim1_welch_test.R
# Alternative differential dependency analysis using Welch's t-test
 
02e_aim1_wilcoxon_test.R
# Alternative differential dependency analysis using the Wilcoxon rank-sum test
 
02f_aim1_continuous_score.R
# Secondary exploratory analysis using the continuous ADRN-MES score.
# This analysis was added after the primary binary comparison
# identified no FDR-significant TFs. Comparison of limma,
# Welch, and Wilcoxon results suggested that modelling cell
# state as a binary ADRN/MES classification may have reduced
# sensitivity by excluding HYBRID samples and the underlying
# ADRN-MES continuum.
 
03_aim2_survival_analysis.R
# Survival analysis in the KOCAK cohort
 
04_crossdataset_integration.R
# Integration of dependency and survival evidence
 
outputs/ # Generated figures and result tables (not committed;
# regenerated from scripts)
 
renv.lock # Locked package versions for reproducibility
 
README.md
```

## Reproducing this analysis

1. Clone this repository.
2. Open R (or RStudio) with the repository as the working directory.
3. Install [renv](https://rstudio.github.io/renv/) if not already installed:
   ```r
   install.packages("renv")
   ```
4. Restore the exact package versions used in this project:
   ```r
   renv::restore()
   ```
5. Place the supervisor-provided data files into `data_raw/` (see Data Availability
   below) — filenames must match exactly what each script expects (see each
   script's header comments).
6. Run the scripts in order:
   ```r
   source("scripts/01_reconcile_data.R")
   source("scripts/02_aim1_differential_dependency.R")
   source("scripts/02b_aim1_sensitivity_nonstringent.R")
   source("scripts/03_aim2_survival_analysis.R")
   source("scripts/04_crossdataset_integration.R")
   ```
   Each script writes its results and figures to `outputs/`.

