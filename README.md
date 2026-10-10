# EASIX-IPI (E-IPI) in diffuse large B-cell lymphoma — analysis code

R code reproducing all analyses in:

> Tan JY, *et al.* The EASIX-IPI (E-IPI): a refined prognostic index for survival and early relapse in DLBCL — a multicentre Asian study. *International Journal of Cancer* (2026).

## The E-IPI score

| Step | Rule |
|---|---|
| Modified EASIX (mEASIX) | LDH ratio (observed / upper limit of normal) × creatinine (mg/dL) × 100 / platelets (×10⁹/L), **rounded to two decimal places** |
| mEASIX tier | low < 0.70 = 0 points; intermediate 0.70–1.50 = 1 point; high > 1.50 = 2 points |
| E-IPI | IPI (0–5) + mEASIX tier points (0–2), range 0–7 |
| Risk groups | low 0–1; low-intermediate 2–3; high-intermediate 4; high 5–7 |

## Data

Patient-level data are **not** included in this repository. They are available from the corresponding author on reasonable request, subject to institutional approvals. To run the code, place the two files in a folder named `data/`:

- `data/development.xlsx` — development cohort (n = 753)
- `data/validation.xlsx` — external validation cohort (n = 533)

The required columns are listed in `DATA_DICTIONARY.md`.

## How to run

Requires R (≥ 4.1) and the packages `readxl`, `survival`, `maxstat` and `pROC`:

```r
install.packages(c("readxl", "survival", "maxstat", "pROC"))
setwd("path/to/this/folder")
source("run_all.R")      # about 4 minutes; tables are written to output/
```

Each script can also be run on its own (`source("05_discrimination.R")`); every script starts by sourcing `00_setup.R`, which reads the data, applies the scoring rule, checks it against the coding in the files and sets the random seed (2026) used for all bootstrap analyses (2,000 resamples).

## Scripts and the results they produce

| Script | Manuscript output |
|---|---|
| `00_setup.R` | Data import, scoring rule, cohort definitions, shared functions |
| `01_table1.R` | Table 1 (P values for development vs validation) |
| `02_measix_prognostic.R` | mEASIX per doubling (Results); Supplementary Tables S1A, S2A, S2B |
| `03_score_construction.R` | Derivation of cut-offs and risk groups (Results); Supplementary Tables S3A, S3B, S3C |
| `04_eipi_survival.R` | Supplementary Table S1B; Figure 3 legend; Supplementary Tables S4, S7 |
| `05_discrimination.R` | Table 2; Supplementary Tables S8, S9A, S9B |
| `06_reclassification.R` | Reclassification (Results); Supplementary Tables S5, S6, S10 |
| `07_calibration.R` | Calibration slope and observed/expected ratio; data for Supplementary Figure 3 |
| `08_secondary.R` | Secondary analysis (patients alive with ≥6 months of follow-up); Supplementary Tables S11A, S11B |
| `09_pod12.R` | Early relapse (POD12); Supplementary Tables S12A, S12B |
| `10_km_validation.R` | Figure 3 (Kaplan–Meier curves, validation cohort) |
| `11_calibration_plot.R` | Supplementary Figure 3 (calibration plot) |
| `run_all.R` | Runs scripts 01–09 and records `sessionInfo()` |

## Licence

MIT
