# Data dictionary

Column names expected in `data/development.xlsx` and `data/validation.xlsx` (one row per patient). Other columns are ignored.

| Column name | Description | Coding |
|---|---|---|
| `SN` | Serial number | integer |
| `Country` | Centre (validation file only) | China, Hong Kong, Korea, Malaysia |
| `OS (months)` | Time from diagnosis to death or last follow-up | months |
| `OS outcome` | Death from any cause | 1 = died, 0 = censored |
| `PFS (months)` | Time from diagnosis to progression, relapse, death or last follow-up | months |
| `PFS outcome` | Progression, relapse or death | 1 = event, 0 = censored |
| `POD12` | Progression or relapse within 12 months of the end of first-line treatment | 1, 0, or `Exclude` (not evaluable) |
| `IPI raw value` (development) / `IPI raw` (validation) | IPI score | 0–5 |
| `IPI-categorised` | IPI risk group | 1 low (0–1), 2 low-intermediate (2), 3 high-intermediate (3), 4 high (4–5) |
| `NCCN raw` | NCCN-IPI score | 0–8 |
| `NCCN category` | NCCN-IPI risk group | 1 low (0–1), 2 low-intermediate (2–3), 3 high-intermediate (4–5), 4 high (6–8) |
| `Modified EASIX` | mEASIX (unrounded) | continuous |
| `3 tier classification (<0.7, 0.7-1.5, >1.5)` | mEASIX tier (after rounding to 2 dp) | 0, 1, 2 |
| `E-IPI raw` | E-IPI score | 0–7 |
| `E-IPI categorised` | E-IPI risk group | 1 low, 2 low-intermediate, 3 high-intermediate, 4 high |
| `LDH ratio` | Serum LDH / institutional upper limit of normal | ratio |
| `Creatinine (mg/dl)` | Serum creatinine | mg/dL |
| `PLT` | Platelet count | ×10⁹/L |
| `Gender code (Male = 1, Female =2)` | Sex | 1 male, 2 female |
| `Age` | Age at diagnosis | years |
| `Age>60=1`, `Stage 3 or 4=1`, `ECOG2-4=1`, `LDHraised=1`, `>1Extranodalsite=1` | IPI factors | 1 present, 0 absent |
| `COO (0=GCB)` (development) / `COO (0=GCB, 1=non-GCB, NR= not reported)` (validation) | Cell of origin (Hans) | 0 GCB, 1 non-GCB, NR |
| `Treatment (RCHOP =1, REPOCH=2)` | First-line regimen | 1 R-CHOP, 2 R-EPOCH |

`00_setup.R` checks that mEASIX equals LDH ratio × creatinine × 100 / platelets, and that the tier, E-IPI score and E-IPI group in the file follow the scoring rule; it stops with an error if not.
