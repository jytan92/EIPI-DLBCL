# =============================================================================
# 05_discrimination.R : discrimination and relative model fit of IPI, NCCN-IPI and E-IPI
#   Development cohort               -> Table 2
#   Validation cohort                -> Supplementary Table S8
#   Validation cohort, by country    -> Supplementary Tables S9A (OS) and S9B (PFS)
# Harrell's C (bootstrap 95% CI), paired bootstrap Delta C vs E-IPI (2000 resamples),
# AIC with each score entered as an ordinal variable. Primary cohorts.
# Run time: several minutes (bootstrap). Outputs: output/table2_S8_discrimination.csv,
#           output/S9_by_country.csv
# =============================================================================
source("00_setup.R")

disc <- rbind(cbind(Cohort = "Development", discrimination(COHORTS$primary$dev)),
              cbind(Cohort = "Validation",  discrimination(COHORTS$primary$val)))

val <- COHORTS$primary$val
by_country <- do.call(rbind, lapply(c("Hong Kong", "Malaysia", "Korea", "China"), function(ctry) {
  cbind(Country = ctry, discrimination(val[val$Country == ctry, ]))
}))

cat("\n== Table 2 (development) and S8 (validation) ==\n")
print(disc[, c("Cohort", "Endpoint", "Score", "C-index", "95% CI", "Delta C vs E-IPI", "P", "AIC")], row.names = FALSE)
cat("\n== S9A/B: validation cohort by country ==\n")
print(by_country[, c("Country", "N", "Events", "Endpoint", "Score", "C-index", "Delta C vs E-IPI",
                     "Delta C 95% CI", "P", "AIC")], row.names = FALSE)

save_table(disc, "table2_S8_discrimination")
save_table(by_country, "S9_by_country")
