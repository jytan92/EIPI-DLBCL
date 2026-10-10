# =============================================================================
# 08_secondary.R : secondary analysis in patients alive with at least 6 months of follow-up
# (development n=710, validation n=495). The score, cut-offs and risk groups are not
# re-derived; the analyses of 05 and 04 are repeated in this population.
#   (a) Discrimination and relative model fit (as in 05)  -> Supplementary Table S11A
#   (b) N (%) and 4-year OS/PFS by E-IPI risk group        -> Supplementary Table S11B
# Run time: about 1 minute (bootstrap). Outputs: output/S11A_secondary_discrimination.csv,
#           output/S11B_secondary_survival.csv
# =============================================================================
source("00_setup.R")

sec <- COHORTS$secondary

# ---- (a) discrimination ---------------------------------------------------------
disc <- rbind(cbind(Cohort = sprintf("Development (n=%d)", nrow(sec$dev)), discrimination(sec$dev)),
              cbind(Cohort = sprintf("Validation (n=%d)",  nrow(sec$val)), discrimination(sec$val)))

# ---- (b) 4-year survival by E-IPI risk group --------------------------------------
surv <- do.call(rbind, lapply(names(sec), function(ck) {
  x <- sec[[ck]]
  do.call(rbind, lapply(1:4, function(g) {
    xs <- x[x$eipi_grp == g, ]
    data.frame(Cohort = sprintf("%s (n=%d)", if (ck == "dev") "Development" else "Validation", nrow(x)),
               `E-IPI group` = GROUP_LABELS[g],
               `N (%)` = sprintf("%d (%.1f%%)", nrow(xs), 100 * nrow(xs) / nrow(x)),
               `4-y OS` = sprintf("%.1f%%", km_at(xs, "OS")["est"]),
               `4-y PFS` = sprintf("%.1f%%", km_at(xs, "PFS")["est"]), check.names = FALSE)
  }))
}))

cat("\n== S11A: discrimination, patients alive with >=6 months of follow-up ==\n")
print(disc[, c("Cohort", "Endpoint", "Score", "C-index", "95% CI", "Delta C vs E-IPI", "P", "AIC")], row.names = FALSE)
cat("\n== S11B: 4-year survival by E-IPI group, patients alive with >=6 months of follow-up ==\n")
print(surv, row.names = FALSE)

save_table(disc, "S11A_secondary_discrimination")
save_table(surv, "S11B_secondary_survival")
