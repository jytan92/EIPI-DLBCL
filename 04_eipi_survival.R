# =============================================================================
# 04_eipi_survival.R : risk-group hazard ratios and 4-year survival
#   (a) Cox HRs by E-IPI risk group, low as reference
#       development = Supplementary Table S1B; validation = Figure 3 legend
#   (b) N (%) and 4-year OS/PFS (Kaplan-Meier) by IPI, NCCN-IPI and E-IPI risk group
#       development = Supplementary Table S4; validation = Supplementary Table S7
# Primary cohorts. Outputs: output/S1B_eipi_hr.csv, output/S4_S7_survival_by_group.csv
# =============================================================================
source("00_setup.R")

SCORES <- c("IPI" = "ipi_grp", "NCCN-IPI" = "nccn_grp", "E-IPI" = "eipi_grp")

# ---- (a) E-IPI risk-group hazard ratios ----------------------------------------
hr_by_group <- function(cohorts) {
  do.call(rbind, lapply(names(cohorts), function(ck) {
    x <- cohorts[[ck]]
    x$grp <- factor(x$eipi_grp, levels = 1:4)
    do.call(rbind, lapply(c("OS", "PFS"), function(ep) {
      s  <- summary(coxph(surv_obj(x, ep) ~ grp, data = x, ties = "efron"))
      lr <- survdiff(surv_obj(x, ep) ~ grp, data = x)
      data.frame(Cohort = ck, Endpoint = ep, `E-IPI group` = GROUP_LABELS,
                 `HR (95% CI)` = c("1.00 (reference)", fmt_hr(s$conf.int[, 1], s$conf.int[, 3], s$conf.int[, 4])),
                 P = c("", fmt_p(s$coefficients[, "Pr(>|z|)"])),
                 `Log-rank P` = c(fmt_p(pchisq(lr$chisq, 3, lower.tail = FALSE)), "", "", ""),
                 check.names = FALSE)
    }))
  }))
}

# ---- (b) N (%) and 4-year survival by risk group, for each score ----------------
survival_by_group <- function(cohorts) {
  do.call(rbind, lapply(names(cohorts), function(ck) {
    x <- cohorts[[ck]]
    do.call(rbind, lapply(1:4, function(g) {
      row <- data.frame(Cohort = ck, `Risk group` = GROUP_LABELS[g], check.names = FALSE)
      for (sc in names(SCORES)) {
        xs <- x[x[[SCORES[sc]]] == g, ]
        row[[paste(sc, "N (%)")]]    <- sprintf("%d (%.1f%%)", nrow(xs), 100 * nrow(xs) / nrow(x))
        row[[paste(sc, "4-y OS")]]   <- sprintf("%.1f%%", km_at(xs, "OS")["est"])
        row[[paste(sc, "4-y PFS")]]  <- sprintf("%.1f%%", km_at(xs, "PFS")["est"])
      }
      row
    }))
  }))
}

hr   <- hr_by_group(COHORTS$primary)
surv <- survival_by_group(COHORTS$primary)

cat("\n== E-IPI hazard ratios (S1B = development; Figure 3 legend = validation) ==\n"); print(hr, row.names = FALSE)
cat("\n== Survival by risk group (S4 = development; S7 = validation) ==\n");           print(surv, row.names = FALSE)

save_table(hr,   "S1B_eipi_hr")
save_table(surv, "S4_S7_survival_by_group")
