# =============================================================================
# 02_measix_prognostic.R : prognostic value of mEASIX
#   (a) mEASIX per doubling (log2), unadjusted and adjusted for IPI, NCCN-IPI or sex
#   (b) mEASIX tiers: 4-year OS/PFS, log-rank test and HRs   (Supplementary Table S1A)
#   (c) mEASIX tier within each IPI stratum                   (Supplementary Table S2)
# Primary cohorts. Outputs: output/measix_continuous.csv, output/S1A_measix_tiers.csv,
#                           output/S2_measix_within_ipi.csv
# =============================================================================
source("00_setup.R")

ENDPOINTS <- c("OS", "PFS")

# Cox model returning HR (95% CI) and Wald P for the first term
cox_first_term <- function(x, endpoint, rhs) {
  f <- as.formula(paste("surv_obj(x, endpoint) ~", rhs))
  s <- summary(coxph(f, data = x, ties = "efron"))
  c(hr = s$conf.int[1, "exp(coef)"], lo = s$conf.int[1, "lower .95"],
    hi = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"])
}

# ---- (a) mEASIX as a continuous variable (per doubling) ----------------------
# IPI and NCCN-IPI entered as their four-level risk group, coded 1-4 (linear)
models <- c("Unadjusted"               = "log2_measix",
            "Adjusted for IPI"         = "log2_measix + ipi_grp",
            "Adjusted for NCCN-IPI"    = "log2_measix + nccn_grp",
            "Adjusted for IPI and sex" = "log2_measix + ipi_grp + male")

continuous <- do.call(rbind, lapply(names(COHORTS$primary), function(ck) {
  x <- COHORTS$primary[[ck]]
  do.call(rbind, lapply(ENDPOINTS, function(ep) {
    do.call(rbind, lapply(names(models), function(m) {
      r <- cox_first_term(x, ep, models[[m]])
      data.frame(Cohort = ck, Endpoint = ep, Model = m,
                 `HR per doubling (95% CI)` = fmt_hr(r["hr"], r["lo"], r["hi"]),
                 P = fmt_p(r["p"]), check.names = FALSE)
    }))
  }))
}))

# ---- (b) mEASIX tiers ------------------------------------------------------------
tier_labels <- c("Low (<0.70)", "Intermediate (0.70-1.50)", "High (>1.50)")

tiers <- do.call(rbind, lapply(names(COHORTS$primary), function(ck) {
  x <- COHORTS$primary[[ck]]
  x$tier_f <- factor(x$tier, levels = 0:2)
  do.call(rbind, lapply(ENDPOINTS, function(ep) {
    fit <- summary(coxph(surv_obj(x, ep) ~ tier_f, data = x, ties = "efron"))
    lr  <- survdiff(surv_obj(x, ep) ~ tier_f, data = x)
    p_lr <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
    do.call(rbind, lapply(0:2, function(k) {
      km <- km_at(x[x$tier == k, ], ep)
      hr <- if (k == 0) "1.00 (reference)" else
        fmt_hr(fit$conf.int[k, 1], fit$conf.int[k, 3], fit$conf.int[k, 4])
      data.frame(Cohort = ck, Endpoint = ep, `mEASIX tier` = tier_labels[k + 1],
                 N = sum(x$tier == k), `4-year % (95% CI)` = sprintf("%.1f (%.1f-%.1f)", km[1], km[2], km[3]),
                 `HR (95% CI)` = hr, P = if (k == 0) "" else fmt_p(fit$coefficients[k, "Pr(>|z|)"]),
                 `Log-rank P` = if (k == 0) fmt_p(p_lr) else "", check.names = FALSE)
    }))
  }))
}))

# ---- (c) mEASIX tier within each IPI stratum (tier as ordinal 0/1/2) -------------
ipi_labels <- c("Low (0-1)", "Low-intermediate (2)", "High-intermediate (3)", "High (4-5)")

within_ipi <- do.call(rbind, lapply(ENDPOINTS, function(ep) {
  x <- COHORTS$primary$dev
  do.call(rbind, lapply(1:4, function(g) {
    xs <- x[x$ipi_grp == g, ]
    r <- cox_first_term(xs, ep, "tier")
    ev <- if (ep == "OS") sum(xs$os_event) else sum(xs$pfs_event)
    data.frame(Endpoint = ep, `IPI stratum` = ipi_labels[g], N = nrow(xs), Events = ev,
               `HR per tier (95% CI)` = fmt_hr(r["hr"], r["lo"], r["hi"]), P = fmt_p(r["p"]),
               check.names = FALSE)
  }))
}))

cat("\n== mEASIX per doubling ==\n");            print(continuous, row.names = FALSE)
cat("\n== mEASIX tiers (S1A) ==\n");             print(tiers, row.names = FALSE)
cat("\n== mEASIX within IPI strata, development (S2) ==\n"); print(within_ipi, row.names = FALSE)

save_table(continuous, "measix_continuous")
save_table(tiers,      "S1A_measix_tiers")
save_table(within_ipi, "S2_measix_within_ipi")
