# =============================================================================
# 07_calibration.R : calibration of the E-IPI in the external validation cohort
# The development cohort provides the predictions; nothing is re-estimated in validation.
#   (a) Predicted vs observed 4-year OS/PFS by E-IPI risk group
#       predicted = development Kaplan-Meier estimate for that group
#       observed  = validation Kaplan-Meier estimate for that group (95% CI, log-log)
#       -> data for Supplementary Figure 3 (plotted in MedCalc)
#   (b) Observed / expected events within 48 months
#       expected = sum over validation patients of the development group-specific
#                  Nelson-Aalen cumulative hazard at min(follow-up, 48 months)
#   (c) Calibration slope: Cox model in validation with the development linear predictor
#       (E-IPI group coefficients from the development cohort) as the only covariate
# Primary cohorts. Outputs: output/calibration_by_group.csv, output/calibration_summary.csv
# =============================================================================
source("00_setup.R")

dev <- COHORTS$primary$dev
val <- COHORTS$primary$val
dev$grp <- factor(dev$eipi_grp, levels = 1:4)
val$grp <- factor(val$eipi_grp, levels = 1:4)
fmt_km <- function(k) sprintf("%.1f (%.1f-%.1f)", k["est"], k["lo"], k["hi"])

# development Nelson-Aalen cumulative hazard as a step function of time
cumhaz_fun <- function(x, ep) {
  f <- survfit(surv_obj(x, ep) ~ 1, ctype = 1)
  stepfun(f$time, c(0, f$cumhaz))
}

by_group <- NULL; summary_rows <- NULL
for (ep in c("OS", "PFS")) {
  tt <- if (ep == "OS") "os_time" else "pfs_time"
  ee <- if (ep == "OS") "os_event" else "pfs_event"
  tot_o <- 0; tot_e <- 0
  for (g in 1:4) {
    a <- dev[dev$eipi_grp == g, ]; b <- val[val$eipi_grp == g, ]
    pred <- km_at(a, ep); obs <- km_at(b, ep)
    H <- cumhaz_fun(a, ep)
    o <- sum(b[[ee]] == 1 & b[[tt]] <= HORIZON)
    e <- sum(H(pmin(b[[tt]], HORIZON)))
    tot_o <- tot_o + o; tot_e <- tot_e + e
    by_group <- rbind(by_group, data.frame(
      Endpoint = ep, `E-IPI group` = GROUP_LABELS[g], `Validation n` = nrow(b),
      `Predicted 4-y % (development KM, 95% CI)` = fmt_km(pred),
      `Observed 4-y % (validation KM, 95% CI)` = fmt_km(obs),
      `Observed - predicted (points)` = sprintf("%+.1f", obs["est"] - pred["est"]),
      `Observed / expected events <=48 mo` = sprintf("%d / %.1f", o, e),
      `O/E ratio` = sprintf("%.2f", o / e), check.names = FALSE))
  }
  # calibration slope
  m_dev <- coxph(surv_obj(dev, ep) ~ grp, data = dev, ties = "efron")
  val$lp <- as.vector(model.matrix(~ grp, data = val)[, -1] %*% coef(m_dev))
  s <- summary(coxph(surv_obj(val, ep) ~ lp, data = val, ties = "efron"))
  summary_rows <- rbind(summary_rows, data.frame(
    Endpoint = ep,
    `Calibration slope (95% CI)` = fmt_hr(s$coefficients[1, "coef"],
                                          s$coefficients[1, "coef"] - 1.96 * s$coefficients[1, "se(coef)"],
                                          s$coefficients[1, "coef"] + 1.96 * s$coefficients[1, "se(coef)"]),
    `Observed / expected events <=48 mo` = sprintf("%d / %.1f", tot_o, tot_e),
    `Overall O/E ratio` = sprintf("%.2f", tot_o / tot_e), check.names = FALSE))
}

cat("\n== Predicted vs observed 4-year survival by E-IPI group (Supplementary Figure 3 data) ==\n")
print(by_group, row.names = FALSE)
cat("\n== Calibration slope and overall observed/expected ratio ==\n")
print(summary_rows, row.names = FALSE)

save_table(by_group, "calibration_by_group")
save_table(summary_rows, "calibration_summary")
