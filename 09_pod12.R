# =============================================================================
# 09_pod12.R : early relapse (POD12) — progression or relapse within 12 months of the end
# of first-line treatment. POD12 status is coded in the data files (0, 1 or "Exclude").
# Not evaluable ("Exclude"): death without documented progression, or alive without
# progression with <12 months of follow-up after the end of treatment.
#   (a) Patients and early-relapse rate by E-IPI score, with cumulative sensitivity,
#       specificity, PPV and NPV at each cut-point (>= score)        -> S12A
#   (b) AUC of continuous E-IPI, IPI and NCCN-IPI scores; DeLong test vs E-IPI
#   (c) E-IPI >= 4 cut-point: sensitivity, specificity, PPV, NPV (Wilson 95% CI)
#   (d) Sensitivity analysis: not-evaluable patients who died counted as failures
#   (b)-(d) -> S12B. Primary cohorts. Outputs: output/S12A_pod12_by_score.csv,
#              output/S12B_pod12_discrimination.csv
# =============================================================================
source("00_setup.R")

fmt_ci <- function(x) sprintf("%.3f (%.3f-%.3f)", x[2], x[1], x[3])

auc_compare <- function(y, x) {
  rocs <- lapply(c(`E-IPI` = "eipi", IPI = "ipi", `NCCN-IPI` = "nccn"),
                 function(s) roc(y, x[[s]], direction = "<", levels = c(0, 1), quiet = TRUE))
  do.call(rbind, lapply(names(rocs), function(nm) {
    p <- if (nm == "E-IPI") "Ref" else
      fmt_p(roc.test(rocs[["E-IPI"]], rocs[[nm]], method = "delong", paired = TRUE)$p.value)
    data.frame(Score = nm, `AUC (95% CI)` = fmt_ci(as.numeric(ci.auc(rocs[[nm]], method = "delong"))),
               `DeLong P vs E-IPI` = p, check.names = FALSE)
  }))
}

wilson <- function(k, n) {
  ci <- prop.test(k, n, correct = FALSE)$conf.int
  sprintf("%d/%d (%.1f%%; %.1f-%.1f%%)", k, n, 100 * k / n, 100 * ci[1], 100 * ci[2])
}

S12A <- NULL; S12B <- NULL
for (ck in c("dev", "val")) {
  x_all <- COHORTS$primary[[ck]]
  label <- if (ck == "dev") "Development" else "Validation"
  x <- x_all[!is.na(x_all$pod12), ]                 # evaluable patients
  y <- x$pod12

  # ---- (a) by E-IPI score --------------------------------------------------------
  for (s in sort(unique(x$eipi))) {
    n_s <- sum(x$eipi == s); e_s <- sum(y[x$eipi == s])
    flag <- x$eipi >= s
    tp <- sum(flag & y == 1); fp <- sum(flag & y == 0); fn <- sum(!flag & y == 1); tn <- sum(!flag & y == 0)
    op <- if (s == min(x$eipi)) rep("-", 4) else
      sprintf("%.2f", c(tp / (tp + fn), tn / (tn + fp), tp / (tp + fp), tn / (tn + fn)))
    S12A <- rbind(S12A, data.frame(Cohort = label, `E-IPI score` = s,
      `Patients, n (%)` = sprintf("%d (%.1f)", n_s, 100 * n_s / nrow(x)),
      `Early relapse, n (%)` = sprintf("%d (%.1f)", e_s, 100 * e_s / n_s),
      `Cut-point` = paste0(">=", s), Sensitivity = op[1], Specificity = op[2], PPV = op[3], NPV = op[4],
      check.names = FALSE))
  }

  # ---- (b) AUC and DeLong --------------------------------------------------------
  auc <- auc_compare(y, x)

  # ---- (c) E-IPI >= 4 -------------------------------------------------------------
  hi <- x$eipi >= 4
  tp <- sum(hi & y == 1); fp <- sum(hi & y == 0); fn <- sum(!hi & y == 1); tn <- sum(!hi & y == 0)

  # ---- (d) sensitivity analysis: deaths without progression counted as failures ---
  died_ne <- x_all$pod12_status == "Exclude" & x_all$os_event == 1
  xs <- x_all[!is.na(x_all$pod12) | died_ne, ]
  ys <- as.integer(xs$pod12_status == "1" | (xs$pod12_status == "Exclude" & xs$os_event == 1))
  auc_s <- auc_compare(ys, xs)

  S12B <- rbind(S12B, data.frame(Cohort = label, Measure = c(
      "Evaluable patients / POD12 events", paste("AUC,", auc$Score),
      "E-IPI >=4: sensitivity", "E-IPI >=4: specificity", "E-IPI >=4: PPV", "E-IPI >=4: NPV",
      "Early relapse rate, E-IPI >=4", "Early relapse rate, E-IPI <=3",
      "Sensitivity analysis: patients / events", paste("Sensitivity analysis: AUC,", auc_s$Score)),
    Value = c(sprintf("%d / %d", nrow(x), sum(y)), auc$`AUC (95% CI)`,
              wilson(tp, tp + fn), wilson(tn, tn + fp), wilson(tp, tp + fp), wilson(tn, tn + fn),
              sprintf("%d/%d (%.1f%%)", tp, tp + fp, 100 * tp / (tp + fp)),
              sprintf("%d/%d (%.1f%%)", fn, fn + tn, 100 * fn / (fn + tn)),
              sprintf("%d / %d", nrow(xs), sum(ys)), auc_s$`AUC (95% CI)`),
    `DeLong P vs E-IPI` = c("", auc$`DeLong P vs E-IPI`, rep("", 6), "", auc_s$`DeLong P vs E-IPI`),
    check.names = FALSE))
}

cat("\n== S12A: POD12 by E-IPI score (evaluable patients) ==\n"); print(S12A, row.names = FALSE)
cat("\n== S12B: discrimination and E-IPI >=4 cut-point ==\n");   print(S12B, row.names = FALSE)

save_table(S12A, "S12A_pod12_by_score")
save_table(S12B, "S12B_pod12_discrimination")
