# =============================================================================
# 11_calibration_plot.R : Supplementary Figure 3 — calibration of the E-IPI in the
# external validation cohort. Predicted 4-year survival (development Kaplan-Meier,
# applied without re-estimation) vs observed 4-year survival (validation Kaplan-Meier,
# 95% CI) for each E-IPI risk group. A = OS, B = PFS.
# Outputs: output/Calibration_EIPI_validation.png (600 dpi), .pdf and .tiff
# =============================================================================
source("00_setup.R")

dev <- COHORTS$primary$dev
val <- COHORTS$primary$val
COLS <- c("#1F4E9A", "#2E9E3E", "#E07B00", "#C62828")   # low, low-int, high-int, high

panel <- function(ep, title) {
  pts <- t(sapply(1:4, function(g) c(km_at(dev[dev$eipi_grp == g, ], ep)["est"],
                                     km_at(val[val$eipi_grp == g, ], ep))))
  colnames(pts) <- c("pred", "obs", "lo", "hi")
  # calibration slope and O/E (as in 07_calibration.R)
  d2 <- dev; v2 <- val
  d2$grp <- factor(d2$eipi_grp, levels = 1:4); v2$grp <- factor(v2$eipi_grp, levels = 1:4)
  b  <- coef(coxph(surv_obj(d2, ep) ~ grp, data = d2, ties = "efron"))
  v2$lp <- as.vector(model.matrix(~ grp, data = v2)[, -1] %*% b)
  sl <- summary(coxph(surv_obj(v2, ep) ~ lp, data = v2))$coefficients
  tt <- if (ep == "OS") "os_time" else "pfs_time"; ee <- if (ep == "OS") "os_event" else "pfs_event"
  o <- 0; e <- 0
  for (g in 1:4) {
    a <- dev[dev$eipi_grp == g, ]; bb <- val[val$eipi_grp == g, ]
    f <- survfit(surv_obj(a, ep) ~ 1, ctype = 1); H <- stepfun(f$time, c(0, f$cumhaz))
    o <- o + sum(bb[[ee]] == 1 & bb[[tt]] <= HORIZON); e <- e + sum(H(pmin(bb[[tt]], HORIZON)))
  }

  par(mar = c(4.5, 4.8, 2.5, 1), mgp = c(2.8, 0.7, 0), las = 1)
  plot(NA, xlim = c(10, 100), ylim = c(10, 100), xaxs = "i", yaxs = "i",
       xlab = "Predicted 4-year survival (%), development cohort",
       ylab = "Observed 4-year survival (%), validation cohort", main = title, font.main = 1)
  abline(h = seq(20, 100, 20), v = seq(20, 100, 20), col = "grey92")
  abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
  arrows(pts[, "pred"], pts[, "lo"], pts[, "pred"], pts[, "hi"], angle = 90, code = 3,
         length = 0.05, col = COLS, lwd = 1.6)
  lines(pts[, "pred"], pts[, "obs"], col = "grey55", lwd = 1)
  points(pts[, "pred"], pts[, "obs"], pch = 21, bg = COLS, col = "black", cex = 1.6)
  legend("topleft", legend = c(GROUP_LABELS, "Perfect calibration"), pch = c(rep(21, 4), NA),
         pt.bg = c(COLS, NA), lty = c(rep(NA, 4), 2), col = c(rep("black", 4), "grey40"),
         bty = "n", cex = 0.85, title = "E-IPI risk group", title.adj = 0)
  text(98, 15, adj = 1, cex = 0.85,
       sprintf("Calibration slope %.2f (%.2f-%.2f)\nObserved/expected events %.2f",
               sl[1, "coef"], sl[1, "coef"] - 1.96 * sl[1, "se(coef)"], sl[1, "coef"] + 1.96 * sl[1, "se(coef)"], o / e))
}

draw <- function() {
  layout(matrix(1:2, 1, 2))
  panel("OS",  "A  Overall survival")
  panel("PFS", "B  Progression-free survival")
}
png(file.path(OUT_DIR, "Calibration_EIPI_validation.png"), width = 11, height = 5.4, units = "in", res = 600); draw(); dev.off()
pdf(file.path(OUT_DIR, "Calibration_EIPI_validation.pdf"), width = 11, height = 5.4); draw(); dev.off()
tiff(file.path(OUT_DIR, "Calibration_EIPI_validation.tiff"), width = 11, height = 5.4, units = "in", res = 600,
     compression = "lzw"); draw(); dev.off()
cat("Saved output/Calibration_EIPI_validation.png, .pdf and .tiff\n")
