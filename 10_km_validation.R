# =============================================================================
# 10_km_validation.R : Kaplan-Meier curves of OS and PFS by E-IPI risk group,
# external validation cohort (n=533), with censoring marks and numbers at risk.
# Base R graphics only. Outputs: output/KM_EIPI_validation.png (600 dpi), .pdf and .tiff
# =============================================================================
source("00_setup.R")

x <- COHORTS$primary$val
x$grp <- factor(x$eipi_grp, levels = 1:4, labels = GROUP_LABELS)
COLS   <- c("#1F4E9A", "#2E9E3E", "#E07B00", "#C62828")   # low, low-int, high-int, high
YEARS  <- 0:6                                              # x-axis, years

km_panel <- function(ep, ylab) {
  tt <- (if (ep == "OS") x$os_time else x$pfs_time) / 12   # months -> years
  ee <-  if (ep == "OS") x$os_event else x$pfs_event
  fit <- survfit(Surv(tt, ee) ~ grp, data = x)
  lr  <- survdiff(Surv(tt, ee) ~ grp, data = x)
  p   <- pchisq(lr$chisq, 3, lower.tail = FALSE)

  par(mar = c(10.5, 7.5, 1.5, 1.5), mgp = c(2.6, 0.7, 0), las = 1, cex.axis = 0.95)
  plot(fit, col = COLS, lwd = 2.2, mark.time = TRUE, mark = "|", cex = 0.7,
       xlim = c(0, max(YEARS)), ylim = c(0, 1), xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
       xlab = "Time (years)", ylab = ylab, bty = "l")
  axis(1, at = YEARS)
  axis(2, at = seq(0, 1, 0.2), labels = seq(0, 100, 20))
  abline(h = seq(0.2, 1, 0.2), col = "grey90", lty = 3)
  text(max(YEARS) * 0.97, 0.06, adj = 1, cex = 1,
       labels = if (p < 0.001) "Log-rank P < 0.001" else sprintf("Log-rank P = %.3f", p))

  # numbers at risk
  s <- summary(fit, times = YEARS, extend = TRUE)
  risk <- matrix(s$n.risk, nrow = 4, byrow = TRUE)
  short <- c("Low", "Low-int", "High-int", "High")
  mtext("Number at risk", side = 1, line = 4.2, adj = 0, at = -1.45, cex = 0.8, font = 2)
  for (g in 1:4) {
    mtext(short[g], side = 1, line = 4.2 + 1.2 * g, adj = 0, at = -1.45, cex = 0.72, col = COLS[g], font = 2)
    mtext(risk[g, ], side = 1, line = 4.2 + 1.2 * g, at = YEARS, cex = 0.72)
  }
  invisible(fit)
}

draw <- function() {
  layout(matrix(1:2, 1, 2))
  par(oma = c(0, 0, 2.2, 0))
  km_panel("OS",  "Overall survival (%)")
  mtext("A", side = 3, line = 0.2, adj = -0.25, font = 2, cex = 1.3)
  km_panel("PFS", "Progression-free survival (%)")
  mtext("B", side = 3, line = 0.2, adj = -0.25, font = 2, cex = 1.3)
  # shared legend across the top
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("top", legend = paste("E-IPI", c("low (0-1)", "low-intermediate (2-3)", "high-intermediate (4)", "high (5-7)")),
         col = COLS, lwd = 2.5, horiz = TRUE, bty = "n", cex = 0.95, inset = 0.005)
}

png(file.path(OUT_DIR, "KM_EIPI_validation.png"), width = 12, height = 5.6, units = "in", res = 600); draw(); dev.off()
pdf(file.path(OUT_DIR, "KM_EIPI_validation.pdf"), width = 12, height = 5.6); draw(); dev.off()
tiff(file.path(OUT_DIR, "KM_EIPI_validation.tiff"), width = 12, height = 5.6, units = "in", res = 600,
     compression = "lzw"); draw(); dev.off()
cat("Saved output/KM_EIPI_validation.png, .pdf and .tiff\n")
