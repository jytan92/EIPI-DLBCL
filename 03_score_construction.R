# =============================================================================
# 03_score_construction.R : derivation and construction of the E-IPI (development cohort)
#   (a) Lower mEASIX cut-off: maximally selected log-rank statistic (OS), full cohort and
#       patients alive with >=6 months follow-up; discrimination with lower cut-off 0.60 vs 0.70
#   (b) Upper mEASIX cut-off: grid search, C-index of IPI + tier (lower cut-off fixed at 0.70)
#   (c) Four-group partition of the E-IPI score (each group >=5% of the cohort)
#   (d) Multivariable Cox model: five IPI factors + mEASIX tier         (Supplementary Table S3A)
#   (e) C-index of alternative point allocations                        (Supplementary Table S3B)
#   (f) mEASIX tier added to IPI vs NCCN-IPI (bootstrap comparison)      (Supplementary Table S3C)
# Validation-cohort columns apply development-derived scores/weights without re-estimation.
# Outputs: output/derivation.csv, output/S3A_multivariable.csv, output/S3B_point_allocation.csv,
#          output/S3C_ipi_vs_nccn_base.csv
# =============================================================================
source("00_setup.R")

d <- COHORTS$primary$dev
v <- COHORTS$primary$val
IPI_FACTORS <- c(age60 = "Age >60", stage34 = "Stage III-IV", ecog2 = "ECOG >=2",
                 ldh_high = "LDH > ULN", en2 = ">1 extranodal site")
for (nm in c("d", "v")) {
  x <- get(nm)
  x$t1 <- as.integer(x$tier == 1); x$t2 <- as.integer(x$tier == 2)
  assign(nm, x)
}
time_event <- function(ep) if (ep == "OS") c("os_time", "os_event") else c("pfs_time", "pfs_event")
cidx <- function(x, ep, score) { te <- time_event(ep); harrell_c(x[[te[1]]], x[[te[2]]], score) }
logrank_chisq <- function(x, ep, group) survdiff(surv_obj(x, ep) ~ group)$chisq

# ---- (a) lower cut-off ---------------------------------------------------------
# Maximally selected log-rank statistic (Lausen & Schumacher; maxstat package), OS:
# candidate cut-points = every observed mEASIX value between the 10th and 90th percentiles,
# split mEASIX <= c vs > c. Run in the full development cohort (primary) and in patients
# alive with at least 6 months of follow-up (secondary population).
lower_cutoff <- function(x) {
  ms <- maxstat.test(Surv(os_time, os_event) ~ measix, data = x, smethod = "LogRank",
                     pmethod = "Lau94", minprop = 0.10, maxprop = 0.90)
  list(cut = unname(ms$estimate), M = max(ms$stats), p = ms$p.value, n = nrow(x))
}
lc_full <- lower_cutoff(d)
lc_6mo  <- lower_cutoff(COHORTS$secondary$dev)

# discrimination with alternative lower cut-offs (upper fixed at 1.50; mEASIX rounded to 2 dp)
lower_compare <- do.call(rbind, lapply(c(0.60, 0.70), function(lo) {
  do.call(rbind, lapply(list(Development = d, Validation = v), function(x) {
    r  <- round(x$measix, 2)
    tr <- ifelse(r < lo, 0, ifelse(r <= 1.50, 1, 2))
    sc <- x$ipi + tr; gp <- eipi_group(sc)
    data.frame(`Lower cut-off` = sprintf("%.2f", lo), n = nrow(x),
               `Tier n (low/int/high)` = paste(tabulate(tr + 1, 3), collapse = "/"),
               `C OS, IPI + tier` = sprintf("%.3f", cidx(x, "OS", sc)),
               `C PFS, IPI + tier` = sprintf("%.3f", cidx(x, "PFS", sc)),
               `C OS, E-IPI groups` = sprintf("%.3f", cidx(x, "OS", gp)),
               `C PFS, E-IPI groups` = sprintf("%.3f", cidx(x, "PFS", gp)), check.names = FALSE)
  }))
}))
lower_compare <- cbind(Cohort = rep(c("Development", "Validation"), 2), lower_compare)

# ---- (b) upper cut-off ---------------------------------------------------------
grid <- round(seq(0.80, 5.00, by = 0.05), 2)
grid_c <- sapply(grid, function(u) {
  tier_u <- ifelse(d$measix < 0.70, 0, ifelse(d$measix <= u, 1, 2))
  cidx(d, "OS", d$ipi + tier_u)
})
grid_rank <- rank(-grid_c, ties.method = "min")

# ---- (c) four-group partitions of the E-IPI score (0-7) ------------------------
cuts <- combn(1:7, 3)                       # group boundaries a < b < c
parts <- do.call(rbind, lapply(seq_len(ncol(cuts)), function(k) {
  b <- cuts[, k]
  g <- findInterval(d$eipi, b)              # 0..3
  n <- tabulate(g + 1, nbins = 4)
  if (min(n) < 0.05 * nrow(d)) return(NULL)
  data.frame(Grouping = sprintf("%d-%d / %d-%d / %d-%d / %d-7", 0, b[1] - 1, b[1], b[2] - 1, b[2], b[3] - 1, b[3]),
             chisq = logrank_chisq(d, "OS", g) + logrank_chisq(d, "PFS", g),
             n = paste(n, collapse = "/"))
}))
parts <- parts[order(-parts$chisq), ]
parts$Rank <- seq_len(nrow(parts))
published <- "0-1 / 2-3 / 4-4 / 5-7"

derivation <- rbind(
  data.frame(Step = sprintf("Lower cut-off, patients alive with >=6 months follow-up (n=%d)", lc_6mo$n),
             Result = sprintf("%.3f", lc_6mo$cut),
             Statistic = sprintf("M = %.3f; Lausen-Schumacher P %s", lc_6mo$M, fmt_p(lc_6mo$p))),
  data.frame(Step = sprintf("Lower cut-off, full development cohort (n=%d)", lc_full$n),
             Result = sprintf("%.3f", lc_full$cut),
             Statistic = sprintf("M = %.3f; Lausen-Schumacher P %s", lc_full$M, fmt_p(lc_full$p))),
  data.frame(Step = "Upper cut-off (lower fixed 0.70; grid 0.80-5.00 by 0.05; OS C-index of IPI + tier)",
             Result = sprintf("1.50 ranked %d of %d (C %.4f); best %.2f (C %.4f)",
                              grid_rank[grid == 1.50], length(grid), grid_c[grid == 1.50],
                              grid[which.max(grid_c)], max(grid_c)),
             Statistic = ""),
  data.frame(Step = paste0("Partition rank ", parts$Rank[1:5], " (OS + PFS log-rank chi-square; each group >=5%)"),
             Result = parts$Grouping[1:5],
             Statistic = sprintf("chi-square %.1f; n = %s", parts$chisq[1:5], parts$n[1:5])),
  data.frame(Step = "Published grouping 0-1 / 2-3 / 4 / 5-7",
             Result = sprintf("rank %d of %d", parts$Rank[parts$Grouping == published], nrow(parts)),
             Statistic = "")
)

# ---- (d) multivariable model (S3A) and spacing/additivity tests -----------------
mv_terms <- c(names(IPI_FACTORS), "t1", "t2")
fit_cox <- function(x, ep, terms) coxph(as.formula(paste("surv_obj(x, ep) ~", paste(terms, collapse = " + "))),
                                         data = x, ties = "efron")
lr_test <- function(big, small) { s <- 2 * (big$loglik[2] - small$loglik[2]); c(chisq = s, p = pchisq(s, 1, lower.tail = FALSE)) }

S3A <- NULL; tests <- NULL; weights <- list()
for (ep in c("OS", "PFS")) {
  m <- fit_cox(d, ep, mv_terms); s <- summary(m); b <- coef(m)
  weights[[ep]] <- b
  ref <- mean(b[names(IPI_FACTORS)])
  S3A <- rbind(S3A, data.frame(Endpoint = ep,
    Variable = c(IPI_FACTORS, "mEASIX intermediate (0.70-1.50)", "mEASIX high (>1.50)"),
    beta = sprintf("%.3f", b),
    `HR (95% CI)` = fmt_hr(s$conf.int[, 1], s$conf.int[, 3], s$conf.int[, 4]),
    P = fmt_p(s$coefficients[, "Pr(>|z|)"]),
    `Equivalent IPI points` = sprintf("%.2f", b / ref), check.names = FALSE))
  lin  <- lr_test(m, fit_cox(d, ep, c(names(IPI_FACTORS), "tier")))
  d$tier_x_ipi <- d$tier * d$ipi
  int  <- lr_test(fit_cox(d, ep, c("ipi", "tier", "tier_x_ipi")), fit_cox(d, ep, c("ipi", "tier")))
  tests <- rbind(tests, data.frame(Endpoint = ep, `Mean IPI-factor beta` = sprintf("%.3f", ref),
    `beta high / beta intermediate` = sprintf("%.2f", b["t2"] / b["t1"]),
    `Linear vs categorical tier (LR)` = sprintf("chi-square %.2f, P=%s", lin["chisq"], fmt_p(lin["p"])),
    `Tier x IPI interaction (LR)` = sprintf("chi-square %.2f, P=%s", int["chisq"], fmt_p(int["p"])),
    check.names = FALSE))
}

# ---- (e) alternative point allocations (S3B) ----------------------------------
allocations <- list(
  "IPI + 0/1/2 (E-IPI)" = function(x, ep) x$ipi + x$tier,
  "Coefficient-weighted (development beta)" = function(x, ep) as.vector(as.matrix(x[, mv_terms]) %*% weights[[ep]]),
  "IPI + 0/1/3" = function(x, ep) x$ipi + ifelse(x$tier == 2, 3, x$tier),
  "IPI + 0/2/4" = function(x, ep) x$ipi + 2 * x$tier,
  "IPI + 0/1/1" = function(x, ep) x$ipi + (x$tier > 0),
  "IPI alone"   = function(x, ep) x$ipi)
S3B <- do.call(rbind, lapply(names(allocations), function(a) {
  f <- allocations[[a]]
  data.frame(Allocation = a,
             `OS dev` = sprintf("%.3f", cidx(d, "OS", f(d, "OS"))),  `OS val` = sprintf("%.3f", cidx(v, "OS", f(v, "OS"))),
             `PFS dev` = sprintf("%.3f", cidx(d, "PFS", f(d, "PFS"))), `PFS val` = sprintf("%.3f", cidx(v, "PFS", f(v, "PFS"))),
             check.names = FALSE)
}))

# ---- (f) IPI-based vs NCCN-IPI-based score (S3C); paired bootstrap ---------------
S3C <- NULL
for (ck in c("Development", "Validation")) {
  x <- if (ck == "Development") d else v
  a <- x$ipi + x$tier; n <- x$nccn + x$tier
  for (ep in c("OS", "PFS")) {
    te <- time_event(ep); tt <- x[[te[1]]]; ee <- x[[te[2]]]
    diffs <- replicate(N_BOOT, {
      i <- sample.int(nrow(x), replace = TRUE)
      harrell_c(tt[i], ee[i], a[i]) - harrell_c(tt[i], ee[i], n[i])
    })
    ca <- harrell_c(tt, ee, a); cn <- harrell_c(tt, ee, n)
    ci <- quantile(diffs, c(0.025, 0.975))
    S3C <- rbind(S3C, data.frame(Cohort = sprintf("%s (n=%d)", ck, nrow(x)), Endpoint = ep,
      `IPI + tier` = sprintf("%.3f", ca), `NCCN-IPI + tier` = sprintf("%.3f", cn),
      `Difference (95% CI)` = sprintf("%+.3f (%+.3f, %+.3f)", ca - cn, ci[1], ci[2]),
      P = fmt_p(2 * min(mean(diffs <= 0), mean(diffs >= 0))),
      `IPI alone` = sprintf("%.3f", harrell_c(tt, ee, x$ipi)),
      `NCCN-IPI alone` = sprintf("%.3f", harrell_c(tt, ee, x$nccn)), check.names = FALSE))
  }
}

cat("\n== Derivation of cut-offs and risk groups ==\n"); print(derivation, row.names = FALSE)
cat("\n== Lower cut-off 0.60 vs 0.70 (upper 1.50) ==\n"); print(lower_compare, row.names = FALSE)
cat("\n== S3A: multivariable Cox model ==\n");         print(S3A, row.names = FALSE)
cat("\n== Spacing and additivity tests ==\n");           print(tests, row.names = FALSE)
cat("\n== S3B: alternative point allocations ==\n");     print(S3B, row.names = FALSE)
cat("\n== S3C: IPI vs NCCN-IPI base ==\n");              print(S3C, row.names = FALSE)

save_table(derivation, "derivation"); save_table(lower_compare, "lower_cutoff_comparison")
save_table(rbind(S3A), "S3A_multivariable"); save_table(tests, "S3A_tests")
save_table(S3B, "S3B_point_allocation"); save_table(S3C, "S3C_ipi_vs_nccn_base")
