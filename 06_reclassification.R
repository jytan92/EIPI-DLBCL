# =============================================================================
# 06_reclassification.R : reclassification by E-IPI relative to IPI and NCCN-IPI
#   (a) Reclassification counts                     -> S5 (development), S10 (validation, by country)
#       Downward : comparator high-int/high -> E-IPI low/low-int
#       Upward i : comparator high-int      -> E-IPI high
#       Upward ii: comparator low/low-int   -> E-IPI high-int/high
#   (b) 4-year OS/PFS of reclassified vs non-reclassified patients (Results text)
#   (c) (b) restricted to R-CHOP-treated patients                    -> S6
# Primary cohorts. Outputs: output/S5_S10_reclassification.csv,
#   output/reclassification_outcomes.csv, output/S6_rchop_reclassification.csv
# =============================================================================
source("00_setup.R")

COMPARATORS <- c("IPI" = "ipi_grp", "NCCN-IPI" = "nccn_grp")
n_pct <- function(k, n) sprintf("%d/%d (%.1f%%)", k, n, if (n > 0) 100 * k / n else 0)

# ---- (a) counts -------------------------------------------------------------------
reclass_counts <- function(x, label) {
  row <- data.frame(Cohort = sprintf("%s (n=%d)", label, nrow(x)), check.names = FALSE)
  for (cp in names(COMPARATORS)) {
    g <- x[[COMPARATORS[cp]]]
    row[[paste("Down vs", cp)]] <- n_pct(sum(g >= 3 & x$eipi_grp <= 2), sum(g >= 3))
  }
  for (cp in names(COMPARATORS)) {
    g <- x[[COMPARATORS[cp]]]
    row[[paste("Up HI->H vs", cp)]] <- n_pct(sum(g == 3 & x$eipi_grp == 4), sum(g == 3))
  }
  for (cp in names(COMPARATORS)) {
    g <- x[[COMPARATORS[cp]]]
    row[[paste("Up L/LI->HI/H vs", cp)]] <- n_pct(sum(g <= 2 & x$eipi_grp >= 3), sum(g <= 2))
  }
  row
}
dev <- COHORTS$primary$dev
val <- COHORTS$primary$val
counts <- rbind(
  reclass_counts(dev, "Development"),
  reclass_counts(val, "Validation"),
  do.call(rbind, lapply(c("Hong Kong", "Malaysia", "Korea", "China"),
                        function(ctry) reclass_counts(val[val$Country == ctry, ], ctry))))

# ---- (b) outcomes of reclassified patients ---------------------------------------
# moved = reclassified by E-IPI; compared with patients in the same comparator band who were not
compare_moved <- function(g, moved, ep) {
  a <- km_at(g[moved, ], ep)["est"]; b <- km_at(g[!moved, ], ep)["est"]
  lr <- survdiff(surv_obj(g, ep) ~ moved)
  c(sprintf("%.1f%% vs %.1f%%", a, b), fmt_p(pchisq(lr$chisq, 1, lower.tail = FALSE)))
}
reclass_outcomes <- function(x, label, directions = c("Downward", "Upward")) {
  do.call(rbind, lapply(names(COMPARATORS), function(cp) {
    grp <- x[[COMPARATORS[cp]]]
    do.call(rbind, lapply(directions, function(dir) {
      if (dir == "Downward") {
        g <- x[grp >= 3, ]; moved <- g$eipi_grp <= 2
        desc <- "comparator HI/H: reclassified to E-IPI L/LI vs remained E-IPI HI/H"
      } else {
        g <- x[grp == 3, ];  moved <- g$eipi_grp == 4
        desc <- "comparator HI: upgraded to E-IPI H vs remained below E-IPI H"
      }
      os <- compare_moved(g, moved, "OS"); pfs <- compare_moved(g, moved, "PFS")
      data.frame(Cohort = label, Comparator = cp, Direction = dir, Comparison = desc,
                 `Reclassified/total` = sprintf("%d/%d", sum(moved), nrow(g)),
                 `4-y OS, reclassified vs not` = os[1], `OS log-rank P` = os[2],
                 `4-y PFS, reclassified vs not` = pfs[1], `PFS log-rank P` = pfs[2], check.names = FALSE)
    }))
  }))
}
outcomes <- rbind(reclass_outcomes(dev, sprintf("Development (n=%d)", nrow(dev))),
                  reclass_outcomes(val, sprintf("Validation (n=%d)", nrow(val))))

# ---- (c) reclassification outcomes, R-CHOP only (S6) --------------------------------
rchop <- function(x) x[x$regimen == 1, ]
S6 <- rbind(reclass_outcomes(rchop(dev), sprintf("Development, R-CHOP (n=%d)", nrow(rchop(dev)))),
            reclass_outcomes(rchop(val), sprintf("Validation, R-CHOP (n=%d)", nrow(rchop(val)))))

# R-EPOCH-treated patients among those reclassified downward (all patients)
repoch <- do.call(rbind, lapply(list(Development = dev, Validation = val), function(x) {
  do.call(rbind, lapply(names(COMPARATORS), function(cp) {
    down <- x[x[[COMPARATORS[cp]]] >= 3 & x$eipi_grp <= 2, ]
    data.frame(Comparator = cp, `Downward reclassified` = nrow(down),
               `of whom R-EPOCH` = sum(down$regimen == 2), check.names = FALSE)
  }))
}))
repoch <- cbind(Cohort = rep(c("Development", "Validation"), each = 2), repoch)

cat("\n== Reclassification counts (S5 = development; S10 = validation and by country) ==\n")
print(counts, row.names = FALSE)
cat("\n== Survival of reclassified vs non-reclassified patients (Results text) ==\n")
print(outcomes[, -4], row.names = FALSE)
cat("\n== S6: reclassification outcomes, R-CHOP-treated patients only ==\n")
print(S6[, -4], row.names = FALSE)
cat("\n== R-EPOCH-treated patients among those reclassified downward ==\n")
print(repoch, row.names = FALSE)

save_table(counts,   "S5_S10_reclassification")
save_table(outcomes, "reclassification_outcomes")
save_table(S6,       "S6_rchop_reclassification")
