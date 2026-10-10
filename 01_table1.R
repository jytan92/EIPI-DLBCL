# =============================================================================
# 01_table1.R : Table 1 — P values for the development vs validation comparison
# Primary cohorts (all patients). Output: output/table1_pvalues.csv
#   Continuous variables : Mann-Whitney U test
#   Categorical variables: chi-square test (Yates continuity correction for 2x2 tables; Fisher's exact test if any expected count < 5)
# Frequencies, medians and IQRs are not produced here.
# =============================================================================
source("00_setup.R")

d <- COHORTS$primary$dev
v <- COHORTS$primary$val

# continuous variable: Mann-Whitney U test
p_continuous <- function(label, var) {
  data.frame(Characteristic = label, Test = "Mann-Whitney U",
             P = fmt_p(wilcox.test(d[[var]], v[[var]], exact = FALSE)$p.value))
}

# categorical variable (binary or several levels); missing values excluded
p_categorical <- function(label, x_dev, x_val) {
  x_dev <- x_dev[!is.na(x_dev)]; x_val <- x_val[!is.na(x_val)]
  tab <- table(c(x_dev, x_val), rep(c("dev", "val"), c(length(x_dev), length(x_val))))
  ct  <- suppressWarnings(chisq.test(tab))
  if (any(ct$expected < 5)) {
    p <- if (nrow(tab) == 2) fisher.test(tab)$p.value
         else fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
    test <- "Fisher's exact"
  } else {
    p <- ct$p.value; test <- "Chi-square"
  }
  data.frame(Characteristic = label, Test = test, P = fmt_p(p))
}

# cell of origin: the column name differs between files, so find it by its prefix
coo <- function(x) suppressWarnings(as.numeric(x[[grep("^COO", names(x))[1]]]))

table1_p <- rbind(
  p_continuous("Age, years", "Age"),
  p_categorical("Age >60 years", d$age60, v$age60),
  p_categorical("Male sex", d$male, v$male),
  p_categorical("Ann Arbor stage III-IV", d$stage34, v$stage34),
  p_categorical("ECOG performance status >=2", d$ecog2, v$ecog2),
  p_categorical("LDH > upper limit of normal", d$ldh_high, v$ldh_high),
  p_categorical(">1 extranodal site", d$en2, v$en2),
  p_categorical("Non-GCB (among available)", coo(d), coo(v)),
  p_continuous("Platelet count", "plt"),
  p_continuous("LDH ratio", "ldh_ratio"),
  p_continuous("Creatinine", "creat"),
  p_continuous("Modified EASIX", "measix"),
  p_categorical("mEASIX tier (3 levels)", d$tier, v$tier),
  p_categorical("IPI risk group (4 levels)", d$ipi_grp, v$ipi_grp),
  p_categorical("NCCN-IPI risk group (4 levels)", d$nccn_grp, v$nccn_grp),
  p_categorical("E-IPI risk group (4 levels)", d$eipi_grp, v$eipi_grp),
  p_categorical("Primary therapy (R-CHOP vs R-EPOCH)", d$regimen, v$regimen),
  p_categorical("POD12 (evaluable patients)", d$pod12, v$pod12),
  p_categorical("Deaths", d$os_event, v$os_event),
  p_categorical("Progression or death", d$pfs_event, v$pfs_event),
  p_categorical("Death within 6 months of diagnosis",
                d$os_event == 1 & d$os_time < 6, v$os_event == 1 & v$os_time < 6)
)

print(table1_p, row.names = FALSE)
save_table(table1_p, "table1_pvalues")
