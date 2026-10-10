# =============================================================================
# EASIX-IPI (E-IPI) in DLBCL — analysis code
# 00_setup.R : packages, data import, scoring rule, cohort definitions, helpers
#
# Run this first; every later script starts with source("00_setup.R").
# Patient-level data are not distributed with this repository (see README).
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)     # data import
  library(survival)   # Kaplan-Meier, Cox regression, concordance
  library(maxstat)    # maximally selected log-rank statistics
  library(pROC)       # ROC / AUC and DeLong test (POD12)
})

# ---- File locations (edit if your files are elsewhere) ----------------------
DEV_FILE <- "data/development.xlsx"
VAL_FILE <- "data/validation.xlsx"
OUT_DIR  <- "output"
dir.create(OUT_DIR, showWarnings = FALSE)

options(warnPartialMatchDollar = TRUE)  # flag accidental partial name matches
set.seed(2026)        # all bootstrap analyses
N_BOOT  <- 2000       # bootstrap resamples
HORIZON <- 48         # months; 4-year survival estimates

# ---- Import and standardise column names -------------------------------------
read_cohort <- function(path) {
  x <- as.data.frame(read_excel(path, .name_repair = "minimal"))
  x <- x[, names(x) != "" & !grepl("^Unnamed", names(x)), drop = FALSE]  # drop blank columns
  # the files hold "Modified EASIX" twice (identical copies); keep the first
  dup <- which(names(x) == "Modified EASIX")
  if (length(dup) > 1) {
    stopifnot(isTRUE(all.equal(x[[dup[1]]], x[[dup[2]]])))
    x <- x[, -dup[-1], drop = FALSE]
  }
  ren <- c(
    "OS (months)" = "os_time",  "OS outcome" = "os_event",
    "PFS (months)" = "pfs_time", "PFS outcome" = "pfs_event",
    "IPI raw value" = "ipi", "IPI raw" = "ipi",
    "IPI-categorised" = "ipi_grp",
    "NCCN raw" = "nccn", "NCCN category" = "nccn_grp",
    "E-IPI raw" = "eipi_file", "E-IPI categorised" = "eipi_grp_file",
    "Modified EASIX" = "measix",
    "3 tier classification (<0.7, 0.7-1.5, >1.5)" = "tier_file",
    "Gender code (Male = 1, Female =2)" = "sex",
    "Age>60=1" = "age60", "Stage 3 or 4=1" = "stage34", "ECOG2-4=1" = "ecog2",
    "LDHraised=1" = "ldh_high", ">1Extranodalsite=1" = "en2",
    "PLT" = "plt", "LDH ratio" = "ldh_ratio", "Creatinine (mg/dl)" = "creat",
    "Treatment (RCHOP =1, REPOCH=2)" = "regimen"
  )
  hit <- names(x) %in% names(ren)
  names(x)[hit] <- ren[names(x)[hit]]
  need <- c("os_time", "os_event", "pfs_time", "pfs_event", "ipi", "ipi_grp", "nccn",
            "nccn_grp", "eipi_file", "eipi_grp_file", "measix", "tier_file", "sex",
            "plt", "ldh_ratio", "creat", "regimen", "POD12")
  missing <- setdiff(need, names(x))
  if (length(missing)) stop("Columns not found in ", path, ": ", paste(missing, collapse = ", "))
  # numbers stored as text in Excel are converted; anything non-numeric stops the run
  num <- c("os_time", "os_event", "pfs_time", "pfs_event", "ipi", "ipi_grp", "nccn", "nccn_grp",
           "eipi_file", "eipi_grp_file", "measix", "tier_file", "plt", "ldh_ratio", "creat",
           "regimen", "age60", "stage34", "ecog2", "ldh_high", "en2")
  for (v in intersect(num, names(x))) {
    if (!is.numeric(x[[v]])) {
      conv <- suppressWarnings(as.numeric(x[[v]]))
      bad <- is.na(conv) & !is.na(x[[v]])
      if (any(bad)) stop("Non-numeric value(s) in ", v, ": ", paste(unique(x[[v]][bad]), collapse = ", "))
      message("Note: ", v, " in ", basename(path), " contained numbers stored as text; converted.")
      x[[v]] <- conv
    }
  }
  x
}

dev_all <- read_cohort(DEV_FILE); dev_all$cohort <- "Development"
val_all <- read_cohort(VAL_FILE); val_all$cohort <- "Validation"
if (!"Country" %in% names(dev_all)) dev_all$Country <- "Singapore"

# ---- E-IPI scoring rule -------------------------------------------------------
# mEASIX = LDH ratio (observed/ULN) x creatinine (mg/dL) x 100 / platelets (x10^9/L)
# mEASIX is rounded to two decimal places before tier assignment:
#   low < 0.70 (0 points); intermediate 0.70-1.50 (1 point); high > 1.50 (2 points)
# E-IPI = IPI (0-5) + mEASIX points (0-2), range 0-7
# Risk groups: low 0-1, low-intermediate 2-3, high-intermediate 4, high 5-7
measix_tier <- function(m) {
  r <- round(m, 2)
  ifelse(r < 0.70, 0L, ifelse(r <= 1.50, 1L, 2L))
}
eipi_group <- function(score) {
  as.integer(cut(score, breaks = c(-Inf, 1, 3, 4, Inf), labels = 1:4))
}
GROUP_LABELS <- c("Low", "Low-intermediate", "High-intermediate", "High")

score_cohort <- function(x) {
  measix_calc <- x[["ldh_ratio"]] * x[["creat"]] * 100 / x[["plt"]]
  stopifnot(max(abs(measix_calc - x[["measix"]])) < 1e-8)   # mEASIX formula check
  x$tier     <- measix_tier(x[["measix"]])
  x$eipi     <- x$ipi + x$tier
  x$eipi_grp <- eipi_group(x$eipi)
  # the data files should already follow the rule; stop if not
  stopifnot(all(x$tier == x$tier_file), all(x$eipi == x$eipi_file),
            all(x$eipi_grp == x$eipi_grp_file))
  x$male   <- ifelse(x$sex == 1, 1L, ifelse(x$sex == 2, 0L, NA_integer_))
  x$log2_measix <- log2(x$measix)
  x$pod12_status <- as.character(x$POD12)            # "0", "1" or "Exclude"
  x$pod12  <- NA_integer_
  x$pod12[x$pod12_status == "1"] <- 1L
  x$pod12[x$pod12_status == "0"] <- 0L
  stopifnot(all(x$pod12_status %in% c("0", "1", "Exclude")))
  x
}

dev_all <- score_cohort(dev_all)
val_all <- score_cohort(val_all)

# ---- Cohort definitions -------------------------------------------------------
# Primary analysis: all patients (early deaths included).
# Secondary analysis: patients alive with at least 6 months of follow-up (OS time >= 6 months).
induction <- function(x) x[x$os_time >= 6, ]

COHORTS <- list(
  primary   = list(dev = dev_all,            val = val_all),
  secondary = list(dev = induction(dev_all), val = induction(val_all))
)

# ---- Helpers ------------------------------------------------------------------
fmt_p  <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt_hr <- function(est, lo, hi) sprintf("%.2f (%.2f-%.2f)", est, lo, hi)

surv_obj <- function(x, endpoint = c("OS", "PFS")) {
  endpoint <- match.arg(endpoint)
  if (endpoint == "OS") Surv(x$os_time, x$os_event) else Surv(x$pfs_time, x$pfs_event)
}

# Kaplan-Meier survival (%) at HORIZON months with 95% CI (log-log)
km_at <- function(x, endpoint, t = HORIZON) {
  f <- survfit(surv_obj(x, endpoint) ~ 1, conf.type = "log-log")
  s <- summary(f, times = t, extend = TRUE)
  c(est = 100 * s$surv, lo = 100 * s$lower, hi = 100 * s$upper)
}

# Harrell's C: higher score = higher risk (reverse = TRUE)
harrell_c <- function(time, event, score) {
  concordance(Surv(time, event) ~ score, reverse = TRUE)$concordance
}

save_table <- function(df, name) {
  write.csv(df, file.path(OUT_DIR, paste0(name, ".csv")), row.names = FALSE)
  invisible(df)
}

# ---- Discrimination helpers (used by 05_discrimination.R and 08_secondary.R) ---
# Risk groups (1-4) are the scores compared; higher group = higher risk.
SCORE_COLS <- c("IPI" = "ipi_grp", "NCCN-IPI" = "nccn_grp", "E-IPI" = "eipi_grp")

# Harrell's C with a percentile bootstrap CI, paired bootstrap difference vs E-IPI
# (difference = C[E-IPI] - C[comparator]; two-sided P = 2 x the smaller tail proportion),
# and AIC of a Cox model with the score entered as an ordinal (linear) term.
# The same resampled patients are used for all three scores (paired comparison).
discrimination <- function(x, n_boot = N_BOOT) {
  do.call(rbind, lapply(c("OS", "PFS"), function(ep) {
    tt <- if (ep == "OS") x$os_time else x$pfs_time
    ee <- if (ep == "OS") x$os_event else x$pfs_event
    obs <- sapply(SCORE_COLS, function(s) harrell_c(tt, ee, x[[s]]))
    boot <- t(replicate(n_boot, {
      i <- sample.int(nrow(x), replace = TRUE)
      sapply(SCORE_COLS, function(s) harrell_c(tt[i], ee[i], x[[s]][i]))
    }))
    do.call(rbind, lapply(names(SCORE_COLS), function(nm) {
      ci  <- quantile(boot[, nm], c(0.025, 0.975))
      aic <- AIC(coxph(Surv(tt, ee) ~ x[[SCORE_COLS[nm]]], ties = "efron"))
      if (nm == "E-IPI") { dC <- "Ref"; dci <- ""; p <- "" } else {
        dd  <- boot[, "E-IPI"] - boot[, nm]
        dC  <- sprintf("%+.3f", obs["E-IPI"] - obs[nm])
        q   <- quantile(dd, c(0.025, 0.975)); dci <- sprintf("(%+.3f, %+.3f)", q[1], q[2])
        p   <- fmt_p(2 * min(mean(dd <= 0), mean(dd >= 0)))
      }
      data.frame(Endpoint = ep, Score = nm, `C-index` = sprintf("%.3f", obs[nm]),
                 `95% CI` = sprintf("%.3f-%.3f", ci[1], ci[2]), `Delta C vs E-IPI` = dC,
                 `Delta C 95% CI` = dci, P = p, AIC = sprintf("%.1f", aic),
                 Events = sum(ee), N = nrow(x), check.names = FALSE)
    }))
  }))
}

# ---- Summary printed when sourced --------------------------------------------
cat(sprintf("Development: %d patients (%d deaths); Validation: %d patients (%d deaths)\n",
            nrow(dev_all), sum(dev_all$os_event), nrow(val_all), sum(val_all$os_event)))
cat(sprintf("Secondary (alive at 6 months): development %d; validation %d\n",
            nrow(COHORTS$secondary$dev), nrow(COHORTS$secondary$val)))
