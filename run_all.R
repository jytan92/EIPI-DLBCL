# =============================================================================
# run_all.R : runs every analysis script in order. Set the working directory to this
# folder (containing data/development.xlsx and data/validation.xlsx) and run:
#   source("run_all.R")
# Total run time: about 5 minutes. All output tables are written to output/.
# =============================================================================
scripts <- c("01_table1.R", "02_measix_prognostic.R", "03_score_construction.R",
             "04_eipi_survival.R", "05_discrimination.R", "06_reclassification.R",
             "07_calibration.R", "08_secondary.R", "09_pod12.R")
for (s in scripts) {
  cat("\n\n############", s, "############\n")
  source(s, echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path("output", "sessionInfo.txt"))
cat("\nAll analyses complete. Tables are in the output/ folder.\n")
