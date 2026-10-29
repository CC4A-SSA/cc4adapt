# Run the adaptation-factor formulation exploration end to end.
#   Rscript code/explore/af_formulations/run_all.R          # all steps
#   Rscript code/explore/af_formulations/run_all.R 3 6      # steps 3..6 only
# Env: AF_OUT_DIR (default data/outputs/af_formulations), AF_RECORD_DIR
# (optional copy, e.g. the project's _claude/outputs/adaptation_factor),
# AF_FS_DIR (Dixon farming-system rasters). See README.md.

if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
args <- as.integer(commandArgs(trailingOnly = TRUE))
from <- if (length(args) >= 1) args[1] else 1L
to   <- if (length(args) >= 2) args[2] else 6L
steps <- c("01_fetch_inputs.R", "02_check_inputs.R", "03_params.R", "04_units.R", "05_compute.R", "06_diagnostics.R")

message("Tests")
testthat::test_dir("code/explore/af_formulations/tests/testthat", reporter = "summary", stop_on_failure = TRUE)
for (i in from:to) {
  message("\n==== ", steps[i])
  status <- system2("Rscript", file.path("code/explore/af_formulations", steps[i]))
  if (status != 0) stop(steps[i], " failed with status ", status)
}
message("\nDone. Outputs in ", Sys.getenv("AF_OUT_DIR", "data/outputs/af_formulations"))
