# Run the adaptation-factor formulation exploration end to end.
#   Rscript code/explore/af_formulations/run_all.R            # core: tests, then steps 01-06 and 12
#   Rscript code/explore/af_formulations/run_all.R 3 6        # core steps 3..6 only
#   Rscript code/explore/af_formulations/run_all.R --pilot    # ERA season pilot: 03b, 03c, (07 by hand), 08, 09, 11, 13, 12
# Env: AF_OUT_DIR (default data/outputs/af_formulations), AF_RECORD_DIR
# (optional copy), AF_FS_DIR (Dixon farming-system rasters), GCP_PROJECT_ID
# (step 07 and 10, run by hand: they need Earth Engine and take hours).
# See README.md.

if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
args <- commandArgs(trailingOnly = TRUE)
pilot <- "--pilot" %in% args
nums <- suppressWarnings(as.integer(setdiff(args, "--pilot")))
from <- if (length(nums) >= 1) nums[1] else 1L
to   <- if (length(nums) >= 2) nums[2] else 6L

core  <- c("01_fetch_inputs.R", "02_check_inputs.R", "03_params.R", "04_units.R", "05_compute.R", "06_diagnostics.R")
pilot_steps <- c("03b_era_site_years.R", "03c_era_site_windows.R", "08_era_split.R", "09_obs_vs_atlas.R",
                 "11_gcm_summary.R", "13_kenya_dotplot.R")
run <- function(f) {
  message("\n==== ", f)
  status <- system2("Rscript", file.path("code/explore/af_formulations", f))
  if (status != 0) stop(f, " failed with status ", status)
}

message("Tests")
testthat::test_dir("code/explore/af_formulations/tests/testthat", reporter = "summary", stop_on_failure = TRUE)
if (pilot) {
  haz <- "data/staging/af_formulations/era_site_season_hazards.csv"
  for (f in pilot_steps[1:2]) run(f)
  if (!file.exists(haz)) {
    stop("Run step 07 by hand first (needs Earth Engine):\n  GCP_PROJECT_ID=<project> PYTHONPATH=../climate-toolkit ",
         "python3 code/explore/af_formulations/07_era_seasons.py --workers 2\nthen re-run with --pilot.")
  }
  for (f in pilot_steps[3:6]) {
    if (f == "11_gcm_summary.R" && !length(list.files("data/staging/af_formulations", "^gcm_vs_obs_sites.*csv$"))) {
      message("skipping 11 and 13: run step 10 by hand first (Earth Engine)"); break
    }
    run(f)
  }
} else {
  for (i in from:to) run(core[i])
}
run("12_export_page_data.R")
message("\nDone. Outputs in ", Sys.getenv("AF_OUT_DIR", "data/outputs/af_formulations"), "; page data in docs/af_formulations/data.js")
