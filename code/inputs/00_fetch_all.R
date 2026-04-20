# ------------------------------------------------------------------------------
# code/inputs/00_fetch_all.R
#
# Driver that sources every fetch / build script in code/inputs/ in the right
# order. Use this to kick off parallel downloads in the background, e.g.
#
#     nohup Rscript code/inputs/00_fetch_all.R > fetch.log 2>&1 &
#
# Design:
#   * Each child script is wrapped in tryCatch so a single failed source (e.g.
#     ERW Zenodo DOI not resolving yet) does not block the rest.
#   * A summary table is printed at the end showing which datasets succeeded,
#     which failed, and where the staged outputs live.
#   * Ordering: dependencies first (template -> prices + constraints + soil-
#     climate-driven suitability -> potentials that may reference those).
#
# This script is NOT a substitute for the step1/step2/step3 pipeline -- it
# only populates data/inputs/ and data/staging/. The mitigation / production
# / adaptation math still lives under code/activity2/ and code/activity3/.
# ------------------------------------------------------------------------------

# ---- Preamble -----------------------------------------------------------
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                      error = function(e) "")
inputs_dir <- if (nzchar(this_file)) dirname(this_file) else "code/inputs"
common_dir <- normalizePath(file.path(inputs_dir, "..", "common"), mustWork = FALSE)
source(file.path(common_dir, "config.R"))
source(file.path(common_dir, "logging.R"))

cca_banner("CCA pipeline orchestrator (inputs -> Activity 2 -> Activity 3)")

# Order: dependencies first. Template must exist for any clipping step.
activity_dir <- normalizePath(file.path(inputs_dir, "..", "activity2"), mustWork = FALSE)
activity3_dir <- normalizePath(file.path(inputs_dir, "..", "activity3"), mustWork = FALSE)

steps <- list(
  # --- Inputs phase (skips silently if outputs already exist) -----------
  # Order matters: naturebase is required by Step 1 (mitigation); mapspam is
  # required by Step 2 + residue proxy + AWD T1 fallback; era is required
  # by Step 2 and build_buffering_coefficients. We fetch them up-front so
  # everything downstream is idempotent.
  list(tag = "naturebase",  path = file.path(inputs_dir, "fetch_naturebase.R")),
  list(tag = "era",         path = file.path(inputs_dir, "fetch_era.R")),
  list(tag = "mapspam",     path = file.path(inputs_dir, "fetch_mapspam.R")),
  list(tag = "template",    path = file.path(inputs_dir, "build_cca_template.R")),
  list(tag = "fao_prices",  path = file.path(inputs_dir, "fetch_fao_prices.R")),
  list(tag = "constraints", path = file.path(inputs_dir, "fetch_constraint_masks.R")),
  list(tag = "erw_pot",     path = file.path(inputs_dir, "fetch_erw_potential.R")),
  list(tag = "basalt",      path = file.path(inputs_dir, "fetch_basalt.R")),
  list(tag = "awd",         path = file.path(inputs_dir, "fetch_cao_awd.R")),
  list(tag = "biochar_suit",path = file.path(inputs_dir, "build_biochar_suitability.R")),
  list(tag = "erw_suit",    path = file.path(inputs_dir, "build_erw_suitability.R")),
  list(tag = "residue",     path = file.path(inputs_dir, "build_residue_proxy.R")),
  # Biochar potential is SYNTHESISED (no download). Must run after both
  # `residue` (build_residue_proxy.R) and `biochar_suit` (build_biochar_suitability.R)
  # because it multiplies their outputs together. Replaces the deprecated
  # fetch_biochar_potential.R, which downloaded Sanderman 2020's MG bundle
  # (cover crops + residue + biochar combined) and mislabelled it biochar-only.
  # See docs/methodology_biochar.md.
  list(tag = "biochar_pot", path = file.path(inputs_dir, "build_biochar_potential.R")),
  # Hazard frequency rasters (drought / heat / waterlogging x moderate /
  # severe / extreme x maiz / rice / sorg) are SYNTHESISED from WorldClim
  # 2.1 monthly climatology + Ecocrop crop optima. Must run after
  # `template` and after `biochar_suit` (which is what triggers the
  # geodata::worldclim_global download that this script reads from disk).
  # See docs/methodology_hazards.md for the derivation.
  list(tag = "hazards",     path = file.path(inputs_dir, "build_hazard_frequencies.R")),
  list(tag = "buffering",   path = file.path(inputs_dir, "build_buffering_coefficients.R")),
  # --- Activity 2 + 3 phase --------------------------------------------
  list(tag = "step1",       path = file.path(activity_dir, "step1_mitigation.R")),
  list(tag = "step2",       path = file.path(activity_dir, "step2_production.R")),
  list(tag = "step3",       path = file.path(activity_dir, "step3_adaptation.R")),
  list(tag = "step4",       path = file.path(activity_dir, "step4_integration.R")),
  list(tag = "step5",       path = file.path(activity_dir, "step5_uncertainty.R")),
  list(tag = "neb",         path = file.path(activity3_dir, "neb.R"))
)

results <- data.frame(
  tag     = character(),
  script  = character(),
  status  = character(),
  seconds = numeric(),
  stringsAsFactors = FALSE
)

# We invoke each child script as its own Rscript subprocess. Child scripts
# intentionally call `quit(status = 0)` when their outputs are already
# staged; running them inside a single R process via source() would end the
# whole orchestrator at the first short-circuit. Subprocess isolation also
# means that a segfault in one heavy step (e.g. a terra + large raster issue)
# does not take the whole run down.
rscript <- file.path(R.home("bin"), "Rscript")
for (step in steps) {
  if (!file.exists(step$path)) {
    cca_log(sprintf("Missing script: %s -- skipping", step$path), level = "WARN")
    results <- rbind(results, data.frame(
      tag = step$tag, script = basename(step$path),
      status = "MISSING", seconds = 0, stringsAsFactors = FALSE
    ))
    next
  }
  cca_banner(sprintf("Run: %s", step$tag))
  t0 <- Sys.time()
  exit <- tryCatch(
    system2(rscript, args = shQuote(step$path), stdout = "", stderr = ""),
    error = function(e) { cca_log(conditionMessage(e), level = "ERROR"); 1L }
  )
  t1 <- Sys.time()
  results <- rbind(results, data.frame(
    tag = step$tag, script = basename(step$path),
    status = if (identical(exit, 0L)) "OK" else "FAIL",
    seconds = as.numeric(difftime(t1, t0, units = "secs")),
    stringsAsFactors = FALSE
  ))
}

cca_banner("Summary")
print(results, row.names = FALSE)

# Exit code reflects whether any step failed. Useful for CI hooks.
if (any(results$status == "FAIL")) {
  cca_log(sprintf("%d step(s) failed -- see log above", sum(results$status == "FAIL")),
          level = "ERROR")
  quit(status = 1)
}
cca_log("All steps OK.")
