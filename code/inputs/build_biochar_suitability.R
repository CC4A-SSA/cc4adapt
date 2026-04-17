# ------------------------------------------------------------------------------
# code/inputs/build_biochar_suitability.R
#
# Builds a 0-1 biochar-suitability surface at 5 arcmin over the CCA region.
# Inputs are static soil + climatological layers (SoilGrids, CHIRPS) pulled
# via the geodata package where it covers them.
#
# Why R (not Python / CropSuiteLite full pipeline):
#   CropSuiteLite is a growing-cycle crop-suitability simulator; biochar is
#   a one-time amendment with NO growing-cycle. All we need is a static fuzzy
#   suitability surface based on four biophysical drivers. Reproducing a
#   single-script implementation here keeps pipeline dependencies tight and
#   avoids adding a Python/conda environment to the CCA container.
#
# Equivalent CropSuiteLite .inf side-car is written alongside the raster for
# teams that do want to run the full simulator.
#
# --- Drivers (project plan values) ---
#
#   * Soil pH (H2O): optimum 5.5-7.5
#       Evidence: biochar benefits are strongest on acidic / moderately-acidic
#       soils; above pH 8 the carbon-stability and pH-correction arguments
#       weaken (Jeffery et al. 2017 meta-analysis). SSA soils trend acidic.
#   * Soil organic carbon (0-30 cm): 0.5-3 percent
#       Evidence: Lehmann 2007 / Woolf 2010 -- biochar yield benefits largest
#       on low-OC soils; very high OC soils derive less benefit.
#   * Clay content (0-30 cm): 10-40 percent
#       Evidence: sandy soils leak biochar; very heavy clays reduce infiltration.
#       15-35 percent is the broad optimum from meta-analyses.
#   * Mean annual rainfall: > 400 mm
#       Evidence: below 400 mm biochar-induced water retention is unlikely to
#       translate into yield benefit (too arid for rainfed cereals); above
#       2500 mm leaching risk erodes stability. We use a broad plateau.
#
# Combination rule:
#   geometric mean of the four memberships -- chosen over arithmetic mean
#   because it zero-outs unsuitable pixels (if any driver is 0, the whole
#   product is 0), which matches the biophysical intuition that biochar
#   fails if ANY driver is out of range.
#
# --- No default thresholding after suitability ---
#
#   Downstream callers can binarise this at whatever cut-off Activity 1
#   workshops agree on (e.g. >=0.3). By default the raster stays continuous.
# ------------------------------------------------------------------------------

# ---- Preamble -----------------------------------------------------------
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                      error = function(e) "")
common_dir <- if (nzchar(this_file)) {
  normalizePath(file.path(dirname(this_file), "..", "common"), mustWork = FALSE)
} else {
  candidates <- c("code/common", "../common", "../../common")
  hit <- candidates[file.exists(file.path(candidates, "config.R"))][1]
  if (is.na(hit)) stop("Cannot locate code/common from getwd() = ", getwd())
  hit
}
source(file.path(common_dir, "config.R"))
source(file.path(common_dir, "logging.R"))
source(file.path(common_dir, "geo_utils.R"))
source(file.path(common_dir, "download.R"))

suppressPackageStartupMessages({
  library(terra)
  library(geodata)
})

cca_banner("Build biochar suitability surface (pH, OC, clay, rainfall)")

staging_dir <- file.path(CCA_PATHS$staging_root, "suitability")
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

staged_tif <- file.path(staging_dir, sprintf("biochar_suitability_%s.tif", CCA_RES_LABEL))
if (cca_skip_if_exists(staged_tif, tag = "biochar_suitability")) quit(status = 0)

# ---- Fuzzy membership helper -------------------------------------------
# Piecewise-linear trapezoid: rises from lo1 to lo2 (0 -> 1), plateau to hi1
# (stays at 1), falls to hi2 (1 -> 0). Values below lo1 or above hi2 are 0.
# This matches the CropSuiteLite .inf encoding exactly (vals/suit pairs).
fuzzy_trapezoid <- function(r, lo1, lo2, hi1, hi2, name = "fuzzy") {
  # Lightweight sanity check; catches typos in the calling config.
  stopifnot(lo1 <= lo2, lo2 <= hi1, hi1 <= hi2)
  # terra arithmetic on rasters -- no cell-by-cell loop needed.
  rising  <- terra::clamp((r - lo1) / (lo2 - lo1), 0, 1)
  falling <- terra::clamp((hi2 - r) / (hi2 - hi1), 0, 1)
  m <- terra::ifel(r <= lo2, rising,
         terra::ifel(r <= hi1, 1, falling))
  m <- terra::clamp(m, 0, 1)
  cca_log(sprintf("%s membership built (trapezoid %g/%g/%g/%g)",
                  name, lo1, lo2, hi1, hi2))
  m
}

# ---- Get soil and climate drivers via geodata --------------------------
# SoilGrids via geodata::soil_world returns 1 km rasters at ~250 m re-sampled
# to 1 km. We pull pH, OC, clay at 0-30 cm depth.
# geodata caches downloads in `path = `; the cca_download sidecar is not
# strictly needed because geodata itself is authoritative.
soil_cache <- file.path(CCA_PATHS$inputs_root, "soilgrids")
dir.create(soil_cache, recursive = TRUE, showWarnings = FALSE)

cca_log("Pulling SoilGrids pH(H2O), OC, clay (0-30 cm)...")
# `var` names map to SoilGrids ISRIC conventions. `depth = 30` is the most
# recent layer (averaged 0-30 cm). Empirical unit check against the
# geodata-served files (2026-04-17):
#   phh2o : delivered in pH units directly (range 3.5-10.2). NO rescale.
#   soc   : delivered in g/kg (range 0.2-508). Divide by 10 for percent.
#   clay  : delivered in percent already (range 0-98). NO rescale.
# An earlier revision of this script divided all three by 10 (carried over
# from the raw ISRIC SoilGrids WCS convention). That drove pH memberships to
# zero and collapsed suitability across the CCA region -- caught by the
# staging raster having mean ~1e-4 instead of ~0.4. Fixed below.
ph   <- geodata::soil_world(var = "phh2o", depth = 30, stat = "mean", path = soil_cache)
oc   <- geodata::soil_world(var = "soc",   depth = 30, stat = "mean", path = soil_cache) / 10   # g/kg -> percent
clay <- geodata::soil_world(var = "clay",  depth = 30, stat = "mean", path = soil_cache)

# Mean annual rainfall: geodata::worldclim_global gives WorldClim 2.1 with
# 12 monthly bands at 10 arcmin; sum to annual mean. CHIRPS would be more
# precise but requires the atlas seasonality pipeline run -- we use WorldClim
# as the always-available baseline; callers can swap in CHIRPS annual mean
# if they have it.
cca_log("Pulling WorldClim 2.1 monthly precipitation (10 arcmin)...")
prec_monthly <- geodata::worldclim_global(var = "prec", res = 10, path = soil_cache)
mar <- terra::app(prec_monthly, sum)     # mm/yr

# ---- Reproject everything to CCA template ------------------------------
template <- terra::rast(cca_template_path())

ph_r   <- terra::project(ph,   template, method = "average")
oc_r   <- terra::project(oc,   template, method = "average")
clay_r <- terra::project(clay, template, method = "average")
mar_r  <- terra::project(mar,  template, method = "average")

# ---- Apply membership functions ----------------------------------------
m_ph    <- fuzzy_trapezoid(ph_r,   lo1 = 4.5, lo2 = 5.5, hi1 = 7.5, hi2 = 8.5, name = "pH")
m_oc    <- fuzzy_trapezoid(oc_r,   lo1 = 0.2, lo2 = 0.5, hi1 = 3.0, hi2 = 5.0, name = "OC(%)")
m_clay  <- fuzzy_trapezoid(clay_r, lo1 = 5,   lo2 = 10,  hi1 = 40,  hi2 = 60,  name = "clay(%)")
m_rain  <- fuzzy_trapezoid(mar_r,  lo1 = 300, lo2 = 400, hi1 = 2500, hi2 = 3500, name = "rainfall(mm)")

# Geometric mean (product^(1/n)) => any zero driver zeroes the pixel.
# Small epsilon avoids NaN when a driver is exactly 0.
eps <- 1e-9
suit <- (m_ph + eps) * (m_oc + eps) * (m_clay + eps) * (m_rain + eps)
suit <- suit ^ (1/4)
suit <- terra::clamp(suit, 0, 1)

cca_assert_range(suit, 0, 1, name = "biochar suitability", strict = TRUE)
cca_assert_not_allna(suit, name = "biochar suitability")

out_tif <- file.path(staging_dir, sprintf("biochar_suitability_%s.tif", CCA_RES_LABEL))
cca_write_cog(suit, out_tif)

# ---- Side-car: CropSuiteLite-compatible .inf --------------------------
# For analysts who want to rerun with the full CropSuiteLite growing-cycle
# machinery. This file is NOT used by this script.
inf_path <- file.path(staging_dir, "biochar.inf")
writeLines(c(
  "name = \t\tbiochar",
  "growing_cycle = \t0",
  "ph_vals = \t\t4.5,5.5,7.5,8.5",
  "ph_suit = \t\t0,1,1,0",
  "organic_carbon_vals = \t0.2,0.5,3.0,5.0",
  "organic_carbon_suit = \t0,1,1,0",
  "clay_vals = \t\t5,10,40,60",
  "clay_suit = \t\t0,1,1,0",
  "prec_vals = \t\t300,400,2500,3500",
  "prec_suit = \t\t0,1,1,0"
), inf_path)
cca_log("Wrote CropSuiteLite .inf side-car: ", inf_path)

cca_banner("Done: biochar suitability built")
