# ------------------------------------------------------------------------------
# code/inputs/build_erw_suitability.R
#
# Builds a 0-1 suitability surface for Enhanced Rock Weathering (ERW) at 5
# arcmin over the CCA region.
#
# Drivers (documented rationale):
#
#   * Soil pH: acidic soils weather basalt fastest; optimum below ~6.5.
#       Goll et al. 2021 / Beerling et al. 2020 show near-linear acceleration
#       of silicate dissolution with decreasing pH.
#   * Clay content: moderate-to-high clays increase particle-water contact
#       time and cation exchange capacity; optimum >= 15 percent.
#       Evidence: Dietzen & Rosing 2023 weathering-kinetics review.
#   * Mean annual temperature: warmer soils weather faster via Arrhenius;
#       tropical optimum ~20-30 C.
#       Evidence: Kantola et al. 2017 "Potential of global croplands and
#       bioenergy crops for climate change mitigation through deployment
#       for enhanced weathering."
#   * Mean annual rainfall: >= 600 mm to keep pores wet enough for
#       dissolution; >2500 mm increases leaching loss.
#
# This script does NOT apply the basalt-proximity mask. That lives in
# cca_mask_by_proximity(..., distance_km = NULL) and is only applied if the
# caller provides a stakeholder-decided distance parameter.
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

cca_banner("Build ERW suitability surface (pH, clay, temperature, rainfall)")

staging_dir <- file.path(CCA_PATHS$staging_root, "suitability")
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

staged_tif <- file.path(staging_dir, sprintf("erw_suitability_%s.tif", CCA_RES_LABEL))
if (cca_skip_if_exists(staged_tif, tag = "erw_suitability")) quit(status = 0)

# Reuse the trapezoid helper from build_biochar_suitability.R by re-defining
# it locally (both scripts run as standalone entry points).
fuzzy_trapezoid <- function(r, lo1, lo2, hi1, hi2, name = "fuzzy") {
  stopifnot(lo1 <= lo2, lo2 <= hi1, hi1 <= hi2)
  rising  <- terra::clamp((r - lo1) / (lo2 - lo1), 0, 1)
  falling <- terra::clamp((hi2 - r) / (hi2 - hi1), 0, 1)
  m <- terra::ifel(r <= lo2, rising,
         terra::ifel(r <= hi1, 1, falling))
  m <- terra::clamp(m, 0, 1)
  cca_log(sprintf("%s membership built (trapezoid %g/%g/%g/%g)",
                  name, lo1, lo2, hi1, hi2))
  m
}

# ---- Pull drivers via geodata ------------------------------------------
soil_cache <- file.path(CCA_PATHS$inputs_root, "soilgrids")
dir.create(soil_cache, recursive = TRUE, showWarnings = FALSE)

cca_log("Pulling SoilGrids pH(H2O) and clay (0-30 cm)...")
# Empirical unit check (see build_biochar_suitability.R):
#   phh2o : pH units (no rescale); clay : percent (no rescale).
ph   <- geodata::soil_world(var = "phh2o", depth = 30, stat = "mean", path = soil_cache)
clay <- geodata::soil_world(var = "clay",  depth = 30, stat = "mean", path = soil_cache)

cca_log("Pulling WorldClim 2.1 mean annual temperature (bio_1, 10 arcmin)...")
# WorldClim bio_1 is annual mean temperature in degrees C times 10 in older
# releases. The 2.1 release uses degrees C directly. geodata handles the
# scaling for us.
worldclim_cache <- soil_cache
tmean <- geodata::worldclim_global(var = "bio", res = 10, path = worldclim_cache)[[1]]

cca_log("Pulling WorldClim 2.1 monthly precipitation...")
prec_monthly <- geodata::worldclim_global(var = "prec", res = 10, path = worldclim_cache)
mar <- terra::app(prec_monthly, sum)

# ---- Reproject to CCA template ----------------------------------------
template <- terra::rast(cca_template_path())

ph_r    <- terra::project(ph,    template, method = "average")
clay_r  <- terra::project(clay,  template, method = "average")
tmean_r <- terra::project(tmean, template, method = "average")
mar_r   <- terra::project(mar,   template, method = "average")

# ---- Memberships -------------------------------------------------------
# pH: full suitability at <=6.5, dropping to 0 at pH 8.
m_ph    <- fuzzy_trapezoid(ph_r,    lo1 = 3.5, lo2 = 4.5, hi1 = 6.5, hi2 = 8.0, name = "pH")
# clay: suitable at >=15%, with gentle taper above ~40% (very heavy clays
# reduce aeration).
m_clay  <- fuzzy_trapezoid(clay_r,  lo1 = 5,   lo2 = 15,  hi1 = 45,  hi2 = 65,  name = "clay(%)")
# temperature: tropical optimum 20-30 C.
m_temp  <- fuzzy_trapezoid(tmean_r, lo1 = 10,  lo2 = 20,  hi1 = 30,  hi2 = 35,  name = "MAT(C)")
# rainfall: same as biochar driver (the pore-water argument is similar).
m_rain  <- fuzzy_trapezoid(mar_r,   lo1 = 400, lo2 = 600, hi1 = 2500, hi2 = 3500, name = "rainfall(mm)")

# Geometric mean, same rationale as biochar suitability
eps <- 1e-9
suit <- (m_ph + eps) * (m_clay + eps) * (m_temp + eps) * (m_rain + eps)
suit <- suit ^ (1/4)
suit <- terra::clamp(suit, 0, 1)

cca_assert_range(suit, 0, 1, name = "ERW suitability", strict = TRUE)
cca_assert_not_allna(suit, name = "ERW suitability")

out_tif <- file.path(staging_dir, sprintf("erw_suitability_%s.tif", CCA_RES_LABEL))
cca_write_cog(suit, out_tif)

inf_path <- file.path(staging_dir, "erw.inf")
writeLines(c(
  "name = \t\terw",
  "growing_cycle = \t0",
  "ph_vals = \t\t3.5,4.5,6.5,8.0",
  "ph_suit = \t\t0,1,1,0",
  "clay_vals = \t\t5,15,45,65",
  "clay_suit = \t\t0,1,1,0",
  "temp_vals = \t\t10,20,30,35",
  "temp_suit = \t\t0,1,1,0",
  "prec_vals = \t\t400,600,2500,3500",
  "prec_suit = \t\t0,1,1,0"
), inf_path)
cca_log("Wrote CropSuiteLite .inf side-car: ", inf_path)

cca_banner("Done: ERW suitability built")
