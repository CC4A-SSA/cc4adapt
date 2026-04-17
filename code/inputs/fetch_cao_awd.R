# ------------------------------------------------------------------------------
# code/inputs/fetch_cao_awd.R
#
# Fetches / constructs an AWD rice methane-mitigation raster for Africa.
#
# --- The project plan cites Cao et al. 2026 ---
#
# The investment document (Activity 2 Step 1) refers to "Cao et al. (2026) for
# alternate wetting and drying in rice" -- a paper that is listed as
# anticipated at the time of writing (April 2026). A published Nature-Food
# raster may or may not be available yet. This script tries two things in
# sequence:
#
#   1. Attempt to resolve a Cao 2026 Zenodo / Dataverse DOI via the
#      respective APIs. If a record exists and publishes a global raster,
#      download and stage it.
#
#   2. Fallback: construct an IPCC Tier-1 AWD potential surface from first
#      principles, using:
#        * MapSPAM 2020 irrigated + rainfed rice harvested area (already
#          fetched by atlas/atlas_exposure/)
#        * Regional IPCC 2019 Refinement Tier-1 emission factors for
#          continuous-flood (EFc) vs multiple-drainage (EFm) rice
#        * GWP20 for CH4 (IPCC AR6) to match the Naturebase crp_rec
#          accounting (which also uses GWP20 -- documented in their metadata)
#
#      The delta (EFc - EFm) in kg CH4 / ha / season, multiplied by the
#      number of seasons per year, by MapSPAM area, and by GWP20 * (1/1000)
#      gives tCO2e/ha/yr. This is a coarser but fully reproducible fallback
#      that we clearly flag in outputs (`method_tag = "IPCC_T1_fallback"`).
#
# --- Explicit caveats ---
#
#   * Tier-1 factors are regional averages; real mitigation depends on
#     soil, cultivar and drainage duration. Where Cao 2026 becomes available
#     it should supersede this.
#   * This script DOES NOT mask to a minimum-irrigation threshold or a
#     minimum rice-area threshold. Those are stakeholder-decided in
#     Activity 1 and are wrapped as optional functions downstream.
# ------------------------------------------------------------------------------

# ---- Preamble ----------------------------------------------------------
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
})

cca_banner("Fetch AWD rice mitigation (Cao 2026 primary; IPCC T1 fallback)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "awd")
staging_dir <- file.path(CCA_PATHS$staging_root, "awd")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

awd_primary  <- file.path(staging_dir, "awd_potential_tco2e_ha_yr_Cao2026.tif")
awd_fallback <- file.path(staging_dir, "awd_potential_tco2e_ha_yr_T1fallback.tif")
if (cca_skip_if_exists(awd_primary, tag = "awd_potential (Cao 2026)")) quit(status = 0)
if (cca_skip_if_exists(awd_fallback, tag = "awd_potential (T1 fallback)")) quit(status = 0)

# ---- 1. Try Cao 2026 ---------------------------------------------------
# The Nature Food norm is to deposit data on Zenodo. We search Zenodo for
# recent records tagged with "AWD" and Cao. This deliberately over-queries
# so that any matching dataset is caught; manual curation still needed.
zenodo_search <- function(query, size = 5L) {
  api <- sprintf('https://zenodo.org/api/records?q=%s&size=%d&sort=mostrecent',
                 utils::URLencode(query, reserved = TRUE), size)
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({ cca_download(api, tmp, force = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok) return(list())
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("[CCA] jsonlite required. install.packages('jsonlite').")
  }
  j <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  j$hits$hits %||% list()
}

cao_hits <- zenodo_search('"alternate wetting and drying" creators.name:"Cao"', 5L)
if (length(cao_hits)) {
  cca_log(sprintf("Zenodo returned %d candidate Cao AWD record(s); manual review needed",
                  length(cao_hits)))
  for (h in cao_hits) {
    title <- h$metadata$title %||% "(untitled)"
    doi   <- h$doi %||% h$metadata$doi %||% ""
    cca_log(sprintf("  - '%s' (DOI %s)", title, doi))
  }
  cca_log("Auto-download of Cao 2026 not performed: verify record identity manually ",
          "and add the confirmed DOI to fetch_cao_awd.R before re-running.",
          level = "WARN")
} else {
  cca_log("No Cao 2026 record found on Zenodo. Proceeding to IPCC Tier-1 fallback.",
          level = "WARN")
}

# ---- 2. IPCC Tier-1 fallback ------------------------------------------
# Reference: IPCC 2019 Refinement, Vol 4 Ch 5.5 (Rice Cultivations), Table 5.11.
# We use the Africa-regional value set as the Tier-1 default for Sub-Saharan
# Africa. These are KG CH4 / HA / DAY values; we convert to a seasonal
# emission by multiplying by growing-season length.
#
# NOTE on EF units: IPCC Tier-1 baseline EFc is 1.19 kg CH4/ha/day for
# continuous flooding, regionally adjusted. Multiple-drainage (AWD) scaling
# factor SFw = 0.71 for single aeration; 0.55 for multiple drainage. We use
# SFw = 0.55 to align with Naturebase crp_rec's "multiple drainage" framing.
#
# The mitigation delta per pixel is therefore:
#   Delta_CH4_kg_ha_yr = EFc * day_count * (1 - SFw) * n_seasons
#   Mitigation_tCO2e_ha_yr = Delta_CH4_kg_ha_yr * GWP20_CH4 / 1000
#   area-weighted = Mitigation * rice_area_frac_per_pixel
#
# With regional, literature-grounded constants:
IPCC_EFc_kg_ch4_ha_day <- 1.19   # Tier-1 continuous flood baseline (Africa)
IPCC_SFw_multiple_drain <- 0.55  # multiple drainage scaling factor
GROWING_DAYS   <- 120            # single-season rice; Tier 1 default
N_SEASONS      <- 1              # NGA/KEN/ETH: predominantly single-season
GWP20_CH4      <- 80             # IPCC AR6 (for consistency with Naturebase)

# MapSPAM rice harvested area (ha/pixel): we expect this to live under
# atlas_exposure outputs; if not yet available we warn and skip construction.
mapspam_rice_candidates <- c(
  file.path(CCA_PATHS$inputs_root, "mapspam", "SPAM2020_rice_H.tif"),     # preferred
  file.path(CCA_PATHS$staging_root, "mapspam", "rice_harvested_ha.tif")
)
mapspam_rice <- mapspam_rice_candidates[file.exists(mapspam_rice_candidates)][1]
if (is.na(mapspam_rice)) {
  cca_log("No MapSPAM rice raster found. IPCC T1 fallback needs it; see ",
          "atlas/atlas_exposure/0_download_data.r. Skipping raster construction for now.",
          level = "WARN")
  cca_banner("Done (partial): AWD raster pending MapSPAM rice input")
  quit(status = 0)
}

cca_log("Using MapSPAM rice area: ", mapspam_rice)
rice_ha <- terra::rast(mapspam_rice)

# Compute per-pixel AWD mitigation
delta_tco2e_per_ha <- IPCC_EFc_kg_ch4_ha_day * GROWING_DAYS *
                      (1 - IPCC_SFw_multiple_drain) * N_SEASONS *
                      GWP20_CH4 / 1000
cca_log(sprintf("IPCC T1 delta = %.3f tCO2e/ha/yr (uniform across pixels)",
                delta_tco2e_per_ha))

# This is a uniform density multiplied by rice-presence. We express as a
# density raster (tCO2e/ha/yr) where pixels with rice area > 0 get the
# regional EF, and pixels without rice are NA. Downstream step1_mitigation.R
# multiplies by MapSPAM rice fraction if a caller opts in.
awd_density <- terra::ifel(rice_ha > 0, delta_tco2e_per_ha, NA)

# Clip to CCA region and write
clipped <- cca_clip_to_cca_region(awd_density, cca_template_path())
cca_assert_not_allna(clipped, name = "AWD T1 fallback")
cca_assert_range(clipped, lo = 0, hi = 50, name = "AWD T1 fallback (tCO2e/ha/yr)")

out_tif <- file.path(staging_dir, "awd_potential_tco2e_ha_yr_T1fallback.tif")
cca_write_cog(clipped, out_tif)

# Tag the provenance sidecar with the method so downstream code can route
# correctly once Cao 2026 becomes available.
writeLines(
  c("{",
    '  "method": "IPCC_Tier1_fallback",',
    '  "EFc_kg_ch4_ha_day": ', IPCC_EFc_kg_ch4_ha_day, ",",
    '  "SFw": ',               IPCC_SFw_multiple_drain, ",",
    '  "GWP20_CH4": ',         GWP20_CH4, ",",
    '  "growing_days": ',      GROWING_DAYS, ",",
    '  "n_seasons": ',         N_SEASONS,
    "}"),
  paste0(out_tif, ".method.json")
)

cca_banner("Done: AWD Tier-1 fallback staged")
