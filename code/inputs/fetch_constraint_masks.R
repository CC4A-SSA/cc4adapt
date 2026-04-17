# ------------------------------------------------------------------------------
# code/inputs/fetch_constraint_masks.R
#
# Fetches layers used for OPTIONAL constraint masks in Step 4 integration:
#   * Protected areas (WDPA)
#   * Slope (derived from SRTM 30 s via geodata::elevation_global)
#   * Urban extent (WorldPop / GHSL SMOD fallback)
#   * Irrigation extent (FAO GMIA or GAEZ)
#
# --- No default thresholds applied ---
#
# Per project-lead rule, this script does NOT apply a slope cut-off, a
# protected-area exclusion, or an urban mask. It produces raw rasters and
# vectors; cca_mask_by_threshold / cca_mask_by_polygon wrappers apply them
# downstream only when callers pass a non-NULL parameter.
#
# --- Sources ---
#
#   * WDPA: https://www.protectedplanet.net/en (monthly CSV/SHP release,
#     stable URL pattern). Alternative: `wdpar` R package. We use direct
#     download because `wdpar` depends on a headless browser not installed
#     in atlas containers.
#   * Elevation: geodata::elevation_global(res = 0.5 * 60) returns SRTM at
#     30 arcsec which we reproject/resample to derive slope at CCA_RES.
#   * Urban: GHSL SMOD 2020 from JRC (stable URLs). WorldPop is an
#     alternative; GHSL is smaller and has built-in urban-rural classes.
#   * Irrigation: FAO AQUASTAT GMIA v5 (stable URL from FAO).
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

cca_banner("Fetch OPTIONAL constraint masks (WDPA, slope, urban, irrigation)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "constraints")
staging_dir <- file.path(CCA_PATHS$staging_root, "constraints")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

# Skip-if-all-present guard. The three ingredients that we can reliably
# auto-produce are per-country slope rasters and the single GHSL SMOD COG.
# WDPA + GMIA are manual-drop so they are NOT in the skip list -- this
# script should re-check them on every run (cheap).
auto_outputs <- c(
  file.path(staging_dir, sprintf("slope_deg_%s_%s.tif", CCA_COUNTRIES, CCA_RES_LABEL)),
  file.path(staging_dir, sprintf("ghsl_smod_%s.tif", CCA_RES_LABEL))
)
if (cca_skip_if_exists(auto_outputs, tag = "constraint masks (slope + GHSL SMOD)")) quit(status = 0)

template <- terra::rast(cca_template_path())

# ---- 1. Protected areas (WDPA) -----------------------------------------
# Protected Planet publishes monthly archives at a predictable URL:
#   https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_<MonYYYY>_Public_shp.zip
# We resolve "current" via redirect; if that fails we fall back to the
# country-specific subsets (per-ISO3, much smaller).
wdpa_for_iso3 <- function(iso3) {
  # Protected Planet requires a date-stamped slug plus a session cookie for
  # the public `WDPA_<MonYYYY>_<ISO3>` package. The "current" alias does
  # not work outside the browser flow (returns 404). Easiest reliable path:
  # wait for the `wdpar` R package (adds headless browser) OR require a
  # manual download of the per-country zip into data/inputs/constraints/wdpa/.
  # We check for a manually-dropped zip first; only warn if absent.
  local_candidates <- list.files(
    file.path(raw_dir, "wdpa"),
    pattern = sprintf("WDPA_.*_%s\\.zip$", iso3),
    full.names = TRUE, ignore.case = TRUE
  )
  if (!length(local_candidates)) {
    cca_log(sprintf(
      "WDPA for %s not available. Manual step: download ",
      iso3),
      sprintf("https://www.protectedplanet.net/country/%s and drop the zip into %s",
              iso3, file.path(raw_dir, "wdpa")),
      level = "WARN")
    return(invisible(NULL))
  }
  dest <- local_candidates[1]
  exdir <- file.path(raw_dir, sprintf("wdpa_%s_extracted", iso3))
  if (!dir.exists(exdir) || !length(list.files(exdir))) utils::unzip(dest, exdir = exdir)
  shps <- list.files(exdir, pattern = "polygons.*\\.shp$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  if (!length(shps)) { cca_log("No polygons shp in WDPA zip for ", iso3, level = "WARN"); return(NULL) }
  v <- terra::vect(shps[1])
  out_v <- file.path(staging_dir, sprintf("wdpa_%s.gpkg", iso3))
  cca_write_vector(v, out_v)
  invisible(out_v)
}
for (iso3 in CCA_COUNTRIES) wdpa_for_iso3(iso3)

# ---- 2. Slope from SRTM (via geodata) ----------------------------------
# geodata::elevation_global pulls a global mosaic; elevation_30s returns
# per-country SRTM. For Africa-wide we need the per-country pulls or the
# global low-res variant. We use the country-level pulls at 30 arcsec for
# accuracy and then derive slope via terra::terrain(..., v="slope").
for (iso3 in CCA_COUNTRIES) {
  dem <- tryCatch(
    geodata::elevation_30s(country = iso3, path = file.path(raw_dir, "dem")),
    error = function(e) { cca_log(conditionMessage(e), level = "WARN"); NULL }
  )
  if (is.null(dem)) next
  # terra::terrain expects projected meters; for lat/lon we pass
  # `unit = "degrees"` and convert downstream. Since slope is only used
  # optionally as a threshold mask (no default), we store as degrees.
  slope_deg <- terra::terrain(dem, v = "slope", unit = "degrees")
  # Resample to CCA_RES using "average" (slope is a density-like quantity).
  template_iso <- terra::rast(cca_template_path(iso3))
  slope_cca <- terra::project(slope_deg, template_iso, method = "average")
  out_tif <- file.path(staging_dir, sprintf("slope_deg_%s_%s.tif", iso3, CCA_RES_LABEL))
  cca_write_cog(slope_cca, out_tif)
}

# ---- 3. Urban extent (GHSL SMOD) --------------------------------------
# GHSL SMOD 2020 at 1 km, global. JRC publishes ZIPs per release. Exact URL
# includes year + epoch; use the 2023A release 1 km product.
ghsl_url  <- "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000/V1-0/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V1_0.zip"
ghsl_zip  <- file.path(raw_dir, "GHS_SMOD_2020_1km.zip")
ok <- tryCatch({ cca_download(ghsl_url, ghsl_zip); TRUE },
               error = function(e) { cca_log(conditionMessage(e), level = "WARN"); FALSE })
if (ok) {
  exdir <- file.path(raw_dir, "ghs_smod_extracted")
  if (!dir.exists(exdir) || !length(list.files(exdir))) utils::unzip(ghsl_zip, exdir = exdir)
  tifs <- list.files(exdir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (length(tifs)) {
    smod <- terra::rast(tifs[1])
    # GHSL ships in Mollweide (EPSG:54009). Reproject to EPSG:4326 first,
    # then resample to template. Use "near" because SMOD is categorical.
    smod_ll <- terra::project(smod, "EPSG:4326", method = "near")
    smod_cca <- terra::project(smod_ll, template, method = "near")
    out_tif <- file.path(staging_dir, sprintf("ghsl_smod_%s.tif", CCA_RES_LABEL))
    cca_write_cog(smod_cca, out_tif)
    cca_log("Wrote GHSL SMOD 2020 (categorical): ", out_tif)
  }
}

# ---- 4. Irrigation (FAO GMIA) -----------------------------------------
# FAO AQUASTAT GMIA v5 (Siebert et al. 2013) is the standard global
# irrigation extent at 5 arcmin. FAO has moved GMIA between hosting
# platforms several times (aquastat -> earthengine -> data.apps.fao.org);
# the previous ASC .gz endpoint now serves an HTML landing page. Until a
# stable URL is agreed with the team we require a manual drop into
# data/inputs/constraints/gmia/ of either:
#   * gmia_v5_aei_pct_aei.asc  (raw ASC, no .prj sidecar)
#   * gmia_v5_aei_pct_aei.tif  (GeoTIFF export)
# The atlas has a working GMIA ingest at atlas/atlas_exposure/; we defer
# to that for the canonical version.
gmia_dir <- file.path(raw_dir, "gmia")
dir.create(gmia_dir, recursive = TRUE, showWarnings = FALSE)
gmia_candidates <- list.files(
  gmia_dir,
  pattern = "gmia_v5_aei_pct_aei\\.(asc|tif)$",
  full.names = TRUE, ignore.case = TRUE
)
if (!length(gmia_candidates)) {
  cca_log("GMIA v5 irrigation raster not available. Manual step: download from ",
          "FAO AQUASTAT GMIA v5 (Siebert et al. 2013) and place the file in ",
          gmia_dir, ". Alternative: reuse atlas/atlas_exposure/ output.",
          level = "WARN")
} else {
  gmia <- terra::rast(gmia_candidates[1])
  if (as.character(terra::crs(gmia, describe = TRUE)$code) != "4326") {
    terra::crs(gmia) <- CCA_CRS
  }
  gmia_cca <- terra::project(gmia, template, method = "average")
  out_tif <- file.path(staging_dir, sprintf("gmia_aei_pct_%s.tif", CCA_RES_LABEL))
  cca_write_cog(gmia_cca, out_tif)
  cca_log("Wrote FAO GMIA irrigation extent: ", out_tif)
}

cca_banner("Done: constraint layers staged")
cca_log("Remember: these are INPUT LAYERS ONLY. Apply as masks explicitly via ",
        "cca_mask_by_polygon(...) or cca_mask_by_threshold(...) with non-NULL args.")
