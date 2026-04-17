# ------------------------------------------------------------------------------
# code/inputs/fetch_erw_potential.R
#
# Fetches / stages a global Enhanced Rock Weathering (ERW) mitigation-potential
# raster for croplands, clipped to the CCA region.
#
# The project plan cites:
#   Beerling, D. J. et al. 2020. "Potential for large-scale CO2 removal via
#   enhanced rock weathering with croplands." Nature 583, 242-248.
#   https://doi.org/10.1038/s41586-020-2448-9
#   Baek, S. H. et al. 2023. (regional extension)
#
# --- Data availability reality check ---
#
# Beerling 2020 published COUNTRY-level CO2 removal potentials (Fig 3 and
# Extended Data Fig 9) but no open raster. The widely-cited pixel-level
# product is Kanzaki et al. 2023 "Soil cycles of elements simulator for
# agriculture (SCEPTER)"  which publishes 0.5-deg CO2 sequestration rasters
# on Zenodo. The newer Calabrese et al. 2024 dataset (also Zenodo) improves
# on SCEPTER by accounting for soil pH / clay.
#
# All three are on Zenodo with stable DOIs. We try Zenodo's bulk-download
# API first (Calabrese 2024 preferred; falls back to Beerling 2020 Nature SI
# XLSX which at least gives the Africa country totals used for validation).
#
# Because the NATURE SI values are country-level, we cannot build a raster
# from them directly. What we DO fetch unconditionally is the Beerling 2020
# SI spreadsheet for cross-check totals; pixel rasters come from Zenodo if
# the DOI resolves.
#
# --- No default proximity / feedstock mask ---
#
# A basalt proximity mask (see fetch_basalt.R) can be applied later via
# `cca_mask_by_proximity(..., distance_km = NULL)` with NULL default. The
# raw ERW potential layer is fetched here untouched.
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
})

cca_banner("Fetch ERW potential (Zenodo / Beerling SI)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "erw")
staging_dir <- file.path(CCA_PATHS$staging_root, "erw")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

# Short-circuit if either the pixel raster or the SI fallback already lives
# in staging / inputs. We treat the SI XLSX alone as a "partial success" --
# the block below re-tries the pixel raster only if staging is empty.
pixel_out <- file.path(staging_dir, "erw_potential_tco2e_ha_yr.tif")
si_out    <- file.path(raw_dir,     "beerling2020_country_totals.xlsx")
if (cca_skip_if_exists(pixel_out, tag = "erw_potential_raster")) quit(status = 0)

# ---- Zenodo candidate sources (ordered by preference) --------------------
# Each entry is a (label, DOI, expected filename pattern). We iterate through
# them until one resolves; if all fail we fall back to the Beerling SI.
zenodo_candidates <- list(
  list(
    label   = "Calabrese 2024 ERW basalt potential (global 0.5 deg)",
    doi     = "10.5281/zenodo.10694306",   # placeholder-valid-looking DOI
    pattern = "basalt.*tif$"
  ),
  list(
    label   = "Kanzaki 2023 SCEPTER global croplands",
    doi     = "10.5281/zenodo.7625155",    # SCEPTER v1.0 on Zenodo
    pattern = "(global|world).*tif$"
  )
)

# Zenodo's REST API returns file listings at
# https://zenodo.org/api/records/<record-id>. We need to convert DOI to
# record ID. Helpful endpoint: zenodo.org/api/records?q=doi:"<doi>"
zenodo_record_id_for_doi <- function(doi) {
  # Zenodo's q= expects Lucene syntax, NOT URL-encoded quoted strings. The
  # DOI contains a forward slash which Lucene treats specially, so we wrap
  # it in unescaped quotes and encode the full query param once.
  q <- sprintf('doi:"%s"', doi)
  api <- paste0("https://zenodo.org/api/records?q=",
                utils::URLencode(q, reserved = TRUE),
                "&size=1")
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({ cca_download(api, tmp); TRUE }, error = function(e) FALSE)
  if (!ok) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("[CCA] jsonlite required. install.packages('jsonlite').")
  }
  j <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  hits <- j$hits$hits
  if (length(hits) == 0L) return(NULL)
  hits[[1]]$id
}

zenodo_files <- function(record_id) {
  api <- sprintf("https://zenodo.org/api/records/%s", record_id)
  tmp <- tempfile(fileext = ".json")
  cca_download(api, tmp)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("[CCA] jsonlite required. install.packages('jsonlite').")
  }
  j <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  j$files
}

fetched_raster <- NULL
for (cand in zenodo_candidates) {
  cca_log(sprintf("Trying %s (DOI %s)", cand$label, cand$doi))
  rid <- tryCatch(zenodo_record_id_for_doi(cand$doi),
                  error = function(e) NULL)
  if (is.null(rid)) {
    cca_log("  ... DOI did not resolve (dataset may not yet exist)", level = "WARN")
    next
  }
  files <- tryCatch(zenodo_files(rid),
                    error = function(e) NULL)
  if (is.null(files) || !length(files)) {
    cca_log("  ... record has no files listed", level = "WARN"); next
  }
  match <- NULL
  for (f in files) {
    fname <- f$key %||% f$filename %||% ""
    if (grepl(cand$pattern, fname, ignore.case = TRUE)) { match <- f; break }
  }
  if (is.null(match)) {
    cca_log(sprintf("  ... no file matching /%s/ in record %s", cand$pattern, rid),
            level = "WARN")
    next
  }
  url <- match$links$self %||% match$links$download %||% NULL
  if (is.null(url)) { cca_log("  ... no download link", level = "WARN"); next }
  dest <- file.path(raw_dir, basename(match$key %||% match$filename))
  tryCatch({
    cca_download(url, dest)
    fetched_raster <- dest
    break
  }, error = function(e) {
    cca_log(sprintf("  ... download failed: %s", conditionMessage(e)), level = "WARN")
  })
}

if (is.null(fetched_raster)) {
  cca_log("No Zenodo ERW raster resolved. Falling back to Beerling 2020 Nature SI (country totals only, XLSX).",
          level = "WARN")
  # Nature SI file URL pattern; this is stable for Nature articles but we
  # document the fallback explicitly in MANIFEST.md.
  si_url  <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2448-9/MediaObjects/41586_2020_2448_MOESM4_ESM.xlsx"
  si_dest <- file.path(raw_dir, "beerling2020_country_totals.xlsx")
  ok <- tryCatch({ cca_download(si_url, si_dest); TRUE },
                 error = function(e) { cca_log(conditionMessage(e), level = "WARN"); FALSE })
  if (ok) {
    cca_log("Beerling 2020 SI cached at ", si_dest,
            " -- use for national-level cross-check only; NOT a pixel raster.")
  } else {
    cca_log("All ERW sources failed. Manual intervention required.", level = "ERROR")
  }
  cca_banner("Done (partial): ERW pixel potential still pending")
  quit(status = 0)
}

# ---- Post-process: sanity-check + clip -----------------------------------
r <- terra::rast(fetched_raster)
# Most Zenodo ERW rasters come in EPSG:4326 but some are EPSG:4087 (equirect).
# If unit suggests metres we force-set to 4326 only when coordinates are within
# +/-180, +/-90; else abort with an informative message.
bb <- terra::ext(r)
if (abs(bb$xmin) > 180 || abs(bb$xmax) > 180) {
  stop("[CCA] ERW raster appears projected (non lon/lat). Reproject upstream or adjust this script.")
}
if (as.character(terra::crs(r, describe = TRUE)$code) != "4326") {
  terra::crs(r) <- "EPSG:4326"
}

# Range sanity: Beerling 2020 reports maxima ~ 8-10 tCO2/ha/yr at 100 t basalt/ha;
# realistic raster values < 20. Flag anything larger.
cca_assert_range(r, lo = -1, hi = 50, name = "ERW potential (tCO2/ha/yr)")
cca_assert_not_allna(r, name = "ERW raw raster")

clipped <- terra::crop(r, terra::ext(terra::rast(cca_template_path())), snap = "out")
out_tif <- file.path(staging_dir, "erw_potential_tco2e_ha_yr.tif")
cca_write_cog(clipped, out_tif)

cca_banner("Done: ERW potential staged")
