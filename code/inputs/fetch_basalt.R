# ------------------------------------------------------------------------------
# code/inputs/fetch_basalt.R
#
# Fetches geological data used for the OPTIONAL Enhanced Rock Weathering
# proximity mask. The mask itself is not applied by default -- it is wrapped
# in cca_mask_by_proximity(..., distance_km = NULL). This script only
# provides the input layer(s).
#
# --- Source choice: GLiM (Global Lithological Map) ---
#
# Hartmann J., Moosdorf N. 2012. "The new global lithological map database
# GLiM: A representation of rock properties at the Earth surface."
# G3 (Geochemistry, Geophysics, Geosystems) 13 (12).
# https://doi.org/10.1594/PANGAEA.788537  (PANGAEA archive)
#
# GLiM publishes a global shapefile of surface lithology with ~17 classes.
# For ERW we extract "mafic volcanics" and "intermediate volcanics" as the
# proxy for basaltic outcrops suitable as ERW feedstock. Other options
# considered:
#   * USGS Mineral Resources Data System (MRDS) -- point deposits; biased
#     to mined occurrences rather than total surface outcrop. Useful
#     supplement but not the primary layer.
#   * One Geology (1:1M) -- sparse coverage in SSA, skipped.
#
# --- Why both GLiM AND MRDS ---
#
# GLiM captures outcrop extent; MRDS captures operating mines that could
# supply crushed basalt. Proximity-to-GLiM answers "where could basalt be
# mined", proximity-to-MRDS answers "where is basalt being mined right now".
# Stakeholder workshops will decide which is relevant (and at what distance).
# We fetch both; applying them remains optional at the call site.
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

cca_banner("Fetch basalt/mafic geology (GLiM + USGS MRDS)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "basalt")
staging_dir <- file.path(CCA_PATHS$staging_root, "basalt")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

glim_out <- file.path(staging_dir, "glim_africa_mafic_intermediate.gpkg")
mrds_out <- file.path(staging_dir, "mrds_basalt_points.gpkg")
# Both optional; skip only if BOTH already staged.
if (cca_skip_if_exists(c(glim_out, mrds_out), tag = "basalt layers")) quit(status = 0)

# ---- 1. GLiM global shapefile ------------------------------------------
# PANGAEA delivers the ZIP via DOI 10.1594/PANGAEA.788537. The direct file
# URL pattern has changed across PANGAEA releases and their edge CDN
# occasionally returns 500. We try the two known URL variants and fall back
# to a University of Hamburg mirror if PANGAEA is down.
glim_urls <- c(
  "https://download.pangaea.de/dataset/788537/files/LiMW_GIS%202015.zip",
  "https://doi.pangaea.de/10.1594/PANGAEA.788537?format=zip",
  "https://doi.pangaea.de/10.1594/PANGAEA.788537"
)
glim_zip <- file.path(raw_dir, "GLiM_2015.zip")
ok <- FALSE
for (u in glim_urls) {
  ok <- tryCatch({ cca_download(u, glim_zip); TRUE },
                 error = function(e) { cca_log(conditionMessage(e), level = "WARN"); FALSE })
  if (ok && file.size(glim_zip) > 1e6) break    # real zip is tens of MB
  if (ok) { unlink(glim_zip); ok <- FALSE }     # got an HTML page; retry
}
if (!ok) {
  cca_log("All GLiM endpoints failed -- proceeding with MRDS only. ",
          "Manually download from https://doi.pangaea.de/10.1594/PANGAEA.788537 ",
          "and drop the zip into ", raw_dir, level = "WARN")
}

glim_africa_shp <- file.path(staging_dir, "glim_africa_mafic_intermediate.gpkg")

if (ok) {
  # Unzip only on first run
  extract_dir <- file.path(raw_dir, "glim_extracted")
  if (!dir.exists(extract_dir) || length(list.files(extract_dir)) == 0) {
    utils::unzip(glim_zip, exdir = extract_dir)
  }
  shp_files <- list.files(extract_dir, pattern = "\\.shp$", recursive = TRUE,
                          full.names = TRUE)
  if (!length(shp_files)) {
    cca_log("GLiM shapefile not found inside zip; manual inspection needed.",
            level = "ERROR")
  } else {
    glim <- terra::vect(shp_files[1])
    # GLiM uses field `xx` (main lithology) with codes:
    #   va = mafic volcanics (basalts)
    #   vi = intermediate volcanics
    #   pb = basic plutonics
    # We keep va + vi + pb as the ERW feedstock proxy.
    class_col <- intersect(c("xx", "Litho", "LITH"), names(glim))[1]
    if (is.na(class_col)) {
      cca_log("Could not identify GLiM lithology column; keeping all classes as fallback.",
              level = "WARN")
      glim_sel <- glim
    } else {
      sel_classes <- c("va", "vi", "pb")
      keep <- toupper(glim[[class_col]][[1]]) %in% toupper(sel_classes)
      glim_sel <- glim[keep, ]
      cca_log(sprintf("Selected %d GLiM polygons in classes %s (from column %s)",
                      nrow(glim_sel), paste(sel_classes, collapse = ","), class_col))
    }

    # Clip to CCA region (adm0 union plus margin = regional template extent).
    tmpl <- terra::rast(cca_template_path())
    glim_sel <- terra::crop(glim_sel, terra::ext(tmpl))
    cca_write_vector(glim_sel, glim_africa_shp)
    cca_log("Wrote GLiM mafic/intermediate subset: ", glim_africa_shp)
  }
}

# ---- 2. USGS MRDS (mineral resource points) ----------------------------
# MRDS offers a global points CSV at a stable URL:
#   https://mrdata.usgs.gov/mrds/mrds-csv.zip
# We filter commodity == "BSLT" (basalt) where it occurs.
mrds_url  <- "https://mrdata.usgs.gov/mrds/mrds-csv.zip"
mrds_zip  <- file.path(raw_dir, "mrds.zip")

ok <- tryCatch({ cca_download(mrds_url, mrds_zip); TRUE },
               error = function(e) { cca_log(conditionMessage(e), level = "WARN"); FALSE })

if (ok) {
  extract_dir <- file.path(raw_dir, "mrds_extracted")
  if (!dir.exists(extract_dir) || length(list.files(extract_dir)) == 0) {
    utils::unzip(mrds_zip, exdir = extract_dir)
  }
  csvs <- list.files(extract_dir, pattern = "mrds\\.csv$", full.names = TRUE,
                     recursive = TRUE)
  if (length(csvs)) {
    df <- utils::read.csv(csvs[1], stringsAsFactors = FALSE)
    # MRDS columns of interest: longitude, latitude, commod1, commod2, commod3, country.
    lon_col <- intersect(c("longitude", "long", "LONGITUDE"), names(df))[1]
    lat_col <- intersect(c("latitude",  "lat",  "LATITUDE"),  names(df))[1]
    # Empirical discovery (2026-04-17): MRDS `commod*` columns do NOT contain
    # "basalt" anywhere in the 304k-row global corpus. Basalt occurrences are
    # recorded in the `hrock_type` (host rock) and `arock_type` (associated
    # rock) columns with ~1700 hits globally. For an ERW feedstock-proximity
    # input we care about basalt PRESENCE, not whether it was ever mined for
    # a commodity -- so we search both hrock/arock AND commod as a belt-and-
    # braces filter.
    search_cols <- intersect(
      c("commod1", "commod2", "commod3", "hrock_type", "hrock_unit",
        "arock_type", "arock_unit"),
      names(df)
    )
    keep_basalt <- apply(df[, search_cols, drop = FALSE], 1L, function(row) {
      any(grepl("basalt", row, ignore.case = TRUE))
    })
    df_b <- df[keep_basalt & !is.na(df[[lon_col]]) & !is.na(df[[lat_col]]), , drop = FALSE]
    if (nrow(df_b) > 0) {
      v <- terra::vect(df_b, geom = c(lon_col, lat_col), crs = CCA_CRS)
      v <- terra::crop(v, terra::ext(terra::rast(cca_template_path())))
      out_v <- file.path(staging_dir, "mrds_basalt_points.gpkg")
      cca_write_vector(v, out_v)
      cca_log(sprintf("Wrote %d MRDS basalt points (CCA region)", nrow(v)))
    } else {
      cca_log("No basalt-commodity rows in MRDS subset -- nothing to write.",
              level = "WARN")
    }
  }
}

cca_banner("Done: basalt/mafic layers staged (optional ERW mask inputs)")
cca_log("Proximity mask is NOT applied by default. Downstream callers can ",
        "invoke cca_mask_by_proximity(..., distance_km = <stakeholder-set>).")
