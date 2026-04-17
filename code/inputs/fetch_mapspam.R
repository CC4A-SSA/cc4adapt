# ------------------------------------------------------------------------------
# code/inputs/fetch_mapspam.R
#
# Fetches MapSPAM 2020 V1.0 (Yu et al. 2020) crop distribution rasters from
# the MAPSPAM Harvard Dataverse and stages country-clipped COGs for the CCA
# pipeline.
#
# Why MapSPAM:
#   * The investment document calls MapSPAM the baseline-yield and dominant-
#     crop source for Step 2 (production benefit). The atlas pipeline does
#     the same (atlas/atlas_exposure/).
#   * AWD IPCC-Tier-1 fallback (fetch_cao_awd.R) multiplies by the MapSPAM
#     rice harvested-area raster.
#   * The biochar residue-availability proxy (build_residue_proxy.R) sums
#     MapSPAM PRODUCTION x residue-ratio per crop.
#
# --- Source ---
#
# IFPRI MAPSPAM 2020 V2r0 (current release as of this pipeline):
#   https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/SWPENT
#
# Empirical discovery (2026-04-17): the Dataverse package ships NOT per-crop
# TIFs but four mega-zips, one per variable:
#   spam2020V2r0_global_harvested_area.geotiff.zip   (H)
#   spam2020V2r0_global_physical_area.geotiff.zip    (A)
#   spam2020V2r0_global_production.geotiff.zip       (P)
#   spam2020V2r0_global_yield.geotiff.zip            (Y)
# Each zip unpacks to per-crop, per-technology TIFs with the file-name pattern
#   spam2020V2r0_global_<V>_<CROP>_<TECH>.tif   (TECH in {A, I, R, S, L, H}).
# We therefore:
#   1. download the 3 zips we need (H, P, Y);
#   2. unzip each to a cache directory;
#   3. copy the TIFs for the crops we care about (filtered by regex) into
#      data/inputs/mapspam/;
#   4. crop each to the CCA regional extent and write a COG in
#      data/staging/mapspam/.
#
# --- No thresholding here ---
#
# Everything is written at native ~5 arcmin. Reprojection to the CCA 5-arcmin
# template is handled downstream (Step 2 does the final resample; MapSPAM is
# already at 5 arcmin so this is effectively a no-op).
# ------------------------------------------------------------------------------

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

cca_banner("Fetch MapSPAM 2020 V1.0 (subset of crops x variables for CCA)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "mapspam")
staging_dir <- file.path(CCA_PATHS$staging_root, "mapspam")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

# ---- What to fetch -------------------------------------------------------
# Crops: a conservative list covering the practices in config/practices.yaml
# (rice for AWD; maize / wheat / sorghum / millet / cassava / groundnut for
# agroforestry + cover crops; soybean / cotton for biochar feedstock mix).
MAPSPAM_CROPS <- c("RICE","MAIZ","WHEA","SORG","MILL","BARL",
                   "SOYB","GROU","CASS","YAMS","POTA","SWPO",
                   "COTT","SUGC","BANA","PLNT","OOIL","COFF",
                   "COCO","OCER")

# Variables: H (harvested area), P (production), Y (yield).
# Technology: A (all technologies bundled). Step 2 may later pull _R (rainfed)
# and _I (irrigated) when needed; for now A covers the main use cases.
MAPSPAM_VARS <- c("H", "P", "Y")
MAPSPAM_TECH <- "A"

# Empirical: the MapSPAM 2020 V2r0 Dataverse record
# (DOI 10.7910/DVN/SWPENT) requires a Guestbook-response POST before files
# can be downloaded, which blocks unattended API access. We therefore fall
# back to MapSPAM 2010 V2r0 (DOI 10.7910/DVN/PRFF8V), which is open-access
# and is the version the main atlas pipeline (atlas/atlas_exposure/) uses.
# To upgrade later, update MAPSPAM_DOI + MAPSPAM_VERSION + MAPSPAM_ZIP_FOR_VAR
# together.
MAPSPAM_DOI <- "10.7910/DVN/PRFF8V"
MAPSPAM_VERSION <- "spam2010v2r0"

# Mapping from our short variable tag (H/P/Y) to the zip file name stem
# present in the Dataverse record.
MAPSPAM_ZIP_FOR_VAR <- c(
  H = "harv_area",
  P = "prod",
  Y = "yield"
)

# ---- Skip-if-complete guard --------------------------------------------
expected <- unlist(lapply(MAPSPAM_CROPS, function(cr) {
  file.path(staging_dir,
            sprintf("%s_global_%s_%s_%s_cca.tif", MAPSPAM_VERSION,
                    MAPSPAM_VARS, cr, MAPSPAM_TECH))
}))
if (cca_skip_if_exists(expected, tag = "mapspam (all crops x vars)")) quit(status = 0)

# ---- Resolve Dataverse file IDs ----------------------------------------
# The Dataverse API returns a latestVersion.files[] list with each file's
# `label` and `dataFile.id`. We build the ID map up-front so we can download
# only the subset we need.
resolve_dv_files <- function(doi) {
  api <- sprintf(
    "https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:%s",
    doi
  )
  tmp <- tempfile(fileext = ".json")
  cca_download(api, tmp)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("[CCA] jsonlite required. install.packages('jsonlite').")
  }
  meta <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  files <- meta$data$latestVersion$files
  out <- data.frame(
    label    = vapply(files, function(f) f$label %||% "", character(1)),
    file_id  = vapply(files, function(f) as.integer(f$dataFile$id %||% NA), integer(1)),
    filesize = vapply(files, function(f) as.numeric(f$dataFile$filesize %||% NA), numeric(1)),
    stringsAsFactors = FALSE
  )
  out
}

dv_files <- tryCatch(resolve_dv_files(MAPSPAM_DOI),
                     error = function(e) { cca_log(conditionMessage(e), level = "ERROR"); NULL })
if (is.null(dv_files) || !nrow(dv_files)) {
  cca_log("MAPSPAM DOI did not resolve. The DOI may have changed; update ",
          "MAPSPAM_DOI in this script. Aborting MAPSPAM fetch.",
          level = "ERROR")
  quit(status = 0)
}
cca_log(sprintf("MAPSPAM package has %d files on Dataverse", nrow(dv_files)))

# ---- Download + unzip each per-variable mega-zip -----------------------
template <- terra::rast(cca_template_path())
bbox_cca <- terra::ext(template)

# Each zip is ~500 MB and holds ~1200 TIFs (42 crops x 6 techs x 5 variants).
# We unzip once into a variable-scoped cache, then iterate over the subset we
# actually need.
for (v in MAPSPAM_VARS) {
  zip_stem <- MAPSPAM_ZIP_FOR_VAR[[v]]
  pattern <- sprintf("%s_global_%s\\.geotiff\\.zip$", MAPSPAM_VERSION, zip_stem)
  hits <- dv_files[grepl(pattern, dv_files$label, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(hits)) {
    cca_log(sprintf("MAPSPAM zip label not found for %s (%s)", v, pattern),
            level = "WARN")
    next
  }
  hit <- hits[1, ]
  zip_path <- file.path(raw_dir, hit$label)
  url <- sprintf("https://dataverse.harvard.edu/api/access/datafile/%d", hit$file_id)
  ok <- tryCatch({ cca_download(url, zip_path); TRUE },
                 error = function(e) { cca_log(conditionMessage(e), level = "WARN"); FALSE })
  if (!ok) next

  exdir <- file.path(raw_dir, sprintf("%s_extracted", v))
  if (!dir.exists(exdir) || !length(list.files(exdir, pattern = "\\.tif$"))) {
    cca_log("Unzipping ", zip_path, " -> ", exdir)
    utils::unzip(zip_path, exdir = exdir)
  }

  # The extracted TIFs use the MAPSPAM naming convention
  #   <prefix>_<V>_<CROP>_<TECH>.tif
  # where prefix includes the version/stem. We iterate our wanted crops and
  # copy + CCA-clip the matching files.
  for (cr in MAPSPAM_CROPS) {
    pat <- sprintf("_%s_%s_%s\\.tif$", v, cr, MAPSPAM_TECH)
    tifs <- list.files(exdir, pattern = pat, full.names = TRUE,
                       recursive = TRUE, ignore.case = TRUE)
    if (!length(tifs)) {
      cca_log(sprintf("MAPSPAM crop %s (%s) not in zip (skipping)", cr, v),
              level = "WARN")
      next
    }
    raw_tif <- tifs[1]
    staged_path <- file.path(staging_dir,
                             sprintf("%s_global_%s_%s_%s_cca.tif",
                                     MAPSPAM_VERSION, v, cr, MAPSPAM_TECH))
    if (cca_skip_if_exists(staged_path, tag = sprintf("mapspam %s_%s_%s", v, cr, MAPSPAM_TECH))) next
    r <- terra::rast(raw_tif)
    if (as.character(terra::crs(r, describe = TRUE)$code) != "4326") {
      terra::crs(r) <- "EPSG:4326"
    }
    r_cca <- terra::crop(r, bbox_cca, snap = "out")
    if (isTRUE(all(is.na(terra::minmax(r_cca, compute = TRUE))))) {
      cca_log(sprintf("MAPSPAM %s_%s_%s: no pixels in CCA extent -- skipping staged write",
                      v, cr, MAPSPAM_TECH), level = "WARN")
      next
    }
    cca_write_cog(r_cca, staged_path)
  }
}

# ---- Convenience aliases -----------------------------------------------
# fetch_cao_awd.R expects SPAM2020_rice_H.tif in a predictable path; Step 2
# (step2_production.R) expects SPAM2020_Y_<crop>_A.tif and the residue proxy
# expects SPAM2020_P_<crop>_A.tif. Create aliases from our CCA-clipped COGs
# so those consumers just work.
alias_dir <- file.path(CCA_PATHS$inputs_root, "mapspam")
dir.create(alias_dir, recursive = TRUE, showWarnings = FALSE)

for (v in MAPSPAM_VARS) {
  for (cr in MAPSPAM_CROPS) {
    src <- file.path(staging_dir,
                     sprintf("%s_global_%s_%s_%s_cca.tif", MAPSPAM_VERSION, v, cr, MAPSPAM_TECH))
    if (!file.exists(src)) next
    alias <- file.path(alias_dir, sprintf("SPAM2020_%s_%s_A.tif", v, tolower(cr)))
    if (!file.exists(alias)) {
      file.copy(src, alias, overwrite = FALSE)
      cca_log("Alias: ", alias)
    }
  }
}
# AWD-specific convenience name used by fetch_cao_awd.R
rice_H <- file.path(alias_dir, "SPAM2020_H_rice_A.tif")
rice_H_alias <- file.path(alias_dir, "SPAM2020_rice_H.tif")
if (file.exists(rice_H) && !file.exists(rice_H_alias)) {
  file.copy(rice_H, rice_H_alias, overwrite = FALSE)
  cca_log("Copied rice H alias for AWD fallback: ", rice_H_alias)
}

cca_banner("Done: MapSPAM subset staged")
