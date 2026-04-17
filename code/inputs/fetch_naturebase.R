# ------------------------------------------------------------------------------
# code/inputs/fetch_naturebase.R
#
# Downloads and unpacks the Naturebase Natural Climate Solutions (NCS)
# "all pathways" bundle so the pipeline is standalone on any machine.
#
# Source:
#   The Nature Conservancy / Naturebase.  Direct URL:
#     https://naturebase-shared-docs.s3.amazonaws.com/Data-download/all_pathways.zip
#   The zip unpacks into a folder named `all_pathways/` containing:
#     * 16 NCS pathway rasters (*_total.tif, 1-km, tCO2e per pixel per year)
#     * pathway-group rollups (total_crop.tif, total_forest.tif, etc.)
#     * all_pathways_tco2e_adm{0,1,2}.csv admin summary tables
#     * Naturebase_metadata.csv (per-pathway methods + citations)
#     * ReadMe_all_pathways.txt + License.txt
#
# Design notes:
#   * Idempotent: if the adm0 CSV (a cheap sentinel) is already present we
#     skip the 500-ish MB download.
#   * Uses cca_download() for the atomic write + provenance sidecar.
#   * The downstream target directory MATCHES what Step 1 expects:
#         data/naturebase/all_pathways/  (not data/inputs/naturebase/)
#     because Step 1 and the Naturebase adm0 CSV lookup were written against
#     that path. Moving the unpacked files into data/naturebase/ keeps the
#     rest of the pipeline unchanged on a fresh clone.
# ------------------------------------------------------------------------------

# ---- Preamble ------------------------------------------------------------
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

cca_banner("Fetch Naturebase NCS bundle (all_pathways.zip)")

NATUREBASE_URL <- "https://naturebase-shared-docs.s3.amazonaws.com/Data-download/all_pathways.zip"

# Destination is data/naturebase/ (not data/inputs/naturebase/) because the
# rest of the pipeline already references CCA_PATHS$naturebase with that path
# — see config.R. Keeping that invariant avoids touching every consumer.
dest_root <- normalizePath(file.path(CCA_ROOT, "data", "naturebase"), mustWork = FALSE)
pathways_dir <- file.path(dest_root, "all_pathways")
dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)

zip_path <- file.path(dest_root, "all_pathways.zip")
sentinel <- file.path(pathways_dir, "all_pathways_tco2e_adm0.csv")

# Skip-if-present guard (a few of the 16 TIFs + the adm0 CSV all have to
# exist; if anything is missing we re-extract from the cached zip).
expected_tifs <- c(
  "crp_cba_total.tif",   # agroforestry
  "crp_scc_cc.tif",      # cover crops (subcomponent)
  "crp_scc_total.tif",   # soil C in croplands (bundle)
  "grs_scg_total.tif",   # grassland soil carbon
  "grs_reg_total.tif",   # reduced grassland emissions
  "grs_grr_total.tif",   # grassland restoration
  "grs_sba_total.tif"    # savanna agroforestry
)
expected <- c(sentinel,
              file.path(pathways_dir, expected_tifs),
              file.path(pathways_dir, "Naturebase_metadata.csv"),
              file.path(pathways_dir, "ReadMe_all_pathways.txt"))

if (cca_skip_if_exists(expected, tag = "naturebase NCS bundle")) quit(status = 0)

# ---- Download (once) -----------------------------------------------------
# The zip is ~500 MB. cca_download() is atomic (writes .part, renames on
# success) and drops a .source.json sidecar next to the zip for provenance.
cca_download(NATUREBASE_URL, zip_path)

# ---- Unpack --------------------------------------------------------------
# Extract directly into data/naturebase/ so the zip's own top-level
# `all_pathways/` directory lands at the expected path.
# unzip() writes files with their intra-zip paths preserved; on first run
# this creates data/naturebase/all_pathways/<files>. On a re-run where a
# partial extract exists, unzip() overwrites atomically (overwrite = TRUE).
cca_log("Unzipping Naturebase bundle to ", dest_root)
extracted <- utils::unzip(zip_path, exdir = dest_root, overwrite = TRUE)
cca_log(sprintf("Unzipped %d file(s)", length(extracted)))

# ---- Validate ------------------------------------------------------------
missing <- expected[!file.exists(expected)]
if (length(missing)) {
  cca_log(sprintf("Naturebase unpack incomplete, missing %d file(s): %s",
                  length(missing), paste(head(basename(missing), 5), collapse = ", ")),
          level = "ERROR")
  quit(status = 0)
}

# Quick range sanity-check on one pathway to confirm the rasters are not
# zero-filled placeholders. Cropland agroforestry dSOC typically peaks
# around 0-50 tCO2e/pixel/yr at 1 km.
suppressPackageStartupMessages(library(terra))
smp <- terra::rast(file.path(pathways_dir, "crp_cba_total.tif"))
mm  <- terra::minmax(smp, compute = TRUE)
cca_log(sprintf("crp_cba_total.tif sanity: min=%.4f max=%.4f", mm[1], mm[2]))
cca_assert_range(smp, lo = -5, hi = 500, name = "Naturebase crp_cba_total")

cca_banner("Done: Naturebase bundle staged at data/naturebase/all_pathways/")
