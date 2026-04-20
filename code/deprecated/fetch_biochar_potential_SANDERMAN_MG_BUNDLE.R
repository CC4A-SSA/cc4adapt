# ------------------------------------------------------------------------------
# code/deprecated/fetch_biochar_potential_SANDERMAN_MG_BUNDLE.R
#
# *** DEPRECATED 2026-04-18 — DO NOT WIRE INTO THE PIPELINE ***
#
# This script downloaded `scenario_crop_MG_dSOC_Y20.tif` from Sanderman 2020
# (Harvard Dataverse DOI 10.7910/DVN/HA17D3) and labelled it the biochar
# mitigation-potential layer. That file is the *Management bundle* —
# cover crops + crop-residue retention + biochar combined into a single
# delta-SOC surface — NOT a biochar-only layer. Sanderman 2020 publishes no
# biochar-only raster. Using the MG bundle as biochar triple-counts cover
# crops and residue retention, which propagated through the entire CCA
# pipeline (Step 1 mitigation, Step 4 opportunity zones, Activity 3 NEB).
#
# Replacement: code/inputs/build_biochar_potential.R synthesises a
# biochar-only layer from MapSPAM residues and pyrolysis parameters.
# Methodology: docs/methodology_biochar.md (postmortem in section 7).
#
# This file is kept under code/deprecated/ for two reasons:
#   1. Reproducibility of pre-2026-04-18 numbers if anyone needs to rerun
#      the old (incorrect) pipeline for comparison.
#   2. The MG bundle remains a useful upper-bound sanity check on the
#      *combined* (cover + residue + biochar) cropland soil-carbon practices.
#
# It is intentionally NOT in the orchestrator (code/inputs/00_fetch_all.R).
# To run it manually for the upper-bound sanity check, invoke directly:
#   Rscript code/deprecated/fetch_biochar_potential_SANDERMAN_MG_BUNDLE.R
# The output is written to data/staging/sanderman_mg_reference/ (not the
# biochar staging dir) so no downstream consumer can pick it up by mistake.
# ------------------------------------------------------------------------------
#
# Original header preserved below for institutional memory.
# ------------------------------------------------------------------------------
#
# Fetches a global biochar mitigation-potential raster and clips to the CCA
# region (NGA/KEN/ETH + 5-arcmin snapped bbox).
#
# Source (cross-checked against Naturebase_metadata.csv, which cites the same
# dataset as the basis for its biochar subcomponent):
#
#   Sanderman, Jonathan; Woolf, Dominic; Lehmann, Johannes; Rivard, Charlotte;
#   Poggio, Laura; Heuvelink, Gerard; Bossio, Deborah, 2020,
#   "Soils Revealed soil carbon futures",
#   https://doi.org/10.7910/DVN/HA17D3,
#   Harvard Dataverse, V3.
#
# The Dataverse package contains per-scenario dSOC rasters at 20- and 50-year
# horizons. The biochar-specific layer we want is
# `scenario_crop_biochar_dSOC_Y20.tif` (delta SOC in tC/ha over 20 y with full
# biochar adoption). We multiply by 44/12 to get tCO2e/ha, divide by the
# horizon to get tCO2e/ha/yr, and clip to NGA/KEN/ETH.
#
# --- Why this dataset and not Woolf et al. 2021 ---
#
# The project plan originally listed Woolf et al. 2021 (Nature Comms). That
# paper publishes a regional *country-level* inventory model, not a
# reproducible global raster. Sanderman et al. 2020 is the canonical global
# raster and is the same source Naturebase cites in their v2 biochar method
# note ("Naturebase_metadata.csv": "Sanderman et al. 2020 ... Harvard
# Dataverse"). Using the same upstream keeps our biochar numbers comparable
# to the Naturebase crp_scc layer.
#
# --- Why we do NOT apply a suitability mask here ---
#
# Biochar suitability is a fuzzy function of soil pH, OC, clay, and rainfall
# (CropSuiteLite). That is computed in a separate step
# (`build_biochar_suitability.py`) and applied in Step-1 mitigation. This
# script only fetches the raw potential; downstream decisions (suitability
# cut-off, residue availability, feedstock cost cut-off) remain explicit.
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

suppressPackageStartupMessages({
  library(terra)
})

cca_banner("Fetch biochar potential (Sanderman 2020, Harvard Dataverse)")

# ---- Target paths --------------------------------------------------------
# DEPRECATED 2026-04-18: outputs now go to data/staging/sanderman_mg_reference/
# (not data/staging/biochar/) so the live pipeline can never accidentally
# pick this up as the biochar mitigation layer.
raw_dir      <- file.path(CCA_PATHS$inputs_root, "sanderman_mg_reference")
staging_dir  <- file.path(CCA_PATHS$staging_root, "sanderman_mg_reference")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

# Short-circuit on successful prior run
staged_out <- file.path(staging_dir, "sanderman2020_crop_MG_bundle_tco2e_ha_yr.tif")
if (cca_skip_if_exists(staged_out, tag = "biochar_potential")) quit(status = 0)

# Harvard Dataverse file IDs for the Sanderman 2020 package are enumerated in
# the Dataverse metadata endpoint. We keep the filename we want and the DOI;
# the resolver API gives us the file ID at runtime. Storing the file ID here
# directly would break if Dataverse reingests the package.
#
# Empirical discovery (2026-04-17, against the live DOI): Sanderman 2020 does
# NOT publish a standalone biochar layer. The closest published raster is
# `scenario_crop_MG_dSOC_Y20.tif` -- the "management bundle" dSOC surface
# which combines cover crops, residue retention, and biochar at IPCC Tier-1
# emission factors (as described in the Soils Revealed methodology doc).
# Naturebase cites this same file in its crp_scc method note. We therefore:
#   (a) download the MG bundle as the REFERENCE layer for biochar-adjacent
#       soil-carbon change, flagging it clearly as a bundle (not biochar-only);
#   (b) expect downstream construction of a biochar-only layer from
#       MapSPAM residue availability (build_residue_proxy.R) combined with a
#       biochar-specific pyrolysis / mineralisation coefficient in Step 1.
DATAVERSE_DOI <- "10.7910/DVN/HA17D3"
WANTED_FILE   <- "scenario_crop_MG_dSOC_Y20.tif"        # management bundle
HORIZON_YEARS <- 20                                     # from filename suffix
BUNDLE_LABEL  <- "sanderman2020_crop_MG_bundle"         # cover+residue+biochar

# ---- 1. Resolve file ID via Dataverse API --------------------------------
# GET https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:...
# returns a JSON doc listing latestVersion.files[] with each file's label
# and dataFile.id. We parse that to find our wanted filename.
resolve_dataverse_id <- function(doi, wanted_label) {
  api <- sprintf(
    "https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:%s",
    doi
  )
  tmp <- tempfile(fileext = ".json")
  cca_download(api, tmp)
  # `jsonlite` is in Imports of arrow already installed; cheap dependency.
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("[CCA] jsonlite is required for Dataverse API resolution. install.packages('jsonlite').")
  }
  meta <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  files <- meta$data$latestVersion$files
  hit <- NULL
  for (f in files) {
    if (identical(f$label, wanted_label)) { hit <- f; break }
  }
  if (is.null(hit)) {
    labels <- vapply(files, function(f) f$label %||% "", character(1))
    stop(sprintf(
      "[CCA] File '%s' not found in Dataverse DOI %s. Available: %s",
      wanted_label, doi, paste(head(labels, 20), collapse = ", ")
    ))
  }
  list(id = hit$dataFile$id, size = hit$dataFile$filesize)
}

resolve <- tryCatch(resolve_dataverse_id(DATAVERSE_DOI, WANTED_FILE),
                    error = function(e) { cca_log(conditionMessage(e), level = "ERROR"); NULL })

if (is.null(resolve)) {
  cca_log("Could not resolve Sanderman biochar file ID; aborting fetch.",
          level = "ERROR")
  quit(status = 0)  # non-fatal: other pipeline stages can still run
}

# ---- 2. Download raw global TIF ------------------------------------------
raw_path <- file.path(raw_dir, WANTED_FILE)
file_url <- sprintf("https://dataverse.harvard.edu/api/access/datafile/%d", resolve$id)
cca_download(file_url, raw_path)

# ---- 3. Sanity-check raw raster ------------------------------------------
# Sanderman 2020 publishes rasters in Mollweide (EPSG:54009) at ~1 km. We
# DO NOT force the CRS -- terra picks up the projection from the TIFF tags.
# cca_clip_to_cca_region handles the Mollweide-to-4326 transform via a
# template-extent reprojection (cheap; avoids reprojecting the global raster).
raw <- terra::rast(raw_path)
cca_assert_not_allna(raw, name = "Sanderman MG bundle raw")
cca_log(sprintf("Sanderman raw CRS: EPSG:%s, dims %s",
                terra::crs(raw, describe = TRUE)$code,
                paste(dim(raw), collapse = " x ")))
cca_assert_range(raw, lo = -50, hi = 200, name = "Sanderman MG dSOC (tC/ha, Y20)")

# ---- 4. Clip (reproject on the way out), then convert units --------------
# Order matters: reproject-then-divide is equivalent to divide-then-reproject
# for this linear transform, but the reproject-first path produces a ~5 MB
# CCA-region raster instead of muliplying a 488 MB global raster.
clipped_tC <- cca_clip_to_cca_region(raw, cca_template_path(), method = "bilinear")
clipped    <- clipped_tC * (44 / 12) / HORIZON_YEARS   # tC/ha/Y20 -> tCO2e/ha/yr

# ---- 5. Write staging COG ------------------------------------------------
# The file staged here is the MANAGEMENT BUNDLE -- it is NOT biochar-only.
# Downstream Step 1 does the biochar attribution from MapSPAM residues
# (build_residue_proxy.R) x biochar-specific pyrolysis yield. This bundle is
# useful as a reference for cross-checking: biochar's contribution to crp_scc
# should be a fraction of the MG bundle.
out_tif <- file.path(staging_dir, sprintf("%s_tco2e_ha_yr.tif", BUNDLE_LABEL))
cca_write_cog(clipped, out_tif)

# Provenance sidecar documenting the bundle caveat.
writeLines(c(
  "{",
  sprintf('  "source": "%s",', "Sanderman 2020 Harvard Dataverse DOI 10.7910/DVN/HA17D3"),
  sprintf('  "file": "%s",', WANTED_FILE),
  sprintf('  "horizon_years": %d,', HORIZON_YEARS),
  '  "unit_output": "tCO2e/ha/yr (horizon-averaged)",',
  '  "bundle": "cover_crops + residue_retention + biochar (Tier-1)",',
  '  "note": "NOT a biochar-only layer. Biochar-only must be built from MapSPAM residues.",',
  '  "conversion": "dSOC (tC/ha over Y20) * 44/12 / 20"',
  "}"),
  paste0(out_tif, ".method.json")
)

cca_banner("Done: Sanderman MG bundle staged at native ~1 km")
cca_log("This is a CROPLAND MANAGEMENT BUNDLE, not biochar-only. Downstream")
cca_log("Step 1 constructs biochar-only from MapSPAM residues x pyrolysis yield.")
