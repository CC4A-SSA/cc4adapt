# ------------------------------------------------------------------------------
# code/inputs/build_residue_proxy.R
#
# Builds a per-pixel crop-residue biomass proxy used to downweight biochar
# mitigation potential by local feedstock availability.
#
# Why we need this:
#   The raw Sanderman 2020 biochar dSOC raster assumes full residue
#   availability everywhere. In reality, biochar potential is capped by
#   how much crop residue a pixel actually produces. Naturebase's own
#   biochar subcomponent does this weighting too -- their method note reads
#   "grid-cell level estimates of crop production, using MapSPAM, to develop
#   a residue production model of crop residues". We follow the same logic
#   so our biochar layer is comparable to theirs.
#
# Formula (matches the atlas convention in atlas/atlas_exposure/):
#   residue_t_pix = sum_crop( production_t_pix[crop] * residue_index[crop] )
#
# where `residue_index` is the harvest-index-derived residue:grain ratio
# drawn from IPCC Table 11.2 (2019 Refinement) -- same table Naturebase
# cites. We expose it as a CSV so values are auditable.
#
# --- No thresholding here ---
#
# This script produces the residue biomass raster (t/pixel/yr). Whether
# biochar benefits are allowed only above some residue threshold is an
# optional downstream decision via cca_mask_by_threshold.
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

suppressPackageStartupMessages({
  library(terra)
})

cca_banner("Build crop-residue biomass proxy (MapSPAM x IPCC harvest indices)")

staging_dir  <- file.path(CCA_PATHS$staging_root, "biochar")
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

staged_tif <- file.path(staging_dir, "residue_biomass_tonnes_per_pixel.tif")
if (cca_skip_if_exists(staged_tif, tag = "residue_proxy")) quit(status = 0)

# ---- Residue index table -------------------------------------------------
# Source: IPCC 2019 Refinement Vol 4 Ch 11 Table 11.2, residue-to-crop ratio
# Rag = above-ground residue dry matter / grain dry matter. We use the
# SSA-applicable values where available and global means otherwise.
# This table is deliberately explicit so the stakeholder workshop can
# override it by editing a single file.
residue_table <- data.frame(
  mapspam_code = c("maiz","whea","rice","barl","sorg","mill","othc", "pota","cass","yams",
                   "swpo","othr","soyb","grou","cnut","cott","sunf","rape","othe","sugc","bana"),
  crop_label   = c("Maize","Wheat","Rice","Barley","Sorghum","Millet","Other cereals",
                   "Potato","Cassava","Yams","Sweet potato","Other roots","Soybean",
                   "Groundnut","Coconut","Cotton","Sunflower","Rapeseed","Other oilseeds",
                   "Sugarcane","Banana"),
  residue_ratio = c(1.03, 1.61, 1.05, 1.69, 1.24, 1.34, 1.30,
                    0.57, 0.59, 0.40, 0.57, 0.60, 1.50,
                    1.10, 0.30, 2.10, 2.00, 1.60, 1.60, 0.16, 0.35),
  notes = c(
    "IPCC T11.2 Africa",
    "IPCC T11.2 global mean",
    "IPCC T11.2 tropical",
    "IPCC T11.2 global mean",
    "IPCC T11.2 Africa",
    "IPCC T11.2 Africa",
    "fallback: mean of cereals",
    "IPCC T11.2 global",
    "IPCC T11.2 tropical",
    "IPCC T11.2 tropical",
    "IPCC T11.2 tropical",
    "IPCC T11.2 roots fallback",
    "IPCC T11.2 global",
    "IPCC T11.2 Africa",
    "IPCC T11.2 tropical",
    "IPCC T11.2 global",
    "IPCC T11.2 global",
    "IPCC T11.2 temperate",
    "IPCC T11.2 oilseeds mean",
    "IPCC T11.2 (bagasse)",
    "fallback: FAO IPCC"
  ),
  stringsAsFactors = FALSE
)
out_csv <- file.path(staging_dir, "residue_indices_ipcc_t11.2.csv")
utils::write.csv(residue_table, out_csv, row.names = FALSE)
cca_log("Wrote residue index table: ", out_csv)

# ---- Find MapSPAM production rasters ------------------------------------
# Expected location after atlas_exposure pipeline run:
#   <atlas>/atlas_exposure/output/SPAM2020_P_*.tif
# We allow a local override under data/inputs/mapspam/.
mapspam_dir <- file.path(CCA_PATHS$inputs_root, "mapspam")
if (!dir.exists(mapspam_dir)) {
  mapspam_dir <- file.path(CCA_ROOT, "..", "atlas_exposure", "output")
}
if (!dir.exists(mapspam_dir)) {
  cca_log("MapSPAM production rasters not found. Run atlas/atlas_exposure/ or ",
          "drop files into data/inputs/mapspam/. Exiting.", level = "ERROR")
  quit(status = 0)
}

# Per-crop tif pattern: SPAM2020_P_<code>_A.tif  (A = all technologies)
tifs <- list.files(mapspam_dir, pattern = "SPAM2020_P_.*_A\\.tif$",
                   full.names = TRUE, ignore.case = TRUE)
if (!length(tifs)) {
  cca_log("No SPAM2020_P_*_A.tif files in ", mapspam_dir, level = "ERROR")
  quit(status = 0)
}

# ---- Compute residue per pixel ------------------------------------------
template <- terra::rast(cca_template_path())
accum <- terra::rast(template); terra::values(accum) <- 0

for (tif in tifs) {
  code <- regmatches(basename(tif),
                     regexpr("(?<=SPAM2020_P_)[a-z]{4}", basename(tif), perl = TRUE))
  ratio <- residue_table$residue_ratio[residue_table$mapspam_code == code]
  if (!length(ratio)) {
    cca_log("No residue ratio for MapSPAM code '", code, "' -- using 0.", level = "WARN")
    next
  }
  r <- terra::rast(tif)
  r <- terra::project(r, template, method = "sum")
  accum <- accum + (r * ratio)
}
# Values above 0 only; NA outside cropland (implicit from template NA mask)
accum <- terra::ifel(accum > 0, accum, NA)

cca_assert_not_allna(accum, name = "residue biomass proxy")
out_tif <- file.path(staging_dir, "residue_biomass_tonnes_per_pixel.tif")
cca_write_cog(accum, out_tif)

cca_banner("Done: residue biomass proxy built")
