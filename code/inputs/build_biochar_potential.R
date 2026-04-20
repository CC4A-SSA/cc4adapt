# ------------------------------------------------------------------------------
# code/inputs/build_biochar_potential.R
#
# (1) WHAT THIS SCRIPT DOES
# Synthesises a per-pixel biochar mitigation-potential raster (tCO2e/ha/yr)
# for the CCA region (NGA + KEN + ETH at 5 arc-minutes) from already-staged
# crop-residue biomass and pyrolysis parameters. No external HTTP download.
#
# (2) DATA IT READS, AND WHERE EACH INPUT CAME FROM
#   - data/staging/biochar/residue_biomass_tonnes_per_pixel.tif
#       built by code/inputs/build_residue_proxy.R from MapSPAM 2010 V2r0
#       production rasters x IPCC 2019 Refinement Vol 4 Ch 11 Table 11.2
#       residue:product ratios. Units: tonnes dry residue / pixel / yr.
#   - data/staging/suitability/biochar_suitability_5arcmin.tif
#       built by code/inputs/build_biochar_suitability.R from SoilGrids
#       (pH, organic C, clay) + WorldClim rainfall, fuzzy trapezoid combine.
#       Units: dimensionless 0..1.
#   - code/inputs/data/biochar_pyrolysis_params.yaml
#       hand-curated scalar parameters with citations. Read at runtime so
#       Step 5 Monte Carlo can perturb them without code edits.
#   - data/staging/template/cca_template_3country_5arcmin.tif
#       defines the output extent + resolution.
#
# (3) METHODOLOGY / FORMULA  (matches docs/methodology_biochar.md exactly)
#
#   biochar_potential_tco2e_ha_yr(p) =
#         residue_biomass_t_per_pixel(p)
#       * R_available           # default 0.50 (Woolf 2010)
#       * Y_biochar              # default 0.32 (Woolf 2010)
#       * C_biochar              # default 0.74 (Lehmann & Joseph 2015)
#       * P_stable_100yr         # default 0.80 (Wang et al. 2016)
#       * (44 / 12)              # C -> CO2e
#       / cellSize_ha(p)         # tonnes/pixel -> tonnes/ha
#       * Suitability(p)         # 0..1 fuzzy multiplier
#
# Because residue_biomass_t_per_pixel is already an annual quantity, NO
# horizon division is applied. (The previous Sanderman-MG-bundle code
# divided by 20 because the source was a 20-yr cumulative dSOC surface;
# that is a unit bug if carried over here.)
#
# (4) PROVENANCE SIDECAR
#   data/staging/biochar/biochar_potential_tco2e_ha_yr.tif.method.json
#   records every parameter value, the formula, the input filenames, and the
#   per-country cross-check results. We do NOT write a .source.json because
#   this is a synthesised layer, not a downloaded one.
#
# (5) CROSS-CHECK RULE
#   The cropland-area-weighted country mean of the output must fall inside
#   the Woolf 2010 SI Tab S6 sustainable-residue envelope:
#       NGA 0.6-1.4   KEN 0.3-0.9   ETH 0.4-1.1   tCO2e/ha/yr.
#   Inside envelope -> INFO. Outside -> WARN. Never aborts (envelopes are
#   themselves uncertain; WARN is enough to trigger human review).
#
# Replaces the deprecated code/inputs/fetch_biochar_potential.R, which
# downloaded Sanderman 2020's MG management bundle (cover crops + residue +
# biochar combined) and incorrectly labelled it biochar-only. See
# docs/methodology_biochar.md S7 for the postmortem.
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

suppressPackageStartupMessages({
  library(terra)
  library(yaml)
})

cca_banner("Build biochar mitigation potential (synthesised from MapSPAM residues)")

# ---- Paths ---------------------------------------------------------------
staging_dir <- file.path(CCA_PATHS$staging_root, "biochar")
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

residue_tif    <- file.path(staging_dir, "residue_biomass_tonnes_per_pixel.tif")
suit_tif       <- file.path(CCA_PATHS$staging_root, "suitability",
                            "biochar_suitability_5arcmin.tif")
params_yaml    <- file.path(CCA_ROOT, "code", "inputs", "data",
                            "biochar_pyrolysis_params.yaml")
out_tif        <- file.path(staging_dir, "biochar_potential_tco2e_ha_yr.tif")
sidecar_json   <- paste0(out_tif, ".method.json")

# Skip if already produced. (CCA_FORCE=1 to override.)
if (cca_skip_if_exists(out_tif, tag = "biochar_potential")) quit(status = 0)

# ---- Sanity-check dependencies exist ------------------------------------
missing <- c(
  if (!file.exists(residue_tif)) residue_tif,
  if (!file.exists(suit_tif))    suit_tif,
  if (!file.exists(params_yaml)) params_yaml
)
if (length(missing)) {
  cca_log("Missing required inputs:", level = "ERROR")
  for (m in missing) cca_log("  ", m, level = "ERROR")
  cca_log("Run build_residue_proxy.R + build_biochar_suitability.R first.",
          level = "ERROR")
  quit(status = 0)  # non-fatal so the orchestrator can continue
}

# ---- Load parameters -----------------------------------------------------
params <- yaml::read_yaml(params_yaml)
R_available    <- as.numeric(params$R_available$value)
Y_biochar      <- as.numeric(params$Y_biochar$value)
C_biochar      <- as.numeric(params$C_biochar$value)
P_stable_100yr <- as.numeric(params$P_stable_100yr$value)

# Defensive parameter sanity-checks. Each is a fraction in [0, 1] except
# Y_biochar which is mass-yield ratio (also bounded [0, 1] for slow pyrolysis).
for (nm in c("R_available", "Y_biochar", "C_biochar", "P_stable_100yr")) {
  v <- get(nm)
  if (!is.finite(v) || v < 0 || v > 1) {
    stop(sprintf("[CCA] Parameter %s out of [0,1]: %s", nm, format(v)))
  }
}

cca_log(sprintf("Parameters: R_avail=%.2f Y_biochar=%.2f C_biochar=%.2f P_stable=%.2f",
                R_available, Y_biochar, C_biochar, P_stable_100yr))

# Combined multiplier: residue (t dry/pixel/yr) -> tCO2e/pixel/yr.
#   biochar_C_yield = R_avail * Y_biochar * C_biochar * P_stable
#   tCO2e/pixel/yr  = residue_t * biochar_C_yield * (44/12)
biochar_c_yield <- R_available * Y_biochar * C_biochar * P_stable_100yr
multiplier      <- biochar_c_yield * (44 / 12)
cca_log(sprintf("Combined biochar-C yield (mass fraction of dry residue) = %.4f",
                biochar_c_yield))
cca_log(sprintf("Multiplier (residue tonnes -> tCO2e tonnes) = %.4f", multiplier))

# ---- Read inputs ---------------------------------------------------------
residue <- terra::rast(residue_tif)        # t dry residue / pixel / yr
suit    <- terra::rast(suit_tif)           # 0..1 dimensionless

cca_assert_crs(residue)
cca_assert_crs(suit)
cca_assert_not_allna(residue, name = "residue biomass proxy")
cca_assert_not_allna(suit,    name = "biochar suitability surface")

# Suitability is built on the 5-arcmin template; the residue proxy was also
# built on the 5-arcmin template. Belt-and-braces: align the two to the same
# grid before multiplication.
if (!terra::compareGeom(residue, suit, stopOnError = FALSE)) {
  cca_log("Residue and suitability grids differ; resampling suitability to residue grid.",
          level = "WARN")
  suit <- terra::resample(suit, residue, method = "bilinear")
}

# ---- Compute biochar potential ------------------------------------------
# Step A: per-pixel tCO2e per year (still totals per pixel).
biochar_per_pixel <- residue * multiplier

# Step B: convert tonnes/pixel/yr to tonnes/ha/yr by dividing by pixel area.
# terra::cellSize(unit="ha") computes the actual area of each cell in
# hectares accounting for latitude (cells shrink toward the poles).
area_ha <- terra::cellSize(biochar_per_pixel, unit = "ha")
biochar_density <- biochar_per_pixel / area_ha

# Step C: apply the fuzzy suitability multiplier (0..1). Pixels with
# unsuitable soils contribute zero; pixels with ideal soils contribute the
# full residue-derived potential.
biochar_density <- biochar_density * suit

# Drop zero-residue pixels back to NA so the layer carries the same
# "no cropland here" mask as the residue proxy.
biochar_density <- terra::ifel(biochar_density > 0, biochar_density, NA)

cca_assert_crs(biochar_density)
cca_assert_res(biochar_density, CCA_RES)
cca_assert_not_allna(biochar_density, name = "biochar potential (synth)")
# Sanity bound: cropland biochar potential is universally < ~10 tCO2e/ha/yr
# (Woolf 2010 Tab S6 max is ~7); we use 20 as a generous safety bound to
# catch unit errors without rejecting outliers.
cca_assert_range(biochar_density, lo = 0, hi = 20,
                 name = "biochar potential (tCO2e/ha/yr)")

# ---- Cross-check against Woolf 2010 country envelope --------------------
# We compute the per-country area-weighted mean over pixels where
# biochar_density > 0 (i.e. cropland with non-zero residue + non-zero
# suitability). The envelope is itself uncertain; we WARN rather than ERROR
# if a country falls outside.
envelope <- params$woolf2010_envelope_tco2e_ha_yr
cross_check_results <- list()

for (iso3 in CCA_COUNTRIES) {
  adm0_path <- cca_admin_path(iso3, "adm0")
  if (!file.exists(adm0_path)) {
    cca_log(sprintf("Admin polygon missing for %s; cross-check skipped", iso3),
            level = "WARN")
    next
  }
  adm0   <- cca_read_vector(adm0_path)
  masked <- terra::mask(biochar_density, adm0)
  weights <- terra::cellSize(masked, unit = "ha")
  # weighted mean = sum(value * area) / sum(area), restricted to non-NA cells.
  num <- as.numeric(terra::global(masked * weights, "sum", na.rm = TRUE))
  den <- as.numeric(terra::global(terra::ifel(is.na(masked), NA, weights),
                                  "sum", na.rm = TRUE))
  wmean <- if (is.finite(den) && den > 0) num / den else NA_real_

  env <- envelope[[iso3]]
  inside <- !is.na(wmean) && wmean >= env$lo && wmean <= env$hi
  lvl    <- if (inside) "INFO" else "WARN"
  cca_log(sprintf(
    "Cross-check %s: cropland-area-weighted mean = %.3f tCO2e/ha/yr; Woolf envelope [%.2f, %.2f] -- %s",
    iso3, wmean, env$lo, env$hi, if (inside) "OK" else "OUT OF ENVELOPE"
  ), level = lvl)
  cross_check_results[[iso3]] <- list(
    weighted_mean = wmean,
    envelope_lo   = env$lo,
    envelope_hi   = env$hi,
    inside        = inside
  )
}

# ---- Write COG ----------------------------------------------------------
cca_write_cog(biochar_density, out_tif)

# ---- Provenance sidecar (.method.json) ----------------------------------
# Build the JSON by hand to keep the dependency footprint minimal (no
# jsonlite required for this script). Values are quoted as needed.
quote_str <- function(s) sprintf('"%s"', gsub('"', '\\\\"', s))
fmt_num   <- function(x) if (is.finite(x)) sprintf("%.6f", x) else "null"

cc_lines <- vapply(names(cross_check_results), function(iso3) {
  r <- cross_check_results[[iso3]]
  sprintf('    "%s": { "weighted_mean": %s, "envelope_lo": %s, "envelope_hi": %s, "inside": %s }',
          iso3, fmt_num(r$weighted_mean), fmt_num(r$envelope_lo),
          fmt_num(r$envelope_hi), tolower(as.character(r$inside)))
}, character(1))

writeLines(c(
  "{",
  paste0('  "layer": ', quote_str("biochar_potential_tco2e_ha_yr"), ","),
  paste0('  "produced_by": ', quote_str("code/inputs/build_biochar_potential.R"), ","),
  paste0('  "doc": ', quote_str("docs/methodology_biochar.md"), ","),
  paste0('  "produced_at_utc": ', quote_str(format(Sys.time(), tz = "UTC", usetz = TRUE)), ","),
  '  "inputs": {',
  paste0('    "residue_biomass_t_per_pixel": ', quote_str(residue_tif), ","),
  paste0('    "biochar_suitability_5arcmin":  ', quote_str(suit_tif), ","),
  paste0('    "params_yaml":                  ', quote_str(params_yaml)),
  '  },',
  '  "parameters": {',
  paste0('    "R_available":    ', fmt_num(R_available), ","),
  paste0('    "Y_biochar":      ', fmt_num(Y_biochar), ","),
  paste0('    "C_biochar":      ', fmt_num(C_biochar), ","),
  paste0('    "P_stable_100yr": ', fmt_num(P_stable_100yr), ","),
  paste0('    "C_to_CO2e":      ', fmt_num(44/12)),
  '  },',
  paste0('  "formula": ',
         quote_str("biochar [tCO2e/ha/yr] = residue [t/pixel/yr] * R_available * Y_biochar * C_biochar * P_stable_100yr * (44/12) / cellSize_ha * suitability"),
         ","),
  '  "cross_check_woolf2010": {',
  paste(cc_lines, collapse = ",\n"),
  '  }',
  "}"
), sidecar_json)

cca_banner("Done: biochar potential (synthesised) staged")
cca_log("Wrote: ", out_tif)
cca_log("Wrote: ", sidecar_json)
cca_log("This layer is BIOCHAR-ONLY, derived from MapSPAM residues x ",
        "pyrolysis parameters (Woolf 2010 + Lehmann & Joseph 2015 + Wang 2016).")
cca_log("See docs/methodology_biochar.md for the full derivation.")
