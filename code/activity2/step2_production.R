# ------------------------------------------------------------------------------
# code/activity2/step2_production.R
#
# Activity 2 Step 2: production benefit (delta-yield in t/ha and USD/ha/yr)
# for every priority practice x country, restricted to the Step-1
# mitigation-eligible pixels.
#
# Inputs:
#   * data/outputs/step1_mitigation/mitigation_<practice>_<iso3>_5arcmin.tif
#     -- used as the spatial anchor (where the practice is technically feasible).
#   * MapSPAM 2020 crop rasters for (a) dominant-crop detection per pixel and
#     (b) baseline yield t/ha. Expected under data/inputs/mapspam/ or
#     atlas/atlas_exposure/output/ (produced by atlas_exposure pipeline).
#   * ERA practice x outcome response ratios from
#     atlas/solutions/data/ERA_Derived.rda (rich meta-analysis table).
#   * FAO producer prices (USD/tonne) from
#     data/staging/faostat/fao_producer_prices_usd_per_tonne_2022_2024.csv
#   * code/activity2/config/practices.yaml (practice catalogue + CCA<->ERA map)
#
# Outputs:
#   * data/outputs/step2_production/dyield_<practice>_<crop>_<iso3>_5arcmin.tif
#     (percent yield change applied to baseline yield, in t/ha/yr)
#   * data/outputs/step2_production/dvalue_<practice>_<iso3>_5arcmin.tif
#     (dollar value of Delta-yield, USD/ha/yr, summed across crops per pixel)
#   * data/outputs/step2_production/production_admin_<iso3>.csv
#
# Design notes (cross-checked against solutions/R/ERAAnalyze.R docstring and
# atlas_exposure conventions):
#
#   [A] Step 2 runs ONLY on pixels where Step 1 mitigation > 0. This is the
#       "mitigation-eligible anchor" demanded by the investment document.
#       Outside those pixels, production benefit is undefined (NA).
#
#   [B] ERA response ratio (RR) is log(MeanT/MeanC). We convert to percent
#       change `pc = exp(RR) - 1` before applying to MapSPAM baseline yield.
#       This mirrors atlas/solutions/R/ERAAnalyze.R back-transform.
#
#   [C] Per-practice mapping from CCA practice -> ERA practice label lives in
#       a lookup table at the top of the script (editable by analysts without
#       touching code). Unmapped practices log a WARN and are skipped.
#
#   [D] FAO producer prices are joined at the Area x Item level. MapSPAM
#       crop codes (4-letter) are mapped to FAO Item strings in a second
#       lookup at the top. Prices are averaged across the available years
#       per (country, item) because FAOSTAT reporting is uneven.
#
#   [E] No implicit masking. Price gaps produce NA dvalue, never a silently
#       imputed value. Downstream NEB propagates NAs cleanly.
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
  library(yaml)
})

cca_banner("Activity 2 | Step 2 | Production benefit (ERA x MapSPAM x prices)")

# ---- Load catalogue + lookup tables ------------------------------------
practices_yaml <- yaml::read_yaml(
  file.path(CCA_ROOT, "code", "activity2", "config", "practices.yaml")
)
practice_defs     <- practices_yaml$practices
country_practices <- practices_yaml$country_practices

# CCA practice key -> ERA PrName (as they appear in ERA_Derived.rda).
# Compound ERA labels (e.g. "Alleycropping-Intercropping") are handled via
# substring match.
cca_to_era <- list(
  agroforestry_cropland       = c("Agroforestry", "Alleycropping", "Hedgerow", "Parklands", "Silvopasture"),
  cover_crops                 = c("Green Manure", "Cover Crop", "Crop Residue Incorporation"),
  biochar                     = c("Biochar"),
  awd_rice                    = c("Improved Rice Management", "Alternate Wetting and Drying"),
  enhanced_rock_weathering    = c("Silicate", "Enhanced Weathering", "Rock Amendment"),
  grassland_soil_carbon       = c("Grazing Management", "Improved Grazing"),
  grassland_reduced_emissions = c("Grazing Management"),
  grassland_restoration       = c("Grassland Restoration", "Pasture Restoration"),
  savanna_agroforestry        = c("Agroforestry", "Parklands")
)

# MapSPAM 4-letter code -> FAO Item string (as FAOSTAT spells them).
spam_to_fao <- c(
  maiz = "Maize (corn)",
  whea = "Wheat",
  rice = "Rice",
  barl = "Barley",
  sorg = "Sorghum",
  mill = "Millet",
  soyb = "Soya beans",
  grou = "Groundnuts, excluding shelled",
  cnut = "Coconuts, in shell",
  cott = "Cotton lint, ginned",
  sunf = "Sunflower seed",
  rape = "Rape or colza seed",
  pota = "Potatoes",
  cass = "Cassava, fresh",
  yams = "Yams",
  swpo = "Sweet potatoes",
  bana = "Bananas",
  sugc = "Sugar cane"
)

# ---- Prices join -------------------------------------------------------
price_csv <- file.path(CCA_PATHS$staging_root, "faostat",
                       "fao_producer_prices_usd_per_tonne_2022_2024.csv")
prices <- if (file.exists(price_csv)) {
  utils::read.csv(price_csv, stringsAsFactors = FALSE)
} else {
  cca_log("FAO prices CSV missing -- run code/inputs/fetch_fao_prices.R. ",
          "Step 2 will produce dyield rasters but dvalue will be NA.", level = "WARN")
  NULL
}

mean_price <- function(area, item) {
  # Return mean USD/tonne for (area, item) across available years, NA if absent.
  if (is.null(prices)) return(NA_real_)
  hits <- prices[prices$Area == area & prices$Item == item, ]
  if (!nrow(hits)) return(NA_real_)
  mean(as.numeric(hits$Value), na.rm = TRUE)
}

# ---- ERA effect sizes --------------------------------------------------
# Prefer the locally-staged copy so a standalone checkout works without a
# sibling atlas/solutions repo. fetch_era.R populates data/inputs/era/ from
# the sibling or Zenodo.
era_rda_candidates <- c(
  file.path(CCA_PATHS$inputs_root, "era", "ERA_Derived.rda"),
  normalizePath(file.path(CCA_ROOT, "..", "solutions", "data", "ERA_Derived.rda"),
                mustWork = FALSE)
)
era_rda <- era_rda_candidates[file.exists(era_rda_candidates)][1]
if (is.na(era_rda)) {
  cca_log("ERA_Derived.rda not found. Run code/inputs/fetch_era.R or drop the ",
          "file into data/inputs/era/. Step 2 aborted.", level = "ERROR")
  quit(status = 0)
}
env <- new.env(); load(era_rda, envir = env)
era <- NULL
for (nm in ls(env)) {
  x <- get(nm, envir = env)
  if (inherits(x, "data.frame") || inherits(x, "data.table")) { era <- x; break }
}
cca_log(sprintf("ERA table: %d rows", nrow(era)))

# Derive aggregate RR per (practice family, crop) using only yield outcomes.
yield_mask <- grepl("yield|productivity", era$Out.SubInd, ignore.case = TRUE)
era_y <- era[yield_mask &
             !is.na(era$MeanT) & !is.na(era$MeanC) &
             era$MeanT > 0 & era$MeanC > 0, ,drop = FALSE]
era_y$RR <- log(era_y$MeanT / era_y$MeanC)

era_effect <- function(cca_practice, spam_code) {
  # Match CCA practice -> ERA labels; crop -> ERA Product via SPAM->FAO
  # (close-enough substring match on Product). Return a list with pc
  # (percent change), n_obs, ci_low, ci_high.
  era_labels <- cca_to_era[[cca_practice]]
  if (is.null(era_labels)) return(NULL)
  fao_name <- spam_to_fao[spam_code]
  if (is.na(fao_name)) fao_name <- spam_code
  rows <- era_y[vapply(era_y$PrName, function(x) {
                  any(vapply(era_labels, function(lbl) grepl(lbl, x, fixed = TRUE),
                             logical(1)))
                }, logical(1)) &
                (grepl(spam_code, era_y$Product, ignore.case = TRUE) |
                 grepl(fao_name,   era_y$Product, ignore.case = TRUE)), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  mu  <- mean(rows$RR, na.rm = TRUE)
  sdv <- stats::sd(rows$RR, na.rm = TRUE)
  n   <- sum(!is.na(rows$RR))
  ci  <- 1.96 * (sdv / sqrt(max(n, 1)))
  list(
    pc = exp(mu) - 1,
    n_obs = n,
    ci_low  = exp(mu - ci) - 1,
    ci_high = exp(mu + ci) - 1
  )
}

# ---- MapSPAM inputs ----------------------------------------------------
mapspam_candidates <- c(
  file.path(CCA_PATHS$inputs_root, "mapspam"),
  file.path(CCA_ROOT, "..", "atlas_exposure", "output")
)
mapspam_dir <- mapspam_candidates[dir.exists(mapspam_candidates)][1]
if (is.na(mapspam_dir)) {
  cca_log("MapSPAM rasters not found. Run atlas/atlas_exposure/ or drop into ",
          file.path(CCA_PATHS$inputs_root, "mapspam"), ". Step 2 aborted.",
          level = "ERROR")
  quit(status = 0)
}
# Yield rasters: file name `SPAM2020_Y_<code>_A.tif` (Y=yield t/ha, A=all tech)
yield_tifs <- list.files(mapspam_dir, pattern = "SPAM2020_Y_.*_A\\.tif$",
                         full.names = TRUE, ignore.case = TRUE)

# ---- Driver --------------------------------------------------------------
for (iso3 in CCA_COUNTRIES) {
  cca_banner(sprintf("Country: %s", iso3))
  template <- terra::rast(cca_template_path(iso3))
  adm0     <- cca_read_vector(cca_admin_path(iso3, "adm0"))
  wanted   <- country_practices[[iso3]] %||% character(0)

  admin_rows <- list()
  for (practice in wanted) {
    step1_tif <- file.path(CCA_PATHS$step1_dir,
                           sprintf("mitigation_%s_%s_%s.tif",
                                   practice, iso3, CCA_RES_LABEL))
    if (!file.exists(step1_tif)) {
      cca_log(sprintf("Step 1 output missing for %s/%s -- skipping Step 2",
                      practice, iso3), level = "WARN")
      next
    }

    # Short-circuit
    dvalue_tif <- file.path(CCA_PATHS$step2_dir,
                            sprintf("dvalue_%s_%s_%s.tif",
                                    practice, iso3, CCA_RES_LABEL))
    if (cca_skip_if_exists(dvalue_tif, tag = sprintf("step2 %s/%s", practice, iso3))) next

    cca_log(sprintf("Computing production benefit for %s/%s", practice, iso3))
    anchor <- terra::rast(step1_tif) > 0    # mitigation-eligible mask

    dvalue_accum <- terra::rast(template); terra::values(dvalue_accum) <- 0
    touched_any  <- FALSE

    for (ytif in yield_tifs) {
      code <- regmatches(basename(ytif),
                         regexpr("(?<=SPAM2020_Y_)[a-z]{4}", basename(ytif), perl = TRUE))
      if (!length(code)) next

      eff <- era_effect(practice, code)
      if (is.null(eff)) next       # no ERA evidence for this crop x practice

      # Reproject yield raster to CCA template and mask to anchor pixels.
      # MapSPAM 2010 V2r0 "Y" rasters are in kg/ha; convert to t/ha on read so
      # downstream arithmetic with FAO prices (USD/tonne) is unit-consistent.
      y_kg_ha <- terra::project(terra::rast(ytif), template, method = "bilinear")
      y <- y_kg_ha / 1000          # t/ha
      y <- terra::mask(y, anchor, maskvalues = FALSE)

      dyield <- y * eff$pc         # t/ha/yr (can be negative for penalties)

      # Dollar value (USD/ha/yr) via FAO price. Both sides are now per-tonne.
      country_name <- switch(iso3, NGA = "Nigeria", KEN = "Kenya", ETH = "Ethiopia")
      p_usd <- mean_price(country_name, spam_to_fao[code] %||% code)
      dval  <- if (is.finite(p_usd)) dyield * p_usd else dyield * NA

      # Write per-crop Delta-yield COG for auditability.
      out_dyield <- file.path(CCA_PATHS$step2_dir,
                              sprintf("dyield_%s_%s_%s_%s.tif",
                                      practice, code, iso3, CCA_RES_LABEL))
      cca_write_cog(dyield, out_dyield)

      # Accumulate value across crops (USD/ha/yr).
      add <- terra::ifel(is.na(dval), 0, dval)
      dvalue_accum <- dvalue_accum + add
      touched_any <- TRUE

      admin_rows[[length(admin_rows) + 1L]] <- data.frame(
        iso3 = iso3, practice = practice, crop = code,
        dyield_t_ha_mean = as.numeric(terra::global(dyield, "mean", na.rm = TRUE)),
        dvalue_usd_ha_mean = as.numeric(terra::global(dval,  "mean", na.rm = TRUE)),
        n_obs  = eff$n_obs,
        ci_low = eff$ci_low, ci_high = eff$ci_high,
        price_usd_tonne = p_usd,
        stringsAsFactors = FALSE
      )
    }

    if (!touched_any) {
      cca_log(sprintf("No overlap between ERA evidence and MapSPAM crops for %s/%s -- no output",
                      practice, iso3), level = "WARN")
      next
    }

    # Clip Delta-value to the anchor and write the summed raster.
    dvalue_accum <- terra::mask(dvalue_accum, anchor, maskvalues = FALSE)
    cca_assert_not_allna(dvalue_accum, name = sprintf("dvalue %s/%s", practice, iso3))
    cca_write_cog(dvalue_accum, dvalue_tif)
  }

  if (length(admin_rows)) {
    df <- do.call(rbind, admin_rows)
    out_csv <- file.path(CCA_PATHS$step2_dir,
                         sprintf("production_admin_%s.csv", iso3))
    utils::write.csv(df, out_csv, row.names = FALSE)
    cca_log("Wrote production admin summary: ", out_csv)
  }
}

cca_banner("Done: Step 2 production")
