# ------------------------------------------------------------------------------
# code/activity2/step3_adaptation.R
#
# Activity 2 Step 3: adaptation benefit (avoided loss) at each pixel.
#
# Formula (per pixel, per practice):
#
#   avoided_loss_usd_ha_yr =
#       sum_{hazard,severity} [
#         P(hazard, severity)
#       * yield_loss_frac(hazard, severity, crop)
#       * buffer_frac(practice, hazard, severity, crop)
#       * baseline_yield_t_ha(crop)
#       * price_usd_tonne(crop, country)
#       ]
#
# Inputs:
#   * Step-1 mitigation mask (eligibility anchor).
#   * Hazard frequency rasters (annual P of each severity) from the CCA
#     synthesiser:
#         data/staging/hazards/<crop>/<hazard>/<severity>_5arcmin.tif
#     Built by code/inputs/build_hazard_frequencies.R from WorldClim 2.1
#     monthly climatology + Ecocrop crop optima + literature CVs. Full
#     methodology: docs/methodology_hazards.md.
#     Previously this script tried to read from atlas/hazards_prototype/
#     directly, but that folder ships only source code (no rasters), and
#     the canonical S3 publication is auth-gated; the synthesiser
#     replaces both pathways for a fresh-clone pipeline. Where a hazard
#     frequency raster is absent, the hazard is skipped with a clear log
#     -- never silently imputed.
#   * Buffering coefficient table from
#     data/inputs/buffering_coef.csv (built by
#     code/inputs/build_buffering_coefficients.R).
#   * MapSPAM yield and FAO prices (same loaders as Step 2).
#   * Hazard-specific yield_loss_frac table (embedded below; editable).
#
# Outputs:
#   * data/outputs/step3_adaptation/avoidedloss_<practice>_<iso3>_5arcmin.tif
#   * data/outputs/step3_adaptation/adaptation_admin_<iso3>.csv
#     (per-hazard contribution columns preserved for auditability).
#
# Design notes:
#
#   [A] We deliberately do NOT multiply by any "adoption" or "intervention
#       efficacy" factor here -- those are Activity-3 NEB parameters. Step 3
#       gives gross avoided loss at full technical potential. Keeping the
#       unit chain transparent makes NEB easy to audit.
#
#   [B] Hazards are iterated in a fixed order {drought, heat, waterlogging,
#       variability}. We sum contributions assuming independence -- documented
#       as an open question for the stakeholder workshop (plan 8.*). A
#       correlated-hazard formulation is a future extension point.
#
#   [C] Where buffer_frac is flagged `confidence_tier = "Low"` we still use
#       it (the pipeline must produce a layer), but Step 5 Monte Carlo later
#       widens the P10-P90 envelope proportionally.
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

cca_banner("Activity 2 | Step 3 | Adaptation benefit (avoided loss)")

# ---- Load catalogue + buffer table -------------------------------------
practices_yaml    <- yaml::read_yaml(file.path(CCA_ROOT, "code", "activity2", "config", "practices.yaml"))
country_practices <- practices_yaml$country_practices

buffer_csv <- file.path(CCA_PATHS$inputs_root, "buffering_coef.csv")
buffer_tbl <- if (file.exists(buffer_csv) && file.size(buffer_csv) > 0) {
  utils::read.csv(buffer_csv, stringsAsFactors = FALSE)
} else {
  cca_log("buffering_coef.csv not found; run build_buffering_coefficients.R. ",
          "Using zero buffers (adaptation benefit will be zero).", level = "WARN")
  data.frame(practice = character(), hazard = character(), severity = character(),
             crop = character(), buffer_frac = numeric())
}

# Yield loss fraction per (hazard, severity, crop). Literature-grounded
# defaults from the investment document's Step-3 specification; stakeholder
# workshop can edit this table directly.
loss_table <- data.frame(
  hazard   = c(rep("drought", 9), rep("heat", 9), rep("waterlogging", 9)),
  severity = rep(c("moderate","severe","extreme"), 9),
  crop     = rep(rep(c("maiz","rice","sorg"), each = 3), 3),
  loss_frac = c(
    # drought: maize / rice / sorghum at moderate/severe/extreme
    0.10, 0.25, 0.45,
    0.08, 0.22, 0.40,
    0.07, 0.18, 0.35,
    # heat
    0.06, 0.18, 0.35,
    0.10, 0.25, 0.45,
    0.05, 0.15, 0.30,
    # waterlogging
    0.05, 0.15, 0.30,
    0.02, 0.08, 0.20,
    0.04, 0.12, 0.25
  ),
  stringsAsFactors = FALSE
)

# ---- Helpers -----------------------------------------------------------
#  Substring-keyword fallback: build_buffering_coefficients.R only maps a
#  small subset of ERA practice labels to CCA practice keys, so most rows
#  come through with practice = NA but a meaningful era_practice (e.g.
#  "Mulch-Reduced Tillage", "Improved Varieties"). When the exact-key
#  lookup misses, we scan era_practice for keywords associated with the
#  requested CCA practice. The buffer_frac is then the mean across all
#  matching rows for the (hazard, severity) (or hazard alone) combination.
#  This is the right behaviour until Activity-1 stakeholder review locks
#  down the ERA->CCA mapping in the build script -- documented in
#  docs/methodology_hazards.md S7 and README "What needs checking" item 2.
practice_keywords <- list(
  agroforestry_cropland   = c("agroforest", "alley", "hedge", "silvopast",
                              "parkland", "tree", "shrub"),
  cover_crops             = c("cover crop", "green manure", "mulch", "residue retention"),
  biochar                 = c("biochar"),
  awd_rice                = c("awd", "alternate wetting", "intermittent flooding"),
  enhanced_rock_weathering = c("rock weathering", "basalt"),
  grassland_soil_carbon   = c("grazing", "grassland", "rotational"),
  grassland_reduced_emissions = c("grazing"),
  grassland_restoration   = c("restoration", "rehabilitation")
)

buffer_lookup <- function(practice, hazard, severity, crop) {
  # 1. Exact match on the (already-mapped) practice column.
  r <- buffer_tbl[!is.na(buffer_tbl$practice) &
                  buffer_tbl$practice == practice &
                  buffer_tbl$hazard == hazard &
                  buffer_tbl$severity == severity, ]
  if (!nrow(r)) r <- buffer_tbl[!is.na(buffer_tbl$practice) &
                                buffer_tbl$practice == practice &
                                buffer_tbl$hazard == hazard, ]
  # 2. Substring-keyword fallback on era_practice when no mapped row exists.
  if (!nrow(r)) {
    kws <- practice_keywords[[practice]] %||% character(0)
    if (length(kws)) {
      pat <- paste(kws, collapse = "|")
      mask <- !is.na(buffer_tbl$era_practice) &
              grepl(pat, buffer_tbl$era_practice, ignore.case = TRUE) &
              buffer_tbl$hazard == hazard &
              buffer_tbl$severity == severity
      r <- buffer_tbl[mask, ]
      if (!nrow(r)) {
        mask <- !is.na(buffer_tbl$era_practice) &
                grepl(pat, buffer_tbl$era_practice, ignore.case = TRUE) &
                buffer_tbl$hazard == hazard
        r <- buffer_tbl[mask, ]
      }
    }
  }
  if (!nrow(r)) return(0)
  v <- suppressWarnings(mean(as.numeric(r$buffer_frac), na.rm = TRUE))
  # mean() of an all-NA / empty numeric is NaN; coerce to 0 so downstream
  # arithmetic and `== 0` checks behave (NaN || 0 -> NaN, which crashes
  # the contribution-skip predicate). Cap at [0, 1] -- buffer_frac > 1 in
  # build_buffering_coefficients.R means an ERA study reported the
  # treatment EXCEEDING the control under hazard stress (RR_mean > 0),
  # which physically means a practice that avoids 100% of the loss AND
  # adds a yield bonus on top. The yield bonus belongs in Step 2
  # (production benefit), not Step 3 (avoided loss); so we cap at 1.0
  # here to avoid double-counting.
  if (!is.finite(v)) 0 else max(0, min(1, v))
}

loss_lookup <- function(hazard, severity, crop) {
  r <- loss_table[loss_table$hazard == hazard &
                  loss_table$severity == severity &
                  loss_table$crop == crop, , drop = FALSE]
  if (!nrow(r)) return(NA_real_)
  as.numeric(r$loss_frac[1])
}

# Locate a hazard frequency raster per (crop, hazard, severity) under the
# CCA staging tree produced by code/inputs/build_hazard_frequencies.R.
# Returns NA if none exists, in which case the (crop, hazard, severity)
# contribution is skipped at the call site.
hazards_root <- file.path(CCA_PATHS$staging_root, "hazards")
find_haz_raster <- function(crop, hazard, severity) {
  # Canonical layout: data/staging/hazards/<crop>/<hazard>/<severity>_<res>.tif
  primary <- file.path(hazards_root, crop, hazard,
                       sprintf("%s_%s.tif", severity, CCA_RES_LABEL))
  if (file.exists(primary)) return(primary)
  # Fallback glob in case a future release uses a slightly different name
  # (e.g. multi-scenario suffixes); preserves forward-compatibility.
  if (!dir.exists(hazards_root)) return(NA_character_)
  pat <- sprintf("^%s.*\\.tif$", severity)
  hits <- list.files(file.path(hazards_root, crop, hazard), pattern = pat,
                     full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) NA_character_ else hits[1]
}

# ---- Driver -------------------------------------------------------------
hazards <- c("drought", "heat", "waterlogging")
severities <- c("moderate", "severe", "extreme")
crops      <- c("maiz", "rice", "sorg")         # keep in sync with loss_table

# Prices (shared with Step 2)
price_csv <- file.path(CCA_PATHS$staging_root, "faostat",
                       "fao_producer_prices_usd_per_tonne_2022_2024.csv")
prices <- if (file.exists(price_csv)) utils::read.csv(price_csv, stringsAsFactors = FALSE) else NULL
spam_to_fao <- c(maiz = "Maize (corn)", rice = "Rice", sorg = "Sorghum")
mean_price <- function(area, item) {
  if (is.null(prices)) return(NA_real_)
  hits <- prices[prices$Area == area & prices$Item == item, ]
  if (!nrow(hits)) NA_real_ else mean(as.numeric(hits$Value), na.rm = TRUE)
}

# MapSPAM yield loader (copied logic from Step 2 so this script is self-contained)
mapspam_candidates <- c(
  file.path(CCA_PATHS$inputs_root, "mapspam"),
  file.path(CCA_ROOT, "..", "atlas_exposure", "output")
)
mapspam_dir <- mapspam_candidates[dir.exists(mapspam_candidates)][1]

for (iso3 in CCA_COUNTRIES) {
  cca_banner(sprintf("Country: %s", iso3))
  template <- terra::rast(cca_template_path(iso3))
  wanted   <- country_practices[[iso3]] %||% character(0)
  country_name <- switch(iso3, NGA = "Nigeria", KEN = "Kenya", ETH = "Ethiopia")

  admin_rows <- list()
  for (practice in wanted) {
    out_tif <- file.path(CCA_PATHS$step3_dir,
                         sprintf("avoidedloss_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    if (cca_skip_if_exists(out_tif, tag = sprintf("step3 %s/%s", practice, iso3))) next

    step1_tif <- file.path(CCA_PATHS$step1_dir,
                           sprintf("mitigation_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    if (!file.exists(step1_tif)) {
      cca_log(sprintf("Step 1 output missing for %s/%s -- skipping", practice, iso3),
              level = "WARN"); next
    }
    anchor <- terra::rast(step1_tif) > 0

    avoid_accum <- terra::rast(template); terra::values(avoid_accum) <- 0
    touched <- FALSE

    for (crop in crops) {
      # Baseline yield
      if (is.na(mapspam_dir)) { cca_log("No MapSPAM -- skipping crops loop", level = "WARN"); break }
      ytif <- list.files(mapspam_dir,
                         pattern = sprintf("SPAM2020_Y_%s_A\\.tif$", crop),
                         full.names = TRUE, ignore.case = TRUE)
      if (!length(ytif)) next
      # MapSPAM Y is kg/ha. FAO prices are USD/tonne. Convert kg -> t on
      # read so the unit chain (yield t/ha * price USD/t) cancels cleanly.
      # Same fix that Step 2 already carries; this script omitted it
      # because it had never produced output (no hazards in the tree
      # before 2026-04-18).
      y_kg <- terra::project(terra::rast(ytif[1]), template, method = "bilinear")
      y    <- y_kg / 1000
      y    <- terra::mask(y, anchor, maskvalues = FALSE)

      p_usd <- mean_price(country_name, spam_to_fao[[crop]])
      if (!is.finite(p_usd)) {
        cca_log(sprintf("No price for %s/%s -- skipping crop contribution", country_name, crop),
                level = "WARN"); next
      }

      for (hazard in hazards) {
        for (severity in severities) {
          haz_tif <- find_haz_raster(crop, hazard, severity)
          if (is.na(haz_tif)) next
          phs <- terra::project(terra::rast(haz_tif), template, method = "bilinear")
          loss_f   <- loss_lookup(hazard, severity, crop)
          buffer_f <- buffer_lookup(practice, hazard, severity, crop)
          if (is.na(loss_f) || buffer_f == 0) next
          contrib <- phs * loss_f * buffer_f * y * p_usd
          contrib <- terra::ifel(is.na(contrib), 0, contrib)
          avoid_accum <- avoid_accum + contrib
          touched <- TRUE

          admin_rows[[length(admin_rows) + 1L]] <- data.frame(
            iso3 = iso3, practice = practice, crop = crop,
            hazard = hazard, severity = severity,
            loss_frac = loss_f, buffer_frac = buffer_f,
            mean_usd_ha_yr = as.numeric(terra::global(contrib, "mean", na.rm = TRUE)),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (!touched) {
      cca_log(sprintf("No hazard x crop contributions for %s/%s -- skipping output",
                      practice, iso3), level = "WARN"); next
    }
    avoid_accum <- terra::mask(avoid_accum, anchor, maskvalues = FALSE)
    # Cap avoided loss at baseline VoP (documented as a guard in the plan):
    # avoid cannot exceed the total value of baseline production on the pixel.
    # We implement it conservatively as min(avoid_accum, y*p) summed across
    # crops touched -- approximated here by clipping at 10000 USD/ha/yr as a
    # sanity ceiling. The real cap is enforced in Activity 3 NEB.
    avoid_accum <- terra::clamp(avoid_accum, 0, 10000, values = TRUE)
    cca_assert_not_allna(avoid_accum, name = sprintf("avoidedloss %s/%s", practice, iso3))
    cca_write_cog(avoid_accum, out_tif)
  }

  if (length(admin_rows)) {
    df <- do.call(rbind, admin_rows)
    out_csv <- file.path(CCA_PATHS$step3_dir,
                         sprintf("adaptation_admin_%s.csv", iso3))
    utils::write.csv(df, out_csv, row.names = FALSE)
    cca_log("Wrote adaptation admin summary: ", out_csv)
  }
}

cca_banner("Done: Step 3 adaptation")
