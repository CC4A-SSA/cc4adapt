# ------------------------------------------------------------------------------
# code/activity2/step3_adaptation.R
#
# Activity 2 Step 3: adaptation benefit (avoided loss) at each pixel.
#
# Formula (per pixel, per practice) -- "same share lost" (F3) split, adopted
# Nov 2026 (issue CC4A-SSA/cc4adapt#1; af_formulations page section 5):
#
#   avoided_loss_usd_ha_yr =
#       sum_{crop,hazard,severity} [
#         P(hazard, severity)
#       * yield_loss_frac(hazard, severity, crop)            L
#       * buffer_frac(practice, hazard, severity, crop)      m (SIGNED)
#       * practice_yield_t_ha(crop)                          A = B * (1 + g)
#       * price_usd_tonne(crop, country)
#       ]
#
#   production_adj_usd_ha_yr =                              (new layer)
#       dvalue (Step 2, = sum_crop (A - B) * price)
#     - sum_crop (A - B) * E[L_crop] * price
#     = sum_crop (A - B) * (1 - E[L_crop]) * price
#   with E[L_crop] = 1 - prod_hazard (1 - sum_severity P * L)   (see [E]).
#
# What changed from the previous version (which the exploration calls F5):
#   * avoided loss used the CONTROL's yield B; it now uses the practice's own
#     yield A = B * (1 + g), with g from the same ERA source and rule as
#     Step 2 (code/common/era_yield_effects.R). Avoided loss rises by (1 + g).
#   * buffer_frac was clamped to [0, 1]; it is now signed (negative = the
#     practice makes bad years worse) with a sanity bound of [-1, 1] that
#     logs a warning when hit (cca_signed_buffer).
#   * the accumulated avoided loss was floored at 0; the floor is removed, so
#     negative avoided loss (maladaptation) now reaches the output. A
#     symmetric +/-10000 USD/ha/yr sanity bound remains.
#   * Step 2 counts the full normal-year gain every year; under F3 the gain
#     is bitten by the same share L as the control. That correction is
#     written here as production_adj_* because this step owns the hazard
#     inputs; Steps 4-5 and Activity 3 prefer it over dvalue when present.
#   The old total over-counted by L (A - B)(1 - m) per stress year; now
#   production + avoided loss equals the stress-year yield gap exactly
#   (single hazard; see [E] for several).
#   Algebra and tests: code/common/adaptation_split.R,
#   tests/testthat/test_adaptation_split.R.
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
#   * ERA_Derived.rda (same loader as Step 2) for the practice gain g.
#   * Step 2 dvalue_<practice>_<iso3>_5arcmin.tif (only for production_adj).
#
# Outputs (file names and units unchanged; one new layer):
#   * data/outputs/step3_adaptation/avoidedloss_<practice>_<iso3>_5arcmin.tif
#     USD/ha/yr; may now be negative.
#   * data/outputs/step3_adaptation/production_adj_<practice>_<iso3>_5arcmin.tif
#     NEW. USD/ha/yr, expected annual production benefit after the hazard
#     bite. Written only when Step 2 dvalue exists and at least one hazard
#     raster was found. If Step 2 is re-run later, re-run Step 3 too
#     (CCA_FORCE=1) or this layer is stale.
#   * data/outputs/step3_adaptation/adaptation_admin_<iso3>.csv
#     (per-hazard contribution columns preserved for auditability; new
#     column `g` = ERA relative yield gain used for A).
#
# Design notes:
#
#   [A] We deliberately do NOT multiply by any "adoption" or "intervention
#       efficacy" factor here -- those are Activity-3 NEB parameters. Step 3
#       gives gross avoided loss at full technical potential. Keeping the
#       unit chain transparent makes NEB easy to audit.
#
#   [B] Hazards are iterated in a fixed order {drought, heat, waterlogging,
#       variability}. We sum avoided-loss contributions over hazards (no
#       joint-year correction) -- documented as an open question for the
#       stakeholder workshop (plan 8.*). Hazard co-occurrence is unresolved
#       project-wide; a correlated-hazard formulation is a future extension
#       point.
#
#   [C] Where buffer_frac is flagged `confidence_tier = "Low"` we still use
#       it (the pipeline must produce a layer), but Step 5 Monte Carlo later
#       widens the P10-P90 envelope proportionally.
#
#   [D] g comes from cca_era_yield_effect(), the function Step 2 uses, so A
#       here and A - B in Step 2 are the same quantity. Where ERA has no
#       evidence for a (practice, crop) pair Step 2 contributes no production
#       for that crop, i.e. g = 0; Step 3 then uses g = 0 (A = B) for that
#       crop too. If ERA_Derived.rda is missing altogether, g = 0 for every
#       crop (avoided loss falls back to the control yield) and a WARN is
#       logged. g is pooled across all ERA seasons (stress and normal); a
#       season-split g is NOT implemented here.
#
#   [E] Production adjustment across hazards. Per hazard, the expected loss
#       share is E[L_h] = sum_s P(h, s) L(h, s) (severity bands are disjoint,
#       see build_hazard_frequencies.R). Hazards are combined as
#       E[L] = 1 - prod_h (1 - E[L_h]), which ASSUMES hazards occur
#       independently and compound multiplicatively within a year. Set
#       SPLIT_HAZARD_COMBINE <- "additive" below to use min(1, sum_h E[L_h])
#       instead, the rule implied by summing avoided loss over hazards in [B]
#       (with it, production + avoided loss equals the expected F3 gap
#       exactly). The two differ only by second-order joint-hazard terms. The
#       bite applies whatever m is, including where buffer_frac = 0.
#
#   [F] MapSPAM yield is used as B. The exploration found it is a mean-year
#       (not normal-year) yield and suggests B_norm = B / (1 - E[L]); that
#       correction is NOT implemented here.
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
source(file.path(common_dir, "era_yield_effects.R"))
source(file.path(common_dir, "adaptation_split.R"))

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
  # mean() of an all-NA / empty numeric is NaN; coerce to 0 (no evidence =
  # no buffering) so downstream arithmetic and `== 0` checks behave.
  # m is SIGNED since Nov 2026: a negative value means the practice makes
  # stress years worse and must reach the output as negative avoided loss.
  # The old [0, 1] clamp is replaced by a [-1, 1] sanity bound that logs a
  # warning (cca_signed_buffer). The old clamp also guarded against counting
  # a yield bonus twice; under the F3 split the bonus sits in production as
  # (A - B)(1 - L) and avoided loss is A L m, which together equal the
  # stress-year gap exactly, so no clamp at 1 is needed for that reason.
  if (!is.finite(v)) return(0)
  cca_signed_buffer(v, label = sprintf("buffer_frac %s/%s/%s", practice, hazard, severity))
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

# ERA yield effect (same source and rule as Step 2) -> practice gain g.
era_y <- cca_load_era_yield(c(
  file.path(CCA_PATHS$inputs_root, "era", "ERA_Derived.rda"),
  normalizePath(file.path(CCA_ROOT, "..", "solutions", "data", "ERA_Derived.rda"),
                mustWork = FALSE)
))
if (is.null(era_y)) {
  cca_log("ERA_Derived.rda not found -- practice gain g = 0 for every crop, so ",
          "avoided loss uses the control yield and no production adjustment is ",
          "made. Run code/inputs/fetch_era.R.", level = "WARN")
}
gain_lookup <- function(practice, crop) {
  eff <- cca_era_yield_effect(era_y, practice, crop)
  if (is.null(eff) || !is.finite(eff$pc)) 0 else eff$pc
}

# How per-hazard expected loss shares are combined for the production
# adjustment; see design note [E]. "independent" or "additive".
SPLIT_HAZARD_COMBINE <- "independent"

# ---- Driver -------------------------------------------------------------
hazards <- c("drought", "heat", "waterlogging")
severities <- c("moderate", "severe", "extreme")
crops      <- c("maiz", "rice", "sorg")         # keep in sync with loss_table

# Prices (shared with Step 2)
price_csv <- file.path(CCA_PATHS$staging_root, "faostat",
                       "fao_producer_prices_usd_per_tonne_2022_2024.csv")
prices <- if (file.exists(price_csv)) utils::read.csv(price_csv, stringsAsFactors = FALSE) else NULL
spam_to_fao <- SPAM_TO_FAO   # shared with Step 2 (code/common/era_yield_effects.R)
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
    padj_tif <- file.path(CCA_PATHS$step3_dir,
                          sprintf("production_adj_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    dvalue_tif <- file.path(CCA_PATHS$step2_dir,
                            sprintf("dvalue_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    # Re-run when production_adj is missing but could be built (Step 2 done).
    skip_paths <- if (file.exists(dvalue_tif)) c(out_tif, padj_tif) else out_tif
    if (cca_skip_if_exists(skip_paths, tag = sprintf("step3 %s/%s", practice, iso3))) next

    step1_tif <- file.path(CCA_PATHS$step1_dir,
                           sprintf("mitigation_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    if (!file.exists(step1_tif)) {
      cca_log(sprintf("Step 1 output missing for %s/%s -- skipping", practice, iso3),
              level = "WARN"); next
    }
    anchor <- terra::rast(step1_tif) > 0

    avoid_accum <- terra::rast(template); terra::values(avoid_accum) <- 0
    touched <- FALSE
    # Production lost to hazards, sum_crop (A - B) E[L_crop] price (USD/ha/yr).
    bite_accum <- terra::rast(template); terra::values(bite_accum) <- 0
    hazard_seen <- FALSE

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

      # Practice yield A = B (1 + g); g as in Step 2 (design note [D]).
      g_c  <- gain_lookup(practice, crop)
      y_pr <- cca_practice_yield(y, g_c)

      el_by_hazard <- list()   # E[L_h] = sum_s P(h, s) L(h, s), per hazard
      for (hazard in hazards) {
        el_h <- NULL
        for (severity in severities) {
          haz_tif <- find_haz_raster(crop, hazard, severity)
          if (is.na(haz_tif)) next
          phs <- terra::project(terra::rast(haz_tif), template, method = "bilinear")
          loss_f <- loss_lookup(hazard, severity, crop)
          if (is.na(loss_f)) next
          # Expected loss share for the production bite: accumulated for every
          # band with a loss fraction, whatever the practice's buffering.
          pl   <- terra::ifel(is.na(phs), 0, phs * loss_f)
          el_h <- if (is.null(el_h)) pl else el_h + pl

          buffer_f <- buffer_lookup(practice, hazard, severity, crop)
          if (buffer_f == 0) next
          # F3 avoided loss on the practice's own yield A; signed.
          contrib <- phs * cca_f3_avoided_loss(y_pr, loss_f, buffer_f) * p_usd
          contrib <- terra::ifel(is.na(contrib), 0, contrib)
          avoid_accum <- avoid_accum + contrib
          touched <- TRUE

          admin_rows[[length(admin_rows) + 1L]] <- data.frame(
            iso3 = iso3, practice = practice, crop = crop,
            hazard = hazard, severity = severity,
            loss_frac = loss_f, buffer_frac = buffer_f, g = g_c,
            mean_usd_ha_yr = as.numeric(terra::global(contrib, "mean", na.rm = TRUE)),
            stringsAsFactors = FALSE
          )
        }
        if (!is.null(el_h)) el_by_hazard[[hazard]] <- el_h
      }

      # Production bite for this crop:
      #   (A - B) - cca_f3_production(A, B, E[L]) = (A - B) E[L].
      if (length(el_by_hazard)) {
        hazard_seen <- TRUE
        if (g_c != 0) {
          el_c <- cca_combine_expected_loss(el_by_hazard, method = SPLIT_HAZARD_COMBINE)
          bite <- (y_pr - y) * el_c * p_usd
          bite <- terra::ifel(is.na(bite), 0, bite)
          bite_accum <- bite_accum + bite
        }
      }
    }

    # ---- Stress-adjusted production layer (design note [E]) ----------------
    # production_adj = dvalue - bite = sum_crop (A - B)(1 - E[L]) price.
    # Crops Step 2 covers but this step does not (no loss_table rows) keep
    # their full dvalue contribution, i.e. E[L] = 0 for them.
    if (hazard_seen && file.exists(dvalue_tif)) {
      dv <- terra::rast(dvalue_tif)
      if (!isTRUE(terra::compareGeom(dv, template, stopOnError = FALSE))) {
        dv <- terra::project(dv, template, method = "bilinear")
      }
      dv   <- terra::mask(dv, anchor, maskvalues = FALSE)
      padj <- dv - bite_accum
      cca_assert_not_allna(padj, name = sprintf("production_adj %s/%s", practice, iso3))
      cca_write_cog(padj, padj_tif)
      cca_log(sprintf("production_adj %s/%s: mean %.2f USD/ha/yr (dvalue mean %.2f)",
                      practice, iso3,
                      as.numeric(terra::global(padj, "mean", na.rm = TRUE)),
                      as.numeric(terra::global(dv,   "mean", na.rm = TRUE))))
    } else if (!file.exists(dvalue_tif)) {
      cca_log(sprintf("No Step 2 dvalue for %s/%s -- production_adj not written",
                      practice, iso3), level = "WARN")
    } else {
      cca_log(sprintf("No hazard rasters for %s/%s -- production_adj not written",
                      practice, iso3), level = "WARN")
    }

    if (!touched) {
      cca_log(sprintf("No hazard x crop contributions for %s/%s -- skipping output",
                      practice, iso3), level = "WARN"); next
    }
    avoid_accum <- terra::mask(avoid_accum, anchor, maskvalues = FALSE)
    # Sanity bound only. The old lower bound of 0 is removed: avoided loss is
    # negative where the practice makes stress years worse (m < 0). The
    # +/-10000 USD/ha/yr magnitude ceiling stands in for the plan's "cannot
    # exceed baseline value of production" guard; the real cap is enforced in
    # Activity 3 NEB. Values hitting it are logged.
    cca_assert_range(avoid_accum, -10000, 10000,
                     name = sprintf("avoidedloss %s/%s", practice, iso3))
    avoid_accum <- terra::clamp(avoid_accum, -10000, 10000, values = TRUE)
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
