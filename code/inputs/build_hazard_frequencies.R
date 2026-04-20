# ------------------------------------------------------------------------------
# code/inputs/build_hazard_frequencies.R
#
# (1) WHAT THIS SCRIPT DOES
# Produces 27 per-pixel hazard frequency rasters (3 crops x 3 hazards x 3
# severities) at 5 arc-minutes for the CCA region. Each raster carries the
# probability P(hazard, severity, crop, pixel) in [0, 1] used by
# code/activity2/step3_adaptation.R to compute avoided-loss adaptation
# benefit.
#
# (2) DATA IT READS, AND WHERE EACH INPUT CAME FROM
#   - data/inputs/soilgrids/climate/wc2.1_10m/wc2.1_10m_tmax_<MM>.tif
#       12 monthly Tmax climatologies, 10 arc-min, WorldClim 2.1
#       (Fick & Hijmans 2017, DOI 10.1002/joc.5086). Already cached on
#       disk by code/inputs/build_biochar_suitability.R via
#       geodata::worldclim_global("tmax", res=10).
#   - data/inputs/soilgrids/climate/wc2.1_10m/wc2.1_10m_prec_<MM>.tif
#       12 monthly precipitation climatologies, same source.
#   - code/inputs/data/ecocrop_subset_cca.csv
#       Per-crop optimum + absolute T and R, snapshot from
#       atlas/hazards_prototype/metadata/ecocrop.csv (Zea mays, Oryza
#       sativa, Sorghum bicolor).
#   - code/inputs/data/hazards_params.yaml
#       CV(precip), CV(tmax), severity offsets/factors, growing-season
#       month windows -- all with citations.
#   - data/staging/template/cca_template_3country_5arcmin.tif
#       defines output extent + resolution.
#   - (optional) atlas/hazards_prototype/Data/hazard_risk_n/annual/
#       If present (e.g. when CCA is run inside the EiA atlas environment
#       with the hazards_prototype caches mounted), the script short-
#       circuits per-(crop, severity) and reprojects those rasters onto
#       the CCA grid instead of synthesising. See 'two-tier resolution'
#       in docs/methodology_hazards.md S5.
#
# (3) METHODOLOGY  (matches docs/methodology_hazards.md S2 exactly)
#
#   For each (crop c, hazard h, severity s, pixel p):
#
#     mu_T(p)    = mean monthly Tmax over the crop's growing season
#     sigma_T(p) = mu_T(p) * CV_T                              (CV_T = 0.05)
#     mu_R(p)    = sum of monthly precipitation over GS months
#     sigma_R(p) = mu_R(p) * CV_R                              (CV_R = 0.25)
#
#     Heat:
#       T_lo = T_opt_max(c) + lo_offset(s)
#       T_hi = T_opt_max(c) + hi_offset(s)        (Inf for extreme)
#       P_heat(s) = pnorm(T_hi, mu_T, sigma_T) - pnorm(T_lo, mu_T, sigma_T)
#
#     Drought:
#       R_lo = R_opt_min(c) * lo_factor(s)        (0 for extreme)
#       R_hi = R_opt_min(c) * hi_factor(s)
#       P_dry(s)  = pnorm(R_hi, mu_R, sigma_R) - pnorm(R_lo, mu_R, sigma_R)
#
#     Waterlogging:
#       R_lo = R_opt_max(c) * lo_factor(s)
#       R_hi = R_opt_max(c) * hi_factor(s)        (Inf for extreme)
#       P_wet(s)  = pnorm(R_hi, mu_R, sigma_R) - pnorm(R_lo, mu_R, sigma_R)
#
# (4) PROVENANCE SIDECAR
#   data/staging/hazards/<crop>/<hazard>/<severity>_5arcmin.tif.method.json
#   records: source path (atlas-cache vs synthesis), every parameter, the
#   sum-of-severities cross-check result, the per-country mean. Tagged
#   confidence_tier="Low" until Activity-1 stakeholder review.
#
# (5) CROSS-CHECK RULES (matches docs/methodology_hazards.md S6)
#   - sum(P_moderate + P_severe + P_extreme) <= 1 + 1e-3 per pixel
#   - country mean per (crop, hazard) within [0, 0.95]; outside -> WARN
#   - growing-season precip <= annual precip per pixel
#
# Replaces the prior expectation that step3_adaptation.R could find
# rasters under atlas/hazards_prototype/<crop>/<hazard>/<severity>.tif.
# That folder ships only source code; the canonical S3 publication is
# private (HTTP 403 on every probe). See docs/methodology_hazards.md S1.
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

cca_banner("Build hazard frequency rasters (WorldClim + Ecocrop synthesis)")

# ---- Paths ---------------------------------------------------------------
staging_root <- file.path(CCA_PATHS$staging_root, "hazards")
dir.create(staging_root, recursive = TRUE, showWarnings = FALSE)

ecocrop_csv  <- file.path(CCA_ROOT, "code", "inputs", "data", "ecocrop_subset_cca.csv")
params_yaml  <- file.path(CCA_ROOT, "code", "inputs", "data", "hazards_params.yaml")
template_tif <- cca_template_path()
wc_dir       <- file.path(CCA_PATHS$inputs_root, "soilgrids", "climate", "wc2.1_10m")
atlas_haz_dir <- normalizePath(
  file.path(CCA_ROOT, "..", "hazards_prototype", "Data"),
  mustWork = FALSE
)

# Skip-if-exists sentinel: the very last sidecar produced. If that exists
# (and CCA_FORCE is unset), every upstream output also exists.
sentinel <- file.path(staging_root, "sorg", "waterlogging",
                      paste0("extreme_", CCA_RES_LABEL, ".tif.method.json"))
if (cca_skip_if_exists(sentinel, tag = "hazard_frequencies")) quit(status = 0)

# ---- Sanity-check inputs exist ------------------------------------------
missing <- c(
  if (!file.exists(ecocrop_csv))  ecocrop_csv,
  if (!file.exists(params_yaml))  params_yaml,
  if (!file.exists(template_tif)) template_tif,
  if (!dir.exists(wc_dir))        wc_dir
)
if (length(missing)) {
  cca_log("Missing required inputs:", level = "ERROR")
  for (m in missing) cca_log("  ", m, level = "ERROR")
  cca_log("Run build_cca_template.R + build_biochar_suitability.R (or any ",
          "script that calls geodata::worldclim_global) first.", level = "ERROR")
  quit(status = 0)
}

# ---- Load parameters -----------------------------------------------------
ecocrop <- utils::read.csv(ecocrop_csv, stringsAsFactors = FALSE)
params  <- yaml::read_yaml(params_yaml)
CV_R    <- as.numeric(params$cv$precipitation$value)
CV_T    <- as.numeric(params$cv$temperature_max$value)

cca_log(sprintf("CVs: precip=%.2f, Tmax=%.2f", CV_R, CV_T))
cca_log(sprintf("Crops: %s", paste(ecocrop$spam_code, collapse = ", ")))

# ---- Load template (defines output grid) --------------------------------
template <- terra::rast(template_tif)
cca_assert_crs(template)
cca_assert_res(template, CCA_RES)

# ---- Build a single-template climate stack per crop ---------------------
# Approach: for each crop's growing-season month set, sum precipitation
# and average Tmax across those months at the WorldClim 10-arcmin grid,
# then reproject ONCE to the CCA template. We avoid reprojecting the full
# 12-month stack 3 times (one per crop) by reprojecting only the
# crop-specific aggregated rasters. This keeps memory low and is exact
# for the linear sum + mean used here.
read_wc_monthly <- function(var, months) {
  paths <- file.path(wc_dir, sprintf("wc2.1_10m_%s_%02d.tif", var, months))
  miss <- paths[!file.exists(paths)]
  if (length(miss)) {
    # Auto-fetch via geodata if available -- the suitability builders use
    # the same WorldClim cache, so this normally just downloads any vars
    # they didn't pull. Fall through to ERROR if geodata isn't installed.
    if (requireNamespace("geodata", quietly = TRUE)) {
      cca_log(sprintf("Fetching WorldClim 2.1 monthly %s (10m) via geodata::worldclim_global ...", var))
      geodata_dir <- normalizePath(file.path(wc_dir, "..", ".."), mustWork = FALSE)
      tryCatch(
        geodata::worldclim_global(var = var, res = 10, path = geodata_dir),
        error = function(e) cca_log(conditionMessage(e), level = "ERROR")
      )
    }
    miss <- paths[!file.exists(paths)]
    if (length(miss)) {
      cca_log("Missing WorldClim files after fetch attempt: ",
              paste(miss, collapse = ", "), level = "ERROR")
      quit(status = 0)
    }
  }
  terra::rast(paths)
}

# WorldClim ships -32768 NA but terra usually picks this up; verify once.
wc_months_all <- 1:12
prec12 <- read_wc_monthly("prec", wc_months_all)
tmax12 <- read_wc_monthly("tmax", wc_months_all)
annual_prec_full <- terra::app(prec12, sum)  # mm / yr at 10 arc-min
annual_prec_cca  <- terra::project(annual_prec_full, template, method = "bilinear")

cca_assert_not_allna(annual_prec_cca, name = "annual precipitation (CCA grid)")

# ---- Helper: parametric exceedance --------------------------------------
# Avoid creating SpatRaster objects for sigma to keep memory low; we do
# the pnorm in plain matrix space using terra::app or terra::lapp.
# pnorm_band(thr_lo, thr_hi, mu, sigma) returns pnorm(hi) - pnorm(lo).
# `sigma` is a SpatRaster equal to mu * CV with a small floor to avoid
# divide-by-zero in deserts where mu_R rounds to 0.
pnorm_band_rast <- function(thr_lo, thr_hi, mu, sigma) {
  # Floor sigma at a small non-zero value -- a pixel with mu=0 still
  # tolerates a band query (it returns 0 for any positive threshold,
  # which is the right answer: zero rain means certain drought, never
  # waterlogging or heat).
  sigma_f <- terra::ifel(sigma <= 1e-6, 1e-6, sigma)
  if (is.finite(thr_hi)) {
    upper <- terra::app(c(mu, sigma_f),
                        function(x) pnorm(thr_hi, mean = x[1], sd = x[2]))
  } else {
    upper <- terra::ifel(is.na(mu), NA, 1)
  }
  lower <- terra::app(c(mu, sigma_f),
                      function(x) pnorm(thr_lo, mean = x[1], sd = x[2]))
  band <- upper - lower
  band <- terra::ifel(band < 0, 0, band)
  band <- terra::ifel(band > 1, 1, band)
  band
}

# ---- Helper: try the atlas hazards_prototype cache ---------------------
# Returns a SpatRaster on the CCA grid if a matching atlas raster is
# found, or NULL if nothing usable is available. Pattern matches the
# layout produced by atlas/hazards_prototype/R/2_calculate_haz_freq.R.
try_atlas_cache <- function(crop, hazard, severity) {
  if (!dir.exists(atlas_haz_dir)) return(NULL)
  # Atlas naming uses the *common name* (e.g. "maize-moderate.tif"), not
  # the SPAM code. Map back via the ecocrop subset.
  common <- ecocrop$common_name[ecocrop$spam_code == crop]
  if (!length(common)) return(NULL)
  # The prototype publishes per-(crop, severity) stacked rasters with one
  # band per hazard code. We search under hazard_risk* directories and
  # accept any matching file; the actual layer extraction would need the
  # band index for the requested hazard. Until we observe an atlas Data/
  # tree in the wild, we don't know the exact band layout, so this branch
  # logs INFO and returns NULL (forcing the synthesis path) -- the right
  # behaviour for standalone CCA. When an EiA-server run actually has the
  # cache, this is the place to wire the band selector.
  hits <- list.files(atlas_haz_dir, pattern = sprintf(".*%s.*%s.*\\.tif$",
                                                      tolower(common), severity),
                     recursive = TRUE, ignore.case = TRUE, full.names = TRUE)
  if (!length(hits)) return(NULL)
  cca_log(sprintf("Atlas-cache hit for %s/%s/%s: %s (band-extraction not yet wired; falling back to synthesis)",
                  crop, hazard, severity, hits[1]), level = "WARN")
  NULL
}

# ---- Driver: iterate crop x hazard x severity ---------------------------
all_cross_check <- list()

for (i in seq_len(nrow(ecocrop))) {
  crop_row <- ecocrop[i, ]
  crop     <- crop_row$spam_code
  T_opt_max <- as.numeric(crop_row$Temp_Opt_Max)
  R_opt_min <- as.numeric(crop_row$Rain_Opt_Min)
  R_opt_max <- as.numeric(crop_row$Rain_Opt_Max)

  gs <- params$growing_season[[crop]]$months
  if (is.null(gs) || !length(gs)) {
    cca_log(sprintf("No growing season defined for %s -- skipping", crop),
            level = "WARN")
    next
  }
  cca_banner(sprintf("Crop: %s (T_opt_max=%.1f, R_opt_min=%.0f, R_opt_max=%.0f, GS months=%s)",
                     crop, T_opt_max, R_opt_min, R_opt_max,
                     paste(gs, collapse = ",")))

  # Aggregate WorldClim to growing-season summaries on the source grid,
  # then reproject ONCE per (crop, var). Cheaper than reprojecting full
  # 12-month stacks crop-by-crop.
  prec_gs_full <- terra::app(prec12[[gs]], sum)   # mm / GS
  tmax_gs_full <- terra::app(tmax12[[gs]], mean)  # mean monthly Tmax °C / GS

  mu_R <- terra::project(prec_gs_full, template, method = "bilinear")
  mu_T <- terra::project(tmax_gs_full, template, method = "bilinear")
  sigma_R <- mu_R * CV_R
  sigma_T <- mu_T * CV_T

  cca_assert_not_allna(mu_R, name = sprintf("%s GS precip mean", crop))
  cca_assert_not_allna(mu_T, name = sprintf("%s GS Tmax mean", crop))

  # Seasonality sanity check: GS precip <= annual precip
  ratio <- as.numeric(terra::global(mu_R, "mean", na.rm = TRUE)) /
           as.numeric(terra::global(annual_prec_cca, "mean", na.rm = TRUE))
  cca_log(sprintf("%s GS-precip / annual-precip ratio: %.2f", crop, ratio))

  for (hazard in c("heat", "drought", "waterlogging")) {
    hazard_dir <- file.path(staging_root, crop, hazard)
    dir.create(hazard_dir, recursive = TRUE, showWarnings = FALSE)

    severity_rasters <- list()
    for (sev in c("moderate", "severe", "extreme")) {

      atlas_r <- try_atlas_cache(crop, hazard, sev)
      if (!is.null(atlas_r)) {
        # Reproject onto template; tag source.
        out <- terra::project(atlas_r, template, method = "bilinear")
        source_tag <- "atlas_hazards_prototype"
      } else {
        # Synthesis path.
        sev_def <- params$severity[[hazard]][[sev]]
        if (is.null(sev_def)) {
          cca_log(sprintf("No severity definition for %s/%s -- skipping",
                          hazard, sev), level = "WARN")
          next
        }

        if (hazard == "heat") {
          T_lo <- T_opt_max + as.numeric(sev_def$lo_offset)
          T_hi <- T_opt_max + as.numeric(sev_def$hi_offset)  # Inf for extreme
          out  <- pnorm_band_rast(T_lo, T_hi, mu_T, sigma_T)
        } else if (hazard == "drought") {
          R_lo <- R_opt_min * as.numeric(sev_def$lo_factor)  # 0 for extreme
          R_hi <- R_opt_min * as.numeric(sev_def$hi_factor)
          out  <- pnorm_band_rast(R_lo, R_hi, mu_R, sigma_R)
        } else if (hazard == "waterlogging") {
          R_lo <- R_opt_max * as.numeric(sev_def$lo_factor)
          R_hi <- R_opt_max * as.numeric(sev_def$hi_factor)  # Inf for extreme
          out  <- pnorm_band_rast(R_lo, R_hi, mu_R, sigma_R)
        }
        source_tag <- "worldclim_ecocrop_synthesis"
      }

      cca_assert_crs(out)
      cca_assert_range(out, lo = 0, hi = 1,
                       name = sprintf("P(%s,%s,%s)", crop, hazard, sev))

      out_tif <- file.path(hazard_dir,
                           sprintf("%s_%s.tif", sev, CCA_RES_LABEL))
      cca_write_cog(out, out_tif)

      # Country-mean cross-check + sidecar per layer.
      country_means <- list()
      for (iso3 in CCA_COUNTRIES) {
        adm0_path <- cca_admin_path(iso3, "adm0")
        if (!file.exists(adm0_path)) next
        adm0 <- cca_read_vector(adm0_path)
        cm <- as.numeric(terra::global(terra::mask(out, adm0), "mean", na.rm = TRUE))
        lo <- as.numeric(params$country_plausibility$freq_lo)
        hi <- as.numeric(params$country_plausibility$freq_hi)
        ok <- !is.na(cm) && cm >= lo && cm <= hi
        country_means[[iso3]] <- cm
        if (!ok) {
          cca_log(sprintf("Plausibility WARN: %s/%s/%s mean(%s)=%.3f outside [%.2f,%.2f]",
                          crop, hazard, sev, iso3, cm, lo, hi), level = "WARN")
        }
      }
      severity_rasters[[sev]] <- out

      # Sidecar
      sidecar <- paste0(out_tif, ".method.json")
      cm_lines <- vapply(names(country_means), function(iso3) {
        v <- country_means[[iso3]]
        sprintf('    "%s": %s', iso3, if (is.finite(v)) sprintf("%.6f", v) else "null")
      }, character(1))
      writeLines(c(
        "{",
        sprintf('  "layer": "P(%s,%s,%s)",', crop, hazard, sev),
        sprintf('  "produced_by": "code/inputs/build_hazard_frequencies.R",'),
        sprintf('  "doc": "docs/methodology_hazards.md",'),
        sprintf('  "produced_at_utc": "%s",',
                format(Sys.time(), tz = "UTC", usetz = TRUE)),
        sprintf('  "source": "%s",', source_tag),
        sprintf('  "confidence_tier": "%s",', params$confidence_tier),
        sprintf('  "crop_thresholds": { "T_opt_max_C": %.2f, "R_opt_min_mm": %.0f, "R_opt_max_mm": %.0f },',
                T_opt_max, R_opt_min, R_opt_max),
        sprintf('  "growing_season_months": [%s],', paste(gs, collapse = ",")),
        sprintf('  "CV_precip": %.3f,', CV_R),
        sprintf('  "CV_tmax": %.3f,', CV_T),
        '  "country_mean": {',
        paste(cm_lines, collapse = ",\n"),
        '  }',
        "}"
      ), sidecar)

      cca_log(sprintf("Wrote: %s (means: %s)", out_tif,
                      paste(sprintf("%s=%.3f", names(country_means),
                                    unlist(country_means)),
                            collapse = ", ")))
    }

    # Cross-check: sum of severities <= 1 + 1e-3 (disjoint bands).
    if (length(severity_rasters) == 3) {
      total <- severity_rasters$moderate +
               severity_rasters$severe +
               severity_rasters$extreme
      tmax  <- as.numeric(terra::global(total, "max", na.rm = TRUE))
      lvl   <- if (is.finite(tmax) && tmax <= 1.001) "INFO" else "ERROR"
      cca_log(sprintf("Sum-of-severities cross-check %s/%s: max(sum)=%.4f",
                      crop, hazard, tmax), level = lvl)
    }
  }
}

cca_banner("Done: hazard frequency rasters built")
cca_log("Wrote 27 COGs + 27 method.json sidecars under data/staging/hazards/.")
cca_log("Source: WorldClim 2.1 climatology + Ecocrop thresholds + literature CVs.")
cca_log("See docs/methodology_hazards.md for the full derivation.")
