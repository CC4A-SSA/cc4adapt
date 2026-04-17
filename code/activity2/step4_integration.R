# ------------------------------------------------------------------------------
# code/activity2/step4_integration.R
#
# Activity 2 Step 4: pixel-level classification + contiguous opportunity
# zones.
#
# Inputs  (per practice x country):
#   M = Step 1 mitigation tCO2e/ha/yr
#   P = Step 2 dvalue USD/ha/yr  (production benefit, summed across crops)
#   A = Step 3 avoided-loss USD/ha/yr (adaptation benefit)
#
# Method:
#   1. Normalise each layer to [0, 1] within country by rank, to remove
#      unit-scale effects (tCO2e vs USD are not comparable on a fixed cutoff).
#   2. Compare each layer to the country-specific quantile cutoff
#      (default 0.6; set in config/thresholds.yaml). Above cutoff = "high".
#   3. Classify each pixel:
#        1 = Triple Win       (M high,  P high, A high)
#        2 = Mitigation+Prod  (M high,  P high, A low)
#        3 = Mitigation+Adapt (M high,  P low,  A high)
#        4 = Constrained      (M high,  P low,  A low)
#        5 = Non-priority     (M low)
#   4. Extract contiguous patches (terra::patches) at configurable
#      connectivity; polygonise; compute zone stats.
#
# Outputs:
#   * data/outputs/step4_integration/classification_<practice>_<iso3>.tif (cat)
#   * data/outputs/step4_integration/opportunity_zones_<practice>_<iso3>.gpkg
#     with attributes area_ha, mean_M, mean_P, mean_A, class.
# ------------------------------------------------------------------------------

# ---- Preamble ---------------------------------------------------------
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

cca_banner("Activity 2 | Step 4 | Integration + opportunity zones")

practices_yaml    <- yaml::read_yaml(file.path(CCA_ROOT, "code", "activity2", "config", "practices.yaml"))
country_practices <- practices_yaml$country_practices
thresholds <- yaml::read_yaml(file.path(CCA_ROOT, "code", "activity2", "config", "thresholds.yaml"))
QCUT     <- thresholds$classification$quantile_cutoff  %||% 0.6
MIN_HA   <- thresholds$opportunity_zones$min_area_ha   %||% 0
CONN     <- thresholds$opportunity_zones$connectivity  %||% 8

# Score a raster to [0, 1] by quantile rank inside the country.
rank_normalise <- function(r) {
  v <- terra::values(r)
  keep <- is.finite(v)
  if (!any(keep)) return(r)
  ranks <- rep(NA_real_, length(v))
  ranks[keep] <- (rank(v[keep], ties.method = "average") - 1) / (sum(keep) - 1)
  terra::setValues(r, ranks)
}

for (iso3 in CCA_COUNTRIES) {
  cca_banner(sprintf("Country: %s", iso3))
  template <- terra::rast(cca_template_path(iso3))
  wanted   <- country_practices[[iso3]] %||% character(0)

  for (practice in wanted) {
    cls_tif <- file.path(CCA_PATHS$step4_dir,
                         sprintf("classification_%s_%s.tif", practice, iso3))
    zones_gpkg <- file.path(CCA_PATHS$step4_dir,
                            sprintf("opportunity_zones_%s_%s.gpkg", practice, iso3))
    if (cca_skip_if_exists(c(cls_tif, zones_gpkg),
                           tag = sprintf("step4 %s/%s", practice, iso3))) next

    M_tif <- file.path(CCA_PATHS$step1_dir,
                       sprintf("mitigation_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    P_tif <- file.path(CCA_PATHS$step2_dir,
                       sprintf("dvalue_%s_%s_%s.tif",     practice, iso3, CCA_RES_LABEL))
    A_tif <- file.path(CCA_PATHS$step3_dir,
                       sprintf("avoidedloss_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))

    if (!file.exists(M_tif)) {
      cca_log(sprintf("No Step 1 output for %s/%s -- skipping Step 4", practice, iso3),
              level = "WARN"); next
    }
    M <- terra::rast(M_tif)
    P <- if (file.exists(P_tif)) terra::rast(P_tif) else {
           cca_log(sprintf("No Step 2 P for %s/%s; treating as zero", practice, iso3),
                   level = "WARN"); M * 0 }
    A <- if (file.exists(A_tif)) terra::rast(A_tif) else {
           cca_log(sprintf("No Step 3 A for %s/%s; treating as zero", practice, iso3),
                   level = "WARN"); M * 0 }

    # Align P, A to M (all should already be same grid but guard anyway)
    P <- terra::project(P, M, method = "bilinear")
    A <- terra::project(A, M, method = "bilinear")

    Mn <- rank_normalise(M)
    Pn <- rank_normalise(P)
    An <- rank_normalise(A)

    mhigh <- Mn >= QCUT
    phigh <- Pn >= QCUT
    ahigh <- An >= QCUT

    cls <- terra::rast(M); terra::values(cls) <- 5L   # default: Non-priority
    cls[mhigh & phigh & ahigh] <- 1L
    cls[mhigh & phigh & !ahigh] <- 2L
    cls[mhigh & !phigh & ahigh] <- 3L
    cls[mhigh & !phigh & !ahigh] <- 4L
    # class 5 = Non-priority (already set)
    cca_write_cog(cls, cls_tif)

    # Contiguous opportunity zones: take any "high-mitigation" pixel (classes 1-4)
    # as a seed; patches() groups contiguous seeds.
    seed <- cls <= 4
    patches <- terra::patches(seed, directions = CONN, allowGaps = FALSE)
    if (is.null(patches) || all(is.na(terra::values(patches)))) {
      cca_log(sprintf("No contiguous high-potential patches for %s/%s", practice, iso3),
              level = "INFO")
      next
    }

    # Polygonise and compute attributes
    zones <- terra::as.polygons(patches, dissolve = TRUE, values = TRUE)
    names(zones)[1] <- "patch_id"

    # Attributes: area_ha, zonal means of M/P/A, and dominant class
    area_ha <- terra::expanse(zones, unit = "ha")
    zone_stats <- data.frame(
      patch_id = zones$patch_id,
      area_ha  = area_ha,
      mean_M   = as.numeric(terra::extract(M, zones, fun = mean, na.rm = TRUE)[,2]),
      mean_P   = as.numeric(terra::extract(P, zones, fun = mean, na.rm = TRUE)[,2]),
      mean_A   = as.numeric(terra::extract(A, zones, fun = mean, na.rm = TRUE)[,2]),
      # terra::extract rejects na.rm when a custom fun is supplied, so we
      # drop NAs inside the function and omit na.rm from the outer call.
      dominant_class = as.integer(terra::extract(cls, zones,
                                                  fun = function(v) {
                                                    v <- v[!is.na(v)]
                                                    if (!length(v)) return(NA_integer_)
                                                    as.integer(names(sort(table(v), decreasing = TRUE))[1])
                                                  })[,2]),
      stringsAsFactors = FALSE
    )

    # Attach zone_stats as SpatVector attributes (by patch_id order)
    zones$area_ha   <- zone_stats$area_ha
    zones$mean_M    <- zone_stats$mean_M
    zones$mean_P    <- zone_stats$mean_P
    zones$mean_A    <- zone_stats$mean_A
    zones$class_dom <- zone_stats$dominant_class

    # Filter by minimum area if set
    if (MIN_HA > 0) {
      zones <- zones[zones$area_ha >= MIN_HA, ]
    }
    cca_write_vector(zones, zones_gpkg)
    cca_log(sprintf("Step 4 %s/%s: %d opportunity zones, total %.0f ha",
                    practice, iso3, nrow(zones), sum(zones$area_ha, na.rm = TRUE)))
  }
}

cca_banner("Done: Step 4 integration")
