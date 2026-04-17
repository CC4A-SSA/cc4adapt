# ------------------------------------------------------------------------------
# code/activity3/neb.R
#
# Activity 3 Net Economic Benefit. For every (country, practice, carbon_price,
# adoption) combination produce:
#
#   NEB_usd_ha_yr = cp * M + A + P - (transition_annuity + implementation + mrv*M)
#   NPV_usd_ha    = sum_{t=1..T} NEB_t / (1+r)^t    (transition paid at t=0)
#   NPV_national  = NPV_usd_ha * area_ha * adoption_frac
#
# where M, P, A are the Step-1 / Step-2 / Step-3 rasters, `cp` is the
# scenario carbon price in USD/tCO2e, and costs come from neb_config.yaml.
#
# Implementation notes:
#
#   [A] Transition cost is annualised over the horizon: transition_annuity =
#       transition * r / (1 - (1+r)^-T). This keeps NEB comparable year-to-year
#       (otherwise year 1 is hugely negative and later years positive).
#   [B] MRV cost is charged per tCO2e per year, so it scales with M.
#   [C] Adoption scales the national aggregation, NOT per-hectare NEB. The
#       per-hectare raster is unchanged; only the national total uses
#       adoption fraction.
#   [D] Negative NPV pixels ARE kept (useful evidence that a practice is
#       uneconomic in a given location at a given price).
#
# Outputs:
#   * data/outputs/activity3_neb/npv_<practice>_<iso3>_cp=<cp>_adopt=<a>.tif
#   * data/outputs/activity3_neb/neb_admin_<iso3>.csv
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

cca_banner("Activity 3 | Net Economic Benefit")

cfg <- yaml::read_yaml(file.path(CCA_ROOT, "code", "activity3", "neb_config.yaml"))
practices_yaml    <- yaml::read_yaml(file.path(CCA_ROOT, "code", "activity2", "config", "practices.yaml"))
country_practices <- practices_yaml$country_practices

T_HORIZON <- cfg$horizon_years
R_DISC    <- cfg$discount_rate

annuity_factor <- function(T_years, r) r / (1 - (1 + r)^(-T_years))
# e.g. T=20, r=0.08 -> 0.1019; transition*0.1019 becomes annualised cost.
trans_annuity <- annuity_factor(T_HORIZON, R_DISC)
npv_mult      <- (1 - (1 + R_DISC)^(-T_HORIZON)) / R_DISC   # sum_{t=1..T} (1+r)^-t

admin_rows <- list()

for (iso3 in CCA_COUNTRIES) {
  cca_banner(sprintf("Country: %s", iso3))
  wanted <- country_practices[[iso3]] %||% character(0)

  for (practice in wanted) {
    costs <- cfg$costs[[practice]]
    if (is.null(costs)) {
      cca_log(sprintf("No cost entry for '%s' in neb_config.yaml -- skipping", practice),
              level = "WARN"); next
    }

    M_tif <- file.path(CCA_PATHS$step1_dir,
                       sprintf("mitigation_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
    if (!file.exists(M_tif)) {
      cca_log(sprintf("Step 1 missing for %s/%s -- skipping NEB", practice, iso3),
              level = "WARN"); next
    }
    M <- terra::rast(M_tif)
    P <- {
      pt <- file.path(CCA_PATHS$step2_dir,
                      sprintf("dvalue_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
      if (file.exists(pt)) terra::rast(pt) else M * 0
    }
    A <- {
      at <- file.path(CCA_PATHS$step3_dir,
                      sprintf("avoidedloss_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
      if (file.exists(at)) terra::rast(at) else M * 0
    }

    # Align P, A to M if needed
    P <- terra::project(P, M, method = "bilinear")
    A <- terra::project(A, M, method = "bilinear")

    cell_area_ha <- terra::cellSize(M, unit = "ha")

    for (cp in cfg$carbon_price_usd_per_tco2e) {
      # Per-hectare annual NEB (independent of adoption)
      annual_cost <- costs$transition * trans_annuity + costs$implementation + costs$mrv * M
      annual_benefit <- cp * M + A + P
      neb_ha_yr <- annual_benefit - annual_cost
      npv_ha    <- neb_ha_yr * npv_mult

      for (ad in cfg$adoption_scenarios) {
        out_tif <- file.path(CCA_PATHS$activity3_dir,
                             sprintf("npv_%s_%s_cp=%s_adopt=%s.tif",
                                     practice, iso3,
                                     format(cp, nsmall = 0),
                                     format(ad,  nsmall = 2)))
        if (cca_skip_if_exists(out_tif,
                               tag = sprintf("neb %s/%s cp=%s ad=%s",
                                             practice, iso3, cp, ad))) next
        cca_write_cog(npv_ha, out_tif)   # per-hectare NPV is adoption-agnostic

        # National NPV summed over pixels, scaled by adoption.
        national_npv <- as.numeric(terra::global(npv_ha * cell_area_ha, "sum", na.rm = TRUE)) * ad
        admin_rows[[length(admin_rows) + 1L]] <- data.frame(
          iso3 = iso3, practice = practice,
          carbon_price = cp, adoption = ad,
          npv_per_ha_mean = as.numeric(terra::global(npv_ha, "mean", na.rm = TRUE)),
          national_npv_usd = national_npv,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(admin_rows)) {
  df <- do.call(rbind, admin_rows)
  for (iso3 in unique(df$iso3)) {
    out_csv <- file.path(CCA_PATHS$activity3_dir,
                         sprintf("neb_admin_%s.csv", iso3))
    utils::write.csv(df[df$iso3 == iso3, , drop = FALSE], out_csv, row.names = FALSE)
    cca_log("Wrote NEB admin summary: ", out_csv)
  }
}

cca_banner("Done: Activity 3 NEB")
