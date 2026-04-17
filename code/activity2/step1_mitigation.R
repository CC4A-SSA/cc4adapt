# ------------------------------------------------------------------------------
# code/activity2/step1_mitigation.R
#
# Activity 2 Step 1: construct the spatial mitigation-potential layer
# (tCO2e/ha/yr) for each (country, practice) combination the project plan
# prioritises.
#
# Inputs consumed:
#   * data/naturebase/all_pathways/*_total.tif          (Naturebase 1-km rasters)
#   * data/staging/template/cca_template_<iso3>_5arcmin.tif  (built by
#     code/inputs/build_cca_template.R)
#   * data/staging/admin/gadm_<iso3>_adm0_<RES>.tif     (admin mask)
#   * code/activity2/config/practices.yaml              (practice catalogue)
#
# Outputs produced:
#   * data/outputs/step1_mitigation/mitigation_<practice>_<iso3>_5arcmin.tif
#     (tCO2e/ha/yr density, COG)
#   * data/outputs/step1_mitigation/mitigation_admin_<iso3>.parquet
#     (per-admin-unit totals; columns adm_id, practice, tco2e_yr, area_ha,
#      tco2e_ha_yr)
#
# --- Key design decisions (cross-checked against Naturebase ReadMe +
#     Naturebase_metadata.csv, investment document Activity 2 Step 1, and the
#     spatial-foundation / hazards_prototype atlas conventions) ---
#
# [A] Resampling strategy: total -> coarse sum -> density
#     Naturebase rasters store totals per 1-km pixel (filenames end in
#     "_total.tif" and the project ReadMe gives `tco2e/ha = tco2e / ha` for
#     aggregation). To preserve totals when moving from 1 km to 5 arcmin
#     (~100 km^2 at the equator), we reproject with method = "sum". We then
#     divide by the area of each coarse cell (terra::cellSize in hectares) to
#     get a stable density in tCO2e/ha/yr -- the unit Activity-3 NEB expects.
#     An alternative (resample with "average" on a pre-computed density
#     surface) would double-count along partial cells at country borders.
#
# [B] No implicit masking
#     Project lead: "Do not make any decisions such as buffer or spatial
#     masking by distance from road. Wrap such operations in a function with
#     masking distance option but with no masking applied."
#     => `apply_crop_mask` and `apply_country_clip` are parameters with
#        defaults such that NO masking runs unless explicitly requested.
#        Naturebase already bakes in cropland agreement (Tubiello 2022), so
#        re-masking to MapSPAM would drop legitimate pixels; cropland mask
#        defaults to FALSE.
#
# [C] Cross-check validation
#     For every Naturebase-sourced practice we compare country totals from
#     the reprojected raster against the Naturebase adm0 CSV values. Expected
#     tolerance: +/- 5 percent (the plan's cross-check criterion). A larger
#     deviation is logged as WARN for triage but does not abort.
#
# [D] External-pending practices
#     Biochar, AWD, ERW are not in the local Naturebase package. The script
#     skips them with a visible WARN so downstream steps know why there is no
#     raster. When the download scripts under code/inputs/ land those layers,
#     this script reads them through the same code path.
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
  library(arrow)
})

cca_banner("Activity 2 | Step 1 | Mitigation potential (Naturebase reprojection)")

# ---- Load practice catalogue -------------------------------------------
practices_yaml <- yaml::read_yaml(
  file.path(CCA_ROOT, "code", "activity2", "config", "practices.yaml")
)
practice_defs       <- practices_yaml$practices
country_practices   <- practices_yaml$country_practices

# Naturebase adm0 totals -- used for cross-check validation ([C] above).
# Columns in the file follow the pattern <pathway_code>_tco2e and
# <pathway_code>_ha (see Naturebase ReadMe_all_pathways.txt).
nb_adm0_csv <- file.path(CCA_PATHS$naturebase, "all_pathways_tco2e_adm0.csv")
nb_adm0 <- if (file.exists(nb_adm0_csv)) {
  utils::read.csv(nb_adm0_csv, stringsAsFactors = FALSE)
} else {
  cca_log(sprintf("Naturebase adm0 CSV missing (%s); cross-check skipped", nb_adm0_csv),
          level = "WARN")
  NULL
}

# ---- Core Step-1 worker --------------------------------------------------
#' Reproject a Naturebase total raster onto the CCA template and convert to
#' density (tCO2e/ha/yr).
#'
#' @param src_tif path to the Naturebase 1-km _total raster
#' @param template SpatRaster (country or 3-country template)
#' @param mask_vect optional SpatVector to clip to (NULL = no clip)
#' @param mask_crop optional cropland-proportion raster to multiply through
#'   (NULL = no additional mask). Passed in only when the caller deliberately
#'   wants to intersect with MapSPAM cropland; defaults to NULL per [B].
reproject_and_densify <- function(src_tif, template,
                                  src_is_total = TRUE,
                                  mask_vect = NULL,
                                  mask_crop = NULL) {
  cca_log("Reading: ", src_tif)
  src <- terra::rast(src_tif)
  # Naturebase metadata labels the CRS as "WGS84 +no_defs" but some files
  # have a minor typo ("WHA83") in their PROJ string. We force to EPSG:4326
  # because the underlying datum is unambiguously WGS84.
  terra::crs(src) <- CCA_CRS

  # [A] Normalise to density before reprojection, then average onto template.
  # Rationale: if the source layer is a TOTAL (tCO2e per pixel at fine res),
  # we first divide by the source pixel area in ha to get tCO2e/ha/yr at the
  # source resolution. We then `terra::project(..., y=template, method="average")`
  # which guarantees output extent == template extent AND computes a properly
  # area-weighted mean of contributing source cells. This is mathematically
  # equivalent to "sum the totals, divide by coarse cell area" but without the
  # extent-alignment quirks of terra::project(..., method="sum", align=TRUE).
  # Staged layers (Sanderman MG bundle, AWD T1) already carry density units;
  # callers pass src_is_total = FALSE so we skip the divide-by-area step.
  if (src_is_total) {
    src_area_ha <- terra::cellSize(src, unit = "ha")
    src <- src / src_area_ha
  }
  density <- terra::project(src, template, method = "average")
  # Explicit crop+extend as belt-and-braces: force density to the template grid
  # even if terra::project surprisingly keeps a pre-existing extent.
  if (!terra::compareGeom(density, template, stopOnError = FALSE)) {
    density <- terra::resample(density, template, method = "average")
  }
  area_ha <- terra::cellSize(density, unit = "ha")

  # [B] Optional clip + mask, NO-OP by default.
  density <- cca_mask_by_polygon(density, polygons = mask_vect, name = "country clip")
  if (!is.null(mask_crop)) {
    density <- density * mask_crop
    cca_log("Applied MapSPAM cropland proportion mask (explicit opt-in).")
  }

  # Reconstruct per-cell total (density * area_ha) for the admin aggregation
  # block downstream. With the normalise-then-average path above, totals are
  # density * coarse cell area; this matches what the old sum-then-density
  # approach produced.
  totals_per_pixel <- density * area_ha
  list(density = density, totals_per_pixel = totals_per_pixel, area_ha = area_ha)
}

# ---- Driver: iterate over country x practice ----------------------------
admin_rows <- list()

for (iso3 in CCA_COUNTRIES) {
  cca_banner(sprintf("Country: %s", iso3))

  template_path <- cca_template_path(iso3)
  if (!file.exists(template_path)) {
    cca_log(sprintf("Template missing for %s -- run build_cca_template.R first. Skipping.",
                    iso3), level = "ERROR")
    next
  }
  template <- terra::rast(template_path)

  adm0 <- cca_read_vector(cca_admin_path(iso3, "adm0"))

  wanted <- country_practices[[iso3]] %||% character(0)
  cca_log(sprintf("%s practices: %s", iso3, paste(wanted, collapse = ", ")))

  for (practice in wanted) {
    def <- practice_defs[[practice]]
    if (is.null(def)) {
      cca_log(sprintf("Unknown practice '%s' -- skipping", practice), level = "WARN")
      next
    }

    # [D] Route by source.
    src <- def$mitigation$source
    src_tif <- NULL
    if (src == "naturebase") {
      src_tif <- file.path(CCA_PATHS$naturebase, def$mitigation$file)
    } else if (src == "staged") {
      # Prefer an "alt_file" (e.g. Cao 2026) if present; else fallback file.
      candidates <- c(
        if (!is.null(def$mitigation$alt_file)) file.path(CCA_PATHS$staging_root, def$mitigation$alt_file),
        file.path(CCA_PATHS$staging_root, def$mitigation$file)
      )
      src_tif <- candidates[file.exists(candidates)][1]
    } else if (src == "staged_optional") {
      staged_tif <- file.path(CCA_PATHS$staging_root, def$mitigation$file)
      if (!file.exists(staged_tif)) {
        cca_log(sprintf("Practice '%s' optional and not yet staged (%s) -- skipping %s",
                        practice, staged_tif, iso3), level = "WARN")
        next
      }
      src_tif <- staged_tif
    } else if (src == "external_pending") {
      cca_log(sprintf("Practice '%s' is external_pending (%s) -- skipping Step 1 for %s.",
                      practice, def$mitigation$citation %||% "", iso3),
              level = "WARN")
      next
    } else {
      cca_log(sprintf("Unsupported source '%s' for '%s' -- skipping.", src, practice),
              level = "WARN")
      next
    }

    if (is.null(src_tif) || is.na(src_tif) || !file.exists(src_tif)) {
      cca_log(sprintf("Source raster missing for practice '%s' (%s) -- skipping %s",
                      practice, src, iso3), level = "WARN")
      next
    }
    cca_log(sprintf("Using mitigation source for %s: %s", practice, src_tif))

    # Optional suitability multiplier (biochar / ERW) from practices.yaml.
    suitability_r <- NULL
    if (!is.null(def$mitigation$suitability_file)) {
      suit_path <- file.path(CCA_PATHS$staging_root, def$mitigation$suitability_file)
      if (file.exists(suit_path)) {
        suitability_r <- terra::rast(suit_path)
        cca_log("Applying suitability scaling from ", suit_path)
      } else {
        cca_log(sprintf("Suitability file listed but missing: %s", suit_path),
                level = "WARN")
      }
    }

    # Naturebase layers are totals per source pixel (file names end "_total.tif",
    # per Naturebase ReadMe). Staged layers (Sanderman MG bundle, AWD T1
    # fallback, ERW potential) are already density (tCO2e/ha/yr per their
    # filename / method.json). We switch the normalisation step based on the
    # practice source tag.
    is_total <- src == "naturebase"
    result <- reproject_and_densify(
      src_tif      = src_tif,
      template     = template,
      src_is_total = is_total,
      mask_vect    = adm0,
      mask_crop    = NULL
    )

    # Suitability scaling: density * suitability (0..1). Only done if a
    # suitability raster was supplied in the practice YAML; never implicit.
    if (!is.null(suitability_r)) {
      suit_aligned <- terra::project(suitability_r, result$density, method = "bilinear")
      result$density <- result$density * suit_aligned
    }

    # ---- Validate output --------------------------------------------------
    cca_assert_crs(result$density)
    cca_assert_res(result$density, CCA_RES)
    cca_assert_not_allna(result$density, name = sprintf("mitigation %s/%s", practice, iso3))
    cca_assert_range(result$density, lo = -5, hi = 100,
                     name = sprintf("mitigation %s/%s (tCO2e/ha/yr)", practice, iso3))

    # ---- Write COG --------------------------------------------------------
    out_tif <- file.path(
      CCA_PATHS$step1_dir,
      sprintf("mitigation_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL)
    )
    cca_write_cog(result$density, out_tif)

    # ---- Admin aggregation ------------------------------------------------
    # total tCO2e/yr for the country = sum of (density * area_ha).
    country_total <- as.numeric(terra::global(result$totals_per_pixel, "sum", na.rm = TRUE))
    country_area  <- as.numeric(terra::global(result$area_ha, "sum", na.rm = TRUE))

    admin_rows[[length(admin_rows) + 1L]] <- data.frame(
      iso3        = iso3,
      practice    = practice,
      tco2e_yr    = country_total,
      area_ha     = country_area,
      tco2e_ha_yr = country_total / country_area,
      stringsAsFactors = FALSE
    )

    # ---- [C] Cross-check vs Naturebase adm0 CSV --------------------------
    # Only applicable for naturebase-sourced practices; for staged layers
    # (Sanderman MG bundle, Cao AWD, ERW) we rely on the per-layer method.json
    # sidecars for provenance instead.
    if (!is.null(nb_adm0) && src == "naturebase") {
      col_tco2e <- paste0(sub("_total$", "", tools::file_path_sans_ext(def$mitigation$file)),
                          "_tco2e")
      iso3_col <- grep("ISO|iso3|GID_0|gaul0_code", names(nb_adm0),
                       value = TRUE, ignore.case = TRUE)[1]
      if (!is.na(iso3_col) && col_tco2e %in% names(nb_adm0)) {
        ref_row <- nb_adm0[nb_adm0[[iso3_col]] == iso3, , drop = FALSE]
        if (nrow(ref_row) == 1L) {
          ref_val <- suppressWarnings(as.numeric(ref_row[[col_tco2e]]))
          if (!is.na(ref_val) && ref_val > 0) {
            dev <- abs(country_total - ref_val) / ref_val
            lvl <- if (dev <= 0.05) "INFO" else if (dev <= 0.15) "WARN" else "ERROR"
            cca_log(sprintf("%s / %s cross-check: raster=%.2e reference=%.2e deviation=%.1f%%",
                            iso3, practice, country_total, ref_val, dev * 100),
                    level = lvl)
          }
        }
      }
    }
  }
}

# ---- Persist admin-level summary (CSV; small, human-readable) ----------
if (length(admin_rows)) {
  admin_df <- do.call(rbind, admin_rows)
  for (iso3 in unique(admin_df$iso3)) {
    out_csv <- file.path(
      CCA_PATHS$step1_dir,
      sprintf("mitigation_admin_%s.csv", iso3)
    )
    utils::write.csv(admin_df[admin_df$iso3 == iso3, , drop = FALSE], out_csv, row.names = FALSE)
    cca_log("Wrote admin summary: ", out_csv)
  }
}

cca_banner("Done: Step 1 mitigation")
