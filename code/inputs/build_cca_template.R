# ------------------------------------------------------------------------------
# code/inputs/build_cca_template.R
#
# Builds the foundational spatial template for the CCA pipeline:
#   * a 5-arcmin reference SpatRaster covering Nigeria + Kenya + Ethiopia;
#   * country-specific templates (one per ISO3);
#   * GADM adm0/adm1/adm2 vectors (GeoParquet + GeoPackage dual write);
#   * rasterized admin masks (one integer layer per admin level) at 5 arcmin.
#
# Why this script exists:
#   Every downstream script (Step 1 mitigation, Step 3 adaptation, NEB) will
#   reproject or mask to ONE canonical template. Centralising the template
#   here prevents subtle cell-alignment drift (terra reports "extents do not
#   match" if origins differ by <1e-8 from different float-literal sources).
#
# Why `geodata::gadm` instead of the atlas GAUL-2024 pipeline:
#   The project lead specified `geodata` as the admin source. GADM uses ISO3
#   country codes directly (matching our CCA_COUNTRIES constant), ships with
#   geodata on CRAN, and includes levels 0-3. We store GADM's adm0/1/2 so
#   downstream admin extraction (hazards_prototype convention) keeps working.
#   GAUL-2024 remains the atlas-wide default; CCA is specifically GADM so
#   country partners (who use GADM in their national MoA workflows) can
#   cross-reference our admin IDs directly.
#
# Why we snap the extent to a 5/60 grid:
#   GADM country polygons have arbitrary decimal extents. If we rasterize to
#   those extents directly, the resulting grid is offset from MapSPAM and
#   from the atlas hazard rasters, forcing a second reprojection on every
#   merge. Snapping to a grid anchored at (-180, -90) with step = 5/60 makes
#   our template interoperable with atlas outputs at the same resolution.
# ------------------------------------------------------------------------------

# ---- Preamble: source the common/ modules in order -----------------------
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                      error = function(e) "")
common_dir <- if (nzchar(this_file)) {
  normalizePath(file.path(dirname(this_file), "..", "common"), mustWork = FALSE)
} else {
  # Interactive fallback: assume CWD is ccadapt/ or ccadapt/code/inputs.
  candidates <- c("code/common", "../common", "../../code/common")
  hit <- candidates[file.exists(file.path(candidates, "config.R"))][1]
  if (is.na(hit)) stop("Cannot locate code/common from getwd() = ", getwd())
  hit
}
source(file.path(common_dir, "config.R"))
source(file.path(common_dir, "logging.R"))
source(file.path(common_dir, "geo_utils.R"))

suppressPackageStartupMessages({
  library(terra)
  library(geodata)       # gadm() admin boundaries
})

cca_banner("Build CCA spatial template (NGA/KEN/ETH at 5 arcmin)")

# Short-circuit: if the regional template + all per-country templates +
# admin rasters already exist, skip the whole build. CCA_FORCE=1 overrides.
expected_outputs <- c(
  cca_template_path(),
  vapply(CCA_COUNTRIES, cca_template_path, character(1)),
  as.vector(outer(
    CCA_COUNTRIES,
    c("adm0", "adm1", "adm2"),
    FUN = function(iso, lv) file.path(CCA_PATHS$admin_dir,
      sprintf("gadm_%s_%s_%s.tif", iso, lv, CCA_RES_LABEL))
  ))
)
if (cca_skip_if_exists(expected_outputs, tag = "cca_template")) {
  cca_log("Template + admin masks already built. Done.")
  quit(status = 0)
}

# ---- 1. Download GADM admin boundaries for each country -------------------
# geodata::gadm downloads the ZIP to `path = ` and caches; re-running is a
# no-op after the first call. We request levels 0, 1, 2 (the levels that
# Activity-3 NEB aggregation and the country profiles will use).
gadm_cache <- file.path(CCA_PATHS$inputs_root, "gadm")
dir.create(gadm_cache, recursive = TRUE, showWarnings = FALSE)

gadm_levels <- c(0L, 1L, 2L)

fetch_gadm <- function(iso3, level) {
  cca_log(sprintf("GADM pull iso3=%s level=%d", iso3, level))
  # geodata::gadm returns a SpatVector; version defaults to the latest packaged
  # in the installed geodata build. We pin resolution = 1 (detailed) because
  # the CCA explorer needs polygonal detail at the sub-national level for
  # country-profile maps.
  geodata::gadm(country = iso3, level = level, resolution = 1, path = gadm_cache)
}

gadm_vects <- list()
for (iso3 in CCA_COUNTRIES) {
  gadm_vects[[iso3]] <- lapply(stats::setNames(gadm_levels, paste0("adm", gadm_levels)),
                               function(lv) fetch_gadm(iso3, lv))
}

# ---- 2. Persist GADM vectors as GeoParquet + GeoPackage -------------------
# Dual format: GeoParquet for the atlas data stack (arrow/DuckDB consumers),
# GeoPackage as a human-inspectable fallback (QGIS, ogr2ogr).
for (iso3 in CCA_COUNTRIES) {
  for (lv in paste0("adm", gadm_levels)) {
    v <- gadm_vects[[iso3]][[lv]]
    cca_write_vector(v, cca_admin_path(iso3, lv, ext = "gpkg"))
  }
}

# ---- 3. Compute the regional 5-arcmin template extent ---------------------
# Union the adm0 polygons of the three countries, then snap the bounding box
# outwards to the nearest 5-arcmin cell edge anchored at (-180, -90). This is
# the same snap strategy used in hazards_prototype to keep outputs aligned
# with MapSPAM (see spatial-foundation/build_base-raster.r: template origin is
# the global grid, not a country-specific one).
# NOTE: do.call(rbind, list_of_SpatVector) does not dispatch to terra's
# rbind method, so we fold the list explicitly. Likewise terra::ext() on a
# list does not work, so we combine bboxes arithmetically.
bbs <- lapply(gadm_vects, function(x) terra::ext(x$adm0))
bb <- terra::ext(
  min(vapply(bbs, function(e) e$xmin, numeric(1))),
  max(vapply(bbs, function(e) e$xmax, numeric(1))),
  min(vapply(bbs, function(e) e$ymin, numeric(1))),
  max(vapply(bbs, function(e) e$ymax, numeric(1)))
)

snap_out <- function(v, step) {
  # Floor lower edge, ceiling upper edge -- both snapped to `step`. The anchor
  # is implicitly (0,0) (global lon/lat origin), which matches the atlas
  # convention that rasters align with the -180..180 / -90..90 grid.
  c(floor(v[1] / step) * step,
    ceiling(v[2] / step) * step,
    floor(v[3] / step) * step,
    ceiling(v[4] / step) * step)
}

snapped <- snap_out(c(bb$xmin, bb$xmax, bb$ymin, bb$ymax), CCA_RES)
template_ext <- terra::ext(snapped[1], snapped[2], snapped[3], snapped[4])
cca_log(sprintf("Regional template extent (snapped): xmin=%.6f xmax=%.6f ymin=%.6f ymax=%.6f",
                snapped[1], snapped[2], snapped[3], snapped[4]))

template_3c <- terra::rast(
  extent     = template_ext,
  resolution = CCA_RES,
  crs        = CCA_CRS
)
terra::values(template_3c) <- NA_real_   # empty skeleton; downstream fills

# ---- 4. Write the regional template and country-specific variants --------
# The regional template is the primary reference grid. Country-specific
# templates are just the regional grid clipped to each country's adm0 extent;
# they exist so memory-lean country-by-country runs don't need to load the
# full 3-country raster.
cca_write_cog(template_3c, cca_template_path())

for (iso3 in CCA_COUNTRIES) {
  bb_c  <- terra::ext(gadm_vects[[iso3]]$adm0)
  snap_c <- snap_out(c(bb_c$xmin, bb_c$xmax, bb_c$ymin, bb_c$ymax), CCA_RES)
  tc <- terra::rast(
    extent     = terra::ext(snap_c[1], snap_c[2], snap_c[3], snap_c[4]),
    resolution = CCA_RES,
    crs        = CCA_CRS
  )
  terra::values(tc) <- NA_real_
  cca_write_cog(tc, cca_template_path(iso3))
}

# ---- 5. Rasterize admin masks at 5 arcmin -------------------------------
# Per country, per admin level: burn the GID_{level} attribute into the
# template grid. These rasters drive fast zonal statistics in downstream
# admin parquet generation (see hazards_prototype/R/3.1_risk_x_exposure_x_ac.R
# for the pattern: integer raster keyed to a lookup table).
for (iso3 in CCA_COUNTRIES) {
  tc <- terra::rast(cca_template_path(iso3))
  for (lv in gadm_levels) {
    lv_tag <- paste0("adm", lv)
    v <- gadm_vects[[iso3]][[lv_tag]]
    # GADM uses column names like GID_0, GID_1, GID_2. These are character
    # (e.g. "NGA.1.3_1") so rasterize needs a numeric key. We build a stable
    # integer ID from row order and write the CSV lookup alongside.
    v$adm_id <- seq_len(nrow(v))
    r <- terra::rasterize(v, tc, field = "adm_id")
    out_tif <- file.path(CCA_PATHS$admin_dir,
                         sprintf("gadm_%s_%s_%s.tif", iso3, lv_tag, CCA_RES_LABEL))
    cca_write_cog(r, out_tif)

    # Lookup table: adm_id -> GID_{lv} and NAME_{lv}. Used downstream to join
    # admin rasters with the vector attributes without round-tripping through
    # the SpatVector for every zonal call.
    lookup_cols <- c("adm_id",
                     paste0("GID_", lv),
                     paste0("NAME_", lv))
    present <- lookup_cols %in% names(v)
    lookup <- as.data.frame(v)[, lookup_cols[present], drop = FALSE]
    out_csv <- file.path(CCA_PATHS$admin_dir,
                         sprintf("gadm_%s_%s_lookup.csv", iso3, lv_tag))
    utils::write.csv(lookup, out_csv, row.names = FALSE)
    cca_log("Wrote lookup: ", out_csv)
  }
}

# ---- 6. Final validation --------------------------------------------------
# Sanity-check the regional template: CRS, resolution, non-empty extent.
regional <- terra::rast(cca_template_path())
cca_assert_crs(regional)
cca_assert_res(regional, CCA_RES)
cca_assert(terra::ncell(regional) > 0, "Regional template has zero cells")
cca_log(sprintf("Regional template OK: %d cols x %d rows (%d cells)",
                terra::ncol(regional), terra::nrow(regional), terra::ncell(regional)))

cca_banner("Done: spatial template built")
