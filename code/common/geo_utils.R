# ------------------------------------------------------------------------------
# code/common/geo_utils.R
#
# Thin geospatial helpers shared by every Activity-2 / Activity-3 script. The
# rule is: if two scripts need the same transformation, it lives here. That
# keeps reprojection parameters, COG write options, and optional masking
# semantics consistent across the pipeline.
#
# Non-negotiable design calls (per project lead guidance):
#   * No default spatial masking. All constraint-mask functions take a
#     `distance_km=` / `threshold=` / `buffer_km=` argument that defaults to
#     NULL, which means the function is a no-op pass-through. Real values
#     get set by Activity-1 stakeholder workshops, not by the coder.
#   * Area-conserving reprojection for any "per-area" quantity (tCO2e,
#     hectares, USD/ha). We use terra::project(..., method = "average") when
#     going from fine -> coarse (Naturebase 1 km -> CCA 5 arcmin) because
#     bilinear would introduce edge artefacts on the cropland boundary.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
})

# geo_utils.R assumes config.R and logging.R have already been sourced by the
# entry-point script. This is deliberate: `sys.frame(1)$ofile` is unreliable
# in interactive / reticulate / testthat contexts. Entry-points under
# code/activity2/* and code/inputs/* source the three common/ files explicitly.
if (!exists("CCA_ROOT") || !exists("cca_log")) {
  stop("[CCA] geo_utils.R requires config.R and logging.R to be sourced first. ",
       "See code/inputs/build_cca_template.R for the expected preamble.")
}

# null-coalesce operator, used below
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- COG writer -----------------------------------------------------------
cca_write_cog <- function(r, path, overwrite = TRUE) {
  # Wrapper around terra::writeRaster that enforces the project-wide COG
  # conventions. Writing to a temp path first and renaming avoids leaving a
  # half-written file on disk if R is killed mid-write.
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  terra::writeRaster(
    r, tmp,
    filetype  = "COG",
    gdal      = CCA_COG_OPTS,
    overwrite = overwrite
  )
  file.rename(tmp, path)
  cca_log("Wrote COG: ", path, " (", file.size(path), " bytes)")
  invisible(path)
}

# Primary vector format for the CCA pipeline is GeoPackage (.gpkg). Every
# GDAL build ships the GPKG driver, it supports layers, indexes and large
# geometries, and QGIS / ogr2ogr / terra / sf all open it natively. For
# web consumption we also accept GeoJSON via `also_geojson = TRUE`
# (small files only -- bloats rapidly for admin2 at full detail).
cca_write_vector <- function(v, path, also_geojson = FALSE, overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # Auto-route by extension so callers can keep using their preferred stem.
  ext <- tolower(tools::file_ext(path))
  if (ext == "" || ext == "parquet") {
    # parquet kept as a *logical* stem because earlier code referenced it; we
    # always write GPKG now.
    path <- sub("\\.parquet$", ".gpkg", path)
    if (tolower(tools::file_ext(path)) == "")
      path <- paste0(path, ".gpkg")
  }
  terra::writeVector(v, path, filetype = "gpkg", overwrite = overwrite)
  cca_log("Wrote GPKG: ", path)

  if (also_geojson) {
    geo <- sub("\\.gpkg$", ".geojson", path)
    terra::writeVector(v, geo, filetype = "GeoJSON", overwrite = overwrite)
    cca_log("Wrote GeoJSON: ", geo)
  }
  invisible(path)
}

# Back-compat alias for earlier code that called cca_write_geoparquet; now
# just writes GPKG.
cca_write_geoparquet <- function(v, path, overwrite = TRUE) {
  cca_write_vector(v, path, overwrite = overwrite)
}

# Tolerant reader: accepts .gpkg, .geojson, or a bare stem.
cca_read_vector <- function(path) {
  if (file.exists(path)) return(terra::vect(path))
  for (ext in c(".gpkg", ".geojson", ".parquet")) {
    alt <- paste0(tools::file_path_sans_ext(path), ext)
    if (file.exists(alt)) return(terra::vect(alt))
  }
  stop("[CCA] No vector file found for stem ", path)
}

# ---- Skip-if-output-exists helper ---------------------------------------
# Every fetcher / builder guards its heavy work with
#     if (cca_skip_if_exists(output_paths, tag = "biochar_potential")) return()
# which returns TRUE (and logs) if ALL listed paths exist with non-zero size.
# Callers can force a rebuild by setting env var CCA_FORCE=1.
cca_skip_if_exists <- function(paths, tag = "output") {
  if (isTRUE(Sys.getenv("CCA_FORCE") %in% c("1", "true", "TRUE"))) return(FALSE)
  if (!length(paths)) return(FALSE)
  ok <- all(vapply(paths, function(p) file.exists(p) && file.size(p) > 0L, logical(1)))
  if (ok) {
    cca_log(sprintf("SKIP %s: all output files present (set CCA_FORCE=1 to rebuild)", tag))
  }
  ok
}

# ---- Reprojection ---------------------------------------------------------
cca_reproject_to_template <- function(
  src,
  template,
  method = c("average", "bilinear", "near", "sum", "mode"),
  name   = "raster"
) {
  # `method = "average"` is the default for area-conserving quantities --
  # it takes the mean of contributing fine-resolution cells, weighted
  # implicitly by the number of cells overlapping each coarse cell. This is
  # the terra equivalent of gdalwarp -r average, which is what the
  # spatial-foundation and hazards_prototype pipelines use for MapSPAM-to-hazard
  # alignment (see atlas/hazards_prototype/R/3_freq_x_exposure.R).
  method <- match.arg(method)
  cca_log("Reprojecting ", name, " -> template (method=", method, ")")
  out <- terra::project(src, template, method = method, align = TRUE)
  cca_assert_crs(out)
  cca_assert_res(out, terra::res(template))
  out
}

# ---- Optional spatial masks ----------------------------------------------
# Every function below is a NO-OP unless the caller passes a concrete
# parameter value. This is a deliberate project-lead decision: the CCA team
# will not silently encode masking defaults (protected areas, slope cut-offs,
# basalt proximity) that belong in the stakeholder-workshop-driven config.

cca_mask_by_proximity <- function(
  r,
  points,
  distance_km = NULL,
  name = "proximity"
) {
  # r        : SpatRaster to optionally mask
  # points   : SpatVector of feature locations (e.g. basalt quarries)
  # distance_km : if NULL, returns r unchanged. If numeric, keeps only pixels
  #               whose cell centre lies within distance_km of any point.
  if (is.null(distance_km)) {
    cca_log(sprintf("%s proximity mask: SKIPPED (distance_km = NULL)", name))
    return(r)
  }
  cca_assert(is.numeric(distance_km) && distance_km > 0,
             "distance_km must be a positive number or NULL")
  cca_log(sprintf("%s proximity mask: applying at %.1f km", name, distance_km))
  # terra::distance returns metres for lon/lat rasters
  d <- terra::distance(r, points)
  terra::mask(r, d <= distance_km * 1000, maskvalues = FALSE)
}

cca_mask_by_threshold <- function(
  r,
  layer,
  threshold = NULL,
  direction = c(">=", "<=", ">", "<"),
  name = "threshold"
) {
  # Generic pixel-level mask. `layer` is another SpatRaster (e.g. slope, OC);
  # `threshold` is the cut-off. When threshold is NULL the input is returned
  # unchanged.
  direction <- match.arg(direction)
  if (is.null(threshold)) {
    cca_log(sprintf("%s mask: SKIPPED (threshold = NULL)", name))
    return(r)
  }
  cca_assert(is.numeric(threshold),
             "threshold must be numeric or NULL")
  cca_log(sprintf("%s mask: keeping layer %s %g", name, direction, threshold))
  keep <- switch(direction,
    ">=" = layer >= threshold,
    "<=" = layer <= threshold,
    ">"  = layer >  threshold,
    "<"  = layer <  threshold
  )
  terra::mask(r, keep, maskvalues = FALSE)
}

cca_mask_by_polygon <- function(
  r,
  polygons = NULL,
  inverse = FALSE,
  name = "polygon"
) {
  # Apply a vector mask only if `polygons` is non-NULL. `inverse = TRUE`
  # excludes pixels inside the polygon (useful for WDPA protected-area
  # exclusion); `inverse = FALSE` keeps only pixels inside the polygon.
  if (is.null(polygons)) {
    cca_log(sprintf("%s polygon mask: SKIPPED (polygons = NULL)", name))
    return(r)
  }
  cca_log(sprintf("%s polygon mask: applying (inverse=%s)", name, inverse))
  terra::mask(r, polygons, inverse = inverse)
}
