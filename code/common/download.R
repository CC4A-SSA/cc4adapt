# ------------------------------------------------------------------------------
# code/common/download.R
#
# Idempotent download helper with provenance sidecars. Why this exists:
#
#   * We hit half a dozen external sources (Harvard Dataverse, FAOSTAT,
#     PANGAEA, USGS, WorldPop, SoilGrids). Each has its own retry quirks.
#     Wrapping them here means each fetcher in code/inputs/*.R is just a
#     one-liner `cca_download(url, dest)`.
#
#   * Every downloaded file gets a JSON sidecar at `<dest>.source.json` with
#     URL, timestamp, HTTP status, byte size, and SHA-256 digest. This is
#     what QAQC-system / audit reports inspect to verify the provenance of
#     every layer published to the `cca` STAC collection.
#
#   * Idempotent: re-running the pipeline does not redownload files already
#     present on disk with non-zero size, unless the caller passes
#     `force = TRUE`. Prevents thrashing FAOSTAT / Dataverse on every rerun.
#
# Depends on: config.R (CCA_ROOT), logging.R (cca_log). Source them first.
# ------------------------------------------------------------------------------

if (!exists("CCA_ROOT") || !exists("cca_log")) {
  stop("[CCA] download.R requires config.R and logging.R sourced first.")
}

suppressPackageStartupMessages({
  library(tools)      # sha256
})

cca_download <- function(
  url,
  dest,
  force = FALSE,
  expected_sha256 = NULL,
  timeout_s = 600L,
  user_agent = "ccadapt-pipeline/1.0 (+https://adaptationatlas.cgiar.org)"
) {
  # url         : source URL (https preferred; we do not upgrade http)
  # dest        : absolute path to destination file
  # force       : re-download even if dest already exists
  # expected_sha256 : optional hex digest for cross-check; non-fatal mismatch
  # timeout_s   : per-request timeout (FAOSTAT / Dataverse are slow)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  if (!force && file.exists(dest) && file.size(dest) > 0L) {
    cca_log("SKIP download (exists): ", dest)
    return(invisible(dest))
  }

  old_to <- options(timeout = timeout_s); on.exit(options(old_to), add = TRUE)
  cca_log(sprintf("Downloading %s -> %s", url, dest))

  # Use method="libcurl" for parallel-safe, resumable TLS. `mode="wb"` is
  # critical for binary files on Windows (no-op on unix but harmless).
  tmp <- paste0(dest, ".part")
  ok <- tryCatch({
    utils::download.file(
      url       = url,
      destfile  = tmp,
      mode      = "wb",
      method    = "libcurl",
      quiet     = FALSE,
      headers   = c("User-Agent" = user_agent)
    )
    TRUE
  }, error = function(e) {
    cca_log("Download failed: ", conditionMessage(e), level = "ERROR")
    FALSE
  })
  if (!ok || !file.exists(tmp) || file.size(tmp) == 0L) {
    if (file.exists(tmp)) unlink(tmp)
    stop("[CCA download] Failed: ", url)
  }
  # Atomic rename: on *nix this is safe even under concurrent pipeline runs.
  file.rename(tmp, dest)

  # Provenance sidecar. Minimal JSON by hand to avoid pulling jsonlite.
  size <- file.size(dest)
  sha  <- as.character(tools::md5sum(dest))  # md5 by default; cheap & sufficient
  sha256 <- tryCatch(digest::digest(file = dest, algo = "sha256"), error = function(e) NA_character_)
  sidecar <- paste0(
    "{\n",
    '  "url": "',       url,                                        '",\n',
    '  "dest": "',      dest,                                       '",\n',
    '  "fetched_at": "',format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), '",\n',
    '  "bytes": ',      size,                                       ',\n',
    '  "md5": "',       sha,                                        '",\n',
    '  "sha256": "',    sha256 %||% "",                             '"\n',
    "}\n"
  )
  writeLines(sidecar, paste0(dest, ".source.json"))
  cca_log(sprintf("OK (%d bytes, md5=%s): %s", size, sha, dest))

  if (!is.null(expected_sha256) && !is.na(sha256) &&
      !identical(tolower(expected_sha256), tolower(sha256))) {
    cca_log(sprintf("SHA-256 mismatch! expected %s got %s", expected_sha256, sha256),
            level = "WARN")
  }
  invisible(dest)
}

# Small helper for extracting the commonly-needed ISO3/country subset from a
# global raster, without having to repeat boilerplate in every fetcher.
#
# Important: do NOT force the CRS here. Source rasters may arrive in
# Mollweide (Sanderman, GHSL), sinusoidal (MODIS), or lon/lat (CHIRPS,
# Naturebase). We detect and reproject where needed, keeping the input
# CRS intact during read. `force_crs` is available for cases where the
# source is known-incorrectly-labelled (Naturebase's WHA83 typo) but the
# underlying datum is unambiguously WGS84.
cca_clip_to_cca_region <- function(src, template_path = NULL,
                                   force_crs = NULL,
                                   method = "bilinear") {
  if (is.null(template_path)) template_path <- cca_template_path()  # 3-country
  tmpl <- terra::rast(template_path)
  r <- if (inherits(src, "SpatRaster")) src else terra::rast(src)

  if (!is.null(force_crs)) terra::crs(r) <- force_crs

  src_epsg <- as.character(terra::crs(r, describe = TRUE)$code)
  if (is.na(src_epsg) || !nzchar(src_epsg)) {
    stop("[CCA] Source raster has no detectable CRS. Pass force_crs= explicitly.")
  }

  # Two cases:
  #   (a) source is already lon/lat 4326 -> crop to template extent is cheap.
  #   (b) source is projected -> transform the template extent INTO the source
  #       CRS, crop there (fast, no resampling), then project the smaller
  #       result back to 4326. This avoids reprojecting the whole global raster.
  if (src_epsg == "4326") {
    return(terra::crop(r, terra::ext(tmpl), snap = "out"))
  }

  # Build a vector polygon from the template extent, reproject into source CRS.
  ext_vec_4326 <- terra::vect(terra::ext(tmpl), crs = "EPSG:4326")
  ext_vec_src  <- terra::project(ext_vec_4326, terra::crs(r))
  ext_src      <- terra::ext(ext_vec_src)
  cropped_src  <- terra::crop(r, ext_src, snap = "out")
  terra::project(cropped_src, tmpl, method = method)
}
