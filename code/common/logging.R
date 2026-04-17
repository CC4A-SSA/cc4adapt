# ------------------------------------------------------------------------------
# code/common/logging.R
#
# Structured logging + assertion helpers. Two goals:
#   1. Every step prints a clear banner, input/output counts, and timing --
#      so a log tail tells the reader exactly what happened and what stage
#      failed without having to read the R trace.
#   2. Validation checks (CRS match, no all-NA raster, P10 <= P50 <= P90)
#      are expressed as helper functions that either pass silently or stop
#      with an informative message -- preventing bad artefacts from being
#      published to S3 under the `cca` STAC collection.
#
# We deliberately do NOT depend on `logger` / `futile.logger` -- adding a
# package for a handful of messages is overkill, and keeping base-R means the
# scripts run in any R session the team has (including the atlas server
# containers which pin tidyverse-only).
# ------------------------------------------------------------------------------

cca_log <- function(..., level = "INFO") {
  msg <- paste0(..., collapse = "")
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %-5s %s\n", ts, level, msg))
}

cca_banner <- function(title) {
  bar <- strrep("-", max(nchar(title) + 4, 60))
  cat("\n", bar, "\n", sprintf("-- %s", title), "\n", bar, "\n", sep = "")
}

# ---- Assertions -----------------------------------------------------------

cca_assert <- function(condition, msg) {
  if (!isTRUE(condition)) stop("[CCA assertion failed] ", msg, call. = FALSE)
  invisible(TRUE)
}

cca_assert_crs <- function(r, expected_crs = "EPSG:4326") {
  # terra stores CRS as a WKT string; compare using the canonical EPSG code.
  got <- terra::crs(r, describe = TRUE)$code
  exp <- sub("^EPSG:", "", expected_crs)
  cca_assert(identical(as.character(got), as.character(exp)),
             sprintf("CRS mismatch: expected EPSG:%s, got EPSG:%s", exp, got))
}

cca_assert_res <- function(r, expected_res, tol = 1e-9) {
  got <- terra::res(r)
  cca_assert(all(abs(got - expected_res) < tol),
             sprintf("Resolution mismatch: expected %.10f, got %s",
                     expected_res, paste(got, collapse = ", ")))
}

cca_assert_not_allna <- function(r, name = "raster") {
  # minmax returns NA if the raster is entirely NA.
  mm <- terra::minmax(r, compute = TRUE)
  cca_assert(!all(is.na(mm)),
             sprintf("%s is entirely NA -- upstream masking / input is likely broken", name))
}

cca_assert_range <- function(r, lo, hi, name = "raster", strict = FALSE) {
  # Range check. `strict = FALSE` (default) only warns on out-of-range values
  # because many upstream layers legitimately have a handful of sentinel
  # values we may want to inspect rather than abort on.
  mm <- terra::minmax(r, compute = TRUE)
  if (is.na(mm[1]) || is.na(mm[2])) return(invisible(FALSE))
  bad <- mm[1] < lo - 1e-9 || mm[2] > hi + 1e-9
  if (bad) {
    msg <- sprintf("%s range [%g, %g] outside expected [%g, %g]",
                   name, mm[1], mm[2], lo, hi)
    if (strict) stop("[CCA assertion failed] ", msg, call. = FALSE)
    cca_log(msg, level = "WARN")
  }
  invisible(!bad)
}

cca_assert_percentile_order <- function(p10, p50, p90, name = "uncertainty") {
  # Used by Step 5 outputs: every pixel must satisfy P10 <= P50 <= P90.
  # We check on a sample for speed; full-raster check would be >100 M comparisons.
  idx <- sample(terra::ncell(p50), size = min(1e5, terra::ncell(p50)))
  v10 <- terra::values(p10)[idx]
  v50 <- terra::values(p50)[idx]
  v90 <- terra::values(p90)[idx]
  keep <- stats::complete.cases(v10, v50, v90)
  cca_assert(all(v10[keep] <= v50[keep] + 1e-9) && all(v50[keep] <= v90[keep] + 1e-9),
             sprintf("%s percentile ordering violated on sampled pixels", name))
}
