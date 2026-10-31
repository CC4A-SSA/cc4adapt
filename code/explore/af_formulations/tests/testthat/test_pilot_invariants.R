# Invariants on the ERA season-pilot staging files (steps 03b, 03c, 07, 08).
# Skipped when the files are absent (fresh clone before the pilot has run).
# Run: Rscript -e 'testthat::test_dir("code/explore/af_formulations/tests/testthat")'

root <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", "..", ".."))
IN <- file.path(root, "data", "staging", "af_formulations")
need <- function(f) { p <- file.path(IN, f); if (!file.exists(p)) testthat::skip(paste("no", f)); p }
rd <- function(f) utils::read.csv(need(f), stringsAsFactors = FALSE)

test_that("03b: yields in one unit, keys complete, irrigation flagged", {
  cmp <- rd("era_comparisons.csv")
  expect_true(all(cmp$Units == "Mg/ha"))
  expect_true(all(is.finite(cmp$MeanC) & is.finite(cmp$MeanT) & cmp$MeanC > 0 & cmp$MeanT > 0))
  expect_true(all(cmp$practice %in% c("cover_crops", "reduced_tillage")))
  expect_true(is.logical(cmp$irrigated))
  expect_true(max(cmp$MeanC) < 30)   # Mg/ha; a kg/ha row would show as hundreds or thousands
})

test_that("03c: one window per site-slot, windows parse, climatology windows capped", {
  w <- rd("era_site_windows.csv")
  expect_equal(anyDuplicated(w[c("Site.Key", "slot")]), 0)
  expect_true(all(grepl("^\\d{2}-\\d{2}:\\d{2}-\\d{2}$", w$fixed_season)))
  days <- function(tok) { a <- as.Date(paste0("2001-", substr(tok, 1, 5))); b <- as.Date(paste0("2001-", substr(tok, 7, 11)))
    as.integer(ifelse(b >= a, b - a, b + 365 - a)) + 1 }
  clim <- w[w$source == "climatology", ]
  expect_true(all(days(clim$fixed_season) <= 31 * 5))
  expect_true(all(w$window_status %in% c("ok", "slot_missing")))
})

test_that("07: bands follow the Atlas thresholds and every row has its parameters", {
  h <- rd("era_site_season_hazards.csv")
  cut <- function(v, c) ifelse(v > c[3], "extreme", ifelse(v > c[2], "severe", ifelse(v > c[1], "moderate", "none")))
  expect_equal(h$NDWS_band,  cut(h$NDWS_monthly_mean,  c(15, 20, 25)))
  expect_equal(h$NDWL0_band, cut(h$NDWL0_monthly_mean, c(2, 5, 8)))
  expect_equal(h$NTx35_band, cut(h$NTx35_monthly_mean, c(7, 14, 21)))
  expect_equal(anyDuplicated(h[c("Site.Key", "season_slot", "year")]), 0)
  expect_true(all(c("soilcp_mm", "root_depth_m") %in% names(h)))
  expect_true(all(is.finite(h$soilcp_mm)))
})

test_that("08: every matched comparison is rainfed and joins one classified season", {
  s <- utils::read.csv(need(file.path("..", "..", "outputs", "af_formulations", "era_split", "era_split_sites.csv")))
  expect_false(any(s$irrigated %in% TRUE))
  expect_true(all(s$NDWS_band %in% c("none", "moderate", "severe", "extreme")))
  expect_true(all(abs(s$lrr - log(s$MeanT / s$MeanC)) < 1e-9))
  expect_true(all(abs(s$diff - (s$MeanT - s$MeanC)) < 1e-9))
})
