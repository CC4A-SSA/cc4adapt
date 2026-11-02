# Tests for code/common/era_yield_effects.R (ERA yield effect shared by
# Steps 2 and 3). Run from the repo root:
#   Rscript -e 'testthat::test_dir("tests/testthat")'

source(normalizePath(file.path(testthat::test_path(), "..", "..",
                               "code", "common", "era_yield_effects.R")))

era <- data.frame(
  PrName     = c("Agroforestry", "Alleycropping-Intercropping", "Agroforestry",
                 "Biochar", "Agroforestry", "Agroforestry"),
  Product    = c("Maize", "maize", "Sorghum", "Maize", "Maize", "Maize"),
  Out.SubInd = c("Crop Yield", "Crop Yield", "Crop Yield", "Crop Yield",
                 "Soil Organic Carbon", "Crop Yield"),
  MeanT      = c(1.2, 2.0, 1.0, 3.0, 9.0, 0),
  MeanC      = c(1.0, 1.6, 1.0, 2.0, 1.0, 1.0),
  stringsAsFactors = FALSE
)

test_that("yield rows keep yield outcomes with positive means", {
  y <- cca_era_yield_rows(era)
  expect_equal(nrow(y), 4)
  expect_equal(y$RR, log(c(1.2, 2.0, 1.0, 3.0) / c(1.0, 1.6, 1.0, 2.0)))
})

test_that("g = exp(mean RR) - 1 over matching practice x crop rows", {
  y <- cca_era_yield_rows(era)
  eff <- cca_era_yield_effect(y, "agroforestry_cropland", "maiz")
  rr <- log(c(1.2, 2.0 / 1.6))
  expect_equal(eff$pc, exp(mean(rr)) - 1)
  expect_equal(eff$n_obs, 2)
  expect_lt(eff$ci_low, eff$pc)
  expect_gt(eff$ci_high, eff$pc)
})

test_that("no evidence returns NULL", {
  y <- cca_era_yield_rows(era)
  expect_null(cca_era_yield_effect(y, "agroforestry_cropland", "rice"))
  expect_null(cca_era_yield_effect(y, "not_a_practice", "maiz"))
  expect_null(cca_era_yield_effect(NULL, "biochar", "maiz"))
})
