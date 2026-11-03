# Tests for code/common/adaptation_split.R (F3 production/adaptation split)
# Run from the repo root:
#   Rscript -e 'testthat::test_dir("tests/testthat")'

source(normalizePath(file.path(testthat::test_path(), "..", "..",
                               "code", "common", "adaptation_split.R")))

grid <- expand.grid(B = c(0.8, 2, 10), g = c(-0.1, 0, 0.1, 0.25, 0.5, 1),
                    L = c(0, 0.1, 0.25, 0.5, 0.8, 1), m = c(-0.5, -0.2, 0, 0.25, 0.5, 1))

test_that("production + avoided loss equals the stress-year gap", {
  s <- cca_f3_split(grid$B, grid$g, grid$L, grid$m)
  A <- grid$B * (1 + grid$g)
  gap <- A * (1 - grid$L * (1 - grid$m)) - grid$B * (1 - grid$L)
  expect_equal(s$production + s$avoided_loss, gap)
  expect_equal(s$total, s$gap)
  expect_equal(s$production, (A - grid$B) * (1 - grid$L))
  expect_equal(s$avoided_loss, A * grid$L * grid$m)
})

test_that("g = 0 gives avoided loss B L m and zero production", {
  s <- cca_f3_split(grid$B, 0, grid$L, grid$m)
  expect_equal(s$avoided_loss, grid$B * grid$L * grid$m)
  expect_equal(s$production, rep(0, nrow(grid)))
})

test_that("avoided loss is larger than the old control-yield rule by (1 + g)", {
  s <- cca_f3_split(grid$B, grid$g, grid$L, grid$m)
  old <- grid$B * grid$L * grid$m
  expect_equal(s$avoided_loss, old * (1 + grid$g))
})

test_that("m < 0 gives negative avoided loss (maladaptation is recorded)", {
  s <- cca_f3_split(B = 2, g = 0.2, L = 0.4, m = c(-0.5, -0.1))
  expect_true(all(s$avoided_loss < 0))
  expect_equal(s$avoided_loss, 2 * 1.2 * 0.4 * c(-0.5, -0.1))
})

test_that("no hazard: production A - B, avoided loss 0", {
  s <- cca_f3_split(B = c(1, 3), g = c(0.3, -0.05), L = 0, m = 0.7)
  expect_equal(s$production, c(1, 3) * c(0.3, -0.05))
  expect_equal(s$avoided_loss, c(0, 0))
  e <- cca_expected_split(B = 3, g = 0.3, p_band = c(0, 0, 0),
                          L = c(0.1, 0.25, 0.45), m = c(0.2, 0.3, 0.4))
  expect_equal(e$production, 3 * 0.3)
  expect_equal(e$avoided_loss, 0)
  expect_equal(e$expected_L, 0)
})

test_that("expected split equals the probability-weighted F3 split", {
  B <- 2.5; g <- 0.2
  p <- c(0.15, 0.08, 0.03); L <- c(0.10, 0.25, 0.45); m <- c(0.3, 0.2, -0.1)
  e <- cca_expected_split(B, g, p, L, m)
  s <- cca_f3_split(B, g, L, m)
  dn <- B * g
  expect_equal(e$expected_L, sum(p * L))
  expect_equal(e$avoided_loss, sum(p * s$avoided_loss))
  expect_equal(e$production, (1 - sum(p)) * dn + sum(p * s$production))
  expect_equal(e$production, dn * (1 - sum(p * L)))
  # total = expected gap over normal and stress seasons
  expect_equal(e$total, (1 - sum(p)) * dn + sum(p * s$gap))
})

test_that("expected split is vectorised over units (matrix bands)", {
  P <- rbind(c(0.1, 0.05, 0.01), c(0.3, 0.1, 0.05))
  e <- cca_expected_split(B = c(1, 2), g = c(0.1, 0.4), p_band = P,
                          L = c(0.1, 0.25, 0.45), m = c(0.2, 0.2, 0.2))
  e1 <- cca_expected_split(1, 0.1, P[1, ], c(0.1, 0.25, 0.45), 0.2 + 0 * 1:3)
  e2 <- cca_expected_split(2, 0.4, P[2, ], c(0.1, 0.25, 0.45), 0.2 + 0 * 1:3)
  expect_equal(e, rbind(e1, e2))
})

test_that("matches the exploration's F3 reference implementation", {
  ref <- normalizePath(file.path(testthat::test_path(), "..", "..", "code", "explore",
                                 "af_formulations", "R", "formulations.R"),
                       mustWork = FALSE)
  skip_if_not(file.exists(ref), "exploration code not present")
  env <- new.env(); sys.source(ref, envir = env)
  sy <- env$af_stress_year(grid$B, grid$g, grid$L, grid$m)
  s  <- cca_f3_split(grid$B, grid$g, grid$L, grid$m)
  expect_equal(s$production, sy$F3_prod)
  expect_equal(s$avoided_loss, sy$F3_al)
})

test_that("combining hazards", {
  el <- list(c(0.1, 0), c(0.2, 0.5))
  expect_equal(cca_combine_expected_loss(el), c(1 - 0.9 * 0.8, 0.5))
  expect_equal(cca_combine_expected_loss(el, "additive"), c(0.3, 0.5))
  expect_equal(cca_combine_expected_loss(list(0.7, 0.6), "additive"), 1)
  expect_equal(cca_combine_expected_loss(list(c(0.2, 0.3))), c(0.2, 0.3))
  expect_error(cca_combine_expected_loss(list()))
})

test_that("arithmetic helpers work on SpatRasters", {
  skip_if_not_installed("terra")
  r <- function(v) terra::rast(matrix(v, 1, 2))
  A <- r(c(2, 3)); B <- r(c(1.5, 3))
  expect_equal(as.vector(terra::values(cca_f3_avoided_loss(A, 0.2, -0.5))), c(-0.2, -0.3))
  expect_equal(as.vector(terra::values(cca_f3_production(A, B, r(c(0.1, 0.5))))), c(0.45, 0))
  comb <- cca_combine_expected_loss(list(r(c(0.1, 0.6)), r(c(0.2, 0.6))), "additive")
  expect_equal(as.vector(terra::values(comb)), c(0.3, 1))
})

test_that("signed buffer: keeps negatives, bounds to [-1, 1] with a warning", {
  expect_equal(cca_signed_buffer(c(-0.4, 0, 0.7)), c(-0.4, 0, 0.7))
  expect_warning(v <- cca_signed_buffer(c(1.3, -2, 0.5), label = "test"),
                 "test: 2 value")
  expect_equal(v, c(1, -1, 0.5))
  expect_equal(cca_signed_buffer(c(NA, NaN, 0.2)), c(NA, NA, 0.2))
  expect_error(cca_signed_buffer("a"))
})

test_that("input checks: NA, lengths, ranges", {
  expect_error(cca_f3_split(B = c(1, 2, 3), g = c(0.1, 0.2), L = 0.1, m = 0.1),
               "length 1 or 3")
  expect_error(cca_f3_split(B = NA_real_, g = 0.1, L = 0.1, m = 0.1), "NA in B")
  expect_error(cca_f3_split(B = 1, g = 0.1, L = 0.1, m = NA_real_), "NA in m")
  expect_error(cca_f3_split(B = 1, g = 0.1, L = 1.2, m = 0.1), "L must lie")
  expect_error(cca_f3_split(B = -1, g = 0.1, L = 0.1, m = 0.1), "B must be")
  expect_error(cca_f3_split(B = 1, g = -1, L = 0.1, m = 0.1), "g must be")
  expect_error(cca_f3_split(B = "1", g = 0.1, L = 0.1, m = 0.1), "numeric")
  expect_error(cca_expected_split(1, 0.1, c(0.6, 0.5), c(0.1, 0.2), c(0, 0)),
               "sum above 1")
  expect_error(cca_expected_split(1, 0.1, c(-0.1, 0.5), c(0.1, 0.2), c(0, 0)),
               "Negative p_band")
  expect_error(cca_expected_split(1, 0.1, c(0.1, 0.2), c(0.1, 0.2, 0.3), c(0, 0)),
               "same number of bands")
  expect_error(cca_expected_split(1, 0.1, c(0.1, NA), c(0.1, 0.2), c(0, 0)),
               "NA in p_band")
  expect_error(cca_expected_split(c(1, 2, 3), 0.1, rbind(c(0.1, 0.2), c(0.1, 0.2)),
                                  c(0.1, 0.2), c(0, 0)), "rows")
})
