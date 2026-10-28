# Tests for code/explore/af_formulations/R/formulations.R
# Run: Rscript -e 'testthat::test_dir("code/explore/af_formulations/tests/testthat")'

src <- normalizePath(file.path(testthat::test_path(), "..", "..", "R", "formulations.R"))
source(src)

grid <- expand.grid(B = c(0.8, 2, 10), g = c(-0.1, 0, 0.1, 0.25, 0.5, 1),
                    L = c(0, 0.1, 0.25, 0.5, 0.8), m = c(-0.2, 0, 0.25, 0.5, 0.75, 1))

test_that("reproduces the 26 Oct worked table", {
  wt <- read.csv(testthat::test_path("fixtures", "worked_table_2026-10-26.csv"))
  sy <- af_stress_year(wt$Bn, wt$An / wt$Bn - 1, wt$L, wt$m)
  expect_equal(sy$Bs, wt$Bs)
  expect_equal(sy$As, wt$As)
  expect_equal(sy$dn, wt$delta_normal)
  expect_equal(sy$ds, wt$delta_stress)
  expect_equal(sy$F1_al, wt$DiD_levels)
  expect_equal(sy$F2_prod, wt$floor_prod)
  expect_equal(sy$F2_al, wt$floor_adapt)
  expect_equal(sy$F3_prod, wt$prop_prod)
  expect_equal(sy$F3_al, wt$prop_adapt)
})

test_that("F0, F1, F2, F3, F5 agree at g = 0; F4 = F3 / (1 - L)", {
  sy <- af_stress_year(grid$B, 0, grid$L, grid$m)
  pos <- sy$m >= 0  # F2 floors negatives, so it only agrees for m >= 0
  expect_equal(sy$F1_al, sy$F3_al)
  expect_equal(sy$F5_al, sy$F3_al)
  expect_equal(sy$F0_al, sy$F3_al)
  expect_equal(sy$F2_al[pos], sy$F3_al[pos])
  ok <- sy$L < 1
  expect_equal(sy$F4_al[ok], sy$F3_al[ok] / (1 - sy$L[ok]))
})

test_that("levels DiD is zero at m = m*", {
  g  <- c(0.1, 0.25, 0.5, 1)
  sy <- af_stress_year(2, g, 0.5, af_m_star(g))
  expect_equal(sy$F1_al, rep(0, length(g)))
  expect_equal(af_m_star(c(0.5, 1)), c(1 / 3, 0.5))
})

test_that("F1, F2, F3 sum exactly to the stress-year gap", {
  sy <- af_stress_year(grid$B, grid$g, grid$L, grid$m)
  for (f in c("F1", "F2", "F3")) expect_equal(sy[[paste0(f, "_total")]], sy$ds)
})

test_that("F0, F4, F5 over-counts have the closed forms in the decisions note", {
  sy <- af_stress_year(grid$B, grid$g, grid$L, grid$m)
  expect_equal(sy$F0_total - sy$ds, sy$dn)
  expect_equal(sy$F5_total - sy$ds, sy$L * (sy$A - sy$B) * (1 - sy$m))
  ok <- sy$L < 1
  expect_equal((sy$F4_total - sy$ds)[ok],
               (sy$A * sy$L^2 * sy$m / (1 - sy$L))[ok])
})

test_that("F1 - F3 gap is -L B g, independent of m; F5 - F3 is m times that", {
  sy <- af_stress_year(grid$B, grid$g, grid$L, grid$m)
  expect_equal(sy$F1_al - sy$F3_al, -sy$L * sy$B * sy$g)
  expect_equal(sy$F5_al - sy$F3_al, -sy$L * sy$m * sy$B * sy$g)
  ok <- sy$F3_al != 0
  expect_equal(((sy$F1_al - sy$F3_al) / sy$F3_al)[ok],
               (-sy$g / ((1 + sy$g) * sy$m))[ok])
})

test_that("sign behaviour: F2 never negative, F3 negative only when m < 0", {
  sy <- af_stress_year(grid$B, grid$g, grid$L, grid$m)
  expect_true(all(sy$F2_al >= 0))
  nz <- sy$L > 0
  expect_identical(sy$F3_al[nz] < 0, sy$m[nz] < 0)
})

test_that("cumulative tiers convert to bands", {
  b <- af_bands_from_cumulative(c(0.5, 0.2), c(0.3, 0.2), c(0.1, 0.05))
  expect_equal(b$moderate, c(0.2, 0))
  expect_equal(b$severe, c(0.2, 0.15))
  expect_equal(rowSums(b), c(0.5, 0.2))
  expect_equal(af_bands_from_cumulative(0.300, 0.305, 0.1)$moderate, 0)  # rounding clip
  expect_error(af_bands_from_cumulative(0.2, 0.4, 0.1), "not monotone")
})

test_that("expected annual values", {
  bands <- data.frame(
    unit = rep(c("u1", "u2"), each = 3),
    severity = rep(c("moderate", "severe", "extreme"), 2),
    p_band = c(0.2, 0.1, 0.05, 0, 0, 0),
    B = 2, g = 0.5, L = rep(c(0.10, 0.25, 0.45), 2), m = rep(c(0.75, 0.5, 0.25), 2)
  )
  ea <- af_expected_annual(bands, "unit")

  # No hazard: everything is production, nothing is avoided loss.
  u2 <- ea[ea$unit == "u2", ]
  for (f in AF_FORMULATIONS) {
    expect_equal(u2[[paste0(f, "_prod")]], 1)
    expect_equal(u2[[paste0(f, "_al")]], 0)
  }

  # F1, F2, F3 give the same expected total; only the split differs.
  u1 <- ea[ea$unit == "u1", ]
  expect_equal(u1$F1_total, u1$F3_total)
  expect_equal(u1$F2_total, u1$F3_total)
  sy <- af_stress_year(2, 0.5, bands$L[1:3], bands$m[1:3])
  expect_equal(u1$F3_al, sum(bands$p_band[1:3] * sy$F3_al))
  expect_equal(u1$p_stress, 0.35)

  expect_error(af_expected_annual(transform(bands, p_band = 0.5), "unit"), "above 1")
  expect_error(af_expected_annual(transform(bands, g = seq_len(6) / 10), "unit"), "varies")
})
