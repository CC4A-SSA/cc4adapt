# ------------------------------------------------------------------------------
# code/explore/af_formulations/R/formulations.R
#
# EXPLORATORY. Alternative ways of splitting a practice's stress-year benefit
# into a production part and an adaptation (avoided-loss) part, so the choice
# of formulation can be compared on Atlas hazard regimes.
# Brief: _claude/notes/2026-10-26 - BRIEF - Adaptation factor formulations ...
# Decisions: _claude/notes/2026-10-27 - DECISIONS - AF formulations session.md
#
# Notation (per hazard h, severity s, crop c, unit u):
#   B  baseline normal-year yield (t/ha)
#   g  practice's normal-year relative gain, A = B * (1 + g)
#   L  share of yield lost by the control under (h, s)
#   m  share of that loss the practice avoids (signed; < 0 = maladaptive)
#
#   Bs = B * (1 - L)                 control, stress year
#   As = A * (1 - L * (1 - m))       practice, stress year
#   dn = A - B                       normal-year gap
#   ds = As - Bs                     stress-year gap
#
# Formulations (stress-year production / avoided loss):
#   F0  v0.1, no split        A - B            / As - Bs      (double counts)
#   F1  levels DiD            A - B            / L(B - (1-m)A)
#   F2  DiD floored at 0      ds - max(0, F1)  / max(0, F1)
#   F3  proportional split    (A - B)(1 - L)   / A L m
#   F4  relative g x B        (A - B)(1 - L)   / (g_s - g_n) B = A L m / (1-L)
#   F5  current cc4adapt      A - B            / B L m
#       (step2 adds the normal-year gain every year; step3 avoided loss on
#        the control's yield)
#
# F1, F2 and F3 sum exactly to ds. F0, F4 and F5 do not; af_stress_year()
# returns the totals so the over-count is visible.
# ------------------------------------------------------------------------------

AF_FORMULATIONS <- c("F0", "F1", "F2", "F3", "F4", "F5")

# Gain at which levels DiD avoided loss is exactly zero: m* = g / (1 + g).
af_m_star <- function(g) g / (1 + g)

# Stress-year production and avoided loss for every formulation.
# All arguments are recycled; returns a data.frame with one row per input row.
af_stress_year <- function(B, g, L, m) {
  n  <- max(length(B), length(g), length(L), length(m))
  B  <- rep_len(B, n); g <- rep_len(g, n); L <- rep_len(L, n); m <- rep_len(m, n)
  if (any(L < 0 | L > 1, na.rm = TRUE)) stop("L must lie in [0, 1]")

  A  <- B * (1 + g)
  Bs <- B * (1 - L)
  As <- A * (1 - L * (1 - m))
  dn <- A - B
  ds <- As - Bs

  f1_al <- L * (B - (1 - m) * A)
  f2_al <- pmax(0, f1_al)
  f3_al <- A * L * m
  # F4 is undefined when the control loses everything (Bs = 0).
  f4_al <- ifelse(L < 1, A * L * m / (1 - L), NA_real_)

  out <- data.frame(
    B = B, g = g, L = L, m = m, A = A, Bs = Bs, As = As, dn = dn, ds = ds,
    F0_prod = dn,             F0_al = ds,
    F1_prod = dn,             F1_al = f1_al,
    F2_prod = ds - f2_al,     F2_al = f2_al,
    F3_prod = dn * (1 - L),   F3_al = f3_al,
    F4_prod = dn * (1 - L),   F4_al = f4_al,
    F5_prod = dn,             F5_al = B * L * m
  )
  for (f in AF_FORMULATIONS) {
    out[[paste0(f, "_total")]] <- out[[paste0(f, "_prod")]] + out[[paste0(f, "_al")]]
  }
  out
}

# Convert the Atlas's cumulative severity tiers (P(>= moderate) >= P(>= severe)
# >= P(>= extreme)) to exclusive bands. Valid for single-hazard marginals only.
# Small negative differences from the 3-decimal rounding of published
# frequencies are clipped to 0; larger ones are an error.
af_bands_from_cumulative <- function(p_mod, p_sev, p_ext, tol = 0.01) {
  b_mod <- p_mod - p_sev
  b_sev <- p_sev - p_ext
  bad <- (b_mod < -tol) | (b_sev < -tol)
  if (any(bad, na.rm = TRUE)) {
    stop("Cumulative tiers not monotone beyond tol in ", sum(bad, na.rm = TRUE), " rows")
  }
  data.frame(moderate = pmax(0, b_mod), severe = pmax(0, b_sev), extreme = p_ext)
}

# Expected annual production and avoided loss per hazard (brief §2, work plan
# Eq. 2). `bands` is a long data.frame with one row per (id, severity) and
# columns p_band, L, m; `id_cols` identify one expected value (e.g. unit,
# practice, hazard, parameter set). B and g must be constant within an id.
#   E[AL]   = sum_s P(s) AL(s)
#   E[prod] = (1 - sum_s P(s)) (A - B) + sum_s P(s) prod(s)
af_expected_annual <- function(bands, id_cols) {
  need <- c(id_cols, "B", "g", "p_band", "L", "m")
  miss <- setdiff(need, names(bands))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  sy  <- af_stress_year(bands$B, bands$g, bands$L, bands$m)
  key <- interaction(bands[id_cols], drop = TRUE, lex.order = TRUE)

  p_stress <- tapply(bands$p_band, key, sum)
  if (any(p_stress > 1 + 1e-9)) stop("Band probabilities sum above 1 for some ids")
  dn <- tapply(sy$dn, key, function(x) {
    if (diff(range(x)) > 1e-9) stop("B or g varies within an id")
    x[1]
  })

  out <- unique(bands[id_cols])
  out <- out[order(interaction(out[id_cols], drop = TRUE, lex.order = TRUE)), , drop = FALSE]
  k   <- as.character(interaction(out[id_cols], drop = TRUE, lex.order = TRUE))
  out$p_stress <- as.numeric(p_stress[k])
  out$dn       <- as.numeric(dn[k])

  for (f in AF_FORMULATIONS) {
    e_al <- tapply(bands$p_band * sy[[paste0(f, "_al")]], key, sum)
    e_pr <- tapply(bands$p_band * sy[[paste0(f, "_prod")]], key, sum)
    out[[paste0(f, "_al")]]    <- as.numeric(e_al[k])
    out[[paste0(f, "_prod")]]  <- (1 - out$p_stress) * out$dn + as.numeric(e_pr[k])
    out[[paste0(f, "_total")]] <- out[[paste0(f, "_prod")]] + out[[paste0(f, "_al")]]
    out[[paste0(f, "_ashare")]] <- ifelse(out[[paste0(f, "_total")]] != 0,
                                          out[[paste0(f, "_al")]] / out[[paste0(f, "_total")]],
                                          NA_real_)
  }
  rownames(out) <- NULL
  out
}
