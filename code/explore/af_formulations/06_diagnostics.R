# ------------------------------------------------------------------------------
# code/explore/af_formulations/06_diagnostics.R
#
# EXPLORATORY. Step 5.5 diagnostics on the F0-F5 results.
#
# Run from the repo root:  Rscript code/explore/af_formulations/06_diagnostics.R
#
# Closed forms (per stress year; tested in tests/testthat):
#   F1_al - F3_al = -L B g                (independent of m)
#   F5_al - F3_al = -L m B g
#   (F1_al - F3_al) / F3_al = -g / ((1 + g) m)
# so the expected annual gap is -g B E[L], E[L] = sum_s P(s) L(s), and for a
# constant m the relative gap is the same in every unit. Location (B, hazard
# frequency, severity mix) scales the absolute gap; with m varying by
# severity it also moves the relative gap through the effective m
#   m_eff = sum_s P(s) L(s) m(s) / E[L].
#
# Shares are the robust F*_share from 05 (al / (|prod| + |al|), in [-1, 1]).
# m_eff and ext_share are NA where E[L] = 0 (no hazard); strata medians use
# na.rm. Strata tables carry p_source so resolved and fallback units can be
# told apart, and exclude p_stress >= 0.99 pairs from the "identified" rows.
#
# Outputs (<AF_OUT_DIR>/diagnostics/):
#   d1_sign.csv          F1 sign problems by practice x hazard x g x m profile
#   d2_gap_units.csv     per-unit gap, B, E[L], m_eff (ERA g; reference profiles)
#   d2_gap_variance.csv  share of spatial variance in log|gap| from B vs E[L]
#   d3_label_flips.csv   units whose dominant label differs from F3
#   d4_strata.csv        FS1 x AEZ, FS1 and AEZ summaries (ERA g)
#   d5_severity_mix.csv  relative gap vs extreme share of E[L], declining profiles
#   d6_bimodal.csv       bimodal vs other units
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}

source("code/explore/af_formulations/R/paths.R")
IN  <- "data/staging/af_formulations"
OUT <- file.path(af_out_dir(), "diagnostics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REC <- af_record_dir()
write_diag <- function(d, name) {
  fwrite(d, file.path(OUT, name))
  if (!is.null(REC)) { dir.create(file.path(REC, "diagnostics"), showWarnings = FALSE); fwrite(d, file.path(REC, "diagnostics", name)) }
}
REF_PROFILES <- c("const_0.50", "decl_mid")
FS <- c("F1", "F2", "F3", "F5")

cfg <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")
r <- as.data.table(arrow::read_parquet(file.path(af_out_dir(), "af_results_maize_historic.parquet")))
r[, B := B_norm]
grid <- fread(file.path(IN, "params", "param_grid.csv"))

# ---- E[L] and m_eff per unit x parameter set -----------------------------
pl <- melt(unique(r[, .(gaul2_code, hazard, p_moderate, p_severe, p_extreme)]),
           id.vars = c("gaul2_code", "hazard"), variable.name = "severity", value.name = "p")
pl[, severity := sub("^p_", "", severity)]
# Join L and m on (hazard, m_profile), not param_id: the era_aez sets built in
# 05_compute carry ids outside param_grid.csv.
lm_ <- unique(grid[, .(hazard, m_profile, severity, L, m)])
el <- merge(merge(r[, .(gaul2_code, param_id, hazard, m_profile)], lm_, by = c("hazard", "m_profile"),
                  allow.cartesian = TRUE), pl, by = c("gaul2_code", "hazard", "severity"))
stopifnot(nrow(el) == 3 * nrow(r))
el <- el[, .(EL_chk = sum(p * L),
             m_eff = fifelse(sum(p * L) > 0, sum(p * L * m) / sum(p * L), NA_real_),
             ext_share = fifelse(sum(p * L) > 0, sum((p * L)[severity == "extreme"]) / sum(p * L), NA_real_)),
         by = .(gaul2_code, param_id)]
r <- merge(r, el, by = c("gaul2_code", "param_id"))
if (r[, max(abs(EL_chk - EL))] > 1e-9) stop("E[L] recomputed here differs from 05_compute")
r[, gap := F1_al - F3_al]
r[, gap_pred := -g * B * EL]
r[, rel_gap := fifelse(F3_al != 0, gap / F3_al, NA_real_)]
chk <- r[, max(abs(gap - gap_pred))]
if (chk > 1e-9) stop("Expected gap differs from -g B E[L] by ", signif(chk, 3))

# ---- D1 sign problems ----------------------------------------------------
d1 <- r[, .(
  units = .N,
  m_star = round(g[1] / (1 + g[1]), 3),
  F1_negative = mean(F1_al < -1e-12),
  F1_positive_without_buffering = mean(F1_al > 1e-12 & m_eff <= 0, na.rm = TRUE),
  F3_negative = mean(F3_al < -1e-12),
  p_stress_ge_099 = mean(p_stress_ge_099)
), by = .(practice, hazard, g_source, g, m_profile)]
setorder(d1, practice, hazard, g, m_profile)
write_diag(d1, "d1_sign.csv")

# ---- D2 gap scaling ------------------------------------------------------
ref <- r[g_source == "era" & m_profile %in% REF_PROFILES]
d2 <- ref[, .(iso3, admin1_name, admin2_name, gaul2_code, fs1, aez, bimodal, p_source, p_stress_ge_099,
              practice, hazard, m_profile, g, B_obs, B_norm, p_stress, EL, m_eff, ext_share,
              F1_al, F3_al, F5_al, gap, rel_gap)]
write_diag(d2, "d2_gap_units.csv")

# log|gap| = log|g| + log B + log E[L]; within a parameter set g is fixed, so
# the spatial variance splits between log B and log E[L] (plus covariance).
var_split <- function(d) {
  d <- d[B > 0 & EL > 0]
  if (nrow(d) < 2) return(data.table(units = nrow(d), var_log_gap = NA_real_, share_B = NA_real_, share_EL = NA_real_, share_cov = NA_real_))
  lb <- log(d$B); le <- log(d$EL); v <- var(lb + le)
  if (!is.finite(v) || v == 0) return(data.table(units = nrow(d), var_log_gap = v, share_B = NA_real_, share_EL = NA_real_, share_cov = NA_real_))
  data.table(units = nrow(d), var_log_gap = v,
             share_B = var(lb) / v, share_EL = var(le) / v,
             share_cov = 2 * cov(lb, le) / v)
}
d2v <- unique(r[g_source == "era", .(gaul2_code, practice, hazard, B, EL, param_id)])[
  , var_split(.SD), by = .(practice, hazard, param_id)]
d2v <- merge(unique(grid[, .(param_id, m_profile)]), d2v, by = "param_id")[m_profile %in% REF_PROFILES]
write_diag(d2v, "d2_gap_variance.csv")

# ---- D3 label flips ------------------------------------------------------
# Dominant label: "adaptation" if |al| > |prod|, else "production".
lab <- function(al, prod) fifelse(abs(al) > abs(prod), "adaptation", "production")
for (f in FS) r[, (paste0(f, "_label")) := lab(get(paste0(f, "_al")), get(paste0(f, "_prod")))]
d3 <- r[, .(
  units = .N,
  F3_adaptation = mean(F3_label == "adaptation"),
  flip_F1_vs_F3 = mean(F1_label != F3_label),
  flip_F2_vs_F3 = mean(F2_label != F3_label),
  flip_F5_vs_F3 = mean(F5_label != F3_label)
), by = .(practice, hazard, g_source, g, m_profile)]
setorder(d3, practice, hazard, g, m_profile)
write_diag(d3, "d3_label_flips.csv")

# ---- D4 strata -----------------------------------------------------------
strata_sum <- function(d) d[, .(
  units = .N, units_adm1_fallback = sum(p_source == "adm1"),
  units_p_stress_ge_099 = sum(p_stress_ge_099), maize_kt = sum(maize_prod_t) / 1e3,
  B_obs_med = median(B_obs), B_norm_med = median(B_norm),
  p_stress_med = median(p_stress), EL_med = median(EL),
  F3_al_med = median(F3_al), F1_al_med = median(F1_al), F5_al_med = median(F5_al),
  gap_med = median(gap), gap_wmean = weighted.mean(gap, maize_prod_t),
  F1_share_med = median(F1_share, na.rm = TRUE), F3_share_med = median(F3_share, na.rm = TRUE),
  F5_share_med = median(F5_share, na.rm = TRUE)
)]
era <- r[g_source %in% c("era", "era_aez", "era_norm") & m_profile %in% REF_PROFILES]
by_vars <- c("practice", "hazard", "g_source", "m_profile")
d4 <- rbind(
  era[, strata_sum(.SD), by = c(by_vars, "fs1", "aez")][, level := "fs1_x_aez"],
  era[, strata_sum(.SD), by = c(by_vars, "fs1")][, level := "fs1"],
  era[, strata_sum(.SD), by = c(by_vars, "aez")][, level := "aez"],
  era[, strata_sum(.SD), by = c(by_vars, "iso3")][, level := "country"],
  fill = TRUE
)
setcolorder(d4, c("level", by_vars, "iso3", "fs1", "aez"))
write_diag(d4, "d4_strata.csv")

# ---- D5 severity mix (declining profiles only) ---------------------------
d5 <- r[g_source == "era" & grepl("^decl", m_profile) & !is.na(m_eff),
        .(units = .N,
          cor_ext_share_m_eff = cor(ext_share, m_eff, use = "complete.obs"),
          m_eff_q10 = quantile(m_eff, 0.1, na.rm = TRUE), m_eff_q90 = quantile(m_eff, 0.9, na.rm = TRUE),
          rel_gap_q10 = quantile(rel_gap, 0.1, na.rm = TRUE),
          rel_gap_q90 = quantile(rel_gap, 0.9, na.rm = TRUE)),
        by = .(practice, hazard, m_profile, g)]
write_diag(d5, "d5_severity_mix.csv")

# ---- D6 bimodal ----------------------------------------------------------
d6 <- era[, strata_sum(.SD), by = c(by_vars, "bimodal")]
write_diag(d6, "d6_bimodal.csv")

# ---- Console summary -----------------------------------------------------
cat("\nGap identity holds: max |gap - (-g B E[L])| =", signif(chk, 2), "\n")
cat("\nD1 (ERA g): share of units with F1 < 0 / F1 > 0 without buffering\n")
print(d1[g_source == "era", .(practice, hazard, m_profile, m_star, F1_negative,
                              F1_positive_without_buffering, p_stress_ge_099)], digits = 2)
cat("\nD2: spatial variance of log|gap| (ERA g, reference profiles)\n")
print(d2v[, .(practice, hazard, m_profile, units, share_B, share_EL, share_cov)], digits = 2)
cat("\nD3 (ERA g): label flips vs F3\n")
print(d3[g_source == "era" & m_profile %in% c(REF_PROFILES, "const_0.25", "decl_lo"),
         .(practice, hazard, m_profile, F3_adaptation, flip_F1_vs_F3, flip_F2_vs_F3, flip_F5_vs_F3)],
      digits = 2)
cat("\nD5: severity mix vs effective m (declining profiles, ERA g)\n")
print(d5, digits = 2)
message("Wrote diagnostics to ", OUT)
