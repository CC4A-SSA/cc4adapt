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
# Shares use the absolute total: share = al / (|prod| + |al|), which stays in
# [-1, 1] when the total is near zero or negative (practices with g < 0).
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

IN  <- "data/staging/af_formulations"
OUT <- file.path(path.expand(Sys.getenv("AF_OUT_DIR", file.path(
  "~/Library/CloudStorage/OneDrive-CGIAR/Gates Carbon Credits for Adaptation",
  "_claude/outputs/adaptation_factor"))), "diagnostics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
REF_PROFILES <- c("const_0.50", "decl_mid")
FS <- c("F1", "F2", "F3", "F5")

cfg <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")
r <- as.data.table(arrow::read_parquet(file.path(IN, "af_results_maize_historic.parquet")))
grid <- fread(file.path(IN, "params", "param_grid.csv"))

# ---- E[L] and m_eff per unit x parameter set -----------------------------
pl <- melt(unique(r[, .(gaul2_code, hazard, p_moderate, p_severe, p_extreme)]),
           id.vars = c("gaul2_code", "hazard"), variable.name = "severity", value.name = "p")
pl[, severity := sub("^p_", "", severity)]
lm_ <- unique(grid[, .(param_id, hazard, severity, L, m)])
el <- merge(merge(r[, .(gaul2_code, param_id, hazard)], lm_, by = c("param_id", "hazard"),
                  allow.cartesian = TRUE), pl, by = c("gaul2_code", "hazard", "severity"))
el <- el[, .(EL = sum(p * L), m_eff = sum(p * L * m) / sum(p * L),
             ext_share = sum((p * L)[severity == "extreme"]) / sum(p * L)),
         by = .(gaul2_code, param_id)]
r <- merge(r, el, by = c("gaul2_code", "param_id"))

for (f in c("F0", FS)) {
  r[, (paste0(f, "_share")) := get(paste0(f, "_al")) /
      (abs(get(paste0(f, "_prod"))) + abs(get(paste0(f, "_al"))))]
}
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
  F2_zeroed = mean(F1_al < -1e-12),
  F3_negative = mean(F3_al < -1e-12)
), by = .(practice, hazard, g_source, g, m_profile)]
setorder(d1, practice, hazard, g, m_profile)
fwrite(d1, file.path(OUT, "d1_sign.csv"))

# ---- D2 gap scaling ------------------------------------------------------
ref <- r[g_source == "era" & m_profile %in% REF_PROFILES]
d2 <- ref[, .(iso3, admin1_name, admin2_name, gaul2_code, fs1, aez, bimodal,
              practice, hazard, m_profile, g, B, p_stress, EL, m_eff, ext_share,
              F1_al, F3_al, F5_al, gap, rel_gap)]
fwrite(d2, file.path(OUT, "d2_gap_units.csv"))

# log|gap| = log|g| + log B + log E[L]; within a parameter set g is fixed, so
# the spatial variance splits between log B and log E[L] (plus covariance).
var_split <- function(d) {
  d <- d[B > 0 & EL > 0]
  lb <- log(d$B); le <- log(d$EL); v <- var(lb + le)
  data.table(units = nrow(d), var_log_gap = v,
             share_B = var(lb) / v, share_EL = var(le) / v,
             share_cov = 2 * cov(lb, le) / v)
}
d2v <- unique(r[g_source == "era", .(gaul2_code, practice, hazard, B, EL, param_id)])[
  , var_split(.SD), by = .(practice, hazard, param_id)]
d2v <- merge(unique(grid[, .(param_id, m_profile)]), d2v, by = "param_id")[m_profile %in% REF_PROFILES]
fwrite(d2v, file.path(OUT, "d2_gap_variance.csv"))

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
fwrite(d3, file.path(OUT, "d3_label_flips.csv"))

# ---- D4 strata -----------------------------------------------------------
strata_sum <- function(d) d[, .(
  units = .N, maize_kt = sum(maize_prod_t) / 1e3,
  B_med = median(B), p_stress_med = median(p_stress), EL_med = median(EL),
  F3_al_med = median(F3_al), F1_al_med = median(F1_al), F5_al_med = median(F5_al),
  gap_med = median(gap), gap_wmean = weighted.mean(gap, maize_prod_t),
  F1_share_med = median(F1_share), F3_share_med = median(F3_share),
  F5_share_med = median(F5_share)
)]
era <- r[g_source == "era" & m_profile %in% REF_PROFILES]
by_vars <- c("practice", "hazard", "m_profile")
d4 <- rbind(
  era[, strata_sum(.SD), by = c(by_vars, "fs1", "aez")][, level := "fs1_x_aez"],
  era[, strata_sum(.SD), by = c(by_vars, "fs1")][, level := "fs1"],
  era[, strata_sum(.SD), by = c(by_vars, "aez")][, level := "aez"],
  era[, strata_sum(.SD), by = c(by_vars, "iso3")][, level := "country"],
  fill = TRUE
)
setcolorder(d4, c("level", by_vars, "iso3", "fs1", "aez"))
fwrite(d4, file.path(OUT, "d4_strata.csv"))

# ---- D5 severity mix (declining profiles only) ---------------------------
d5 <- r[g_source == "era" & grepl("^decl", m_profile),
        .(units = .N,
          cor_ext_share_m_eff = cor(ext_share, m_eff, use = "complete.obs"),
          m_eff_q10 = quantile(m_eff, 0.1, na.rm = TRUE), m_eff_q90 = quantile(m_eff, 0.9, na.rm = TRUE),
          rel_gap_q10 = quantile(rel_gap, 0.1, na.rm = TRUE),
          rel_gap_q90 = quantile(rel_gap, 0.9, na.rm = TRUE)),
        by = .(practice, hazard, m_profile, g)]
fwrite(d5, file.path(OUT, "d5_severity_mix.csv"))

# ---- D6 bimodal ----------------------------------------------------------
d6 <- era[, strata_sum(.SD), by = c(by_vars, "bimodal")]
fwrite(d6, file.path(OUT, "d6_bimodal.csv"))

# ---- Console summary -----------------------------------------------------
cat("\nGap identity holds: max |gap - (-g B E[L])| =", signif(chk, 2), "\n")
cat("\nD1 (ERA g): share of units with F1 < 0 / F1 > 0 without buffering\n")
print(d1[g_source == "era", .(practice, hazard, m_profile, m_star, F1_negative,
                              F1_positive_without_buffering)], digits = 2)
cat("\nD2: spatial variance of log|gap| (ERA g, reference profiles)\n")
print(d2v[, .(practice, hazard, m_profile, units, share_B, share_EL, share_cov)], digits = 2)
cat("\nD3 (ERA g): label flips vs F3\n")
print(d3[g_source == "era" & m_profile %in% c(REF_PROFILES, "const_0.25", "decl_lo"),
         .(practice, hazard, m_profile, F3_adaptation, flip_F1_vs_F3, flip_F2_vs_F3, flip_F5_vs_F3)],
      digits = 2)
cat("\nD5: severity mix vs effective m (declining profiles, ERA g)\n")
print(d5, digits = 2)
message("Wrote diagnostics to ", OUT)
