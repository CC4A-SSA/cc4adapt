# ------------------------------------------------------------------------------
# code/explore/af_formulations/05_compute.R
#
# EXPLORATORY. Step 5.4: expected annual production and avoided loss under
# formulations F0-F5 for every maize admin-2 unit x parameter set, per hazard.
# Units: t/ha/yr. Price is deferred (it scales totals, not the split).
#
# Run from the repo root:  Rscript code/explore/af_formulations/05_compute.R
#
# P(h, s) per unit:
#   * Single-hazard marginal per cumulative tier from the Atlas exclusive
#     classes, e.g. P(dry >= s) = [dry + dry+heat + dry+wet + dry+heat+wet]
#     / [any + none]  (hazard briefing S1.7; production-weighted).
#   * Tiers forced monotone (P(>=sev) <= P(>=mod), P(>=ext) <= P(>=sev)) with
#     a flag when the adjustment exceeds 0.01, then differenced to bands.
#   * Admin-2 row used when present, not NaN and any + none >= MIN_TIER_T;
#     otherwise the unit's admin-1 row (p_source records which).
# B per unit from 04_units.R (MapSPAM maize production / harvested area).
#
# Outputs:
#   data/staging/af_formulations/af_results_maize_historic.parquet
#   <AF_OUT_DIR>/af_results_maize_historic.parquet (copy for the _claude record)
#   AF_OUT_DIR defaults to the project's _claude/outputs/adaptation_factor.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}
source("code/explore/af_formulations/R/formulations.R")

IN  <- "data/staging/af_formulations"
OUT <- Sys.getenv("AF_OUT_DIR", file.path(
  "~/Library/CloudStorage/OneDrive-CGIAR/Gates Carbon Credits for Adaptation",
  "_claude/outputs/adaptation_factor"))
MIN_TIER_T <- 100
SEVS <- c("moderate", "severe", "extreme")
HAZ_CLASSES <- list(
  dry  = c("dry", "dry+heat", "dry+wet", "dry+heat+wet"),
  heat = c("heat", "dry+heat", "heat+wet", "dry+heat+wet"),
  wet  = c("wet", "dry+wet", "heat+wet", "dry+heat+wet")
)

cfg <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")

# ---- P(h, >= s) per admin unit -------------------------------------------
tiers <- rbindlist(lapply(SEVS, function(s) {
  as.data.table(arrow::read_parquet(file.path(IN, sprintf("atlas_prod_t_jagermeyr_%s.parquet", s))))
}))
tiers <- tiers[crop == cfg$crop & hazard_vars == cfg$hazard_vars & scenario == cfg$scenario]
tiers[, level := fifelse(!is.na(admin2_name), "adm2", fifelse(!is.na(admin1_name), "adm1", "adm0"))]
tiers <- tiers[level != "adm0"]
tiers[, code := fifelse(level == "adm2", gaul2_code, gaul1_code)]

w <- dcast(tiers, level + code + severity ~ hazard, value.var = "value")
w[, tier_t := any + none]
cum <- rbindlist(lapply(names(HAZ_CLASSES), function(h) {
  w[, .(level, code, severity, tier_t, hazard = h,
        p_cum = rowSums(.SD) / tier_t), .SDcols = HAZ_CLASSES[[h]]]
}))
cum <- cum[is.finite(p_cum)]
# Per-pixel integer rounding in the Atlas tiers lets a few admin-2 ratios
# exceed 1 (by < 1 % in the 2026-10-19 tiers); clip to [0, 1].
cum[, p_cum := pmin(pmax(p_cum, 0), 1)]
cum <- dcast(cum, level + code + hazard ~ severity, value.var = c("p_cum", "tier_t"))
cum[, tier_t := tier_t_moderate]
cum[, mono_adj := pmax(0, p_cum_severe - p_cum_moderate, p_cum_extreme - pmin(p_cum_severe, p_cum_moderate))]
cum[, p_cum_severe  := pmin(p_cum_severe, p_cum_moderate)]
cum[, p_cum_extreme := pmin(p_cum_extreme, p_cum_severe)]
bands <- af_bands_from_cumulative(cum$p_cum_moderate, cum$p_cum_severe, cum$p_cum_extreme)
cum[, c("p_moderate", "p_severe", "p_extreme") := bands]

# ---- Units and their P source --------------------------------------------
units <- as.data.table(arrow::read_parquet(file.path(IN, "units.parquet")))
units <- units[maize_prod_t > 0 & is.finite(B)]

a2 <- cum[level == "adm2" & tier_t >= MIN_TIER_T]
a1 <- cum[level == "adm1"]
p_keys <- c("hazard", "p_moderate", "p_severe", "p_extreme", "mono_adj")
up <- rbindlist(lapply(names(HAZ_CLASSES), function(h) {
  x <- merge(units[, .(gaul2_code, gaul1_code)], a2[hazard == h, c("code", p_keys), with = FALSE],
             by.x = "gaul2_code", by.y = "code", all.x = TRUE)
  x[, p_source := fifelse(is.na(hazard), NA_character_, "adm2")]
  miss <- x[is.na(hazard), .(gaul2_code, gaul1_code)]
  y <- merge(miss, a1[hazard == h, c("code", p_keys), with = FALSE],
             by.x = "gaul1_code", by.y = "code", all.x = TRUE)
  y[, p_source := fifelse(is.na(hazard), NA_character_, "adm1")]
  rbind(x[!is.na(hazard)], y[!is.na(hazard)], fill = TRUE)
}))
message(sprintf("P source (unit x hazard): %s; units with no P for some hazard: %d",
                paste(names(table(up$p_source)), table(up$p_source), collapse = ", "),
                uniqueN(setdiff(units$gaul2_code, up[, .N, by = gaul2_code][N == 3, gaul2_code]))))

# ---- Expand to unit x parameter set x severity ---------------------------
grid <- fread(file.path(IN, "params", "param_grid.csv"))
pl <- melt(up, id.vars = c("gaul2_code", "hazard", "p_source", "mono_adj"),
           measure.vars = c("p_moderate", "p_severe", "p_extreme"),
           variable.name = "severity", value.name = "p_band")
pl[, severity := sub("^p_", "", severity)]

long <- merge(grid, pl, by = c("hazard", "severity"), allow.cartesian = TRUE)
long <- merge(long, units[, .(gaul2_code, B)], by = "gaul2_code")
setorder(long, gaul2_code, param_id, severity)

ea <- as.data.table(af_expected_annual(as.data.frame(long), c("gaul2_code", "param_id")))

# Identity check: F1, F2, F3 must give the same expected total.
d12 <- max(abs(ea$F1_total - ea$F3_total)); d23 <- max(abs(ea$F2_total - ea$F3_total))
if (max(d12, d23) > 1e-9) stop("F1/F2/F3 totals differ: ", signif(max(d12, d23), 3))

# ---- Attach parameters, P and unit attributes -----------------------------
pw <- dcast(pl, gaul2_code + hazard + p_source + mono_adj ~ severity, value.var = "p_band")
setnames(pw, SEVS, paste0("p_", SEVS))
pars <- unique(grid[, .(param_id, practice, hazard, g, g_source, m_profile)])
res <- merge(ea, pars, by = "param_id")
res <- merge(res, pw, by = c("gaul2_code", "hazard"))
res <- merge(units[, .(gaul2_code, gaul1_code, iso3, admin1_name, admin2_name,
                       maize_prod_t, maize_harv_ha, B, fs1, fs2, aez, bimodal, bimodal_share)],
             res, by = "gaul2_code")
res[, `:=`(F1_minus_F3_al = F1_al - F3_al, F1_negative = F1_al < 0)]
setcolorder(res, c("iso3", "admin1_name", "admin2_name", "gaul1_code", "gaul2_code",
                   "practice", "hazard", "g_source", "g", "m_profile", "param_id"))

f <- "af_results_maize_historic.parquet"
arrow::write_parquet(res, file.path(IN, f))
dir.create(path.expand(OUT), recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(res, file.path(path.expand(OUT), f))
message(sprintf("%s: %d rows (%d units x %d parameter sets); copied to %s",
                f, nrow(res), uniqueN(res$gaul2_code), uniqueN(res$param_id), path.expand(OUT)))
