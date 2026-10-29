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
#   * Tiers clipped to [0, 1], forced monotone (P(>=sev) <= P(>=mod),
#     P(>=ext) <= P(>=sev)); mono_adj records the size of that adjustment.
#     Then differenced to bands.
#   * Admin-2 row used when present, complete across tiers, not NaN and
#     any + none >= MIN_TIER_T; otherwise the unit's admin-1 row under the same
#     rules (p_source records which). Units with no usable P for a hazard the
#     grid needs are written to units_no_p.csv and excluded.
#   * p_none = P(no hazard at all at the moderate tier) and p_dry_wet =
#     P(dry and wet co-occur, moderate tier), from the exclusive classes, are
#     carried per unit for the co-occurrence discussion. They are not used in
#     the per-hazard formulas, whose production term treats seasons without
#     hazard h as normal years even if another hazard occurred (reported per
#     hazard; not summable across hazards).
#
# Baseline yield. B_obs = MapSPAM production / harvested area is a mean over
# stressed and normal seasons. Under the model's own assumptions
# B_obs = B_norm (1 - E[L]), E[L] = sum_s P(s) L(s), so the formulations use
# B_norm = B_obs / (1 - E[L]) (science review, 29 Oct 2026). E[L] depends on
# the placeholder loss table, so this correction inherits its uncertainty;
# both B_obs and B_norm are kept.
#
# g sources: "grid" (scenario values), "era" (ERA all-SSA median per practice,
# pooled across seasons), "era_aez" (ERA median for the unit's AEZ16simple
# stratum where it rests on at least MIN_STRATUM_STUDIES studies, else the
# all-SSA value; g_stratum_n says which), and "era_norm" (normal-season g from
# the ERA season split, 08_era_split.R: dry hazard, moderate threshold,
# study-averaged mean; present only when era_split_g.csv exists).
#
# Outputs (AF_OUT_DIR, default data/outputs/af_formulations; copy to
# AF_RECORD_DIR when set):
#   af_results_maize_historic.parquet, units_no_p.csv
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}
source("code/explore/af_formulations/R/formulations.R")
source("code/explore/af_formulations/R/paths.R")

IN  <- "data/staging/af_formulations"
MIN_TIER_T <- 100
MIN_STRATUM_STUDIES <- 3
P_ONE <- 0.99
SEVS <- c("moderate", "severe", "extreme")
HAZ_CLASSES <- list(
  dry  = c("dry", "dry+heat", "dry+wet", "dry+heat+wet"),
  heat = c("heat", "dry+heat", "heat+wet", "dry+heat+wet"),
  wet  = c("wet", "dry+wet", "heat+wet", "dry+heat+wet")
)

cfg  <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")
grid <- fread(file.path(IN, "params", "param_grid.csv"))
era_g <- fread(file.path(IN, "params", "era_g.csv"))
HAZ_USED <- sort(unique(grid$hazard))

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
w[, p_none := none / tier_t]
w[, p_dry_wet := (`dry+wet` + `dry+heat+wet`) / tier_t]
cum <- rbindlist(lapply(names(HAZ_CLASSES), function(h) {
  w[, .(level, code, severity, tier_t, p_none, p_dry_wet, hazard = h,
        p_cum = rowSums(.SD) / tier_t), .SDcols = HAZ_CLASSES[[h]]]
}))
# Per-pixel integer rounding in the Atlas tiers lets a few admin-2 ratios
# exceed 1 (by < 1 % in the 2026-10-19 tiers); clip to [0, 1].
cum[, p_cum := pmin(pmax(p_cum, 0), 1)]
cum <- dcast(cum, level + code + hazard ~ severity, value.var = c("p_cum", "tier_t", "p_none", "p_dry_wet"))
cum <- cum[is.finite(p_cum_moderate) & is.finite(p_cum_severe) & is.finite(p_cum_extreme) &
             is.finite(tier_t_moderate) & tier_t_moderate >= MIN_TIER_T]
cum[, `:=`(tier_t = tier_t_moderate, p_none = p_none_moderate, p_dry_wet = p_dry_wet_moderate)]
cum[, mono_adj := pmax(0, p_cum_severe - p_cum_moderate, p_cum_extreme - pmin(p_cum_severe, p_cum_moderate))]
cum[, p_cum_severe  := pmin(p_cum_severe, p_cum_moderate)]
cum[, p_cum_extreme := pmin(p_cum_extreme, p_cum_severe)]
bands <- af_bands_from_cumulative(cum$p_cum_moderate, cum$p_cum_severe, cum$p_cum_extreme)
cum[, c("p_moderate", "p_severe", "p_extreme") := bands]
stopifnot(!anyNA(cum[, .(p_moderate, p_severe, p_extreme)]))

# ---- Units and their P source --------------------------------------------
units <- as.data.table(arrow::read_parquet(file.path(IN, "units.parquet")))
units <- units[maize_prod_t > 0 & is.finite(B)]

p_cols <- c("p_moderate", "p_severe", "p_extreme", "mono_adj", "p_none", "p_dry_wet", "tier_t")
pick <- function(h) {
  a2 <- cum[level == "adm2" & hazard == h, c("code", p_cols), with = FALSE]
  a1 <- cum[level == "adm1" & hazard == h, c("code", p_cols), with = FALSE]
  x <- merge(units[, .(gaul2_code, gaul1_code)], a2, by.x = "gaul2_code", by.y = "code", all.x = TRUE)
  x[, p_source := fifelse(is.na(tier_t), NA_character_, "adm2")]
  y <- merge(x[is.na(p_source), .(gaul2_code, gaul1_code)], a1, by.x = "gaul1_code", by.y = "code", all.x = TRUE)
  y[, p_source := fifelse(is.na(tier_t), NA_character_, "adm1")]
  out <- rbind(x[!is.na(p_source)], y, fill = TRUE)
  out[, hazard := h]
  out
}
up <- rbindlist(lapply(HAZ_USED, pick))
stopifnot(!anyDuplicated(up[, .(gaul2_code, hazard)]))

no_p <- up[is.na(p_source), .(gaul2_code, hazard)]
no_p <- merge(no_p, units[, .(gaul2_code, iso3, admin1_name, admin2_name, maize_prod_t)], by = "gaul2_code")
no_p[, reason := "no admin-2 or admin-1 row with complete tiers and tier_t >= MIN_TIER_T"]
af_write(no_p, "units_no_p.csv", writer = fwrite)
up <- up[!is.na(p_source)]
message(sprintf("P source (unit x hazard): %s; unit-hazard pairs without P: %d (%d units)",
                paste(names(table(up$p_source)), table(up$p_source), collapse = ", "),
                nrow(no_p), uniqueN(no_p$gaul2_code)))

# ---- E[L] per unit x hazard, B_norm ---------------------------------------
Ltab <- unique(grid[, .(hazard, severity, L)])
pl <- melt(up[, .(gaul2_code, hazard, p_moderate, p_severe, p_extreme)],
           id.vars = c("gaul2_code", "hazard"), variable.name = "severity", value.name = "p_band")
pl[, severity := sub("^p_", "", severity)]
el <- merge(pl, Ltab, by = c("hazard", "severity"))[, .(EL = sum(p_band * L)), by = .(gaul2_code, hazard)]
up <- merge(up, el, by = c("gaul2_code", "hazard"))
up <- merge(up, units[, .(gaul2_code, B_obs = B, aez16simple)], by = "gaul2_code")
up[, B_norm := B_obs / (1 - EL)]

# ---- Parameter sets: grid + ERA all-SSA + ERA stratum ---------------------
era_strat <- era_g[stratum != "all", .(practice, aez16simple = stratum, g_strat = g, g_stratum_n = n_studies)]
era_all   <- era_g[stratum == "all", .(practice, g_all = g)]
sets <- unique(grid[, .(param_id, practice, hazard, g_source, g, m_profile)])
# era_aez sets reuse the era rows' (practice, hazard, m_profile) with a unit-specific g.
aez_sets <- sets[g_source == "era"][, `:=`(g_source = "era_aez", param_id = param_id + max(sets$param_id))]
sets <- rbind(sets, aez_sets)
split_file <- file.path(af_out_dir(), "era_split", "era_split_g.csv")
if (file.exists(split_file)) {
  gn <- fread(split_file)[hazard == "dry" & threshold == "moderate", .(practice, g_norm)]
  norm_sets <- sets[g_source == "era"][, `:=`(g_source = "era_norm", param_id = param_id + 2 * max(grid$param_id))]
  norm_sets <- merge(norm_sets, gn, by = "practice")[, g := g_norm][, g_norm := NULL]
  sets <- rbind(sets, norm_sets)
  message("era_norm sets added from ", split_file)
}

ug <- merge(up[, .(gaul2_code, hazard, aez16simple, B_norm)], sets, by = "hazard", allow.cartesian = TRUE)
ug <- merge(ug, era_strat, by = c("practice", "aez16simple"), all.x = TRUE)
ug <- merge(ug, era_all, by = "practice", all.x = TRUE)
ug[g_source == "era_aez", g := fifelse(!is.na(g_strat) & g_stratum_n >= MIN_STRATUM_STUDIES, g_strat, g_all)]
ug[g_source != "era_aez", g_stratum_n := NA_integer_]
ug <- ug[!is.na(g), .(gaul2_code, hazard, param_id, practice, g_source, g, g_stratum_n, m_profile, B = B_norm)]

# ---- Expand to severity and compute ---------------------------------------
Lm <- unique(grid[, .(hazard, severity, m_profile, L, m)])
long <- merge(ug, Lm, by = c("hazard", "m_profile"), allow.cartesian = TRUE)
long <- merge(long, pl, by = c("gaul2_code", "hazard", "severity"))
stopifnot(nrow(long) == 3 * nrow(ug))
setorder(long, gaul2_code, param_id, severity)

ea <- as.data.table(af_expected_annual(as.data.frame(long), c("gaul2_code", "param_id")))
d12 <- max(abs(ea$F1_total - ea$F3_total)); d23 <- max(abs(ea$F2_total - ea$F3_total))
if (max(d12, d23) > 1e-9) stop("F1/F2/F3 totals differ: ", signif(max(d12, d23), 3))

# ---- Assemble ------------------------------------------------------------
res <- merge(ea, unique(ug[, .(gaul2_code, param_id, hazard, practice, g_source, g, g_stratum_n, m_profile)]),
             by = c("gaul2_code", "param_id"))
res <- merge(res, up[, .(gaul2_code, hazard, p_source, mono_adj, p_moderate, p_severe, p_extreme,
                         p_none, p_dry_wet, tier_t, EL, B_obs, B_norm)], by = c("gaul2_code", "hazard"))
res <- merge(units[, .(gaul2_code, gaul1_code, iso3, admin1_name, admin2_name, maize_prod_t,
                       maize_harv_ha, fs1, fs2, aez, aez16simple, bimodal, bimodal_share)],
             res, by = "gaul2_code")
res[, `:=`(F1_minus_F3_al = F1_al - F3_al, F1_negative = F1_al < 0, p_stress_ge_099 = p_stress >= P_ONE)]
for (f in AF_FORMULATIONS) {  # robust share, in [-1, 1]; _ashare (al / total) kept for continuity
  res[, (paste0(f, "_share")) := get(paste0(f, "_al")) / (abs(get(paste0(f, "_prod"))) + abs(get(paste0(f, "_al"))))]
}
setcolorder(res, c("iso3", "admin1_name", "admin2_name", "gaul1_code", "gaul2_code",
                   "practice", "hazard", "g_source", "g", "g_stratum_n", "m_profile", "param_id",
                   "p_source", "B_obs", "B_norm", "EL"))

path <- af_write(res, "af_results_maize_historic.parquet")
message(sprintf("%s: %d rows (%d units x %d parameter sets); p_stress >= %.2f in %d unit-hazard pairs",
                path, nrow(res), uniqueN(res$gaul2_code), uniqueN(res$param_id), P_ONE,
                uniqueN(res[p_stress_ge_099 == TRUE, .(gaul2_code, hazard)])))
