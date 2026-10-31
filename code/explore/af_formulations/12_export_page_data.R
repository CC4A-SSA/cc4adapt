# ------------------------------------------------------------------------------
# code/explore/af_formulations/12_export_page_data.R
#
# EXPLORATORY. Builds docs/af_formulations/data.js, the data file behind the
# interactive page, from the outputs of 05, 06, 08, 09 and 11.
#
# Run from the repo root:  Rscript code/explore/af_formulations/12_export_page_data.R
# Reads AF_OUT_DIR (default data/outputs/af_formulations) and
# data/staging/af_formulations/. Writes window.AF_DATA = {...} as JS.
# ------------------------------------------------------------------------------
suppressMessages({library(data.table); library(jsonlite)})
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
source("code/explore/af_formulations/R/paths.R")
r <- as.data.table(arrow::read_parquet(file.path(af_out_dir(), "af_results_maize_historic.parquet")))
D <- file.path(af_out_dir(), "diagnostics")
d3 <- fread(file.path(D, "d3_label_flips.csv")); d4 <- fread(file.path(D, "d4_strata.csv"))
units <- unique(r[, .(gaul2_code, iso3, admin1_name, admin2_name, hazard, fs1, aez, aez16simple, bimodal, p_source, p_stress_ge_099,
                      maize_kt = round(maize_prod_t/1e3, 1), B_obs = round(B_obs, 2), B_norm = round(B_norm, 2),
                      p_stress = round(p_stress, 3), EL = round(EL, 4), p_mod = round(p_moderate,3), p_sev = round(p_severe,3), p_ext = round(p_extreme,3))])
# per-unit results for the reference sets (era g; const_0.50 and decl_mid), all practices/hazards
ref <- r[g_source == "era" & m_profile %in% c("const_0.50","decl_mid"),
         .(gaul2_code, hazard, practice, g_source, m_profile, g = round(g,3),
           F0_al = round(F0_al,4), F1_al = round(F1_al,4), F2_al = round(F2_al,4), F3_al = round(F3_al,4), F4_al = round(F4_al,4), F5_al = round(F5_al,4),
           F0_prod = round(F0_prod,4), F1_prod = round(F1_prod,4), F2_prod = round(F2_prod,4), F3_prod = round(F3_prod,4), F4_prod = round(F4_prod,4), F5_prod = round(F5_prod,4),
           F0_share = round(F0_share,3), F1_share = round(F1_share,3), F2_share = round(F2_share,3), F3_share = round(F3_share,3), F4_share = round(F4_share,3), F5_share = round(F5_share,3))]
# severity-mix: ext_share & m_eff for decl profiles (cover_crops dry, era)
pl <- melt(unique(r[, .(gaul2_code, hazard, p_moderate, p_severe, p_extreme)]), id.vars=c("gaul2_code","hazard"), variable.name="severity", value.name="p")
pl[, severity := sub("^p_","",severity)]
grid <- fread("data/staging/af_formulations/params/param_grid.csv")
lm_ <- unique(grid[, .(hazard, m_profile, severity, L, m)])
sm <- merge(pl, lm_[m_profile == "decl_mid"], by=c("hazard","severity"), allow.cartesian=TRUE)[
  , .(EL = sum(p*L), m_eff = if (sum(p*L)>0) sum(p*L*m)/sum(p*L) else NA_real_, ext_share = if (sum(p*L)>0) sum((p*L)[severity=="extreme"])/sum(p*L) else NA_real_), by=.(gaul2_code, hazard, m_profile)]
sm <- sm[!is.na(m_eff), .(gaul2_code, hazard, m_profile, m_eff=round(m_eff,3), ext_share=round(ext_share,3))]
era_g <- fread("data/staging/af_formulations/params/era_g.csv")
m_prof <- fread("data/staging/af_formulations/params/param_grid.csv")[, .(m_moderate = m[severity=="moderate"][1], m_severe = m[severity=="severe"][1], m_extreme = m[severity=="extreme"][1]), by=m_profile]
L <- unique(grid[, .(hazard, severity, L)])
strata <- d4[g_source %in% c("era","era_aez","era_norm") & level %in% c("fs1","aez","country"),
             .(level, practice, hazard, g_source, m_profile, iso3, fs1, aez, units, units_adm1_fallback, units_p_stress_ge_099, maize_kt=round(maize_kt),
               B_obs_med=round(B_obs_med,2), B_norm_med=round(B_norm_med,2), p_stress_med=round(p_stress_med,3), EL_med=round(EL_med,4), gap_med=round(gap_med,4),
               F0_al_med=round(F0_al_med,4), F1_al_med=round(F1_al_med,4), F2_al_med=round(F2_al_med,4), F3_al_med=round(F3_al_med,4), F4_al_med=round(F4_al_med,4), F5_al_med=round(F5_al_med,4),
               F0_share_med=round(F0_share_med,3), F1_share_med=round(F1_share_med,3), F2_share_med=round(F2_share_med,3), F3_share_med=round(F3_share_med,3), F4_share_med=round(F4_share_med,3), F5_share_med=round(F5_share_med,3))]
flips <- d3[g_source=="era", .(practice, hazard, g=round(g,3), m_profile, units, F3_adaptation=round(F3_adaptation,3), flip_F1_vs_F3=round(flip_F1_vs_F3,3), flip_F2_vs_F3=round(flip_F2_vs_F3,3), flip_F5_vs_F3=round(flip_F5_vs_F3,3))]
E <- file.path(af_out_dir(), "era_split")
obs <- fread(file.path(E, "obs_vs_atlas_sites.csv"))[, .(sites=.N, seasons=sum(n_seasons), obs_dry=round(weighted.mean(p_dry_obs, n_seasons),2), atlas_dry=round(mean(p_stress_dry),2),
  obs_dry_sev=round(weighted.mean(p_dry_sev_obs, n_seasons),2), atlas_dry_sev=round(mean(p_sev_dry),2), obs_wet=round(weighted.mean(p_wet_obs, n_seasons),2), atlas_wet=round(mean(p_stress_wet),2)), by=iso3]
a <- fread("data/staging/af_formulations/era_site_season_hazards.csv")[season_slot == 1]
b <- fread("data/staging/af_formulations/era_site_season_hazards_ggcmi.csv")
syi <- unique(fread("data/staging/af_formulations/era_site_years.csv")[, .(Site.Key, iso)])
pw <- merge(merge(a[, .(Site.Key, year, days_site = counted_days, dry_site = NDWS_band != "none")],
                  b[, .(Site.Key, year, days_ggcmi = counted_days, dry_ggcmi = NDWS_band != "none")], by = c("Site.Key", "year")), syi, by = "Site.Key")
win <- pw[, .(sites = uniqueN(Site.Key), seasons = .N, days_site = round(median(days_site)), days_ggcmi = round(median(days_ggcmi)),
              p_dry_site = round(mean(dry_site), 2), p_dry_ggcmi = round(mean(dry_ggcmi), 2)), by = iso][order(-seasons)][seasons >= 8]
obs_g <- fread(file.path(E, "obs_vs_atlas_sites_ggcmi.csv"))[, .(obs_dry_ggcmi = round(weighted.mean(p_dry_obs, n_seasons), 2)), by = iso3]
obs <- merge(obs, obs_g, by = "iso3", all.x = TRUE)
gcm_file <- file.path(E, "gcm_vs_obs_summary.csv")
gcm <- if (file.exists(gcm_file)) {
  gs <- fread(gcm_file); srcs <- sub("^p_dry_", "", grep("^p_dry_", names(gs), value = TRUE))
  rbindlist(lapply(srcs, function(sv) data.table(src = sv, sites = sum(!is.na(gs[[paste0("p_dry_", sv)]])),
    seasons = sum(gs[[paste0("n_", sv)]], na.rm = TRUE), p_dry = round(mean(gs[[paste0("p_dry_", sv)]], na.rm = TRUE), 2),
    p_sev = round(mean(gs[[paste0("p_sev_", sv)]], na.rm = TRUE), 2), ndws = round(mean(gs[[paste0("ndws_", sv)]], na.rm = TRUE), 1),
    r_atlas = round(cor(gs[[paste0("p_dry_", sv)]], gs$atlas_p_dry, use = "complete.obs"), 2))))[
    , rbind(.SD, data.table(src = "atlas", sites = nrow(gs), seasons = NA, p_dry = round(mean(gs$atlas_p_dry), 2), p_sev = round(mean(gs$atlas_p_sev), 2), ndws = NA, r_atlas = NA))]
} else NULL
gcm_sites <- if (file.exists(gcm_file)) {
  gs <- fread(gcm_file)
  gs[, .(site = Site.Key, admin2 = admin2_name, atlas = round(atlas_p_dry, 2), atlas_sev = round(atlas_p_sev, 2), p_source = atlas_p_source,
         obs_all = round(p_dry_obs_all, 2), obs_era = round(p_dry_obs_era, 2), n_era = n_obs_era,
         access = round(p_dry_gcm_accesscm2, 2), gfdl = round(p_dry_gcm_gfdlesm4, 2))][order(-atlas)]
} else NULL
era_split <- list(g = fread(file.path(E, "era_split_g.csv")), test = fread(file.path(E, "era_split_test.csv")), obs = obs, window = win, gcm = gcm, gcm_sites = gcm_sites)
out <- list(generated = format(Sys.time(), "%Y-%m-%d"), era_split = era_split, n_units = uniqueN(r$gaul2_code), n_param_sets = uniqueN(r$param_id),
            era_g = era_g, m_profiles = m_prof, loss = L, units = units, ref = ref, sevmix = sm, strata = strata, flips = flips)
writeLines(paste0("window.AF_DATA = ", toJSON(out, dataframe="rows", na="null", auto_unbox=TRUE, digits=NA), ";"), "docs/af_formulations/data.js")
cat("rows: units", nrow(units), "ref", nrow(ref), "sevmix", nrow(sm), "strata", nrow(strata), "flips", nrow(flips), "\n")
