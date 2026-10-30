# ------------------------------------------------------------------------------
# code/explore/af_formulations/09_obs_vs_atlas.R
#
# EXPLORATORY. Observed stress frequency at ERA sites (07: CHIRPS v3 + AgERA5
# site-seasons on site windows) against the Atlas P(h, >= moderate) of the
# admin-2 unit containing the site (05). Speaks to the benchmark question
# (science review finding 3): the two are not the same quantity.
#
# Run from the repo root:  Rscript code/explore/af_formulations/09_obs_vs_atlas.R
# Output: AF_OUT_DIR/era_split/obs_vs_atlas_sites.csv
# ------------------------------------------------------------------------------
suppressMessages({library(data.table); library(sf)})
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
source("code/explore/af_formulations/R/paths.R")
# AF_HAZ_FILE selects the classification (default: site windows; the
# --window ggcmi output gives the Atlas-window variant). AF_HAZ_TAG names the output.
HAZ <- Sys.getenv("AF_HAZ_FILE", "data/staging/af_formulations/era_site_season_hazards.csv")
TAG <- Sys.getenv("AF_HAZ_TAG", "")
h <- fread(HAZ)
obs <- h[, .(n_seasons = .N, p_dry_obs = mean(NDWS_band != "none"), p_dry_sev_obs = mean(NDWS_band %in% c("severe","extreme")),
             p_wet_obs = mean(NDWL0_band != "none"), lat = lat[1], lon = lon[1]), by = Site.Key]
g <- st_as_sf(arrow::read_parquet("data/staging/af_formulations/gaul24_adm2_ken_eth_nga.parquet")); if (is.na(st_crs(g))) st_crs(g) <- 4326
pts <- st_as_sf(obs, coords = c("lon","lat"), crs = 4326)
j <- st_join(pts, g[, c("iso3","admin2_name","gaul2_code")], left = FALSE)
r <- as.data.table(arrow::read_parquet("data/outputs/af_formulations/af_results_maize_historic.parquet"))
atl <- unique(r[hazard %in% c("dry","wet"), .(gaul2_code, hazard, p_stress, p_sev = p_severe + p_extreme, p_source)])
jd <- merge(as.data.table(st_drop_geometry(j)), dcast(atl, gaul2_code + p_source ~ hazard, value.var = c("p_stress","p_sev")), by = "gaul2_code")
cat("ERA sites inside KEN/ETH/NGA maize units:", nrow(jd), "| seasons:", sum(jd$n_seasons), "\n")
print(jd[, .(sites = .N, seasons = sum(n_seasons), obs_dry = round(weighted.mean(p_dry_obs, n_seasons),2), atlas_dry = round(mean(p_stress_dry),2),
             obs_dry_sev = round(weighted.mean(p_dry_sev_obs, n_seasons),2), atlas_dry_sev = round(mean(p_sev_dry),2),
             obs_wet = round(weighted.mean(p_wet_obs, n_seasons),2), atlas_wet = round(mean(p_stress_wet),2)), by = iso3])
cat("site-level correlation obs vs Atlas P(dry>=mod):", round(cor(jd$p_dry_obs, jd$p_stress_dry), 2), "\n")
fwrite(jd, file.path(af_out_dir(), "era_split", paste0("obs_vs_atlas_sites", TAG, ".csv")))
