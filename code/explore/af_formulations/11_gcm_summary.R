# ------------------------------------------------------------------------------
# code/explore/af_formulations/11_gcm_summary.R
#
# EXPLORATORY. Where does the Atlas-vs-observed stress-frequency gap come from?
# Per ERA site (Kenya by default), P(dry >= moderate) on the same GGCMI window:
#   obs_all      observed climate, every season 1995-2014 (10, source paired)
#   obs_era      observed climate, ERA trial seasons only (07 --window ggcmi)
#   gcm_<model>  NEX-GDDP historical, 1995-2014 (10, source nex_gddp)
#   atlas        Atlas tier P for the containing admin-2 unit (05)
# Differences: obs_era - obs_all = trial-year selection; gcm - obs_all = GCM
# climate through the same water balance; atlas - gcm = the Atlas's own
# water balance / thresholds / aggregation.
#
# Run from the repo root:  Rscript code/explore/af_formulations/11_gcm_summary.R
# Output: AF_OUT_DIR/era_split/gcm_vs_obs_summary.csv (+ console summary)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
source("code/explore/af_formulations/R/paths.R")
IN <- "data/staging/af_formulations"

files <- list.files(IN, pattern = "^gcm_vs_obs_sites.*\\.csv$", full.names = TRUE)
files <- files[!grepl("dryrun", files)]
g <- rbindlist(lapply(files, fread), fill = TRUE)[year <= 2013]   # NEX historical ends 2014; last complete window is 2013
g[, src := fifelse(source == "paired", "obs_all", paste0("gcm_", tolower(gsub("-", "", model))))]
g[, dry := NDWS_band != "none"]
g[, sev := NDWS_band %in% c("severe", "extreme")]
per <- g[, .(n = .N, p_dry = mean(dry), p_sev = mean(sev), ndws = mean(NDWS_monthly_mean, na.rm = TRUE)), by = .(Site.Key, src)]

era <- fread(file.path(IN, "era_site_season_hazards_ggcmi.csv"))[Site.Key %in% unique(g$Site.Key)]
era_s <- era[, .(n = .N, p_dry = mean(NDWS_band != "none"), p_sev = mean(NDWS_band %in% c("severe", "extreme")),
                 ndws = mean(NDWS_monthly_mean, na.rm = TRUE)), by = Site.Key][, src := "obs_era"]
per <- rbind(per, era_s, fill = TRUE)

atlas <- fread(file.path(af_out_dir(), "era_split", "obs_vs_atlas_sites.csv"))[, .(Site.Key, iso3, admin2_name, gaul2_code,
            atlas_p_dry = p_stress_dry, atlas_p_sev = p_sev_dry, atlas_p_source = p_source)]
w <- dcast(per, Site.Key ~ src, value.var = c("p_dry", "p_sev", "ndws", "n"))
w <- merge(atlas, w, by = "Site.Key")
fwrite(w, file.path(af_out_dir(), "era_split", "gcm_vs_obs_summary.csv"))

srcs <- setdiff(unique(per$src), character(0))
cat(sprintf("%d sites; sources: %s\n", nrow(w), paste(srcs, collapse = ", ")))
summ <- per[, .(sites = uniqueN(Site.Key), seasons = sum(n), p_dry = round(mean(p_dry), 2), p_sev = round(mean(p_sev), 2),
                ndws_per_month = round(mean(ndws), 1)), by = src]
summ <- rbind(summ, data.table(src = "atlas", sites = nrow(w), seasons = NA, p_dry = round(mean(w$atlas_p_dry), 2),
                               p_sev = round(mean(w$atlas_p_sev), 2), ndws_per_month = NA))
print(summ)
cat("\nsite-level correlations with Atlas P(dry):\n")
for (s in srcs) { v <- w[[paste0("p_dry_", s)]]; if (!is.null(v)) cat(sprintf("  %-14s r = %.2f (n = %d)\n", s, cor(v, w$atlas_p_dry, use = "complete.obs"), sum(!is.na(v)))) }
