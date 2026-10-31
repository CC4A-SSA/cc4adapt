# ------------------------------------------------------------------------------
# code/explore/af_formulations/13_kenya_dotplot.R
#
# EXPLORATORY. Static PNG of the Kenyan per-site comparison (11_gcm_summary.R
# output): observed vs two climate models vs the Atlas tier.
# Run from the repo root:  Rscript code/explore/af_formulations/13_kenya_dotplot.R
# Writes AF_OUT_DIR/era_split/kenya_gcm_vs_obs.png
# ------------------------------------------------------------------------------
suppressMessages({library(data.table); library(ggplot2)})
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
source("code/explore/af_formulations/R/paths.R")
g <- fread(file.path(af_out_dir(), "era_split", "gcm_vs_obs_summary.csv"))
g[, label := paste0(admin2_name, "  (", sub(" B.*", "", Site.Key), ")")]
g <- g[order(atlas_p_dry, p_dry_obs_all)]; g[, label := factor(label, levels = unique(label))]
long <- melt(g[, .(label, `Observed (CHIRPS v3 + AgERA5)` = p_dry_obs_all, `NEX-GDDP ACCESS-CM2` = p_dry_gcm_accesscm2,
                   `NEX-GDDP GFDL-ESM4` = p_dry_gcm_gfdlesm4, `Atlas tier (admin-2)` = atlas_p_dry)], id.vars = "label")
seg <- g[, .(label, lo = pmin(p_dry_obs_all, p_dry_gcm_accesscm2, p_dry_gcm_gfdlesm4, atlas_p_dry), hi = pmax(p_dry_obs_all, p_dry_gcm_accesscm2, p_dry_gcm_gfdlesm4, atlas_p_dry))]
cols <- c(`Observed (CHIRPS v3 + AgERA5)` = "#eb6834", `NEX-GDDP ACCESS-CM2` = "#1baf7a", `NEX-GDDP GFDL-ESM4` = "#0d7a52", `Atlas tier (admin-2)` = "#0b0b0b")
shp  <- c(`Observed (CHIRPS v3 + AgERA5)` = 16, `NEX-GDDP ACCESS-CM2` = 16, `NEX-GDDP GFDL-ESM4` = 1, `Atlas tier (admin-2)` = 15)
p <- ggplot() + geom_segment(data = seg, aes(x = lo, xend = hi, y = label, yend = label), colour = "#c3c2b7", linewidth = 1) +
  geom_point(data = long, aes(value, label, colour = variable, shape = variable), size = 3.2, stroke = 1.3) +
  scale_colour_manual(values = cols, name = NULL) + scale_shape_manual(values = shp, name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Share of seasons with NDWS at or above moderate, 1995-2013, GGCMI maize window", y = NULL,
       title = "Kenyan ERA maize sites: dry-stress frequency by climate source vs the Atlas tier",
       subtitle = "Same sites, window, seasons and water balance; Atlas value is the containing admin-2 unit. Exploratory.") +
  theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title.position = "plot")
ggsave(file.path(af_out_dir(), "era_split", "kenya_gcm_vs_obs.png"), p, width = 10, height = 6.5, dpi = 150, bg = "white")
cat("ok\n")
