# ------------------------------------------------------------------------------
# code/explore/af_formulations/08_era_split.R
#
# EXPLORATORY. Split the ERA maize effect of each practice into normal-season
# and stress-season estimates using the site-season hazard classes from
# 07_era_seasons.py, and test whether the stress-season effect scales with
# control yield (proportional) or not (additive).
#
# Run from the repo root:  Rscript code/explore/af_formulations/08_era_split.R
#
# Inputs:  data/staging/af_formulations/era_comparisons.csv (03b)
#          data/staging/af_formulations/era_site_season_hazards.csv (07)
# Outputs (AF_OUT_DIR/era_split/):
#   era_split_g.csv        g_norm, g_stress per practice x hazard x severity
#                          threshold, study-averaged, with site-season counts
#   era_split_test.csv     proportional-vs-additive regression per practice x
#                          hazard
#   era_split_sites.csv    per comparison: hazard bands and derived quantities
#
# Definitions. A site-season is "stress" for hazard h when its Atlas-style
# band is at or above the chosen threshold (default: moderate). g = exp(mean
# study-level log(MeanT/MeanC)) - 1 within each class. The control's
# stress-season loss L_c is not observable from ERA (no same-site normal-year
# control), so the proportional-vs-additive test uses cross-site variation:
#   proportional: log(MeanT/MeanC) in stress seasons independent of MeanC
#   additive:     (MeanT - MeanC), Mg/ha, in stress seasons independent of MeanC
# Both regressions are reported with study-clustered standard errors; which
# slope is closer to zero is the evidence. Low power is expected.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")
source("code/explore/af_formulations/R/paths.R")

IN  <- "data/staging/af_formulations"
OUT <- file.path(af_out_dir(), paste0("era_split", Sys.getenv("AF_HAZ_TAG", ""))); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
THRESH <- c("moderate", "severe")
HAZ <- c(dry = "NDWS", wet = "NDWL0", heat = "NTx35")
ORD <- c(none = 0, moderate = 1, severe = 2, extreme = 3)

cmp <- fread(file.path(IN, "era_comparisons.csv"))
stopifnot(all(cmp$Units == "Mg/ha"))          # 03b harmonises units; the additive test needs one unit
n_irr <- sum(cmp$irrigated %in% TRUE)
cmp <- cmp[!(irrigated %in% TRUE)]           # water stress is not observable for irrigated trials
HAZ_FILE <- Sys.getenv("AF_HAZ_FILE", file.path(IN, "era_site_season_hazards.csv"))
haz <- fread(HAZ_FILE)
cmp[, slot := fifelse(is.na(season_idx), 1L, as.integer(season_idx))]
haz[, slot := as.integer(season_slot)]

d <- merge(cmp, haz[, c("Site.Key", "slot", "year", "calendar_status", "counted_days", "total_mm", "wrsi",
                        paste0(HAZ, "_monthly_mean"), paste0(HAZ, "_band")), with = FALSE],
           by.x = c("Site.Key", "slot", "yr"), by.y = c("Site.Key", "slot", "year"))
d[, lrr := log(MeanT / MeanC)]
d[, diff := MeanT - MeanC]
message(sprintf("matched %d of %d rainfed comparisons (%d irrigated excluded; %d site-seasons)", nrow(d), nrow(cmp), n_irr, uniqueN(d[, .(Site.Key, slot, yr)])))

# ---- g_norm and g_stress ---------------------------------------------------
gsplit <- rbindlist(lapply(names(HAZ), function(h) rbindlist(lapply(THRESH, function(th) {
  b <- d[[paste0(HAZ[[h]], "_band")]]
  x <- copy(d)[, stress := ORD[b] >= ORD[[th]]]
  x <- x[!is.na(stress)]
  st <- x[, .(lrr = mean(lrr), n = .N), by = .(practice, Code, stress)]
  ss <- x[, .(n_site_seasons = uniqueN(paste(Site.Key, slot, yr))), by = .(practice, stress)]
  g  <- st[, .(g = exp(mean(lrr)) - 1, g_median = exp(median(lrr)) - 1, n_studies = .N, n_obs = sum(n)),
           by = .(practice, stress)]
  merge(g, ss, by = c("practice", "stress"))[, `:=`(hazard = h, threshold = th)]
}))))
gw <- dcast(gsplit, practice + hazard + threshold ~ stress, value.var = c("g", "g_median", "n_studies", "n_obs", "n_site_seasons"))
setnames(gw, c("g_FALSE", "g_TRUE", "g_median_FALSE", "g_median_TRUE", "n_studies_FALSE", "n_studies_TRUE",
               "n_obs_FALSE", "n_obs_TRUE", "n_site_seasons_FALSE", "n_site_seasons_TRUE"),
         c("g_norm", "g_stress", "g_norm_median", "g_stress_median", "studies_norm", "studies_stress",
           "obs_norm", "obs_stress", "site_seasons_norm", "site_seasons_stress"), skip_absent = TRUE)
gw[, g_stress_minus_norm := g_stress - g_norm]
fwrite(gw, file.path(OUT, "era_split_g.csv"))
print(gw[, .(practice, hazard, threshold, g_norm = round(g_norm, 3), g_stress = round(g_stress, 3),
             studies_norm, studies_stress, site_seasons_norm, site_seasons_stress)])

# ---- proportional vs additive -----------------------------------------------
cluster_se <- function(fit, cl) {
  X <- model.matrix(fit); u <- residuals(fit); cl <- as.factor(cl)
  meat <- Reduce(`+`, lapply(split(seq_along(u), cl), function(i) { s <- crossprod(X[i, , drop = FALSE], u[i]); s %*% t(s) }))
  bread <- solve(crossprod(X)); G <- nlevels(cl); n <- nrow(X); k <- ncol(X)
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
  sqrt(diag(V))
}
tests <- rbindlist(lapply(names(HAZ), function(h) {
  b <- d[[paste0(HAZ[[h]], "_band")]]
  x <- copy(d)[, stress := ORD[b] >= ORD[["moderate"]]][stress %in% TRUE]
  rbindlist(lapply(unique(x$practice), function(p) {
    y <- x[practice == p]
    if (nrow(y) < 10 || uniqueN(y$Code) < 3) return(NULL)
    f1 <- lm(lrr ~ log(MeanC), data = y); f2 <- lm(diff ~ MeanC, data = y)
    s1 <- cluster_se(f1, y$Code); s2 <- cluster_se(f2, y$Code)
    data.table(practice = p, hazard = h, n_obs = nrow(y), n_studies = uniqueN(y$Code),
               prop_slope = coef(f1)[[2]], prop_se = s1[[2]], prop_t = coef(f1)[[2]] / s1[[2]],
               add_slope = coef(f2)[[2]], add_se = s2[[2]], add_t = coef(f2)[[2]] / s2[[2]])
  }))
}))
fwrite(tests, file.path(OUT, "era_split_test.csv"))
print(tests[, .(practice, hazard, n_obs, n_studies, prop_t = round(prop_t, 2), add_t = round(add_t, 2))])
fwrite(d, file.path(OUT, "era_split_sites.csv"))
message("wrote ", OUT)
