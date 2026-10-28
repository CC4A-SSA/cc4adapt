# ------------------------------------------------------------------------------
# code/explore/af_formulations/03_params.R
#
# EXPLORATORY. Step 5.3: normal-year gain g per practice from ERA, and the
# expanded parameter grid (practice x hazard x severity x g x m profile).
#
# Run from the repo root:  Rscript code/explore/af_formulations/03_params.R
#
# Inputs:
#   * code/explore/af_formulations/config/params.yaml
#   * data/inputs/era/ERA_Derived.rda (sibling solutions checkout, Apr 2026)
# Outputs (data/staging/af_formulations/params/):
#   * era_g.csv       g per practice, overall and by ERA AEZ16simple class
#   * param_grid.csv  one row per practice x hazard x severity x g x m profile
#
# ERA notes:
#   [A] Single-practice comparisons only (PrName exactly as configured), maize
#       Crop Yield, Invalid.Comparison != "Y", MeanT and MeanC > 0.
#   [B] Effect = log(MeanT / MeanC), computed here. In this ERA_Derived vintage
#       `yi` is the plain ratio MeanT / MeanC and `pc` equals 100 * yi - 200,
#       i.e. 100 points below the percent change, so neither is used.
#   [C] Studies (Code) are averaged first so heavily replicated studies do not
#       dominate; g = exp(median of study means) - 1. Mean and IQR reported.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}

cfg <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")
OUT <- "data/staging/af_formulations/params"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- ERA g --------------------------------------------------------------
env <- new.env()
nm  <- load("data/inputs/era/ERA_Derived.rda", envir = env)
era <- as.data.table(get(nm[1], envir = env))

pr  <- vapply(cfg$practices, `[[`, "", "era_prname")
sel <- era[Product.Simple == "Maize" & Out.SubInd == "Crop Yield" &
             PrName %in% pr & (is.na(Invalid.Comparison) | Invalid.Comparison != "Y")]
sel[, `:=`(MeanT = as.numeric(MeanT), MeanC = as.numeric(MeanC))]
sel <- sel[is.finite(MeanT) & is.finite(MeanC) & MeanT > 0 & MeanC > 0]
sel[, lrr := log(MeanT / MeanC)]
sel[, practice := names(pr)[match(PrName, pr)]]

summarise_g <- function(d) {
  st <- d[, .(lrr = mean(lrr)), by = Code]
  q  <- quantile(st$lrr, c(0.25, 0.75))
  data.table(n_obs = nrow(d), n_studies = nrow(st),
             g = exp(median(st$lrr)) - 1, g_mean = exp(mean(st$lrr)) - 1,
             g_q25 = exp(q[[1]]) - 1, g_q75 = exp(q[[2]]) - 1)
}
era_g <- rbind(
  sel[, summarise_g(.SD), by = practice][, stratum := "all"],
  sel[!is.na(AEZ16simple), summarise_g(.SD), by = .(practice, stratum = AEZ16simple)],
  fill = TRUE
)
setcolorder(era_g, c("practice", "stratum"))
fwrite(era_g, file.path(OUT, "era_g.csv"))
print(era_g, digits = 3)

# ---- Parameter grid ------------------------------------------------------
g_era <- era_g[stratum == "all", .(practice, g, g_source = "era")]
rows <- list()
for (p in names(cfg$practices)) {
  g_vals <- rbind(data.table(g = unlist(cfg$g_grid), g_source = "grid"),
                  g_era[practice == p, .(g, g_source)])
  for (h in cfg$practices[[p]]$hazards) {
    lh <- cfg$loss_frac[[cfg$hazard_map[[h]]]]
    for (mp in names(cfg$m_profiles)) {
      for (s in names(lh)) {
        rows[[length(rows) + 1]] <- data.table(
          practice = p, hazard = h, severity = s, L = lh[[s]],
          m_profile = mp, m = cfg$m_profiles[[mp]][[s]], g_vals
        )
      }
    }
  }
}
grid <- rbindlist(rows)
grid[, param_id := .GRP, by = .(practice, hazard, g_source, g, m_profile)]
fwrite(grid, file.path(OUT, "param_grid.csv"))
message(sprintf("param_grid.csv: %d rows, %d parameter sets", nrow(grid), uniqueN(grid$param_id)))
