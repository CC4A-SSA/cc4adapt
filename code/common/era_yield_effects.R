# ------------------------------------------------------------------------------
# code/common/era_yield_effects.R
#
# ERA yield effect of a CCA practice on a MapSPAM crop, shared by
# Step 2 (production benefit) and Step 3 (adaptation benefit).
#
# Until 2026-11 this logic lived only in step2_production.R. Step 3 now needs
# the same normal-year relative gain g to build the practice's own yield
# A = B * (1 + g) (see code/common/adaptation_split.R), so the lookup tables
# and the effect function were moved here unchanged. Step 2's output is
# numerically identical to before.
#
# Method (unchanged from Step 2 design note [B]):
#   * keep ERA rows whose Out.SubInd mentions yield or productivity and whose
#     MeanT and MeanC are both positive;
#   * RR = log(MeanT / MeanC) per row;
#   * for a (practice, crop) pair, take the rows whose PrName contains any of
#     the practice's ERA labels and whose Product matches the SPAM code or the
#     FAO item name; g = exp(mean(RR)) - 1.
#   The pooled RR mixes all seasons (stress and normal years alike). Splitting
#   it by season is discussed in code/explore/af_formulations/ and is NOT done
#   here.
#
# Pure base R; sourceable without config.R (tests rely on that).
# ------------------------------------------------------------------------------

# CCA practice key -> ERA PrName labels (as they appear in ERA_Derived.rda).
# Compound ERA labels (e.g. "Alleycropping-Intercropping") are handled via
# substring match.
CCA_TO_ERA <- list(
  agroforestry_cropland       = c("Agroforestry", "Alleycropping", "Hedgerow", "Parklands", "Silvopasture"),
  cover_crops                 = c("Green Manure", "Cover Crop", "Crop Residue Incorporation"),
  biochar                     = c("Biochar"),
  awd_rice                    = c("Improved Rice Management", "Alternate Wetting and Drying"),
  enhanced_rock_weathering    = c("Silicate", "Enhanced Weathering", "Rock Amendment"),
  grassland_soil_carbon       = c("Grazing Management", "Improved Grazing"),
  grassland_reduced_emissions = c("Grazing Management"),
  grassland_restoration       = c("Grassland Restoration", "Pasture Restoration"),
  savanna_agroforestry        = c("Agroforestry", "Parklands")
)

# MapSPAM 4-letter code -> FAO Item string (as FAOSTAT spells them).
SPAM_TO_FAO <- c(
  maiz = "Maize (corn)",
  whea = "Wheat",
  rice = "Rice",
  barl = "Barley",
  sorg = "Sorghum",
  mill = "Millet",
  soyb = "Soya beans",
  grou = "Groundnuts, excluding shelled",
  cnut = "Coconuts, in shell",
  cott = "Cotton lint, ginned",
  sunf = "Sunflower seed",
  rape = "Rape or colza seed",
  pota = "Potatoes",
  cass = "Cassava, fresh",
  yams = "Yams",
  swpo = "Sweet potatoes",
  bana = "Bananas",
  sugc = "Sugar cane"
)

# Load ERA_Derived.rda from the first existing candidate path and return the
# yield-outcome subset with an RR column, or NULL when no file is found.
cca_load_era_yield <- function(candidates) {
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) return(NULL)
  env <- new.env(); load(hit, envir = env)
  era <- NULL
  for (nm in ls(env)) {
    x <- get(nm, envir = env)
    if (inherits(x, "data.frame") || inherits(x, "data.table")) { era <- x; break }
  }
  if (is.null(era)) return(NULL)
  if (exists("cca_log", mode = "function")) {
    cca_log(sprintf("ERA table (%s): %d rows", basename(hit), nrow(era)))
  }
  cca_era_yield_rows(era)
}

# Yield-outcome rows with positive means, plus RR = log(MeanT / MeanC).
cca_era_yield_rows <- function(era) {
  era <- as.data.frame(era)
  yield_mask <- grepl("yield|productivity", era$Out.SubInd, ignore.case = TRUE)
  era_y <- era[yield_mask &
               !is.na(era$MeanT) & !is.na(era$MeanC) &
               era$MeanT > 0 & era$MeanC > 0, , drop = FALSE]
  era_y$RR <- log(era_y$MeanT / era_y$MeanC)
  era_y
}

# Pooled ERA yield effect for one (CCA practice, SPAM crop) pair.
# Returns NULL when the practice is unmapped or no ERA rows match; otherwise a
# list with pc (= g, relative yield change), n_obs, ci_low, ci_high.
cca_era_yield_effect <- function(era_y, cca_practice, spam_code,
                                 cca_to_era = CCA_TO_ERA,
                                 spam_to_fao = SPAM_TO_FAO) {
  era_labels <- cca_to_era[[cca_practice]]
  if (is.null(era_labels) || is.null(era_y) || !nrow(era_y)) return(NULL)
  fao_name <- spam_to_fao[spam_code]
  if (is.na(fao_name)) fao_name <- spam_code
  rows <- era_y[vapply(era_y$PrName, function(x) {
                  any(vapply(era_labels, function(lbl) grepl(lbl, x, fixed = TRUE),
                             logical(1)))
                }, logical(1)) &
                (grepl(spam_code, era_y$Product, ignore.case = TRUE) |
                 grepl(fao_name,   era_y$Product, ignore.case = TRUE)), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  mu  <- mean(rows$RR, na.rm = TRUE)
  sdv <- stats::sd(rows$RR, na.rm = TRUE)
  n   <- sum(!is.na(rows$RR))
  ci  <- 1.96 * (sdv / sqrt(max(n, 1)))
  list(
    pc = exp(mu) - 1,
    n_obs = n,
    ci_low  = exp(mu - ci) - 1,
    ci_high = exp(mu + ci) - 1
  )
}
