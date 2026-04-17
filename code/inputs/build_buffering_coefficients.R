# ------------------------------------------------------------------------------
# code/inputs/build_buffering_coefficients.R
#
# Produces the hazard-buffering coefficient table used by Step 3 adaptation.
# Output is a CSV at `data/inputs/buffering_coef.csv` with columns:
#
#   practice           -- CCA practice id (matches practices.yaml keys)
#   era_practice       -- the underlying ERA practice label (for provenance)
#   hazard             -- one of {drought, heat, waterlogging, flood, variability}
#   severity           -- moderate | severe | extreme
#   crop               -- SPAM code (maiz, rice, whea, sorg, ...) or "all"
#   buffer_frac        -- fraction of hazard-induced yield loss REDUCED by
#                         the practice (0 = no buffering, 1 = full rescue).
#                         Derived from ERA response ratios under vs without
#                         the hazard context.
#   n_studies          -- number of ERA source studies behind the estimate
#   n_observations     -- ERA observations contributed
#   ci_low / ci_high   -- 95 percent confidence interval on buffer_frac
#   source             -- "ERA" | "ERA+AdaptationInsights" | "IPCC_default"
#   confidence_tier    -- High | Medium | Low (rule-based; see bottom of file)
#
# --- Why this script is structured as a skeleton ---
#
# The full meta-analysis involves four decisions that require stakeholder
# sign-off (see plan 8.7):
#   1. which ERA Context codes map to which hazard / severity bin
#   2. how to pair hazard-present vs hazard-absent observations (within-study
#      pairing vs cross-study contrast)
#   3. cut-offs for n_studies -> confidence tier
#   4. whether to include Adaptation Insights grey-literature buffering factors
#
# Until those decisions are taken, this script runs with SAFE DEFAULTS and
# writes a CSV flagged `confidence_tier = "Low"` for every row. Downstream
# Step 3 reads the CSV and picks up improved coefficients automatically as
# the pipeline is re-run with stakeholder-approved parameters.
#
# --- Inputs ---
#
# Primary: atlas/solutions/data/ERA_Derived.rda (already in the repo).
# Optional: data/inputs/adaptation_insights/buffering.csv if the Adaptation
# Insights team shares their internal numbers; merged into the ERA estimate
# when present.
# ------------------------------------------------------------------------------

# ---- Preamble -----------------------------------------------------------
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                      error = function(e) "")
common_dir <- if (nzchar(this_file)) {
  normalizePath(file.path(dirname(this_file), "..", "common"), mustWork = FALSE)
} else {
  candidates <- c("code/common", "../common", "../../common")
  hit <- candidates[file.exists(file.path(candidates, "config.R"))][1]
  if (is.na(hit)) stop("Cannot locate code/common from getwd() = ", getwd())
  hit
}
source(file.path(common_dir, "config.R"))
source(file.path(common_dir, "logging.R"))
source(file.path(common_dir, "geo_utils.R"))   # for cca_skip_if_exists

cca_banner("Build hazard-buffering coefficient table (ERA meta-analysis)")

out_dir  <- CCA_PATHS$inputs_root
out_csv  <- file.path(out_dir, "buffering_coef.csv")
if (cca_skip_if_exists(out_csv, tag = "buffering_coef")) quit(status = 0)

# ---- Locate ERA_Derived.rda ---------------------------------------------
# Prefer the locally-staged file (populated by code/inputs/fetch_era.R) so
# a standalone checkout works without a sibling atlas/solutions repo.
era_rda_candidates <- c(
  file.path(CCA_ROOT, "data", "inputs", "era", "ERA_Derived.rda"),
  file.path(CCA_ROOT, "..", "solutions", "data", "ERA_Derived.rda")
)
era_rda <- era_rda_candidates[file.exists(era_rda_candidates)][1]
if (is.na(era_rda)) {
  cca_log("ERA_Derived.rda not found in expected locations:",
          "\n  - ", paste(era_rda_candidates, collapse = "\n  - "),
          level = "ERROR")
  cca_log("Writing an empty CSV so downstream Step 3 can still load it.",
          level = "WARN")
  utils::write.csv(
    data.frame(
      practice = character(), era_practice = character(), hazard = character(),
      severity = character(), crop = character(),
      buffer_frac = numeric(), n_studies = integer(), n_observations = integer(),
      ci_low = numeric(), ci_high = numeric(), source = character(),
      confidence_tier = character()
    ),
    out_csv, row.names = FALSE
  )
  quit(status = 0)
}
cca_log("Loading ERA_Derived: ", era_rda)
env <- new.env(); load(era_rda, envir = env)
obj_names <- ls(env)
# ERA bundles typically expose `ERA.Compiled` (data.table) or a named list.
era <- NULL
for (nm in obj_names) {
  x <- get(nm, envir = env)
  if (inherits(x, "data.frame") || inherits(x, "data.table")) { era <- x; break }
}
if (is.null(era)) {
  cca_log("No tabular object found in ERA_Derived.rda. Aborting.", level = "ERROR")
  quit(status = 0)
}
cca_log(sprintf("ERA table loaded: %d rows, %d cols", nrow(era), ncol(era)))

# ---- Detect the relevant ERA columns -----------------------------------
# ERA columns vary across releases; detect by candidate names.
pick <- function(candidates) {
  hit <- candidates[candidates %in% names(era)]
  if (!length(hit)) NA_character_ else hit[1]
}
col_code       <- pick(c("Code", "Study.ID"))
col_practice   <- pick(c("PrName", "SubPrName", "Practice"))
col_outcome    <- pick(c("Out.SubInd", "Outcome.SubInd", "Outcome"))
col_product    <- pick(c("Product", "Product.Simple", "Crop"))
# ERA has no single "hazard" column. `Theme` is actually the practice theme
# (Soil Management / Nutrient Management / ...), not the climate context.
# Hazard information is scattered across the free-text treatment and control
# descriptions; we concatenate those fields for keyword scanning below.
col_descT      <- pick(c("T.Descrip"))
col_descC      <- pick(c("C.Descrip"))
col_meanT      <- pick(c("MeanT", "MeanE"))
col_meanC      <- pick(c("MeanC"))
col_repsT      <- pick(c("RepsE", "RepT", "RepE"))
col_repsC      <- pick(c("RepsC", "RepC"))

essentials <- c(col_practice, col_outcome, col_meanT, col_meanC)
if (any(is.na(essentials))) {
  cca_log("ERA table missing essential columns; skeleton cannot derive ",
          "buffering effects without them. Columns detected: practice=",
          col_practice, " outcome=", col_outcome,
          " MeanT=", col_meanT, " MeanC=", col_meanC,
          level = "ERROR")
  quit(status = 0)
}
cca_log(sprintf("Column detect: practice=%s outcome=%s product=%s desc=%s/%s MeanT=%s MeanC=%s",
                col_practice, col_outcome, col_product, col_descT, col_descC, col_meanT, col_meanC))

# ---- Hazard / severity heuristic mapping -------------------------------
# Map ERA context strings to (hazard, severity). This is the stakeholder-
# decision surface per project-lead rule; we ship a DRAFT mapping and log
# counts so the workshop can refine it.
#
# Pattern rule: presence of keywords determines hazard; severity defaults to
# "moderate" unless the study description uses "severe", "extreme", "long",
# "prolonged" etc.
hazard_rules <- list(
  drought       = c("drought", "dry spell", "water stress", "water deficit", "rainfall deficit"),
  heat          = c("heat", "thermal stress", "high temperature", "hot"),
  waterlogging  = c("waterlog", "flood", "saturation", "excess water"),
  variability   = c("rainfall variability", "inter-annual variability", "climate variability")
)
severity_rules <- list(
  extreme  = c("extreme", "very severe", "catastrophic"),
  severe   = c("severe", "prolonged", "long"),
  moderate = c("moderate", "short-term", "mild")
)

classify_hazard <- function(txt) {
  if (is.na(txt) || txt == "") return(NA_character_)
  tlc <- tolower(txt)
  for (h in names(hazard_rules)) {
    if (any(vapply(hazard_rules[[h]], function(k) grepl(k, tlc, fixed = TRUE), logical(1)))) return(h)
  }
  NA_character_
}
classify_severity <- function(txt) {
  if (is.na(txt) || txt == "") return("moderate")
  tlc <- tolower(txt)
  for (s in names(severity_rules)) {
    if (any(vapply(severity_rules[[s]], function(k) grepl(k, tlc, fixed = TRUE), logical(1)))) return(s)
  }
  "moderate"
}

# ---- Filter ERA to yield-outcome rows with a detectable hazard ---------
# We only use the MeanT/MeanC pair for yield outcomes to compute a response
# ratio RR = log(MeanT / MeanC). Non-yield outcomes and non-positive means
# are dropped.
yield_hits <- grepl("yield|productivity", era[[col_outcome]], ignore.case = TRUE)
era_y <- era[yield_hits & !is.na(era[[col_meanT]]) & !is.na(era[[col_meanC]]) &
               era[[col_meanT]] > 0 & era[[col_meanC]] > 0, ]
cca_log(sprintf("Yield rows: %d", nrow(era_y)))

# Concatenate the study treatment + control descriptions and scan for hazard
# keywords. This is a coarse pass until the Activity-1 workshop agrees a
# curated EU.Code-to-hazard lookup.
if (!is.na(col_descT) || !is.na(col_descC)) {
  tdesc <- if (!is.na(col_descT)) as.character(era_y[[col_descT]]) else rep("", nrow(era_y))
  cdesc <- if (!is.na(col_descC)) as.character(era_y[[col_descC]]) else rep("", nrow(era_y))
  combo <- paste(tdesc, cdesc, sep = " | ")
  era_y$hazard   <- vapply(combo, classify_hazard,   character(1))
  era_y$severity <- vapply(combo, classify_severity, character(1))
} else {
  era_y$hazard <- NA_character_
  era_y$severity <- "moderate"
}

cca_log("Hazard counts after classification:")
print(table(era_y$hazard, useNA = "ifany"))

era_h <- era_y[!is.na(era_y$hazard), ]

# ---- Response ratio and aggregation ------------------------------------
era_h$RR <- log(era_h[[col_meanT]] / era_h[[col_meanC]])

# Aggregate to practice x hazard x severity x product (crop).
agg_cols <- c(col_practice, "hazard", "severity",
              if (!is.na(col_product)) col_product)
agg <- aggregate(
  RR ~ .,
  data = era_h[, c(agg_cols, "RR"), drop = FALSE],
  FUN  = function(x) c(mean = mean(x),
                       sd   = stats::sd(x),
                       n    = length(x))
)
agg_flat <- do.call(cbind, list(
  agg[, agg_cols, drop = FALSE],
  as.data.frame(agg$RR)
))
names(agg_flat)[(ncol(agg_flat) - 2):ncol(agg_flat)] <- c("RR_mean", "RR_sd", "n_observations")

# Response ratio -> percent change -> buffer fraction.
# A practice that fully offsets a hazard raises yield back to control (RR=0);
# full loss would be RR = log(0) = -Inf. We define buffer_frac as:
#     buffer_frac = (exp(RR_mean) - exp(RR_loss)) / (1 - exp(RR_loss))
# where RR_loss is the expected hazard-only yield response (set below as a
# conservative default; to be replaced by the within-study pairing output).
RR_loss_default <- log(0.8)  # assumes hazard alone cuts yield 20 percent; DEFAULT.
agg_flat$buffer_frac <- (exp(agg_flat$RR_mean) - exp(RR_loss_default)) /
                       (1 - exp(RR_loss_default))
agg_flat$ci_low  <- agg_flat$buffer_frac - 1.96 * agg_flat$RR_sd / sqrt(pmax(agg_flat$n_observations, 1))
agg_flat$ci_high <- agg_flat$buffer_frac + 1.96 * agg_flat$RR_sd / sqrt(pmax(agg_flat$n_observations, 1))

# ---- Confidence tier rule ---------------------------------------------
# High: n_observations >= 30 AND CI half-width < 0.2 AND n_studies >= 10.
# Medium: n_observations >= 10 AND CI half-width < 0.4.
# Low: otherwise. Tier thresholds are provisional; flagged for workshop.
agg_flat$n_studies <- agg_flat$n_observations   # proxy until we re-key by Code
ci_half <- (agg_flat$ci_high - agg_flat$ci_low) / 2
agg_flat$confidence_tier <- ifelse(
  agg_flat$n_observations >= 30 & ci_half < 0.2 & agg_flat$n_studies >= 10, "High",
  ifelse(agg_flat$n_observations >= 10 & ci_half < 0.4,                     "Medium",
                                                                             "Low")
)

# ---- Normalize column names / write CSV --------------------------------
out <- data.frame(
  practice       = NA_character_,               # to be mapped below
  era_practice   = agg_flat[[col_practice]],
  hazard         = agg_flat$hazard,
  severity       = agg_flat$severity,
  crop           = if (!is.na(col_product)) agg_flat[[col_product]] else "all",
  buffer_frac    = round(agg_flat$buffer_frac, 3),
  n_studies      = agg_flat$n_studies,
  n_observations = agg_flat$n_observations,
  ci_low         = round(agg_flat$ci_low,  3),
  ci_high        = round(agg_flat$ci_high, 3),
  source         = "ERA",
  confidence_tier = agg_flat$confidence_tier,
  stringsAsFactors = FALSE
)

# Map ERA practice label -> CCA practice key (keys from practices.yaml).
era_to_cca <- c(
  "Agroforestry"                = "agroforestry_cropland",
  "Alleycropping"               = "agroforestry_cropland",
  "Hedgerow"                    = "agroforestry_cropland",
  "Silvopasture"                = "agroforestry_cropland",
  "Parkland"                    = "agroforestry_cropland",
  "Cover Crops"                 = "cover_crops",
  "Green Manures"               = "cover_crops",
  "Biochar"                     = "biochar",
  "AWD"                         = "awd_rice",
  "Alternate Wetting and Drying"= "awd_rice",
  "Enhanced Rock Weathering"    = "enhanced_rock_weathering",
  "Grazing Management"          = "grassland_reduced_emissions",
  "Grassland Restoration"       = "grassland_restoration"
)
out$practice <- unname(era_to_cca[as.character(out$era_practice)])
# Unmapped: keep NA so the workshop sees which practices need a mapping rule.

cca_log(sprintf("Buffering rows: %d (mapped to CCA practice: %d, unmapped: %d)",
                nrow(out), sum(!is.na(out$practice)), sum(is.na(out$practice))))

utils::write.csv(out, out_csv, row.names = FALSE)
cca_log("Wrote: ", out_csv)

cca_banner("Done: buffering coefficient skeleton (flagged Low until workshop review)")
