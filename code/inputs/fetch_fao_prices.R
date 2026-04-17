# ------------------------------------------------------------------------------
# code/inputs/fetch_fao_prices.R
#
# Fetches FAOSTAT Producer Prices annual table and produces a tidy CSV of
# national farm-gate prices (USD/t) for the commodities that MapSPAM and
# ERA reference in Steps 2 and 3.
#
# --- Source ---
#
# FAOSTAT "Producer Prices - Annual" bulk download:
#   https://bulks-faostat.fao.org/production/Prices_E_All_Data_(Normalized).zip
# (Endpoint is the "normalized" flavour: long-form, one row per
# country x item x year x value.)
#
# Atlas convention (atlas/atlas_exposure/2_spam_convert_intdollar.r) uses a
# MapSPAM-compatible crop list and converts Int$2015 to USD using the LCU
# producer price. We reproduce that list here in a stable CSV so Step 2 /
# Step 3 downstream scripts can join without pulling the FAO bulk every
# time.
#
# --- Why refresh this (vs atlas defaults) ---
#
# The atlas pipeline was parameterised for 2017 prices. The project plan
# asks for 2022-2024 prices (the most recent 3-year window available in
# FAOSTAT) so NEB numbers are not a decade old.
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
# geo_utils.R defines cca_skip_if_exists() + the masking helpers. We need the
# skip helper to honour the pipeline-wide idempotency contract.
source(file.path(common_dir, "geo_utils.R"))
source(file.path(common_dir, "download.R"))

cca_banner("Fetch FAOSTAT producer prices (2022-2024 annual)")

raw_dir     <- file.path(CCA_PATHS$inputs_root, "faostat")
staging_dir <- file.path(CCA_PATHS$staging_root, "faostat")
dir.create(raw_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

staged_csv <- file.path(staging_dir, "fao_producer_prices_usd_per_tonne_2022_2024.csv")
if (cca_skip_if_exists(staged_csv, tag = "fao_prices")) quit(status = 0)

FAOSTAT_ZIP_URL <- "https://bulks-faostat.fao.org/production/Prices_E_All_Data_(Normalized).zip"
zip_path <- file.path(raw_dir, "Prices_E_All_Data_Normalized.zip")

cca_download(FAOSTAT_ZIP_URL, zip_path)

extract_dir <- file.path(raw_dir, "faostat_prices_extracted")
if (!dir.exists(extract_dir) || !length(list.files(extract_dir))) {
  utils::unzip(zip_path, exdir = extract_dir)
}
csvs <- list.files(extract_dir, pattern = "Prices_E_All_Data_\\(Normalized\\)\\.csv$",
                   full.names = TRUE)
if (!length(csvs)) {
  stop("[CCA] FAOSTAT price CSV not found after unzip.")
}

# The normalized CSV is ~300 MB; we filter on read to only keep the 3
# countries and the 2022-2024 year window. Using data.table::fread is 10x
# faster than read.csv; fall back to base if data.table missing.
use_dt <- requireNamespace("data.table", quietly = TRUE)
# FAOSTAT reporting lags vary by country -- Nigeria and Ethiopia rarely
# publish price data within 2 years of the current year, so we widen the
# window to the last 6 years and downstream code selects "most recent
# available per country x item" rather than a fixed year.
YEAR_WINDOW <- (as.integer(format(Sys.Date(), "%Y")) - 6L):as.integer(format(Sys.Date(), "%Y"))

if (use_dt) {
  prices <- data.table::fread(csvs[1], showProgress = FALSE)
  # FAOSTAT columns: Area Code, Area, Item Code, Item, Element Code, Element,
  # Year Code, Year, Unit, Value, Flag. We keep only Area in the 3 countries,
  # the USD/tonne element, and the widened year window.
  sel <- prices[
    Area %in% c("Nigeria", "Kenya", "Ethiopia") &
    Element == "Producer Price (USD/tonne)" &
    Year %in% YEAR_WINDOW,
    .(Area, Item, Year, Value, Unit)
  ]
} else {
  df <- utils::read.csv(csvs[1], stringsAsFactors = FALSE)
  keep <- df$Area %in% c("Nigeria", "Kenya", "Ethiopia") &
          df$Element == "Producer Price (USD/tonne)" &
          df$Year %in% YEAR_WINDOW
  sel <- df[keep, c("Area", "Item", "Year", "Value", "Unit")]
}

cca_log(sprintf("Filtered FAO prices rows: %d", nrow(sel)))

out_csv <- file.path(staging_dir, "fao_producer_prices_usd_per_tonne_2022_2024.csv")
utils::write.csv(sel, out_csv, row.names = FALSE)
cca_log("Wrote: ", out_csv)

# ---- Cross-check --------------------------------------------------------
# Spot-check: Kenya maize price should be in the range 200-500 USD/t; Nigeria
# rice 300-800 USD/t; Ethiopia teff 600-1200 USD/t. Large deviations likely
# indicate unit confusion (LCU/tonne rather than USD/tonne).
spot <- list(
  list(country = "Kenya",    item = "Maize (corn)", lo = 150, hi = 700),
  list(country = "Nigeria",  item = "Rice",         lo = 200, hi = 1200),
  list(country = "Ethiopia", item = "Maize (corn)", lo = 150, hi = 700)
)
for (s in spot) {
  # fixed = TRUE because FAO item labels contain regex metacharacters like
  # "Maize (corn)" -- without fixed the parens would be treated as a capture
  # group and never match. `ignore.case` is ignored with fixed=TRUE, so we
  # normalise both sides to tolower for a case-insensitive match.
  vals <- as.numeric(sel$Value[sel$Area == s$country &
                               grepl(tolower(s$item), tolower(sel$Item), fixed = TRUE)])
  vals <- vals[!is.na(vals)]
  if (!length(vals)) {
    cca_log(sprintf("Spot-check %s/%s: no values found", s$country, s$item),
            level = "WARN")
    next
  }
  mu <- mean(vals)
  lvl <- if (mu >= s$lo && mu <= s$hi) "INFO" else "WARN"
  cca_log(sprintf("Spot-check %s %s: mean USD/t = %.0f (expected %d-%d)",
                  s$country, s$item, mu, s$lo, s$hi),
          level = lvl)
}

cca_banner("Done: FAO producer prices staged")
