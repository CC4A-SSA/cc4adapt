# ------------------------------------------------------------------------------
# code/explore/af_formulations/04_units.R
#
# EXPLORATORY. Unit table for the formulation comparison: one row per GAUL 2024
# admin-2 unit in KEN, ETH, NGA with maize baseline and stratum tags.
#
# Run from the repo root:  Rscript code/explore/af_formulations/04_units.R
#
# Per unit (polygon-cell overlap weights from terra::extract(exact = TRUE)):
#   maize_prod_t, maize_harv_ha   MapSPAM 2020 (Atlas harmonised), 0.05 deg
#   B                             maize_prod_t / maize_harv_ha (t/ha)
#   fs1, fs2                      dominant Dixon farming system (level 1 / 2)
#   aez                           dominant AEZ class: HarvestChoice AEZ16
#                                 collapsed to aridity (rule in hazards_prototype
#                                 R/misc/0.7_create_aez_temp_zones.R) x mean
#                                 annual temperature (WorldClim 2.1 10' bio_1:
#                                 <15 low, 15-25 medium, >25 high). The 0.7
#                                 script uses 2.5' bio_1; 10' is used here.
#   bimodal_share                 share of maize production in cells with a
#                                 bimodal rainfall climatology (see below)
# "Dominant" = largest share of maize production in the unit; units with no
# mapped maize fall back to largest area share and are flagged.
#
# Bimodal rule (crude, for flagging only): on the WorldClim 10' monthly
# precipitation climatology a month is wet if P >= 60 mm; a cell is bimodal
# if the circular year has >= 2 separate runs of wet months, each >= 2 months.
# This catches East African long/short rains and Ethiopian belg/kiremt; it
# does not catch a relative mid-season dip that stays above 60 mm.
#
# Output: data/staging/af_formulations/units.parquet
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}

IN <- "data/staging/af_formulations"
WC <- "data/inputs/soilgrids/climate/wc2.1_10m"
WET_MM <- 60
MIN_RUN <- 2

gaul <- sf::st_as_sf(arrow::read_parquet(file.path(IN, "gaul24_adm2_ken_eth_nga.parquet")))
if (is.na(sf::st_crs(gaul))) sf::st_crs(gaul) <- 4326
v <- terra::vect(gaul[, c("iso3", "admin1_name", "admin2_name", "gaul1_code", "gaul2_code")])

prod <- rast(file.path(IN, "spam2020_maize_prod_t.tif"))
harv <- rast(file.path(IN, "spam2020_maize_harv_area_ha.tif"))
tmpl <- prod

# ---- Bimodality from monthly rainfall ------------------------------------
prec <- crop(rast(file.path(WC, sprintf("wc2.1_10m_prec_%02d.tif", 1:12))), ext(tmpl) + 1)
is_bimodal <- function(p) {
  if (anyNA(p)) return(NA_integer_)
  wet <- p >= WET_MM
  if (all(wet) || !any(wet)) return(0L)
  # rotate so the series starts on a dry month, then count wet runs
  k <- which(!wet)[1]
  w <- c(wet[k:12], if (k > 1) wet[1:(k - 1)])
  r <- rle(w)
  as.integer(sum(r$values & r$lengths >= MIN_RUN) >= 2)
}
bimodal <- app(prec, is_bimodal)

# ---- AEZ class -----------------------------------------------------------
aez16 <- rast(list.files(file.path(IN, "aez"), pattern = "^AEZ16.*\\.tif$", full.names = TRUE)[1])
arid_rcl <- matrix(c(0, 1, 4, 1, 7, 1, 11, 1,
                     1, 2, 5, 2, 8, 2, 12, 2,
                     2, 3, 6, 3, 9, 3, 13, 3,
                     3, 4, 10, 4, 14, 4), ncol = 2, byrow = TRUE)
aridity <- classify(aez16, arid_rcl, others = NA)
mat <- crop(rast(file.path(WC, "wc2.1_10m_bio_1.tif")), ext(tmpl) + 1)
tclass <- classify(mat, matrix(c(-Inf, 15, 1, 15, 25, 2, 25, Inf, 3), ncol = 3, byrow = TRUE))
ARID_LAB <- c("Arid", "Semiarid", "Subhumid", "Humid")
TEMP_LAB <- c("Low Temp", "Medium Temp", "High Temp")

# ---- Farming systems -----------------------------------------------------
fs1 <- rast(file.path(IN, "farming_systems", "FS1.tif"))
fs2 <- rast(file.path(IN, "farming_systems", "FS2.tif"))
fs1_lab <- cats(fs1)[[1]]; fs2_lab <- cats(fs2)[[1]]
fs1 <- as.int(fs1); fs2 <- as.int(fs2)  # extract codes, relabel after

# ---- Stack on the 0.05 deg MapSPAM grid ----------------------------------
near <- function(r) resample(crop(r, ext(tmpl) + 1), tmpl, method = "near")
stk <- c(prod, harv, near(bimodal), near(aridity), near(tclass), near(fs1), near(fs2))
names(stk) <- c("prod", "harv", "bimodal", "arid", "temp", "fs1", "fs2")

cells <- as.data.table(terra::extract(stk, v, exact = TRUE))
cells[, `:=`(prod = fifelse(is.na(prod), 0, prod * fraction),
             harv = fifelse(is.na(harv), 0, harv * fraction))]

dominant <- function(code, w, a) {
  ok <- !is.na(code)
  if (!any(ok)) return(NA_integer_)
  ww <- if (sum(w[ok]) > 0) w[ok] else a[ok]
  s <- tapply(ww, code[ok], sum)
  as.integer(names(s)[which.max(s)])
}

units <- cells[, .(
  maize_prod_t   = sum(prod),
  maize_harv_ha  = sum(harv),
  n_cells        = .N,
  bimodal_share  = if (sum(prod[!is.na(bimodal)]) > 0)
                     sum((prod * bimodal)[!is.na(bimodal)]) / sum(prod[!is.na(bimodal)])
                   else mean(bimodal, na.rm = TRUE),
  arid           = dominant(arid, prod, fraction),
  temp           = dominant(temp, prod, fraction),
  fs1_code       = dominant(fs1, prod, fraction),
  fs2_code       = dominant(fs2, prod, fraction)
), by = ID]

attrs <- as.data.table(as.data.frame(v))[, ID := .I]
units <- merge(attrs, units, by = "ID", all.x = TRUE)
units[, B := fifelse(maize_harv_ha > 0, maize_prod_t / maize_harv_ha, NA_real_)]
units[, strata_weight := fifelse(maize_prod_t > 0, "maize_prod", "area")]
units[, aez := fifelse(is.na(arid) | is.na(temp), NA_character_,
                       paste(TEMP_LAB[temp], "-", ARID_LAB[arid]))]
units[, fs1 := fs1_lab[[2]][match(fs1_code, fs1_lab[[1]])]]
units[, fs2 := fs2_lab[[2]][match(fs2_code, fs2_lab[[1]])]]
units[, bimodal := bimodal_share >= 0.5]
units[, c("ID", "arid", "temp") := NULL]

arrow::write_parquet(units, file.path(IN, "units.parquet"))
message(sprintf("units.parquet: %d units; %d with maize; %d bimodal",
                nrow(units), sum(units$maize_prod_t > 0, na.rm = TRUE),
                sum(units$bimodal, na.rm = TRUE)))
print(units[maize_prod_t > 0, .(units = .N, prod_kt = round(sum(maize_prod_t) / 1e3),
                                 B_med = round(median(B, na.rm = TRUE), 2),
                                 bimodal = sum(bimodal, na.rm = TRUE)), by = iso3])
