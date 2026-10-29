# ------------------------------------------------------------------------------
# code/explore/af_formulations/03c_era_site_windows.R
#
# EXPLORATORY. One growing-season window per ERA site x season slot, for
# 07_era_seasons.py. Needed because the GGCMI rainfed-maize calendar holds one
# season per pixel and, in bimodal East Africa, often the short rains, so it
# cannot serve ERA's slot 1 / slot 2 distinction.
#
# Run from the repo root:  Rscript code/explore/af_formulations/03c_era_site_windows.R
#
# Rule, in order of preference:
#   1. ERA planting and harvest dates for that site x slot (modal planting
#      month-day across its site-years; harvest likewise). source = "era_dates".
#   2. WorldClim 2.1 10' monthly precipitation climatology at the site: a month
#      is wet if P >= WET_MM; wet runs of >= MIN_RUN months, circular year.
#      Slot k = k-th wet run in calendar order starting from the first dry
#      month after the driest month. Window = first day of the run's first
#      month to the last day of the (run length + 1)-th month, capped at
#      SEASON_MAX_MONTHS, so a 3-month wet run gives a ~4-month maize season.
#      source = "climatology".
#   3. GGCMI nearest-point token (first season). source = "ggcmi".
# Sites whose ERA slot has no matching climatology run are flagged
# window_status = "slot_missing" and get the last available run.
#
# Output: data/staging/af_formulations/era_site_windows.csv
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(terra) })
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")

IN <- "data/staging/af_formulations"
WC <- "data/inputs/soilgrids/climate/wc2.1_10m"
WET_MM <- 60; MIN_RUN <- 2; SEASON_MAX_MONTHS <- 5
mdays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
tok <- function(m1, m2) sprintf("%02d-01:%02d-%02d", m1, m2, mdays[m2])

sy <- fread(file.path(IN, "era_site_years.csv"))
sy[, slot := fifelse(is.na(season_idx), 1L, as.integer(season_idx))]
sites <- unique(sy[, .(Site.Key, lat, lon)])

prec <- rast(file.path(WC, sprintf("wc2.1_10m_prec_%02d.tif", 1:12)))
pm <- terra::extract(prec, sites[, .(lon, lat)])[, -1]
colnames(pm) <- sprintf("p%02d", 1:12)

runs_for <- function(p) {
  if (anyNA(p)) return(list())
  wet <- p >= WET_MM
  if (all(wet)) return(list(c(1, 12)))
  if (!any(wet)) return(list())
  k <- which.min(p)                     # start scan after the driest month
  idx <- ((k:(k + 11)) - 1) %% 12 + 1
  w <- wet[idx]; r <- rle(w); out <- list(); pos <- 1
  for (i in seq_along(r$lengths)) {
    if (r$values[i] && r$lengths[i] >= MIN_RUN) {
      s <- idx[pos]; len <- r$lengths[i]
      e <- idx[min(pos + min(len, SEASON_MAX_MONTHS - 1), 12)]
      out[[length(out) + 1]] <- c(s, e)
    }
    pos <- pos + r$lengths[i]
  }
  # order runs by calendar start so slot 1 is the earlier season in the year
  out[order(vapply(out, `[`, 1, 1))]
}

win <- rbindlist(lapply(seq_len(nrow(sites)), function(i) {
  rs <- runs_for(as.numeric(pm[i, ]))
  slots <- sort(unique(sy[Site.Key == sites$Site.Key[i], slot]))
  rbindlist(lapply(slots, function(s) {
    if (length(rs) == 0) return(data.table(Site.Key = sites$Site.Key[i], slot = s, fixed_season = NA_character_,
                                           source = "none", window_status = "no_wet_run"))
    r <- rs[[min(s, length(rs))]]
    data.table(Site.Key = sites$Site.Key[i], slot = s, fixed_season = tok(r[1], r[2]), source = "climatology",
               window_status = fifelse(s <= length(rs), "ok", "slot_missing"))
  }))
}))
win[, n_wet_runs := vapply(Site.Key, function(k) length(runs_for(as.numeric(pm[match(k, sites$Site.Key), ]))), 1L)]

# ERA dates override where present for that site x slot
ed <- sy[!is.na(plant_start) & !is.na(harvest_end),
         .(era_plant = names(which.max(table(substr(plant_start, 6, 10)))),
           era_harv  = names(which.max(table(substr(harvest_end, 6, 10)))), n_dated = .N), by = .(Site.Key, slot)]
win <- merge(win, ed, by = c("Site.Key", "slot"), all.x = TRUE)
win[!is.na(era_plant), `:=`(fixed_season_clim = fixed_season, fixed_season = paste0(era_plant, ":", era_harv),
                            source = "era_dates", window_status = "ok")]

fwrite(win, file.path(IN, "era_site_windows.csv"))
cat(sprintf("era_site_windows.csv: %d site-slots; by source: %s; slot_missing: %d; no window: %d\n",
            nrow(win), paste(names(table(win$source)), table(win$source), collapse = ", "),
            sum(win$window_status == "slot_missing"), sum(is.na(win$fixed_season))))
print(win[, .N, by = .(slot, source, window_status)][order(slot, source)])
