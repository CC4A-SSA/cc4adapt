# ------------------------------------------------------------------------------
# code/explore/af_formulations/03b_era_site_years.R
#
# EXPLORATORY. ERA maize comparisons for the two practices, listed by site and
# measurement year, as the input to the season classifier (07_era_seasons.py).
#
# Run from the repo root:  Rscript code/explore/af_formulations/03b_era_site_years.R
#
# Output: data/staging/af_formulations/era_site_years.csv, one row per
# (practice, Site.Key, M.Year) with coordinates, year, ERA season slot
# (Season.Start, 1 or 2 where the site is bimodal), buffer radius and whether
# ERA holds planting / harvest dates. era_comparisons.csv keeps the comparison
# rows (MeanT, MeanC, control definition) for 08_era_split.R.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
if (!dir.exists("code/explore/af_formulations")) stop("Run from the cc4adapt repo root")

cfg <- yaml::read_yaml("code/explore/af_formulations/config/params.yaml")
pr  <- vapply(cfg$practices, `[[`, "", "era_prname")

era_file <- "data/inputs/era/ERA_Derived.rda"
env <- new.env(); nm <- load(era_file, envir = env)
# Pin the ERA input: md5 and the fetch record, so a re-run can confirm the same file.
src <- if (file.exists(paste0(era_file, ".source.json"))) jsonlite::read_json(paste0(era_file, ".source.json")) else list()
jsonlite::write_json(list(file = era_file, md5 = unname(tools::md5sum(era_file)), bytes = file.size(era_file),
                          object = nm[1], versions = sort(unique(as.character(get(nm[1], envir = env)$Version))), source = src),
                     "data/staging/af_formulations/era_pin.json", auto_unbox = TRUE, pretty = TRUE)
era <- as.data.table(get(nm[1], envir = env))
m <- era[Product.Simple == "Maize" & Out.SubInd == "Crop Yield" & PrName %in% pr &
           (is.na(Invalid.Comparison) | Invalid.Comparison != "Y")]
m[, `:=`(MeanT = as.numeric(MeanT), MeanC = as.numeric(MeanC))]
m <- m[is.finite(MeanT) & is.finite(MeanC) & MeanT > 0 & MeanC > 0]
m[, practice := names(pr)[match(PrName, pr)]]
m[, yr := M.Year.Start]
# Yields to Mg/ha (ERA holds a minority of rows in kg/ha); irrigation flagged so
# 08 can exclude trials whose water stress the climate record cannot see.
stopifnot(all(m$Units %in% c("Mg/ha", "kg/ha", NA)))
m[Units == "kg/ha", `:=`(MeanC = MeanC / 1000, MeanT = MeanT / 1000, Units = "Mg/ha")]
m[is.na(Units), Units := "Mg/ha"]
m[, irrigated := (Irrigation.C %in% TRUE) | (Irrigation.T %in% TRUE)]

sy <- m[, .(lat = Latitude[1], lon = Longitude[1], buffer = Buffer[1], iso = ISO.3166.1.alpha.3[1],
            yr = yr[1], season_idx = Season.Start[1], site_type = Site.Type[1],
            has_plant = any(!is.na(Plant.Start)), has_harv = any(!is.na(Harvest.End)),
            plant_start = as.character(min(Plant.Start, na.rm = TRUE)),
            harvest_end = as.character(max(Harvest.End, na.rm = TRUE)),
            n_comp = .N, n_studies = uniqueN(Code)),
        by = .(PrName, practice, Site.Key, M.Year)]
sy[plant_start == "Inf" | plant_start == "-Inf", plant_start := NA_character_]
sy[harvest_end == "Inf" | harvest_end == "-Inf", harvest_end := NA_character_]

OUT <- "data/staging/af_formulations"
fwrite(sy, file.path(OUT, "era_site_years.csv"))
fwrite(m[, .(practice, PrName, SubPrName, PrName.Base, Code, Site.Key, M.Year, yr, season_idx = Season.Start,
             Latitude, Longitude, MeanC, MeanT, Units, irrigated, Irrigation.C, Irrigation.T, Site.Type, AEZ16simple)],
       file.path(OUT, "era_comparisons.csv"))
cat(sprintf("era_site_years.csv: %d practice x site x year rows (%d unique site-years, %d from 1981), %d sites, %d studies; era_comparisons.csv: %d rows (%d irrigated)\n",
            nrow(sy), uniqueN(sy[, .(Site.Key, M.Year)]), uniqueN(sy[yr >= 1981, .(Site.Key, M.Year)]), uniqueN(sy$Site.Key), uniqueN(m$Code), nrow(m), sum(m$irrigated)))
