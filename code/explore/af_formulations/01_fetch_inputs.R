# ------------------------------------------------------------------------------
# code/explore/af_formulations/01_fetch_inputs.R
#
# EXPLORATORY. Step 5.2 of the adaptation-factor formulation brief: pull the
# public inputs for KEN, ETH, NGA into data/staging/af_formulations/ and write
# a MANIFEST.csv recording source, pin and checksum of every file.
#
# Run from the repo root:  Rscript code/explore/af_formulations/01_fetch_inputs.R
#
# Inputs (all public):
#   * Atlas hazard x exposure tiers, variable=prod_t, period=jagermeyr,
#     ENSEMBLEmean, 3 severities -- pinned by S3 versionId (pins from
#     hazards_prototype/HANDOVER_2026-10-12_carbon-credits-adaptation-hazard-
#     briefing.md S7). Filtered to the three countries and the crops below; all
#     admin levels, scenarios and hazard_vars kept.
#   * Atlas exposure denominator crop-livestock_all_res-25.parquet (pinned),
#     rows with unit t or ha for the same crops.
#   * GAUL 2024 admin-2 boundaries. A versioned GET on this key returns 403,
#     so the unversioned key is read and the script stops unless its current
#     x-amz-version-id equals the pin.
#   * HarvestChoice AEZ16 r2.0 TIF (Harvard Dataverse doi:10.7910/DVN/M7XIUB,
#     CC0, file id 2723111).
#   * Dixon farming systems FS1/FS2 rasters (Dixon et al. 2019/2021, FAO
#     "Farming Systems and Food Security in Africa", level-1 and level-2
#     subsystems, 5 arcmin). Copied from a local checkout of AICCRA_EFA
#     (AF_FS_DIR, default ../AICCRA_EFA/Data/Farming_systems); not
#     downloadable here. A collaborator needs that folder or the FAO source.
#   * WorldClim 2.1 10' bio_1 and prec_01..12 via geodata::worldclim_global
#     into data/inputs/soilgrids/climate/ (same cache the biochar step uses;
#     skipped when present).
#   * MapSPAM 2020 (Atlas harmonised) maize production and harvested area,
#     windowed /vsicurl read of the SSA COGs. The object's version id at
#     fetch time is stored in <file>.versionId and reused on later runs.
#
# ERA_Derived.rda is fetched separately by code/inputs/fetch_era.R.
#
# Re-running skips files that already exist; delete a file to re-fetch it.
# Downloads go to <file>.part and are renamed on completion.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

if (!dir.exists("code/explore/af_formulations")) {
  stop("Run from the cc4adapt repo root (getwd() = ", getwd(), ")")
}

ISO3  <- c("KEN", "ETH", "NGA")
CROPS <- c("maize", "sorghum", "wheat", "bean", "pearl-millet", "small-millet")
OUT   <- "data/staging/af_formulations"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

S3 <- "https://digital-atlas.s3.amazonaws.com/"
TIER_KEY <- paste0(
  "domain=hazard_exposure/source=nex-gddp-cmip6/region=ssa/",
  "processing=hazard-risk-exposure/variable=prod_t/period=jagermeyr/",
  "model=ENSEMBLEmean/severity=%s/int=multi-hazard.parquet"
)
TIER_PINS <- c(
  moderate = "kPR4BCorjzg86rnT2hPGA2M9j0JHY57s",
  severe   = "ckaDGRQzyM3gbvrk91_.fDu4a72K2rxi",
  extreme  = "oJRJxUNbBU37rRPBCPU8qpwUC0K28qqk"
)
DENOM_KEY <- paste0(
  "domain=exposure/type=combined/source=glw4-2020_spam2020AA/region=ssa/",
  "processing=atlas-harmonized/variable=crop-livestock_all_res-25.parquet"
)
DENOM_PIN <- "IolrzLCVLfGzd5Ryw1eACOeVp8Sy179Q"
GAUL_KEY <- paste0(
  "domain=boundaries/type=admin/source=gaul2024/region=africa/",
  "processing=analysis-ready/level=adm2/atlas_gaul24_a2_africa.parquet"
)
GAUL_PIN <- "66c7MiKvkW8w_4q6wfQXly0wy7FRKEIa"
AEZ_URL  <- "https://dataverse.harvard.edu/api/access/datafile/2723111"
FS_DIR   <- normalizePath(Sys.getenv("AF_FS_DIR", "../AICCRA_EFA/Data/Farming_systems"), mustWork = FALSE)
AEZ_MD5  <- "8e59b8bb5e19f19b6ce536f6f5ea8acb"  # Dataverse file 2723111, checked 2026-10-27
WC_DIR   <- "data/inputs/soilgrids/climate"

manifest <- list()
log_file <- function(item, source, pin, path, note = "") {
  manifest[[length(manifest) + 1]] <<- data.frame(
    item = item, source = source, pin = pin, path = path,
    bytes = file.size(path), md5 = unname(tools::md5sum(path)),
    fetched_at = format(file.mtime(path), "%Y-%m-%dT%H:%M:%S%z"), note = note
  )
  message(sprintf("  %-28s %s (%.1f MB)", item, basename(path), file.size(path) / 1e6))
}

con <- dbConnect(duckdb::duckdb())
invisible(dbExecute(con, "INSTALL httpfs; LOAD httpfs;"))

s3_version <- function(key) {
  hdr <- curlGetHeaders(paste0(S3, key))
  trimws(sub("(?i)^x-amz-version-id:\\s*", "",
             grep("(?i)^x-amz-version-id", hdr, value = TRUE, perl = TRUE), perl = TRUE))
}
download_part <- function(url, path) {
  tmp <- paste0(path, ".part")
  utils::download.file(url, tmp, mode = "wb", quiet = TRUE)
  invisible(file.rename(tmp, path))
}
q   <- function(x) as.character(dbQuoteString(con, x))
inq <- function(x) paste(vapply(x, q, ""), collapse = ", ")

copy_query <- function(sql, path) {
  tmp <- paste0(path, ".part")
  dbExecute(con, sprintf("COPY (%s) TO %s (FORMAT parquet, COMPRESSION zstd)", sql, q(tmp)))
  invisible(file.rename(tmp, path))
}

# ---- 1. Hazard x exposure tiers ------------------------------------------
message("Atlas prod_t tiers (jagermeyr)")
for (sev in names(TIER_PINS)) {
  path <- file.path(OUT, sprintf("atlas_prod_t_jagermeyr_%s.parquet", sev))
  url  <- paste0(S3, sprintf(TIER_KEY, sev), "?versionId=", TIER_PINS[[sev]])
  if (!file.exists(path)) {
    copy_query(sprintf(
      "SELECT * FROM read_parquet(%s, hive_partitioning = false)
       WHERE iso3 IN (%s) AND crop IN (%s)", q(url), inq(ISO3), inq(CROPS)), path)
  }
  log_file(paste0("prod_t_", sev), paste0(S3, sprintf(TIER_KEY, sev)), TIER_PINS[[sev]], path)
}

# ---- 2. Exposure denominator ---------------------------------------------
message("Exposure denominator")
path <- file.path(OUT, "atlas_denominator_res-25_t_ha.parquet")
if (!file.exists(path)) {
  copy_query(sprintf(
    "SELECT * FROM read_parquet(%s, hive_partitioning = false)
     WHERE iso3 IN (%s) AND crop IN (%s) AND unit IN ('t', 'ha')",
    q(paste0(S3, DENOM_KEY, "?versionId=", DENOM_PIN)), inq(ISO3), inq(CROPS)), path)
}
log_file("denominator_t_ha", paste0(S3, DENOM_KEY), DENOM_PIN, path)

# ---- 3. GAUL admin-2 (unversioned read, version checked) -----------------
message("GAUL 2024 admin-2")
path <- file.path(OUT, "gaul24_adm2_ken_eth_nga.parquet")
if (!file.exists(path)) {
  ver <- s3_version(GAUL_KEY)
  if (!identical(ver, GAUL_PIN)) stop("GAUL current version ", ver, " != pin ", GAUL_PIN)
  iso_col <- intersect(c("iso3", "ISO3", "iso3_code"),
                       dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet(%s)", q(paste0(S3, GAUL_KEY))))$column_name)[1]
  if (is.na(iso_col)) stop("No iso3 column in GAUL parquet")
  copy_query(sprintf("SELECT * FROM read_parquet(%s) WHERE %s IN (%s)",
                     q(paste0(S3, GAUL_KEY)), iso_col, inq(ISO3)), path)
}
log_file("gaul24_adm2", paste0(S3, GAUL_KEY), GAUL_PIN, path,
         "versioned GET is 403; current version checked equal to pin")

# ---- 4. HarvestChoice AEZ16 ----------------------------------------------
message("HarvestChoice AEZ")
zip <- file.path(OUT, "aez_r2.0_tif.zip")
if (!file.exists(zip)) download_part(AEZ_URL, zip)
if (unname(tools::md5sum(zip)) != AEZ_MD5) stop("AEZ zip md5 mismatch; delete ", zip, " and re-run")
aez_dir <- file.path(OUT, "aez")
# The Dataverse zip holds one zip per scheme (AEZ5, AEZ8, AEZ16); unpack AEZ16.
if (!length(list.files(aez_dir, pattern = "^AEZ16.*\\.tif$", recursive = TRUE))) {
  utils::unzip(zip, exdir = aez_dir)
  inner <- list.files(aez_dir, pattern = "^AEZ16.*\\.zip$", full.names = TRUE)
  if (length(inner) != 1) stop("Expected one AEZ16 zip inside ", zip)
  utils::unzip(inner, exdir = aez_dir)
  unlink(file.path(aez_dir, "__MACOSX"), recursive = TRUE)
}
log_file("aez_r2.0_zip", "doi:10.7910/DVN/M7XIUB (file 2723111)", "r2.0 / dataset v3.1", zip)

# ---- 5. Dixon farming systems (local copy) -------------------------------
message("Dixon farming systems")
fs_files <- c("FS1.tif", "FS1.tif.aux.xml", "FS2.tif", "FS2.tif.aux.xml", "FS_Metadata.csv")
if (!all(file.exists(file.path(FS_DIR, fs_files)))) stop("Farming-system files missing in ", FS_DIR)
fs_out <- file.path(OUT, "farming_systems")
dir.create(fs_out, showWarnings = FALSE)
for (f in fs_files) {
  dst <- file.path(fs_out, f)
  if (!file.exists(dst)) file.copy(file.path(FS_DIR, f), dst)
  if (!grepl("aux.xml$", f)) log_file(paste0("dixon_", f), file.path(FS_DIR, f), "local copy", dst,
                                      "Dixon et al. FAO farming systems of Africa, 2016 update; not fetchable here")
}

# ---- 5b. WorldClim 2.1 10' (bio_1, monthly prec) --------------------------
message("WorldClim 2.1 10'")
wc_files <- file.path(WC_DIR, "wc2.1_10m", c("wc2.1_10m_bio_1.tif", sprintf("wc2.1_10m_prec_%02d.tif", 1:12)))
if (!all(file.exists(wc_files))) {
  dir.create(WC_DIR, recursive = TRUE, showWarnings = FALSE)
  invisible(geodata::worldclim_global(var = "bio",  res = 10, path = WC_DIR))
  invisible(geodata::worldclim_global(var = "prec", res = 10, path = WC_DIR))
}
if (!all(file.exists(wc_files))) stop("WorldClim files missing under ", WC_DIR)
log_file("worldclim_bio_1", "WorldClim 2.1 via geodata::worldclim_global", "v2.1", wc_files[1])
log_file("worldclim_prec_01", "WorldClim 2.1 via geodata::worldclim_global", "v2.1", wc_files[2],
         "prec_01..12 all present; one logged")

# ---- 6. MapSPAM 2020 (Atlas harmonised) maize production and area ---------
# Windowed /vsicurl read of the maize layer over the three countries; the
# object's current x-amz-version-id is recorded as its pin.
message("MapSPAM maize rasters (Atlas harmonised COGs)")
SPAM_KEYS <- c(
  prod_t       = "domain=exposure/type=crop/source=spam2020v1r2_ssa/region=ssa/processing=atlas-harmonized/variable=prod_t/spam_prod_t_all.tif",
  harv_area_ha = "domain=exposure/type=crop/source=spam2020v1r2_ssa/region=ssa/processing=atlas-harmonized/variable=harv-area_ha/spam_harv-area_ha_all.tif"
)
gaul <- sf::st_as_sf(arrow::read_parquet(file.path(OUT, "gaul24_adm2_ken_eth_nga.parquet")))
ext3 <- terra::ext(terra::vect(gaul)) + 0.5
for (v in names(SPAM_KEYS)) {
  path  <- file.path(OUT, sprintf("spam2020_maize_%s.tif", v))
  vfile <- paste0(path, ".versionId")
  if (!file.exists(path)) {
    ver <- s3_version(SPAM_KEYS[[v]])
    r <- terra::rast(paste0("/vsicurl/", S3, SPAM_KEYS[[v]]))
    if (!"maize" %in% names(r)) stop("No maize layer in ", SPAM_KEYS[[v]])
    tmp <- paste0(path, ".part.tif")
    terra::writeRaster(terra::crop(r[["maize"]], ext3), tmp, overwrite = TRUE)
    invisible(file.rename(tmp, path))
    writeLines(ver, vfile)
  }
  ver <- if (file.exists(vfile)) readLines(vfile, warn = FALSE)[1] else "unknown (fetched before pins were stored)"
  log_file(paste0("spam_maize_", v), paste0(S3, SPAM_KEYS[[v]]), ver, path,
           "maize layer, cropped to KEN/ETH/NGA bbox + 0.5 deg; pin = version id at fetch")
}

# ---- Manifest ------------------------------------------------------------
m <- do.call(rbind, manifest)
utils::write.csv(m, file.path(OUT, "MANIFEST.csv"), row.names = FALSE)
message("Wrote ", file.path(OUT, "MANIFEST.csv"), " (", nrow(m), " files)")
dbDisconnect(con, shutdown = TRUE)
