# ------------------------------------------------------------------------------
# code/common/config.R
#
# Central configuration for the CCA pipeline. Any constant that more than one
# script depends on belongs here, NOT duplicated downstream. This prevents the
# classic bug of "Step 3 re-projected to a slightly different grid than Step 1"
# because someone typed `res = 0.08333` in one place and `res = 5/60` in another.
#
# Design notes:
#   * CRS is EPSG:4326 to match MapSPAM, CHIRPS, Naturebase, CMIP6 hazards --
#     the upstream inputs we will resample FROM. Reprojecting to a projected
#     CRS would force a second resampling for every downstream merge; not worth
#     it for a pan-African / 3-country study at 5 arc-minutes.
#   * Resolution is expressed as the fraction 5/60 so it is numerically exact
#     (0.08333333... decimal literals introduce rounding at ~1e-8 that terra
#     compares as "different extents").
#   * ISO3 list is ordered NGA, KEN, ETH -- matches the investment document
#     sequencing (Activity 1 rollout, rapid assessment first).
# ------------------------------------------------------------------------------

# ---- Project root (resolved relative to this file so scripts can be invoked
# ---- from anywhere: `Rscript code/activity2/step1_mitigation.R` or
# ---- sourced from the repo root both work).
cca_project_root <- function() {
  # Resolve in priority order:
  #   1. Rscript --file=... argument (most reliable under Rscript)
  #   2. A sourced frame's $ofile attribute (source() in interactive R)
  #   3. getSrcDirectory() called inside a sourced file
  #   4. Walk up from cwd looking for a `ccadapt/` folder name
  cargs <- commandArgs(trailingOnly = FALSE)
  fhit <- grep("^--file=", cargs, value = TRUE)
  if (length(fhit)) {
    path <- sub("^--file=", "", fhit[1])
    if (nzchar(path)) {
      # config.R is at <root>/code/common/config.R -- go up two parents from
      # whichever script was actually invoked (it could be build_cca_template.R
      # or 00_fetch_all.R; both live at <root>/code/<subdir>/).
      script_dir <- normalizePath(dirname(path), mustWork = FALSE)
      # Walk up until we hit a directory named "ccadapt" (the repo root we want).
      wd <- script_dir
      while (basename(wd) != "ccadapt" && dirname(wd) != wd) wd <- dirname(wd)
      if (basename(wd) == "ccadapt") return(wd)
    }
  }
  # source()'d: the active top frame may have an $ofile attribute
  frames <- sys.frames()
  for (i in seq_along(frames)) {
    of <- tryCatch(frames[[i]]$ofile, error = function(e) NULL)
    if (!is.null(of) && nzchar(of)) {
      wd <- normalizePath(dirname(of), mustWork = FALSE)
      while (basename(wd) != "ccadapt" && dirname(wd) != wd) wd <- dirname(wd)
      if (basename(wd) == "ccadapt") return(wd)
    }
  }
  # Last resort: walk up from cwd.
  wd <- normalizePath(getwd(), mustWork = FALSE)
  while (basename(wd) != "ccadapt" && dirname(wd) != wd) wd <- dirname(wd)
  if (basename(wd) != "ccadapt") {
    stop("[CCA] Cannot locate the ccadapt project root. Run Rscript from inside ",
         "the ccadapt folder or its subdirectories.")
  }
  wd
}

CCA_ROOT <- cca_project_root()

# ---- Spatial reference ----------------------------------------------------
CCA_CRS   <- "EPSG:4326"
CCA_RES   <- 5 / 60               # 5 arc-minutes, exact rational -- see notes.
CCA_RES_LABEL <- "5arcmin"        # File-naming token; keep in sync with CCA_RES.

# ---- Countries ------------------------------------------------------------
# Priority in the investment doc: Nigeria and Kenya first (rapid assessment),
# Ethiopia after Activity-1 practice confirmation (defaulting to
# "Agroforestry + cover crops" per stakeholder clarification on the plan).
CCA_COUNTRIES <- c("NGA", "KEN", "ETH")

# ---- Paths ----------------------------------------------------------------
CCA_PATHS <- list(
  inputs_root      = file.path(CCA_ROOT, "data", "inputs"),
  staging_root     = file.path(CCA_ROOT, "data", "staging"),
  outputs_root     = file.path(CCA_ROOT, "data", "outputs"),
  naturebase       = file.path(CCA_ROOT, "data", "naturebase", "all_pathways"),
  text_root        = file.path(CCA_ROOT, "text"),
  # Staging subfolders -- pre-created here so downstream scripts can rely on them
  template_dir     = file.path(CCA_ROOT, "data", "staging", "template"),
  admin_dir        = file.path(CCA_ROOT, "data", "staging", "admin"),
  # Outputs subfolders
  step1_dir        = file.path(CCA_ROOT, "data", "outputs", "step1_mitigation"),
  step2_dir        = file.path(CCA_ROOT, "data", "outputs", "step2_production"),
  step3_dir        = file.path(CCA_ROOT, "data", "outputs", "step3_adaptation"),
  step4_dir        = file.path(CCA_ROOT, "data", "outputs", "step4_integration"),
  step5_dir        = file.path(CCA_ROOT, "data", "outputs", "step5_uncertainty"),
  activity3_dir    = file.path(CCA_ROOT, "data", "outputs", "activity3_neb")
)

# Create output directories idempotently -- safer than scattering dir.create()
# calls across scripts where they may run out of order.
invisible(lapply(CCA_PATHS, function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}))

# ---- COG write options ----------------------------------------------------
# Matches the atlas convention used in spatial-foundation and GeoCloudConvert:
# DEFLATE compression (lossless, widely supported by GDAL/web clients),
# 512-pixel internal tiles (good balance between HTTP range-request count
# and per-request payload for a Leaflet/MapLibre viewer).
CCA_COG_OPTS <- c(
  "TILED=YES",
  "BLOCKSIZE=512",
  "COMPRESS=DEFLATE",
  "BIGTIFF=IF_SAFER"
)

# ---- GeoParquet write options --------------------------------------------
# Snappy is the atlas default. WKB geometry encoding because arrow/ogr2ogr
# consumers on the web stack read it natively.
CCA_GEOPARQUET_OPTS <- c(
  "COMPRESSION=SNAPPY",
  "GEOMETRY_ENCODING=WKB"
)

# ---- Helpers --------------------------------------------------------------
cca_template_path <- function(iso3 = NULL) {
  # Country-specific templates live beside the regional one. If iso3 is NULL
  # returns the 3-country bounding-box template.
  stem <- if (is.null(iso3)) "cca_template_3country" else paste0("cca_template_", iso3)
  file.path(CCA_PATHS$template_dir, paste0(stem, "_", CCA_RES_LABEL, ".tif"))
}

cca_admin_path <- function(iso3, level, ext = "gpkg") {
  stopifnot(level %in% c("adm0", "adm1", "adm2"))
  file.path(CCA_PATHS$admin_dir,
            sprintf("gadm_%s_%s.%s", iso3, level, ext))
}
