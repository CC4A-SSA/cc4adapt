# ------------------------------------------------------------------------------
# code/inputs/fetch_era.R
#
# Stages the Evidence for Resilient Agriculture (ERA) derived meta-analysis
# table used by Step 2 (practice yield effects) and build_buffering_coefficients.
#
# Sources tried in order (stop at first success):
#   (1) A copy already present in the CCA data/inputs/era/ folder.
#   (2) A sibling atlas checkout at ../solutions/data/ERA_Derived.rda
#       (the main atlas ships ERA inside the `solutions` repo).
#   (3) The public ERA Zenodo deposit (DOI resolved at runtime).
#
# Why this exists:
#   For a standalone CCA checkout on a fresh machine there is no sibling
#   `atlas/solutions/` folder. Having an explicit fetcher means the pipeline
#   can run end-to-end without requiring the user to clone the whole atlas
#   monorepo -- it just needs either (a) an ERA_Derived.rda drop in
#   data/inputs/era/ or (b) a Zenodo record that `jsonlite` can resolve.
#
# The script does NOT invent values if ERA is unavailable. Downstream
# consumers (Step 2, build_buffering_coefficients) already fall back to
# literature defaults when the table is missing, with clear WARN logs. This
# fetcher's job is to land the RDA if it can and stay silent otherwise.
# ------------------------------------------------------------------------------

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
source(file.path(common_dir, "geo_utils.R"))
source(file.path(common_dir, "download.R"))

cca_banner("Fetch ERA derived meta-analysis table")

era_dir <- file.path(CCA_PATHS$inputs_root, "era")
dir.create(era_dir, recursive = TRUE, showWarnings = FALSE)
era_dest <- file.path(era_dir, "ERA_Derived.rda")

# ---- (1) Already present ? ---------------------------------------------
if (cca_skip_if_exists(era_dest, tag = "ERA derived")) quit(status = 0)

# ---- (2) Sibling atlas/solutions checkout ------------------------------
# Works on the atlas-monorepo dev setup; harmless elsewhere.
sibling <- normalizePath(file.path(CCA_ROOT, "..", "solutions", "data", "ERA_Derived.rda"),
                         mustWork = FALSE)
if (file.exists(sibling)) {
  cca_log("Found ERA in sibling atlas checkout: ", sibling)
  file.copy(sibling, era_dest, overwrite = FALSE)
  cca_log("Copied to: ", era_dest)
  # Minimal provenance sidecar.
  writeLines(c("{",
               sprintf('  "source": "sibling atlas/solutions checkout at %s",', sibling),
               sprintf('  "fetched_at": "%s",', format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
               sprintf('  "bytes": %d', file.size(era_dest)),
               "}"),
             paste0(era_dest, ".source.json"))
  quit(status = 0)
}

# ---- (3) Zenodo lookup -------------------------------------------------
# ERA is intermittently published on Zenodo under the EiA (Excellence in
# Agronomy) initiative. We search for the most recent candidate; any
# matching record still requires MANUAL confirmation of identity before
# auto-download, per the project-lead rule on silent data ingest.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  cca_log("jsonlite not installed; cannot resolve Zenodo. install.packages('jsonlite').",
          level = "ERROR")
  quit(status = 0)
}

zenodo_search <- function(query, size = 5L) {
  api <- sprintf('https://zenodo.org/api/records?q=%s&size=%d&sort=mostrecent',
                 utils::URLencode(query, reserved = TRUE), size)
  tmp <- tempfile(fileext = ".json")
  ok <- tryCatch({ cca_download(api, tmp, force = TRUE); TRUE },
                 error = function(e) FALSE)
  if (!ok) return(list())
  j <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
  j$hits$hits %||% list()
}

hits <- zenodo_search('"Evidence for Resilient Agriculture" ERA_Derived', 5L)
if (!length(hits)) {
  cca_log("ERA not found via sibling or Zenodo. Manual step: drop ERA_Derived.rda into ",
          era_dir, ". Step 2 + buffering coefficients fall back to literature defaults ",
          "without ERA, which is flagged as a WARN.", level = "WARN")
  quit(status = 0)
}

cca_log(sprintf("Zenodo returned %d candidate ERA record(s); manual review needed:",
                length(hits)))
for (h in hits) {
  title <- h$metadata$title %||% "(untitled)"
  doi   <- h$doi %||% h$metadata$doi %||% ""
  cca_log(sprintf("  - '%s' (DOI %s)", title, doi))
}
cca_log("Auto-download not performed: confirm the record identity, then drop the ",
        "file into ", era_dir, " (or set ERA_ZENODO_DOI and re-run).",
        level = "WARN")

cca_banner("Done: ERA fetch (manual confirmation required)")
