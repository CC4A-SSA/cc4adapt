# Shared output locations for the af_formulations scripts.
#   AF_OUT_DIR  canonical results dir (default data/outputs/af_formulations, gitignored)
#   AF_RECORD_DIR  optional second copy for the project's _claude record; used
#                  only when the env var is set and its parent exists.
af_out_dir <- function() {
  d <- Sys.getenv("AF_OUT_DIR", "data/outputs/af_formulations")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
af_record_dir <- function() {
  d <- Sys.getenv("AF_RECORD_DIR", "")
  if (!nzchar(d)) return(NULL)
  d <- path.expand(d)
  if (!dir.exists(dirname(d))) return(NULL)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
af_write <- function(obj, name, writer = arrow::write_parquet) {
  writer(obj, file.path(af_out_dir(), name))
  rec <- af_record_dir()
  if (!is.null(rec)) writer(obj, file.path(rec, name))
  invisible(file.path(af_out_dir(), name))
}
