# ------------------------------------------------------------------------------
# code/common/adaptation_split.R
#
# How a practice's yield benefit is split between PRODUCTION (Step 2 / the
# production layer) and ADAPTATION (Step 3 avoided loss).
#
# Decision (project lead, Nov 2026; tracking issue CC4A-SSA/cc4adapt#1;
# https://cc4a-ssa.github.io/cc4adapt/af_formulations/ section 5): use the
# "same share lost" rule, called F3 in code/explore/af_formulations/.
#
# Notation (per crop c, hazard h, severity band s, pixel):
#   B   baseline (control) yield, t/ha
#   g   practice's relative yield gain from ERA, so A = B * (1 + g)
#   L   share of yield the control loses in band s, in [0, 1]
#   m   share of that loss the practice avoids. SIGNED: m < 0 means the
#       practice makes bad years worse (maladaptation).
#   P   annual probability of band s (bands are disjoint, sum_s P <= 1)
#
#   Control, stress year    Bs = B (1 - L)
#   Practice, stress year   As = A (1 - L (1 - m))
#   Stress-year gap         As - Bs
#
# F3 split of the stress-year gap:
#   production    = (A - B)(1 - L)      the gain, bitten by the same share L
#   avoided loss  = A L m               loss avoided on the PRACTICE's yield
#   production + avoided loss = As - Bs exactly (no double counting).
#   In a normal season production = A - B and avoided loss = 0.
#
# Expected annual values over bands:
#   E[L]          = sum_s P(s) L(s)
#   E[production] = (A - B)(1 - E[L])
#   E[avoided]    = sum_s P(s) A L(s) m(s)
#
# What this replaces (F5 in the exploration): production A - B every year and
# avoided loss B L m on the CONTROL's yield, with m clamped to [0, 1] and the
# total clamped at >= 0. That over-counted the total by L (A - B)(1 - m) per
# stress year and could not record maladaptation.
#
# Reference algebra: code/explore/af_formulations/R/formulations.R
# (af_stress_year, af_expected_annual). Re-implemented here on purpose:
# pipeline code must not source code/explore/.
#
# The arithmetic helpers (cca_f3_production, cca_f3_avoided_loss,
# cca_combine_expected_loss) use only +, -, * so they work on plain numerics
# and on terra SpatRasters alike. The checked wrappers (cca_f3_split,
# cca_expected_split, cca_expected_loss) are for numeric inputs.
#
# Pure base R (terra only touched if a SpatRaster is passed); sourceable
# without config.R so tests can load it alone.
# ------------------------------------------------------------------------------

# Sanity bound on |m|. A share of loss avoided above 1 would mean the practice
# yields MORE under stress than in a normal year; below -1 would mean it more
# than doubles the control's loss. Both are treated as data errors.
CCA_BUFFER_BOUND <- 1

# ---- Internal helpers ----------------------------------------------------

.cca_split_warn <- function(msg) {
  if (exists("cca_log", mode = "function")) cca_log(msg, level = "WARN")
  warning(msg, call. = FALSE)
}

.cca_recycle <- function(args) {
  lens <- vapply(args, length, integer(1))
  n <- max(lens)
  if (any(!lens %in% c(1L, n))) {
    stop("Arguments ", paste(names(args), collapse = ", "),
         " must have length 1 or ", n, " (got ", paste(lens, collapse = ", "), ")",
         call. = FALSE)
  }
  lapply(args, rep_len, length.out = n)
}

.cca_check_numeric <- function(args) {
  for (nm in names(args)) {
    x <- args[[nm]]
    if (!is.numeric(x)) stop(nm, " must be numeric", call. = FALSE)
    if (anyNA(x)) stop("NA in ", nm, call. = FALSE)
    if (any(!is.finite(x))) stop("Non-finite value in ", nm, call. = FALSE)
  }
  invisible(TRUE)
}

# ---- Buffering coefficient ----------------------------------------------

#' Signed buffering coefficient m with a sanity bound.
#'
#' Replaces the old [0, 1] clamp. Negative values are kept (maladaptation).
#' Values outside [-bound, bound] are clamped to the bound and a warning is
#' logged naming `label`. NA / NaN are returned as NA (callers decide what
#' missing evidence means; Step 3 uses 0 = no buffering).
#'
#' @param m numeric vector.
#' @param bound positive number, default CCA_BUFFER_BOUND (1).
#' @param label text used in the warning (e.g. "agroforestry/drought/severe").
#' @return numeric vector, same length as m.
cca_signed_buffer <- function(m, bound = CCA_BUFFER_BOUND, label = "buffer_frac") {
  if (!is.numeric(m)) stop("m must be numeric", call. = FALSE)
  m[is.nan(m)] <- NA_real_
  out <- !is.na(m) & abs(m) > bound
  if (any(out)) {
    .cca_split_warn(sprintf(
      "%s: %d value(s) outside [-%g, %g] (range %s); clamped to the bound",
      label, sum(out), bound, bound,
      paste(signif(range(m[out]), 4), collapse = " to ")))
    m[out] <- sign(m[out]) * bound
  }
  m
}

# ---- F3 arithmetic (numeric or SpatRaster) -------------------------------

#' Practice yield A = B (1 + g).
cca_practice_yield <- function(B, g) B * (1 + g)

#' F3 avoided loss in one stress band: A L m (practice's own yield).
cca_f3_avoided_loss <- function(A, L, m) A * L * m

#' F3 production: (A - B)(1 - L). With L = 0 this is the normal-season gain;
#' with L = E[L] it is the expected annual production.
cca_f3_production <- function(A, B, L) (A - B) * (1 - L)

#' Combine per-hazard expected loss shares into one share for the crop.
#'
#' method = "independent" (default): 1 - prod_h (1 - E[L_h]). Assumes hazards
#'   occur independently and that losses compound multiplicatively in a year
#'   hit by more than one hazard. This is an ASSUMPTION: co-occurrence of
#'   hazards is unresolved project-wide.
#' method = "additive": min(1, sum_h E[L_h]). Matches the way Step 3 sums
#'   avoided loss over hazards (no joint-year correction).
#'
#' @param el list of per-hazard E[L] (numeric vectors of equal length, or
#'   SpatRasters on the same grid). Each must lie in [0, 1].
#' @return same type as the elements of el.
cca_combine_expected_loss <- function(el, method = c("independent", "additive")) {
  method <- match.arg(method)
  if (!is.list(el) || !length(el)) stop("el must be a non-empty list", call. = FALSE)
  if (method == "independent") {
    surv <- Reduce(`*`, lapply(el, function(x) 1 - x))
    return(1 - surv)
  }
  s <- Reduce(`+`, el)
  if (inherits(s, "SpatRaster")) terra::clamp(s, upper = 1, values = TRUE) else pmin(s, 1)
}

# ---- Checked numeric wrappers -------------------------------------------

#' F3 split for one stress band (vectorised over rows).
#'
#' @param B baseline yield (>= 0).
#' @param g relative gain (> -1).
#' @param L control loss share in the band, [0, 1].
#' @param m signed share of loss avoided (any finite value; bound it upstream
#'   with cca_signed_buffer()).
#' All arguments recycled to a common length (each must be length 1 or n).
#' @return data.frame with B, g, L, m, A, Bs, As, gap (= As - Bs),
#'   production, avoided_loss, total (= production + avoided_loss = gap).
cca_f3_split <- function(B, g, L, m) {
  args <- .cca_recycle(list(B = B, g = g, L = L, m = m))
  .cca_check_numeric(args)
  B <- args$B; g <- args$g; L <- args$L; m <- args$m
  if (any(B < 0)) stop("B must be >= 0", call. = FALSE)
  if (any(g <= -1)) stop("g must be > -1", call. = FALSE)
  if (any(L < 0 | L > 1)) stop("L must lie in [0, 1]", call. = FALSE)

  A  <- cca_practice_yield(B, g)
  Bs <- B * (1 - L)
  As <- A * (1 - L * (1 - m))
  prod <- cca_f3_production(A, B, L)
  al   <- cca_f3_avoided_loss(A, L, m)
  data.frame(B = B, g = g, L = L, m = m, A = A, Bs = Bs, As = As,
             gap = As - Bs, production = prod, avoided_loss = al,
             total = prod + al)
}

# Coerce band inputs to an n x k matrix.
.cca_band_matrix <- function(x, nm, n = NULL) {
  if (!is.numeric(x)) stop(nm, " must be numeric", call. = FALSE)
  if (is.matrix(x)) return(x)
  matrix(x, nrow = if (is.null(n)) 1L else n, ncol = length(x), byrow = TRUE)
}

#' Expected annual loss share E[L] = sum_s P(s) L(s).
#'
#' @param p_band band probabilities: length-k vector (one unit) or n x k
#'   matrix. Must be >= 0 with row sums <= 1 (disjoint bands).
#' @param L loss shares, same shape as p_band (or a length-k vector recycled
#'   over rows), in [0, 1].
#' @return numeric vector of length n.
cca_expected_loss <- function(p_band, L) {
  P <- .cca_band_matrix(p_band, "p_band")
  Lm <- .cca_band_matrix(L, "L", n = nrow(P))
  if (!identical(dim(P), dim(Lm))) {
    stop("p_band and L must have the same number of bands (and rows)", call. = FALSE)
  }
  .cca_check_numeric(list(p_band = as.vector(P), L = as.vector(Lm)))
  if (any(P < 0)) stop("Negative p_band", call. = FALSE)
  if (any(rowSums(P) > 1 + 1e-9)) stop("Band probabilities sum above 1", call. = FALSE)
  if (any(Lm < 0 | Lm > 1)) stop("L must lie in [0, 1]", call. = FALSE)
  as.vector(rowSums(P * Lm))
}

#' Expected annual F3 production and avoided loss over severity bands.
#'
#' @param B,g baseline yield and relative gain, length 1 or n.
#' @param p_band band probabilities, length-k vector or n x k matrix.
#' @param L,m band loss shares and signed buffering, same shape as p_band or
#'   length-k vectors recycled over rows.
#' @return data.frame with B, g, A, p_stress (= sum_s P), expected_L,
#'   production (= (A - B)(1 - E[L])), avoided_loss (= sum_s P A L m), total.
cca_expected_split <- function(B, g, p_band, L, m) {
  P <- .cca_band_matrix(p_band, "p_band")
  n <- max(nrow(P), length(B), length(g))
  if (nrow(P) == 1L && n > 1L) P <- P[rep(1L, n), , drop = FALSE]
  if (nrow(P) != n) stop("p_band has ", nrow(P), " rows but B/g imply ", n, call. = FALSE)
  Lm <- .cca_band_matrix(L, "L", n = n)
  Mm <- .cca_band_matrix(m, "m", n = n)
  if (!identical(dim(Lm), dim(P)) || !identical(dim(Mm), dim(P))) {
    stop("p_band, L and m must have the same number of bands (and rows)", call. = FALSE)
  }
  bg <- .cca_recycle(list(B = B, g = g))
  if (length(bg$B) != n) bg <- lapply(bg, rep_len, length.out = n)
  .cca_check_numeric(list(B = bg$B, g = bg$g, m = as.vector(Mm)))
  if (any(bg$B < 0)) stop("B must be >= 0", call. = FALSE)
  if (any(bg$g <= -1)) stop("g must be > -1", call. = FALSE)

  EL <- cca_expected_loss(P, Lm)          # checks P and L
  A  <- cca_practice_yield(bg$B, bg$g)
  prod <- cca_f3_production(A, bg$B, EL)
  al   <- A * rowSums(P * Lm * Mm)
  data.frame(B = bg$B, g = bg$g, A = A, p_stress = rowSums(P),
             expected_L = EL, production = prod, avoided_loss = al,
             total = prod + al)
}

# ---- Pipeline path helper -------------------------------------------------

#' Production layer for a practice x country: the stress-adjusted F3 layer
#' written by Step 3 when it exists, else Step 2's dvalue (normal-year gain).
#' Requires config.R (CCA_PATHS, CCA_RES_LABEL). Returns a path (which may
#' not exist if neither layer was produced).
cca_production_path <- function(practice, iso3, quiet = FALSE) {
  adj <- file.path(CCA_PATHS$step3_dir,
                   sprintf("production_adj_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
  raw <- file.path(CCA_PATHS$step2_dir,
                   sprintf("dvalue_%s_%s_%s.tif", practice, iso3, CCA_RES_LABEL))
  if (file.exists(adj)) return(adj)
  if (!quiet && file.exists(raw) && exists("cca_log", mode = "function")) {
    cca_log(sprintf(paste0("No stress-adjusted production layer for %s/%s; using Step 2 ",
                           "dvalue (normal-year gain, not reduced in stress years)"),
                    practice, iso3), level = "WARN")
  }
  raw
}
