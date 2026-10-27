# ------------------------------------------------------------------------------
# code/explore/af_formulations/02_check_inputs.R
#
# EXPLORATORY. Sanity checks on the step 5.2 inputs before any P(h, s) is
# derived. Each check prints PASS / WARN / FAIL with the measured value.
# Identities follow hazards_prototype HANDOVER_2026-10-12 S1.4 and S1.7.
#
# Run from the repo root:  Rscript code/explore/af_formulations/02_check_inputs.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(terra)
})

IN   <- "data/staging/af_formulations"
SETS <- "NDWS+NTx35+NDWL0"
CROP <- "maize"
CLASSES <- c("dry", "heat", "wet", "dry+heat", "dry+wet", "heat+wet", "dry+heat+wet")

verdict <- function(label, ok, value, warn = FALSE) {
  tag <- if (ok) "PASS" else if (warn) "WARN" else "FAIL"
  cat(sprintf("[%s] %-58s %s\n", tag, label, value))
}

tiers <- rbindlist(lapply(c("moderate", "severe", "extreme"), function(s) {
  as.data.table(read_parquet(file.path(IN, sprintf("atlas_prod_t_jagermeyr_%s.parquet", s))))
}))
tiers[, level := fifelse(!is.na(admin2_name), "adm2", fifelse(!is.na(admin1_name), "adm1", "adm0"))]

cat("\n== Tier contents\n")
print(tiers[, .N, by = .(severity)])
print(dcast(tiers[level == "adm2" & scenario == "historic", .N, by = .(crop, hazard_vars)],
            crop ~ hazard_vars, value.var = "N"))
cat("scenarios:", paste(sort(unique(paste(tiers$scenario, tiers$timeframe))), collapse = " | "), "\n")
cat("exposure:", unique(tiers$exposure_var), unique(tiers$exposure_unit), "\n\n")

h <- tiers[crop == CROP & hazard_vars == SETS & scenario == "historic"]
w <- dcast(h, iso3 + level + gaul0_code + gaul1_code + gaul2_code + severity ~ hazard,
           value.var = "value")

cat("== Identities (maize, ", SETS, ", historic)\n", sep = "")
a2  <- w[level == "adm2"]
nan <- a2[, .(units = uniqueN(gaul2_code), nan = uniqueN(gaul2_code[is.nan(any)])), by = iso3]
for (i in seq_len(nrow(nan))) {
  verdict(sprintf("%s adm2 units NaN in every class", nan$iso3[i]),
          nan$nan[i] / nan$units[i] < 0.4,
          sprintf("%d of %d", nan$nan[i], nan$units[i]), warn = TRUE)
}

ok <- a2[!is.nan(any) & (any + none) > 0]
ok[, cls_sum := rowSums(.SD), .SDcols = CLASSES]
rel <- ok[any > 100, abs(cls_sum - any) / any]
verdict("7 classes sum to `any` (adm2, any > 100 t): share within 1 %",
        mean(rel <= 0.01) > 0.95, sprintf("%.1f %% of %d", 100 * mean(rel <= 0.01), length(rel)))

tot <- ok[, .(spread = diff(range(any + none)) / mean(any + none)), by = gaul2_code]
verdict("any + none identical across severities (rel. spread < 1 %)",
        mean(tot$spread < 0.01) > 0.95, sprintf("%.1f %% of units", 100 * mean(tot$spread < 0.01)))

# Single-hazard marginals must be monotone across cumulative tiers.
ok[, p_dry  := (dry + `dry+heat` + `dry+wet` + `dry+heat+wet`) / (any + none)]
ok[, p_heat := (heat + `dry+heat` + `heat+wet` + `dry+heat+wet`) / (any + none)]
ok[, p_wet  := (wet + `dry+wet` + `heat+wet` + `dry+heat+wet`) / (any + none)]
mono <- dcast(melt(ok[, .(gaul2_code, severity, p_dry, p_heat, p_wet)],
                   id.vars = c("gaul2_code", "severity")),
              gaul2_code + variable ~ severity, value.var = "value")
mono <- mono[complete.cases(mono)]
viol <- mono[, (severe - moderate > 0.01) | (extreme - severe > 0.01)]
verdict("marginals monotone mod >= sev >= ext (tol 0.01)", mean(viol) < 0.01,
        sprintf("%d violations of %d", sum(viol), nrow(mono)))

cat("\n== Denominator\n")
den <- as.data.table(read_parquet(file.path(IN, "atlas_denominator_res-25_t_ha.parquet")))
print(den[, .N, by = .(unit, tech, exposure)][order(unit)])
den2 <- den[crop == CROP & unit == "t" & !is.na(admin2_name) & tech == "all",
            .(den_t = sum(value)), by = gaul2_code]
cmp <- merge(ok[severity == "severe", .(gaul2_code, tier_t = any + none)], den2, by = "gaul2_code")
cmp <- cmp[den_t > 100]
r <- cmp[, abs(tier_t - den_t) / den_t]
verdict("tier any+none vs denominator t (adm2, maize): within 2 %",
        mean(r <= 0.02) > 0.9, sprintf("%.1f %% of %d units; median rel diff %.4f",
                                       100 * mean(r <= 0.02), length(r), median(r)), warn = TRUE)

cat("\n== GAUL admin-2\n")
g <- as.data.table(read_parquet(file.path(IN, "gaul24_adm2_ken_eth_nga.parquet")))
cat("columns:", paste(names(g), collapse = ", "), "\n")
code_col <- grep("gaul2|adm2_code|ADM2_CODE", names(g), value = TRUE, ignore.case = TRUE)[1]
iso_col  <- grep("^iso3$", names(g), value = TRUE, ignore.case = TRUE)[1]
print(g[, .N, by = c(iso_col)])
miss <- setdiff(unique(a2$gaul2_code), g[[code_col]])
verdict("every tier adm2 code present in GAUL", length(miss) == 0,
        sprintf("%d missing (join column %s)", length(miss), code_col))

cat("\n== Rasters\n")
aez <- rast(list.files(file.path(IN, "aez"), pattern = "^AEZ16.*\\.tif$", full.names = TRUE)[1])
fs2 <- rast(file.path(IN, "farming_systems", "FS2.tif"))
cat("AEZ16:", paste(dim(aez), collapse = "x"), "res", paste(signif(res(aez), 4), collapse = ","), "\n")
print(head(cats(aez)[[1]], 20))
cat("FS2 res", paste(signif(res(fs2), 4), collapse = ","), "\n")
