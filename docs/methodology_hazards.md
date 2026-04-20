# Hazard frequency rasters — sourcing + methodology

**Status:** authoritative reference for every script that produces or consumes a "hazard frequency" raster in the CCA pipeline (currently `code/activity2/step3_adaptation.R`).
**Owner:** ccadapt maintainers.
**Last revised:** 2026-04-18.

This document explains, step by step, how the per-pixel × per-(crop, hazard,
severity) hazard frequency rasters consumed by Step 3 are derived. As with
[`methodology_biochar.md`](methodology_biochar.md), this is the source of
truth referenced by the headers of every modified `.R` and `.yaml` file in
the change; if anything below changes, those headers, MANIFEST.md, and the
README must change with it.

---

## 1. Why we cannot just point at the atlas hazards_prototype outputs

The CCA pipeline's `code/activity2/step3_adaptation.R` was originally written
to look for hazard frequency rasters under
`atlas/hazards_prototype/<crop>/<hazard>/<severity>.tif`. Investigation
(2026-04-18) shows three reasons that does not work as a standalone
arrangement:

1. **The atlas/hazards_prototype/ folder ships only source code, not output
   rasters.** The prototype's `R/2_calculate_haz_freq.R` writes COGs into a
   sibling `Data/hazard_timeseries_risk/` directory which is populated only
   when the prototype is *run* — it has not been run inside this checkout.
2. **The canonical published outputs live in a private S3 bucket.** The
   atlas `metadata/data.json` documents that hazard rasters are published at
   `s3://digital-atlas/domain=hazard_impact/source=atlas_cmip6/region=ssa/time={time}/processing=hazard-risk/variable={variable}/model={model}/severity={severity}/`.
   We confirmed empirically (HTTP `403 Forbidden` on every probe) that this
   prefix requires authenticated S3 access. Public consumers — including
   anyone replicating CCA on a fresh laptop — cannot pull from it.
3. **Even if S3 were public, the published rasters multiply hazard frequency
   by exposure (`variable=vop_intld15` etc.).** Step 3 needs *raw* hazard
   frequency P(h, s, c, p), not frequency × value-of-production. The atlas
   does not currently publish a separately-addressable raw-frequency layer.

Conclusion: the only defensible per-pixel hazard frequency for CCA's
standalone-replicable pipeline must be **synthesised locally** from
publicly-available climatological data and crop-specific thresholds. As
with biochar, every term then has an auditable citation.

A second, lower-priority pathway is preserved: if the atlas
`hazards_prototype` outputs *are* available locally (e.g. when CCA is run
inside the EiA server with the prototype caches mounted), the synthesiser
short-circuits and uses them — see `§5 Two-tier resolution` below.

---

## 2. The first-principles formula

For each grid cell *p*, crop *c*, hazard *h*, severity *s*:

```
P(h, s, c, p)  [dimensionless 0..1]
  = empirical_freq(h, s, c, p)                    if a per-year time series exists
  = parametric_exceedance(h, s, c, p)             otherwise (default)
```

The **parametric exceedance** path treats inter-annual climate variability
as approximately normal around the WorldClim 2.1 climatological mean with
a hazard-specific coefficient of variation drawn from the SSA literature.
Per pixel, hazard, and severity:

### Heat (growing-season Tmax > crop optimum)

```
mu_T(p)   = mean monthly Tmax over the crop's growing season         # WorldClim 2.1
sigma_T(p) = mu_T(p) * CV_T                                          # CV_T = 0.05 (literature)
T_thr(s, c) = T_opt_max(c) + delta_T(s)                               # Ecocrop + per-severity offset

P_heat(s, c, p) = pnorm(T_thr(s, c) + delta_band, mu_T, sigma_T)
                - pnorm(T_thr(s, c),               mu_T, sigma_T)     # for moderate, severe
P_heat(extreme, c, p) = 1 - pnorm(T_thr(extreme, c), mu_T, sigma_T)   # for extreme
```

with severity offsets (above each crop's `Temp_Opt_Max`):
- moderate: +1°C  ≤ T < +3°C
- severe:   +3°C  ≤ T < +5°C
- extreme:  T ≥ +5°C

### Drought (growing-season precipitation < crop minimum)

```
mu_R(p)    = total monthly precipitation summed over growing season   # WorldClim 2.1
sigma_R(p) = mu_R(p) * CV_R                                           # CV_R = 0.25 (SSA literature)
R_thr(s, c) = R_opt_min(c) * f_dry(s)                                  # Ecocrop * per-severity multiplier

P_drought(moderate, c, p) = pnorm(0.7 * R_opt_min, mu_R, sigma_R)
                          - pnorm(0.5 * R_opt_min, mu_R, sigma_R)
P_drought(severe,   c, p) = pnorm(0.5 * R_opt_min, mu_R, sigma_R)
                          - pnorm(0.3 * R_opt_min, mu_R, sigma_R)
P_drought(extreme,  c, p) = pnorm(0.3 * R_opt_min, mu_R, sigma_R)
```

### Waterlogging (growing-season precipitation > crop maximum)

```
P_wet(moderate, c, p) = pnorm(1.5 * R_opt_max, mu_R, sigma_R)
                       - pnorm(1.2 * R_opt_max, mu_R, sigma_R)
P_wet(severe,   c, p) = pnorm(2.0 * R_opt_max, mu_R, sigma_R)
                       - pnorm(1.5 * R_opt_max, mu_R, sigma_R)
P_wet(extreme,  c, p) = 1 - pnorm(2.0 * R_opt_max, mu_R, sigma_R)
```

All severity bands are non-overlapping by construction, so the three
severity rasters for a given (hazard, crop) sum to ≤ 1 at every pixel
(equality only when there is no "no-hazard" outcome in the band — usually
< 1 because the central probability mass falls outside any of the three
bands).

---

## 3. Inputs the synthesiser reads

| Input | Path | Origin | Notes |
|---|---|---|---|
| WorldClim 2.1 monthly Tmax | `data/inputs/soilgrids/climate/wc2.1_10m/wc2.1_10m_tmax_<MM>.tif` | `geodata::worldclim_global("tmax", res=10)` (already cached by the suitability builders) | 12 monthly climatological means for 1970–2000 |
| WorldClim 2.1 monthly precipitation | `data/inputs/soilgrids/climate/wc2.1_10m/wc2.1_10m_prec_<MM>.tif` | same | 12 monthly climatological totals |
| Ecocrop thresholds | `code/inputs/data/ecocrop_subset_cca.csv` | imported once from `atlas/hazards_prototype/metadata/ecocrop.csv` (Zea mays, Oryza sativa, Sorghum bicolor) | per-crop optimum + absolute T and R |
| Hazard / severity parameters | `code/inputs/data/hazards_params.yaml` | hand-curated, citations in file | growing season per crop, CVs, severity offsets |
| CCA template | `data/staging/template/cca_template_3country_5arcmin.tif` | `build_cca_template.R` | output extent + resolution |

All inputs are local. No external HTTP download is required by the
synthesiser. `geodata::worldclim_global()` was already invoked by the
suitability builders, so the WorldClim cache is on disk after the
suitability stage.

---

## 4. Output

For each crop ∈ {`maiz`, `rice`, `sorg`}, hazard ∈ {`drought`, `heat`,
`waterlogging`}, severity ∈ {`moderate`, `severe`, `extreme`}, the
synthesiser writes one COG plus one `.method.json` sidecar:

```
data/staging/hazards/<crop>/<hazard>/<severity>_5arcmin.tif
data/staging/hazards/<crop>/<hazard>/<severity>_5arcmin.tif.method.json
```

For 3 × 3 × 3 = **27 hazard frequency COGs** plus 27 sidecars. The
filenames embed the resolution token (`_5arcmin`) so that future ensemble
runs at finer resolution land at distinct paths.

`step3_adaptation.R` reads every existing combination via a glob; missing
files are skipped silently per the existing graceful-skip contract.

---

## 5. Two-tier resolution (atlas-cache → synthesis)

The synthesiser tries inputs in this order; the *first* one that yields a
usable raster wins (per pixel, per crop × hazard × severity):

1. **Atlas cache pass-through.** If
   `atlas/hazards_prototype/Data/hazard_risk_n/annual/<crop>-<severity>.tif`
   (or any of the prototype's documented output trees) is mounted and
   readable, the synthesiser reprojects it onto the CCA template and
   tags the output `.method.json` with `source="atlas_hazards_prototype"`.
2. **Local synthesis.** Otherwise, compute the parametric exceedance
   probabilities of §2 from WorldClim + Ecocrop. Tag with
   `source="worldclim_ecocrop_synthesis"`.

Tier 1 currently never fires in standalone CCA (the prototype's `Data/`
directory is not in the repo and S3 is private). It is wired so that an
EiA-server run gets the more authoritative inputs without a code change.

---

## 6. Cross-check rules

After the synthesiser runs:

- **Sum of severities ≤ 1.** For each (crop, hazard, pixel),
  `P(moderate) + P(severe) + P(extreme) ≤ 1` (severity bands are
  disjoint). Logged at INFO; deviations > 1.001 logged ERROR.
- **Per-country plausibility.** For each (crop, hazard), the country-mean
  drought frequency must be in [0, 0.95] (plausibility band; values near
  1 indicate a crop entirely outside its rainfall envelope, which is
  physically possible — e.g. rice in arid northern Kenya — and so this
  is WARN, not ERROR).
- **Seasonality sanity check.** Growing-season-summed precipitation
  must be ≤ annual precipitation. Logged at INFO with the per-country
  ratio.

---

## 7. Limitations and the upgrade path

This is a **climatological baseline approximation**, not a CMIP6 ensemble.
Specifically:

- **Inter-annual variability is parametric, not empirical.** A literature
  CV of 0.05 (temperature) and 0.25 (precipitation) is applied uniformly.
  Real CVs vary spatially — Sahel rainfall CVs reach 0.40+, Ethiopian
  highland rainfall CVs are nearer 0.15. A future improvement is to
  derive per-pixel CV from CHIRPS 1981–2024 monthly time series (~3 GB
  download) or from CRU TS 4.07 (~150 MB).
- **No climate-change projections.** The frequencies are for the WorldClim
  1970–2000 baseline. SSP-specific projections require pulling the atlas
  CMIP6 ensemble from S3 (currently auth-gated) or running the
  hazards_prototype pipeline end-to-end.
- **Growing season is fixed per crop, not per-pixel.** A future improvement
  is to multiply by the GGCMI Phase 3 crop calendar (Zenodo
  10.5281/zenodo.5062513) to get pixel-specific growing-season windows.
- **Hazards are independent.** Step 3 already sums hazards as if
  independent; correlated-hazard formulations (drought ↔ heat
  co-occurrence) are flagged in [README.md "What needs checking" §6](../README.md).

The output `.method.json` sidecars carry a `confidence_tier = "Low"` flag
until Activity-1 stakeholder review confirms the parameter choices — the
same pattern used by `build_buffering_coefficients.R` for the ERA-derived
table.

---

## 8. Citations

| Element | Source |
|---|---|
| WorldClim 2.1 monthly climatology | Fick & Hijmans 2017 *Int. J. Climatol.* (DOI 10.1002/joc.5086); accessed via `geodata::worldclim_global()` |
| Ecocrop crop optima | FAO Ecocrop database, snapshot from `atlas/hazards_prototype/metadata/ecocrop.csv` (Zea mays, Oryza sativa, Sorghum bicolor) |
| SSA precipitation CV | Nicholson 2018 *Climate of Africa*, rainfall coefficient of variation tables for 1981–2010 — central value 0.20–0.30 |
| SSA monthly Tmax CV | New, Hulme, Jones 2000 *J. Clim.* — coefficient of variation of monthly Tmax climatology 0.04–0.07 in tropical Africa |
| Severity offset convention | Atlas `metadata/haz_classes.csv` (moderate/severe/extreme nomenclature), with band widths chosen to match the equivalent CMIP6-derived bins published by the atlas |
| Hazard naming (drought / heat / waterlogging) | Atlas `metadata/haz_metadata.csv` `type` field aggregations: dry → drought, heat → heat, wet → waterlogging |
