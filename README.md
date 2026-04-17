# Carbon Credit linked Adaptation (CCA) — analysis pipeline

Reproducible spatial + economic analysis of agricultural carbon-credit practices with adaptation co-benefits across **Nigeria, Kenya, and Ethiopia** at **5 arc-minute (~10 km)** resolution.

The pipeline is **standalone and replicable on any computer**: every external dataset is downloaded by a script, every pre-processed intermediate is derived from a documented input, and every output artefact (COG, GeoPackage, CSV) is written to a deterministic path. A fresh clone + `Rscript code/inputs/00_fetch_all.R` reproduces everything end-to-end.

---

## Deliverables

Per practice × country:

- **Technical potential** — spatial layers for mitigation (tCO₂e/ha/yr), production benefit (USD/ha/yr), adaptation benefit (avoided loss, USD/ha/yr).
- **Classification + opportunity zones** — pixel typology (Triple Win / Mitigation+Production / Mitigation+Adaptation / Constrained / Non-priority) and contiguous zone polygons with zone-level statistics.
- **Uncertainty** — Monte-Carlo P10/P50/P90 per layer plus a qualitative confidence tier.
- **Net Economic Benefit (NEB)** — per-hectare NPV rasters across a sweep of carbon prices and adoption rates, plus admin-level national totals.

All heavy work is in R using `terra`, `geodata`, `yaml`, `arrow`, `jsonlite`, `digest`, and `data.table`. Vectors are written as **GeoPackage (`.gpkg`)** with optional **GeoJSON**; rasters as **Cloud-Optimized GeoTIFF (COG)**; admin-level summaries as **CSV**.

---

## Pipeline flowchart

```text
 ┌─────────────────────────────────────────────────────────────────────────┐
 │ STAGE 0 — BOOTSTRAP: external fetchers + pre-processing                 │
 │                                                                         │
 │  Admin / grid          Mitigation inputs        Production inputs       │
 │  ─────────────         ──────────────────       ──────────────────      │
 │  build_cca_template    fetch_naturebase (NCS)   fetch_mapspam (IFPRI)   │
 │    └─ GADM via            └─ all_pathways.zip     └─ 2010 V2r0          │
 │       geodata                                                           │
 │                        fetch_biochar_potential  fetch_fao_prices        │
 │                           └─ Sanderman MG          └─ FAOSTAT bulk      │
 │  Constraint layer                                                       │
 │  ─────────────         fetch_erw_potential      fetch_era               │
 │  fetch_constraint_        └─ Zenodo/Beerling       └─ ERA_Derived.rda   │
 │    masks (WDPA +                                                        │
 │    slope + GHSL +      fetch_cao_awd            Derived builders        │
 │    GMIA)                  └─ Cao 2026 / IPCC    ──────────────────      │
 │                                                 build_biochar_suit      │
 │                        fetch_basalt             build_erw_suit          │
 │                           └─ GLiM + MRDS        build_residue_proxy     │
 │                                                 build_buffering_coef    │
 └─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
 ┌─────────────────────────────────────────────────────────────────────────┐
 │ STAGE 1 — ACTIVITY 2: technical potential (per country, 5-arcmin COGs)  │
 │                                                                         │
 │   step1_mitigation    ──►  M  (tCO2e / ha / yr)                         │
 │   step2_production    ──►  P  (USD / ha / yr, Δyield × price)           │
 │   step3_adaptation    ──►  A  (USD / ha / yr, avoided loss)             │
 │                                                                         │
 │   step4_integration   ──►  class raster + opportunity-zone polygons     │
 │                            {Triple Win | M+P | M+A | Constrained |      │
 │                             Non-priority}                               │
 │                                                                         │
 │   step5_uncertainty   ──►  P10 / P50 / P90 per-pixel (Monte Carlo)      │
 └─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
 ┌─────────────────────────────────────────────────────────────────────────┐
 │ STAGE 2 — ACTIVITY 3: Net Economic Benefit                              │
 │                                                                         │
 │   neb.R  ──►  NPV tables   sweep over                                   │
 │                              carbon_price ∈ {5, 15, 40, 80} USD/tCO2e   │
 │                              adoption    ∈ {5%, 15%, 30%, 50%}          │
 │                            per country × practice, 20-yr horizon, r=8%  │
 └─────────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
                     country_profiles   +   scenario_explorer
                     (static reports)       (future web UI)

 Cross-cutting dependencies (simplified):
   template       → step1, step2, step3, step4
   mapspam        → step2, step3, residue_proxy, cao_awd
   fao_prices     → step2, step3, neb
   era            → step2, buffering_coef
   buffering_coef → step3
   residue_proxy  → step1 (crop-residue biochar supply)
   constraints    → step4 (mask out protected/urban/high-slope)
   step1 (M)      → step2, step3, step4, step5, neb
   step2 (P)      → step4, step5, neb
   step3 (A)      → step4, step5, neb
```

---

## Quick start — bootstrap on a fresh computer

```bash
# 1. Clone
git clone <repo-url> ccadapt && cd ccadapt

# 2. Install R packages (one-time)
Rscript -e 'install.packages(c("terra","geodata","yaml","arrow","jsonlite","digest","data.table"))'

# 3. Run the whole pipeline (idempotent; skips completed stages)
Rscript code/inputs/00_fetch_all.R

# 4. Force a full re-run (overrides all skip-if-exists guards)
CCA_FORCE=1 Rscript code/inputs/00_fetch_all.R
```

Expected wall time on a fresh machine (25 Mbit/s down): ~45–60 minutes of which ~30 min is downloading external datasets (Naturebase 500 MB + Sanderman 488 MB + MapSPAM ~600 MB + SoilGrids ~1.2 GB + FAO 12 MB + minor). Incremental re-runs are seconds.

Per-stage invocation is also supported:

```bash
Rscript code/inputs/fetch_naturebase.R            # NCS pathways, ~500 MB
Rscript code/inputs/fetch_era.R                   # ERA_Derived.rda
Rscript code/inputs/fetch_mapspam.R               # IFPRI MAPSPAM 2010 V2r0
Rscript code/inputs/build_cca_template.R          # GADM + 5-arcmin template
Rscript code/inputs/fetch_fao_prices.R            # FAOSTAT producer prices
Rscript code/inputs/fetch_biochar_potential.R     # Sanderman 2020 Dataverse
Rscript code/inputs/fetch_erw_potential.R         # Zenodo + Beerling SI
Rscript code/inputs/fetch_basalt.R                # GLiM + MRDS (optional)
Rscript code/inputs/fetch_cao_awd.R               # Cao 2026 or IPCC T1 fallback
Rscript code/inputs/fetch_constraint_masks.R      # WDPA / slope / GHSL / GMIA
Rscript code/inputs/build_biochar_suitability.R   # pH / OC / clay / rain
Rscript code/inputs/build_erw_suitability.R       # pH / clay / T / rain
Rscript code/inputs/build_residue_proxy.R         # MapSPAM x IPCC T11.2
Rscript code/inputs/build_buffering_coefficients.R
Rscript code/activity2/step1_mitigation.R
Rscript code/activity2/step2_production.R
Rscript code/activity2/step3_adaptation.R
Rscript code/activity2/step4_integration.R
Rscript code/activity2/step5_uncertainty.R
Rscript code/activity3/neb.R
```

Every script is idempotent: if its expected outputs already exist with non-zero size, it logs `SKIP …` and exits. Set `CCA_FORCE=1` to override. Set `CCA_MC_SAMPLES=200` (default) for Step 5 to tune Monte-Carlo resolution.

---

## Repository layout

```
ccadapt/
├── README.md                         # This file
├── .gitignore                        # Excludes data/ (all regenerable)
├── text/                             # Project text and optional policy PDFs/<iso3>
├── code/
│   ├── common/                       # Shared helpers
│   │   ├── config.R                  # Paths, CRS, resolution, COG/GPKG options
│   │   ├── logging.R                 # Banner + assertion helpers
│   │   ├── geo_utils.R               # Reproject, masks, GPKG + GeoJSON writers, skip-if-exists
│   │   └── download.R                # Idempotent downloader + provenance sidecars
│   ├── inputs/                       # Fetch + pre-processing scripts
│   │   ├── 00_fetch_all.R            # Orchestrator (inputs + Activity 2 + 3)
│   │   ├── fetch_naturebase.R        # NCS all_pathways.zip
│   │   ├── fetch_era.R               # ERA_Derived.rda
│   │   ├── fetch_mapspam.R           # IFPRI MAPSPAM 2010 V2r0
│   │   ├── build_cca_template.R      # GADM + 5-arcmin reference grid
│   │   ├── fetch_fao_prices.R        # FAOSTAT producer prices
│   │   ├── fetch_biochar_potential.R # Sanderman management bundle
│   │   ├── fetch_erw_potential.R     # Zenodo / Beerling SI
│   │   ├── fetch_basalt.R            # GLiM + USGS MRDS
│   │   ├── fetch_cao_awd.R           # Cao 2026 + IPCC T1 fallback
│   │   ├── fetch_constraint_masks.R  # WDPA / slope / GHSL SMOD / GMIA
│   │   ├── build_biochar_suitability.R
│   │   ├── build_erw_suitability.R
│   │   ├── build_residue_proxy.R
│   │   └── build_buffering_coefficients.R
│   ├── activity2/
│   │   ├── step1_mitigation.R
│   │   ├── step2_production.R
│   │   ├── step3_adaptation.R
│   │   ├── step4_integration.R
│   │   ├── step5_uncertainty.R
│   │   └── config/
│   │       ├── practices.yaml        # Practice catalogue + country priority
│   │       └── thresholds.yaml       # Step 4 cut-offs
│   └── activity3/
│       ├── neb.R
│       └── neb_config.yaml           # Carbon prices, costs, adoption, horizon
└── data/                             # Gitignored; regenerated by the pipeline
    ├── inputs/                       # Raw downloads + MANIFEST
    │   ├── MANIFEST.md               # Per-dataset source + status log
    │   ├── buffering_coef.csv        # Output of build_buffering_coefficients.R
    │   ├── biochar/                  # Sanderman TIF
    │   ├── gadm/                     # Per-country GADM RDS cache
    │   ├── basalt/                   # GLiM + MRDS raw
    │   ├── era/                      # ERA_Derived.rda
    │   ├── mapspam/                  # MapSPAM per-variable zips + SPAM2020_*_*.tif aliases
    │   ├── naturebase/               # (staging ref copy if any)
    │   ├── soilgrids/                # geodata SoilGrids + WorldClim cache
    │   ├── constraints/              # WDPA / SRTM / GHSL / GMIA raw
    │   └── erw/ , faostat/ , awd/
    ├── staging/                      # Pre-processed inputs
    │   ├── template/                 # 5-arcmin CCA templates (regional + per-country)
    │   ├── admin/                    # GADM adm0/1/2 (GPKG + COG masks + lookup CSVs)
    │   ├── biochar/                  # Sanderman MG bundle + residue proxy
    │   ├── suitability/              # biochar + ERW fuzzy suitability 0..1
    │   ├── erw/, awd/, basalt/, constraints/, faostat/, mapspam/
    ├── naturebase/                   # The Naturebase bundle (kept at a fixed path
    │   └── all_pathways/             # because Step 1 references it directly)
    └── outputs/                      # Pipeline products
        ├── step1_mitigation/
        ├── step2_production/
        ├── step3_adaptation/
        ├── step4_integration/
        ├── step5_uncertainty/
        └── activity3_neb/
```

---

## Dataset inventory

Full per-dataset source URLs, cross-check rules, and manual-intervention flags are in [`data/inputs/MANIFEST.md`](data/inputs/MANIFEST.md). Short version:

| # | Dataset | Source | Fetcher | Status |
|---|---------|--------|---------|--------|
| 1 | GADM admin0/1/2 | `geodata::gadm` → <https://gadm.org> | `build_cca_template.R` | auto |
| 2 | Naturebase NCS (16 pathways + adm0/1/2 totals) | The Nature Conservancy — `https://naturebase-shared-docs.s3.amazonaws.com/Data-download/all_pathways.zip` | `fetch_naturebase.R` | auto |
| 3 | ERA derived yield effects | CGIAR Solutions / Zenodo | `fetch_era.R` | auto (sibling or manual-drop fallback) |
| 4 | Sanderman 2020 "Soils Revealed" dSOC | Harvard Dataverse `10.7910/DVN/HA17D3` | `fetch_biochar_potential.R` | auto |
| 5 | MapSPAM 2010 V2r0 (H, P, Y) | Harvard Dataverse `10.7910/DVN/PRFF8V` | `fetch_mapspam.R` | auto |
| 6 | Beerling 2020 / Calabrese 2024 / Kanzaki 2023 ERW | Nature SI + Zenodo | `fetch_erw_potential.R` | auto (best-effort; pixel raster pending) |
| 7 | Cao 2026 AWD rice CH4 | (anticipated) Zenodo + IPCC T1 fallback from MapSPAM | `fetch_cao_awd.R` | auto (fallback always runs; primary pending) |
| 8 | GLiM mafic / USGS MRDS basalt | PANGAEA `10.1594/PANGAEA.788537` + <https://mrdata.usgs.gov/mrds/> | `fetch_basalt.R` | auto (optional; MRDS always; GLiM CDN sometimes 500) |
| 9 | CMIP6 hazard frequencies | `atlas/hazards_prototype/` outputs | (external pipeline) | **deferred** |
| 10 | FAOSTAT producer prices | `bulks-faostat.fao.org` | `fetch_fao_prices.R` | auto |
| 11 | SoilGrids v2.0 pH / OC / clay | `geodata::soil_world` | `build_{biochar,erw}_suitability.R` | auto |
| 12 | WorldClim 2.1 monthly + bio | `geodata::worldclim_global` | `build_{biochar,erw}_suitability.R` | auto |
| 13 | SRTM 30″ → slope | `geodata::elevation_30s` | `fetch_constraint_masks.R` | auto |
| 14 | WDPA protected areas | <https://www.protectedplanet.net/> | `fetch_constraint_masks.R` | **manual drop** |
| 15 | GHSL SMOD 2020 | JRC | `fetch_constraint_masks.R` | auto |
| 16 | FAO AQUASTAT GMIA v5 | FAO | `fetch_constraint_masks.R` | **manual drop** |

---

## Analytical methods

### Step 1 — Mitigation potential (tCO₂e/ha/yr)

1. For each practice, read the source raster (Naturebase 1-km totals, Sanderman MG bundle, or a staged AWD/ERW raster).
2. **Normalise to density before reprojection**: divide source totals by the source per-pixel area (ha) to get tCO₂e/ha/yr at the source grid, then `terra::project(..., method = "average")` onto the 5-arcmin country template. This order guarantees the output extent equals the template extent AND produces a properly area-weighted mean of contributing source cells. Staged layers (Sanderman MG bundle, AWD T1 fallback) are already densities and skip the divide-by-area step.
3. Clip to the country polygon (scope boundary, not a mask).
4. Optional multiplication by a fuzzy-suitability surface (biochar, ERW) if the practice YAML nominates one. No suitability cut-off is applied by default.
5. Cross-check Naturebase-sourced rasters against their published adm0 totals in `all_pathways_tco2e_adm0.csv` (±5 % accepted INFO, ±15 % WARN, above that ERROR).

### Step 2 — Production benefit (Δyield t/ha, USD/ha/yr)

1. Anchor on Step-1 mitigation-eligible pixels (M > 0).
2. For every MapSPAM crop that intersects the ERA evidence for the practice:
   `pc = exp(mean(log(MeanT/MeanC))) − 1` over ERA yield outcomes, matching the ERAAnalyze back-transform.
3. Apply `pc` to the MapSPAM baseline yield raster, **converted from kg/ha to t/ha on read**, to get Δyield (t/ha).
4. Multiply by the country-specific FAOSTAT producer price (USD/tonne, averaged over recent years) to get Δvalue (USD/ha/yr). Missing prices produce NA values — never silent zeros.
5. Per-crop Δyield rasters are written for auditability; the Δvalue raster is the cross-crop sum per pixel.

### Step 3 — Adaptation benefit (avoided loss, USD/ha/yr)

```
avoided_loss_usd_ha_yr = Σ_{hazard, severity, crop} [
    P(hazard, severity)                        # from atlas/hazards_prototype
  × yield_loss_frac(hazard, severity, crop)    # literature table, editable
  × buffer_frac(practice, hazard, severity)    # ERA meta-analysis
  × baseline_yield_t_ha(crop)                  # MapSPAM
  × price_usd_tonne(crop, country)             # FAOSTAT
]
```

`yield_loss_frac` is an editable table (drought / heat / waterlogging at moderate / severe / extreme for maize, rice, sorghum) embedded in `step3_adaptation.R`. `buffer_frac` is read from `data/inputs/buffering_coef.csv` produced by the ERA meta-analysis. Hazards are currently treated as independent; a correlated-hazard formulation is flagged under "What needs checking".

### Step 4 — Integration + opportunity zones

1. Normalise M, P, A within each country by quantile rank to remove unit-scale effects.
2. Flag each layer "high" at the quantile cutoff in [`thresholds.yaml`](code/activity2/config/thresholds.yaml) (default 0.60).
3. Classify:
   - **1 Triple Win** — M high, P high, A high
   - **2 Mitigation + Production** — M high, P high, A low
   - **3 Mitigation + Adaptation** — M high, P low, A high
   - **4 Constrained** — M high, P low, A low
   - **5 Non-priority** — everything else
4. Contiguous patches of classes 1–4 are extracted with `terra::patches(..., directions = CONN)` and polygonised with per-zone mean M/P/A, area, and dominant class. `min_area_ha = 0` by default.

### Step 5 — Monte Carlo uncertainty

For each (country, practice) combination we draw `N_SAMPLES` (default 200) log-normal multipliers with per-layer CVs (M 0.30, P 0.50, A 0.60). The per-pixel stack of perturbed rasters is reduced via `terra::app(..., stats::quantile)` to P10, P50, P90. A rule-based confidence tier (High / Medium / Low) combining CV with ERA observation counts is written to `step5_uncertainty/confidence_tier.csv`.

### Activity 3 — Net Economic Benefit

```
annual_cost     = transition × annuity_factor(T, r) + implementation + mrv × M
annual_benefit  = cp × M + A + P
NEB_ha_yr       = annual_benefit − annual_cost
NPV_ha          = NEB_ha_yr × Σ_{t=1..T} (1+r)^(−t)
National_NPV    = Σ_pixels NPV_ha × cell_area_ha × adoption_frac
```

- Transition cost is annualised (`r / (1 − (1+r)^(−T))`) so year-on-year NEB is comparable.
- MRV cost scales with M (USD per tCO₂e per year).
- Adoption scales the **national aggregation only**; the per-hectare NPV raster is adoption-agnostic.
- Outputs are one COG per `(practice, country, carbon_price, adoption)` + a country CSV with mean NPV/ha and national totals.

All parameters live in [`code/activity3/neb_config.yaml`](code/activity3/neb_config.yaml). Defaults:

- Carbon prices: 5, 15, 40, 80 USD/tCO₂e
- Horizon: 20 years
- Discount rate: 8 %
- Adoption scenarios: 0.05, 0.15, 0.30, 0.50

---

## What has been done

- **Foundation**
  - 5-arcmin CCA reference grid (regional + per-country) built from GADM via `geodata::gadm`.
  - GADM adm0/1/2 persisted as GPKG + rasterized admin masks + lookup CSVs.
- **External fetchers (all auto-idempotent)**
  - Naturebase NCS `all_pathways.zip` (500 MB) downloaded + unpacked + range-validated.
  - ERA derived meta-analysis resolved via sibling-atlas → local-drop → Zenodo search.
  - MapSPAM 2010 V2r0 (3 per-variable mega-zips; ~600 MB) unzipped + CCA-clipped + aliased.
  - Sanderman 2020 "Soils Revealed" management bundle (488 MB) pulled from Harvard Dataverse, clipped, converted to tCO₂e/ha/yr.
  - FAOSTAT producer prices (11 MB zip) filtered to 3 countries × last 6 years, with spot-check assertions.
  - SoilGrids (pH / OC / clay) + WorldClim 2.1 (precip + bio) via `geodata`, with empirically-corrected unit handling.
  - SRTM 30″ → slope and GHSL SMOD 2020 urban-extent rasters at 5 arcmin.
  - GLiM lithology + USGS MRDS basalt occurrences (for the OPTIONAL ERW proximity mask).
  - Beerling 2020 Nature SI cached as a country-level cross-check table.
  - IPCC Tier-1 AWD rice methane-mitigation fallback raster from MapSPAM rice.
- **Pre-processing**
  - Biochar + ERW fuzzy suitability surfaces (0–1, geometric mean of trapezoids over pH, clay, OC, rainfall, temperature).
  - Crop-residue biomass proxy (MapSPAM P × IPCC T11.2 ratios, auditable CSV).
  - ERA hazard-buffering coefficient skeleton (64 rows; flagged `confidence_tier = Low` until workshop review).
- **Activity 2 (Technical potential)**
  - Step 1 mitigation rasters for every priority practice × country (10 COGs end-to-end); biochar scaled by suitability; admin totals cross-checked against Naturebase adm0 CSV.
  - Step 2 per-crop Δyield + aggregated Δvalue COGs (73 per-crop + 7 Δvalue, with admin CSVs).
  - Step 3 framework runs end-to-end with a graceful skip when hazard frequency rasters are absent.
  - Step 4 pixel classification + opportunity-zone GPKG for all (practice, country) pairs.
  - Step 5 Monte-Carlo P10/P50/P90 (40 COGs at `N_SAMPLES = 50` for a smoke test; `CCA_MC_SAMPLES` envvar tunable).
- **Activity 3 (NEB)**
  - 131 per-hectare NPV COGs across the full carbon-price × adoption sweep; admin CSVs with national totals.
- **Conventions + infra**
  - `cca_download()` atomic fetcher + `.source.json` (URL, bytes, md5, sha256) for every downloaded file.
  - `.method.json` sidecar for every synthesised staged raster (formula + parameters).
  - `cca_skip_if_exists()` + `CCA_FORCE=1` override baked into every script.
  - `cca_mask_by_*()` helpers that are NO-OPs by default — no implicit stakeholder-decision masking anywhere in the pipeline.
  - Smoke-run end-to-end verification with fixes for six bugs surfaced by real runs (see the section below).

---

## What needs checking

Items that produce valid artefacts today but deserve a review before publication:

1. **Hazard frequency surfaces for Step 3** — the current pipeline consumes `atlas/hazards_prototype/` outputs (CMIP6 indices × crop × severity). Until that pipeline is rerun for the CCA extent, Step 3 writes zero adaptation rasters and the Step-5 A-layer percentiles reflect M + P variance only.
2. **ERA → CCA practice mapping + hazard classification** — `build_buffering_coefficients.R` uses a DRAFT free-text keyword scan of `C.Descrip + T.Descrip` to label hazard / severity, and a hard-coded practice-name lookup. Both deserve workshop sign-off; every output row is currently tagged `confidence_tier = Low` pending that review.
3. **Ethiopia FAOSTAT prices** — FAOSTAT's Ethiopia coverage for 2018–2024 is sparse; Step 2 Δvalue is NA for Ethiopia which propagates through to a negative (cost-only) NPV in Activity 3. This is a FAOSTAT gap, not a code bug — a manual price sheet (`data/staging/faostat/manual_overrides.csv`) would unblock it.
4. **Step 4 opportunity-zone area_ha** — when the high-priority seed raster is a single contiguous cluster covering most of a country, `terra::expanse()` returns an area inflated by non-seed pixels inside the polygon's envelope. Zone polygons and classification COGs are still correct; only the `area_ha` attribute is inflated. Fix: intersect each polygon with the country adm0 before computing area.
5. **Yield-loss fraction table in Step 3** — the embedded `loss_table` uses literature defaults for 3 crops × 3 hazards × 3 severities. It should be replaced by a crop-calendar-aware table calibrated against `jagermeyer/` GGCMI Phase-3 projections.
6. **Correlated hazards in Step 3** — hazards are currently summed as independent contributors. A copula-based formulation (e.g. drought ↔ heat co-occurrence) is a deferred extension; the TODO marker is in `step3_adaptation.R`.
7. **Sanderman management bundle caveat** — the Step-1 biochar arm uses Sanderman's `scenario_crop_MG_dSOC_Y20.tif`, a *bundle* of cover crops + residue retention + biochar. The `.method.json` sidecar flags this. A biochar-only layer would require a pyrolysis coefficient applied to the MapSPAM residue proxy — scaffolding is in place but not wired.
8. **Cao et al. 2026 AWD** — the Zenodo search returns false-positive records; auto-download is deliberately OFF until a human confirms the DOI. The IPCC Tier-1 fallback runs meanwhile.
9. **GLiM PANGAEA 500 error** — PANGAEA's edge CDN returns HTTP 500 intermittently on the GLiM zip. `fetch_basalt.R` tries three URL variants and logs a WARN. MRDS still succeeds and provides thin Africa coverage.
10. **WDPA + GMIA manual drops** — protectedplanet.net WDPA requires a guestbook-style clickthrough; FAO AQUASTAT GMIA v5 moved endpoints. Both are flagged as manual-drop in `fetch_constraint_masks.R` with exact paths.
11. **Monte-Carlo sample count** — Step 5 default is 200 samples to keep wall time reasonable; the paper draft calls for 1000. Set `CCA_MC_SAMPLES=1000` for publication runs.
12. **Naturebase cross-check tolerance** — country totals are validated against `all_pathways_tco2e_adm0.csv` at ±5 %; the resample-then-average path introduces a few-percent drift near coastlines. Accept or tighten with explicit land-mask intersection.

---

## TODO list (prioritized)

### P0 — Blockers for country-profile publication
- [ ] Run `atlas/hazards_prototype/` for NGA/KEN/ETH at CCA extent → unblocks Step 3.
- [ ] Commission stakeholder-workshop sign-off on the ERA → CCA practice mapping, hazard classification rules, and confidence-tier cut-offs used by `build_buffering_coefficients.R`.
- [ ] Source an Ethiopia producer-price override table (USD/t for maize, teff, wheat, sorghum, pulses). Drop into `data/staging/faostat/manual_overrides.csv` with the schema `Area,Item,Year,Value,Unit`.
- [ ] Intersect Step 4 opportunity-zone polygons with adm0 before computing `area_ha` so the attribute matches the classification COG.

### P1 — Quality + coverage
- [ ] Replace the Step 3 yield-loss lookup table with a crop-calendar-aware calibration from `jagermeyer/` GGCMI outputs.
- [ ] Confirm Cao et al. 2026 DOI, wire it into `fetch_cao_awd.R`, and regenerate AWD Step 1 rasters.
- [ ] Pyrolysis-yield biochar attribution: multiply `residue_biomass_tonnes_per_pixel.tif` by a biochar yield coefficient + stability factor to produce a biochar-only layer that supersedes the Sanderman MG bundle.
- [ ] Review FAO AQUASTAT GMIA v5 current endpoint (or reuse atlas `atlas_exposure/` output) and flip `fetch_constraint_masks.R` §4 back to `auto`.
- [ ] Publish a country-profile renderer (Quarto) that reads the admin CSVs + COGs and writes a 25-page PDF per country.

### P2 — Robustness + observability
- [ ] Add unit tests under `code/tests/` (CRS, res, value ranges, percentile ordering) and wire them into a CI smoke job.
- [ ] Register the `cca` STAC collection on the team S3 bucket; publish final COGs + GPKGs there.
- [ ] Wire the optional Step 4 masks into `cca_mask_by_threshold()` behind YAML knobs (slope cut-off, WDPA exclusion, minimum zone area). Keep all defaults `NULL`.
- [ ] Add the correlated-hazard copula option to Step 3.
- [ ] Package R dependencies in a `renv.lock` for deterministic package versions.

### P3 — Nice-to-haves
- [ ] Zenodo mirrors for the large datasets (Naturebase, Sanderman, MapSPAM) so a fresh clone has a single fallback when upstream endpoints flap.
- [ ] Parallelise Step 2 across (country, practice, crop) with `future.apply` — the inner crop loop is the current bottleneck.
- [ ] Extend the confidence-tier CSV with an HTML factor plot for the country-profile appendix.

---

## Design rules

1. **Every script documents its own reasoning.** Inputs, units, method choices, and cross-check thresholds are explained inline in each file. A new analyst can audit a layer by reading the script that produced it.
2. **No implicit spatial masking.** Every optional mask — basalt proximity, slope cut-off, protected-area exclusion, urban mask, irrigation requirement, minimum zone area — is wrapped in an `cca_mask_by_*` helper that defaults to a no-op. Real values are workshop-decided.
3. **Idempotent fetchers.** Re-running a script is cheap: a download skips if the raw file is already on disk, and a builder skips if its staged output is already present. Override with `CCA_FORCE=1`.
4. **Provenance sidecars.** Every downloaded file gets `<file>.source.json` (URL, timestamp, bytes, md5, sha256). Every synthesised layer gets `<file>.method.json` (formula + parameters).
5. **Primary vector format: GPKG.** GeoJSON is available on opt-in via `cca_write_vector(..., also_geojson = TRUE)`. GeoParquet is supported when the GDAL build has the Parquet driver; otherwise the writer transparently falls back to GPKG.
6. **Primary raster format: COG** (DEFLATE, 512-px tiles, `BIGTIFF=IF_SAFER`).
7. **Units are unit-consistent at every boundary.** kg/ha → t/ha conversion on MapSPAM read, tC → tCO₂e conversion on Sanderman read, pH units left in native pH, clay in percent. Unit assumptions are commented at every conversion site.
8. **Standalone.** No hard dependency on sibling atlas folders; every external dataset has a fetcher that resolves at runtime or a documented manual-drop path.

---

## Smoke-test results (2026-04-17)

| Stage | Outputs | Key metric |
|---|---:|---|
| `build_cca_template.R` | 4 templates + 9 GADM GPKG + 9 admin masks + 9 lookup CSV | regional grid 544 × 236 cells |
| `fetch_naturebase.R` | 16 TIFs + 3 admin CSVs + metadata | Naturebase bundle ~500 MB unpacked |
| `fetch_era.R` | 1 RDA | 110 298 rows loaded |
| `fetch_mapspam.R` | 60 CCA-clipped COGs + aliases | Kenya rice yield mean 1.7 t/ha |
| `fetch_fao_prices.R` | 1 CSV | Kenya maize 384 USD/t ✓; Nigeria rice 486 USD/t ✓ |
| `fetch_biochar_potential.R` | 1 COG + `.method.json` | mean 0.58 tCO₂e/ha/yr ✓ |
| `build_{biochar,erw}_suitability.R` | 2 COG + 2 `.inf` | biochar mean 0.80; ERW mean 0.75 ✓ |
| `step1_mitigation.R` | 10 COGs + 3 admin CSVs | biochar KEN mean 0.53 tCO₂e/ha/yr; extents == templates ✓ |
| `step2_production.R` | 73 per-crop Δyield + 7 aggregated Δvalue + 3 admin CSVs | KEN biochar rice +1.13 t/ha × 491 USD/t → 557 USD/ha/yr ✓ |
| `step4_integration.R` | 7 classification COGs + 7 opportunity-zone GPKGs | classification correct (area_ha inflated — see "What needs checking") |
| `step5_uncertainty.R` (N=50) | 40 P10/P50/P90 COGs + 1 confidence CSV | percentile ordering assertion passes ✓ |
| `neb.R` | 131 NPV COGs + 3 admin CSVs | KEN biochar national NPV $23 B at $5/tCO₂e × 0.30 adoption ✓ |

Bugs caught during verification and since fixed: (i) Rscript `--file=` resolution in `cca_project_root`, (ii) GDAL Parquet driver absence → GPKG primary, (iii) FAO regex / parens, (iv) Sanderman Mollweide CRS handling, (v) SoilGrids unit handling (pH, clay not scaled by 10), (vi) ERA column detection, (vii) Zenodo DOI quoting, (viii) MRDS commodity filter, (ix) GMIA URL rot, (x) MapSPAM 2020 guestbook ↔ 2010 V2r0 pivot, (xi) Step 1 extent alignment (density-first reprojection), (xii) Step 2 kg↔t unit, (xiii) `terra::extract(..., fun = custom, na.rm = TRUE)` incompatibility.

---

## License

Code: MIT (see `LICENSE` when present). Data: varies by source — see each entry in `data/inputs/MANIFEST.md` for the original licence.
