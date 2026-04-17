# CCA input data manifest

Every file consumed by the CCA pipeline is logged here. One section per dataset. "Status" values:

- **local** — already present in the repo (no fetch needed).
- **auto** — fetched by a script in `code/inputs/` when the pipeline runs.
- **pending** — not yet sourced; listed so it is obvious what is missing.

All paths below are relative to the repo root `ccadapt/`. Resolution and CRS are noted so cross-checks are straightforward.

---

## Already local

### 1. Naturebase NCS pathways (16 rasters + adm0/1/2 CSVs)

- **Status:** local
- **Location:** `data/naturebase/all_pathways/`
- **Source:** <https://naturebase.org> (Nature Conservancy, retrieved 2025-04-15)
- **Citation:** See `data/naturebase/all_pathways/Naturebase_metadata.csv` per pathway.
- **CRS / resolution:** EPSG:4326 / 0.00898315° (~1 km at equator)
- **Units (rasters):** Total tCO₂e per 1-km pixel per year (filename `*_total.tif`).
- **Used by:** `code/activity2/step1_mitigation.R` for `agroforestry_cropland`, `cover_crops`, `grassland_soil_carbon`, `grassland_reduced_emissions`, `grassland_restoration`, `savanna_agroforestry`.
- **Cross-check:** admin-0 totals compared against `all_pathways_tco2e_adm0.csv` (±5 % acceptance).
- **Notes:** Naturebase embeds a cropland-agreement mask (Tubiello 2022). No further MapSPAM mask is applied by default — see `code/common/geo_utils.R::cca_mask_by_polygon` if a caller opts in.

### 2. Project text

- **Status:** local
- **Location:** `text/`
- **Contents:** Project scoping document (plus any stakeholder-workshop notes and policy PDFs that get added during Activity 1).
- **Used by:** humans (scope reference); not read programmatically.

---

## Fetched automatically by pipeline scripts

### 3. GADM administrative boundaries (adm0, adm1, adm2)

- **Status:** auto
- **Fetched by:** `code/inputs/build_cca_template.R` (wraps `geodata::gadm`)
- **Cache location:** `data/inputs/gadm/`
- **Source:** <https://gadm.org>
- **Countries:** NGA, KEN, ETH
- **CRS:** EPSG:4326
- **Output:** `data/staging/admin/gadm_<iso3>_<adm>.gpkg` (primary vector), `gadm_<iso3>_<adm>_<res>.tif` (rasterized mask), `gadm_<iso3>_<adm>_lookup.csv` (adm_id → GID/NAME lookup).
- **Notes:** `geodata::gadm` is used so country partners (national ministries in NGA/KEN/ETH) can cross-reference our admin IDs directly with their own GADM-based datasets. The atlas-wide GAUL 2024 source (see `atlas/spatial-foundation/build_atlas-boundaries_gaul24.r`) is available as an alternate but not used here.

### 4. 5-arcmin CCA template raster

- **Status:** auto
- **Built by:** `code/inputs/build_cca_template.R`
- **Output:** `data/staging/template/cca_template_{3country,NGA,KEN,ETH}_5arcmin.tif`
- **Resolution:** 5/60 degrees (0.08333…°, exact rational)
- **Extent:** union of NGA+KEN+ETH adm0 snapped outward to the nearest 5-arcmin cell aligned with the global (–180, –90) grid origin.

---

## Fetch scripts ready (download or build on demand)

### 5. AWD rice mitigation raster

- **Status:** auto (with fallback) — primary Cao 2026 source is anticipated; fallback builds an IPCC-Tier-1 AWD raster from MapSPAM rice area.
- **Practice:** `awd_rice`
- **Blocker for:** Step 1 Nigeria AWD arm (full precision until Cao 2026 lands).
- **Fetch script:** [fetch_cao_awd.R](../../code/inputs/fetch_cao_awd.R)
- **Output:** `data/staging/awd/awd_potential_tco2e_ha_yr_T1fallback.tif` (fallback) or Cao raster when a Zenodo DOI resolves.
- **Method tag:** written into `<tif>.method.json` (`IPCC_Tier1_fallback` vs `Cao2026`).
- **Cross-check:** Cao 2026 primary source is searched on Zenodo by creator + "alternate wetting and drying"; auto-download is deliberately OFF until a human confirms the record identity.

### 6. Biochar mitigation potential

- **Status:** auto
- **Source:** Sanderman J. et al. 2020, "Soils Revealed soil carbon futures" — Harvard Dataverse DOI [`10.7910/DVN/HA17D3`](https://doi.org/10.7910/DVN/HA17D3). File `scenario_crop_biochar_dSOC_Y20.tif` (delta-SOC over 20 yr).
- **Why Sanderman (not Woolf 2021):** Woolf 2021 publishes a country-level inventory model, no raster. Naturebase itself cites Sanderman 2020 as the biochar basis; keeping our upstream identical to theirs makes the numbers directly comparable.
- **Fetch script:** [fetch_biochar_potential.R](../../code/inputs/fetch_biochar_potential.R) (resolves file ID via Dataverse API, downloads, converts tC/ha over Y20 → tCO₂e/ha/yr, clips to CCA extent)
- **Output:** `data/staging/biochar/biochar_potential_tco2e_ha_yr.tif` (COG, native ~1 km)
- **Downstream:** multiplied by biochar suitability (§9) and residue-availability proxy (§14) in Step 1.

### 7. Enhanced rock weathering (ERW) mitigation potential

- **Status:** auto (best-effort)
- **Sources (in preference order):**
  1. Calabrese et al. 2024 — Zenodo DOI `10.5281/zenodo.10694306` (preferred)
  2. Kanzaki et al. 2023 SCEPTER v1 — Zenodo DOI `10.5281/zenodo.7625155`
  3. Beerling et al. 2020 Nature SI XLSX (country totals only, for cross-check) <https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2448-9/MediaObjects/41586_2020_2448_MOESM4_ESM.xlsx>
- **Fetch script:** [fetch_erw_potential.R](../../code/inputs/fetch_erw_potential.R)
- **Output:** `data/staging/erw/erw_potential_tco2e_ha_yr.tif` when a raster source resolves; otherwise just the Beerling SI cached for country cross-check.

### 8. Basalt / mafic geology (OPTIONAL ERW proximity input)

- **Status:** auto
- **Primary:** GLiM (Hartmann & Moosdorf 2012) at PANGAEA DOI `10.1594/PANGAEA.788537` — filtered to mafic / intermediate volcanics + basic plutonics (`va`, `vi`, `pb`).
- **Supplementary:** USGS MRDS <https://mrdata.usgs.gov/mrds/> — filtered to basalt-commodity occurrences.
- **Fetch script:** [fetch_basalt.R](../../code/inputs/fetch_basalt.R)
- **Output:** `data/staging/basalt/glim_africa_mafic_intermediate.gpkg` and `data/staging/basalt/mrds_basalt_points.gpkg`.
- **Used by:** `cca_mask_by_proximity()` ONLY when a caller passes a non-NULL `distance_km`. No default.

### 9. Biochar + ERW fuzzy suitability rasters (0–1)

- **Status:** auto (R implementation; CropSuiteLite `.inf` side-car written for optional full-simulator reruns)
- **Produced by:**
  - [build_biochar_suitability.R](../../code/inputs/build_biochar_suitability.R) — pH / OC / clay / rainfall
  - [build_erw_suitability.R](../../code/inputs/build_erw_suitability.R) — pH / clay / temperature / rainfall
- **Dependencies (auto via geodata):** SoilGrids (`soil_world`), WorldClim 2.1 (`worldclim_global`).
- **Output:** `data/staging/suitability/{biochar,erw}_suitability_5arcmin.tif` + `biochar.inf`, `erw.inf` CropSuiteLite side-cars.
- **Note:** membership trapezoids are piecewise-linear; combination is geometric mean so any zero driver zeros the pixel. Thresholding is deferred to Activity-1 stakeholder workshop.

### 10. Hazard-buffering coefficient table

- **Status:** auto (skeleton; flagged `confidence_tier = "Low"` until workshop review)
- **Produced by:** [build_buffering_coefficients.R](../../code/inputs/build_buffering_coefficients.R)
- **Input:** `atlas/solutions/data/ERA_Derived.rda` (already bundled with the atlas); optional merge with Adaptation Insights CSV if placed at `data/inputs/adaptation_insights/buffering.csv`.
- **Output:** `data/inputs/buffering_coef.csv` with columns `practice, era_practice, hazard, severity, crop, buffer_frac, n_studies, n_observations, ci_low, ci_high, source, confidence_tier`.
- **Caveats:** ERA context-to-hazard mapping rules, within-study pairing strategy, confidence-tier cut-offs, and Adaptation Insights merge are all flagged for Activity-1 stakeholder sign-off. The skeleton runs with safe defaults and re-runs automatically when the config is updated.

### 11. FAO farm-gate prices (2022–2024)

- **Status:** auto
- **Source:** FAOSTAT bulk "Producer Prices – Annual (Normalized)" <https://bulks-faostat.fao.org/production/Prices_E_All_Data_(Normalized).zip>
- **Fetch script:** [fetch_fao_prices.R](../../code/inputs/fetch_fao_prices.R)
- **Output:** `data/staging/faostat/fao_producer_prices_usd_per_tonne_2022_2024.csv` (Area / Item / Year / USD_per_t).
- **Cross-check:** spot-checks against Kenya maize, Nigeria rice, Ethiopia maize price bands logged as INFO or WARN.

### 12. Country policy / NDC / NAP PDFs

- **Status:** pending (manual)
- **Source:** UNFCCC NDC registry, country climate strategy portals.
- **Target location:** `text/policy/<iso3>/`
- **Used by:** Activity 1 (non-code) + Activity 4 scenario-D parameterisation.

### 13. Constraint-mask layers (OPTIONAL; default off)

- **Status:** auto
- **Sources:**
  - WDPA — <https://www.protectedplanet.net/downloads/WDPA_current_<ISO3>?type=shapefile> (per-country shapefile)
  - SRTM 30″ — via `geodata::elevation_30s` per country, slope derived with `terra::terrain(..., unit="degrees")`
  - GHSL SMOD 2020 — JRC release R2023A 1 km <https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_SMOD_GLOBE_R2023A/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000/V1-0/GHS_SMOD_E2020_GLOBE_R2023A_54009_1000_V1_0.zip>
  - FAO AQUASTAT GMIA v5 irrigation extent — <https://www.fao.org/nr/water/aquastat/irrigationmap/gmia_v5_aei_pct_aei.asc.gz>
- **Fetch script:** [fetch_constraint_masks.R](../../code/inputs/fetch_constraint_masks.R)
- **Output:** `data/staging/constraints/{wdpa_<iso3>.gpkg,slope_deg_<iso3>_5arcmin.tif,ghsl_smod_5arcmin.tif,gmia_aei_pct_5arcmin.tif}`
- **Applied only when:** a caller of `cca_mask_by_polygon()` / `cca_mask_by_threshold()` explicitly passes a non-NULL argument.

### 14. Crop-residue biomass proxy (for biochar)

- **Status:** auto (requires MapSPAM outputs from §15 or from `atlas/atlas_exposure/`)
- **Produced by:** [build_residue_proxy.R](../../code/inputs/build_residue_proxy.R)
- **Method:** per-pixel sum of `production × residue_ratio` across MapSPAM crops, ratios from IPCC 2019 Refinement Vol 4 Table 11.2 (csv written beside the raster for auditability).
- **Output:** `data/staging/biochar/residue_biomass_tonnes_per_pixel.tif` + `residue_indices_ipcc_t11.2.csv`.

### 15. MapSPAM 2010 V2r0 crop distribution rasters

- **Status:** auto
- **Source:** IFPRI MapSPAM 2010 V2r0, Harvard Dataverse DOI [`10.7910/DVN/PRFF8V`](https://doi.org/10.7910/DVN/PRFF8V). File IDs are resolved at runtime via the Dataverse API; a DOI change requires only updating `MAPSPAM_DOI` in [`fetch_mapspam.R`](../../code/inputs/fetch_mapspam.R).
- **Why 2010 V2r0 (not 2020 V2r0):** the 2020 package (`10.7910/DVN/SWPENT`) requires a Dataverse Guestbook response before any file download, which blocks unattended API access. The 2010 V2r0 package is open-access and is the version `atlas/atlas_exposure/` uses, so CCA numbers remain directly comparable to atlas benchmarks.
- **Fetched by:** [fetch_mapspam.R](../../code/inputs/fetch_mapspam.R)
- **Layout quirk:** files ship as three per-variable mega-zips (`…_harv_area.geotiff.zip`, `…_prod.geotiff.zip`, `…_yield.geotiff.zip`). The fetcher downloads each, unzips once, copies the wanted `_<CROP>_A.tif` files, CCA-clips them to the regional bbox, and writes COGs.
- **Crops pulled:** RICE, MAIZ, WHEA, SORG, MILL, BARL, SOYB, GROU, CASS, YAMS, POTA, SWPO, COTT, SUGC, BANA, PLNT, OOIL, COFF, COCO, OCER.
- **Variables pulled:** H (harvested area, ha), P (production, tonnes), Y (yield, **kg/ha** — Step 2 converts to t/ha on read so arithmetic with FAO USD/tonne prices is unit-consistent).
- **Output:** `data/staging/mapspam/spam2010v2r0_global_{H,P,Y}_<CROP>_A_cca.tif` + `data/inputs/mapspam/SPAM2020_{H,P,Y}_<crop>_A.tif` aliases used by `fetch_cao_awd.R`, `build_residue_proxy.R`, Step 2, and Step 3.
- **Used by:** Step 2 production benefit, Step 3 adaptation benefit, biochar residue proxy, AWD IPCC Tier-1 fallback.

---

## Change log

Entries get dated when a dataset changes status (pending → auto → local). Keep the most recent first.

- 2026-04-17 (late) — added §15 MapSPAM fetcher (IFPRI Dataverse). Originally targeted 2020 V2r0 (`10.7910/DVN/SWPENT`) but that package is guestbook-gated; pivoted to open-access 2010 V2r0 (`10.7910/DVN/PRFF8V`). Wired as the canonical source for Step 2 baseline yield, the biochar residue proxy, and the AWD IPCC T1 fallback.
- 2026-04-17 (late) — fetch/pre-process scripts written for AWD, biochar, ERW, basalt, suitability, FAO prices, constraint masks, residue proxy, buffering coefficients. Status flipped from pending → auto for items 5–11, 13, 14. Item 12 (policy PDFs) remains manual.
- 2026-04-17 — initial manifest created; Naturebase local, GADM auto, all others pending.
