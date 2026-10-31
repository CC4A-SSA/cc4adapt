# Adaptation-factor formulation exploration

**EXPLORATORY. Nothing here is an estimate.** Loss shares (L) are bronze placeholders, buffering (m) is a scenario grid, ERA gains (g) are pooled across seasons. See "Status of inputs" below and `README_EXPLORATORY.md` in any outputs folder. Tracking: CC4A-SSA/cc4adapt#1. Brief and decisions: project `_claude/notes/` (2026-10-26 BRIEF, 2026-10-27 DECISIONS).

## Question

How does the choice of avoided-loss formulation (levels DiD vs proportional split vs the current pipeline) change the adaptation factor and the production/adaptation split, across real Atlas hazard regimes? Totals are identical under every formulation; only the label moves, and only the adaptation share is claimed as a carbon-finance co-benefit.

## Run

From the repo root:

```bash
Rscript code/explore/af_formulations/run_all.R        # tests, then steps 01-06
Rscript code/explore/af_formulations/run_all.R 3 6    # steps 3-6 only
```

| Step | Script | Does | Writes |
|---|---|---|---|
| 01 | `01_fetch_inputs.R` | Pulls pinned Atlas `prod_t` tiers, exposure denominator, GAUL 2024 admin-2, HarvestChoice AEZ16, WorldClim 10′, MapSPAM maize COG windows; copies Dixon FS rasters | `data/staging/af_formulations/` + `MANIFEST.csv` (source, pin, md5) |
| 02 | `02_check_inputs.R` | Identity checks from the hazard briefing (classes sum to `any`; `any+none` stable; marginals monotone; GAUL join) | console |
| 03 | `03_params.R` | ERA g per practice (all-SSA and by AEZ16simple); expands the parameter grid | `params/era_g.csv`, `params/param_grid.csv` |
| 04 | `04_units.R` | One row per admin-2 unit: MapSPAM maize B_obs, dominant Dixon FS1/FS2, AEZ classes, bimodal flag | `units.parquet` |
| 05 | `05_compute.R` | P(h, s) from Atlas exclusive classes (admin-1 fallback); B_norm; F0–F5 per unit × parameter set | `af_results_maize_historic.parquet`, `units_no_p.csv` |
| 06 | `06_diagnostics.R` | Sign, gap-scaling, label-flip, strata, severity-mix, bimodal tables | `diagnostics/d1..d6.csv` |
| 03b | `03b_era_site_years.R` | ERA maize comparisons for the two practices by site × year × season slot | `era_site_years.csv`, `era_comparisons.csv` |
| 07 | `07_era_seasons.py` | Classify each ERA site-season by hazard band from observed daily climate (CHIRPS v3 + AgERA5 via Earth Engine, GGCMI maize window, Atlas bands). Needs `GCP_PROJECT_ID` and `PYTHONPATH=../climate-toolkit`; `--limit N` for a dry run, resumable | `era_site_season_hazards.csv` |
| 08 | `08_era_split.R` | g_norm vs g_stress per practice × hazard; proportional-vs-additive test on stress seasons | `era_split/*.csv` |
| 09 | `09_obs_vs_atlas.R` | Observed stress frequency at ERA sites vs the Atlas P of the containing admin-2 unit | `era_split/obs_vs_atlas_sites.csv` |

`07 --window ggcmi` reruns the classification on the GGCMI (Atlas jagermeyr) window to separate the window effect from the GCM effect; output `era_site_season_hazards_ggcmi.csv`.

| 10 | `10_gcm_vs_obs.py` | Same sites, window and period (1995–2014) scored on observed climate (`--models ""`) and on NEX-GDDP historical runs (`--models ACCESS-CM2 --tag _access-cm2`, one model per run) | `gcm_vs_obs_sites*.csv` |
| 11 | `11_gcm_summary.R` | Per-site P(dry ≥ moderate) by source (observed all years, observed ERA years, each GCM) beside the Atlas P of the unit | `era_split/gcm_vs_obs_summary.csv` |
| 12 | `12_export_page_data.R` | Builds `docs/af_formulations/data.js` for the interactive page from 05, 06, 08, 09, 11 outputs | `docs/af_formulations/data.js` |
| 13 | `13_kenya_dotplot.R` | Static PNG of the Kenyan per-site comparison | `era_split/kenya_gcm_vs_obs.png` |

Library: `R/formulations.R` (formulations, band conversion, expected annual), `R/paths.R`. Tests: `tests/testthat/` (`Rscript -e 'testthat::test_dir("code/explore/af_formulations/tests/testthat")'`). Runtime identity checks that must pass on every run: F1 = F2 = F3 expected totals; expected F1 − F3 gap = −g·B·E[L]; E[L] in 06 equals E[L] in 05.

Steps 03b–13 are the ERA season-classification pilot (science review findings 2 and 4).

**Tools.** Season-by-season stress at ERA sites is computed with the [CGIAR climate-toolkit](https://github.com/CGIAR-Climate-Data-Hub/climate-toolkit) (`evaluate_hazards()`; commit `453c12a`): it fetches daily CHIRPS v3 + AgERA5 (or NEX-GDDP) through Earth Engine, runs its soil-water balance, and returns season counts of NDWS / NDWL0 / NTx35 for the window passed in. This repo supplies the windows (03c), converts counts to Atlas monthly-mean bands with Atlas thresholds (07), matches seasons to ERA comparisons and does the statistics (08), and joins to Atlas district frequencies (09, 11). Two toolkit workarounds live in 10: NEX model/scenario via the source module's `DEFAULT_MODEL` / `DEFAULT_SCENARIO`, and fetching to 30 June of the following year for windows crossing 31 December. **The Atlas's own per-season classes are not used (server-only); all stress seasons here are Atlas-style, computed from weather records or climate-model output.**

`run_all.R` does not run 07; run it by hand:

```bash
GCP_PROJECT_ID=ee-peetmate PYTHONPATH=../climate-toolkit python3 code/explore/af_formulations/07_era_seasons.py --workers 2
```

## Reproduce from scratch

Everything below the data pulls is deterministic given the pinned inputs. Steps 07 and 10 need an Earth Engine project and take hours; everything else runs in minutes.

```bash
git clone https://github.com/CC4A-SSA/cc4adapt.git && cd cc4adapt
git clone https://github.com/CGIAR-Climate-Data-Hub/climate-toolkit.git ../climate-toolkit   # pinned below
Rscript code/inputs/fetch_era.R                                   # ERA_Derived.rda
export AF_FS_DIR=/path/to/Farming_systems                         # Dixon FS1/FS2 rasters (ask Pete)
Rscript code/explore/af_formulations/run_all.R                    # tests, 01-06, page data
Rscript code/explore/af_formulations/run_all.R --pilot            # 03b, 03c, then it stops and asks for 07
GCP_PROJECT_ID=<ee-project> PYTHONPATH=../climate-toolkit \
  python3 code/explore/af_formulations/07_era_seasons.py --workers 2            # ~1 h; resumable
GCP_PROJECT_ID=<ee-project> PYTHONPATH=../climate-toolkit \
  python3 code/explore/af_formulations/07_era_seasons.py --window ggcmi --workers 2
GCP_PROJECT_ID=<ee-project> PYTHONPATH=../climate-toolkit \
  python3 code/explore/af_formulations/10_gcm_vs_obs.py --iso KEN --models ""                       # observed, ~1 h
GCP_PROJECT_ID=<ee-project> PYTHONPATH=../climate-toolkit \
  python3 code/explore/af_formulations/10_gcm_vs_obs.py --iso KEN --models ACCESS-CM2 --tag _access-cm2
GCP_PROJECT_ID=<ee-project> PYTHONPATH=../climate-toolkit \
  python3 code/explore/af_formulations/10_gcm_vs_obs.py --iso KEN --models GFDL-ESM4 --tag _gfdl-esm4
AF_HAZ_FILE=data/staging/af_formulations/era_site_season_hazards_ggcmi.csv AF_HAZ_TAG=_ggcmi Rscript code/explore/af_formulations/08_era_split.R
AF_HAZ_FILE=data/staging/af_formulations/era_site_season_hazards_ggcmi.csv AF_HAZ_TAG=_ggcmi Rscript code/explore/af_formulations/09_obs_vs_atlas.R
Rscript code/explore/af_formulations/run_all.R --pilot            # 08, 09, 11, 13, 12
```

**Pins.** Atlas tiers and denominator by S3 `versionId` (in `01_fetch_inputs.R`); GAUL 2024 and MapSPAM COGs by version id checked at fetch; HarvestChoice AEZ16 r2.0 by md5; climate-toolkit at commit `453c12a` (2026-07-27) with Earth Engine NEX-GDDP `1.1`; ERA_Derived.rda as fetched by `code/inputs/fetch_era.R` (the copy used here was dated 2026-04-20 and has the `yi`/`pc` issue noted in `03_params.R`). Not pinned: CHIRPS v3 and AgERA5 as served by Earth Engine on the run date (Oct 2026), and the Dixon rasters (local copy).

**Software.** R ≥ 4.4 with `data.table`, `arrow`, `duckdb`, `terra`, `sf`, `yaml`, `geodata`, `jsonlite`, `ggplot2`, `testthat`. Python 3.12 with `pandas` and the climate-toolkit's requirements (`earthengine-api`, `xarray`, `xee`). An Earth Engine account with a Cloud project.

**Checks a re-run must pass.** The 70 unit tests; 02's identity checks; 05's F1 = F2 = F3 totals; 06's gap identity and E[L] agreement; 12 reproducing `data.js` without a diff when nothing upstream changed.

### Environment variables

| Var | Default | Meaning |
|---|---|---|
| `AF_OUT_DIR` | `data/outputs/af_formulations` (gitignored) | canonical outputs |
| `AF_RECORD_DIR` | unset | optional second copy (Pete uses the project `_claude/outputs/adaptation_factor`); only used if set and its parent exists |
| `AF_FS_DIR` | `../AICCRA_EFA/Data/Farming_systems` | Dixon farming-system rasters (not downloadable; ask Pete) |
| `GCP_PROJECT_ID` | unset | Earth Engine project for steps 07 and 10 (CGIAR-Climate-Data-Hub/climate-toolkit checkout on `PYTHONPATH`) |
| `AF_HAZ_FILE`, `AF_HAZ_TAG` | site-window file, "" | which 07 classification 08/09 read and how their outputs are suffixed |

### Prerequisites not fetched by 01

- `data/inputs/era/ERA_Derived.rda` — run `code/inputs/fetch_era.R`.
- Dixon FS rasters — see `AF_FS_DIR`.
- R ≥ 4.4 with `data.table`, `arrow`, `duckdb`, `terra`, `sf`, `yaml`, `geodata`, `testthat`.

## Formulations

Per hazard h, severity s: B normal-year baseline yield, A = B(1+g), L control loss share, m share of loss the practice avoids (signed). Stress year: Bs = B(1−L), As = A(1−L(1−m)).

| ID | Production (stress yr) | Avoided loss (stress yr) | Notes |
|---|---|---|---|
| F0 | A − B | As − Bs | v0.1; double counts |
| F1 | A − B | L(B − (1−m)A) | levels DiD; what MELIAF v3.0 defines (g in kg/ha) |
| F2 | Δs − max(0, F1) | max(0, F1) | F1 floored at 0, per band |
| F3 | (A − B)(1 − L) | A·L·m | proportional split; relative g × stress-year control yield |
| F4 | as F3 | A·L·m/(1 − L) | relative in-year ratios × normal-year B |
| F5 | A − B | B·L·m | current cc4adapt step2 + step3; each arm's own loss × B |

Decomposition: Δs = Δn − L·Δn + B·L·m + g·B·L·m. F3 puts the hazard bite on the gain in production and the loss avoided on the extra yield in adaptation; F1 nets the two into adaptation (F1 = F3 − L·B·g); F5 counts neither (F5 = F3 − L·B·g·m). Expected annual gap F1 − F3 = −g·B·E[L]; relative gap = −g/((1+g)·m), independent of location for constant m; F1 changes sign at m\* = g/(1+g).

## Status of inputs

| Input | Status |
|---|---|
| P(h, s) | Atlas NEX-GDDP 1995–2014 ensemble, admin-2 production-weighted, generic NDWS/NTx35/NDWL0 "≥ moderate" tiers. About a third of unit-hazard pairs fall back to admin-1 (`p_source`). P_stress ≥ 0.99 in ~6 % of pairs (`p_stress_ge_099`): no normal season, split unidentifiable there. Irrigated units get rainfed water-stress frequencies. |
| B | MapSPAM 2020 (Atlas SSA) mean-year yield `B_obs`; `B_norm = B_obs/(1 − E[L])` is what the formulations use. Inherits the loss table's uncertainty. |
| L | **Placeholder**, bronze maize table from `step3_adaptation.R`, read per band. AEZ-invariant. |
| m | **Scenario grid**, 10 signed profiles. No evidence-based values exist. |
| g | ERA maize Crop Yield, single-practice comparisons, study-averaged log RR. `era` is pooled across seasons (MELIAF v3.0: not usable as g_norm); `era_aez` uses the unit's AEZ16simple stratum where ≥ 3 studies; `era_norm` is g_norm from the season split (08, dry ≥ moderate). For these two practices the pooled and normal-season values are close. |
| Season | GGCMI rainfed-maize window, one per year; bimodal units flagged. |
| Price | Deferred. Scales totals, not the split. |

## Column dictionary (`af_results_maize_historic.parquet`)

| Column | Meaning |
|---|---|
| `iso3`, `admin1_name`, `admin2_name`, `gaul1_code`, `gaul2_code` | GAUL 2024 unit |
| `practice`, `hazard` | parameter set; hazards are Atlas marginals `dry` (NDWS), `wet` (NDWL0) |
| `g_source` | `grid` scenario, `era` all-SSA ERA median (pooled), `era_aez` stratum median, `era_norm` normal-season g from the ERA season split |
| `g`, `g_stratum_n` | normal-year gain used; number of ERA studies behind a stratum g (NA unless `era_aez`) |
| `m_profile`, `param_id` | buffering profile (values per severity in `params/param_grid.csv`); set id, deterministic for grid/era sets, `era_aez` ids = era ids + 210 |
| `p_source` | `adm2` or `adm1` hazard row used |
| `B_obs`, `B_norm`, `EL` | mean-year yield; normal-year yield; expected control loss share Σ P(s)·L(s) |
| `p_moderate`, `p_severe`, `p_extreme` | band probabilities (differenced tiers) |
| `p_stress` | Σ bands; `p_stress_ge_099` flags ≥ 0.99 |
| `p_none`, `p_dry_wet` | P(no hazard at all), P(dry ∩ wet), moderate tier, from exclusive classes; not used in per-hazard formulas |
| `mono_adj` | size of monotone adjustment applied to tiers (0 = none) |
| `tier_t` | production (t) behind the hazard row |
| `dn` | A − B |
| `F*_prod`, `F*_al`, `F*_total` | expected annual production, avoided loss, total, t/ha/yr |
| `F*_share` | `al / (|prod| + |al|)`, in [−1, 1]; use this |
| `F*_ashare` | `al / total`; unstable near zero totals; kept for continuity |
| `F1_minus_F3_al`, `F1_negative` | the gap; F1 < 0 flag |
| `maize_prod_t`, `maize_harv_ha`, `fs1`, `fs2`, `aez`, `aez16simple`, `bimodal`, `bimodal_share` | unit attributes from 04 |

Diagnostics add `m_eff` (Σ P·L·m / Σ P·L, NA when E[L] = 0), `ext_share` (extreme band's share of E[L]), `gap`, `rel_gap` (= gap / F3_al).

## Extending

- **New practice:** add to `config/params.yaml` (`era_prname`, `hazards`). Yaml only.
- **New crop:** not yet a config change. Maize is hard-coded in 01 (SPAM layer, file names), 03 (`Product.Simple`), 04 (column names), 05/06 (output name). `loss_frac` is keyed by hazard only. Sorghum/wheat/bean tier rows are already staged.
- **Heat:** `hazard_map` has `heat` (NTx35). Crop-specific NTxM tiers exist in the Atlas but are prototype.
