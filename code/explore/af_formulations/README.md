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

Library: `R/formulations.R` (formulations, band conversion, expected annual), `R/paths.R`. Tests: `tests/testthat/` (`Rscript -e 'testthat::test_dir("code/explore/af_formulations/tests/testthat")'`). Runtime identity checks that must pass on every run: F1 = F2 = F3 expected totals; expected F1 − F3 gap = −g·B·E[L]; E[L] in 06 equals E[L] in 05.

### Environment variables

| Var | Default | Meaning |
|---|---|---|
| `AF_OUT_DIR` | `data/outputs/af_formulations` (gitignored) | canonical outputs |
| `AF_RECORD_DIR` | unset | optional second copy (Pete uses the project `_claude/outputs/adaptation_factor`); only used if set and its parent exists |
| `AF_FS_DIR` | `../AICCRA_EFA/Data/Farming_systems` | Dixon farming-system rasters (not downloadable; ask Pete) |

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
| g | ERA maize Crop Yield, single-practice comparisons, study-averaged median log RR. **Pooled across seasons**, so an upper bound on g_norm when m > 0 (MELIAF v3.0: pooled effects cannot be used as g_norm). `era_aez` uses the unit's AEZ16simple stratum where ≥ 3 studies. |
| Season | GGCMI rainfed-maize window, one per year; bimodal units flagged. |
| Price | Deferred. Scales totals, not the split. |

## Column dictionary (`af_results_maize_historic.parquet`)

| Column | Meaning |
|---|---|
| `iso3`, `admin1_name`, `admin2_name`, `gaul1_code`, `gaul2_code` | GAUL 2024 unit |
| `practice`, `hazard` | parameter set; hazards are Atlas marginals `dry` (NDWS), `wet` (NDWL0) |
| `g_source` | `grid` scenario, `era` all-SSA ERA median, `era_aez` stratum median |
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
