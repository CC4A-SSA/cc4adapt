# Biochar mitigation potential — sourcing + methodology

**Status:** authoritative reference for every script that produces or consumes a "biochar" raster in the CCA pipeline.
**Owner:** ccadapt maintainers.
**Last revised:** 2026-04-18.

This document explains, step by step, how the per-pixel biochar mitigation
potential layer for Nigeria, Kenya, and Ethiopia is derived. It is the source
of truth referenced by the headers of every modified `.R` and `.yaml` file in
the fix; if anything below changes, those headers and `data/inputs/MANIFEST.md`
must change with it.

---

## 1. Why a published biochar raster is not used

The first principle is that we never label a layer "biochar" unless it
actually represents only biochar. We surveyed the four most plausible
upstream candidates and rejected each:

1. **Sanderman et al. 2020 — "Soils Revealed soil carbon futures"**
   ([Harvard Dataverse DOI 10.7910/DVN/HA17D3](https://doi.org/10.7910/DVN/HA17D3)).
   This package publishes per-scenario delta-SOC rasters. The scenarios are
   *bundles* (`MG` = improved cropland management = cover crops + residue
   retention + biochar combined; `AF` = agroforestry; etc.). It does **not**
   publish a biochar-only raster. Using `scenario_crop_MG_dSOC_Y20.tif` as a
   biochar layer — as the previous version of this codebase did — silently
   triple-counts cover crops + residue retention as biochar. This was the
   bug that motivated the present rewrite.
2. **Woolf et al. 2010, *Nature Communications* — "Sustainable biochar to
   mitigate global climate change."**
   This is the canonical global biochar-potential study. It publishes
   **country-level totals**, not a raster. We therefore cannot use it as a
   per-pixel input, but we *do* use its country envelopes as a cross-check
   on the synthesised raster (see §6).
3. **Naturebase v3 NCS pathways "Biochar production".**
   The pathway is listed in the `all_pathways.zip` index, but at the time of
   writing the corresponding raster is empty/placeholder; the project team
   has not yet published a spatial layer for it.
4. **Karhu 2011, Lehmann & Joseph 2015, Wang et al. 2016, Smith et al. 2016.**
   These supply *parameters* (pyrolysis biochar yield, biochar carbon
   content, mean residence time, biochar suitability response curves) but
   not spatial layers.

**Conclusion:** the only defensible per-pixel biochar layer for Africa at
5-arcmin must be **synthesised from first principles**. Every term in the
formula then has a single citation.

---

## 2. The first-principles formula

For each grid cell *p* and each MapSPAM crop *c* present at *p*:

```
M_biochar(p, c)  [tCO2e ha-1 yr-1]
   =  Production(p, c)        # tonnes fresh-weight per pixel per year
                              # MapSPAM 2010 V2r0 SPAM2020_P_<crop>_A.tif
    × DM_frac(c)              # dry-matter fraction of harvested product
                              # IPCC 2019 Refinement Vol 4 Ch 11 Table 11.2
    × RPR(c)                  # residue : product ratio (dry mass)
                              # IPCC 2019 V4 T11.2 (already in build_residue_proxy.R)
    × R_available             # fraction of dry residue free for pyrolysis
                              # default 0.50 (Woolf 2010 conservative scenario)
    × Y_biochar               # kg biochar per kg dry residue (slow pyrolysis)
                              # default 0.32 (Woolf 2010)
    × C_biochar               # carbon mass fraction of biochar
                              # default 0.74 (Lehmann & Joseph 2015)
    × P_stable_100yr          # fraction of biochar-C stable over 100 yr
                              # default 0.80 (Wang et al. 2016 meta-analysis)
    × (44 / 12)               # C  ->  CO2e
    × (1 / area_ha(p))        # tonnes-per-pixel  ->  tonnes-per-ha
    × Suitability(p)          # 0..1 fuzzy soil-pH/OC/clay multiplier
                              # build_biochar_suitability.R
```

The **per-pixel biochar potential** is the sum of `M_biochar(p, c)` over all
crops *c* present at *p*. The crop loop and the residue-ratio table are
**already implemented** in `code/inputs/build_residue_proxy.R`, which writes
`data/staging/biochar/residue_biomass_tonnes_per_pixel.tif` (units:
tonnes dry residue per pixel per year). We reuse that raster directly:

```
biochar_potential_tco2e_ha_yr(p)
   =  residue_biomass_t_per_pixel(p)
    × R_available × Y_biochar × C_biochar × P_stable_100yr
    × (44 / 12)
    / cellSize_ha(p)
    × Suitability(p)
```

Because `residue_biomass_t_per_pixel` is already an annual quantity (MapSPAM
production is per year, residue ratios are dry-matter per year), **no
horizon division is applied**. This also fixes a latent unit bug in the
previous Sanderman-MG-bundle code, which divided by 20 because the source
was a 20-yr cumulative dSOC surface.

The four scalar parameters live in
[`code/inputs/data/biochar_pyrolysis_params.yaml`](../code/inputs/data/biochar_pyrolysis_params.yaml)
so they can be sensitivity-analysed by Step 5 (Monte Carlo) without code edits.

---

## 3. Inputs the synthesiser reads

| Input | Path | Produced by | Units |
|---|---|---|---|
| Crop-residue biomass | `data/staging/biochar/residue_biomass_tonnes_per_pixel.tif` | `code/inputs/build_residue_proxy.R` | tonnes dry residue / pixel / yr |
| Biochar fuzzy suitability | `data/staging/suitability/biochar_suitability_5arcmin.tif` | `code/inputs/build_biochar_suitability.R` | dimensionless 0..1 |
| Pyrolysis parameters | `code/inputs/data/biochar_pyrolysis_params.yaml` | hand-curated, citations in file | dimensionless |
| Country template grid | `data/staging/template/cca_template_3country_5arcmin.tif` | `code/inputs/build_cca_template.R` | (defines extent + resolution) |

No external HTTP download. The synthesiser is pure local computation.

---

## 4. Output

| Path | Format | Units | Purpose |
|---|---|---|---|
| `data/staging/biochar/biochar_potential_tco2e_ha_yr.tif` | COG | tCO2e / ha / yr | consumed by `code/activity2/step1_mitigation.R` as the biochar mitigation source |
| `data/staging/biochar/biochar_potential_tco2e_ha_yr.tif.method.json` | JSON sidecar | n/a | records every parameter value, the script git-SHA at run time, the cross-check ranges |

The **`.method.json`** sidecar is mandatory — the layer is synthesised, not
downloaded, so the canonical `.source.json` provenance is meaningless. The
sidecar lets a reviewer reconstruct exactly which parameter values produced
a given staged raster.

---

## 5. What happens to the old Sanderman MG bundle

The old fetcher is **not deleted** — it is moved to
`code/deprecated/fetch_biochar_potential_SANDERMAN_MG_BUNDLE.R` with a
tombstone header pointing here. Two reasons to keep it:

1. **Reproducibility of pre-fix outputs.** Anyone returning to the repo to
   reproduce the pre-2026-04-18 numbers can still rerun the deprecated
   fetcher.
2. **Upper-bound sanity check.** The MG bundle (cover crops + residue
   retention + biochar) is a useful upper bound for the *combined*
   mitigation of those three practices in cropland soils. If the new
   biochar-only layer ever exceeds the MG bundle at the same pixel, that
   is a red flag.

The downloaded raster (if it exists on disk) is moved out of the
`data/staging/biochar/` directory into
`data/staging/sanderman_mg_reference/` so no downstream consumer can pick
it up by accident.

---

## 6. Cross-check rule

After the synthesiser runs, the per-country **cropland-area-weighted mean**
of the output raster must fall inside the Woolf 2010 envelope:

| Country | Lower bound (tCO2e/ha/yr) | Upper bound | Source |
|---|---:|---:|---|
| Nigeria | 0.6 | 1.4 | Woolf 2010 SI Tab S6 (sustainable-residue scenario) |
| Kenya | 0.3 | 0.9 | Woolf 2010 SI Tab S6 |
| Ethiopia | 0.4 | 1.1 | Woolf 2010 SI Tab S6 |

Inside envelope → `INFO`. Outside envelope → `WARN`, run does **not** fail
(the envelopes are themselves uncertain; a WARN is enough to trigger
human review).

---

## 7. Postmortem of the MG-bundle bug (2026-04-18)

A colleague reviewing
[`code/inputs/fetch_biochar_potential.R`](https://github.com/AdaptationAtlas/cc4adapt/blob/main/code/inputs/fetch_biochar_potential.R)
flagged that the file it downloaded
([Harvard Dataverse fileId 4202270](https://dataverse.harvard.edu/file.xhtml?fileId=4202270&version=3.0))
"has nothing to do with biochar". Investigation confirmed:

- The script downloaded `scenario_crop_MG_dSOC_Y20.tif` — the *Management
  bundle* (cover crops + residue retention + biochar combined).
- Sanderman 2020 publishes no biochar-only raster. The script's own
  comments admitted the file was a bundle, but the resulting layer was
  still labelled `biochar` in `practices.yaml`, `step1_mitigation.R`,
  `neb_config.yaml`, and the README.
- `data/inputs/MANIFEST.md` contradicted the script outright, claiming the
  file was `scenario_crop_biochar_dSOC_Y20.tif` — a filename that does not
  exist on Dataverse.
- The `.method.json` sidecar disclosed the bundle status, but the disclosure
  did not propagate up to consumers.

**Net effect of the bug:** every downstream number that named biochar
(opportunity zones, NPV under carbon-price sweeps, country profiles) was
attributable to the wrong intervention. The fix (this document + the
synthesiser + the deprecation of the old fetcher) replaces the mislabelled
bundle with a transparent, parameter-driven biochar-only layer. Numbers
shift downward — biochar's per-hectare mitigation is genuinely lower than
the cover+residue+biochar bundle — so break-even carbon prices in
Activity 3 NEB rise. That shift is documented numerically in the
verification log produced by the synthesiser run.

**Lesson for future layers:** the layer name must match the formula, not
the source filename. A `.method.json` sidecar is necessary but not
sufficient; consumer-facing names (YAML keys, README rows, output
filenames) must also reflect what the formula actually computes.
