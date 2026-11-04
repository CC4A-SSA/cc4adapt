# cc4adapt — Carbon Credits for Adaptation (CC4A), sub-Saharan Africa

Working repository for the CIAT side of CC4A (Gates INV-089367, with TNC). Everything here is exploratory unless a page says otherwise.

## What is here

| Part | Source | Published at |
|---|---|---|
| Practice register | `docs/practice_definitions/` (`practice_definitions.csv` is the data; `versions/` holds frozen releases) | https://cc4a-ssa.github.io/cc4adapt/practice_definitions/ |
| Adaptation-factor formulations (F0–F5) on Atlas data | `code/explore/af_formulations/` (see its `README.md`), page in `docs/af_formulations/` | https://cc4a-ssa.github.io/cc4adapt/af_formulations/ |
| Shared code | `code/common/` (config, logging, download helpers, F3 production/adaptation split, ERA yield effects), `code/inputs/fetch_era.R` | — |
| Tests | `tests/testthat/` (shared code), `code/explore/af_formulations/tests/` | — |
| References | `text/REFERENCES.md` | — |

GitHub Pages serves `main` `/docs`. Discussion threads on both pages are GitHub Discussions in this repository.

Progress on the formulation work is tracked in [issue 1](https://github.com/CC4A-SSA/cc4adapt/issues/1).

## Archived prototype pipeline

The April 2026 CCA prototype pipeline (Steps 1–5 and Activity 3 NEB, input fetchers, methodology notes and walkthrough notebook) was moved on 4 November 2026 to [CC4A-SSA/cc4adapt-pipeline-v0](https://github.com/CC4A-SSA/cc4adapt-pipeline-v0), with its git history. That repository is read-only. The last state of the pipeline in this repository is the parent of the commit that removed it.

## Licence

Data and outputs: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (`LICENSE`). Layers that incorporate TNC/Naturebase inputs follow the consortium data-sharing agreement (to be confirmed with TNC on 5 Nov). Inputs keep their own licences — see `text/REFERENCES.md`.
