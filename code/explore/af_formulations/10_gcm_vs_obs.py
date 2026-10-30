"""
code/explore/af_formulations/10_gcm_vs_obs.py

EXPLORATORY. Why is the Atlas P(dry >= moderate) so much higher than the
observed share of stress seasons at Kenyan ERA sites, even on the Atlas's own
season window? This scores the same sites, window and period (1995-2014, the
Atlas baseline) on (a) observed climate (CHIRPS v3 + AgERA5) and (b) NEX-GDDP
historical runs of Atlas ensemble members, through the same toolkit water
balance. If (b) is near the Atlas value and (a) is not, the gap is the GCM
climate; if both are low, it is the Atlas water balance or thresholds.

Run from the cc4adapt repo root:
    GCP_PROJECT_ID=ee-peetmate PYTHONPATH=../climate-toolkit \
      python3 code/explore/af_formulations/10_gcm_vs_obs.py [--iso KEN] [--models ACCESS-CM2,GFDL-ESM4] [--limit N]

Output: data/staging/af_formulations/gcm_vs_obs_sites.csv
        one row per site x source x year (NDWS, NDWL0, NTx35 counts and bands)

The toolkit's calculate_hazards() does not take model/scenario, but its fetch
helper does; the wrapper below injects them when source == "nex_gddp".
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SITES = ROOT / "data/outputs/af_formulations/era_split/obs_vs_atlas_sites.csv"
OUT = ROOT / "data/staging/af_formulations/gcm_vs_obs_sites.csv"
LOG = ROOT / "data/staging/af_formulations/gcm_vs_obs.log.jsonl"
Y0, Y1 = 1995, 2014
DAYS_PER_MONTH = 365.25 / 12
BANDS = {"NDWS": (15, 20, 25), "NDWL0": (2, 5, 8), "NTx35": (7, 14, 21)}
_MODEL = {"model": None}


def band(v, cuts):
    if v is None or pd.isna(v):
        return None
    m, s, e = cuts
    return "extreme" if v > e else "severe" if v > s else "moderate" if v > m else "none"


def install_nex_wrapper():
    from climate_toolkit.calculate_hazards import hazards as H

    orig = H.get_climate_data_for_season

    def wrapped(*a, **k):
        if k.get("source") == "nex_gddp":
            k.setdefault("model", _MODEL["model"])
            k.setdefault("scenario", "historical")
        return orig(*a, **k)

    H.get_climate_data_for_season = wrapped


def ggcmi_token(lat, lon):
    from climate_toolkit.crop_calendar.ggcmi import extract_point_calendar

    rows = extract_point_calendar(lat, lon, "maize", system="rf")
    rows = rows[rows["fixed_season_token"].notna()].sort_values("season_id")
    return None if rows.empty else rows.iloc[0]["fixed_season_token"]


def score(site, lat, lon, source, model, project):
    import climate_toolkit as ct

    token = ggcmi_token(lat, lon)
    t0 = time.time()
    kw = dict(source=source, fixed_season=token, calendar_source=None if token else "ggcmi", ee_project_id=project)
    if source == "paired":
        kw.update(precip_source="chirps_v3_daily_rnl", temp_source="agera_5")
    _MODEL["model"] = model
    try:
        r = ct.evaluate_hazards("maize", (lat, lon), f"{Y0}-01-01", f"{Y1}-12-31", **kw)
    except Exception as e:  # noqa: BLE001
        return [], {"site": site, "source": source, "model": model, "status": "exception", "error": str(e)[:300]}
    if "error" in r:
        return [], {"site": site, "source": source, "model": model, "status": "error", "error": str(r["error"])[:300]}
    counted = {ps.get("season_number"): ps.get("mean_length_days") for ps in (r.get("baseline_ltm") or {}).get("per_season", [])}
    rows = []
    for c in r.get("baseline_ltm_comparisons", []):
        met = c.get("metrics", {})
        n_days = counted.get(c.get("season_number"))
        row = {"Site.Key": site, "lat": lat, "lon": lon, "source": source, "model": model or "", "year": c.get("year"),
               "fixed_season": token, "counted_days": n_days, "total_mm": met.get("total_mm", {}).get("actual"),
               "wrsi": met.get("wrsi", {}).get("actual")}
        for h in BANDS:
            n = met.get(h, {}).get("actual")
            mm = (n * DAYS_PER_MONTH / n_days) if (n is not None and n_days) else None
            row[f"{h}_days"] = n
            row[f"{h}_monthly_mean"] = round(mm, 2) if mm is not None else None
            row[f"{h}_band"] = band(mm, BANDS[h])
        rows.append(row)
    return rows, {"site": site, "source": source, "model": model, "status": "ok", "years": len(rows), "secs": round(time.time() - t0)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso", default="KEN")
    ap.add_argument("--models", default="ACCESS-CM2,GFDL-ESM4")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--project", default=os.environ.get("GCP_PROJECT_ID", ""))
    a = ap.parse_args()
    if not a.project:
        sys.exit("Set GCP_PROJECT_ID or --project")
    install_nex_wrapper()

    sites = pd.read_csv(SITES)
    sites = sites[sites["iso3"] == a.iso][["Site.Key", "lat", "lon"]].drop_duplicates()
    if a.limit:
        sites = sites.head(a.limit)
    jobs = [(s, "paired", None) for s in sites.itertuples(index=False)]
    jobs += [(s, "nex_gddp", m) for m in a.models.split(",") for s in sites.itertuples(index=False)]

    out = OUT.with_suffix(".dryrun.csv") if a.limit else OUT
    done = set()
    if out.exists() and not a.limit:
        prev = pd.read_csv(out)
        done = set(zip(prev["Site.Key"], prev["source"], prev["model"].fillna("")))
    elif a.limit and out.exists():
        out.unlink()
    todo = [j for j in jobs if (j[0][0], j[1], j[2] or "") not in done]
    print(f"{len(sites)} sites, {len(jobs)} jobs, {len(todo)} to do, workers={a.workers}", flush=True)

    with ThreadPoolExecutor(max_workers=a.workers) as ex, open(LOG, "a") as log:
        futs = [ex.submit(score, s[0], float(s[1]), float(s[2]), src, m, a.project) for s, src, m in todo]
        for i, f in enumerate(as_completed(futs), 1):
            rows, status = f.result()
            log.write(json.dumps(status, default=str) + "\n"); log.flush()
            print(f"[{i}/{len(todo)}] {status}", flush=True)
            if rows:
                pd.DataFrame(rows).to_csv(out, mode="a", header=not out.exists(), index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
