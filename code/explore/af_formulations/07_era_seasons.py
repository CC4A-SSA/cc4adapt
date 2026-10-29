"""
code/explore/af_formulations/07_era_seasons.py

EXPLORATORY. Classify ERA maize site-seasons as stress or normal, per hazard,
using observed daily climate (CHIRPS v3 precipitation + AgERA5 temperature via
Earth Engine) and the Atlas severity bands, so g_norm and g_stress can be
estimated separately (science review findings 2 and 4, 29 Oct 2026).

Run from the cc4adapt repo root, with the climate-toolkit checkout importable:
    GCP_PROJECT_ID=ee-peetmate PYTHONPATH=../climate-toolkit \
      python3 code/explore/af_formulations/07_era_seasons.py [--limit N] [--workers W]

Input:  data/staging/af_formulations/era_site_years.csv  (from 03b, R)
Output: data/staging/af_formulations/era_site_season_hazards.csv
        one row per site x season slot x year with season window, day counts,
        Atlas-style monthly means and severity bands.

Method
  * Sites are grouped by (Site.Key, season slot). The slot is ERA's Season.Start
    (1 or 2) where set, else 1. The window comes from era_site_windows.csv
    (03c: ERA planting/harvest dates, else WorldClim rainfall climatology);
    the GGCMI calendar is the fallback when 03c has no window.
  * One evaluate_hazards() call per group covers all its years (fetch is per
    site; chunks are cached under climate-toolkit/outputs/cache).
  * The toolkit counts NDWS (ERATIO < 0.5), NDWL0 (logging > 0) and NTx35 over
    the whole window. The Atlas classifies the mean monthly count over the
    season months, so monthly_mean = count * 30.44 / counted_days, and bands use
    the Atlas generic thresholds: NDWS 15/20/25, NDWL0 2/5/8, NTx35 7/14/21.
  * Toolkit thresholds for NTx35 (10/20/25) are not used.

Everything here is a first pass; the water balance uses the toolkit's FAO-56
starter crop parameters and a 1 m root zone.
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
IN = ROOT / "data/staging/af_formulations/era_site_years.csv"
WIN = ROOT / "data/staging/af_formulations/era_site_windows.csv"
OUT = ROOT / "data/staging/af_formulations/era_site_season_hazards.csv"
LOG = ROOT / "data/staging/af_formulations/era_site_season_hazards.log.jsonl"
MIN_YEAR = 1981
DAYS_PER_MONTH = 365.25 / 12
BANDS = {"NDWS": (15, 20, 25), "NDWL0": (2, 5, 8), "NTx35": (7, 14, 21)}


def band(v, cuts):
    if v is None or pd.isna(v):
        return None
    m, s, e = cuts
    return "extreme" if v > e else "severe" if v > s else "moderate" if v > m else "none"


WINDOWS = {}


def season_token(lat, lon, slot, site=None):
    w = WINDOWS.get((site, int(slot)))
    if w and isinstance(w.get("fixed_season"), str):
        return w["fixed_season"], f"{w['source']}:{w['window_status']}"
    from climate_toolkit.crop_calendar.ggcmi import extract_point_calendar

    rows = extract_point_calendar(lat, lon, "maize", system="rf")
    rows = rows[rows["fixed_season_token"].notna()]
    if rows.empty:
        return None, "no_calendar"
    have = set(int(x) for x in rows["season_id"])
    if slot in have:
        r = rows[rows["season_id"] == slot].iloc[0]
        return r["fixed_season_token"], "calendar_slot_%d" % slot
    r = rows.sort_values("season_id").iloc[0]
    return r["fixed_season_token"], ("calendar_slot_missing" if slot != int(r["season_id"]) else "calendar_slot_1")


def run_group(key, g, project):
    import climate_toolkit as ct

    site, slot = key[0], int(key[1])
    lat, lon = float(g["lat"].iloc[0]), float(g["lon"].iloc[0])
    years = sorted(int(y) for y in g["yr"].dropna().unique() if y >= MIN_YEAR)
    if not years:
        return [], {"site": site, "slot": slot, "status": "no_years"}
    token, cal_status = season_token(lat, lon, slot, site)
    t0 = time.time()
    try:
        r = ct.evaluate_hazards(
            # Fetch through mid next year so a year-crossing window (e.g. Oct-Jan)
            # of the last year is complete; rows are filtered back to `years`.
            "maize", (lat, lon), f"{years[0]}-01-01", f"{years[-1] + 1}-06-30",
            source="paired", precip_source="chirps_v3_daily_rnl", temp_source="agera_5",
            fixed_season=token, calendar_source=None if token else "ggcmi",
            ee_project_id=project,
        )
    except Exception as e:  # noqa: BLE001
        return [], {"site": site, "slot": slot, "status": "exception", "error": str(e)[:300]}
    if "error" in r:
        return [], {"site": site, "slot": slot, "status": "error", "error": str(r["error"])[:300], "secs": round(time.time() - t0)}
    win = (r.get("baseline_ltm") or {}).get("water_balance_methodology", {}).get("analysis_window", {})
    rows = []
    for c in r.get("baseline_ltm_comparisons", []):
        if c.get("year") not in years:
            continue
        met = c.get("metrics", {})
        counted = None
        for ps in (r.get("baseline_ltm") or {}).get("per_season", []):
            if ps.get("season_number") == c.get("season_number"):
                counted = ps.get("mean_length_days")
        row = {
            "Site.Key": site, "season_slot": slot, "lat": lat, "lon": lon, "year": c.get("year"),
            "season_number": c.get("season_number"), "calendar_status": cal_status, "fixed_season": token,
            "onset": win.get("onset_date"), "cessation": win.get("cessation_date"), "counted_days": counted,
            "total_mm": met.get("total_mm", {}).get("actual"), "mean_tavg": met.get("mean_tavg", {}).get("actual"),
            "wrsi": met.get("wrsi", {}).get("actual"),
        }
        for h in BANDS:
            n = met.get(h, {}).get("actual")
            mm = (n * DAYS_PER_MONTH / counted) if (n is not None and counted) else None
            row[f"{h}_days"] = n
            row[f"{h}_monthly_mean"] = round(mm, 2) if mm is not None else None
            row[f"{h}_band"] = band(mm, BANDS[h])
        rows.append(row)
    return rows, {"site": site, "slot": slot, "status": "ok", "years": len(rows), "secs": round(time.time() - t0), "calendar": cal_status}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="run only the first N site-slot groups (dry run)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--project", default=os.environ.get("GCP_PROJECT_ID", ""))
    ap.add_argument("--window", choices=["site", "ggcmi"], default="site",
                    help="site: 03c windows (ERA dates / climatology); ggcmi: GGCMI rainfed-maize window as the Atlas uses")
    a = ap.parse_args()
    global OUT
    if a.window == "ggcmi":
        OUT = OUT.with_name("era_site_season_hazards_ggcmi.csv")
    if not a.project:
        sys.exit("Set GCP_PROJECT_ID or --project")

    sy = pd.read_csv(IN)
    if WIN.exists() and a.window == "site":
        for r in pd.read_csv(WIN).to_dict("records"):
            WINDOWS[(r["Site.Key"], int(r["slot"]))] = r
        print(f"windows loaded for {len(WINDOWS)} site-slots", flush=True)
    sy["slot"] = sy["season_idx"].fillna(1).astype(int)
    sy = sy[sy["yr"] >= MIN_YEAR]
    groups = list(sy.groupby(["Site.Key", "slot"]))
    if a.window == "ggcmi":
        # one GGCMI window per site: collapse slots so the Atlas-style single season is scored once
        sy["slot"] = 1
        groups = list(sy.groupby(["Site.Key", "slot"]))
    if a.limit:
        groups = groups[: a.limit]
    print(f"{len(groups)} site-slot groups, {sy['yr'].nunique()} distinct years, workers={a.workers}", flush=True)

    out = OUT.with_suffix(".dryrun.csv") if a.limit else OUT
    done = set()
    if out.exists() and not a.limit:
        prev = pd.read_csv(out)
        done = set(zip(prev["Site.Key"], prev["season_slot"]))
        print(f"resuming: {len(done)} groups already in {out.name}", flush=True)
    elif a.limit and out.exists():
        out.unlink()
    todo = [(k, g) for k, g in groups if k not in done]

    with ThreadPoolExecutor(max_workers=a.workers) as ex, open(LOG, "a") as log:
        futs = {ex.submit(run_group, k, g, a.project): k for k, g in todo}
        for i, f in enumerate(as_completed(futs), 1):
            rows, status = f.result()
            log.write(json.dumps(status, default=str) + "\n"); log.flush()
            print(f"[{i}/{len(todo)}] {status}", flush=True)
            if rows:
                pd.DataFrame(rows).to_csv(out, mode="a", header=not out.exists(), index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
