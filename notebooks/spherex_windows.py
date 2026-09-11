#!/usr/bin/env python
"""When SPHEREx observed each survey comet, per phase group.

Reads the SPHEREx catalog's ``data/phase_assignment.csv`` (one row per
exposure with its phase group) and, per comet, the ``data/apphot/<T>.csv``
photometry for the exposure times and heliocentric distances, and writes one
row per (comet, phase) with the JD and r_h range covered.  The Af-rho figures
shade these windows so the ZTF trend can be read against the SPHEREx epochs.
Only the columns needed are read; the apphot tables are large.

Output: results/afrho/spherex_windows.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

import ztfcomet as zc


def main():
    root = zc.SPHEREX_ROOT
    pa_path = root / "data" / "phase_assignment.csv"
    if not pa_path.exists():
        print(f"no SPHEREx catalog at {root} (set ZTFCOMET_SPHEREX)")
        return 1
    pa = pd.read_csv(pa_path, dtype={"target": str})
    ours = {p.stem for p in (zc.result_dir("photometry")).glob("*.csv")}
    rows = []
    for target, grp in pa.groupby("target"):
        if target not in ours:
            continue
        ap_path = root / "data" / "apphot" / f"{target}.csv"
        if not ap_path.exists():
            continue
        ap = pd.read_csv(ap_path, usecols=["filename", "jd_utc", "r_hel", "r_obs"]).drop_duplicates("filename")
        m = grp.merge(ap, on="filename", how="left", suffixes=("", "_ap"))
        for phase, g in m.groupby("phase"):
            rows.append(dict(target=target, phase=int(phase), n_exposures=len(g),
                             jd_min=float(g["jd_utc"].min()), jd_max=float(g["jd_utc"].max()),
                             rh_min=float(g["r_hel"].min()), rh_max=float(g["r_hel"].max()),
                             delta_min=float(g["r_obs"].min()), delta_max=float(g["r_obs"].max()),
                             arc=int(g["arc"].iloc[0]), epoch=int(g["epoch"].iloc[0])))
    out = pd.DataFrame(rows).sort_values(["target", "phase"])
    path = zc.result_path("afrho", "spherex_windows.csv")
    out.to_csv(path, index=False)
    print(f"{len(out)} windows for {out['target'].nunique()} comets -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
