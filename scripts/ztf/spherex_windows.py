#!/usr/bin/env python
"""When SPHEREx observed each survey comet, per phase group.

Reads the SPHEREx stage's ``data/comspec/phase_assignment.csv`` (one row per
exposure with its phase group) and, per comet, the ``results/apphot/photometry/<T>.csv``
photometry for the exposure times and heliocentric distances, and writes one
row per (comet, phase) with the JD and r_h range covered.  The Af-rho figures
shade these windows so the ZTF trend can be read against the SPHEREx epochs.
Only the columns needed are read; the apphot tables are large.

Two further sources qualify each window.  ``results/comspec/phase_map.csv`` gives
the phase's ``arc`` label (``arc_io``): an arc *index* like the assignment's
integer -- ``in`` until a perihelion resolved inside the SPHEREx coverage,
``out`` after it -- so a comet observed only after perihelion reads ``in``
throughout; the orbital direction of an epoch is decided from T_p in
``afrho_trends.py``, not from this label.  ``results/comspec/gas_fit.csv``
gives, for the phases that were fitted, the mean r_h and JD of the channels
that carried the production rate (``rh_fit``, ``jd_fit``); they differ from
the exposure means by under 3 %, but they are the geometry the Q belongs to
and the Af-rho estimate of ``afrho_trends.py`` is evaluated there.

Output: results/ztf/afrho/spherex_windows.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd

import ztfcomet as zc

S = {"target": str}


def main():
    root = zc.SPHEREX_ROOT
    pa_path = zc.SPHEREX_PHASE_CSV
    if not pa_path.exists():
        print(f"no SPHEREx phase assignment at {pa_path} (run scripts/comspec/main.py group, or set ZTFCOMET_SPHEREX)")
        return 1
    pa = pd.read_csv(pa_path, dtype=S)
    ours = {p.stem for p in (zc.result_dir("photometry")).glob("*.csv")}
    rows = []
    for target, grp in pa.groupby("target"):
        if target not in ours:
            continue
        ap_path = zc.SPHEREX_APPHOT_DIR / f"{target}.csv"
        if not ap_path.exists():
            continue
        ap = pd.read_csv(ap_path, usecols=["filename", "jd_utc", "r_hel", "r_obs"]).drop_duplicates("filename")
        m = grp.merge(ap, on="filename", how="left", suffixes=("", "_ap"))
        for phase, g in m.groupby("phase"):
            rows.append(dict(target=target, phase=int(phase), n_exposures=len(g),
                             jd_min=float(g["jd_utc"].min()), jd_max=float(g["jd_utc"].max()),
                             jd_mean=float(g["jd_utc"].mean()),
                             rh_min=float(g["r_hel"].min()), rh_max=float(g["r_hel"].max()),
                             rh_mean=float(g["r_hel"].mean()),
                             delta_min=float(g["r_obs"].min()), delta_max=float(g["r_obs"].max()),
                             arc=int(g["arc"].iloc[0]), epoch=int(g["epoch"].iloc[0])))
    out = pd.DataFrame(rows)
    # orbital direction from the phase map; the fitted geometry from the gas table
    pm_path = zc.SPHEREX_COMSPEC_RESULT_DIR / "phase_map.csv"
    gf_path = zc.SPHEREX_COMSPEC_RESULT_DIR / "gas_fit.csv"
    if pm_path.exists():
        pm = pd.read_csv(pm_path, dtype=S, usecols=["target", "phase", "arc"]).rename(columns={"arc": "arc_io"})
        out = out.merge(pm, on=["target", "phase"], how="left")
    else:
        out["arc_io"] = ""
    if gf_path.exists():
        gf = pd.read_csv(gf_path, dtype=S, usecols=["target", "phase", "r_hel_mean", "jd_utc_mean"])
        gf = gf.rename(columns={"r_hel_mean": "rh_fit", "jd_utc_mean": "jd_fit"}).drop_duplicates(["target", "phase"])
        out = out.merge(gf, on=["target", "phase"], how="left")
    else:
        out["rh_fit"] = np.nan; out["jd_fit"] = np.nan
    out["fitted"] = np.isfinite(out["rh_fit"])
    out = out.sort_values(["target", "phase"])
    path = zc.result_path("afrho", "spherex_windows.csv")
    out.to_csv(path, index=False)
    print(f"{len(out)} windows for {out['target'].nunique()} comets "
          f"({int(out['fitted'].sum())} with a production-rate fit) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
