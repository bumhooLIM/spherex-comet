"""
Attach the ZTF dust context to the catalog summaries.

The ZTF stage (``scripts/ztf/afrho_trends.py``, package ``ztfcomet``) estimates, for every SPHEREx
phase group of a comet, the r-band A(0°)fρ at the group's mean r_h from the ZTF series:
the frames inside the window when ZTF observed then (``direct``), otherwise the fitted
heliocentric law of the orbital phase the epoch falls in (``trend``), otherwise a bounded
extrapolation of that law (``trend_extrap``), otherwise nothing (``none``) with the reason --
at ρ = 10 000 and 20 000 km.  Its table ``results/ztf/afrho/spherex_afrho.csv`` has one row per
comet, phase and aperture.

This script pivots that table to one row per (target, phase) -- ``results/comspec/afrho_ztf.csv`` --
and attaches its summary columns to ``results/comspec/gas_fit.csv`` and ``results/comspec/phase_map.csv``
through :func:`spherex_comspec.dataio.attach_afrho_ztf`, which the pipeline itself calls
whenever it rewrites those tables, so the columns survive a rerun as long as
``afrho_ztf.csv`` is current.  Columns attached:

    afrho_rh_au                               r_h (au) the estimate was evaluated at: the gas fit's
                                              r_hel_mean where the phase was fitted, else the exposure mean
    afrho_10k_cm, afrho_10k_err_cm, afrho_10k_method      ρ = 10 000 km
    afrho_20k_cm, afrho_20k_err_cm, afrho_20k_method      ρ = 20 000 km
    afrho_note                                "10k: ...; 20k: ..." -- how each value was obtained, or why none

The 1σ ranges, frame counts, grades and legs stay in ``afrho_ztf.csv``.

Usage
-----
    python scripts/comspec/attach_afrho_ztf.py                    # results/ztf/afrho/spherex_afrho.csv
    python scripts/comspec/attach_afrho_ztf.py --ztf-afrho PATH   # another spherex_afrho.csv
    COMSPEC_ZTF_AFRHO=PATH python scripts/comspec/attach_afrho_ztf.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from spherex_comspec import directory as _dir  # noqa: E402
from spherex_comspec.dataio import AFRHO_ATTACH_COLUMNS, attach_afrho_ztf  # noqa: E402

S = {"target": str}
APERTURES = ((10000, "10k"), (20000, "20k"))
PER_APERTURE = ["afrho_cm", "afrho_err_cm", "afrho_lo_cm", "afrho_hi_cm", "method", "leg", "n", "grade", "note"]


def widen(long: pd.DataFrame) -> pd.DataFrame:
    """One row per (target, phase) with both apertures side by side."""
    keys = ["target", "phase"]
    head = (long.drop_duplicates(keys)[keys + ["fitted", "rh", "rh_source", "jd", "jd_min", "jd_max"]]
            .rename(columns={"rh": "afrho_rh_au", "jd": "afrho_jd"}))
    wide = head
    for rho, tag in APERTURES:
        s = long[np.isclose(long["rho_km"], rho)][keys + PER_APERTURE]
        s = s.rename(columns={c: f"afrho_{tag}_{c.replace('afrho_', '')}" for c in PER_APERTURE})
        wide = wide.merge(s, on=keys, how="left")
    note = []
    for _, r in wide.iterrows():
        note.append("; ".join(f"{tag}: {r.get(f'afrho_{tag}_note', '')}" for _, tag in APERTURES))
    wide["afrho_note"] = note
    wide["phase"] = wide["phase"].astype(int)
    return wide.sort_values(keys).reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--ztf-afrho", type=Path, default=_dir.ZTF_AFRHO_SRC,
                    help=f"the spherex_afrho.csv of the ZTF stage (default {_dir.ZTF_AFRHO_SRC})")
    args = ap.parse_args(argv)
    src = args.ztf_afrho
    if not src.is_file():
        print(f"no {src}: run scripts/ztf/afrho_trends.py first (or pass --ztf-afrho)")
        return 1
    long = pd.read_csv(src, dtype=S)
    wide = widen(long)
    _dir.AFRHO_ZTF_CSV.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(_dir.AFRHO_ZTF_CSV, index=False)
    n_val = {tag: int(np.isfinite(wide[f"afrho_{tag}_cm"]).sum()) for _, tag in APERTURES}
    print(f"{len(wide)} phases of {wide['target'].nunique()} comets -> {_dir.AFRHO_ZTF_CSV}  "
          f"(values at 10k: {n_val['10k']}, 20k: {n_val['20k']})")
    for name, jd_col, mid in (("gas_fit.csv", "jd_utc_mean", None), ("phase_map.csv", "jd_mid", ("jd_start", "jd_end"))):
        p = _dir.RESULT_DIR / name
        if not p.is_file():
            print(f"  {name}: not found, skipped")
            continue
        df = pd.read_csv(p, dtype=S)
        work = df.assign(jd_mid=0.5 * (df[mid[0]] + df[mid[1]])) if mid else df
        out = attach_afrho_ztf(work, wide, jd_col=jd_col)
        if mid:
            out = out.drop(columns=["jd_mid"])
        out.to_csv(p, index=False)
        has = {tag: int(np.isfinite(out[f"afrho_{tag}_cm"]).sum()) for _, tag in APERTURES}
        stale = int(out["afrho_note"].fillna("").str.startswith("stale").sum())
        print(f"  {name}: {len(out)} rows, Afrho at 10k for {has['10k']}, at 20k for {has['20k']}"
              + (f", {stale} stale" if stale else "") + f"  [{', '.join(AFRHO_ATTACH_COLUMNS)}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
