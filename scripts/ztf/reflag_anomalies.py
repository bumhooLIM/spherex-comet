#!/usr/bin/env python
"""Apply the single-frame anomaly flag to every photometry table on disk.

``flag_anomalies`` needs only the table itself -- Af-rho and time -- so tables
written before the flag existed can be brought up to date without the FITS
files.  ``quality_ok`` and ``flags`` are recomputed.  Idempotent.
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pandas as pd

import ztfcomet as zc
from ztfcomet import phot


def main():
    total = 0
    for f in sorted(glob.glob(str(zc.result_dir("photometry") / "*.csv"))):
        t = pd.read_csv(f, dtype={"target": str})
        if t.empty:
            continue
        before = int(t["quality_ok"].sum())
        t = phot.flag_anomalies(t)
        n = int(t["flag_anomalous_bright"].sum())
        total += n
        t.to_csv(f, index=False)
        if n:
            print(f"{Path(f).stem:10s} {n:3d} anomalous frame(s); clean rows {before} -> {int(t['quality_ok'].sum())}")
    print(f"{total} anomalous frames flagged in total")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
