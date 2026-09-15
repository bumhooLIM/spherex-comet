"""
Before/after report of the 2026-09-15 case revisions (``spherex_comspec.revisions``).

Compares the current ``results/comspec/gas_fit.csv`` / ``continuum_summary.csv`` and the ZTF
``results/ztf/afrho/spherex_afrho.csv`` with a snapshot of the same tables taken before the
revisions were applied, for every (target, phase) of the memo, and prints the markdown tables
that ``doc/comspec/case_revisions.md`` carries.

    python notebooks/comspec/case_revision_report.py BEFORE_DIR [--afrho-before FILE]

``BEFORE_DIR`` holds ``gas_fit.csv``, ``continuum_summary.csv`` and ``spherex_afrho.csv`` of the
previous run (2026-09-14 state).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from spherex_comspec.revisions import CASE_REVISIONS, BAND_OF        # noqa: E402

S = {"target": str}
PRETTY = {"H2O": "H₂O", "CO2": "CO₂", "CO": "CO"}
#: new phase -> old phases, for the comet the memo regrouped
REGROUPED = {("240P", 2): [2, 3], ("240P", 3): [4]}


def q_text(row, s):
    if row is None:
        return "no fit"
    st = str(row[f"Q_{s}_status"])
    nsig = row[f"Q_{s}_nsig"]
    if st == "not_covered":
        return "not covered"
    if st == "rejected":
        return f"rejected ({nsig:+.1f}σ)"
    if st in ("upper_limit", "negative_fit"):
        return f"< {row[f'Q_{s}_upper_limit']:.1e} ({nsig:+.1f}σ)"
    tag = "" if st == "detected" else " marginal"
    return f"{row[f'Q_{s}']:.2e} ± {row[f'Q_{s}_err']:.1e} ({nsig:.1f}σ, n_eff {row[f'Q_{s}_n_eff']:.1f}){tag}"


def verdict_text(c, t, ph, band):
    r = c[(c.target == t) & (c.phase == ph) & (c.band == band)]
    if r.empty:
        return "—"
    x = r.iloc[0]
    o = f"o{x.poly_order_used:.0f}" if np.isfinite(x.poly_order_used) else "—"
    return f"{x.verdict} ({o}, {int(x.n_emission)} em, {x.cont_lo_um:.2f}–{x.cont_hi_um:.2f})"


def afrho_text(a, t, ph, rho=10000):
    r = a[(a.target == t) & (a.phase == ph) & (a.rho_km == rho)]
    if r.empty:
        return "—"
    x = r.iloc[0]
    if not np.isfinite(x.afrho_cm):
        return f"none ({str(x.note)[:60]})"
    return f"{x.afrho_cm:.0f} ± {x.afrho_err_cm:.0f} cm ({x.method}, {x.leg})"


def main(argv=None):
    argv = list(argv if argv is not None else sys.argv[1:])
    before = Path(argv[0])
    f0, f1 = pd.read_csv(before / "gas_fit.csv", dtype=S), pd.read_csv(ROOT / "results/comspec/gas_fit.csv", dtype=S)
    c0, c1 = (pd.read_csv(before / "continuum_summary.csv", dtype=S),
              pd.read_csv(ROOT / "results/comspec/continuum_summary.csv", dtype=S))
    a0 = pd.read_csv(before / "spherex_afrho.csv", dtype=S)
    a1 = pd.read_csv(ROOT / "results/ztf/afrho/spherex_afrho.csv", dtype=S)
    row = lambda f, t, ph: (f[(f.target == t) & (f.phase == ph)].iloc[0]
                            if len(f[(f.target == t) & (f.phase == ph)]) else None)
    print("| comet | phase | directive (memo) | band verdicts before → after | Q before → after |")
    print("|---|---|---|---|---|")
    for r in CASE_REVISIONS:
        t, ph = r.target, r.phase
        olds = REGROUPED.get((t, ph), [ph])
        directive = r.describe() or r.memo
        verd = []
        for b in r.bands:
            verd.append(f"{b}: " + " / ".join(verdict_text(c0, t, o, b) for o in olds) + " → " + verdict_text(c1, t, ph, b))
        qs = []
        for s in ("H2O", "CO2", "CO"):
            b0 = " / ".join(q_text(row(f0, t, o), s) for o in olds)
            b1 = q_text(row(f1, t, ph), s)
            if b0 != b1:
                qs.append(f"{PRETTY[s]}: {b0} → {b1}")
        if r.memo.startswith("ztf_afrho"):
            qs.append("Afρ(10k): " + " / ".join(afrho_text(a0, t, o) for o in olds) + " → " + afrho_text(a1, t, ph))
        print(f"| {t} | {ph} | {directive} | {'; '.join(verd) or '—'} | {'; '.join(qs) or 'unchanged'} |")
    # census
    def census(f):
        out = {}
        for s in ("H2O", "CO2", "CO"):
            det = f[(f[f"Q_{s}_status"] == "detected") & (f[f"Q_{s}_n_eff"] >= 2)]
            out[s] = (len(det), int((f[f"Q_{s}_status"] == "marginal").sum()), int((f[f"Q_{s}_status"] == "rejected").sum()))
        return out
    print("\ncensus (robust >= 3 sigma n_eff >= 2, marginal, rejected):", "before", census(f0), "after", census(f1))
    print("fits before/after:", len(f0), len(f1), "| chi2 median before/after: %.2f %.2f" % (f0.chi2_red.median(), f1.chi2_red.median()))


if __name__ == "__main__":
    main()
