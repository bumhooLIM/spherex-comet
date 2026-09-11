#!/usr/bin/env python
"""Heliocentric dependence of Af-rho for every survey comet.

For each comet, band (r as the main series, g as a separate secondary one)
and aperture (10,000 and 20,000 km): fit ``Afrho = A r_h^-x`` on the whole
series and on each orbital leg, grade the reliability of every slope, and --
where both legs are sampled -- locate the peak of activity in days from
perihelion and refit the rising and fading phases split at that peak rather
than at perihelion.  Science lives in ``ztfcomet.activity``; this script is
selection, bookkeeping and figures.  Reads only the local photometry tables;
Horizons is queried once per comet for the perihelion time and the result is
cached in ``results/activity/elements.csv``.

Outputs
  results/activity/afrho_trends.csv     one row per comet x band x aperture x leg
  results/activity/afrho_peaks.csv      one row per comet x band x aperture (two-sided only)
  results/activity/elements.csv         cached q, e, Tp_jd
  fig/<target>/afrho_trend_<target>.png per-comet trends and peak
  fig/survey/afrho_trends_overview.png  slope and peak distributions
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ztfcomet as zc
from ztfcomet import activity as ac
from ztfcomet import orbit
from ztfcomet import rcparams  # noqa: F401
from survey import DEFAULT_LIST, read_designations

OUT = zc.RESULT_ROOT / "activity"
LEG_COLOUR = {"inbound": "tab:blue", "outbound": "tab:orange", "all": "0.4",
              "rising": "tab:green", "fading": "crimson"}


# ------------------------------------------------------------- elements
def designations():
    """slug -> designation from the survey list, so Horizons sees "2022 E2",
    not the results-directory name "2022E2"."""
    try:
        return {zc.target_slug(d): d for d in read_designations(zc.PROJECT_ROOT / DEFAULT_LIST)}
    except Exception:                                           # noqa: BLE001
        return {}


def perihelion_time(target, designation, phot, cache):
    """Tp_jd for *target*, from the cache or one Horizons elements query."""
    if target in cache.index and np.isfinite(cache.loc[target, "Tp_jd"]):
        return float(cache.loc[target, "Tp_jd"])
    jd = float(phot["obsjd"].median())
    record = ""
    try:
        record = zc.get_target(designation).resolve_orbit_record(jd)
        el = orbit.fetch_elements(record, epoch_jd=jd)
        if el is None:
            print(f"  {target}: Horizons returned no elements for {record!r}")
    except Exception as exc:                                    # noqa: BLE001
        print(f"  {target}: no elements ({str(exc)[:60]})")
        el = None
    row = dict(record=str(record), q=el.q if el else np.nan,
               e=el.e if el else np.nan, Tp_jd=el.Tp_jd if el else np.nan, epoch_jd=jd)
    cache.loc[target] = row
    return float(row["Tp_jd"])


# --------------------------------------------------------------- figures
def fig_target(target, phot, trends, peaks, tp_jd, outpath):
    two_sided = bool(trends["two_sided"].any())
    nrow = 2 if two_sided else 1
    fig, axes = plt.subplots(nrow, 2, figsize=(14, 5.2 * nrow), squeeze=False)
    for col, rho in enumerate((10000, 20000)):
        ax = axes[0, col]
        for band, mk, ms, z in (("g", "s", 4, 1), ("r", "o", 6, 2)):
            pts, _ = ac.select_points(phot, band, rho)
            if pts.empty:
                continue
            for leg, s in pts.groupby("leg"):
                c = LEG_COLOUR[leg] if band == "r" else "0.6"
                ax.errorbar(s["r"], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt=mk, ms=ms, color=c,
                            mfc=c if band == "r" else "none", alpha=0.8 if band == "r" else 0.5,
                            lw=0.8, capsize=0, zorder=z,
                            label=f"{band}  {leg}" if band == "r" else (f"g-band" if leg == "inbound" else None))
        t = trends[(trends["rho_km"] == rho) & (trends["band"] == "r") & (trends["split"] == "perihelion")]
        if t.empty:
            t = trends[(trends["rho_km"] == rho) & (trends["band"] == "r") & (trends["split"] == "all")]
        txt = []
        for _, f in t.iterrows():
            if not np.isfinite(f["x"]):
                continue
            rr = np.geomspace(f["rh_min"], f["rh_max"], 50)
            ax.plot(rr, 10 ** (f["a"] - f["x"] * np.log10(rr)), "-", color=LEG_COLOUR[f["leg"]], lw=2)
            txt.append(f"{f['leg']}: x = {f['x']:.2f} ± {f['x_err_scaled']:.2f}  n={int(f['n'])}  [{f['grade']}]")
        ax.set(xscale="log", yscale="log", xlabel=r"$r_h$ (au)", ylabel=r"$A(0°)f\rho$ (cm)",
               title=rf"{target}   $\rho$ = {rho:,} km")
        # r_h spans less than a decade for most comets: plain numbers, not 1.2x10^0
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.get_major_formatter().set_scientific(False)
        ax.xaxis.get_minor_formatter().set_scientific(False)
        ax.text(0.02, 0.02, "\n".join(txt), transform=ax.transAxes, fontsize=9, va="bottom",
                bbox=dict(fc="white", ec="0.7", alpha=0.85))
        ax.legend(fontsize=8, frameon=False, loc="upper right")
        if two_sided:
            ax2 = axes[1, col]
            pts, _ = ac.select_points(phot, "r", rho)
            tt = pts["obsjd"] - tp_jd
            for leg, s in pts.groupby("leg"):
                ax2.errorbar(tt[s.index], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt="o", ms=5,
                             color=LEG_COLOUR[leg], lw=0.8, capsize=0)
            pk = (peaks[(peaks["rho_km"] == rho) & (peaks["band"] == "r")]
                  if len(peaks) else peaks)
            ax2.axvline(0, color="crimson", ls="--", lw=1, label="perihelion")
            if len(pk) and pk.iloc[0]["bracketed"]:
                p = pk.iloc[0]
                ax2.axvspan(p["t_lo"], p["t_hi"], color="tab:green", alpha=0.15)
                ax2.axvline(p["t_peak"], color="tab:green", lw=2,
                            label=f"peak {p['t_peak']:+.0f} d  [{p['t_lo']:+.0f}, {p['t_hi']:+.0f}]")
                for _, f in trends[(trends["rho_km"] == rho) & (trends["band"] == "r")
                                   & (trends["split"] == "peak")].iterrows():
                    if np.isfinite(f["x"]):
                        ax2.plot([], [], color=LEG_COLOUR[f["leg"]], lw=2,
                                 label=f"{f['leg']}: x = {f['x']:.2f} ± {f['x_err_scaled']:.2f} [{f['grade']}]")
            elif len(pk):
                ax2.text(0.02, 0.95, "peak not bracketed", transform=ax2.transAxes, fontsize=9, va="top")
            ax2.set(yscale="log", xlabel=r"$T - T_p$ (days)", ylabel=r"$A(0°)f\rho$ (cm)")
            ax2.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=120)
    plt.close(fig)


def fig_overview(trends, peaks, outpath):
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    good = trends[trends["grade"].isin(["A", "B"]) & (trends["band"] == "r")]
    for i, rho in enumerate((10000, 20000)):
        for leg, off in (("inbound", -0.15), ("outbound", 0.15), ("all", 0.0)):
            s = good[(good["rho_km"] == rho) & (good["split"].isin(["perihelion", "all"])) & (good["leg"] == leg)]
            if s.empty:
                continue
            ax[0].errorbar(np.full(len(s), i + off) + np.random.default_rng(1).uniform(-.04, .04, len(s)),
                           s["x"], yerr=s["x_err_scaled"], fmt="o", ms=5, color=LEG_COLOUR[leg],
                           alpha=.8, lw=.8, capsize=0, label=leg if i == 0 else None)
    ax[0].set(xticks=[0, 1], xticklabels=["10,000 km", "20,000 km"], ylabel="power-law index x",
              title="r-band slopes, grades A and B")
    ax[0].axhline(0, color="0.6", lw=1); ax[0].legend(fontsize=9, frameon=False)
    for i, rho in enumerate((10000, 20000)):
        s = good[(good["rho_km"] == rho) & (good["split"] == "perihelion")]
        for leg, c in (("inbound", "tab:blue"), ("outbound", "tab:orange")):
            v = s[s["leg"] == leg]["x"].dropna()
            if len(v):
                ax[1].hist(v, bins=np.arange(-2, 8.5, 0.5), histtype="step", lw=2 if rho == 10000 else 1,
                           ls="-" if rho == 10000 else "--", color=c,
                           label=f"{leg} {rho // 1000}k  (median {v.median():.1f}, n={len(v)})")
    ax[1].set(xlabel="power-law index x", ylabel="comets", title="slope distribution by leg")
    ax[1].legend(fontsize=8, frameon=False)
    p = peaks[(peaks["band"] == "r") & peaks["bracketed"]].sort_values("t_peak")
    for i, rho in enumerate((10000, 20000)):
        q = p[p["rho_km"] == rho]
        ax[2].errorbar(q["t_peak"], np.arange(len(q)) + 0.2 * i,
                       xerr=[q["t_peak"] - q["t_lo"], q["t_hi"] - q["t_peak"]], fmt="o", ms=5,
                       color="tab:green" if rho == 10000 else "0.4", label=f"{rho // 1000}k km", lw=.8)
        if rho == 10000:
            ax[2].set_yticks(np.arange(len(q))); ax[2].set_yticklabels(q["target"], fontsize=8)
    ax[2].axvline(0, color="crimson", ls="--", lw=1)
    ax[2].set(xlabel=r"activity peak, $T - T_p$ (days)", title="bracketed peaks (r-band)")
    ax[2].legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--targets", nargs="+")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    cpath = OUT / "elements.csv"
    cache = (pd.read_csv(cpath, dtype={"target": str}).set_index("target") if cpath.exists()
             else pd.DataFrame(columns=["record", "q", "e", "Tp_jd", "epoch_jd"]).rename_axis("target"))
    files = sorted(glob.glob(str(zc.RESULT_ROOT / "*" / "photometry_*.csv")))
    if args.targets:
        want = {zc.target_slug(t) for t in args.targets}
        files = [f for f in files if Path(f).parent.name in want]

    names = designations()
    trends, peaks = [], []
    for f in files:
        target = Path(f).parent.name
        phot = pd.read_csv(f, dtype={"target": str})
        if phot.empty:
            continue
        tp = perihelion_time(target, names.get(target, target), phot, cache)
        tr, pk = ac.analyse(phot, target, tp_jd=tp)
        trends.append(tr); peaks.append(pk)
        r10 = tr[(tr["band"] == "r") & (tr["rho_km"] == 10000) & (tr["split"] != "peak")]
        print(f"{target:10s} Tp={'ok' if np.isfinite(tp) else '--'}  " +
              "  ".join(f"{row['leg']}: x={row['x']:.2f}±{row['x_err_scaled']:.2f} n={int(row['n'])} [{row['grade']}]"
                        for _, row in r10.iterrows() if np.isfinite(row["x"])) +
              ("" if pk.empty else "  peak " + ", ".join(
                  f"{p['band']}{p['rho_km'] // 1000:.0f}k:{p['t_peak']:+.0f}d{'' if p['bracketed'] else '(unbracketed)'}"
                  for _, p in pk.iterrows())), flush=True)
        if not args.no_figures and len(tr):
            fig_target(target, phot, tr, pk, tp, zc.fig_dir(target) / f"afrho_trend_{target}.png")
    cache.to_csv(cpath)
    trends = pd.concat([t for t in trends if len(t)], ignore_index=True)
    peaks = pd.concat([p for p in peaks if len(p)], ignore_index=True) if any(len(p) for p in peaks) else pd.DataFrame()
    trends.to_csv(OUT / "afrho_trends.csv", index=False)
    peaks.to_csv(OUT / "afrho_peaks.csv", index=False)
    if not args.no_figures and len(peaks):
        fig_overview(trends, peaks, zc.fig_dir("survey") / "afrho_trends_overview.png")
    print(f"\n{len(trends)} trend rows, {len(peaks)} peak rows -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
