#!/usr/bin/env python
"""Heliocentric dependence of Af-rho for every survey comet.

For each comet, band (r as the main series, g as a separate secondary one)
and aperture (10,000 and 20,000 km): fit ``Afrho = A r_h^-x`` on each phase of
activity, grade the reliability of every slope, test whether one power law
serves the whole r_h range or a broken one is needed (free break, and a
break imposed at 3 au), and follow the g - r colour of the dust with r_h.

Phases are defined by the **peak of activity, not perihelion**: a comet
sampled on both sides of T_p gets a smoothed peak, and its series splits
into *rising* and *fading* there; when the maximum sits at an edge of the
coverage the whole series is one phase.  A comet sampled on one side only is
one leg named by its orbital direction.  Science lives in
``ztfcomet.activity``; this script is selection, bookkeeping and figures.
Horizons is queried once per comet for T_p and cached.

Outputs
  results/afrho/trends.csv       one row per comet x band x aperture x leg
  results/afrho/peaks.csv        two-sided comets: the activity peak
  results/afrho/breaks.csv       broken-law tests per leg (free break and 3 au)
  results/afrho/colour.csv       g - r excess vs r_h, with the step change point
  results/afrho/elements.csv     cached q, e, Tp_jd
  fig/afrho/<target>_trend.png   per-comet: fits, peak, break, colour
  fig/afrho/survey_overview.png  slope, break and colour distributions
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

OUT = zc.result_dir("afrho")
LEG_COLOUR = {"inbound": "tab:blue", "outbound": "tab:orange", "all": "0.4",
              "rising": "tab:green", "fading": "crimson"}
RHOS = (10000, 20000)


# ------------------------------------------------------------- elements
def designations():
    """slug -> designation from the survey list, so Horizons sees "2022 E2",
    not the results-file name "2022E2"."""
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
def _plain_log_axis(ax):
    for f in (ax.xaxis.set_major_formatter, ax.xaxis.set_minor_formatter):
        f(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_minor_formatter().set_scientific(False)


def _label(f):
    return f"{f['leg']}: x = {f['x']:.2f} ± {f['x_err_scaled']:.2f}  n={int(f['n'])}  [{f['grade']}]"


def fig_target(target, phot, trends, peaks, breaks, colours, tp_jd, outpath):
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))
    two_sided = bool(trends["two_sided"].any())
    for col, rho in enumerate(RHOS):
        tr = trends[(trends["rho_km"] == rho) & (trends["band"] == "r")]
        prim = tr[tr["primary"]]
        # ---- row 0: Afrho vs r_h, primary legs
        ax = axes[0, col]
        pts, _ = ac.select_points(phot, "r", rho)
        gpts, _ = ac.select_points(phot, "g", rho)
        if not gpts.empty:
            ax.errorbar(gpts["r"], gpts["afrho0_cm"], yerr=gpts["afrho0_cm_err"], fmt="s", ms=3.5, mfc="none",
                        color="0.6", alpha=.5, lw=.6, capsize=0, zorder=1, label="g-band")
        if not pts.empty:
            leg_of = pd.Series("all", index=pts.index)
            if tp_jd is not None and np.isfinite(tp_jd) and len(prim) and (prim["split"] == "peak").all():
                t = pts["obsjd"] - tp_jd
                tpk = float(prim["t_peak"].iloc[0])
                if bool(prim["peak_bracketed"].iloc[0]):
                    leg_of = pd.Series(np.where(t < tpk, "rising", "fading"), index=pts.index)
                else:
                    leg_of = pd.Series(prim["leg"].iloc[0], index=pts.index)
            elif len(prim):
                leg_of = pts["leg"]
            for leg, s in pts.groupby(leg_of):
                c = LEG_COLOUR.get(leg, "0.4")
                ax.errorbar(s["r"], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt="o", ms=5.5, color=c,
                            lw=.8, capsize=0, zorder=2, label=f"r  {leg}")
        txt = []
        for _, f in prim.iterrows():
            if np.isfinite(f["x"]):
                rr = np.geomspace(f["rh_min"], f["rh_max"], 50)
                ax.plot(rr, 10 ** (f["a"] - f["x"] * np.log10(rr)), "-", color=LEG_COLOUR[f["leg"]], lw=2)
                txt.append(_label(f))
        bk = (breaks[(breaks["rho_km"] == rho) & (breaks["band"] == "r") & breaks["primary"]]
              if len(breaks) else breaks)
        for _, b in bk.iterrows():
            if b["preferred"]:
                xb = np.log10(b["r_break"])
                lr = np.log10(np.geomspace(b["rh_min"], b["rh_max"], 60))
                sub = pts[pts["r"].between(b["rh_min"], b["rh_max"])]
                # anchor the broken curve on the data at the break
                y0 = np.interp(xb, np.log10(sub["r"].sort_values()), np.log10(sub.sort_values("r")["afrho0_cm"]))
                yy = y0 - b["x_inner"] * (lr - xb) - (b["x_outer"] - b["x_inner"]) * np.maximum(lr - xb, 0)
                ax.plot(10 ** lr, 10 ** yy, "--", color="k", lw=1.5, zorder=3)
                ax.axvline(b["r_break"], color="k", ls=":", lw=1)
                txt.append(f"{b['leg']}: break at {b['r_break']:.2f} au  [{b['r_break_lo']:.2f}, {b['r_break_hi']:.2f}]  "
                           f"x = {b['x_inner']:.1f} → {b['x_outer']:.1f}  (ΔBIC {b['dbic']:.0f})")
            elif b["testable"]:
                txt.append(f"{b['leg']}: single law adequate (ΔBIC {b['dbic']:.0f})")
            if b["at3_testable"]:
                txt.append(f"{b['leg']} at 3 au: x(<3) = {b['x_lt3']:.2f} ± {b['x_lt3_err']:.2f},  "
                           f"x(>3) = {b['x_gt3']:.2f} ± {b['x_gt3_err']:.2f}")
        if not pts.empty and pts["r"].min() < 3.0 < pts["r"].max():
            ax.axvline(3.0, color="0.7", ls="-.", lw=1)
        ax.set(xscale="log", yscale="log", xlabel=r"$r_h$ (au)", ylabel=r"$A(0°)f\rho$ (cm)",
               title=rf"{target}   $\rho$ = {rho:,} km")
        _plain_log_axis(ax)
        ax.text(0.02, 0.02, "\n".join(txt) or "no fit", transform=ax.transAxes, fontsize=8.5, va="bottom",
                bbox=dict(fc="white", ec="0.7", alpha=0.85))
        ax.legend(fontsize=8, frameon=False, loc="upper right")
        # ---- row 1: time from perihelion
        ax2 = axes[1, col]
        if two_sided and tp_jd is not None and np.isfinite(tp_jd) and not pts.empty:
            tt = pts["obsjd"] - tp_jd
            for leg, s in pts.groupby(leg_of):
                ax2.errorbar(tt[s.index], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt="o", ms=5,
                             color=LEG_COLOUR.get(leg, "0.4"), lw=.8, capsize=0)
            ax2.axvline(0, color="crimson", ls="--", lw=1, label="perihelion")
            pk = peaks[(peaks["rho_km"] == rho) & (peaks["band"] == "r")] if len(peaks) else peaks
            if len(pk):
                p = pk.iloc[0]
                if p["bracketed"]:
                    ax2.axvspan(p["t_lo"], p["t_hi"], color="tab:green", alpha=0.15)
                    ax2.axvline(p["t_peak"], color="tab:green", lw=2,
                                label=f"peak {p['t_peak']:+.0f} d  [{p['t_lo']:+.0f}, {p['t_hi']:+.0f}]")
                else:
                    ax2.text(0.02, 0.95, f"maximum at the edge of coverage ({p['t_peak']:+.0f} d): "
                             f"{'still rising' if p['t_peak'] > tt.median() else 'already fading'}",
                             transform=ax2.transAxes, fontsize=9, va="top")
            ax2.set(yscale="log", xlabel=r"$T - T_p$ (days)", ylabel=r"$A(0°)f\rho$ (cm)")
            ax2.legend(fontsize=8, frameon=False, loc="lower right")
        else:
            ax2.text(0.5, 0.5, "one-sided coverage: peak not measurable", ha="center", va="center",
                     transform=ax2.transAxes, fontsize=10, color="0.4")
            ax2.set_axis_off()
        # ---- row 2: colour
        ax3 = axes[2, col]
        pairs = ac.colour_pairs(phot, rho)
        cl = colours[(colours["rho_km"] == rho) & (colours["leg"] == "all")] if len(colours) else colours
        if len(pairs) >= 3:
            ax3.errorbar(pairs["r"], pairs["excess"], yerr=pairs["excess_err"], fmt="o", ms=5, color="tab:purple",
                         lw=.8, capsize=0, label=rf"{len(pairs)} g/r pairs (Δt < 1 d): $-2.5\log(Af\rho_g/Af\rho_r)$")
            ax3.axhline(0, color="0.6", lw=1)
            if len(cl):
                c = cl.iloc[0]
                if c["change_preferred"]:
                    ax3.axvspan(c["r_change_lo"], c["r_change_hi"], color="tab:purple", alpha=.12)
                    ax3.axvline(c["r_change"], color="tab:purple", lw=2)
                    ax3.hlines(c["colour_below"], pairs["r"].min(), c["r_change"], color="k", lw=2)
                    ax3.hlines(c["colour_above"], c["r_change"], pairs["r"].max(), color="k", lw=2)
                    lab = (f"colour changes at {c['r_change']:.2f} au [{c['r_change_lo']:.2f}, {c['r_change_hi']:.2f}]: "
                           f"Δ = {c['delta_colour']:+.3f} ± {c['delta_colour_err']:.3f} mag  (ΔBIC {c['dbic']:.0f})")
                else:
                    lab = f"no colour change preferred (ΔBIC {c['dbic']:.0f}); slope {c['slope']:+.2f} ± {c['slope_err']:.2f} mag/dex"
                ax3.text(0.02, 0.97, lab, transform=ax3.transAxes, fontsize=8.5, va="top",
                         bbox=dict(fc="white", ec="0.7", alpha=0.85))
            ax3.set(xscale="log", xlabel=r"$r_h$ (au)", ylabel=r"colour excess over solar (mag)")
            _plain_log_axis(ax3)
            ax3.legend(fontsize=8, frameon=False, loc="lower right")
        else:
            ax3.text(0.5, 0.5, "fewer than 3 same-night g/r pairs", ha="center", va="center",
                     transform=ax3.transAxes, fontsize=10, color="0.4")
            ax3.set_axis_off()
    fig.tight_layout()
    fig.savefig(outpath, dpi=110)
    plt.close(fig)


def fig_overview(trends, peaks, breaks, colours, outpath):
    fig, ax = plt.subplots(1, 3, figsize=(19, 5.6))
    pr = trends[trends["primary"] & (trends["band"] == "r") & trends["grade"].isin(["A", "B"])]
    order = ["rising", "fading", "inbound", "outbound"]
    rng = np.random.default_rng(1)
    for i, rho in enumerate(RHOS):
        for j, leg in enumerate(order):
            s = pr[(pr["rho_km"] == rho) & (pr["leg"] == leg)]
            if s.empty:
                continue
            xpos = i * 5 + j + rng.uniform(-.12, .12, len(s))
            ax[0].errorbar(xpos, s["x"], yerr=s["x_err_scaled"], fmt="o", ms=5, color=LEG_COLOUR[leg],
                           alpha=.85, lw=.8, capsize=0, label=f"{leg} (n={len(s)})" if i == 0 else None)
    ax[0].set(xticks=[1.5, 6.5], xticklabels=["10,000 km", "20,000 km"], ylabel="power-law index x",
              title="primary-leg indices, r-band, grades A/B")
    ax[0].axhline(0, color="0.6", lw=1); ax[0].legend(fontsize=9, frameon=False)
    b = breaks[(breaks["band"] == "r") & breaks["primary"] & breaks["preferred"]]
    for i, rho in enumerate(RHOS):
        s = b[b["rho_km"] == rho]
        if s.empty:
            continue
        ax[1].errorbar(s["r_break"], s["x_outer"] - s["x_inner"],
                       xerr=[np.clip(s["r_break"] - s["r_break_lo"], 0, None), np.clip(s["r_break_hi"] - s["r_break"], 0, None)],
                       fmt="o" if rho == 10000 else "s", ms=6, color="k" if rho == 10000 else "0.5",
                       lw=.8, capsize=0, label=f"{rho // 1000}k km (n={len(s)})")
        for _, r in s.iterrows():
            ax[1].annotate(r["target"], (r["r_break"], r["x_outer"] - r["x_inner"]), fontsize=7,
                           xytext=(3, 3), textcoords="offset points")
    ax[1].axvline(3.0, color="0.7", ls="-.", lw=1); ax[1].axhline(0, color="0.6", lw=1)
    ax[1].set(xscale="log", xlabel=r"$r_h$ of the break (au)", ylabel=r"$x_{\rm outer} - x_{\rm inner}$",
              title="broken laws preferred by ΔBIC > 6")
    _plain_log_axis(ax[1]); ax[1].legend(fontsize=9, frameon=False)
    c = colours[colours["leg"] == "all"] if len(colours) else colours
    for i, rho in enumerate(RHOS):
        s = c[c["rho_km"] == rho] if len(c) else c
        if len(s) == 0:
            continue
        ax[2].scatter(0.5 * (s["rh_min"] + s["rh_max"]), s["excess_median"], s=12 + s["n_pairs"] / 2,
                      color="tab:purple" if rho == 10000 else "0.5", alpha=.7,
                      label=f"{rho // 1000}k km, median excess per comet (n={len(s)})")
        ch = s[s["change_preferred"]]
        if len(ch):
            ax[2].errorbar(ch["r_change"], ch["colour_above"], xerr=[np.clip(ch["r_change"] - ch["r_change_lo"], 0, None),
                           np.clip(ch["r_change_hi"] - ch["r_change"], 0, None)], fmt="^", ms=7,
                           color="crimson" if rho == 10000 else "0.3", lw=.8, capsize=0,
                           label=f"{rho // 1000}k: preferred change (n={len(ch)})")
            for _, r in ch.iterrows():
                ax[2].annotate(f"{r['target']} {r['delta_colour']:+.2f}", (r["r_change"], r["colour_above"]),
                               fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[2].axhline(0, color="0.6", lw=1); ax[2].axvline(3.0, color="0.7", ls="-.", lw=1)
    ax[2].set(xscale="log", xlabel=r"$r_h$ (au)", ylabel="(g − r) excess over solar (mag)",
              title="dust colour and where it changes")
    _plain_log_axis(ax[2]); ax[2].legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--targets", nargs="+")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    cpath = OUT / "elements.csv"
    cache = (pd.read_csv(cpath, dtype={"target": str}).set_index("target") if cpath.exists()
             else pd.DataFrame(columns=["record", "q", "e", "Tp_jd", "epoch_jd"]).rename_axis("target"))
    files = sorted(glob.glob(str(zc.result_dir("photometry") / "*.csv")))
    if args.targets:
        want = {zc.target_slug(t) for t in args.targets}
        files = [f for f in files if Path(f).stem in want]

    names = designations()
    T, P, B, C = [], [], [], []
    for f in files:
        target = Path(f).stem
        phot = pd.read_csv(f, dtype={"target": str})
        if phot.empty:
            continue
        tp = perihelion_time(target, names.get(target, target), phot, cache)
        tr, pk, bk, cl = ac.analyse(phot, target, tp_jd=tp)
        T.append(tr); P.append(pk); B.append(bk); C.append(cl)
        r10 = tr[(tr["band"] == "r") & (tr["rho_km"] == 10000) & tr["primary"]] if len(tr) else tr
        line = "  ".join(f"{r['leg']}: x={r['x']:.2f}±{r['x_err_scaled']:.2f} n={int(r['n'])} [{r['grade']}]"
                         for _, r in r10.iterrows() if np.isfinite(r["x"]))
        if len(pk):
            p = pk[(pk["band"] == "r") & (pk["rho_km"] == 10000)]
            if len(p):
                p = p.iloc[0]; line += f"  peak {p['t_peak']:+.0f}d{'' if p['bracketed'] else '(edge)'}"
        b10 = bk[(bk["band"] == "r") & (bk["rho_km"] == 10000) & bk["primary"] & bk["preferred"]] if len(bk) else bk
        for _, b in b10.iterrows():
            line += f"  break@{b['r_break']:.2f}au({b['x_inner']:.1f}→{b['x_outer']:.1f})"
        c10 = cl[(cl["rho_km"] == 10000) & (cl["leg"] == "all")] if len(cl) else cl
        for _, c in c10.iterrows():
            line += f"  colour{'✓' if c['change_preferred'] else '–'}@{c['r_change']:.2f}au Δ{c['delta_colour']:+.2f}"
        print(f"{target:10s} Tp={'ok' if np.isfinite(tp) else '--'}  {line}", flush=True)
        if not args.no_figures and len(tr):
            fig_target(target, phot, tr, pk, bk, cl, tp, zc.fig_path("afrho", "trend.png", target))
    cache.to_csv(cpath)
    cat = lambda L: pd.concat([x for x in L if len(x)], ignore_index=True) if any(len(x) for x in L) else pd.DataFrame()
    trends, peaks, breaks, colours = cat(T), cat(P), cat(B), cat(C)
    trends.to_csv(OUT / "trends.csv", index=False); peaks.to_csv(OUT / "peaks.csv", index=False)
    breaks.to_csv(OUT / "breaks.csv", index=False); colours.to_csv(OUT / "colour.csv", index=False)
    if not args.no_figures and len(trends):
        fig_overview(trends, peaks, breaks, colours, zc.fig_path("afrho", "survey_overview.png"))
    print(f"\n{len(trends)} trend rows, {len(peaks)} peak rows, {len(breaks)} break rows, "
          f"{len(colours)} colour rows -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
