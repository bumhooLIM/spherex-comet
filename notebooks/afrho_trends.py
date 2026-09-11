#!/usr/bin/env python
"""Heliocentric dependence of Af-rho for every survey comet.

For each comet, band (r as the main series, g separately) and aperture
(10,000 and 20,000 km): fit ``Afrho = A r_h^-x`` on each phase of activity,
grade every slope, test whether one power law serves the whole r_h range or
a broken one is needed (free break, and a break imposed at 3 au), follow the
g - r colour of the dust with r_h, and find outbursts.

Phases are defined by the **peak of activity, not perihelion**: a comet
sampled on both sides of T_p gets a smoothed peak and splits into *rising*
and *fading* there when the maximum is interior to the coverage; a
one-sided comet is one leg named by its orbital direction.  Outburst windows
and small isolated groups of points at the ends of the r_h range are
removed from the primary fits and fitted separately; where a broken law is
preferred the two sides are fitted as segments -- the "divided" trend.
Science lives in ``ztfcomet.activity``; this script is selection,
bookkeeping and figures.  Horizons is queried once per comet for T_p and
cached; the SPHEREx observing windows come from ``spherex_windows.py``.

Outputs
  results/afrho/trends.csv, peaks.csv, breaks.csv, colour.csv, outbursts.csv
  results/afrho/elements.csv           cached q, e, Tp_jd
  fig/afrho/trend/<target>.png         per comet: fits, peak, break, outbursts, colour, SPHEREx windows
  fig/afrho/survey_overview.png        slope, break and colour distributions
"""

from __future__ import annotations

import argparse
import glob
import subprocess
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
DRAW_GRADES = {"A", "B", "C"}          # a grade-D line is arithmetic, not a trend


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


def spherex_windows():
    """Per (comet, phase) JD and r_h ranges of the SPHEREx exposures, built on demand."""
    path = zc.result_path("afrho", "spherex_windows.csv", create=False)
    if not path.exists():
        subprocess.run([sys.executable, str(Path(__file__).with_name("spherex_windows.py"))], check=False)
    return pd.read_csv(path, dtype={"target": str}) if path.exists() else pd.DataFrame()


# --------------------------------------------------------------- figures
def _plain_log_axis(ax):
    for f in (ax.xaxis.set_major_formatter, ax.xaxis.set_minor_formatter):
        f(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_minor_formatter().set_scientific(False)


def _label(f):
    return f"{f['leg']}: x = {f['x']:.2f} ± {f['x_err_scaled']:.2f}  n={int(f['n'])}  [{f['grade']}]"


def _shade_spherex(ax, win, key, lo_col, hi_col, ymax_frac=0.97):
    """Light bands where SPHEREx observed, labelled by phase number."""
    for _, w in win.iterrows():
        ax.axvspan(w[lo_col], w[hi_col], color="tab:cyan", alpha=0.18, lw=0, zorder=0)
        ax.text(0.5 * (w[lo_col] + w[hi_col]), ymax_frac, f"S{int(w['phase'])}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8, color="tab:cyan", fontweight="bold")


def fig_target(target, phot, trends, peaks, breaks, colours, outbursts, windows, tp_jd, outpath):
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))
    two_sided = bool(trends["two_sided"].any())
    win = windows[windows["target"] == target] if len(windows) else windows
    for col, rho in enumerate(RHOS):
        tr = trends[(trends["rho_km"] == rho) & (trends["band"] == "r")]
        prim = tr[tr["primary"] & tr["split"].isin(["peak", "perihelion", "all"])]
        segs = tr[tr["primary"] & (tr["split"] == "segment")]
        ob = outbursts[(outbursts["rho_km"] == rho) & (outbursts["band"] == "r")] if len(outbursts) else outbursts
        # ---- row 0: Afrho vs r_h
        ax = axes[0, col]
        pts, _ = ac.select_points(phot, "r", rho)
        gpts, _ = ac.select_points(phot, "g", rho)
        if len(win):
            _shade_spherex(ax, win, "phase", "rh_min", "rh_max")
        if not gpts.empty:
            ax.errorbar(gpts["r"], gpts["afrho0_cm"], yerr=gpts["afrho0_cm_err"], fmt="s", ms=3.5, mfc="none",
                        color="0.6", alpha=.5, lw=.6, capsize=0, zorder=1, label="g-band")
        leg_of = pd.Series("all", index=pts.index)
        in_ob = pd.Series(False, index=pts.index)
        if not pts.empty:
            for _, w in ob.iterrows():
                in_ob |= (pts["obsjd"] >= w["jd_start"]) & (pts["obsjd"] <= w["jd_end"])
            if tp_jd is not None and np.isfinite(tp_jd) and len(prim) and (prim["split"] == "peak").all():
                t = pts["obsjd"] - tp_jd
                tpk = float(prim["t_peak"].iloc[0])
                if bool(prim["peak_interior"].iloc[0]) and set(prim["leg"]) >= {"rising", "fading"}:
                    leg_of = pd.Series(np.where(t < tpk, "rising", "fading"), index=pts.index)
                else:
                    leg_of = pd.Series(prim["leg"].iloc[0], index=pts.index)
            elif len(prim):
                leg_of = pts["leg"]
            # tails: recompute per primary leg the same way analyse did
            in_tail = pd.Series(False, index=pts.index)
            for leg, s in pts[~in_ob].groupby(leg_of[~in_ob]):
                m = ac.sparse_tail(s["log_rh"].to_numpy())
                in_tail.loc[s.index[m]] = True
            for leg, s in pts.groupby(leg_of):
                c = LEG_COLOUR.get(leg, "0.4")
                core = s[~in_ob[s.index] & ~in_tail[s.index]]
                ax.errorbar(core["r"], core["afrho0_cm"], yerr=core["afrho0_cm_err"], fmt="o", ms=5.5, color=c,
                            lw=.8, capsize=0, zorder=2, label=f"r  {leg}")
                hollow = s[in_ob[s.index] | in_tail[s.index]]
                if len(hollow):
                    ax.errorbar(hollow["r"], hollow["afrho0_cm"], yerr=hollow["afrho0_cm_err"], fmt="o", ms=5.5,
                                mfc="none", mec=c, ecolor=c, lw=.8, capsize=0, zorder=2)
            quiet = pts[~in_ob]
            g, sm = ac.smooth_trend(quiet["log_rh"], quiet["log_afrho"], quiet["log_afrho_err"])
            if len(g):
                ax.plot(10 ** g, 10 ** sm, "-", color="k", lw=2.2, alpha=.85, zorder=4, label="smoothed, all phases")
        txt = []
        for _, f in prim.iterrows():
            if not np.isfinite(f["x"]):
                continue
            if f["grade"] in DRAW_GRADES:
                rr = np.geomspace(f["rh_min"], f["rh_max"], 50)
                ax.plot(rr, 10 ** (f["a"] - f["x"] * np.log10(rr)), ":" if f["has_segments"] else "-",
                        color=LEG_COLOUR[f["leg"]], lw=1.6 if f["has_segments"] else 2.2, zorder=3)
                txt.append(_label(f))
            else:
                txt.append(f"{f['leg']}: unconstrained  n={int(f['n'])}  [{f['grade']}]")
            if f["n_outburst"]:
                txt[-1] += f"   ({int(f['n_outburst'])} outburst pts excluded)"
            if f["n_tail"]:
                txt[-1] += f"   ({int(f['n_tail'])} isolated pts excluded)"
        for _, f in segs.iterrows():
            if np.isfinite(f["x"]) and f["grade"] in DRAW_GRADES:
                rr = np.geomspace(f["rh_min"], f["rh_max"], 50)
                ax.plot(rr, 10 ** (f["a"] - f["x"] * np.log10(rr)), "-", color=LEG_COLOUR[f["segment_of"]], lw=2.6, zorder=5)
            txt.append(f"  segment {f['leg']}: x = {f['x']:.2f} ± {f['x_err_scaled']:.2f}  n={int(f['n'])}  [{f['grade']}]")
        bk = (breaks[(breaks["rho_km"] == rho) & (breaks["band"] == "r") & breaks["primary"]]
              if len(breaks) else breaks)
        for _, b in bk.iterrows():
            if b["preferred"]:
                ax.axvline(b["r_break"], color="k", ls=":", lw=1)
                txt.append(f"{b['leg']}: break at {b['r_break']:.2f} au [{b['r_break_lo']:.2f}, {b['r_break_hi']:.2f}]  (ΔBIC {b['dbic']:.0f})")
            if b["at3_testable"]:
                txt.append(f"{b['leg']} at 3 au: x(<3) = {b['x_lt3']:.2f} ± {b['x_lt3_err']:.2f},  x(>3) = {b['x_gt3']:.2f} ± {b['x_gt3_err']:.2f}")
        for _, w in ob.iterrows():
            txt.append(f"outburst at {w['rh_start']:.2f} au, +{w['rise_dex']:.2f} dex, {int(w['n_points'])} pts"
                       + ("" if w["ended"] else " (not ended in coverage)"))
        if not pts.empty and pts["r"].min() < 3.0 < pts["r"].max():
            ax.axvline(3.0, color="0.7", ls="-.", lw=1)
        ax.set(xscale="log", yscale="log", xlabel=r"$r_h$ (au)", ylabel=r"$A(0°)f\rho$ (cm)",
               title=rf"{target}   $\rho$ = {rho:,} km")
        _plain_log_axis(ax)
        ax.text(0.02, 0.02, "\n".join(txt) or "no fit", transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(fc="white", ec="0.7", alpha=0.85))
        ax.legend(fontsize=7.5, frameon=False, loc="upper right")
        # ---- row 1: time from perihelion
        ax2 = axes[1, col]
        if tp_jd is not None and np.isfinite(tp_jd) and not pts.empty:
            tt = pts["obsjd"] - tp_jd
            if len(win):
                _shade_spherex(ax2, win.assign(t_lo=win["jd_min"] - tp_jd, t_hi=win["jd_max"] - tp_jd), "phase", "t_lo", "t_hi")
            for leg, s in pts.groupby(leg_of):
                ax2.errorbar(tt[s.index], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt="o", ms=5,
                             color=LEG_COLOUR.get(leg, "0.4"), lw=.8, capsize=0)
            for _, w in ob.iterrows():
                ax2.axvspan(w["t_start"], w["t_end"], color="tab:orange", alpha=0.18, lw=0)
            ax2.axvline(0, color="crimson", ls="--", lw=1, label="perihelion")
            pk = peaks[(peaks["rho_km"] == rho) & (peaks["band"] == "r")] if len(peaks) else peaks
            if len(pk):
                p = pk.iloc[0]
                if p["interior"]:
                    ax2.axvspan(p["t_lo"], p["t_hi"], color="tab:green", alpha=0.15)
                    ax2.axvline(p["t_peak"], color="tab:green", lw=2,
                                label=f"peak {p['t_peak']:+.0f} d  [{p['t_lo']:+.0f}, {p['t_hi']:+.0f}]"
                                      + ("" if p["bracketed"] else "  (plateau before)" if p["drop_lo"] < p["drop_hi"] else "  (plateau after)"))
                else:
                    ax2.text(0.02, 0.95, f"maximum at the edge of coverage ({p['t_peak']:+.0f} d): "
                             f"{'still rising' if p['t_peak'] > tt.median() else 'already fading'}",
                             transform=ax2.transAxes, fontsize=9, va="top")
            elif not two_sided:
                ax2.text(0.02, 0.95, "one-sided coverage: peak not measurable", transform=ax2.transAxes, fontsize=9, va="top", color="0.4")
            ax2.set(yscale="log", xlabel=r"$T - T_p$ (days)", ylabel=r"$A(0°)f\rho$ (cm)")
            ax2.legend(fontsize=8, frameon=False, loc="lower right")
        else:
            ax2.text(0.5, 0.5, "no perihelion time", ha="center", va="center", transform=ax2.transAxes, color="0.4")
            ax2.set_axis_off()
        # ---- row 2: colour
        ax3 = axes[2, col]
        pairs = ac.colour_pairs(phot, rho)
        cl = colours[(colours["rho_km"] == rho) & (colours["leg"] == "all")] if len(colours) else colours
        if len(pairs) >= 3:
            if len(win):
                _shade_spherex(ax3, win, "phase", "rh_min", "rh_max")
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
            ax3.set(xscale="log", xlabel=r"$r_h$ (au)", ylabel="colour excess over solar (mag)")
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
    pr = trends[trends["primary"] & trends["split"].isin(["peak", "perihelion"]) & (trends["band"] == "r")
                & trends["grade"].isin(["A", "B"])]
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
              title="primary-phase indices, r-band, grades A/B")
    ax[0].axhline(0, color="0.6", lw=1); ax[0].legend(fontsize=9, frameon=False)
    b = breaks[(breaks["band"] == "r") & breaks["primary"] & breaks["preferred"]] if len(breaks) else breaks
    for i, rho in enumerate(RHOS):
        s = b[b["rho_km"] == rho] if len(b) else b
        if len(s) == 0:
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
              title="broken laws preferred by ΔBIC > 6 (scaled errors)")
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
    windows = spherex_windows()
    names = designations()
    T, P, B, C, O = [], [], [], [], []
    for f in files:
        target = Path(f).stem
        phot = pd.read_csv(f, dtype={"target": str})
        if phot.empty:
            continue
        tp = perihelion_time(target, names.get(target, target), phot, cache)
        tr, pk, bk, cl, ob = ac.analyse(phot, target, tp_jd=tp)
        T.append(tr); P.append(pk); B.append(bk); C.append(cl); O.append(ob)
        r10 = tr[(tr["band"] == "r") & (tr["rho_km"] == 10000) & tr["primary"]] if len(tr) else tr
        line = "  ".join(f"{r['leg']}: x={r['x']:.2f}±{r['x_err_scaled']:.2f} n={int(r['n'])} [{r['grade']}]"
                         for _, r in r10.iterrows() if np.isfinite(r["x"]))
        if len(pk):
            p = pk[(pk["band"] == "r") & (pk["rho_km"] == 10000)]
            if len(p):
                p = p.iloc[0]; line += f"  peak {p['t_peak']:+.0f}d{'' if p['interior'] else '(edge)'}"
        if len(ob):
            o = ob[(ob["band"] == "r") & (ob["rho_km"] == 10000)]
            line += "".join(f"  OUTBURST@{w['rh_start']:.2f}au(+{w['rise_dex']:.2f}dex,{int(w['n_points'])}pts)" for _, w in o.iterrows())
        print(f"{target:10s} Tp={'ok' if np.isfinite(tp) else '--'}  {line}", flush=True)
        if not args.no_figures and len(tr):
            fig_target(target, phot, tr, pk, bk, cl, ob, windows, tp, zc.fig_kind_path("afrho", "trend", target))
    cache.to_csv(cpath)
    cat = lambda L: pd.concat([x for x in L if len(x)], ignore_index=True) if any(len(x) for x in L) else pd.DataFrame()
    trends, peaks, breaks, colours, outbursts = cat(T), cat(P), cat(B), cat(C), cat(O)
    for name, df in (("trends", trends), ("peaks", peaks), ("breaks", breaks), ("colour", colours), ("outbursts", outbursts)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    if not args.no_figures and len(trends):
        fig_overview(trends, peaks, breaks, colours, zc.fig_path("afrho", "survey_overview.png"))
    print(f"\n{len(trends)} trend rows, {len(peaks)} peak rows, {len(breaks)} break rows, "
          f"{len(colours)} colour rows, {len(outbursts)} outbursts -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
