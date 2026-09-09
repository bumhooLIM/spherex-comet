#!/usr/bin/env python
"""Where does the 1/rho coma law stop holding?  A survey-wide test.

The 24P study certified 1/rho only to ~12,000 km because that is as far as
10 px reaches at Delta ~ 0.8 au -- a floor set by coverage, not a detected
turnover.  The survey's distant targets push 10 px to 30-40,000 km, and about
1,150 clean frames do so with outer-annulus S/N > 5, so the question can be
asked directly.

Two biases are known to steepen a naive slope: the PSF core, which at
Delta = 4.5 au spans ~10,000 km, and a sky error on faint frames.  Both are
modelled per frame with the PSF-convolved coma model and a free sky offset,
fitted twice -- slope free, and slope fixed at 1 -- so a steepening that
survives is one neither bias can produce.

The diagnostic is the ratio of the observed profile to the fitted m = 1 model
(nucleus + PSF-convolved 1/rho coma + sky).  It is 1 wherever 1/rho holds, at
every radius including inside the PSF core, and it is binned in km so that a
departure tied to physical distance separates from one tied to pixels or S/N.

Frames whose comet profile does not peak in the innermost annulus are
excluded and counted: a windowed centroid on an asymmetric coma lands on the
light centre rather than the nucleus, and the model's nucleus term then has
nothing to fit.  Uses the profile tables from ``ztfcomet.profile`` -- no FITS,
no network.

Outputs
  results/profile_survey_fits.csv       per frame: window, S/N, both fits
  results/profile_survey_ratio.csv      per (frame, annulus): obs / m=1 model
  results/profile_survey_targets.csv    per target: rho_max and median slopes
  fig/survey/profile_survey_slope.png   corrected slope vs physical window
  fig/survey/profile_survey_ratio.png   obs / (1/rho model) stacked in km
"""

from __future__ import annotations

import argparse
import glob
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ztfcomet as zc
from ztfcomet import profile as pf
from ztfcomet import rcparams  # noqa: F401
from profile_resolution import load

SNR_BANDS = [-np.inf, 5, 20, 60, np.inf]
SNR_LABELS = ["<5", "5-20", "20-60", ">60"]
KM_EDGES = np.geomspace(500, 50000, 25)


def survey_targets():
    return sorted(Path(p).parent.name
                  for p in glob.glob(str(zc.RESULT_ROOT / "*" / "profile_summary_*.csv")))


# ------------------------------------------------------------------ fits
def fit_target(target, snr_min=5.0, oversample=4):
    """Both model fits for every usable clean frame of one target."""
    prof, summ = load(target)
    peak = prof.pivot(index="file", columns="r_pix", values="mean_clip").fillna(-np.inf).idxmax(axis=1)
    summ = summ.set_index("file")
    summ["peak_r"] = peak
    summ = summ.reset_index()
    use = summ["clean"] & (summ["sb10_over_skyerr"] > snr_min)
    rows, ratio = [], []
    n_offpeak = int((use & (summ["peak_r"] > 0.5)).sum())
    for _, s in summ[use & (summ["peak_r"] <= 0.5)].iterrows():
        q = prof[prof["file"] == s["file"]].sort_values("r_pix")
        obs = q["mean_clip"].to_numpy(float)
        if np.isfinite(obs).sum() < 6:
            continue
        stack = pd.DataFrame({"r_pix": q["r_pix"], "median": q["star_median_norm"],
                              "std": q["star_std_norm"], "n_stars": q["n_stars"]})
        comet = q[["r_pix", "mean_clip", "std_clip", "n_clip"]].reset_index(drop=True)
        rec = dict(target=target, file=s["file"], leg=s["leg"], r=s["r"], delta=s["delta"],
                   km_per_pix=s["km_per_pix"], fwhm_pix=s["fwhm_pix"],
                   rho_fit_rmin_km=s["rho_fit_rmin_km"], rho_rmax_km=s["rho_rmax_km"],
                   snr10=s["sb10_over_skyerr"], slope_naive=s["slope_comet"])
        # Same bound the model uses; a fit that sits on it has not converged
        # to a sky level, it has run out of room.
        s_bound = 3.0 * float(np.nanmax(np.abs(obs[np.isfinite(obs)][-3:]))) + 1e-3
        try:
            free = pf.fit_coma_model(comet, stack, oversample=oversample, fit_sky=True)
            one = pf.fit_coma_model(comet, stack, oversample=oversample, fit_sky=True, m_fixed=1.0)
        except Exception as exc:                                # noqa: BLE001
            rec["error"] = str(exc)[:60]
            rows.append(rec)
            continue
        n = int(one["n_points"])
        rec.update(m_free=free["m"], sky_free=free["sky_offset"], chi2_free=free["chi2_red"],
                   nuc_free=free["nucleus_fraction"], sky_one=one["sky_offset"],
                   chi2_one=one["chi2_red"], nuc_one=one["nucleus_fraction"], n_points=n,
                   dchi2=n * (one["chi2_red"] - free["chi2_red"]),
                   # n_clip counts oversampled sub-pixels, which are bilinear
                   # interpolates of the same data: the annulus error is ~1/oversample
                   # too small and any chi-square ~oversample**2 too large.  The
                   # effective value is the one to threshold.
                   dchi2_eff=n * (one["chi2_red"] - free["chi2_red"]) / oversample ** 2,
                   at_bound=(abs(free["sky_offset"]) > 0.95 * s_bound
                             or abs(one["sky_offset"]) > 0.95 * s_bound),
                   ok=bool(free["success"] and one["success"]))
        rows.append(rec)
        # Ratio to the 1/rho model, sky removed from both sides.
        model = np.asarray(one["model"], float)
        radii = q["r_pix"].to_numpy(float)
        keep = np.isfinite(obs) & (radii <= 10.0)
        if model.shape[0] != keep.sum():
            model = model[:keep.sum()] if model.shape[0] > keep.sum() else model
        coma = model - one["sky_offset"]
        good = coma > 0
        ratio.append(pd.DataFrame(dict(
            target=target, file=s["file"], r_pix=radii[keep][good],
            rho_km=q["rho_km"].to_numpy(float)[keep][good],
            y=(obs[keep][good] - one["sky_offset"]) / coma[good])))
    return pd.DataFrame(rows), (pd.concat(ratio, ignore_index=True) if ratio else pd.DataFrame()), n_offpeak


# ------------------------------------------------------------- analysis
def bin_ratio(ratio, fits, mask, min_frames=8):
    """Median obs/(1/rho model) in km bins, for frames satisfying *mask*."""
    files = set(fits.loc[mask, "file"])
    r = ratio[ratio["file"].isin(files)].copy()
    r["bin"] = pd.cut(r["rho_km"], KM_EDGES)
    g = r.groupby("bin", observed=True)["y"].agg(["median", "count",
                                                  lambda v: v.quantile(.16), lambda v: v.quantile(.84)])
    g.columns = ["median", "count", "q16", "q84"]
    g = g[g["count"] >= min_frames].reset_index()
    g["rho_km"] = g["bin"].apply(lambda b: np.sqrt(b.left * b.right))
    return g


def rho_max_from(binned, floor=0.8):
    if binned.empty:
        return np.nan
    b = binned.sort_values("rho_km")
    ok = b[b["median"] >= floor]
    return float(ok["rho_km"].max()) if len(ok) else np.nan


def per_target(fits, ratio):
    rows = []
    for t, f in fits.groupby("target"):
        good = f["ok"].fillna(False) & ~f["at_bound"].fillna(True)
        hi = good & (f["snr10"] > 20)
        b = bin_ratio(ratio, f, hi, min_frames=5)
        rows.append(dict(target=t, n_used=len(f), n_good=int(good.sum()), n_snr20=int(hi.sum()),
                         delta_med=f["delta"].median(), r_med=f["r"].median(),
                         rho10_med=f["rho_rmax_km"].median(),
                         m_naive=-f.loc[good, "slope_naive"].median(),
                         m_free=f.loc[good, "m_free"].median(),
                         m_free_snr20=f.loc[hi, "m_free"].median() if hi.any() else np.nan,
                         frac_m1_ok=(f.loc[good, "dchi2_eff"] < 9).mean() if good.any() else np.nan,
                         rho_max_km=rho_max_from(b)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------- figures
def fig_slope(fits, outpath):
    g = fits[fits["ok"].fillna(False) & ~fits["at_bound"].fillna(True)].copy()
    g["band"] = pd.cut(g["snr10"], SNR_BANDS, labels=SNR_LABELS)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    cols = dict(zip(SNR_LABELS, ["0.7", "tab:blue", "tab:orange", "crimson"]))
    xb = np.geomspace(3000, 45000, 10)
    for band, s in g.groupby("band", observed=True):
        ax[0].scatter(s["rho_rmax_km"], s["m_free"], s=8, c=cols[band], alpha=.4)
        med = s.groupby(pd.cut(s["rho_rmax_km"], xb), observed=True)["m_free"].median()
        xc = [np.sqrt(b.left * b.right) for b in med.index]
        ax[0].plot(xc, med.values, "o-", c=cols[band], ms=5, label=f"S/N(10 px) {band}  n={len(s)}")
    ax[0].axhline(1, color="tab:green", ls="-.")
    ax[0].set(xscale="log", xlabel=r"$\rho$ at 10 px (km)", ylabel="coma slope m  (PSF model + free sky)",
              ylim=(-0.5, 3.5), title="does the corrected slope steepen with the physical window?")
    ax[0].legend(fontsize=9, frameon=False)
    acc = g.assign(ok1=g["dchi2_eff"] < 9)
    for band, s in acc.groupby("band", observed=True):
        m = s.groupby(pd.cut(s["rho_rmax_km"], xb), observed=True)["ok1"].agg(["mean", "size"])
        m = m[m["size"] >= 8]
        xc = [np.sqrt(b.left * b.right) for b in m.index]
        ax[1].plot(xc, m["mean"], "o-", c=cols[band], ms=5, label=f"S/N {band}")
    ax[1].set(xscale="log", xlabel=r"$\rho$ at 10 px (km)", ylabel=r"fraction with effective $\Delta\chi^2(m=1) < 9$",
              ylim=(0, 1.05), title="is m = 1 statistically acceptable?")
    ax[1].legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)


def fig_ratio(fits, ratio, outpath):
    good = fits["ok"].fillna(False) & ~fits["at_bound"].fillna(True)
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    cols = dict(zip(SNR_LABELS, ["0.7", "tab:blue", "tab:orange", "crimson"]))
    band = pd.cut(fits["snr10"], SNR_BANDS, labels=SNR_LABELS)
    for lab in SNR_LABELS:
        b = bin_ratio(ratio, fits, good & (band == lab))
        if b.empty:
            continue
        ax[0].plot(b["rho_km"], b["median"], "o-", c=cols[lab], ms=4, label=f"S/N(10 px) {lab}")
    ax[0].axhline(1, color="tab:green", ls="-.", label=r"$1/\rho$")
    ax[0].axhline(0.8, color="0.5", ls=":", lw=1)
    ax[0].set(xscale="log", xlabel=r"$\rho$ (km)", ylabel=r"observed / ($1/\rho$ model)", ylim=(0.3, 1.6),
              title="departure tied to km, or to S/N?")
    ax[0].legend(fontsize=9, frameon=False)
    hi = good & (fits["snr10"] > 20)
    rt = pd.qcut(fits["r"], 3, labels=["inner", "mid", "outer"])
    for lab, c in zip(["inner", "mid", "outer"], ["tab:blue", "tab:orange", "crimson"]):
        m = hi & (rt == lab)
        if m.sum() < 8:
            continue
        b = bin_ratio(ratio, fits, m)
        rr = fits.loc[m, "r"]
        ax[1].plot(b["rho_km"], b["median"], "o-", c=c, ms=4,
                   label=f"$r_h$ {rr.min():.1f}-{rr.max():.1f} au  n={int(m.sum())}")
        ax[1].fill_between(b["rho_km"], b["q16"], b["q84"], color=c, alpha=.12)
    ax[1].axhline(1, color="tab:green", ls="-.")
    ax[1].axhline(0.8, color="0.5", ls=":", lw=1)
    ax[1].set(xscale="log", xlabel=r"$\rho$ (km)", ylabel=r"observed / ($1/\rho$ model)", ylim=(0.3, 1.6),
              title=r"S/N > 20 only: does the break move with $r_h$?")
    ax[1].legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(outpath); plt.close(fig)


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--targets", nargs="+", help="subset (default: every target with profiles)")
    ap.add_argument("--snr-min", type=float, default=5.0)
    args = ap.parse_args(argv)
    targets = [zc.target_slug(t) for t in args.targets] if args.targets else survey_targets()

    fits, ratios, offpeak = [], [], {}
    t0 = time.time()
    for i, t in enumerate(targets, 1):
        f, r, n_off = fit_target(t, snr_min=args.snr_min)
        fits.append(f); ratios.append(r); offpeak[t] = n_off
        print(f"[{i}/{len(targets)}] {t:10s} fitted {len(f):4d}  off-peak excluded {n_off:3d}  "
              f"({time.time() - t0:.0f} s)", flush=True)
    fits = pd.concat([f for f in fits if len(f)], ignore_index=True)
    ratio = pd.concat([r for r in ratios if len(r)], ignore_index=True)
    fits.to_csv(zc.RESULT_ROOT / "profile_survey_fits.csv", index=False)
    ratio.to_csv(zc.RESULT_ROOT / "profile_survey_ratio.csv", index=False)

    good = fits["ok"].fillna(False) & ~fits["at_bound"].fillna(True)
    print(f"\nframes fitted {len(fits)}, converged off the sky bound {int(good.sum())}, "
          f"off-peak excluded {sum(offpeak.values())}")
    band = pd.cut(fits["snr10"], SNR_BANDS, labels=SNR_LABELS)
    print("\n=== corrected slope by S/N band (median [16-84%]) and m=1 acceptance ===")
    for lab in SNR_LABELS:
        s = fits[good & (band == lab)]
        if len(s) < 5:
            continue
        m = s["m_free"]
        print(f"  S/N {lab:6s} n={len(s):4d}  m={m.median():.2f} [{m.quantile(.16):.2f},{m.quantile(.84):.2f}]"
              f"  naive {-s['slope_naive'].median():.2f}  m=1 accepted {100 * (s['dchi2_eff'] < 9).mean():.0f}%")
    print("\n=== obs / (1/rho model) in km, by S/N band ===")
    for lab in SNR_LABELS:
        b = bin_ratio(ratio, fits, good & (band == lab))
        if b.empty:
            continue
        print(f"  {lab:6s} " + "  ".join(f"{r / 1000:.0f}k:{m:.2f}" for r, m in zip(b["rho_km"], b["median"])))
        print(f"         rho_max (median >= 0.8): {rho_max_from(b):.0f} km")
    hi = good & (fits["snr10"] > 20)
    rt = pd.qcut(fits["r"], 3, labels=["inner", "mid", "outer"])
    print("\n=== S/N > 20, by r_h tercile ===")
    for lab in ["inner", "mid", "outer"]:
        m = hi & (rt == lab)
        b = bin_ratio(ratio, fits, m)
        if b.empty:
            continue
        rr = fits.loc[m, "r"]
        print(f"  r_h {rr.min():.1f}-{rr.max():.1f} au n={int(m.sum()):4d}: rho_max {rho_max_from(b):.0f} km   "
              + "  ".join(f"{r / 1000:.0f}k:{v:.2f}" for r, v in zip(b["rho_km"], b["median"])))

    tab = per_target(fits, ratio)
    tab["n_offpeak"] = tab["target"].map(offpeak)
    tab.to_csv(zc.RESULT_ROOT / "profile_survey_targets.csv", index=False)
    print("\n=== per target (sorted by rho at 10 px) ===")
    print(tab.sort_values("rho10_med", ascending=False)
          .to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    figdir = zc.fig_dir("survey")
    fig_slope(fits, figdir / "profile_survey_slope.png")
    fig_ratio(fits, ratio, figdir / "profile_survey_ratio.png")
    print("\nfigures written to", figdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
