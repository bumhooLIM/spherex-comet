#!/usr/bin/env python
"""Is the steep inbound coma profile physical, or a pixel-scale effect?

Tests the hypothesis that a comet's radial profile looks steeper at large
observer distance because a fixed 10 px window then spans a much larger
cometocentric distance, reaching beyond where the 1/rho law holds.  Uses the
profile tables written by ``ztfcomet.profile`` -- no FITS access, no network.

Outputs (for TARGET, default 24P)
  results/<T>/profile_resolution_<T>.csv     per-frame: scales, slopes by method
  results/profile_feasibility_survey.csv     per survey target: usable window
  fig/<T>/profile_resolution_<T>.png         hypothesis tests
  fig/<T>/profile_validity_<T>.png           where 1/rho holds, km vs pixels
  fig/<T>/profile_correction_<T>.png         naive vs fixed-km vs PSF-model slopes
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import ztfcomet as zc
from ztfcomet import profile as pf
from ztfcomet import rcparams  # noqa: F401


def load(target):
    paths = zc.profile_paths(target, create=False)
    S = {"target": str}
    prof = pd.read_csv(paths["profile"], dtype=S)
    summ = pd.read_csv(paths["summary"], dtype=S)
    phot = pd.read_csv(zc.photometry_path(target, create=False), dtype=S)
    phot = phot[np.isclose(phot["rho_km"], phot["rho_km"].min())].drop_duplicates("file")
    summ = summ.merge(phot[["file", "ssky", "nsky", "snr", "msky"]], on="file", how="left")
    summ["leg"] = np.where(summ["r_rate"] < 0, "pre", "post")
    summ["clean"] = summ["quality_ok"].astype(bool)
    # Outer-annulus S/N: SB at 10 px against the sky-level uncertainty.
    sb10 = prof[prof["r_pix"] == 10.0].set_index("file")["mean_clip"]
    summ["sb10"] = summ["file"].map(sb10)
    summ["sky_err"] = summ["ssky"] / np.sqrt(summ["nsky"].clip(lower=1))
    summ["sb10_over_skyerr"] = summ["sb10"] / summ["sky_err"]
    return prof, summ


# ---------------------------------------------------------------- A: scales
def hypothesis_tests(summ):
    out = {}
    g = summ[summ["clean"] & np.isfinite(summ["slope_comet"])]
    for leg in ("pre", "post", "all"):
        s = g if leg == "all" else g[g["leg"] == leg]
        if len(s) < 8:
            continue
        for col, lab in (("rho_rmax_km", "rho10km"), ("delta", "delta"),
                         ("sb10_over_skyerr", "snr10"), ("fwhm_km", "fwhm_km")):
            rho, p = spearmanr(s[col], s["slope_comet"], nan_policy="omit")
            out[f"{leg}:{lab}"] = (float(rho), float(p), len(s))
    return out


# --------------------------------------------------- B: where 1/rho holds
def rho_sb_diagnostic(prof, summ, min_frames=8):
    """rho*SB normalised per frame; flat == 1/rho.  Binned in km and in px."""
    # prof already carries delta (and r, r_rate) per frame; merge only the
    # per-frame fields that live in the summary.
    p = prof.merge(summ[["file", "leg", "clean", "fit_rmin_pix"]], on="file")
    p = p[p["clean"] & np.isfinite(p["mean_clip_norm"]) & (p["mean_clip_norm"] > 0)]
    p["y"] = p["rho_km"] * p["mean_clip_norm"]
    infit = (p["r_pix"] >= p["fit_rmin_pix"]) & (p["r_pix"] <= 10.0)
    norm = p[infit].groupby("file")["y"].median()
    p["y_norm"] = p["y"] / p["file"].map(norm)
    p = p[infit & np.isfinite(p["y_norm"])]

    km_edges = np.geomspace(800, 20000, 17)
    p["km_bin"] = pd.cut(p["rho_km"], km_edges)
    by_km = (p.groupby(["leg", "km_bin"], observed=True)["y_norm"]
             .agg(["median", "count"]).reset_index())
    by_km["rho_km"] = by_km["km_bin"].apply(lambda b: np.sqrt(b.left * b.right))
    by_km = by_km[by_km["count"] >= min_frames]

    # Same thing against pixel radius, split by observer-distance tercile.
    p["delta_bin"] = pd.qcut(p["delta"], 3, labels=["near", "mid", "far"])
    by_px = (p.groupby(["delta_bin", "r_pix"], observed=True)["y_norm"]
             .agg(["median", "count"]).reset_index())
    by_px = by_px[by_px["count"] >= min_frames]

    # rho_max: largest km bin (post-perihelion, best-resolved) where the
    # median rho*SB is still within 20% of the plateau (slope within ~-1.3).
    post = by_km[by_km["leg"] == "post"].sort_values("rho_km")
    plateau = float(post["median"].max()) if len(post) else np.nan
    ok = post[post["median"] >= 0.8 * plateau]
    rho_max = float(ok["rho_km"].max()) if len(ok) else np.nan
    return p, by_km, by_px, rho_max


# ------------------------------------------------------------ C/D/E: refits
def refit_methods(prof, summ, rho_max_km, oversample=4):
    """Per frame: naive px-window slope, fixed-km-window slope, native-pixel
    slope, and the PSF-convolved model slope."""
    rows = []
    for _, s in summ.iterrows():
        q = prof[prof["file"] == s["file"]].sort_values("r_pix")
        if q.empty:
            continue
        r, rk = q["r_pix"].to_numpy(), q["rho_km"].to_numpy()
        rmin = float(s["fit_rmin_pix"])
        rec = dict(file=s["file"], leg=s["leg"], clean=s["clean"], delta=s["delta"],
                   r=s["r"], km_per_pix=s["km_per_pix"], fwhm_pix=s["fwhm_pix"],
                   rho_rmax_km=s["rho_rmax_km"], sb10_over_skyerr=s["sb10_over_skyerr"],
                   slope_naive=s["slope_comet"], slope_star=s["slope_star"])
        rec["slope_native"], _, _ = pf.fit_powerlaw(r, q["mean_clip_native_norm"], rmin, 10.0)

        # Fixed physical window: from the PSF core out to rho_max, if that
        # leaves at least four annuli.
        lo_km = max(rmin * s["km_per_pix"], 1500.0)
        hi_km = min(rho_max_km, s["rho_rmax_km"])
        n_ann = int(np.sum((rk >= lo_km) & (rk <= hi_km)))
        rec["km_window_lo"], rec["km_window_hi"], rec["km_window_n"] = lo_km, hi_km, n_ann
        rec["slope_kmwin"] = (pf.fit_powerlaw(rk, q["mean_clip_norm"], lo_km, hi_km)[0]
                              if n_ann >= 4 else np.nan)

        # PSF-convolved forward model, from the stored comet + star profiles.
        stack = pd.DataFrame({"r_pix": r, "median": q["star_median_norm"],
                              "std": q["star_std_norm"], "n_stars": q["n_stars"]})
        comet = q[["r_pix", "mean_clip", "std_clip", "n_clip"]].reset_index(drop=True)
        try:
            fit = pf.fit_coma_model(comet, stack, oversample=oversample)
            rec.update(slope_model=-fit["m"], nucleus_fraction=fit["nucleus_fraction"],
                       model_chi2=fit["chi2_red"])
            fixed = pf.fit_coma_model(comet, stack, oversample=oversample, m_fixed=1.0)
            rec["chi2_fixed_m1"] = fixed["chi2_red"]
            # Free sky offset: the test of the sky-subtraction explanation.
            sky = pf.fit_coma_model(comet, stack, oversample=oversample, fit_sky=True)
            rec.update(slope_model_sky=-sky["m"], sky_offset=sky["sky_offset"],
                       sky_offset_over_err=sky["sky_offset"] / s["sky_err"] if s["sky_err"] > 0 else np.nan,
                       model_sky_chi2=sky["chi2_red"])
        except Exception as exc:                                # noqa: BLE001
            rec.update(slope_model=np.nan, nucleus_fraction=np.nan, model_chi2=np.nan,
                       chi2_fixed_m1=np.nan, model_error=str(exc)[:60])
        rows.append(rec)
    return pd.DataFrame(rows)


# --------------------------------------------------- F: survey feasibility
def survey_feasibility(rho_max_km, fit_rmin_fwhm=1.5, min_annuli=4, half_width=0.25):
    rows = []
    for path in sorted(glob.glob(str(zc.result_dir("photometry") / "*.csv"))):
        t = pd.read_csv(path, dtype={"target": str})
        if t.empty or "delta" not in t:
            continue
        t = t[np.isclose(t["rho_km"], t["rho_km"].min())].drop_duplicates("file")
        kpp = t["pixscale"] * t["delta"] * pf.KM_PER_ARCSEC_AU
        core_px = fit_rmin_fwhm * t["fwhm_pix"]
        rho_max_px = rho_max_km / kpp
        n_valid_px = (np.minimum(rho_max_px, 10.0) - core_px)
        feasible = n_valid_px >= min_annuli * 2 * half_width
        rows.append(dict(
            target=t["target"].iloc[0] if "target" in t else Path(path).stem,
            frames=len(t), delta_median=float(t["delta"].median()),
            km_per_pix_median=float(kpp.median()), fwhm_px_median=float(t["fwhm_pix"].median()),
            rho_max_px_median=float(rho_max_px.median()),
            valid_px_median=float(n_valid_px.median()),
            frac_feasible=float(feasible.mean()),
            frac_rho10_beyond_rhomax=float((10.0 * kpp > rho_max_km).mean()),
        ))
    return pd.DataFrame(rows)


# --------------------------------------------------------------- figures
def fig_hypothesis(summ, corr, outpath):
    g = summ[summ["clean"] & np.isfinite(summ["slope_comet"])]
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    col = {"pre": "tab:blue", "post": "tab:orange"}
    for leg, s in g.groupby("leg"):
        ax[0, 0].scatter(s["rho_rmax_km"], s["slope_comet"], s=18, c=col[leg], label=f"{leg}-perihelion")
        ax[0, 1].scatter(s["delta"], s["slope_comet"], s=18, c=col[leg])
        ax[1, 0].scatter(s["sb10_over_skyerr"], s["slope_comet"], s=18, c=col[leg])
    for a in ax.flat[:3]:
        a.axhline(-1, color="tab:green", ls="-.")
        a.set_ylim(-3.5, 0.5)
    ax[0, 0].set(xlabel=r"$\rho$ at 10 px (km)", ylabel="naive slope")
    r, p, n = corr.get("all:rho10km", (np.nan, np.nan, 0))
    ax[0, 0].set_title(f"slope vs outer radius in km   (Spearman {r:+.2f}, p={p:.1e}, n={n})", fontsize=12)
    ax[0, 0].legend(fontsize=10, frameon=False)
    ax[0, 1].set(xlabel=r"$\Delta$ (au)", ylabel="naive slope")
    r, p, n = corr.get("all:delta", (np.nan, np.nan, 0))
    ax[0, 1].set_title(f"slope vs observer distance   (Spearman {r:+.2f})", fontsize=12)
    ax[1, 0].set(xlabel="SB(10 px) / sky uncertainty", ylabel="naive slope", xscale="log")
    r, p, n = corr.get("all:snr10", (np.nan, np.nan, 0))
    ax[1, 0].set_title(f"slope vs outer-annulus S/N   (Spearman {r:+.2f})", fontsize=12)
    # coverage bars: km range each frame's fit spans
    s = g.sort_values("delta")
    ax[1, 1].fill_between(s["delta"], s["rho_fit_rmin_km"], s["rho_rmax_km"], color="0.8", step="mid",
                          label="fit window [1.5 FWHM, 10 px] in km")
    ax[1, 1].plot(s["delta"], s["fwhm_km"], color="k", lw=1, label="FWHM (km)")
    ax[1, 1].set(xlabel=r"$\Delta$ (au)", ylabel=r"$\rho$ (km)", yscale="log")
    ax[1, 1].legend(fontsize=10, frameon=False)
    fig.suptitle("Does the naive slope track the physical window?  (clean frames)")
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)


def fig_validity(by_km, by_px, rho_max, outpath):
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    for leg, s in by_km.groupby("leg"):
        ax[0].plot(s["rho_km"], s["median"], "o-", label=f"{leg}-perihelion", ms=5)
    ax[0].axhline(1, color="tab:green", ls="-.", label=r"$1/\rho$")
    ax[0].axhline(0.8, color="0.5", ls=":", lw=1)
    if np.isfinite(rho_max):
        ax[0].axvline(rho_max, color="crimson", ls="--", label=rf"$\rho_\max$ ≈ {rho_max:.0f} km")
    ax[0].set(xscale="log", xlabel=r"$\rho$ (km)", ylabel=r"$\rho\,$SB / plateau",
              title="in physical units: departure tied to km?", ylim=(0.2, 1.5))
    ax[0].legend(fontsize=10, frameon=False)
    for db, s in by_px.groupby("delta_bin", observed=True):
        ax[1].plot(s["r_pix"], s["median"], "o-", ms=4, label=f"$\\Delta$ {db}")
    ax[1].axhline(1, color="tab:green", ls="-.")
    ax[1].axhline(0.8, color="0.5", ls=":", lw=1)
    ax[1].set(xscale="log", xlabel="radius (px)", ylabel=r"$\rho\,$SB / plateau",
              title="in pixels: departure tied to pixel radius?", ylim=(0.2, 1.5))
    ax[1].legend(fontsize=10, frameon=False)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)


def fig_correction(refit, outpath):
    g = refit[refit["clean"]]
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))
    x = g["delta"]
    ax[0].scatter(x, g["slope_naive"], s=14, c="0.5", label="naive [1.5 FWHM, 10 px]")
    ax[0].scatter(x, g["slope_kmwin"], s=14, c="tab:blue", label="fixed km window")
    ax[0].scatter(x, g["slope_model"], s=14, c="crimson", label="PSF-convolved model")
    ax[0].scatter(x, g["slope_model_sky"], s=14, marker="^", c="tab:green",
                  label="PSF model + free sky")
    ax[0].axhline(-1, color="tab:green", ls="-.")
    ax[0].set(xlabel=r"$\Delta$ (au)", ylabel="coma slope", ylim=(-3.5, 0.5),
              title="three estimates vs observer distance")
    ax[0].legend(fontsize=9, frameon=False)
    ax[1].scatter(g["slope_naive"], g["slope_native"], s=12, c="k")
    lim = (-4, 1); ax[1].plot(lim, lim, "--", color="0.6")
    ax[1].set(xlabel="slope, oversampled", ylabel="slope, native pixels", xlim=lim, ylim=lim,
              title="oversampling changes nothing")
    ax[2].scatter(g["sb10_over_skyerr"], g["slope_naive"], s=12, c="0.5", label="naive")
    ax[2].scatter(g["sb10_over_skyerr"], g["slope_model_sky"], s=12, c="tab:green", marker="^",
                  label="model + free sky")
    ax[2].axhline(-1, color="tab:green", ls="-.")
    ax[2].set(xlabel="SB(10 px) / sky uncertainty", ylabel="coma slope", xscale="log",
              ylim=(-3.5, 0.5), title="steepness is an outer-annulus S/N effect")
    ax[2].legend(fontsize=9, frameon=False)
    fig.tight_layout(); fig.savefig(outpath, dpi=150); plt.close(fig)


# ----------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--target", default="24P")
    args = ap.parse_args(argv)
    target = zc.target_slug(args.target)
    figdir = zc.fig_dir("profile"); resdir = zc.result_dir("profile")

    prof, summ = load(target)
    corr = hypothesis_tests(summ)
    print("=== A. Spearman correlations of naive slope (clean frames) ===")
    for k, (r, p, n) in corr.items():
        print(f"  {k:14s} rho={r:+.2f}  p={p:.1e}  n={n}")

    p, by_km, by_px, rho_max = rho_sb_diagnostic(prof, summ)
    print(f"\n=== B. 1/rho validity ===\n  rho_max (post-perihelion, rho*SB >= 0.8 plateau): {rho_max:.0f} km")
    print("  binned rho*SB by leg (km):")
    for leg, s in by_km.groupby("leg"):
        print(f"    {leg:4s} " + "  ".join(f"{r:.0f}:{m:.2f}" for r, m in zip(s["rho_km"], s["median"])))
    print("  binned rho*SB by delta tercile (px):")
    for db, s in by_px.groupby("delta_bin", observed=True):
        print(f"    {str(db):4s} " + "  ".join(f"{r:.1f}:{m:.2f}" for r, m in zip(s["r_pix"], s["median"])))

    refit = refit_methods(prof, summ, rho_max)
    refit.to_csv(zc.result_path("profile", "resolution.csv", target), index=False)
    g = refit[refit["clean"]]
    print("\n=== C/D/E. slope by method (clean frames, median [16-84%]) ===")
    for c in ("slope_naive", "slope_model", "slope_model_sky"):
        pre = g[g["leg"] == "pre"]
        print(f"  pre-perihelion frames within [-1.3,-0.7] by {c:16s}: "
              f"{int(pre[c].between(-1.3, -0.7).sum())}/{len(pre)}")
    for leg in ("pre", "post"):
        s = g[g["leg"] == leg]
        def q(c): 
            v = s[c].dropna(); return (f"{v.median():+.2f} [{v.quantile(.16):+.2f},{v.quantile(.84):+.2f}] n={len(v)}"
                                       if len(v) else "n=0")
        print(f"  {leg}: naive {q('slope_naive')} | km-window {q('slope_kmwin')} | "
              f"model {q('slope_model')} | model+sky {q('slope_model_sky')} | native {q('slope_native')}")
        so = s["sky_offset_over_err"].dropna()
        if len(so):
            print(f"        fitted sky offset / sky uncertainty: median {so.median():+.2f}, "
                  f"16-84% [{so.quantile(.16):+.2f}, {so.quantile(.84):+.2f}]")
        nf = s["nucleus_fraction"].dropna()
        if len(nf): print(f"        nucleus fraction median {nf.median():.2f}; "
                          f"frames with chi2(m=1) > 2x chi2(free): {int((s['chi2_fixed_m1'] > 2*s['model_chi2']).sum())}/{len(s)}")
    far = g[g["delta"] > g["delta"].quantile(0.67)]
    print(f"  far third (delta > {g['delta'].quantile(.67):.2f} au): naive {far['slope_naive'].median():+.2f}, "
          f"km-window {far['slope_kmwin'].median():+.2f} (n={far['slope_kmwin'].notna().sum()}), "
          f"model {far['slope_model'].median():+.2f}")
    d = (g["slope_native"] - g["slope_naive"]).dropna()
    print(f"  native - oversampled slope: median {d.median():+.3f}, 16-84% [{d.quantile(.16):+.2f}, {d.quantile(.84):+.2f}]")

    feas = survey_feasibility(rho_max)
    feas.to_csv(zc.result_path("profile", "feasibility.csv"), index=False)
    print(f"\n=== F. survey feasibility with rho_max = {rho_max:.0f} km ===")
    print(feas.sort_values("delta_median").to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    fig_hypothesis(summ, corr, zc.fig_path("profile", "resolution.png", target))
    fig_validity(by_km, by_px, rho_max, zc.fig_path("profile", "validity.png", target))
    fig_correction(refit, zc.fig_path("profile", "correction.png", target))
    print("\nfigures written to", figdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
