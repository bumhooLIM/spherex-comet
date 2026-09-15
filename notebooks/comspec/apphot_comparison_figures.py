"""Figures for the previous-vs-revised photometry comparison (reads results/apphot_comparison/)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent           # notebooks/comspec/
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))     # rcparams.py
import rcparams  # noqa: E402,F401
import matplotlib.pyplot as plt  # noqa: E402

RES = ROOT / "results" / "comspec" / "studies" / "apphot_comparison"
FIG = ROOT / "fig" / "comspec" / "studies" / "apphot_comparison"
FIG.mkdir(parents=True, exist_ok=True)
R = pd.read_csv(RES / "matched_rows.csv")
B = pd.read_csv(RES / "band_emission_runs.csv")
Q = pd.read_csv(RES / "Q_runs.csv")
T = pd.read_csv(RES / "per_target.csv")
R["ratio"] = R.flux_new / R.flux_old
R["dsig"] = (R.flux_new - R.flux_old) / R.err_old
CLS = ["clean", "star-masked (old only)", "bad pixels (both)", "bad (new only)"]
COL = {"clean": "tab:blue", "star-masked (old only)": "tab:red", "bad pixels (both)": "tab:orange", "bad (new only)": "tab:green"}
BANDS = ["2.7um", "4.3um", "4.7um"]
BCOL = {"2.7um": "tab:blue", "4.3um": "tab:orange", "4.7um": "tab:green"}
ok = R.flux_old > 3 * R.err_old


def running_median(x, y, nb=8):
    x, y = np.asarray(x), np.asarray(y)
    edges = np.quantile(x, np.linspace(0, 1, nb + 1))
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x <= hi)
        if m.sum() >= 5:
            xs.append(np.median(x[m])); ys.append(np.median(y[m]))
    return np.array(xs), np.array(ys)


# ---------------------------------------------------------------- figure 1: row level
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
e = R[ok & (R.window == "emission")]
for c in CLS:
    s = e[e.cls == c]
    ax[0].scatter(s.r_ap_pix, s.ratio, s=14, alpha=0.5, color=COL[c], label=f"{c} (n={len(s)})")
xs, ys = running_median(e[e.cls == "clean"].r_ap_pix, e[e.cls == "clean"].ratio)
ax[0].plot(xs, ys, "k-", lw=3, label="clean: running median")
ax[0].axhline(1, color="k", ls="--", lw=1)
ax[0].set(xscale="log", yscale="log", ylim=(0.5, 5), xlabel="aperture radius [px]", ylabel="flux revised / previous",
          title="emission windows, previous flux > 3σ")
ax[0].legend(fontsize=13, loc="upper left")
s = R[R.window == "emission"]
ax[1].scatter(s.r_ap_pix, (s.sky_new - s.sky_old) * s.area_old / s.err_old, s=10, alpha=0.4, color="tab:gray")
xs, ys = running_median(s.r_ap_pix, (s.sky_new - s.sky_old) * s.area_old / s.err_old)
ax[1].plot(xs, ys, "k-", lw=3, label="running median")
ax[1].axhline(0, color="k", ls="--", lw=1)
ax[1].set(xscale="log", ylim=(-2, 4), xlabel="aperture radius [px]", ylabel="(sky$_{new}$ − sky$_{old}$) × area / σ$_{old}$",
          title="sky-level shift, revised − previous")
ax[1].legend(fontsize=14)
ax[2].hist(np.clip(s.err_new / s.err_old, 0.3, 1.5), bins=60, range=(0.3, 1.5), histtype="step", lw=3, label="formal revised / previous")
ax[2].hist(np.clip(s.err_emp / s.err_old, 0.3, 1.5), bins=60, range=(0.3, 1.5), histtype="step", lw=3, label="empirical revised / previous")
ax[2].axvline(1, color="k", ls="--", lw=1)
ax[2].set(xlabel="error ratio", ylabel="channels", title="photometric error")
ax[2].legend(fontsize=14)
fig.suptitle("Previous vs revised photometry at the rule aperture, channel level (annulus 150–300 × 10³ km → 15–20 px)", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "row_level.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- figure 2: windows per band
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
for a, b in zip(ax, BANDS):
    for w, ls in [("emission", "-"), ("continuum", "--")]:
        for c in ["clean", "star-masked (old only)"]:
            s = R[ok & (R.window == w) & (R.band == b) & (R.cls == c)]
            if len(s) < 5:
                continue
            a.hist(np.clip(np.log10(s.ratio), -0.5, 1.0), bins=45, range=(-0.5, 1.0), histtype="step", lw=3, ls=ls,
                   color=COL[c], label=f"{w}, {c} (n={len(s)}, med {np.median(s.ratio):.2f})")
    a.axvline(0, color="k", ls=":", lw=1)
    a.set(xlabel="log$_{10}$(flux revised / previous)", ylabel="channels", title=f"{b} band")
    a.legend(fontsize=12)
fig.suptitle("Emission versus continuum windows: the change is a pixel-treatment effect, not a band effect", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "windows_per_band.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- figure 3: band emission, same rows
key = ["target", "phase", "band"]
a_ = B[B.run == "a_old_matched"].set_index(key)
b_ = B[B.run == "b_new_matched"].set_index(key)
m = a_.join(b_, lsuffix="_1", rsuffix="_2", how="inner").reset_index()
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
for a, b in zip(ax, BANDS):
    s = m[m.band == b]
    acc = s.verdict_1.isin(["PASS", "WARN"]) & s.verdict_2.isin(["PASS", "WARN"])
    for sel, col, lab, mk in [(acc, BCOL[b], "accepted in both", "o"), (~acc, "tab:gray", "rejected in one", "x")]:
        ss = s[sel]
        a.errorbar(ss.band_flux_1, ss.band_flux_2, xerr=ss.band_flux_err_1, yerr=ss.band_flux_err_2, fmt=mk, color=col,
                   ms=6, alpha=0.6, lw=0.8, label=f"{lab} (n={len(ss)})")
    lim = np.nanpercentile(np.abs(s[["band_flux_1", "band_flux_2"]].values), [2, 99.5])
    lo, hi = max(lim[0], 1e-19), lim[1] * 2
    a.plot([lo, hi], [lo, hi], "k--", lw=1)
    sig = s[acc & (s.band_flux_1 > 2 * s.band_flux_err_1)]
    a.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi), xlabel="band flux, previous photometry [W m$^{-2}$]",
          ylabel="band flux, revised photometry [W m$^{-2}$]",
          title=f"{b}\nmedian ratio {np.nanmedian(sig.band_flux_2 / sig.band_flux_1):.3f} (n={len(sig)}, previous > 2σ)")
    a.legend(fontsize=13, loc="upper left")
fig.suptitle("Continuum-subtracted band emission from identical channel sets (pipeline method, physical space)", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "band_emission_same_rows.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- figure 4: Q
kq = ["target", "phase", "species"]
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
for a, sp in zip(ax, ["H2O", "CO2", "CO"]):
    for (r1, r2), lab, mk, fc in [(("a_old_matched", "b_new_matched"), "same rows", "o", None),
                                  (("c_old_policy", "d_new_policy"), "each pipeline's own selection", "s", "none")]:
        x = Q[(Q.run == r1) & (Q.species == sp)].set_index(kq)
        y = Q[(Q.run == r2) & (Q.species == sp)].set_index(kq)
        j = x.join(y, lsuffix="_1", rsuffix="_2", how="inner")
        j = j[(j.status_1 == "detected") & (j.status_2 == "detected") & (j.n_eff_1 >= 2) & (j.n_eff_2 >= 2)]
        a.errorbar(j.Q_fit_1, j.Q_fit_2, xerr=j.Q_err_1, yerr=j.Q_err_2, fmt=mk, mfc=fc, ms=7, alpha=0.7, lw=0.8,
                   label=f"{lab} (n={len(j)}, median ratio {np.median(j.Q_fit_2 / j.Q_fit_1):.2f})")
    lim = a.get_xlim(); lo = min(lim[0], a.get_ylim()[0]); hi = max(lim[1], a.get_ylim()[1])
    a.plot([lo, hi], [lo, hi], "k--", lw=1)
    a.set(xscale="log", yscale="log", xlabel="Q previous photometry [s$^{-1}$]", ylabel="Q revised photometry [s$^{-1}$]",
          title=f"{sp}: robust in both (detected, n$_{{eff}}$ ≥ 2)")
    a.legend(fontsize=13, loc="upper left")
fig.suptitle("Production rates from the two photometry sets", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "Q_comparison.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- figure 5: per target
acc = m[m.verdict_1.isin(["PASS", "WARN"]) & m.verdict_2.isin(["PASS", "WARN"]) & (m.band_flux_1 > 2 * m.band_flux_err_1)].copy()
acc["fr"] = acc.band_flux_2 / acc.band_flux_1
pt = acc.groupby("target").agg(n=("fr", "size"), fr=("fr", "median")).join(T.set_index("target"))
fig, ax = plt.subplots(1, 3, figsize=(24, 7))
sc = ax[0].scatter(pt.r_hel, pt.fr, s=40 + 600 * pt.frac_star_masked, c=pt.r_ap_pix, cmap="viridis", norm=plt.matplotlib.colors.LogNorm(), alpha=0.8, edgecolor="k")
for t, r in pt.iterrows():
    if abs(np.log(r.fr)) > 0.15:
        ax[0].annotate(t, (r.r_hel, r.fr), fontsize=12, xytext=(4, 4), textcoords="offset points")
ax[0].axhline(1, color="k", ls="--", lw=1)
ax[0].set(xscale="log", yscale="log", xlabel="r$_h$ [au]", ylabel="median band-flux ratio revised / previous",
          title="per target\n(size ∝ star-masked fraction)")
plt.colorbar(sc, ax=ax[0], label="aperture radius [px]")
ax[1].scatter(T.frac_star_masked, T.frac_bad_both, s=40, alpha=0.7)
for _, r in T.iterrows():
    if r.frac_star_masked > 0.3 or r.frac_bad_both > 0.3:
        ax[1].annotate(r.target, (r.frac_star_masked, r.frac_bad_both), fontsize=12, xytext=(4, 4), textcoords="offset points")
ax[1].set(xlabel="star-masked fraction, previous set", ylabel="bad-pixel fraction, both sets",
          title="why targets differ")
ax[2].scatter(T.n_old_policy, T.n_new_policy, s=40, alpha=0.7)
lim = (0, max(T.n_old_policy.max(), T.n_new_policy.max()) * 1.05)
ax[2].plot(lim, lim, "k--", lw=1)
ax[2].set(xlim=lim, ylim=lim, xlabel="previous set, previous policy", ylabel="revised set, dc_main policy",
          title="channels kept per target\n(each pipeline's own selection)")
fig.suptitle("Per-target view", y=1.02)
fig.tight_layout(); fig.savefig(FIG / "per_target.png", bbox_inches="tight"); plt.close(fig)
print("figures written:", sorted(p.name for p in FIG.glob("*.png")))
