#!/usr/bin/env python
"""Is an Af-rho trend that reverses with distance real, or an aperture effect?

C/2024 E1 rises from 4.7 to ~3.4 au and then fades toward perihelion in the
10,000 km aperture, while the 40,000 km aperture declines monotonically.  The
ratio Afrho(10k)/Afrho(40k) climbs from 1.3 at 2 au to 2.2 at 4.5 au.  For a
1/rho coma that ratio is 1 at every distance, so something in the photometry
depends on Delta.  Four candidates are tested from the tables alone:

1. sky over-subtraction of coma light in the annulus -- the annulus runs
   from 3 rho to 4 rho + 20 px with a median estimator, so for a 1/rho coma
   it removes 7.8% of the coma flux at 10k and 11.8% at 40k: a 4%
   differential, in the direction of lowering the ratio;
2. the point-source aperture correction applied to an extended source --
   tabulated per row (`apcor`), a few per cent at most;
3. resolution -- rho/FWHM per aperture, and the PSF-corrected profile slope
   from the profile study for comparison;
4. the sky *level*: the aperture sum is sky-dominated at large aperture and
   large Delta, so a systematic error of a per cent in the sky level removes
   tens of per cent of the 40k flux and almost nothing from 10k.

Outputs
  fig/ztf/afrho/systematics_<target>.png
  results/ztf/afrho/systematics_<target>.csv     per r_h bin: medians of every quantity above
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import ztfcomet as zc
from ztfcomet import rcparams  # noqa: F401

KM_PER_ARCSEC_AU = 725.27


def curve_of_growth_slope(w, e):
    """Per frame: m_eff = 1 - dlog Afrho / dlog rho across the apertures present."""
    rows = {}
    for f, row in w.iterrows():
        ok = row.notna() & (row > 0)
        if ok.sum() < 3:
            continue
        x = np.log10(row.index[ok].astype(float)); y = np.log10(row[ok].values)
        wt = 1.0 / np.maximum(e.loc[f][ok].values / (row[ok].values * np.log(10)), 1e-3) ** 2
        A = np.column_stack([np.ones_like(x), x]); sw = np.sqrt(wt)
        c, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
        rows[f] = 1.0 - c[1]
    return pd.Series(rows, name="m_eff")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--target", default="2024E1")
    ap.add_argument("--band", default="ZTF_r")
    args = ap.parse_args(argv)
    t = zc.target_slug(args.target)
    d = pd.read_csv(zc.photometry_path(t, create=False), dtype={"target": str})
    d = d[(d["filter"] == args.band) & d["quality_ok"].astype(bool) & np.isfinite(d["afrho0_cm"]) & (d["afrho0_cm"] > 0)]
    rhos = sorted(d["rho_km"].unique())
    edges = np.quantile(d["r"], np.linspace(0, 1, 7))
    d["bin"] = pd.cut(d["r"], edges, include_lowest=True)
    w = d.pivot_table(index="file", columns="rho_km", values="afrho0_cm")
    e = d.pivot_table(index="file", columns="rho_km", values="afrho0_cm_err")
    meta = d.drop_duplicates("file").set_index("file")[["obsjd", "r", "delta", "fwhm_pix", "pixscale", "bin"]]
    meff = curve_of_growth_slope(w, e)
    fr = meta.join(meff).join(w[rhos[0]].rename("af_small")).join(w[rhos[-1]].rename("af_large"))
    fr["ratio"] = fr["af_small"] / fr["af_large"]

    # per-aperture sky and aperture-correction bookkeeping
    per = []
    for rho in rhos:
        s = d[np.isclose(d["rho_km"], rho)].copy()
        npix = np.pi * s["rho_pix"] ** 2
        s["sky_over_flux"] = s["msky"] * npix / s["source_sum"]
        s["skyerr_over_flux"] = (s["ssky"] / np.sqrt(s["nsky"].clip(lower=1))) * npix / s["source_sum"]
        s["rho_fwhm"] = s["rho_pix"] / s["fwhm_pix"]
        g = s.groupby("bin", observed=True)
        per.append(pd.DataFrame({"rho_km": rho, "afrho": g["afrho0_cm"].median(), "sky_over_flux": g["sky_over_flux"].median(),
                                 "skyerr_over_flux": g["skyerr_over_flux"].median(), "apcor_mag": g["apcor"].median(),
                                 "rho_fwhm": g["rho_fwhm"].median(), "n": g.size()}).reset_index())
    per = pd.concat(per, ignore_index=True)
    binned = fr.groupby("bin", observed=True).agg(r=("r", "median"), delta=("delta", "median"), ratio=("ratio", "median"),
                                                 m_eff=("m_eff", "median"), m16=("m_eff", lambda v: v.quantile(.16)),
                                                 m84=("m_eff", lambda v: v.quantile(.84)), n=("m_eff", "size")).reset_index()
    # PSF-model profile slope from the profile study, if present
    prof = pd.DataFrame()
    sf = zc.result_path("profile", "survey_fits.csv", create=False)
    if sf.exists():
        p = pd.read_csv(sf, dtype={"target": str}); p = p[p["target"] == t]
        p = p[p["ok"].astype(str).str.lower().eq("true") & ~p["at_bound"].astype(str).str.lower().eq("true") & (p["m_free"] < 2.95)]
        if len(p):
            p["bin"] = pd.cut(p["r"], edges, include_lowest=True)
            prof = p.groupby("bin", observed=True).agg(m_psf=("m_free", "median"), n_psf=("m_free", "size"),
                                                       snr10=("snr10", "median")).reset_index()
            binned = binned.merge(prof, on="bin", how="left")
    out = binned.copy(); out["bin"] = out["bin"].astype(str)
    out.to_csv(zc.result_path("afrho", f"systematics_{t}.csv"), index=False)

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    cmap = plt.get_cmap("viridis")
    for k, rho in enumerate(rhos):
        s = d[np.isclose(d["rho_km"], rho)].sort_values("r")
        ax[0, 0].errorbar(s["r"], s["afrho0_cm"], yerr=s["afrho0_cm_err"], fmt="o", ms=4, lw=.6, capsize=0,
                          color=cmap(k / max(len(rhos) - 1, 1)), alpha=.8, label=f"{int(rho) // 1000}k km")
    ax[0, 0].set(xscale="log", yscale="log", xlabel=r"$r_h$ (au)", ylabel=r"$A(0°)f\rho$ (cm)",
                 title=f"{t}: the same comet through five apertures")
    ax[0, 0].legend(fontsize=9, frameon=False)
    ax[0, 1].plot(binned["r"], binned["ratio"], "ko-", label=rf"$Af\rho$({int(rhos[0]) // 1000}k) / $Af\rho$({int(rhos[-1]) // 1000}k)")
    ax[0, 1].axhline(1, color="tab:green", ls="-.", label=r"$1/\rho$ coma")
    axm = ax[0, 1].twinx()
    axm.plot(binned["r"], binned["m_eff"], "s-", color="tab:red", label=r"$m_{\rm eff}$ from the curve of growth")
    axm.fill_between(binned["r"], binned["m16"], binned["m84"], color="tab:red", alpha=.15)
    if "m_psf" in binned:
        axm.plot(binned["r"], binned["m_psf"], "^--", color="tab:blue", label=r"$m$ from the PSF-modelled profile (≤10 px)")
    axm.set_ylabel("coma slope m"); axm.set_ylim(0.5, 2.0)
    h1, l1 = ax[0, 1].get_legend_handles_labels(); h2, l2 = axm.get_legend_handles_labels()
    ax[0, 1].legend(h1 + h2, l1 + l2, fontsize=8, frameon=False, loc="upper left")
    ax[0, 1].set(xscale="log", xlabel=r"$r_h$ (au)", ylabel="aperture ratio", title="is the ratio the coma, or the photometry?")
    for k, rho in enumerate(rhos):
        s = per[per["rho_km"] == rho]
        rr = binned.set_index("bin").loc[s["bin"], "r"].values
        ax[1, 0].plot(rr, s["sky_over_flux"], "o-", color=cmap(k / max(len(rhos) - 1, 1)), label=f"{int(rho) // 1000}k km")
    ax[1, 0].axhline(10, color="0.6", ls=":", lw=1)
    ax[1, 0].set(xscale="log", yscale="log", xlabel=r"$r_h$ (au)", ylabel="sky inside the aperture / source flux",
                 title="a 1% sky error costs (this) % of the flux")
    ax[1, 0].legend(fontsize=9, frameon=False)
    for k, rho in enumerate((rhos[0], rhos[-1])):
        s = per[per["rho_km"] == rho]
        rr = binned.set_index("bin").loc[s["bin"], "r"].values
        ax[1, 1].plot(rr, -s["apcor_mag"], "o-" if k == 0 else "s--", color="tab:purple" if k == 0 else "0.5",
                      label=f"aperture correction, {int(rho) // 1000}k (mag)")
        ax[1, 1].plot(rr, s["rho_fwhm"] / 10, "^:" if k == 0 else "v:", color="tab:brown" if k == 0 else "0.5",
                      label=f"ρ / FWHM ÷ 10, {int(rho) // 1000}k")
    ax[1, 1].set(xscale="log", xlabel=r"$r_h$ (au)", title="resolution and aperture correction (small effects)")
    ax[1, 1].legend(fontsize=8, frameon=False)
    for a in ax.flat:
        for f in (a.xaxis.set_major_formatter, a.xaxis.set_minor_formatter):
            f(matplotlib.ticker.ScalarFormatter())
        a.xaxis.get_major_formatter().set_scientific(False); a.xaxis.get_minor_formatter().set_scientific(False)
    fig.tight_layout(); fig.savefig(zc.fig_path("afrho", f"systematics_{t}.png"), dpi=130); plt.close(fig)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))
    print("\nper-aperture sky/flux by bin:")
    print(per.pivot(index="bin", columns="rho_km", values="sky_over_flux").round(1).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
