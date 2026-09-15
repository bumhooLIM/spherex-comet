"""
Previous (`_archive/data_apphot_previous`) versus revised (`data/apphot`) aperture photometry
inside the gas-emission windows, and what the difference does to the continuum
subtraction and the production-rate fit.

Why this exists
---------------
The catalog's Q values were first derived from `data/apphot`; the revised photometry
changed the sky annulus (fixed 150 000-300 000 km -> fixed 15-20 px), the bad-pixel
treatment (zeroed pixels -> masked pixels with an effective-area correction), the
Gaia handling (pixels near stars masked -> rows flagged, pixels kept) and the error
model.  Every exposure exists in both sets, so the rows can be matched one to one
and each change measured where it matters: the channels inside the 2.7 / 4.3 /
4.7 um emission windows.

Four runs of the pipeline's own continuum subtraction (physical flux space) are
compared per (target, phase, band):

    a  previous values, rows kept by the revised baseline policy   } same rows:
    b  revised values,  the same rows                              } photometry only
    c  previous values, the previous policy (badphot = flux <= 0)  } each pipeline's
    d  revised values,  the revised baseline policy (dc_main)      } own selection

Outputs: results/studies/apphot_comparison/*.csv and fig/studies/apphot_comparison/*.png.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent           # notebooks/comspec/
ROOT = HERE.parents[1]                            # the project root (the packages)
sys.path.insert(0, str(ROOT / "notebooks"))     # rcparams.py
sys.path.insert(0, str(ROOT))
import rcparams  # noqa: E402,F401  (project figure style)
import matplotlib.pyplot as plt  # noqa: E402

from spherex_comspec.config import (BAND_WINDOWS, VARIANTS, ApertureConfig, ContinuumConfig,  # noqa: E402
                                    FitConfig, ModelParams)
from spherex_comspec.continuum import insufficient, process_group  # noqa: E402
from spherex_comspec.dataio import PhaseAssignment, aperture_for, load_apphot, select_spectrum  # noqa: E402
from spherex_comspec.fitting import fit_production_rates  # noqa: E402

warnings.filterwarnings("ignore")
from spherex_comspec import directory as _dir  # noqa: E402

# The previous photometry this study compared against (spherex_apphot config 5502194856bc,
# fixed 15-20 px annulus) was deleted in the 2026-09-15 cleanup; the study's products are
# kept in results/comspec/studies/apphot_comparison/.  Point OLD_DIR at a copy to rerun.
OLD_DIR = ROOT / "_archive" / "data_apphot_previous"
RES = _dir.STUDY_RESULT_DIR / "apphot_comparison"
FIG = _dir.STUDY_FIG_DIR / "apphot_comparison"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
MAIN = VARIANTS["dc_main"]
CCFG, FCFG, ACFG = ContinuumConfig(), FitConfig(), ApertureConfig()
EM = {b: w["em"] for b, w in BAND_WINDOWS.items()}
CONT = {b: w["cont"] for b, w in BAND_WINDOWS.items()}
SPECIES = ("H2O", "CO2", "CO")


def classify(wl: float):
    for b, (lo, hi) in EM.items():
        if lo <= wl <= hi:
            return b, "emission"
    for b, (lo, hi) in CONT.items():
        if lo <= wl <= hi:
            return b, "continuum"
    return "", "other"


def old_frame(o: pd.DataFrame, flags: pd.Series | None) -> pd.DataFrame:
    """A `select_spectrum`-shaped frame from previous-photometry rows (physical flux space)."""
    d = pd.DataFrame(dict(
        filename=o.filename.values, wl=o.wl.values, wlwidth=o.wlwidth.values,
        flux=o.source_sum_mjy.values, err=o.source_sum_err_mjy.values,
        flux_raw=o.source_sum_mjy.values, err_raw=o.source_sum_err_mjy.values,
        distcorr_factor=(o.r_hel ** 2 * o.r_obs ** 2).values, snr=o.snr.values,
        r_hel=o.r_hel.values, r_obs=o.r_obs.values, jd_utc=o.jd_utc.values,
        sourceflag=(flags.values if flags is not None else np.full(len(o), "0")),
        badphot=o.badphot.astype(bool).values,
        frac_badpix_ap=(o.nbadpix / o.aperture_area_pixel2).values,
        n_gaia=np.nan, gmag_eff=np.nan))
    return d


def run_groups(raw: pd.DataFrame, target: str, r_ap_km: float, run: str, band_rows: list, q_rows: list):
    """Continuum-subtract and fit every phase of one spectrum; append the band and Q rows."""
    raw = raw.copy()
    raw.attrs["rejection"] = dict(n_total=len(raw), n_nonfinite=0, n_badphot=0, n_flag=0, n_kept=len(raw))
    for ph, g in raw.groupby("phase"):
        g = g.sort_values("wl").reset_index(drop=True)
        g.attrs["rejection"] = raw.attrs["rejection"]
        if insufficient(g, CCFG) is not None:
            continue
        out = process_group(g, target, r_ap_km, int(ph), int(g.epoch.iloc[0]), CCFG, "physical")
        summ = pd.DataFrame(out["summary"]) if not isinstance(out["summary"], pd.DataFrame) else out["summary"]
        for r in summ.to_dict("records"):
            band_rows.append(dict(run=run, target=target, phase=int(ph), band=r["band"], verdict=r["verdict"],
                                  n_cont_used=r.get("n_cont_used"), n_emission=r.get("n_emission"),
                                  band_flux=r.get("band_flux_W_m2"), band_flux_err=r.get("band_flux_err_W_m2"),
                                  snr_peak=r.get("snr_peak"), f_peak=r.get("f_peak_mjy"),
                                  cont_center=r.get("cont_at_center_mjy"), chi2_red=r.get("chi2_red"),
                                  n_flag_a=r.get("n_flag_a"), n_flag_b=r.get("n_flag_b"),
                                  r_hel=r.get("r_hel_mean")))
        pts = out["points"]
        ok = summ[summ.verdict.isin(FCFG.accept_verdicts) & (summ.n_emission > 0)]
        fit_in = pts[(pts.role == "emission") & pts.band.isin(ok.band)].copy()
        fit_in = fit_in[np.isfinite(fit_in.emis_raw_mjy) & np.isfinite(fit_in.emis_raw_err_mjy) & (fit_in.emis_raw_err_mjy > 0)]
        if fit_in.empty:
            continue
        try:
            f = fit_production_rates(fit_in.reset_index(drop=True), ModelParams(rho_ap_km=r_ap_km), FCFG, "physical")
        except ValueError:
            continue
        for i, s in enumerate(f.species):
            q_rows.append(dict(run=run, target=target, phase=int(ph), species=s, status=f.status[i],
                               Q_fit=f.Q_fit[i], Q_err=f.Q_err[i], n_eff=f.n_eff[i], n_points=f.n_points,
                               chi2_red=f.chi2 / f.dof if f.dof else np.nan, r_hel=f.geometry.get("r_hel_mean", np.nan)
                               if hasattr(f, "geometry") else np.nan))


def main() -> None:
    assign = PhaseAssignment.load()
    rows, band_rows, q_rows, target_rows = [], [], [], []
    targets = sorted(p.stem for p in OLD_DIR.glob("*.csv"))
    for t in targets:
        new = load_apphot(t)
        r_ap_km, lab, cov = aperture_for(t, ACFG)
        o = pd.read_csv(OLD_DIR / f"{t}.csv")
        o = o[o.r_ap_km == r_ap_km].drop(columns=["phase", "phase_update"], errors="ignore")
        # load_apphot reads a column subset; the row-level comparison needs the full record
        n = pd.read_csv(_dir.APPHOT_DIR / f"{t}.csv")
        n = n[n.ap_label == lab].drop(columns=["epoch"], errors="ignore") \
             .merge(new[["filename", "ap_label", "jd_utc_exact"]], on=["filename", "ap_label"], how="left")
        if o.empty or n.empty:
            target_rows.append(dict(target=t, r_ap_km=r_ap_km, n_old=len(o), n_new=len(n), note="no common aperture"))
            continue
        a_t = assign[assign.target == t][["filename", "phase", "epoch"]]
        m = o.merge(n, on="filename", suffixes=("_o", "_n")).merge(a_t, on="filename", how="left")
        bw = m.wl_o.map(classify)
        m["band"] = [b for b, _ in bw]
        m["window"] = [w for _, w in bw]
        # ---- row-level record (emission + continuum windows only) -----------------------------
        keep_new = (~(m.badphot_n.astype(bool) & (m.frac_badpix_ap > MAIN.flags.max_frac_badpix))
                    & ~m.sourceflag.astype(str).isin(MAIN.flags.drop_flags)
                    & np.isfinite(m.source_sum_mjy_o) & np.isfinite(m.source_sum_mjy_n)
                    & (m.source_sum_err_mjy_o > 0) & (m.source_sum_err_mjy_n > 0))
        star = (m.nbadpix > 0) & (m.n_badpix_ap == 0)
        both = (m.nbadpix > 0) & (m.n_badpix_ap > 0)
        cls = np.where(star, "star-masked (old only)", np.where(both, "bad pixels (both)",
                       np.where(m.nbadpix == 0, "clean", "bad (new only)")))
        w = m[m.window != "other"]
        rows.append(pd.DataFrame(dict(
            target=t, r_ap_km=r_ap_km, r_ap_pix=w.r_ap_pix, r_hel=w.r_hel_o, band=w.band, window=w.window,
            wl=w.wl_o, flux_old=w.source_sum_mjy_o, flux_new=w.source_sum_mjy_n,
            err_old=w.source_sum_err_mjy_o, err_new=w.source_sum_err_mjy_n, err_emp=w.source_sum_err_empirical_mjy,
            sky_old=w.annulus_median_mjy_per_pix, sky_new=w.sky_median_mjy_per_pix,
            skystd_old=w.bkg_std_mjy_per_pix, skystd_new=w.sky_std_mjy_per_pix,
            area_old=w.aperture_area_pixel2, area_new=w.aperture_area_pix2, area_eff_new=w.aperture_area_eff_pix2,
            nbadpix_old=w.nbadpix, nbadpix_new=w.n_badpix_ap, frac_badpix_new=w.frac_badpix_ap,
            badphot_old=w.badphot_o.astype(bool), badphot_new=w.badphot_n.astype(bool),
            sourceflag=w.sourceflag.astype(str), star_dist_px=w.neargaia_dist_pixel, star_gmag=w.neargaia_gmag,
            r_in_old_km=w.r_in_pixel / w.r_ap_pixel * r_ap_km, r_in_new_km=15.0 / w.r_ap_pix * r_ap_km,
            cls=cls[m.window != "other"], kept_new_policy=keep_new[m.window != "other"], phase=w.phase)))
        # ---- the four continuum runs ---------------------------------------------------------
        mm = m[keep_new & m.phase.notna()]
        run_groups(old_frame(mm.rename(columns={c: c[:-2] for c in mm.columns if c.endswith("_o")}), mm.sourceflag)
                   .assign(phase=mm.phase.values, epoch=mm.epoch.values), t, r_ap_km, "a_old_matched", band_rows, q_rows)
        nb = pd.DataFrame(dict(filename=mm.filename.values, wl=mm.wl_n.values, wlwidth=mm.wlwidth_n.values,
                               flux=mm.source_sum_mjy_n.values, err=mm.source_sum_err_mjy_n.values,
                               flux_raw=mm.source_sum_mjy_n.values, err_raw=mm.source_sum_err_mjy_n.values,
                               distcorr_factor=mm.distcorr_factor.values, snr=mm.snr_n.values,
                               r_hel=mm.r_hel_n.values, r_obs=mm.r_obs_n.values, jd_utc=mm.jd_utc_exact.values,
                               sourceflag=mm.sourceflag.astype(str).values, badphot=mm.badphot_n.astype(bool).values,
                               frac_badpix_ap=mm.frac_badpix_ap.values, n_gaia=mm.n_gaia.values,
                               gmag_eff=mm.gmag_eff.values, phase=mm.phase.values, epoch=mm.epoch.values))
        run_groups(nb, t, r_ap_km, "b_new_matched", band_rows, q_rows)
        oc = o.merge(a_t, on="filename", how="left")
        oc = oc[~oc.badphot.astype(bool) & oc.phase.notna() & np.isfinite(oc.source_sum_mjy) & (oc.source_sum_err_mjy > 0)]
        run_groups(old_frame(oc, None).assign(phase=oc.phase.values, epoch=oc.epoch.values), t, r_ap_km,
                   "c_old_policy", band_rows, q_rows)
        nd = select_spectrum(new, lab, MAIN, a_t)
        nd = nd[nd.phase.notna()].copy()
        # dc_main fits its continuum in distance-corrected space; this comparison is in physical
        # space throughout, so take the physical columns select_spectrum carries alongside
        nd["flux"], nd["err"] = nd.flux_raw, nd.err_raw
        run_groups(nd, t, r_ap_km, "d_new_policy", band_rows, q_rows)
        target_rows.append(dict(target=t, r_ap_km=r_ap_km, ap_coverage=cov, n_old=len(o), n_new=len(n), n_matched=len(m),
                                n_kept_new_policy=int(keep_new.sum()), n_old_policy=len(oc), n_new_policy=len(nd),
                                r_hel=float(o.r_hel.median()), r_ap_pix=float(n.r_ap_pix.median()),
                                frac_star_masked=float(star.mean()), frac_bad_both=float(both.mean())))
        print(f"{t:10s} r_ap {r_ap_km:6.0f} km  matched {len(m):5d}  kept {int(keep_new.sum()):5d}  old-policy {len(oc):5d}  new-policy {len(nd):5d}")
    R = pd.concat(rows, ignore_index=True)
    R.to_csv(RES / "matched_rows.csv", index=False)
    B = pd.DataFrame(band_rows)
    B.to_csv(RES / "band_emission_runs.csv", index=False)
    Q = pd.DataFrame(q_rows)
    Q.to_csv(RES / "Q_runs.csv", index=False)
    pd.DataFrame(target_rows).to_csv(RES / "per_target.csv", index=False)
    print("rows", len(R), "| band rows", len(B), "| Q rows", len(Q))


if __name__ == "__main__":
    main()
