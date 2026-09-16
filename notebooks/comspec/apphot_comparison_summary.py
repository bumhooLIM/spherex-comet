"""Summary tables for the previous-vs-revised photometry comparison (reads results/apphot_comparison/)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "comspec" / "studies" / "apphot_comparison"
R = pd.read_csv(RES / "matched_rows.csv")
B = pd.read_csv(RES / "band_emission_runs.csv")
Q = pd.read_csv(RES / "Q_runs.csv")
T = pd.read_csv(RES / "per_target.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)


def q(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return dict(n=0, med=np.nan, p16=np.nan, p84=np.nan)
    return dict(n=len(x), med=np.median(x), p16=np.percentile(x, 16), p84=np.percentile(x, 84))


R["ratio"] = R.flux_new / R.flux_old
R["dsig"] = (R.flux_new - R.flux_old) / R.err_old
ok = R.flux_old > 3 * R.err_old
CLS = ["clean", "star-masked (old only)", "bad pixels (both)"]

# ---- 1. channel level -------------------------------------------------------------------------
rows = []
for w in ["emission", "continuum"]:
    for c in CLS:
        s = R[ok & (R.window == w) & (R.cls == c)]
        rows.append(dict(window=w, cls=c, frac_of_rows=(R[R.window == w].cls == c).mean(),
                         **{f"ratio_{k}": v for k, v in q(s.ratio).items()},
                         **{f"dsig_{k}": v for k, v in q(s.dsig).items() if k != "n"}))
S1 = pd.DataFrame(rows).round(3); S1.to_csv(RES / "summary_row_level.csv", index=False)
print("=== channel level (rule aperture, previous flux > 3 sigma)\n", S1.to_string(index=False))
e = R[R.window == "emission"]
rows = [dict(quantity="sky shift (new-old) x area / err_old", **q((e.sky_new - e.sky_old) * e.area_old / e.err_old)),
        dict(quantity="err_new / err_old", **q(e.err_new / e.err_old)),
        dict(quantity="err_empirical_new / err_old", **q(e.err_emp / e.err_old)),
        dict(quantity="skystd_new / skystd_old", **q(e.skystd_new / e.skystd_old))]
cl = R[ok & (R.window == "emission") & (R.cls == "clean")]
for lo, hi in [(0.8, 1.5), (1.5, 3), (3, 6), (6, 13)]:
    s = cl[(cl.r_ap_pix >= lo) & (cl.r_ap_pix < hi)]
    rows.append(dict(quantity=f"clean ratio, r_ap {lo}-{hi} px (r_in old {s.r_in_old_km.median():.0f} km, new {s.r_in_new_km.median():.0f} km)", **q(s.ratio)))
S1b = pd.DataFrame(rows).round(3); S1b.to_csv(RES / "summary_row_level_extra.csv", index=False)
print("\n", S1b.to_string(index=False))

# ---- 2. band level ----------------------------------------------------------------------------
key = ["target", "phase", "band"]


def pair(r1, r2):
    a = B[B.run == r1].set_index(key); b = B[B.run == r2].set_index(key)
    return a.join(b, lsuffix="_1", rsuffix="_2", how="inner")


ab, cd = pair("a_old_matched", "b_new_matched"), pair("c_old_policy", "d_new_policy")
ct = pd.crosstab(ab.verdict_1, ab.verdict_2, rownames=["previous"], colnames=["revised"]); ct.to_csv(RES / "verdict_transitions_same_rows.csv")
print("\n=== verdicts, same rows\n", ct)
ct2 = pd.crosstab(cd.verdict_1, cd.verdict_2, rownames=["previous policy"], colnames=["revised policy"]); ct2.to_csv(RES / "verdict_transitions_own_policy.csv")
print("\n=== verdicts, each pipeline's own selection\n", ct2)
acc = ab[ab.verdict_1.isin(["PASS", "WARN"]) & ab.verdict_2.isin(["PASS", "WARN"])]
rows = []
for b in ["2.7um", "4.3um", "4.7um"]:
    s = acc[acc.index.get_level_values("band") == b]
    ds = (s.band_flux_2 - s.band_flux_1) / s.band_flux_err_1
    sig = s[s.band_flux_1 > 2 * s.band_flux_err_1]
    rows.append(dict(band=b, n_both_accepted=len(s), n_prev_flux_gt2sig=len(sig),
                     **{f"fluxratio_{k}": v for k, v in q(sig.band_flux_2 / sig.band_flux_1).items() if k != "n"},
                     **{f"dsig_{k}": v for k, v in q(ds).items() if k != "n"},
                     frac_abs_dsig_gt1=float((np.abs(ds) > 1).mean()), frac_abs_dsig_gt3=float((np.abs(ds) > 3).mean()),
                     snr_peak_ratio_med=float(np.nanmedian(s.snr_peak_2 / s.snr_peak_1))))
S2 = pd.DataFrame(rows).round(3); S2.to_csv(RES / "summary_band_same_rows.csv", index=False)
print("\n=== band emission, same rows, accepted in both\n", S2.to_string(index=False))
rows = []
for r in ["a_old_matched", "b_new_matched", "c_old_policy", "d_new_policy"]:
    s = B[B.run == r]
    rows.append(dict(run=r, groups=s.groupby(["target", "phase"]).ngroups, band_rows=len(s), PASS=int((s.verdict == "PASS").sum()),
                     WARN=int((s.verdict == "WARN").sum()), FAIL=int((s.verdict == "FAIL").sum()),
                     accepted_with_emission=int((s.verdict.isin(["PASS", "WARN"]) & (s.n_emission > 0)).sum()),
                     snr_peak_gt3=int((s.snr_peak > 3).sum())))
S2b = pd.DataFrame(rows); S2b.to_csv(RES / "summary_runs.csv", index=False); print("\n=== the four runs\n", S2b.to_string(index=False))

# ---- 3. Q level -------------------------------------------------------------------------------
kq = ["target", "phase", "species"]


def pairq(r1, r2):
    a = Q[Q.run == r1].set_index(kq); b = Q[Q.run == r2].set_index(kq)
    return a.join(b, lsuffix="_1", rsuffix="_2", how="inner")


rows, pairs = [], []
for lab, (r1, r2) in {"same_rows": ("a_old_matched", "b_new_matched"), "own_policy": ("c_old_policy", "d_new_policy")}.items():
    m = pairq(r1, r2)
    for sp in ["H2O", "CO2", "CO"]:
        s = m[m.index.get_level_values("species") == sp]
        d1 = (s.status_1 == "detected") & (s.n_eff_1 >= 2); d2 = (s.status_2 == "detected") & (s.n_eff_2 >= 2)
        both = s[d1 & d2]; ratio = both.Q_fit_2 / both.Q_fit_1
        ns = (both.Q_fit_2 - both.Q_fit_1) / np.maximum(both.Q_err_1, both.Q_err_2)
        n1 = int(((Q.run == r1) & (Q.species == sp) & (Q.status == "detected") & (Q.n_eff >= 2)).sum())
        n2 = int(((Q.run == r2) & (Q.species == sp) & (Q.status == "detected") & (Q.n_eff >= 2)).sum())
        rows.append(dict(comparison=lab, species=sp, robust_prev=n1, robust_rev=n2, n_common=len(s), n_both_robust=len(both),
                         lost=int((d1 & ~d2).sum()), gained=int((~d1 & d2).sum()),
                         **{f"Qratio_{k}": v for k, v in q(ratio).items() if k != "n"},
                         **{f"nsig_{k}": v for k, v in q(ns).items() if k != "n"},
                         max_abs_nsig=float(np.nanmax(np.abs(ns))) if len(ns) else np.nan,
                         frac_gt3sig=float((np.abs(ns) > 3).mean()) if len(ns) else np.nan))
        pp = both.reset_index()[["target", "phase", "species", "Q_fit_1", "Q_err_1", "Q_fit_2", "Q_err_2"]].copy()
        pp["comparison"] = lab; pp["ratio"] = ratio.values; pp["nsig"] = ns.values; pairs.append(pp)
S3 = pd.DataFrame(rows).round(3); S3.to_csv(RES / "summary_Q.csv", index=False); print("\n=== Q\n", S3.to_string(index=False))
P = pd.concat(pairs, ignore_index=True)
# why does Q change on identical rows?  join the group's bad-pixel / star-mask fractions (emission windows)
g = R[(R.window == "emission") & R.kept_new_policy].groupby(["target", "phase"]).agg(
    n_em=("cls", "size"), frac_star=("cls", lambda x: (x == "star-masked (old only)").mean()),
    frac_badboth=("cls", lambda x: (x == "bad pixels (both)").mean()), r_ap_pix=("r_ap_pix", "median")).reset_index()
P = P.merge(g, on=["target", "phase"], how="left"); P.round(4).to_csv(RES / "Q_pairs.csv", index=False)
s = P[P.comparison == "same_rows"]
print("\n=== same rows: Q ratio vs the group's affected-channel fractions (emission windows)")
for sp in ["H2O", "CO2", "CO"]:
    x = s[s.species == sp]
    if len(x) < 4:
        continue
    lo = x[(x.frac_star + x.frac_badboth) < 0.2]; hi = x[(x.frac_star + x.frac_badboth) >= 0.2]
    print(f"  {sp}: affected < 20 % of channels: n={len(lo)} ratio med {np.median(lo.ratio) if len(lo) else np.nan:.3f} | >= 20 %: n={len(hi)} ratio med {np.median(hi.ratio) if len(hi) else np.nan:.3f}"
          f" | Spearman(ratio, affected frac) = {x[['ratio']].assign(f=x.frac_star + x.frac_badboth).corr(method='spearman').iloc[0, 1]:.2f}")
    print("    largest changes:", x.reindex(x.nsig.abs().sort_values(ascending=False).index).head(4)[["target", "phase", "ratio", "nsig", "frac_star", "frac_badboth", "r_ap_pix"]].round(2).to_dict("records"))
