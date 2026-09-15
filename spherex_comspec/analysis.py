"""
Cross-variant studies: source-flag contamination and the distance correction.

Every function here reads the products of :func:`pipeline.run_variant` and
never refits.  Two questions are answered:

1. **Contamination.**  What does excluding flag ``a`` (bright blends) or flag
   ``b`` (flux-ratio blends) do to the continuum verdicts, to the number of
   fittable groups, and to the retrieved Q of the groups that survive both
   policies?  A policy that changes Q for the groups it *keeps* is telling you
   the contamination was there; one that only shrinks the sample is telling
   you the cut was uninformative.
2. **Distance correction.**  Fitting the continuum in ``F x r_h^2 Delta^2``
   space rather than physical space changes the polynomial and therefore the
   emission; the retrieved Q are compared group by group.  Q itself is always
   retrieved in physical space (see :mod:`fitting`), so any difference is
   entirely a continuum-placement effect.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import directory as _dir
from .config import EMISSION_DTYPES, PLACEHOLDERS, SPECIES, Variant
from .logging_utils import get_logger

__all__ = ["load_fits", "load_continuum", "variant_census", "paired_Q", "compare_variants",
           "distcorr_effect", "placeholder_table", "write_analysis"]

log = get_logger("analysis")

_ROBUST_NEFF = 2.0


def load_fits(variant: str) -> pd.DataFrame:
    p = _dir.variant_dirs(variant)["results"] / "gas_fit.csv"
    return pd.read_csv(p, dtype={"target": str}) if p.is_file() else pd.DataFrame()


def load_continuum(variant: str) -> pd.DataFrame:
    p = _dir.variant_dirs(variant)["results"] / "continuum_summary.csv"
    return pd.read_csv(p, dtype=EMISSION_DTYPES) if p.is_file() else pd.DataFrame()


def variant_census(variant: str) -> dict:
    """Sample size, verdicts and detection counts of one variant."""
    f, c = load_fits(variant), load_continuum(variant)
    out = dict(variant=variant, n_groups=int(c.groupby(["target", "phase"]).ngroups) if len(c) else 0,
               n_band_rows=len(c), n_fitted=len(f))
    if len(c):
        vc = c.verdict.value_counts()
        for v in ("PASS", "WARN", "FAIL"):
            out[f"n_{v}"] = int(vc.get(v, 0))
        out["frac_saved"] = float((vc.get("PASS", 0) + vc.get("WARN", 0)) / len(c))
        out["median_err_scale"] = float(c.err_scale.median())
    for s in SPECIES:
        if len(f):
            det = f[f[f"Q_{s}_status"] == "detected"]
            out[f"n_det_{s}"] = len(det)
            out[f"n_robust_{s}"] = int((det[f"Q_{s}_n_eff"] >= _ROBUST_NEFF).sum())
            out[f"n_marginal_{s}"] = int((f[f"Q_{s}_status"] == "marginal").sum())
            out[f"n_rejected_{s}"] = int((f[f"Q_{s}_status"] == "rejected").sum())
            out[f"n_ul_{s}"] = int(f[f"Q_{s}_status"].isin(["upper_limit", "negative_fit"]).sum())
            out[f"median_relerr_{s}"] = float((det[f"Q_{s}_err"] / det[f"Q_{s}"]).median()) if len(det) else np.nan
            if s == "H2O" and "h2o_source" in det:
                out["n_det_H2O_hot"] = int((det.h2o_source == "hot").sum())
                out["n_robust_H2O_hot"] = int(((det.h2o_source == "hot") & (det.Q_H2O_n_eff >= _ROBUST_NEFF)).sum())
        else:
            out.update({f"n_det_{s}": 0, f"n_robust_{s}": 0, f"n_ul_{s}": 0, f"median_relerr_{s}": np.nan})
    return out


def paired_Q(ref: str, other: str, species: Sequence[str] = SPECIES) -> pd.DataFrame:
    """
    Q of the groups fitted under *both* variants, side by side.

    Returns one row per (target, phase, species) with the two values, their
    errors, the ratio ``Q_other / Q_ref`` and its significance
    ``(Q_other - Q_ref) / sqrt(err_ref^2 + err_other^2)``.
    """
    a, b = load_fits(ref), load_fits(other)
    if a.empty or b.empty:
        return pd.DataFrame()
    m = a.merge(b, on=["target", "phase"], suffixes=("_ref", "_oth"))
    rows = []
    for _, r in m.iterrows():
        for s in species:
            qa, qb = r[f"Q_{s}_fit_ref"], r[f"Q_{s}_fit_oth"]
            ea, eb = r[f"Q_{s}_err_ref"], r[f"Q_{s}_err_oth"]
            if not (np.isfinite(qa) and np.isfinite(qb)):
                continue
            rows.append(dict(
                target=r.target, phase=int(r.phase), species=s,
                Q_ref=qa, Q_oth=qb, err_ref=ea, err_oth=eb,
                status_ref=r[f"Q_{s}_status_ref"], status_oth=r[f"Q_{s}_status_oth"],
                ratio=qb / qa if qa != 0 else np.nan,
                nsig=(qb - qa) / np.hypot(ea, eb) if np.isfinite(ea) and np.isfinite(eb) else np.nan,
                n_eff_ref=r[f"Q_{s}_n_eff_ref"], n_eff_oth=r[f"Q_{s}_n_eff_oth"],
                r_hel=r.r_hel_mean_ref, n_flag_a=r.get("n_flag_a_ref", np.nan),
                n_flag_b=r.get("n_flag_b_ref", np.nan)))
    return pd.DataFrame(rows)


def _pair_summary(pq: pd.DataFrame, label: str) -> List[dict]:
    """Median ratio and significance per species, restricted to genuine detections in both."""
    out = []
    for s in SPECIES:
        d = pq[(pq.species == s) & (pq.status_ref == "detected") & (pq.status_oth == "detected")]
        both = pq[pq.species == s]
        out.append(dict(comparison=label, species=s, n_common=len(both), n_both_detected=len(d),
                        n_lost=int(((both.status_ref == "detected") & (both.status_oth != "detected")).sum()),
                        n_gained=int(((both.status_ref != "detected") & (both.status_oth == "detected")).sum()),
                        median_ratio=float(d.ratio.median()) if len(d) else np.nan,
                        p16_ratio=float(d.ratio.quantile(0.16)) if len(d) else np.nan,
                        p84_ratio=float(d.ratio.quantile(0.84)) if len(d) else np.nan,
                        frac_changed_gt1sig=float((d.nsig.abs() > 1).mean()) if len(d) else np.nan,
                        frac_changed_gt3sig=float((d.nsig.abs() > 3).mean()) if len(d) else np.nan,
                        max_abs_nsig=float(d.nsig.abs().max()) if len(d) else np.nan))
    return out


def compare_variants(ref: str, others: Sequence[str]) -> Dict[str, pd.DataFrame]:
    """
    Flag-policy comparison against a reference variant.

    Returns
    -------
    dict
        ``census`` -- one row per variant; ``pairs`` -- per-species medians per
        comparison; ``paired`` -- the full per-group table.
    """
    census = pd.DataFrame([variant_census(v) for v in [ref, *others]])
    pairs, paired = [], []
    for o in others:
        pq = paired_Q(ref, o)
        if pq.empty:
            continue
        pq.insert(0, "comparison", f"{o} vs {ref}")
        paired.append(pq)
        pairs.extend(_pair_summary(pq, f"{o} vs {ref}"))
    return dict(census=census, pairs=pd.DataFrame(pairs),
                paired=pd.concat(paired, ignore_index=True) if paired else pd.DataFrame())


def distcorr_effect(dc: str, raw: str) -> Dict[str, pd.DataFrame]:
    """
    What the distance-corrected continuum changes, relative to the physical-space run.

    Two levels: the continuum verdicts of every band-row present in both runs,
    and the retrieved Q of every group fitted in both.  Also relates the Q
    change to the spread of ``r_h^2 Delta^2`` inside the group, which is the
    quantity the correction removes.
    """
    ca, cb = load_continuum(raw), load_continuum(dc)
    verdicts = pd.DataFrame()
    if len(ca) and len(cb):
        m = ca.merge(cb, on=["target", "phase", "band"], suffixes=("_raw", "_dc"))
        verdicts = pd.crosstab(m.verdict_raw, m.verdict_dc).reindex(
            index=["PASS", "WARN", "FAIL"], columns=["PASS", "WARN", "FAIL"], fill_value=0)
        m["snr_ratio"] = m.snr_peak_dc / m.snr_peak_raw
        m["dc_spread"] = m.distcorr_factor_mean_dc  # placeholder column for the join below
    pq = paired_Q(raw, dc)
    summary = pd.DataFrame(_pair_summary(pq, f"{dc} vs {raw}")) if len(pq) else pd.DataFrame()
    if len(pq):
        # Spread of r_h^2 Delta^2 across the *channels* of each group -- the quantity the
        # correction removes -- read from the fitted-channel files rather than per-band means.
        pq = pq.merge(_group_dc_spread(dc), on=["target", "phase"], how="left")
    return dict(verdicts=verdicts, summary=summary, paired=pq,
                band_rows=m[["target", "phase", "band", "verdict_raw", "verdict_dc",
                             "snr_peak_raw", "snr_peak_dc", "snr_ratio"]] if len(ca) and len(cb) else pd.DataFrame())


def _group_dc_spread(variant: str) -> pd.DataFrame:
    """``max/min - 1`` of ``distcorr_factor`` over the channels that entered each group's fit."""
    d = _dir.variant_dirs(variant)["lines"]
    rows = []
    for f in sorted(d.glob("*_points.csv")):
        p = pd.read_csv(f, usecols=["target", "phase", "distcorr_factor"], dtype={"target": str})
        if len(p) and p.distcorr_factor.min() > 0:
            rows.append(dict(target=str(p.target.iloc[0]), phase=int(p.phase.iloc[0]),
                             dc_spread=float(p.distcorr_factor.max() / p.distcorr_factor.min() - 1)))
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=["target", "phase", "dc_spread"])


def placeholder_table() -> pd.DataFrame:
    return pd.DataFrame(PLACEHOLDERS).sort_values("priority").reset_index(drop=True)


def write_analysis(main: str, flag_variants: Sequence[str], raw_variant: Optional[str] = None,
                   extra_variants: Sequence[str] = ()) -> Dict[str, object]:
    """
    Run every study and write the tables under ``results/comspec/``.

    ``extra_variants`` -- the badphot-policy run and any previous baseline -- join the
    flag census and the paired-Q comparison against ``main``; every comparison is the
    same pairwise machinery, so the studies differ only in which run is paired.
    """
    _dir.ensure_dirs()
    out = {}
    others = [v for v in flag_variants if v != main] + [v for v in extra_variants if v != main]
    cmp = compare_variants(main, others)
    cmp["census"].to_csv(_dir.STUDY_RESULT_DIR / "flag_policy_census.csv", index=False)
    cmp["pairs"].to_csv(_dir.STUDY_RESULT_DIR / "flag_policy_pairs.csv", index=False)
    cmp["paired"].to_csv(_dir.STUDY_RESULT_DIR / "flag_policy_paired_Q.csv", index=False)
    out["flags"] = cmp
    if raw_variant:
        dce = distcorr_effect(main, raw_variant)
        dce["summary"].to_csv(_dir.STUDY_RESULT_DIR / "distcorr_effect_summary.csv", index=False)
        dce["paired"].to_csv(_dir.STUDY_RESULT_DIR / "distcorr_effect_paired_Q.csv", index=False)
        dce["band_rows"].to_csv(_dir.STUDY_RESULT_DIR / "distcorr_effect_band_rows.csv", index=False)
        dce["verdicts"].to_csv(_dir.STUDY_RESULT_DIR / "distcorr_effect_verdicts.csv")
        out["distcorr"] = dce
    pt = placeholder_table()
    pt.to_csv(_dir.RESULT_DIR / "placeholders.csv", index=False)
    out["placeholders"] = pt
    log.info("analysis tables written to %s", _dir.RESULT_DIR)
    return out
