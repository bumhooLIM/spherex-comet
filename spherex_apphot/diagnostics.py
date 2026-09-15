"""
Growth-curve diagnostics and source-flag effectiveness.

The question this module answers is: *do the Gaia contamination flags actually
catch the measurements that a star has corrupted?*  A flag that is never raised
is useless; one that is raised on everything is equally useless.  Both failure
modes are invisible unless the flags are tested against an independent symptom
of contamination, which is what the growth curve provides.

The symptom
-----------
Enclosed flux against aperture radius is the growth curve.  Differentiate it and
you get the mean surface brightness of each annulus,

.. math:: \\Sigma_i = \\frac{F_i - F_{i-1}}{A_i - A_{i-1}}

For a steady-state coma the surface brightness falls monotonically with
cometocentric distance (canonically as :math:`1/\\rho`), so :math:`\\Sigma_i`
should *decrease* outward.  A field star entering the aperture at some radius
produces a local **rise** instead -- a signature that no plausible coma model
produces and that is therefore diagnostic of contamination rather than of comet
physics.

:func:`growth_metrics` measures the most significant such rise per exposure, in
units of its own uncertainty, and :func:`flag_effectiveness` scores each flag
against it.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from .logging_utils import get_logger

__all__ = ["FLAG_PRIORITY", "worst_flag", "growth_metrics",
           "flag_effectiveness", "contamination_summary",
           "survey_flag_effectiveness"]

log = get_logger("diagnostics")

#: Flag precedence, most severe first -- the same order the photometry applies.
FLAG_PRIORITY = ("a", "b", "c", "d", "0")
_RANK = {f: i for i, f in enumerate(FLAG_PRIORITY)}


def worst_flag(flags: Sequence[str]) -> str:
    """Highest-priority flag in a group (``'a'`` beats ``'b'`` beats ...)."""
    best = "0"
    for f in flags:
        f = str(f)
        if _RANK.get(f, 99) < _RANK.get(best, 99):
            best = f
    return best


def growth_metrics(
    phot: pd.DataFrame,
    ap_kind: str = "km",
    min_steps: int = 4,
    exclude_badphot: bool = True,
) -> pd.DataFrame:
    """
    Per-exposure growth-curve statistics and a contamination indicator.

    Parameters
    ----------
    phot : pandas.DataFrame
        Photometry table from the pipeline.
    ap_kind : str
        Which aperture family defines the curve.  ``"km"`` gives the densest
        radial sampling (19 radii) and is the right choice here.
    min_steps : int
        Exposures with fewer usable annuli than this are dropped -- a curve of
        two points cannot show a rise.
    exclude_badphot : bool
        Drop apertures containing bad pixels, whose shrunken effective area
        would otherwise mimic a surface-brightness step.

    Returns
    -------
    pandas.DataFrame
        One row per exposure:

        ``sb_rise_sigma``
            Significance of the largest outward *rise* in annulus surface
            brightness.  Negative or small positive values are consistent with a
            monotonically declining coma; large positive values indicate that
            something was added to the aperture from outside.
        ``sb_rise_r_pix``
            Radius at which that rise occurs -- for a star, roughly its
            separation from the comet.
        ``growth_ratio``
            :math:`F_{\\rm max}/F_{\\rm min}` over the curve, a crude measure of
            how much flux the outer apertures add.
        ``flag_worst``, ``flag_outer``
            Highest-priority flag over all apertures, and the flag at the
            largest radius.
        ``badphot_any``, ``n_steps``, ``snr_max``, ``wl``, ``epoch``,
        ``gmag_brightest_outer``, ``dist_gmag_nearest``
            Context columns.
    """
    need = {"filename", "ap_kind", "r_ap_pix", "source_sum_mjy",
            "source_sum_err_mjy", "aperture_area_eff_pix2", "sourceflag"}
    missing = need - set(phot.columns)
    if missing:
        raise KeyError(f"growth_metrics needs column(s) {sorted(missing)}")

    all_ap = phot[phot["ap_kind"] == ap_kind]
    sub = all_ap[~all_ap["badphot"].astype(bool)] if exclude_badphot else all_ap
    # `badphot_any` has to come from the *unfiltered* rows: asking it of `sub`
    # after badphot rows were removed would answer False by construction.
    badphot_by_file = all_ap.groupby("filename")["badphot"].any()

    rows = []
    for fname, grp in sub.groupby("filename", sort=False):
        g = grp.sort_values("r_ap_pix")
        r = g["r_ap_pix"].to_numpy(float)
        F = g["source_sum_mjy"].to_numpy(float)
        E = g["source_sum_err_mjy"].to_numpy(float)
        A = g["aperture_area_eff_pix2"].to_numpy(float)

        ok = np.isfinite(F) & np.isfinite(A) & np.isfinite(r)
        r, F, E, A = r[ok], F[ok], E[ok], A[ok]
        if r.size < min_steps + 1:
            continue

        dA = np.diff(A)
        good = dA > 0
        if good.sum() < min_steps:
            continue
        dF = np.diff(F)[good]
        dA = dA[good]
        # Consecutive enclosed fluxes share their inner pixels, so their errors
        # are correlated; adding in quadrature over-states the annulus error and
        # therefore makes `sb_rise_sigma` conservative.  Preferring a
        # conservative indicator is the right way round here: it under-reports
        # contamination rather than inventing it.
        dE = np.sqrt(E[1:] ** 2 + E[:-1] ** 2)[good]
        r_mid = (0.5 * (r[1:] + r[:-1]))[good]

        sb = dF / dA
        sb_err = dE / dA

        if sb.size < 2:
            continue
        rise = np.diff(sb)
        rise_err = np.sqrt(sb_err[1:] ** 2 + sb_err[:-1] ** 2)
        with np.errstate(invalid="ignore", divide="ignore"):
            sig = np.where(rise_err > 0, rise / rise_err, np.nan)
        if not np.isfinite(sig).any():
            continue
        k = int(np.nanargmax(sig))

        flags = g["sourceflag"].astype(str).tolist()
        outer = g.iloc[-1]
        rows.append({
            "filename": fname,
            "epoch": g["epoch"].iloc[0] if "epoch" in g else np.nan,
            "wl": g["wl"].iloc[0] if "wl" in g else np.nan,
            "n_steps": int(sb.size),
            "snr_max": float(np.nanmax(g["snr"].to_numpy(float))) if "snr" in g else np.nan,
            "sb_rise_sigma": float(sig[k]),
            "sb_rise_r_pix": float(r_mid[1:][k]) if k < r_mid.size - 1 else float(r_mid[-1]),
            "growth_ratio": float(np.nanmax(F) / np.nanmin(F)) if np.nanmin(F) > 0 else np.nan,
            "flag_worst": worst_flag(flags),
            "flag_outer": str(outer["sourceflag"]),
            "badphot_any": bool(badphot_by_file.get(fname, False)),
            "gmag_brightest_outer": float(outer.get("gmag_brightest", np.nan)),
            "dist_gmag_nearest": float(outer.get("dist_gmag_nearest", np.nan)),
            "n_gaia_outer": float(outer.get("n_gaia", np.nan)),
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        log.info("growth metrics for %d exposures; median sb_rise_sigma = %.2f",
                 len(out), out["sb_rise_sigma"].median())
    return out


def flag_effectiveness(
    metrics: pd.DataFrame,
    threshold: float = 3.0,
    flag_col: str = "flag_worst",
) -> pd.DataFrame:
    """
    Score each flag against the growth-curve contamination indicator.

    Parameters
    ----------
    metrics : pandas.DataFrame
        Output of :func:`growth_metrics`.
    threshold : float
        ``sb_rise_sigma`` above which an exposure is treated as contaminated.
        3 sigma is the conventional choice and is what the notebook uses; the
        conclusion should not be sensitive to it, which is worth checking.
    flag_col : str
        ``"flag_worst"`` or ``"flag_outer"``.

    Returns
    -------
    pandas.DataFrame
        One row per flag and per cumulative cut (``a``, ``a+b``, ``a+b+c``),
        with:

        ``n_flagged``      exposures the cut removes
        ``recall``         fraction of contaminated exposures it removes
        ``precision``      fraction of what it removes that is contaminated
        ``purity_kept``    fraction of the surviving sample that is clean
        ``n_kept``         exposures surviving the cut

        ``precision`` above the contaminated base rate means the flag carries
        real information; equal to it means the flag is uncorrelated with
        contamination.
    """
    if metrics.empty:
        return pd.DataFrame()

    bad = metrics["sb_rise_sigma"] > float(threshold)
    n_bad, n_tot = int(bad.sum()), len(metrics)
    base = n_bad / n_tot if n_tot else np.nan
    flags = metrics[flag_col].astype(str)

    def _score(name: str, sel: pd.Series) -> dict:
        n_sel = int(sel.sum())
        kept = ~sel
        n_kept = int(kept.sum())
        return {
            "cut": name,
            "n_flagged": n_sel,
            "frac_flagged": n_sel / n_tot if n_tot else np.nan,
            "recall": float((sel & bad).sum() / n_bad) if n_bad else np.nan,
            "precision": float((sel & bad).sum() / n_sel) if n_sel else np.nan,
            "purity_kept": float((kept & ~bad).sum() / n_kept) if n_kept else np.nan,
            "n_kept": n_kept,
        }

    rows = [_score(f"flag == '{f}'", flags == f) for f in FLAG_PRIORITY]
    rows += [
        _score("cut a", flags == "a"),
        _score("cut a+b", flags.isin(["a", "b"])),
        _score("cut a+b+c", flags.isin(["a", "b", "c"])),
        _score("cut badphot", metrics["badphot_any"]),
    ]
    out = pd.DataFrame(rows)
    out.attrs["base_rate"] = base
    out.attrs["n_contaminated"] = n_bad
    out.attrs["n_total"] = n_tot
    out.attrs["threshold"] = float(threshold)
    return out


def contamination_summary(metrics: pd.DataFrame, threshold: float = 3.0) -> str:
    """One-paragraph plain-language reading of :func:`flag_effectiveness`."""
    if metrics.empty:
        return "no exposures with a usable growth curve"
    eff = flag_effectiveness(metrics, threshold=threshold)
    base = eff.attrs["base_rate"]
    lines = [
        f"{eff.attrs['n_total']} exposures with a usable growth curve; "
        f"{eff.attrs['n_contaminated']} ({100*base:.1f} %) show an outward "
        f"surface-brightness rise above {threshold:g} sigma.",
        "",
        f"{'cut':<14s} {'n':>6s} {'recall':>8s} {'precis.':>8s} {'lift':>6s} {'kept pure':>10s}",
    ]
    for _, r in eff.iterrows():
        if not str(r["cut"]).startswith("cut"):
            continue
        lift = r["precision"] / base if base and np.isfinite(r["precision"]) else np.nan
        lines.append(f"{r['cut']:<14s} {r['n_flagged']:>6d} {r['recall']:>8.2f} "
                     f"{r['precision']:>8.2f} {lift:>6.2f} {r['purity_kept']:>10.3f}")
    lines += ["",
              "lift = precision / base rate.  > 1 means the cut preferentially "
              "removes contaminated exposures; ~1 means it is uncorrelated."]
    return "\n".join(lines)


def survey_flag_effectiveness(
    apphot_dir,
    threshold: float = 3.0,
    min_exposures: int = 30,
    cuts: Sequence[str] = ("cut a", "cut a+b", "cut a+b+c"),
) -> pd.DataFrame:
    """
    Score the source flags across every processed target.

    A single comet answers "did the flags work *here*"; the survey answers
    whether the flag definitions generalise -- which is the question that
    matters before they are used to cut a science sample.

    Parameters
    ----------
    apphot_dir : Path
        Directory of ``<slug>.csv`` photometry tables from the pipeline.
    threshold : float
        ``sb_rise_sigma`` above which an exposure counts as contaminated.
    min_exposures : int
        Skip targets with fewer usable growth curves than this; the
        precision/recall of a handful of exposures is noise.
    cuts : sequence of str
        Which cumulative cuts from :func:`flag_effectiveness` to tabulate.

    Returns
    -------
    pandas.DataFrame
        One row per target: ``base_rate`` (contaminated fraction) plus
        ``recall_``/``lift_``/``kept_`` for each cut.  ``lift`` is precision
        divided by the base rate, so 1 means the cut is uncorrelated with
        contamination and anything above it is informative.
    """
    from pathlib import Path as _Path

    rows = []
    for f in sorted(_Path(apphot_dir).glob("*.csv")):
        d = pd.read_csv(f, dtype={"sourceflag": str})
        if "sourceflag" not in d.columns:
            continue
        m = growth_metrics(d, ap_kind="km")
        if m.empty or len(m) < min_exposures:
            log.info("%s: only %d usable growth curve(s); skipped", f.stem, len(m))
            continue
        eff = flag_effectiveness(m, threshold=threshold)
        base = eff.attrs["base_rate"]
        row = {"objdesig": str(d["objdesig"].iloc[0]), "n_exposures": len(m),
               "base_rate": base, "n_gaia_max": int(d["n_gaia"].max())}
        for cut in cuts:
            hit = eff[eff["cut"] == cut]
            if hit.empty:
                continue
            e = hit.iloc[0]
            tag = cut.split()[-1]
            row[f"recall_{tag}"] = e["recall"]
            row[f"lift_{tag}"] = e["precision"] / base if base else np.nan
            row[f"kept_{tag}"] = 1.0 - e["frac_flagged"]
        rows.append(row)

    out = pd.DataFrame(rows)
    if not out.empty:
        log.info("survey: %d targets, median contaminated fraction %.3f",
                 len(out), out["base_rate"].median())
    return out
