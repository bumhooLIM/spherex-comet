"""
Reading the revised photometry and writing every intermediate and final product.

Input schema
------------
``data/apphot_revised/<target>.csv`` from ``spherex_apphot`` -- one row per
(exposure, aperture).  The columns this package relies on:

``filename, epoch, wl, wlwidth, ap_label, ap_kind, r_ap_km, r_ap_pix,
source_sum_mjy, source_sum_err_mjy, source_sum_err_empirical_mjy,
flux_distcorr_mjy, flux_distcorr_err_mjy, distcorr_factor,
sourceflag, badphot, frac_badpix_ap, r_hel, r_obs, jd_utc, detector``

Two facts about that table drive the design here.  Rows are duplicated once per
aperture, so any per-exposure quantity must be read at one aperture or after
``drop_duplicates("filename")``.  And the aperture set is *not* the same for
every exposure: the photometry refuses an aperture below the PSF or beyond the
sky annulus, so :func:`aperture_for` has to check coverage rather than assume a
radius exists.

The photometry tables are never written to.  The group assignment lives in
``data/comspec/phase_assignment.csv`` and is merged on read.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import re

import numpy as np
import pandas as pd

from . import directory as _dir
from .config import (ALL_EM_WINDOWS, EMISSION_DTYPES, EMISSION_WINDOWS, ApertureConfig,
                     FitConfig, FlagPolicy, Variant)
from .logging_utils import get_logger

__all__ = [
    "slug", "aperture_label", "list_targets", "load_apphot", "exposure_table",
    "aperture_for", "aperture_choice", "select_spectrum", "PhaseAssignment", "heliocentric_velocity",
    "emission_paths", "save_emission", "load_summary", "load_points", "list_catalog",
    "load_fit_input", "save_fit_table", "save_fit_lines", "REQUIRED_COLUMNS",
]

log = get_logger("dataio")

#: columns read from the photometry table
USECOLS = [
    "filename", "obsid", "date_obs", "epoch", "wl", "wlwidth", "ap_label", "ap_kind", "r_ap_km", "r_ap_pix",
    "source_sum_mjy", "source_sum_err_mjy", "source_sum_err_empirical_mjy",
    "flux_distcorr_mjy", "flux_distcorr_err_mjy", "distcorr_factor",
    "sourceflag", "badphot", "frac_badpix_ap", "n_gaia", "gmag_eff",
    "r_hel", "r_obs", "jd_utc", "detector", "snr",
    "psf_fwhm_pix", "r_in_pix", "pixel_scale_km",
]

REQUIRED_COLUMNS = ("target", "r_ap_km", "phase", "band", "role", "wl", "wlwidth",
                    "emis_mjy", "emis_err_mjy", "emis_raw_mjy", "emis_raw_err_mjy",
                    "distcorr_factor", "r_hel", "r_obs")

_CACHE: Dict[str, pd.DataFrame] = {}
_CACHE_MAX = 3


def slug(target: str) -> str:
    """``"2025 W2" -> "2025W2"``; ``"73P/Schwassmann" -> "73P"``."""
    return str(target).split("/")[0].replace(" ", "").strip()


def aperture_label(r_km) -> str:
    """``20000 -> "2e4"``, ``60000 -> "6e4"``, ``100000 -> "1e5"`` -- the file-name convention."""
    e = int(np.floor(np.log10(float(r_km))))
    return f"{float(r_km) / 10 ** e:g}e{e}"


def list_targets(apphot_dir: Optional[Path] = None) -> List[str]:
    d = Path(apphot_dir) if apphot_dir else _dir.APPHOT_DIR
    return sorted(p.stem for p in d.glob("*.csv"))


def load_apphot(target: str, apphot_dir: Optional[Path] = None) -> pd.DataFrame:
    """Read one photometry table (needed columns only, small LRU cache)."""
    key = slug(target)
    if key in _CACHE:
        return _CACHE[key]
    d = Path(apphot_dir) if apphot_dir else _dir.APPHOT_DIR
    path = d / f"{key}.csv"
    if not path.is_file():
        raise FileNotFoundError(f"{path} (available: {list_targets(d)[:8]} ...)")
    cols = pd.read_csv(path, nrows=0).columns
    use = [c for c in USECOLS if c in cols]
    df = pd.read_csv(path, usecols=use, dtype={"sourceflag": str, "obsid": str, "date_obs": str})
    # The photometry CSVs were written with float_format="%.8g", which rounds a Julian date of
    # ~2.46e6 to 0.1 day.  date_obs is a string and keeps full precision, so the exposure time
    # is rebuilt from it; jd_utc is retained only where date_obs is missing.
    if "date_obs" in df and df.date_obs.notna().any():
        from astropy.time import Time
        t = pd.to_datetime(df.date_obs, errors="coerce", utc=True)
        jd = Time(t.dt.tz_convert(None).to_numpy(), format="datetime64", scale="utc").jd
        df["jd_utc_exact"] = np.where(np.isfinite(jd), jd, df.jd_utc)
    else:
        df["jd_utc_exact"] = df.jd_utc
    # heliocentric radial velocity of every pointing from the ephemeris r_hel(t) -- the CO
    # g-factor depends on it (Swings effect, see fluorescence.swings_factor)
    df["v_hel_kms"] = _velocity_per_pointing(df)
    if "source_sum_err_empirical_mjy" in df and "flux_distcorr_err_empirical_mjy" not in df:
        # the empirical error in distance-corrected space is the same factor away
        df["flux_distcorr_err_empirical_mjy"] = (
            df["source_sum_err_empirical_mjy"] * df["distcorr_factor"])
    while len(_CACHE) >= _CACHE_MAX:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = df
    return df


# ------------------------------------------------------------------- heliocentric velocity
AU_KM_PER_DAY_TO_KMS = 1.495978707e8 / 86400.0


def heliocentric_velocity(jd, r_hel, sigma_days: float = 5.0, degree: int = 4,
                          min_effective_points: float = 6.0) -> np.ndarray:
    """
    ``d r_h / dt`` [km s^-1] at every exposure from the ephemeris ``r_hel(t)`` of one target.

    The photometry carries ``r_hel`` to 1e-6 au and the exact exposure time from ``date_obs``, so
    the radial velocity follows from a local polynomial fit of ``r_hel`` against time.  At each
    distinct time a polynomial of ``degree`` is fitted to *all* the target's times with Gaussian
    weights ``exp(-(t - t0)^2 / 2 sigma^2)``; ``sigma`` is widened (x3, x10, then flat) until the
    effective number of points reaches ``min_effective_points``, and the degree is lowered when
    even that is small.  Against JPL Horizons heliocentric rates this recovers v_h to 0.002 km/s
    (95th percentile) and 0.05 km/s at worst over 118 pointings of 10 comets, edge pointings of
    sparse epochs near perihelion included -- the curvature of r_h(t) there defeats a wide
    quadratic, which the quartic follows.  Positive = receding from the Sun.  A target observed at
    a single instant gets NaN, which the model treats as a Swings factor of 1 (and the fit flags).
    """
    jd = np.asarray(jd, dtype=float)
    r = np.asarray(r_hel, dtype=float)
    out = np.full(len(jd), np.nan)
    ok = np.isfinite(jd) & np.isfinite(r)
    if ok.sum() < 2:
        return out
    t_u, idx = np.unique(jd[ok], return_inverse=True)
    if len(t_u) < 2:
        return out
    r_u = np.zeros(len(t_u))
    r_u[idx] = r[ok]                                       # one r_hel per distinct time
    v_u = np.full(len(t_u), np.nan)
    for i, t0 in enumerate(t_u):
        x = t_u - t0
        for s in (sigma_days, 3 * sigma_days, 10 * sigma_days, np.inf):
            w = np.ones_like(x) if not np.isfinite(s) else np.exp(-0.5 * (x / s) ** 2)
            n_eff = w.sum() ** 2 / (w ** 2).sum()
            if n_eff >= min_effective_points:
                break
        deg = int(min(degree, max(1, round(n_eff) - 2)))
        m = w > 1e-6
        if m.sum() < 2 or np.ptp(x[m]) <= 0:
            continue
        c = np.polyfit(x[m], r_u[m], deg, w=np.sqrt(w[m]))
        v_u[i] = c[-2] * AU_KM_PER_DAY_TO_KMS               # d r / d t at t0 [au/day -> km/s]
    out[ok] = v_u[idx]
    return out


def _velocity_per_pointing(df: pd.DataFrame) -> np.ndarray:
    """Velocity per pointing (``obsid`` shares one ephemeris) mapped back onto every row."""
    if "obsid" in df:
        key = df["obsid"].astype(str)
    else:
        key = df["filename"].astype(str)
    ptg = (pd.DataFrame({"key": key, "jd": df["jd_utc_exact"], "r": df["r_hel"]})
             .groupby("key", sort=False).agg(jd=("jd", "min"), r=("r", "first")))
    ptg["v"] = heliocentric_velocity(ptg["jd"].to_numpy(), ptg["r"].to_numpy())
    return key.map(ptg["v"]).to_numpy(dtype=float)


# ---------------------------------------------------------------------------- exposures
def exposure_table(target: str) -> pd.DataFrame:
    """
    One row per exposure, time-ordered, with the band-sampling flags the grouping needs.

    An exposure "samples" a band when at least one of its aperture rows inside the
    window is not ``badphot`` and has a finite flux -- the same criterion the
    original notebook used, so the regrouping is comparable.
    """
    df = load_apphot(target)
    obs = df.drop_duplicates("filename")
    # The atomic unit is one pointing (obsid): every detector of a pointing shares its ephemeris,
    # so r_hel is identical within it and a cut can never fall inside one.  This is the unit the
    # original notebook used (its jd_utc was per pointing).  The exposure time is rebuilt from
    # date_obs because the stored jd_utc is rounded to 0.1 day.
    ep = (obs.groupby("obsid")
             .agg(jd_utc=("jd_utc_exact", "min"), r_hel=("r_hel", "first"),
                  r_hel_lo=("r_hel", "min"), r_hel_hi=("r_hel", "max"),
                  r_obs=("r_obs", "first"), epoch=("epoch", "first"),
                  n_epoch=("epoch", "nunique"), weight=("filename", "size"),
                  filename=("filename", "first"))
             .reset_index().sort_values("jd_utc").reset_index(drop=True))
    assert (ep.n_epoch == 1).all(), f"{target}: one exposure spans two epochs"
    assert np.allclose(ep.r_hel_lo, ep.r_hel_hi), f"{target}: r_hel not unique within an exposure"
    ep = ep.drop(columns=["n_epoch", "r_hel_lo", "r_hel_hi"])

    good = df[~df.badphot.astype(bool) & np.isfinite(df.source_sum_mjy)]
    for band, (lo, hi) in EMISSION_WINDOWS.items():
        hit = set(good.loc[good.wl.between(lo, hi), "obsid"])
        ep[band] = ep.obsid.isin(hit)
    return ep


# ----------------------------------------------------------------------------- aperture
def aperture_for(target: str, cfg: ApertureConfig) -> Tuple[float, str, float]:
    """The single aperture a target is analysed at: ``(r_ap_km, ap_label, coverage)``."""
    c = aperture_choice(target, cfg)
    return float(c["r_ap_km"]), str(c["ap_label"]), float(c["coverage"])


def _coverage_table(df: pd.DataFrame, cfg: ApertureConfig) -> pd.DataFrame:
    n_exp = df.filename.nunique()
    km = df[df.ap_kind == cfg.kind]
    cov = (km.groupby(["ap_label", "r_ap_km"])
             .agg(coverage=("filename", "nunique"), r_ap_pix=("r_ap_pix", "median"),
                  psf_fwhm_pix=("psf_fwhm_pix", "median") if "psf_fwhm_pix" in km else ("r_ap_pix", "size"))
             .reset_index())
    cov["coverage"] = cov["coverage"] / n_exp
    return cov.sort_values("r_ap_km").reset_index(drop=True)


def _rule_rh(target: str, df: pd.DataFrame, cov: pd.DataFrame, cfg: ApertureConfig) -> dict:
    """The previous rule: ``near_km`` inside ``rh_split_au``, ``far_km`` beyond, promoted for coverage."""
    rh_mean = float(df.drop_duplicates("filename").r_hel.mean())
    want = cfg.near_km if rh_mean < cfg.rh_split_au else cfg.far_km
    hit = cov[np.isclose(cov.r_ap_km, want)]
    if len(hit) and float(hit.coverage.iloc[0]) >= cfg.min_coverage:
        return dict(r_ap_km=float(want), ap_label=str(hit.ap_label.iloc[0]),
                    coverage=float(hit.coverage.iloc[0]), rule="rh", promoted=False)
    bigger = cov[(cov.r_ap_km >= want) & (cov.coverage >= cfg.min_coverage)]
    if bigger.empty:
        # last resort: the best-covered aperture at or above the rule radius
        bigger = cov[cov.r_ap_km >= want].sort_values("coverage", ascending=False)
        if bigger.empty:
            raise ValueError(f"{target}: no {cfg.kind} aperture >= {want:g} km in the table")
    row = bigger.iloc[0]
    have = float(hit.coverage.iloc[0]) if len(hit) else 0.0
    log.warning("%s: rule aperture %g km covers %.0f %% of exposures (< %.0f %%); "
                "promoted to %g km (%.0f %%)", target, want, 100 * have,
                100 * cfg.min_coverage, row.r_ap_km, 100 * row.coverage)
    return dict(r_ap_km=float(row.r_ap_km), ap_label=str(row.ap_label),
                coverage=float(row.coverage), rule="rh", promoted=True)


def aperture_choice(target: str, cfg: ApertureConfig) -> dict:
    """
    The single aperture a target is analysed at, with the evidence for the choice.

    ``cfg.rule == "snr"`` (default since 2026-09-12): every ``km`` aperture present for at
    least ``cfg.min_coverage`` of the exposures, at least ``cfg.min_psf_mult`` PSF FWHM in
    radius, and no larger than ``cfg.max_ap_frac_annulus`` times the sky annulus' inner
    radius (so the background is measured outside the coma the aperture integrates) is a
    candidate; its score is the median S/N of the *star-free* channels inside the emission
    windows (finite flux, ``frac_badpix_ap <= 0.05``, source flag ``0``), and it needs at
    least ``cfg.min_clean_channels`` of them.  The chosen aperture is the smallest candidate
    within ``cfg.snr_tol`` of the best score -- the S/N of a coma is nearly flat with radius,
    and the smaller aperture keeps the annulus farther away in units of its own radius.
    When no aperture satisfies the PSF or annulus bound the bound is relaxed in that order
    (logged, ``relaxed`` in the result); in a field so crowded that no aperture can be
    scored the smallest bounded aperture is taken (``relaxed = "contaminated"``); the
    ``"rh"`` rule is the last resort.  ``cfg.rule == "rh"``: the previous rule.

    Returns
    -------
    dict
        ``r_ap_km, ap_label, coverage, rule`` and, for ``"snr"``, ``snr_median, snr_best,
        r_ap_km_best, r_in_km, n_candidates, relaxed``.
    """
    df = load_apphot(target)
    cov = _coverage_table(df, cfg)
    if cfg.rule == "rh":
        return _rule_rh(target, df, cov, cfg)
    if cfg.rule != "snr":
        raise ValueError(f"unknown aperture rule {cfg.rule!r} (snr | rh)")

    # annulus inner radius at the comet [km], the same for every aperture of an exposure
    r_in_km = float(np.nanmedian(df.r_in_pix * df.pixel_scale_km)) \
        if "r_in_pix" in df and "pixel_scale_km" in df else np.nan
    km = df[df.ap_kind == cfg.kind]
    em = np.zeros(len(km), bool)
    for lo, hi in EMISSION_WINDOWS.values():
        em |= km.wl.between(lo, hi).to_numpy()
    ecol = "source_sum_err_empirical_mjy" if "source_sum_err_empirical_mjy" in km else "source_sum_err_mjy"
    # The score uses only channels without any Gaia source inside r_ap + PSF FWHM (flag "0"):
    # star flux inflates the "S/N" of a large aperture on a faint comet (2024 N1 at 80 000 km:
    # 46 % of its emission-window channels carry flag b and two of them are 90 mJy spikes), so a
    # score over flagged channels would drive the rule to the most contaminated aperture.
    clean = (em & np.isfinite(km.source_sum_mjy) & np.isfinite(km[ecol]) & (km[ecol] > 0)
             & ~(km.badphot.astype(bool) & (km.frac_badpix_ap > 0.05))
             & (km.sourceflag.astype(str) == "0"))
    d = km[clean]
    score = (d.source_sum_mjy / d[ecol]).groupby(d.r_ap_km).agg(["median", "size"])
    cov["snr_median"] = cov.r_ap_km.map(score["median"])
    cov["n_snr"] = cov.r_ap_km.map(score["size"]).fillna(0).astype(int)
    flagged = km[em & (km.sourceflag.astype(str).isin(["a", "b"]))].groupby("r_ap_km").size()
    n_em = km[em].groupby("r_ap_km").size()
    cov["frac_flag_ab"] = cov.r_ap_km.map((flagged / n_em).fillna(0.0)).fillna(0.0)

    geom = (cov.coverage >= cfg.min_coverage)
    psf_ok = cov.r_ap_pix >= cfg.min_psf_mult * cov.psf_fwhm_pix.fillna(0)
    ann_ok = ~np.isfinite(r_in_km) | (cov.r_ap_km <= cfg.max_ap_frac_annulus * r_in_km)
    scored = (cov.n_snr >= cfg.min_clean_channels) & np.isfinite(cov.snr_median)
    relaxed = ""
    cand = cov[geom & psf_ok & ann_ok & scored]
    if cand.empty:
        cand, relaxed = cov[geom & (cov.r_ap_pix >= 1.0 * cov.psf_fwhm_pix.fillna(0)) & ann_ok & scored], "psf"
    if cand.empty:
        cand, relaxed = cov[geom & psf_ok & scored], "annulus"
    if cand.empty:
        # a crowded field: no aperture has enough star-free channels to be scored.  Take the
        # smallest aperture the bounds allow -- the least contaminated -- rather than any S/N.
        pool = cov[geom & psf_ok & ann_ok]
        if pool.empty:
            pool = cov[geom & (cov.r_ap_pix >= 1.0 * cov.psf_fwhm_pix.fillna(0))]
        if pool.empty:
            out = _rule_rh(target, df, cov, cfg)
            out.update(rule="snr->rh", relaxed="all", r_in_km=r_in_km, n_candidates=0,
                       snr_median=np.nan, snr_best=np.nan, r_ap_km_best=np.nan, frac_flag_ab=np.nan)
            log.warning("%s: no aperture satisfies the S/N rule; fell back to the r_h rule", target)
            return out
        row = pool.sort_values("r_ap_km").iloc[0]
        log.warning("%s: fewer than %d star-free emission channels at every aperture; smallest "
                    "bounded aperture %g km chosen", target, cfg.min_clean_channels, row.r_ap_km)
        return dict(r_ap_km=float(row.r_ap_km), ap_label=str(row.ap_label), coverage=float(row.coverage),
                    rule="snr", snr_median=float(row.snr_median), snr_best=np.nan, r_ap_km_best=np.nan,
                    r_in_km=r_in_km, n_candidates=int(len(pool)), relaxed="contaminated",
                    frac_flag_ab=float(row.frac_flag_ab))
    best = float(cand.snr_median.max())
    # "within snr_tol of the best" must also work when every score is negative (faint targets)
    ok = cand[cand.snr_median >= best - cfg.snr_tol * abs(best)].sort_values("r_ap_km")
    if ok.empty:
        ok = cand.sort_values("snr_median", ascending=False)
    row = ok.iloc[0]
    if relaxed:
        log.warning("%s: aperture rule relaxed (%s): %g km chosen", target, relaxed, row.r_ap_km)
    return dict(r_ap_km=float(row.r_ap_km), ap_label=str(row.ap_label), coverage=float(row.coverage),
                rule="snr", snr_median=float(row.snr_median), snr_best=best,
                r_ap_km_best=float(cand.loc[cand.snr_median.idxmax(), "r_ap_km"]),
                r_in_km=r_in_km, n_candidates=int(len(cand)), relaxed=relaxed,
                frac_flag_ab=float(row.frac_flag_ab))


# ------------------------------------------------------------------------- flag policy
def select_spectrum(df: pd.DataFrame, ap_label: str, variant: Variant,
                    assignment: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    Rows of one aperture that pass the variant's flag policy, as a spectrum table.

    Parameters
    ----------
    df : DataFrame
        A photometry table from :func:`load_apphot`.
    ap_label : str
        Aperture to keep, e.g. ``"km20000"``.
    variant : Variant
        Supplies the flag policy and the flux/error columns.
    assignment : DataFrame, optional
        ``(filename, phase)`` from the regrouping; merged on when given.

    Returns
    -------
    DataFrame
        Sorted by wavelength, with ``flux``/``err`` in the variant's flux space,
        ``flux_raw``/``err_raw`` in physical mJy, ``distcorr_factor``, the flag
        columns, and a ``rejection`` dict in ``.attrs``.
    """
    pol = variant.flags
    d = df[df.ap_label == ap_label]
    n0 = len(d)

    ok_fin = (np.isfinite(d[variant.flux_column]) & np.isfinite(d[variant.err_column])
              & (d[variant.err_column] > 0) & np.isfinite(d.distcorr_factor)
              & (d.distcorr_factor > 0))
    if pol.drop_badphot:
        if pol.max_frac_badpix is None:
            ok_bp = ~d.badphot.astype(bool)
        else:
            ok_bp = ~(d.badphot.astype(bool) & (d.frac_badpix_ap > pol.max_frac_badpix))
    else:
        ok_bp = pd.Series(True, index=d.index)
    ok_flag = ~d.sourceflag.astype(str).isin(pol.drop_flags) if pol.drop_flags \
        else pd.Series(True, index=d.index)

    keep = ok_fin & ok_bp & ok_flag
    out = d[keep].copy()
    out["flux"] = out[variant.flux_column]
    out["err"] = out[variant.err_column]
    out["flux_raw"] = out["source_sum_mjy"]
    out["err_raw"] = out[variant.error_column]
    # the stored jd_utc is rounded to 0.1 day (see load_apphot); every product carries the
    # exact time rebuilt from date_obs instead
    out["jd_utc"] = out["jd_utc_exact"]
    if assignment is not None:
        out = out.merge(assignment[["filename", "phase"]], on="filename", how="left")
        if out.phase.isna().any():
            raise ValueError(f"{int(out.phase.isna().sum())} rows have no phase assignment")
        out["phase"] = out.phase.astype(int)
    cols = ["filename", "epoch", "phase", "wl", "wlwidth", "flux", "err", "flux_raw", "err_raw",
            "distcorr_factor", "snr", "r_hel", "r_obs", "v_hel_kms", "jd_utc", "detector",
            "sourceflag", "badphot", "frac_badpix_ap", "n_gaia", "gmag_eff"]
    out = out[[c for c in cols if c in out.columns]].sort_values("wl").reset_index(drop=True)
    out.attrs["rejection"] = dict(
        n_total=n0, n_nonfinite=int((~ok_fin).sum()),
        n_badphot=int((ok_fin & ~ok_bp).sum()),
        n_flag=int((ok_fin & ok_bp & ~ok_flag).sum()), n_kept=int(keep.sum()))
    return out


# ---------------------------------------------------------------------- phase assignment
class PhaseAssignment:
    """The per-exposure ``phase`` labels, persisted as ``data/comspec/phase_assignment.csv``."""

    PATH = "phase_assignment.csv"
    MAP = "phase_map.csv"

    @classmethod
    def path(cls) -> Path:
        return _dir.DATA_DIR / cls.PATH

    @classmethod
    def save(cls, assignment: pd.DataFrame, group_map: pd.DataFrame) -> Tuple[Path, Path]:
        _dir.ensure_dirs()
        p1 = cls.path()
        p2 = _dir.RESULT_DIR / cls.MAP
        assignment.to_csv(p1, index=False)
        group_map.to_csv(p2, index=False)
        return p1, p2

    @classmethod
    def load(cls) -> pd.DataFrame:
        p = cls.path()
        if not p.is_file():
            raise FileNotFoundError(f"{p}: run the grouping step first")
        return pd.read_csv(p, dtype={"target": str})

    @classmethod
    def for_target(cls, target: str, assignment: Optional[pd.DataFrame] = None) -> pd.DataFrame:
        a = assignment if assignment is not None else cls.load()
        return a[a.target == slug(target)]


# ---------------------------------------------------------------------------- emission
def emission_paths(variant: str, target: str, r_ap_km) -> Tuple[Path, Path]:
    d = _dir.variant_dirs(variant)["emission"]
    stem = f"{slug(target)}_{aperture_label(r_ap_km)}km"
    return d / f"{stem}.csv", d / f"{stem}_points.csv"


def save_emission(variant: str, summary: pd.DataFrame, points: pd.DataFrame) -> Tuple[Path, Path]:
    """
    Write (merge) the per-band summary and the point-level spectrum of one target.

    Rows are keyed on ``(phase, band)`` and replaced in place, so re-running one
    phase updates that phase without touching the others.
    """
    if summary.empty:
        return None, None
    target = str(summary.target.iloc[0])
    r_ap = float(summary.r_ap_km.iloc[0])
    p_sum, p_pts = emission_paths(variant, target, r_ap)
    p_sum.parent.mkdir(parents=True, exist_ok=True)

    def _merge(new: pd.DataFrame, path: Path) -> pd.DataFrame:
        key = ["phase", "band"]
        if path.exists():
            old = pd.read_csv(path, dtype=EMISSION_DTYPES)
            common = new[key].drop_duplicates()
            hit = old[key].merge(common.assign(_h=1), on=key, how="left")["_h"].notna().to_numpy()
            keep_old = old[~hit]
            if len(keep_old):
                new = pd.concat([keep_old, new], ignore_index=True)
        return new.sort_values(key + (["wl"] if "wl" in new.columns else [])).reset_index(drop=True)

    _merge(summary, p_sum).to_csv(p_sum, index=False)
    _merge(points, p_pts).to_csv(p_pts, index=False)
    return p_sum, p_pts


def load_summary(variant: str, target: str, r_ap_km) -> pd.DataFrame:
    p, _ = emission_paths(variant, target, r_ap_km)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_csv(p, dtype=EMISSION_DTYPES)


def load_points(variant: str, target: str, r_ap_km) -> pd.DataFrame:
    _, p = emission_paths(variant, target, r_ap_km)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_csv(p, dtype=EMISSION_DTYPES)


def list_catalog(variant: str) -> pd.DataFrame:
    """Every (target, aperture, phase) in a variant's emission directory, with band verdicts."""
    d = _dir.variant_dirs(variant)["emission"]
    rows = []
    for f in sorted(d.glob("*km.csv")):
        if f.stem.endswith("_points"):
            continue
        s = pd.read_csv(f, dtype=EMISSION_DTYPES)
        for (tgt, ap, ph), g in s.groupby(["target", "r_ap_km", "phase"]):
            ok = g[g.verdict.isin(("PASS", "WARN")) & (g.n_emission > 0)]
            rows.append(dict(target=tgt, r_ap_km=float(ap), phase=int(ph),
                             epoch=int(g.epoch.iloc[0]) if "epoch" in g else -1,
                             n_bands_ok=len(ok), bands_ok="+".join(sorted(ok.band)),
                             n_emission_ok=int(ok.n_emission.sum()) if len(ok) else 0))
    out = pd.DataFrame(rows)
    return out.sort_values(["target", "phase"]).reset_index(drop=True) if len(out) else out


def load_fit_input(variant: str, target: str, r_ap_km, phase: int,
                   cfg: Optional[FitConfig] = None) -> pd.DataFrame:
    """
    The channels that will be fitted for one (target, aperture, phase).

    Only ``role == "emission"`` channels from bands whose continuum earned an
    accepted verdict.  A FAIL band's residual is the residual of a rejected
    continuum and measures nothing.
    """
    cfg = cfg or FitConfig()
    summ = load_summary(variant, target, r_ap_km)
    pts = load_points(variant, target, r_ap_km)
    summ = summ[summ.phase == phase]
    pts = pts[(pts.phase == phase) & (pts.role == "emission")]
    ok = summ[summ.verdict.isin(cfg.accept_verdicts) & (summ.n_emission > 0)]
    if cfg.bands is not None:
        ok = ok[ok.band.isin(cfg.bands)]
    if ok.empty:
        raise ValueError(f"{target} phase {phase} @ {r_ap_km:g} km: no usable band "
                         f"({dict(zip(summ.band, summ.verdict))})")
    cont_cols = [c for c in ok.columns if re.match(r"^(cont_lam_ref_um|cont_c\d|cont_cov_\d\d)$", c)]
    out = pts[pts.band.isin(ok.band)].merge(
        ok[["band", "verdict", "err_scale", "chi2_red", *cont_cols]].rename(columns={"chi2_red": "cont_chi2_red"}),
        on="band", how="left")
    missing = [c for c in REQUIRED_COLUMNS if c not in out.columns]
    if missing:
        raise ValueError(f"fit input missing columns: {missing}")
    out = out[np.isfinite(out.emis_raw_mjy) & np.isfinite(out.emis_raw_err_mjy)
              & (out.emis_raw_err_mjy > 0)]
    if out.empty:
        raise ValueError(f"{target} phase {phase}: no finite channels")
    return out.sort_values("wl").reset_index(drop=True)


# ------------------------------------------------------------------------------ results
def save_fit_table(variant: str, rows, filename: str = "gas_fit.csv") -> Path:
    d = _dir.variant_dirs(variant)["results"]
    d.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows) if not isinstance(rows, pd.DataFrame) else rows
    p = d / filename
    df.to_csv(p, index=False)
    return p


def save_fit_lines(variant: str, fit, points: pd.DataFrame, curves: dict) -> Tuple[Path, Path]:
    d = _dir.variant_dirs(variant)["lines"]
    d.mkdir(parents=True, exist_ok=True)
    stem = f"{fit.target}_{aperture_label(fit.r_ap_km)}km_ph{fit.phase}"
    curve = pd.DataFrame({"lam_um": curves["lam"], "model_intrinsic_mjy": curves["intrinsic"],
                          "model_convolved_mjy": curves["convolved"]})
    for s, v in curves["per_species"].items():
        curve[f"intrinsic_{s}_mjy"] = v
    curve.insert(0, "target", fit.target)
    curve.insert(1, "r_ap_km", fit.r_ap_km)
    curve.insert(2, "phase", fit.phase)
    curve["r_hel_repr"] = curves["r_hel_repr"]
    curve["r_obs_repr"] = curves["r_obs_repr"]
    curve["v_hel_repr_kms"] = curves.get("v_hel_repr", np.nan)
    p_curve = d / f"{stem}_curve.csv"
    curve.to_csv(p_curve, index=False)
    obs = points.copy()
    obs["model_mjy"] = curves["at_points"]
    obs["resid_mjy"] = obs["emis_raw_mjy"] - obs["model_mjy"]
    obs["resid_sigma"] = obs["resid_mjy"] / obs["emis_raw_err_mjy"]
    keep = ["target", "r_ap_km", "phase", "band", "wl", "wlwidth", "emis_raw_mjy",
            "emis_raw_err_mjy", "emis_mjy", "emis_err_mjy", "distcorr_factor",
            "model_mjy", "resid_mjy", "resid_sigma", "r_hel", "r_obs", "v_hel_kms", "sourceflag"]
    p_pts = d / f"{stem}_points.csv"
    obs[[c for c in keep if c in obs.columns]].to_csv(p_pts, index=False)
    return p_curve, p_pts
