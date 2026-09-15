"""
Reflectance spectra (review item S11).

The primitive notebook computed ``source_sum_mjy / sun_jy`` and normalised by a
median.  Three defects followed from that, and all three are fixed here.

1. **No distance scaling.**  A comet's reflected flux scales as
   ``1 / (r_hel^2 * r_obs^2)``: the illuminating sunlight falls off with
   heliocentric distance and the reflected light again with observer distance.
   Without that factor, spectra from epochs at different distances are not on a
   common scale, and combining them mixes a geometric trend into what is read as
   a compositional one.
2. **Fragile error propagation.**  ``refl_norm * (flux_err / flux)`` diverges as
   the flux approaches zero and changes sign with it, so the faintest and most
   interesting points got the least trustworthy error bars.  The relation is
   linear in flux, so the uncertainty simply scales -- and the uncertainty of
   the normalisation itself is now included.
3. **Non-robust normalisation.**  A plain median over every wavelength included
   the water band it was meant to be independent of.  Normalisation now uses an
   inverse-variance weighted, sigma-clipped mean over a configurable continuum
   window.

Definition
----------
Up to a constant that cancels in the normalisation::

    R(lambda) = F_obs(lambda) / F_sun(lambda) * r_hel^2 * r_obs^2

``sun_jy`` is the solar spectrum convolved to the same pixel bandpass, so the
ratio is taken band by band and the instrument's spectral response divides out.
The absolute normalisation of ``sun_jy`` is irrelevant provided it is consistent
across wavelength, which is why only *normalised* reflectance is reported.
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .logging_utils import get_logger

__all__ = ["add_reflectance", "normalise_reflectance", "bin_spectrum",
           "REFLECTANCE_COLUMNS"]

log = get_logger("reflectance")

REFLECTANCE_COLUMNS = ["refl", "refl_err", "refl_norm", "refl_norm_err"]


def _weighted_clipped_mean(
    values: np.ndarray,
    errors: np.ndarray,
    sigma: float = 3.0,
    maxiters: int = 5,
) -> Tuple[float, float, int]:
    """
    Inverse-variance weighted mean with iterative sigma clipping.

    Parameters
    ----------
    values, errors : ndarray
        Measurements and their 1-sigma uncertainties.
    sigma : float
        Clipping threshold, in units of the *scatter* of the sample (not of the
        individual error bars), so a badly under-estimated error set cannot
        defeat the clipping.
    maxiters : int
        Maximum clipping iterations.

    Returns
    -------
    mean, mean_err, n_used : float, float, int
        ``(nan, nan, 0)`` if nothing usable survives.
    """
    v = np.asarray(values, dtype=np.float64)
    e = np.asarray(errors, dtype=np.float64)
    good = np.isfinite(v)
    # Points with a missing or non-positive error still carry information; they
    # are given the median weight rather than discarded.
    e = np.where(np.isfinite(e) & (e > 0), e, np.nan)
    if good.sum() == 0:
        return float("nan"), float("nan"), 0
    if not np.isfinite(e[good]).any():
        e = np.ones_like(v)
    else:
        e = np.where(np.isfinite(e), e, np.nanmedian(e[good]))

    keep = good.copy()
    for _ in range(int(maxiters)):
        if keep.sum() < 2:
            break
        w = 1.0 / e[keep] ** 2
        mu = float(np.sum(w * v[keep]) / np.sum(w))
        scatter = float(np.std(v[keep], ddof=1))
        if not np.isfinite(scatter) or scatter == 0:
            break
        new = good & (np.abs(v - mu) <= sigma * scatter)
        if new.sum() == keep.sum() or new.sum() < 2:
            keep = new if new.sum() >= 2 else keep
            break
        keep = new

    if keep.sum() == 0:
        return float("nan"), float("nan"), 0
    w = 1.0 / e[keep] ** 2
    mu = float(np.sum(w * v[keep]) / np.sum(w))
    mu_err = float(np.sqrt(1.0 / np.sum(w)))
    return mu, mu_err, int(keep.sum())


def add_reflectance(df: pd.DataFrame, cfg) -> pd.DataFrame:
    """
    Add the un-normalised reflectance and its uncertainty.

    Parameters
    ----------
    df : pandas.DataFrame
        Photometry table; needs ``source_sum_mjy``, ``source_sum_err_mjy``,
        ``sun_jy`` and, when ``cfg.refl_apply_distance``, ``r_hel`` and
        ``r_obs``.
    cfg : Config

    Returns
    -------
    pandas.DataFrame
        A copy with ``refl`` and ``refl_err`` added.

    Notes
    -----
    Negative fluxes propagate to negative reflectances and are kept (S5).
    """
    out = df.copy()
    missing = {"source_sum_mjy", "sun_jy"} - set(out.columns)
    if missing:
        raise KeyError(f"reflectance needs column(s) {sorted(missing)}")

    flux_jy = out["source_sum_mjy"].to_numpy(dtype=float) * 1e-3
    err_jy = out.get("source_sum_err_mjy", pd.Series(np.nan, index=out.index)) \
                .to_numpy(dtype=float) * 1e-3
    sun = out["sun_jy"].to_numpy(dtype=float)

    with np.errstate(invalid="ignore", divide="ignore"):
        scale = np.where(np.isfinite(sun) & (sun != 0), 1.0 / sun, np.nan)

    if cfg.refl_apply_distance:
        for col in ("r_hel", "r_obs"):
            if col not in out.columns:
                raise KeyError(f"refl_apply_distance needs the {col!r} column")
        scale = scale * out["r_hel"].to_numpy(float) ** 2 * out["r_obs"].to_numpy(float) ** 2

    if cfg.refl_apply_phase:
        if "alpha" not in out.columns:
            raise KeyError("refl_apply_phase needs the 'alpha' column")
        scale = scale * np.power(10.0, 0.4 * cfg.refl_phase_beta * out["alpha"].to_numpy(float))

    out["refl"] = flux_jy * scale
    out["refl_err"] = np.abs(err_jy * scale)   # linear in flux, so errors scale
    return out


def normalise_reflectance(
    df: pd.DataFrame,
    cfg,
    group_cols: Sequence[str] = ("epoch", "ap_label"),
) -> pd.DataFrame:
    """
    Normalise reflectance to unity over a continuum wavelength window.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :func:`add_reflectance`; must also carry ``wl``.
    cfg : Config
        ``refl_norm_window_um`` sets the continuum window.
    group_cols : sequence of str
        Columns defining an independently normalised spectrum.  The default
        normalises each epoch and each aperture size separately, which is what
        makes spectral shapes comparable between them.

    Returns
    -------
    pandas.DataFrame
        A copy with ``refl_norm`` and ``refl_norm_err``.  Groups with no usable
        continuum point get NaN rather than a silently wrong scale.

    Notes
    -----
    The uncertainty on the normalisation constant is propagated, so a spectrum
    normalised on few or noisy continuum points honestly reports larger errors
    everywhere::

        sigma(R/R0)^2 = (sigma_R / R0)^2 + (R * sigma_R0 / R0^2)^2
    """
    out = df.copy()
    if "refl" not in out.columns:
        out = add_reflectance(out, cfg)
    if "wl" not in out.columns:
        raise KeyError("normalise_reflectance needs the 'wl' column")

    lo, hi = cfg.refl_norm_window_um
    out["refl_norm"] = np.nan
    out["refl_norm_err"] = np.nan

    group_cols = [c for c in group_cols if c in out.columns]
    groups = out.groupby(group_cols, dropna=False) if group_cols else [((), out)]

    for key, sub in (groups if group_cols else groups):
        in_win = sub["wl"].between(lo, hi)
        r0, r0_err, n_used = _weighted_clipped_mean(
            sub.loc[in_win, "refl"].to_numpy(float),
            sub.loc[in_win, "refl_err"].to_numpy(float),
        )
        if not np.isfinite(r0) or r0 == 0:
            log.warning("group %s: no usable continuum in %.2f-%.2f um "
                        "(%d candidate points); reflectance left un-normalised",
                        key, lo, hi, int(in_win.sum()))
            continue
        r = sub["refl"].to_numpy(float)
        rerr = sub["refl_err"].to_numpy(float)
        out.loc[sub.index, "refl_norm"] = r / r0
        out.loc[sub.index, "refl_norm_err"] = np.sqrt(
            (rerr / r0) ** 2 + (r * r0_err / r0 ** 2) ** 2)
        log.debug("group %s: R0=%.4g +/- %.2g from %d points", key, r0, r0_err, n_used)

    return out


def bin_spectrum(
    df: pd.DataFrame,
    value_col: str = "refl_norm",
    err_col: str = "refl_norm_err",
    wl_bins: Optional[np.ndarray] = None,
    n_bins: int = 60,
    sigma: float = 3.0,
) -> pd.DataFrame:
    """
    Collapse a scatter of single-exposure points into a binned spectrum.

    Each SPHEREx exposure contributes one wavelength, so a target's spectrum
    arrives as a cloud of hundreds of points with very unequal errors.  Binning
    them with inverse-variance weights and sigma clipping is a far better
    estimator of the underlying spectrum than the raw scatter, and it is the
    form in which band depths should be measured.

    Parameters
    ----------
    df : pandas.DataFrame
        Must carry ``wl``, ``value_col`` and ``err_col``.
    value_col, err_col : str
        Columns to combine.
    wl_bins : ndarray, optional
        Bin edges [um].  Defaults to ``n_bins`` equal bins spanning the data.
    n_bins : int
        Number of bins when ``wl_bins`` is not given.
    sigma : float
        Clipping threshold passed to the weighted mean.

    Returns
    -------
    pandas.DataFrame
        One row per populated bin: ``wl`` (weighted centre), ``wl_lo``,
        ``wl_hi``, ``value``, ``error``, ``n``.
    """
    need = {"wl", value_col, err_col} - set(df.columns)
    if need:
        raise KeyError(f"bin_spectrum needs column(s) {sorted(need)}")

    d = df[np.isfinite(df["wl"]) & np.isfinite(df[value_col])]
    if d.empty:
        return pd.DataFrame(columns=["wl", "wl_lo", "wl_hi", "value", "error", "n"])

    if wl_bins is None:
        wl_bins = np.linspace(float(d["wl"].min()), float(d["wl"].max()) + 1e-9, int(n_bins) + 1)
    wl_bins = np.asarray(wl_bins, dtype=float)

    idx = np.digitize(d["wl"].to_numpy(float), wl_bins) - 1
    rows = []
    for b in range(len(wl_bins) - 1):
        sel = idx == b
        if not sel.any():
            continue
        sub = d[sel]
        mu, mu_err, n = _weighted_clipped_mean(
            sub[value_col].to_numpy(float), sub[err_col].to_numpy(float), sigma=sigma)
        if n == 0:
            continue
        rows.append({"wl": float(sub["wl"].mean()),
                     "wl_lo": float(wl_bins[b]), "wl_hi": float(wl_bins[b + 1]),
                     "value": mu, "error": mu_err, "n": n})
    return pd.DataFrame(rows)
