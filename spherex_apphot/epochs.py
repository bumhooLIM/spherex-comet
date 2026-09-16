"""
Grouping exposures into observing epochs.

SPHEREx revisits a field in short bursts separated by long gaps, so the natural
unit for combining a comet's exposures is the visit, not the whole mission.

Naming
------
The primitive code called this quantity ``phase``, which collides with the solar
*phase angle* ``alpha`` carried in the same table -- the standard meaning of
"phase" in comet photometry.  The column is called ``epoch`` here (review item
"naming collision"), and the DataFrame also gains ``epoch_n``, the number of
exposures in that epoch, which the stacking step needs.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from .logging_utils import get_logger

__all__ = ["group_epochs", "epoch_table"]

log = get_logger("epochs")


def group_epochs(
    df: pd.DataFrame,
    gap_days: float = 28.0,
    key_dateobs: str = "DATE-OBS",
) -> pd.DataFrame:
    """
    Sort exposures chronologically and label contiguous observing epochs.

    Parameters
    ----------
    df : pandas.DataFrame
        Cutout index rows for one target.
    gap_days : float
        A gap strictly greater than this starts a new epoch.
    key_dateobs : str
        Column holding the observation time.  ``jd_utc`` is used as a fallback
        if this column is absent or entirely unparseable.

    Returns
    -------
    pandas.DataFrame
        A copy sorted by time with a fresh ``RangeIndex`` and two new columns,
        ``epoch`` (1-based) and ``epoch_n``.

    Raises
    ------
    ValueError
        If no usable time column exists.
    """
    if df.empty:
        out = df.copy()
        out["epoch"] = pd.Series(dtype="int64")
        out["epoch_n"] = pd.Series(dtype="int64")
        return out

    out = df.copy()

    times: Optional[pd.Series] = None
    if key_dateobs in out.columns:
        times = pd.to_datetime(out[key_dateobs], errors="coerce", utc=True)
        if times.isna().all():
            log.warning("%r is present but unparseable; falling back to jd_utc", key_dateobs)
            times = None
    if times is None:
        if "jd_utc" not in out.columns:
            raise ValueError(f"neither {key_dateobs!r} nor 'jd_utc' available for epoch grouping")
        times = pd.to_datetime(out["jd_utc"], origin="julian", unit="D", errors="coerce", utc=True)
    if times.isna().all():
        raise ValueError("no parseable observation times; cannot group epochs")

    n_bad = int(times.isna().sum())
    if n_bad:
        # Undated rows sort last and land in the final epoch; they are flagged
        # rather than dropped so nothing disappears silently.
        log.warning("%d row(s) have no parseable observation time", n_bad)

    order = np.argsort(times.values, kind="stable")
    out = out.iloc[order].reset_index(drop=True)
    times = times.iloc[order].reset_index(drop=True)

    gap = times.diff() > pd.Timedelta(days=float(gap_days))
    out["epoch"] = gap.cumsum().astype("int64") + 1
    out["epoch_n"] = out.groupby("epoch")["epoch"].transform("size").astype("int64")

    log.info("%d exposures grouped into %d epoch(s) at a %g-day gap",
             len(out), out["epoch"].nunique(), gap_days)
    return out


def epoch_table(df: pd.DataFrame, key_dateobs: str = "DATE-OBS") -> pd.DataFrame:
    """
    Summarise each epoch: exposure count, date range, wavelength and geometry.

    Returned for logging and for the notebook; not part of the photometry
    output.
    """
    if df.empty or "epoch" not in df.columns:
        return pd.DataFrame()

    t = pd.to_datetime(df[key_dateobs], errors="coerce", utc=True) if key_dateobs in df else pd.NaT
    work = df.assign(_t=t)
    agg = work.groupby("epoch").agg(
        n_exposures=("epoch", "size"),
        date_start=("_t", "min"),
        date_end=("_t", "max"),
        wl_min=("wl", "min"),
        wl_max=("wl", "max"),
        r_hel_mean=("r_hel", "mean"),
        r_obs_mean=("r_obs", "mean"),
        alpha_mean=("alpha", "mean"),
        vmag_mean=("vmag", "mean"),
    )
    return agg.reset_index()
