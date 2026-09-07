"""Orbital geometry: perihelion elements and the r_h <-> time relation.

An Af-rho vs heliocentric distance plot is really a plot against orbital
*phase*, and r_h alone cannot express it: a comet passes every r_h twice, once
inbound and once outbound, and never goes below perihelion at all.  Plotting
signed r_h leaves a hole between -q and +q where no data can exist and puts a
meaningless zero in the middle of it.

Measuring from perihelion instead -- ``r_h - q``, negative inbound -- puts both
legs on one continuous axis that starts at a real physical point.  The time from
perihelion follows from Kepler's equation, so a date axis can be laid over the
same (linear in r_h) ticks even though time is strongly non-linear in r_h.

All three conic cases are handled: elliptic, near-parabolic and hyperbolic.
"""

from __future__ import annotations

import logging

import numpy as np
from astroquery.jplhorizons import Horizons

__all__ = ["GAUSS_K", "fetch_elements", "time_from_perihelion",
           "rh_from_time", "PerihelionInfo"]

log = logging.getLogger(__name__)

#: Gaussian gravitational constant, AU^(3/2) / day.
GAUSS_K = 0.01720209895

#: Within this of 1, treat the orbit as parabolic (Kepler's equation degenerates).
_PARABOLIC_TOL = 1e-3

_CACHE: dict = {}


class PerihelionInfo(dict):
    """Perihelion elements, with attribute access for readability."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    def __repr__(self):
        return (f"PerihelionInfo(q={self.get('q')!r}, e={self.get('e')!r}, "
                f"Tp_jd={self.get('Tp_jd')!r})")


def fetch_elements(target_id, epoch_jd=None, use_cache=True):
    """Heliocentric osculating elements from Horizons.

    Parameters
    ----------
    target_id : str or int
        Whatever :func:`ztfcomet.horizons.resolve_target_id` produced.
    epoch_jd : float, optional
        Epoch at which to osculate; defaults to J2000.  For a comet with
        significant non-gravitational forces the elements drift, so pass an
        epoch near the observations.

    Returns
    -------
    PerihelionInfo or None
        ``q`` (au), ``e``, ``Tp_jd``, ``n`` (deg/day), ``a`` (au), ``P`` (days).
        ``None`` when the query fails -- callers fall back to plain r_h.
    """
    key = (str(target_id), None if epoch_jd is None else round(float(epoch_jd), 1))
    if use_cache and key in _CACHE:
        return _CACHE[key]

    try:
        table = Horizons(id=target_id, location="@sun",
                         epochs=[epoch_jd or 2451545.0]).elements()
        row = table[0]
        info = PerihelionInfo(
            targetname=str(row["targetname"]),
            q=float(row["q"]), e=float(row["e"]), Tp_jd=float(row["Tp_jd"]),
            n=float(row["n"]) if "n" in table.columns else np.nan,
            a=float(row["a"]) if "a" in table.columns else np.nan,
            P=float(row["P"]) if "P" in table.columns else np.nan,
        )
    except Exception as exc:                                    # noqa: BLE001
        log.warning("Could not fetch elements for %r: %s", target_id, exc)
        info = None

    if use_cache:
        _CACHE[key] = info
    return info


def time_from_perihelion(rh, q, e, signed_by=None):
    """Days between perihelion and the epoch at heliocentric distance *rh*.

    Solves Kepler's equation for the eccentric (or hyperbolic) anomaly at *rh*
    and converts to time.  The result is the **magnitude** unless *signed_by* is
    given.

    Parameters
    ----------
    rh : array_like
        Heliocentric distance in au.  Values below perihelion return 0 and, for
        a closed orbit, values beyond aphelion return NaN -- the comet reaches
        neither, and a saturated value would put a wrong date on an axis.
    q : float
        Perihelion distance, au.
    e : float
        Eccentricity.
    signed_by : array_like, optional
        Anything whose sign marks the orbital leg -- typically ``r_rate``.
        Negative (inbound) gives a negative time.

    Returns
    -------
    ndarray
        Days from perihelion.

    Notes
    -----
    Elliptic and hyperbolic cases use the standard Kepler and hyperbolic-Kepler
    relations; ``|e - 1| < 1e-3`` uses Barker's equation, since both forms
    degenerate at the parabolic limit.
    """
    rh = np.atleast_1d(np.asarray(rh, dtype=float))
    out = np.full(rh.shape, np.nan)

    # Inside perihelion is unreachable. Clamp to 0 rather than NaN so an axis
    # tick sitting fractionally below q still gets a label. Do this BEFORE the
    # early return, or an all-inside input comes back as NaN.
    finite = np.isfinite(rh)
    out[finite & (rh < q)] = 0.0

    reachable = finite & (rh >= q)
    if e < 1.0 - _PARABOLIC_TOL:
        # Beyond aphelion is unreachable too; leave those NaN rather than
        # silently saturating every one of them at half a period.
        aphelion = q * (1.0 + e) / (1.0 - e)
        reachable &= rh <= aphelion * (1.0 + 1e-9)

    if not np.any(reachable):
        return out

    r = rh[reachable]
    if abs(e - 1.0) < _PARABOLIC_TOL:
        # Barker's equation: tan(nu/2) = sqrt(r/q - 1)
        d = np.sqrt(np.maximum(r / q - 1.0, 0.0))
        dt = (np.sqrt(2.0) * q ** 1.5 / GAUSS_K) * (d + d ** 3 / 3.0)
    elif e < 1.0:
        a = q / (1.0 - e)
        cos_e = np.clip((1.0 - r / a) / e, -1.0, 1.0)
        ecc_anom = np.arccos(cos_e)
        mean_anom = ecc_anom - e * np.sin(ecc_anom)
        dt = mean_anom / (GAUSS_K / a ** 1.5)
    else:
        a = q / (e - 1.0)
        cosh_h = np.clip((1.0 + r / a) / e, 1.0, None)
        hyp_anom = np.arccosh(cosh_h)
        mean_anom = e * np.sinh(hyp_anom) - hyp_anom
        dt = mean_anom / (GAUSS_K / a ** 1.5)

    out[reachable] = dt

    if signed_by is not None:
        sign = np.where(np.asarray(signed_by, dtype=float) < 0, -1.0, 1.0)
        out = out * sign
    return out


def rh_from_time(days, q, e, n_steps=4096):
    """Heliocentric distance at *days* from perihelion.

    Inverse of :func:`time_from_perihelion`, by monotonic interpolation.  Used
    when a figure wants date ticks rather than distance ticks.
    """
    days = np.atleast_1d(np.asarray(days, dtype=float))
    span = np.nanmax(np.abs(days)) if np.any(np.isfinite(days)) else 1.0
    # Build r -> t over a generous range, then invert.
    r_grid = q + np.linspace(0.0, max(50.0, 10.0 * q), n_steps) ** 1.5
    t_grid = time_from_perihelion(r_grid, q, e)
    good = np.isfinite(t_grid)
    if good.sum() < 2 or t_grid[good].max() < span:
        r_grid = q + np.linspace(0.0, 500.0, n_steps) ** 1.5
        t_grid = time_from_perihelion(r_grid, q, e)
        good = np.isfinite(t_grid)
    return np.interp(np.abs(days), t_grid[good], r_grid[good])
