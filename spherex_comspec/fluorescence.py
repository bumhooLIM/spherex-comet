"""
Layer 3-4 inputs from the reconstructed fluorescence database (``data/fluorescence``).

The database is built by ``notebooks/comspec/fluorescence_gfm/build_fluorescence_db.py`` with the General
Fluorescence Model of Villanueva et al. (2011, 2012) -- HITRAN 2020 line lists pumped by a
Kurucz-continuum x Fraunhofer-line solar spectrum, with level-by-level cascade -- and is validated
against the published GSFC g-factors to ~10 % (``doc/comspec/fluorescence_database.md``).  It gives, per
species, the emission spectral density

    g_lambda(lambda; T_rot)   [photons s^-1 molecule^-1 um^-1]  at 1 au, v_h = 0,

on an R = 5000 grid over 0.7-5.0 um for T_rot = 30, 50, 70, 100, 130 K, and for CO the v(1-0)
g-factor versus heliocentric velocity (the Swings effect: the cometary CO lines sit on the solar
CO Fraunhofer lines at v_h = 0, and g(CO) is ~25 % larger once |v_h| exceeds ~10 km/s; no other
band changes by more than 1 %).

This module reads those tables once and serves three things to :mod:`gasmodel`:

* :func:`g_lambda` -- g_lambda rebinned (flux-conserving) onto any wavelength grid and
  interpolated linearly in T_rot;
* :func:`swings_factor` -- g(v_h) / g(0) for CO, exactly 1 for every other species and for an
  unknown (NaN) velocity;
* :func:`band_total` -- the integral of g_lambda over a wavelength interval, i.e. the band
  g-factor the fit effectively uses.

Positive v_h is motion away from the Sun (post-perihelion).
"""

from __future__ import annotations

import json
import os
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from . import directory as _dir

__all__ = ["FLUOR_DIR", "SWINGS_SPECIES", "available", "describe", "g_lambda", "swings_factor",
           "band_total", "rebin_density", "profile_temperatures"]

#: the database lives with the project (``data/fluorescence/``), not with a run's
#: (overridable) data tree; ``COMSPEC_FLUOR_DIR`` is kept as an alias of
#: ``COMSPEC_FLUORESCENCE_DIR`` for older shell scripts.
FLUOR_DIR: Path = Path(os.environ["COMSPEC_FLUOR_DIR"]).expanduser() if os.environ.get("COMSPEC_FLUOR_DIR") \
    else _dir.FLUORESCENCE_DIR
#: species whose g-factor is scaled with heliocentric velocity
SWINGS_SPECIES = ("CO",)


def available(species=("H2O", "CO2", "CO")) -> bool:
    return all((FLUOR_DIR / "profiles" / f"{s}_gprofile.csv").is_file() for s in species) \
        and (FLUOR_DIR / "co_swings.csv").is_file()


def describe() -> dict:
    p = FLUOR_DIR / "build_meta.json"
    return json.loads(p.read_text()) if p.is_file() else {}


@lru_cache(maxsize=None)
def _profile(species: str):
    p = FLUOR_DIR / "profiles" / f"{species}_gprofile.csv"
    if not p.is_file():
        raise FileNotFoundError(f"{p}: build the fluorescence database first "
                                "(notebooks/comspec/fluorescence_gfm/build_fluorescence_db.py)")
    d = pd.read_csv(p)
    temps = sorted(int(c[7:10]) for c in d.columns if c.startswith("g_lam_T") and c.endswith("K"))
    lam = d["lam_um"].to_numpy(dtype=float)
    g = np.vstack([d[f"g_lam_T{T:03d}K"].to_numpy(dtype=float) for T in temps])
    return lam, np.asarray(temps, dtype=float), g


def profile_temperatures(species: str = "H2O") -> np.ndarray:
    return _profile(species)[1]


def _edges(lam: np.ndarray) -> np.ndarray:
    """Cell edges of a monotonic grid: geometric midpoints, extrapolated at both ends."""
    lam = np.asarray(lam, dtype=float)
    if lam.size == 1:
        return np.array([lam[0] * (1 - 1e-4), lam[0] * (1 + 1e-4)])
    mid = np.sqrt(lam[:-1] * lam[1:])
    return np.concatenate([[lam[0] ** 2 / mid[0]], mid, [lam[-1] ** 2 / mid[-1]]])


def rebin_density(lam_src, y_src, lam_dst) -> np.ndarray:
    """
    Flux-conserving rebin of a spectral density from one grid to another.

    The integral of ``y`` over any interval is preserved: the cumulative integral on the source
    cells is interpolated at the destination cell edges and differenced.  Destination cells
    outside the source range receive zero.
    """
    es, ed = _edges(np.asarray(lam_src, float)), _edges(np.asarray(lam_dst, float))
    cum = np.concatenate([[0.0], np.cumsum(np.asarray(y_src, float) * np.diff(es))])
    c = np.interp(ed, es, cum, left=0.0, right=cum[-1])
    return np.diff(c) / np.diff(ed)


def _at_temperature(species: str, T_rot: float):
    lam, temps, g = _profile(species)
    T = float(np.clip(T_rot, temps[0], temps[-1]))
    j = int(np.searchsorted(temps, T))
    if j == 0 or temps[j - 1] == T:
        return lam, g[max(j - 1, 0) if temps[max(j - 1, 0)] == T else j]
    w = (T - temps[j - 1]) / (temps[j] - temps[j - 1])
    return lam, (1 - w) * g[j - 1] + w * g[j]


def g_lambda(species: str, lam_um, T_rot: float = 70.0) -> np.ndarray:
    """g_lambda [photons s^-1 molecule^-1 um^-1] of ``species`` on the grid ``lam_um`` at T_rot."""
    lam_s, y = _at_temperature(species, T_rot)
    lam = np.atleast_1d(np.asarray(lam_um, dtype=float))
    return rebin_density(lam_s, y, lam)


def band_total(species: str, lo_um: float, hi_um: float, T_rot: float = 70.0) -> float:
    """Integral of g_lambda over [lo, hi] -- the g-factor of everything the species emits there."""
    lam, y = _at_temperature(species, T_rot)
    e = _edges(lam)
    dl = np.diff(e)
    m = (lam >= lo_um) & (lam <= hi_um)
    return float(np.sum(y[m] * dl[m]))


@lru_cache(maxsize=None)
def _swings(species: str):
    p = FLUOR_DIR / "co_swings.csv"
    if not p.is_file():
        raise FileNotFoundError(f"{p}: build the fluorescence database first")
    d = pd.read_csv(p)
    temps = sorted(int(c[7:10]) for c in d.columns if c.startswith("ratio_T") and c.endswith("K"))
    r = np.vstack([d[f"ratio_T{T:03d}K"].to_numpy(dtype=float) for T in temps])
    return d["v_h_kms"].to_numpy(dtype=float), np.asarray(temps, dtype=float), r


def swings_factor(species: str, v_h_kms, T_rot: float = 70.0):
    """
    ``g(v_h) / g(v_h = 0)`` of the species' strongest band, interpolated linearly in ``T_rot``
    between the tabulated 30, 70 and 130 K curves (they differ by up to +-4 % where the ratio
    rises, 2-10 km/s, and by ~1 % on the plateau).  Exactly 1 for species outside
    :data:`SWINGS_SPECIES` and for a non-finite velocity.  Velocities beyond the table
    (+-60 km/s) take the edge value.
    """
    v = np.asarray(v_h_kms, dtype=float)
    scalar = v.ndim == 0
    v = np.atleast_1d(v)
    out = np.ones_like(v)
    if species in SWINGS_SPECIES:
        vg, temps, r = _swings(species)
        T = float(np.clip(T_rot, temps[0], temps[-1]))
        j = int(np.searchsorted(temps, T))
        if j == 0 or temps[j - 1] == T:
            ratio = r[j] if (j < len(temps) and temps[j] == T) else r[max(j - 1, 0)]
        else:
            w = (T - temps[j - 1]) / (temps[j] - temps[j - 1])
            ratio = (1 - w) * r[j - 1] + w * r[j]
        ok = np.isfinite(v)
        out[ok] = np.interp(np.clip(v[ok], vg[0], vg[-1]), vg, ratio)
    return float(out[0]) if scalar else out
