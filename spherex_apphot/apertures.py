"""
Construction and validity testing of the aperture set for one exposure.

The aperture list is not fixed across the survey: a radius specified in km maps
to a different number of pixels for every exposure, because the plate scale in
km/pixel depends on the observer distance.  This module turns the configured
radii into a concrete, validated list for a single exposure and records *why*
an aperture was rejected.

Validity
--------
An aperture is used only when

``PSF_FWHM_pix <= r_ap_pix < r_in_pix``

The lower bound rejects apertures smaller than the point-spread function, whose
flux is neither resolved nor correctable without an encircled-energy model.  The
upper bound keeps the aperture clear of the background annulus, whose inner
radius is resolved per exposure by :func:`annulus_radii_pix` (150 000 km, floored
at 15 px and capped at 40 px inside the 91-pixel cutout -- the geometric failure
the primitive code had, review item S2, cannot recur).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

__all__ = ["ApertureSpec", "pixel_scale_km", "annulus_radii_pix", "build_apertures"]

#: 1 arcsec in radians.
_ARCSEC_RAD = np.pi / (180.0 * 3600.0)
#: 1 astronomical unit in km (IAU 2012 definition).
_AU_KM = 1.495978707e8


@dataclass(frozen=True)
class ApertureSpec:
    """One aperture, resolved for one exposure."""

    label: str            #: stable identifier, e.g. ``"km10000"`` or ``"pix02.0"``
    kind: str             #: ``"km"`` (fixed physical) or ``"pix"`` (fixed angular)
    r_ap_pix: float       #: radius in pixels -- what photutils is given
    r_ap_km: float        #: radius at the comet in km
    r_ap_arcsec: float    #: radius on the sky in arcsec
    valid: bool
    reject_reason: str = ""   #: ``""``, ``"below_psf_fwhm"`` or ``"beyond_r_in"``


def pixel_scale_km(pix_scale_arcsec: float, r_obs_au: float) -> float:
    """
    Plate scale at the comet, in km per pixel.

    Parameters
    ----------
    pix_scale_arcsec : float
        Angular pixel size [arcsec], from the ``pix_scale`` column / ``PIX-SCL``.
    r_obs_au : float
        Observer-target distance [au].

    Returns
    -------
    float
        km per pixel, or ``nan`` if either input is not finite/positive.

    Notes
    -----
    Small-angle approximation, exact to far better than the pixel size for any
    plausible geometry.
    """
    if not (np.isfinite(pix_scale_arcsec) and np.isfinite(r_obs_au)) or r_obs_au <= 0:
        return float("nan")
    return float(pix_scale_arcsec) * _ARCSEC_RAD * float(r_obs_au) * _AU_KM


def annulus_radii_pix(cfg, kmpp: float) -> tuple[float, float]:
    """
    Sky-annulus radii [pixel] for one exposure.

    With ``cfg.annulus_r_in_km`` set, the inner radius is that physical distance
    converted with the plate scale ``kmpp`` [km/pixel], floored at
    ``cfg.annulus_r_in_pix_min`` and capped at ``cfg.annulus_r_in_pix_max`` (the
    91-pixel cutout allows a ring out to ~45 px); the ring is
    ``cfg.annulus_width_pix`` wide, capped at ``cfg.annulus_r_out_pix_max``.
    Without it, or with a non-finite plate scale, the fixed ``r_in_pix`` /
    ``r_out_pix`` ring is returned.
    """
    if cfg.annulus_r_in_km is None or not np.isfinite(kmpp) or kmpp <= 0:
        return float(cfg.r_in_pix), float(cfg.r_out_pix)
    r_in = float(np.clip(cfg.annulus_r_in_km / kmpp, cfg.annulus_r_in_pix_min, cfg.annulus_r_in_pix_max))
    r_out = float(min(r_in + cfg.annulus_width_pix, cfg.annulus_r_out_pix_max))
    return r_in, r_out


def build_apertures(
    cfg,
    *,
    pix_scale_arcsec: float,
    r_obs_au: float,
    psf_fwhm_pix: float,
) -> List[ApertureSpec]:
    """
    Resolve the configured aperture radii for one exposure.

    Parameters
    ----------
    cfg : Config
        Supplies ``aperture_radii_km``, ``aperture_radii_pix``, ``r_in_pix`` and
        the two ``skip_ap_*`` switches.
    pix_scale_arcsec : float
        Angular pixel size for this exposure.
    r_obs_au : float
        Observer-target distance for this exposure.
    psf_fwhm_pix : float
        PSF FWHM in pixels (``PSF_FWHM`` / ``pix_scale``).  Non-finite values
        disable the lower bound rather than rejecting everything.

    Returns
    -------
    list of ApertureSpec
        Every configured aperture, in a stable order: the km-defined radii
        ascending, then the pixel-defined radii.  Invalid entries are included
        with ``valid=False`` so the caller can count and log them; whether they
        reach the output file is decided by ``cfg.drop_invalid_apertures``.
    """
    kmpp = pixel_scale_km(pix_scale_arcsec, r_obs_au)
    r_in_pix, _ = annulus_radii_pix(cfg, kmpp)
    specs: List[ApertureSpec] = []

    def _check(r_pix: float) -> tuple[bool, str]:
        if not np.isfinite(r_pix) or r_pix <= 0:
            return False, "not_finite"
        if cfg.skip_ap_below_psf_fwhm and np.isfinite(psf_fwhm_pix) and r_pix < psf_fwhm_pix:
            return False, "below_psf_fwhm"
        if cfg.skip_ap_beyond_r_in and r_pix >= r_in_pix:
            return False, "beyond_r_in"
        return True, ""

    for r_km in cfg.aperture_radii_km:
        r_pix = float(r_km) / kmpp if np.isfinite(kmpp) and kmpp > 0 else float("nan")
        ok, why = _check(r_pix)
        specs.append(ApertureSpec(
            label=f"km{int(round(float(r_km)))}", kind="km",
            r_ap_pix=r_pix, r_ap_km=float(r_km),
            r_ap_arcsec=r_pix * pix_scale_arcsec if np.isfinite(r_pix) else float("nan"),
            valid=ok, reject_reason=why,
        ))

    for r_pix_cfg in cfg.aperture_radii_pix:
        r_pix = float(r_pix_cfg)
        ok, why = _check(r_pix)
        specs.append(ApertureSpec(
            label=f"pix{r_pix:04.1f}", kind="pix",
            r_ap_pix=r_pix,
            r_ap_km=r_pix * kmpp if np.isfinite(kmpp) else float("nan"),
            r_ap_arcsec=r_pix * pix_scale_arcsec,
            valid=ok, reject_reason=why,
        ))

    return specs
