"""
Gaia contamination bookkeeping.

Gaia stars are deliberately **not** masked.  Painting a disc over every
catalogued source, as the primitive pipeline did, removes real coma flux
wherever a star happens to lie near the comet, and does so with a radius that
grows linearly with stellar brightness -- a large, position-dependent, entirely
uncorrected bias on exactly the measurement the project cares about.

Instead the photometry is performed naively and each measurement carries a
record of what was inside it, so contamination can be filtered (or modelled)
downstream:

============================  ==============================================
``sourceflag``                meaning
============================  ==============================================
``'a'``                       a *bright* blend: the combined G magnitude
                              inside ``r_ap + 2 * PSF_FWHM`` is brighter than
                              ``sourceflag_bright_gmag`` (default 13)
``'b'``                       the blended Gaia flux inside ``r_ap + PSF_FWHM``
                              exceeds ``sourceflag_b_flux_frac`` (default 0.2)
                              of the comet's own predicted flux
``'c'``                       any source lies inside ``r_ap + PSF_FWHM``
``'d'``                       ``SNR < sourceflag_snr_min``
``'0'``                       none of the above
============================  ==============================================

Only the highest-priority flag is stored, in the order a > b > c > d > 0.

Effective magnitude
-------------------
Blends are summed in flux, not in magnitude::

    G_eff = -2.5 * log10( sum_i 10**(-0.4 * G_i) )

so two G = 15 sources give G_eff = 14.25, as intended.

Flag ``'b'`` compares fluxes, not magnitudes::

    10**(-0.4 * G_eff)  >  frac * 10**(-0.4 * vmag)

equivalently ``G_eff < vmag - 2.5*log10(frac)``, i.e. ``G_eff < vmag + 1.75``
for the default ``frac = 0.2``.  It therefore asks the physically meaningful
question -- do the blended stars contribute more than a fifth of what the comet
does? -- rather than comparing two magnitudes on a scale where multiplying one
by 0.2 has no photometric meaning.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

__all__ = ["GaiaPixelField", "SourceFlagResult", "effective_gmag",
           "mag_to_flux", "project_gaia", "evaluate_sourceflag",
           "SOURCEFLAG_COLUMNS"]

SOURCEFLAG_COLUMNS = [
    "n_gaia", "gmag_brightest", "gmag_eff", "gmag_eff_bright",
    "gmag_nearest", "dist_gmag_nearest", "sourceflag",
]


def mag_to_flux(mag) -> np.ndarray:
    """
    Relative flux for a magnitude: ``10**(-0.4 * mag)``.

    The zero point is irrelevant here -- every use is a ratio of two of these --
    but the conversion must happen *before* any comparison, because magnitudes
    are logarithmic and do not scale linearly.
    """
    return np.power(10.0, -0.4 * np.asarray(mag, dtype=np.float64))


def effective_gmag(gmag: np.ndarray) -> float:
    """
    Combined magnitude of a set of sources, summed in flux.

    Parameters
    ----------
    gmag : ndarray
        Individual G magnitudes.  Non-finite entries are ignored.

    Returns
    -------
    float
        ``-2.5 * log10(sum(10**(-0.4 * G)))``, or ``nan`` for an empty set.
    """
    g = np.asarray(gmag, dtype=np.float64)
    g = g[np.isfinite(g)]
    if g.size == 0:
        return float("nan")
    return float(-2.5 * np.log10(np.sum(np.power(10.0, -0.4 * g))))


@dataclass
class GaiaPixelField:
    """Gaia sources of one exposure, already projected into pixel coordinates."""

    x: np.ndarray          #: 0-indexed pixel x
    y: np.ndarray          #: 0-indexed pixel y
    gmag: np.ndarray
    dist: np.ndarray       #: distance from the target centroid [pixel]

    def __len__(self) -> int:
        return int(self.x.size)


@dataclass
class SourceFlagResult:
    """Per-aperture contamination summary; field names match the CSV columns."""

    n_gaia: int
    gmag_brightest: float
    gmag_eff: float
    gmag_eff_bright: float
    gmag_nearest: float
    dist_gmag_nearest: float
    sourceflag: str

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in SOURCEFLAG_COLUMNS}


def project_gaia(
    gaia,
    wcs,
    xcen: float,
    ycen: float,
    shape: Sequence[int],
    pad_pix: float = 30.0,
) -> GaiaPixelField:
    """
    Project a Gaia subset into the pixel frame of one cutout.

    Parameters
    ----------
    gaia : GaiaSubset
        Sources for the target's sky region.
    wcs : astropy.wcs.WCS
        WCS of the science extension (includes the SIP distortion).
    xcen, ycen : float
        Target centroid in 0-indexed pixel coordinates.
    shape : (ny, nx)
        Cutout shape.
    pad_pix : float
        Keep sources up to this far outside the array, so a star just off the
        edge still counts towards a large aperture.

    Returns
    -------
    GaiaPixelField
        Sources sorted by increasing distance from the centroid, which makes
        "nearest source" a simple first-element lookup.

    Notes
    -----
    The catalogue carries no proper motions, so positions are epoch J2016.0.
    At the 6.2 arcsec pixel scale, a typical high-proper-motion star moves well
    under a pixel over the mission, so this is not a limiting error -- but it is
    an approximation.
    """
    ny, nx = int(shape[0]), int(shape[1])
    if gaia is None or len(gaia) == 0:
        e = np.empty(0)
        return GaiaPixelField(e, e.copy(), e.copy(), e.copy())

    # `quiet=True`: the SIP inverse is solved iteratively and does not converge
    # for points far outside the field.  Callers pre-filter to genuine
    # neighbours (see catalog.NeighborIndex), but a stray source must degrade to
    # a non-finite coordinate -- which the mask below drops -- rather than raise
    # and abort the target.
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            x, y = wcs.all_world2pix(gaia.ra, gaia.dec, 0, quiet=True)
        except Exception:                       # noqa: BLE001 - fall back to the linear WCS
            x, y = wcs.wcs_world2pix(gaia.ra, gaia.dec, 0)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    keep = (np.isfinite(x) & np.isfinite(y)
            & (x > -pad_pix) & (x < nx + pad_pix)
            & (y > -pad_pix) & (y < ny + pad_pix))
    x, y, g = x[keep], y[keep], np.asarray(gaia.gmag, dtype=np.float64)[keep]

    dist = np.hypot(x - float(xcen), y - float(ycen))
    order = np.argsort(dist, kind="stable")
    return GaiaPixelField(x[order], y[order], g[order], dist[order])


def evaluate_sourceflag(
    field: GaiaPixelField,
    r_ap_pix: float,
    psf_fwhm_pix: float,
    vmag: float,
    snr: float,
    cfg,
) -> SourceFlagResult:
    """
    Evaluate the contamination flag for one aperture of one exposure.

    Parameters
    ----------
    field : GaiaPixelField
        Projected Gaia sources for this exposure.
    r_ap_pix : float
        Photometric aperture radius [pixel].
    psf_fwhm_pix : float
        PSF FWHM [pixel].  Non-finite values fall back to 0, which reduces both
        test radii to ``r_ap_pix`` rather than discarding the test.
    vmag : float
        Predicted total V magnitude of the comet (JPL/Horizons, via the index).
    snr : float
        Signal-to-noise ratio of this measurement; drives flag ``'d'``.
    cfg : Config
        Supplies the radius multipliers and thresholds.

    Returns
    -------
    SourceFlagResult

    Notes
    -----
    ``'a'`` and ``'b'`` ask different questions and are not nested.  ``'a'`` is
    absolute -- is there a *bright* star nearby at all, within a generous
    ``r_ap + 2*FWHM``, whose wings and ghosts will contaminate the aperture
    regardless of how bright the comet is.  ``'b'`` is relative -- inside the
    tighter ``r_ap + FWHM``, do the blended stars matter compared with *this*
    comet.  A faint comet next to a moderate star trips ``'b'`` but not
    ``'a'``; a bright comet next to a very bright star trips ``'a'`` but might
    not trip ``'b'``.
    """
    fwhm = float(psf_fwhm_pix) if np.isfinite(psf_fwhm_pix) else 0.0
    r_near = float(r_ap_pix) + cfg.sourceflag_near_psf_mult * fwhm
    r_bright = float(r_ap_pix) + cfg.sourceflag_bright_psf_mult * fwhm

    if len(field) == 0:
        near_g = np.empty(0)
        bright_g = np.empty(0)
        g_nearest, d_nearest = float("nan"), float("nan")
    else:
        near_g = field.gmag[field.dist <= r_near]
        bright_g = field.gmag[field.dist <= r_bright]
        # `field` is distance-sorted, so element 0 is the nearest source.
        g_nearest, d_nearest = float(field.gmag[0]), float(field.dist[0])

    n_gaia = int(near_g.size)
    gmag_eff = effective_gmag(near_g)
    gmag_eff_bright = effective_gmag(bright_g)
    gmag_brightest = float(np.nanmin(near_g)) if n_gaia else float("nan")

    # 'b': do the blended stars inside r_near contribute more than
    # `sourceflag_b_flux_frac` of the comet's own predicted flux?  Compared in
    # flux, so the test means what it says.
    flag_b = False
    if n_gaia > 0 and np.isfinite(gmag_eff) and np.isfinite(vmag):
        flag_b = bool(mag_to_flux(gmag_eff)
                      > cfg.sourceflag_b_flux_frac * mag_to_flux(float(vmag)))

    flag = "0"
    if np.isfinite(gmag_eff_bright) and gmag_eff_bright < cfg.sourceflag_bright_gmag:
        flag = "a"
    elif flag_b:
        flag = "b"
    elif n_gaia > 0:
        flag = "c"
    elif np.isfinite(snr) and float(snr) < cfg.sourceflag_snr_min:
        flag = "d"

    return SourceFlagResult(
        n_gaia=n_gaia,
        gmag_brightest=gmag_brightest,
        gmag_eff=gmag_eff,
        gmag_eff_bright=gmag_eff_bright,
        gmag_nearest=g_nearest,
        dist_gmag_nearest=d_nearest,
        sourceflag=flag,
    )
