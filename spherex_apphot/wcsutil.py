"""
World-coordinate helpers.

Two independent WCS paths exist in this project and it is worth being explicit
about which is which:

* **From the file.**  The science extension carries the full, already
  cutout-shifted WCS including the SIP distortion.  Anything that has the FITS
  open should use :func:`wcs_from_header` -- it is the authoritative version.
* **From the index.**  :func:`reconstruct_cutout_wcs` rebuilds an equivalent
  WCS from the Parquet row alone, which lets the notebook draw sky coverage for
  thousands of exposures without opening thousands of files.

The two must agree; :func:`compare_wcs` exists so the notebook can assert that
they do rather than assume it.
"""

from __future__ import annotations

import re
import warnings
from typing import Any, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

from .logging_utils import get_logger

__all__ = ["wcs_from_header", "reconstruct_cutout_wcs", "footprint_radec",
           "north_angle_deg", "compare_wcs"]

log = get_logger("wcsutil")

#: SIP coefficient keywords: ``A_ORDER``, ``A_i_j``, and the ``AP``/``BP``
#: inverse sets.  The index also carries ``PV1_*``/``PV2_*``, which are inert
#: for a ``TAN-SIP`` projection and are deliberately not copied.
_SIP_RE = re.compile(r"^[AB]P?(_ORDER|_\d+_\d+)$")

_BASE_KEYS = ("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
              "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CDELT1", "CDELT2",
              "CUNIT1", "CUNIT2", "LONPOLE", "LATPOLE", "RADESYS", "WCSAXES")


def wcs_from_header(header: fits.Header) -> WCS:
    """Build a WCS from a science-extension header (SIP included)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header, relax=True)


def _get(row: Any, key: str, default=None):
    """Read ``key`` from a namedtuple, Series or mapping."""
    if hasattr(row, "_fields") or hasattr(row, key):
        val = getattr(row, key, default)
        if val is not default:
            return val
    try:
        return row[key]
    except (KeyError, TypeError, IndexError):
        return default


def reconstruct_cutout_wcs(row: Any) -> WCS:
    """
    Rebuild a cutout WCS from one row of the cutout index.

    Parameters
    ----------
    row : namedtuple, Series or mapping
        A row of ``db_filtered.parq``.  Must carry the frame WCS keywords plus
        ``ltv1``/``ltv2`` and ``cutout_size``.

    Returns
    -------
    astropy.wcs.WCS

    Notes
    -----
    The reference pixel is shifted by the cutout origin (``CRPIX + LTV``), which
    is the same transformation the FITS exporter applies, so the result should
    match the header WCS to round-off.  Use :func:`compare_wcs` to check.

    ``itertuples`` renames columns that are not valid Python identifiers (for
    instance ``DATE-OBS`` becomes ``_57``), which is why every lookup goes
    through :func:`_get` and why only identifier-safe keywords are read here.
    """
    hdr = fits.Header()
    hdr["WCSAXES"] = int(_get(row, "WCSAXES", 2) or 2)
    hdr["RADESYS"] = _get(row, "RADESYS", "ICRS") or "ICRS"
    for key in ("CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2"):
        val = _get(row, key)
        if val is not None:
            hdr[key] = val
    for key in ("CRVAL1", "CRVAL2", "LONPOLE", "LATPOLE",
                "CDELT1", "CDELT2", "PC1_1", "PC1_2", "PC2_1", "PC2_2"):
        val = _get(row, key)
        if val is not None and np.isfinite(val):
            hdr[key] = float(val)

    crpix1, crpix2 = _get(row, "CRPIX1"), _get(row, "CRPIX2")
    ltv1, ltv2 = _get(row, "ltv1", 0.0), _get(row, "ltv2", 0.0)
    if crpix1 is None or crpix2 is None:
        raise ValueError("row has no CRPIX1/CRPIX2; cannot reconstruct a WCS")
    hdr["CRPIX1"] = float(crpix1) + float(ltv1 or 0.0)
    hdr["CRPIX2"] = float(crpix2) + float(ltv2 or 0.0)

    size = int(_get(row, "cutout_size", 0) or 0)
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = size
    hdr["NAXIS2"] = size

    fields = getattr(row, "_fields", None)
    if fields is None:
        fields = list(getattr(row, "index", [])) or list(getattr(row, "keys", lambda: [])())
    for name in fields:
        if not _SIP_RE.match(str(name)):
            continue
        val = _get(row, name)
        if val is not None and np.isfinite(val):
            hdr[str(name)] = float(val)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(hdr, relax=True)


def footprint_radec(wcs: WCS, nx: int, ny: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sky coordinates of a cutout's four corners, closed into a polygon.

    Returns
    -------
    ra, dec : ndarray
        Five points (the first repeated last).  RA is unwrapped so that a
        footprint straddling 0h plots as a single polygon instead of a streak
        across the whole axis.
    """
    x = np.array([0, nx, nx, 0, 0], dtype=float)
    y = np.array([0, 0, ny, ny, 0], dtype=float)
    sky = wcs.pixel_to_world(x, y)
    ra = np.asarray(sky.ra.deg, dtype=float)
    dec = np.asarray(sky.dec.deg, dtype=float)
    if np.nanmax(ra) - np.nanmin(ra) > 180.0:
        ra = np.where(ra > 180.0, ra - 360.0, ra)
    return ra, dec


def north_angle_deg(wcs: WCS, x: float, y: float) -> float:
    """
    Position angle of celestial north in the pixel frame, at ``(x, y)``.

    Measured in degrees counter-clockwise from the ``+y`` pixel axis, which is
    the rotation :mod:`stacking` must undo to put north up.

    Notes
    -----
    Evaluated by finite difference through the full WCS (SIP included) rather
    than from the PC matrix, so local distortion is accounted for.
    """
    step = 0.5
    c0 = wcs.pixel_to_world(x, y)
    c1 = wcs.pixel_to_world(x, y + step)
    dra = (c1.ra.deg - c0.ra.deg)
    dra = (dra + 180.0) % 360.0 - 180.0
    dra *= np.cos(np.radians(c0.dec.deg))
    ddec = c1.dec.deg - c0.dec.deg
    return float(np.degrees(np.arctan2(dra, ddec)))


def compare_wcs(wcs_a: WCS, wcs_b: WCS, shape: Sequence[int],
                n: int = 7) -> Tuple[float, float]:
    """
    Compare two WCS solutions on a grid, in pixels.

    Returns
    -------
    max_offset_pix, rms_offset_pix : float
        Displacement obtained by mapping a grid through ``wcs_a`` to the sky and
        back through ``wcs_b``.
    """
    ny, nx = int(shape[0]), int(shape[1])
    gx, gy = np.meshgrid(np.linspace(0, nx - 1, n), np.linspace(0, ny - 1, n))
    sky = wcs_a.pixel_to_world(gx.ravel(), gy.ravel())
    bx, by = wcs_b.world_to_pixel(sky)
    d = np.hypot(bx - gx.ravel(), by - gy.ravel())
    return float(np.nanmax(d)), float(np.sqrt(np.nanmean(d ** 2)))
