"""
Robust band stacking (review item S8).

The primitive stacker took a plain ``nanmedian`` of raw cutouts aligned only to
the nearest whole pixel.  Three things were wrong with that, and all three are
addressed here.

1. **No background subtraction.**  The zodiacal foreground dominates a SPHEREx
   frame and varies strongly between visits, so the median was largely a median
   of sky levels rather than of comet.  Each frame is now background-subtracted
   with its own annulus median before entering the stack.
2. **Whole-pixel registration.**  Rounding the centroid introduces up to half a
   pixel of jitter against a PSF that is only ~1 pixel across.  Frames are now
   shifted to a common sub-pixel centre with NaN-safe bilinear interpolation.
3. **Ignored field rotation.**  SPHEREx revisits a field at different roll
   angles, so stacking in detector coordinates smears any real morphology.
   ``align_north`` rotates every frame to a common north-up orientation, which
   also makes the stack's own WCS meaningful.

The combination itself is a sigma-clipped **median** (``Config.stack_combine``;
``"mean"`` is available), and the result carries a per-pixel count map, so a
pixel built from two frames is distinguishable from one built from forty.  The
median is the right default for these fields: sigma clipping removes obvious
outliers, but in a crowded field a moderately bright star that lands on the same
stamp position in a handful of frames survives clipping and pulls a mean, while
the median ignores it as long as it is a minority.

NaN handling
------------
Interpolating an array that contains NaN spreads the NaN over the interpolation
kernel.  Every resampling step here therefore interpolates a zero-filled copy
and a companion weight image, then divides -- so bad pixels neither contaminate
their neighbours nor bias the result low.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field as _field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clip
from scipy import ndimage

from .logging_utils import get_logger

__all__ = ["StackInput", "StackResult", "shift_subpixel", "rotate_nan_safe",
           "extract_stamp", "stack_frames", "write_stack_fits"]

log = get_logger("stacking")

#: Weight below which a resampled pixel is considered unreliable and set to NaN.
_MIN_WEIGHT = 0.5


@dataclass
class StackInput:
    """One frame's contribution to a stack."""

    sci: np.ndarray
    badpix: np.ndarray
    xcen: float
    ycen: float
    sky_median: float = 0.0
    north_deg: float = 0.0        #: position angle of north, from ``wcsutil``
    label: str = ""


@dataclass
class StackResult:
    """Output of :func:`stack_frames`."""

    image: np.ndarray             #: sigma-clipped median or mean [mJy/pixel]
    count: np.ndarray             #: frames contributing to each pixel
    error: np.ndarray             #: sampling error of that statistic, per pixel
    n_frames: int
    radius: int
    labels: List[str] = _field(default_factory=list)
    north_aligned: bool = False
    combine: str = "median"       #: statistic stored in ``image``
    pix_scale_arcsec: float = float("nan")

    @property
    def shape(self) -> Tuple[int, int]:
        return self.image.shape


def shift_subpixel(img: np.ndarray, dy: float, dx: float, order: int = 1) -> np.ndarray:
    """
    Shift an image by a fractional number of pixels, NaN-safe.

    Parameters
    ----------
    img : ndarray
        2-D image; NaN marks missing data.
    dy, dx : float
        Shift in pixels (``scipy.ndimage`` convention: positive moves content
        toward larger indices).
    order : int
        Spline order; 1 (bilinear) keeps the noise correlation local and is the
        right choice for a ~1-pixel PSF.

    Returns
    -------
    ndarray
        Shifted image with NaN where the interpolated weight fell below 0.5,
        i.e. where the value would have been dominated by missing data.
    """
    good = np.isfinite(img)
    if not good.any():
        return np.full_like(img, np.nan)
    filled = np.where(good, img, 0.0)
    num = ndimage.shift(filled, (dy, dx), order=order, mode="constant", cval=0.0, prefilter=False)
    den = ndimage.shift(good.astype(float), (dy, dx), order=order, mode="constant",
                        cval=0.0, prefilter=False)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > _MIN_WEIGHT, num / den, np.nan)
    return out


def rotate_nan_safe(img: np.ndarray, angle_deg: float, order: int = 1) -> np.ndarray:
    """
    Rotate an image about its centre, NaN-safe, preserving shape.

    Parameters
    ----------
    img : ndarray
    angle_deg : float
        Counter-clockwise rotation applied to the image content.
    order : int
        Spline order.
    """
    if not np.isfinite(angle_deg) or abs(angle_deg) < 1e-6:
        return img
    good = np.isfinite(img)
    filled = np.where(good, img, 0.0)
    kw = dict(reshape=False, order=order, mode="constant", cval=0.0, prefilter=False)
    num = ndimage.rotate(filled, angle_deg, **kw)
    den = ndimage.rotate(good.astype(float), angle_deg, **kw)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > _MIN_WEIGHT, num / den, np.nan)


def extract_stamp(
    img: np.ndarray,
    xcen: float,
    ycen: float,
    radius: int,
    subpixel: bool = True,
) -> np.ndarray:
    """
    Cut a ``(2*radius+1)`` square stamp centred on ``(xcen, ycen)``.

    Parameters
    ----------
    img : ndarray
        Source image; NaN marks missing data.
    xcen, ycen : float
        Requested centre, 0-indexed.  The integer part selects the stamp and the
        fractional part is removed by a sub-pixel shift when ``subpixel``.
    radius : int
        Half-width in pixels.
    subpixel : bool
        Apply the fractional-pixel correction.

    Returns
    -------
    ndarray
        Stamp of shape ``(2*radius+1, 2*radius+1)``, NaN-padded where the
        request runs past the array edge.
    """
    radius = int(radius)
    size = 2 * radius + 1
    out = np.full((size, size), np.nan, dtype=np.float64)
    if not (np.isfinite(xcen) and np.isfinite(ycen)):
        return out

    xi, yi = int(np.floor(xcen + 0.5)), int(np.floor(ycen + 0.5))
    h, w = img.shape
    x0, x1 = xi - radius, xi + radius + 1
    y0, y1 = yi - radius, yi + radius + 1
    sx0, sx1 = max(0, x0), min(w, x1)
    sy0, sy1 = max(0, y0), min(h, y1)
    if sx0 >= sx1 or sy0 >= sy1:
        return out
    out[sy0 - y0:sy0 - y0 + (sy1 - sy0), sx0 - x0:sx0 - x0 + (sx1 - sx0)] = img[sy0:sy1, sx0:sx1]

    if subpixel:
        # Move the true centroid onto the exact stamp centre.
        out = shift_subpixel(out, dy=(yi - float(ycen)), dx=(xi - float(xcen)))
    return out


def stack_frames(
    items: Sequence[StackInput],
    radius: int,
    *,
    subtract_sky: bool = True,
    subpixel: bool = True,
    align_north: bool = False,
    sigma: float = 3.0,
    maxiters: int = 5,
    combine: str = "median",
    pix_scale_arcsec: float = float("nan"),
) -> Optional[StackResult]:
    """
    Combine frames into one robust stack.

    Parameters
    ----------
    items : sequence of StackInput
        Frames to combine.  Bad pixels are set to NaN before resampling, so they
        drop out of both the stack and the count map.
    radius : int
        Stamp half-width [pixel].
    subtract_sky : bool
        Remove each frame's own background level before combining.
    subpixel : bool
        Register to a common sub-pixel centre.
    align_north : bool
        Rotate each frame so celestial north points up.
    sigma, maxiters : float, int
        Sigma-clipping parameters applied before combining.
    combine : {'median', 'mean'}
        Statistic used on the surviving pixels.  ``'median'`` is the default and
        the more robust choice against residual stars.
    pix_scale_arcsec : float
        Recorded on the result for the output WCS.

    Returns
    -------
    StackResult or None
        ``None`` when no frame yielded any usable pixel.
    """
    stamps: List[np.ndarray] = []
    labels: List[str] = []
    for it in items:
        img = np.where(it.badpix, np.nan, np.asarray(it.sci, dtype=np.float64))
        if subtract_sky and np.isfinite(it.sky_median):
            img = img - float(it.sky_median)
        stamp = extract_stamp(img, it.xcen, it.ycen, radius, subpixel=subpixel)
        if align_north:
            # north_deg is measured CCW from +y; rotating by -north_deg puts
            # north back onto +y.
            stamp = rotate_nan_safe(stamp, -float(it.north_deg))
        if np.isfinite(stamp).any():
            stamps.append(stamp)
            labels.append(it.label)

    if not stamps:
        return None

    cube = np.array(stamps, dtype=np.float64)
    with warnings.catch_warnings():
        # Edge padding and bad pixels legitimately leave NaN in the cube;
        # sigma_clip masks them, which is what we want, but says so loudly.
        warnings.filterwarnings("ignore", message=".*invalid values.*")
        clipped = sigma_clip(cube, sigma=sigma, maxiters=maxiters, axis=0,
                             masked=True, copy=False)
    count = (~clipped.mask & np.isfinite(cube)).sum(axis=0).astype(np.int32)

    if combine not in ("median", "mean"):
        raise ValueError(f"unknown combine {combine!r}; use 'median' or 'mean'")

    with np.errstate(invalid="ignore", divide="ignore"):
        if combine == "median":
            image = np.ma.median(clipped, axis=0).filled(np.nan)
        else:
            image = np.ma.mean(clipped, axis=0).filled(np.nan)
        std = np.ma.std(clipped, axis=0, ddof=1).filled(np.nan)
        sem = std / np.sqrt(np.maximum(count, 1))
        # The median's sampling error is larger than the mean's by ~sqrt(pi/2)
        # for Gaussian noise; applying that factor keeps the reported ERROR
        # plane an honest uncertainty on the statistic actually stored.
        if combine == "median":
            sem = sem * np.sqrt(np.pi / 2.0)
        error = np.where(count > 1, sem, np.nan)

    return StackResult(image=np.asarray(image), count=count, error=np.asarray(error),
                       n_frames=len(stamps), radius=int(radius), labels=labels,
                       north_aligned=bool(align_north), combine=str(combine),
                       pix_scale_arcsec=float(pix_scale_arcsec))


def write_stack_fits(
    result: StackResult,
    path: Path,
    *,
    objdesig: str,
    epoch: int,
    band: str,
    wl_range: Tuple[float, float],
    flux_unit: str = "mJy/pixel",
    extra: Optional[Dict[str, Tuple]] = None,
) -> Path:
    """
    Write a stack as a three-extension FITS file.

    Extensions are ``IMAGE`` (the stack), ``COUNT`` (frames per pixel) and
    ``ERROR`` (sampling error of the stored statistic).  The primitive stacker wrote a
    single unnamed plane with no count map, no error and no WCS.

    A tangent-plane WCS centred on the comet is attached only when the stack was
    north-aligned; without that, the frames share no sky orientation and a WCS
    would be actively misleading.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    primary = fits.PrimaryHDU()
    h = primary.header
    h["OBJDESIG"] = (str(objdesig), "Target designation")
    h["EPOCH"] = (int(epoch), "Observing epoch (see epochs.group_epochs)")
    h["BAND"] = (str(band), "Custom SPHEREx spectral band")
    h["WL_MIN"] = (float(wl_range[0]), "[um] band lower edge")
    h["WL_MAX"] = (float(wl_range[1]), "[um] band upper edge")
    h["NCOMBINE"] = (int(result.n_frames), "Frames entering the stack")
    h["BUNIT"] = (flux_unit, "Unit of the IMAGE array")
    h["RAD_PIX"] = (int(result.radius), "Stamp half-width [pixel]")
    h["NORTHUP"] = (bool(result.north_aligned), "Frames rotated to north-up")
    h["COMBINE"] = (f"sigma-clipped {result.combine}", "Combination statistic")
    if np.isfinite(result.pix_scale_arcsec):
        h["PIXSCALE"] = (float(result.pix_scale_arcsec), "[arcsec/pixel]")
    for key, val in (extra or {}).items():
        h[key] = val

    img = fits.ImageHDU(np.asarray(result.image, dtype=np.float32), name="IMAGE")
    if result.north_aligned and np.isfinite(result.pix_scale_arcsec):
        c = float(result.radius)
        ih = img.header
        ih["CTYPE1"], ih["CTYPE2"] = "RA---TAN", "DEC--TAN"
        ih["CRPIX1"] = ih["CRPIX2"] = c + 1.0            # FITS is 1-indexed
        ih["CRVAL1"] = ih["CRVAL2"] = 0.0                # relative frame
        ih["CDELT1"] = -result.pix_scale_arcsec / 3600.0
        ih["CDELT2"] = result.pix_scale_arcsec / 3600.0
        ih["CUNIT1"] = ih["CUNIT2"] = "deg"
    else:
        img.header["COMMENT"] = ("Frames were not rotated to a common sky "
                                 "orientation; axes are detector pixels.")

    hdul = fits.HDUList([
        primary, img,
        fits.ImageHDU(np.asarray(result.count, dtype=np.int32), name="COUNT"),
        fits.ImageHDU(np.asarray(result.error, dtype=np.float32), name="ERROR"),
    ])
    hdul.writeto(path, overwrite=True)
    log.info("stack written: %s (%d frames)", path.name, result.n_frames)
    return path
