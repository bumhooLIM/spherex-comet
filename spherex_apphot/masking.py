"""
Bad-pixel identification.

The project draws a firm distinction between the two kinds of "masking" the
primitive code conflated:

**Flag masking (here).**  A flagged pixel is genuinely unusable.  It is removed
from the photometry entirely -- not zeroed, *excluded* -- so an aperture that
contains one sums the remaining pixels over a correspondingly reduced effective
area.  See :func:`phot.measure_exposure`.

**Source masking (not here).**  Gaia stars are *not* masked.  Painting discs
over them, as the primitive code did, deletes real coma flux and biases every
aperture that happens to sit near a star.  Contamination is instead recorded per
measurement by :mod:`sourceflag`, leaving the decision to the analyst.
"""

from __future__ import annotations

from typing import Iterable, Optional, Tuple

import numpy as np

from .logging_utils import get_logger

__all__ = ["flag_to_mask", "build_badpix_mask", "badpix_report"]

log = get_logger("masking")


def flag_to_mask(flag: np.ndarray, bad_bits: Iterable[int]) -> np.ndarray:
    """
    Convert a bitmask plane into a boolean bad-pixel mask.

    Parameters
    ----------
    flag : ndarray
        Integer flag plane.  Cast to ``uint32`` internally so that a negative
        ``int32`` value (bit 31 set) does not sign-extend during the shift.
    bad_bits : iterable of int
        Bit positions to treat as bad.  See :data:`config.FLAG_BIT_MEANINGS`.

    Returns
    -------
    ndarray of bool
        ``True`` where at least one of ``bad_bits`` is set.

    Notes
    -----
    The bit test is vectorised over the whole plane at once rather than looped
    per bit, which matters because this runs once per exposure over the full
    catalogue.
    """
    bits = [int(b) for b in bad_bits]
    if not bits:
        return np.zeros(flag.shape, dtype=bool)
    invalid = [b for b in bits if not 0 <= b <= 31]
    if invalid:
        raise ValueError(f"flag bits outside 0-31: {invalid}")

    flag_u = np.asarray(flag).astype(np.uint32, copy=False)
    selector = np.uint32(0)
    for b in bits:
        selector |= np.uint32(1) << np.uint32(b)
    return (flag_u & selector) != 0


def build_badpix_mask(
    sci: np.ndarray,
    var: np.ndarray,
    flag: np.ndarray,
    bad_bits: Iterable[int],
    *,
    mask_nonfinite_sci: bool = True,
    mask_bad_variance: bool = True,
) -> Tuple[np.ndarray, dict]:
    """
    Combine every reason a pixel cannot be used into one mask.

    Parameters
    ----------
    sci, var, flag : ndarray
        The three cutout planes.
    bad_bits : iterable of int
        Flag bits to reject.
    mask_nonfinite_sci : bool
        Reject non-finite science pixels.  Required in practice: in the sample
        data roughly 16 % of science pixels are NaN and *none* of them carries
        any flag bit, so relying on the flag plane alone leaks NaN into the
        aperture sums.
    mask_bad_variance : bool
        Reject non-finite or non-positive variance.  A zero variance would
        otherwise contribute zero uncertainty and inflate the reported SNR
        (review item R5).

    Returns
    -------
    mask : ndarray of bool
        ``True`` marks an unusable pixel.
    report : dict
        Pixel counts per rejection reason, for the log and for diagnostics.
    """
    sci = np.asarray(sci)
    m_flag = flag_to_mask(flag, bad_bits)
    m_sci = ~np.isfinite(sci) if mask_nonfinite_sci else np.zeros(sci.shape, bool)
    if mask_bad_variance:
        with np.errstate(invalid="ignore"):
            m_var = ~np.isfinite(var) | (np.asarray(var) <= 0)
    else:
        m_var = np.zeros(sci.shape, bool)

    mask = m_flag | m_sci | m_var
    report = {
        "n_flag": int(m_flag.sum()),
        "n_nonfinite_sci": int(m_sci.sum()),
        "n_bad_var": int(m_var.sum()),
        "n_badpix": int(mask.sum()),
        "frac_badpix": float(mask.mean()),
    }
    return mask, report


def badpix_report(mask: np.ndarray, name: Optional[str] = None) -> str:
    """One-line human summary of a bad-pixel mask."""
    n = int(mask.sum())
    return f"{name or 'cutout'}: {n} bad pixels ({100.0 * mask.mean():.1f} %)"
