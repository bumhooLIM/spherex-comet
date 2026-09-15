"""Layer 5: the SPHEREx instrument response (ported unchanged from ``emission-fitter``).

SPHEREx disperses with linear variable filters, so every measurement in ``data/apphot`` is a flux
density through one filter channel rather than a monochromatic sample. The photometry table records
that channel directly: ``wl`` is its centre and ``wlwidth`` its FWHM.

This module turns those two numbers into a normalised Gaussian throughput and band-averages the
model through it, which is exactly what the detector does::

    F_obs(i) = int F_nu(lam) T_i(lam) dlam / int T_i(lam) dlam

Doing the averaging in F_nu (mJy) rather than F_lambda avoids having to define an effective
wavelength for the conversion afterwards.

A Gaussian is a stand-in for the true SPHEREx line spread function, whose measured shape is not yet
in this repository -- see ``config.PLACEHOLDERS`` and the instrument-model gap noted in the concept
document. It is the right *width* by construction, since ``wlwidth`` is the tabulated FWHM.
"""

from __future__ import annotations

import numpy as np

__all__ = ["FWHM_TO_SIGMA", "gaussian_bandpass", "bandpass_average", "bandpass_matrix"]

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def gaussian_bandpass(lam_um, wl_c: float, wlwidth: float) -> np.ndarray:
    """Normalised Gaussian throughput centred on ``wl_c`` with FWHM ``wlwidth``.

    Returns a profile that integrates to 1 over ``lam_um``; the normalisation is analytic rather
    than numerical so a truncated grid does not bias the average.
    """
    sigma = wlwidth * FWHM_TO_SIGMA
    lam = np.asarray(lam_um, dtype=float)
    return np.exp(-0.5 * ((lam - wl_c) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def bandpass_matrix(lam_hi, wl_c, wlwidth, n_sigma: float = 5.0) -> np.ndarray:
    """Throughput-weighting matrix ``W`` with ``W @ F_hi`` giving the band-averaged flux.

    Parameters
    ----------
    lam_hi : array_like, shape (n_lam,)
        High-resolution wavelength grid [um]; must be sorted ascending.
    wl_c, wlwidth : array_like, shape (n_points,)
        Channel centres and FWHMs [um].
    n_sigma : float
        Truncation radius. Contributions beyond this are set to zero before renormalising, which
        keeps the matrix sparse-ish without biasing the result (5 sigma omits ~6e-7 of the weight).

    Returns
    -------
    ndarray, shape (n_points, n_lam)
        Each row integrates to 1 against ``lam_hi`` under the trapezoid rule, so the matrix product
        is a proper weighted mean rather than an unnormalised integral.
    """
    lam = np.asarray(lam_hi, dtype=float)
    wl_c = np.atleast_1d(np.asarray(wl_c, dtype=float))
    wlw = np.atleast_1d(np.asarray(wlwidth, dtype=float))
    sigma = wlw[:, None] * FWHM_TO_SIGMA

    d = (lam[None, :] - wl_c[:, None]) / sigma
    W = np.exp(-0.5 * d ** 2)
    W[np.abs(d) > n_sigma] = 0.0

    # trapezoid weights of the grid itself, so the row normalisation matches the integration rule
    dl = np.gradient(lam)
    W = W * dl[None, :]
    norm = W.sum(axis=1, keepdims=True)
    bad = norm[:, 0] <= 0
    norm[bad] = 1.0                     # channel entirely off the grid -> row of zeros, not NaN
    W = W / norm
    W[bad] = 0.0
    return W


def bandpass_average(lam_hi, flux_hi, wl_c, wlwidth, n_sigma: float = 5.0) -> np.ndarray:
    """Band-average ``flux_hi`` through Gaussian channels at ``(wl_c, wlwidth)``.

    Convenience wrapper over :func:`bandpass_matrix` for a single spectrum.
    """
    W = bandpass_matrix(lam_hi, wl_c, wlwidth, n_sigma=n_sigma)
    return W @ np.asarray(flux_hi, dtype=float)
