"""Forward model: gas production rates -> monochromatic coma emission spectrum.

Ported from ``emission-fitter/gasmodel.py`` (v0.1.0); every equation was verified against a
re-derivation in ``doc/comspec/fitting_methodology.md`` section 9.  Since 2026-09-11 Layers 3-4 take the
fluorescence efficiencies and band shapes from the reconstructed GSFC-style database
(:mod:`fluorescence`, ``ModelParams.profile_source = "gfm"``) and scale g(CO) with the comet's
heliocentric velocity (Swings effect, ``ModelParams.co_swings``); the previous Gaussian bands
remain available as ``profile_source = "gaussian"``.  Implements Layers 0-4 of the concept design
(``doc/comspec/model_concept.md``):

===== ===================================================================================
Layer  Content
===== ===================================================================================
0      Geometry -- r_h/Delta scalings, r_h-scaled lifetimes, projected aperture
1      Coma density -- Haser profile per species
2      Aperture integration -- Yamamoto (1981) Bessel filling factor
3      Excitation -- g(r_h, v_h) = g(1 au) s(v_h) / r_h^2, optional pump-opacity correction
4      Monochromatic emission -- the species' g_lambda template times the photon energy
===== ===================================================================================

The single most useful structural fact, and the reason the fitter is a linear solve: for one
species every band shares the same ``N_ap``, the same ``r_h^-2`` and the same ``1/(4 pi Delta^2)``,
so the *relative* weights between its bands are fixed and only the amplitude depends on geometry.
Each species therefore has a geometry-independent spectral shape (:func:`species_shape_mjy`) times a
scalar amplitude (:func:`amplitude_per_unit_Q`).  The velocity dependence of g(CO) is a scalar too
(one band dominates), so it multiplies the amplitude and the fit stays linear.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.integrate import quad
from scipy.special import k0, k1, modstruve

from . import fluorescence as _fl
from .config import (BANDS, C_LIGHT, H_PLANCK, KM_M, AU_M, KAPPA_PUMP, RHO_TAU_REF_KM, SPECIES,
                     TAU_1AU, WM2UM_TO_MJY, ModelParams)

__all__ = [
    "filling_factor", "tau_at_rh", "column_per_unit_Q", "amplitude_per_unit_Q", "swings_factor",
    "band_profile", "species_shape_mjy", "hires_grid", "spectrum_mjy", "opacity_factor",
]


# ============================================================================ Layer 2: filling factor
def _ik0_struve(x):
    """int_x^inf K0(y) dy, vectorised.

    Closed form via modified Struve functions below x = 30; the asymptotic series above, where the
    closed form loses precision to cancellation.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    big = x > 30.0
    small = ~big
    xs = x[small]
    out[small] = 0.5 * np.pi - 0.5 * np.pi * xs * (
        k0(xs) * modstruve(-1, xs) + k1(xs) * modstruve(0, xs))
    xb = x[big]
    out[big] = np.sqrt(np.pi / (2 * xb)) * np.exp(-xb) * (1 - 5 / (8 * xb) + 129 / (128 * xb ** 2))
    return out


def _ik0_quad(x):
    """int_x^inf K0(y) dy by direct quadrature (scalar); the independent cross-check."""
    if x > 700.0:
        return 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        val, _ = quad(k0, x, np.inf, limit=200)
    return val


def filling_factor(x, method: str = "struve"):
    r"""Yamamoto (1981) parent-species filling factor :math:`f_1(x)` for a centred circular beam.

    With :math:`x = \rho_{ap} / (v_g \tau)`,

    .. math::
        f_1(x) = x\,g(x), \qquad g(x) = \frac1x - K_1(x) + \int_x^\infty K_0(y)\,dy

    Parameters
    ----------
    x : array_like
        Aperture radius in units of the photodissociation scale length.
    method : {"struve", "quad"}
        Route used for the Bessel integral. ``"quad"`` exists to cross-check ``"struve"``.

    Returns
    -------
    ndarray or float
        :math:`f_1(x)`, tending to :math:`\pi x/2` as :math:`x\to0` and to 1 as
        :math:`x\to\infty`.

    Notes
    -----
    Below ``x = 1e-4`` the term ``1/x - K1(x)`` is a catastrophic cancellation of two ~1/x
    quantities, so the series
    :math:`g(x) = \pi/2 + x[\tfrac12\ln(x/2) + \gamma/2 - \tfrac34] + O(x^3)` is used instead. Its
    leading term is the small-aperture limit and the O(x) correction is the first departure from a
    pure 1/rho column.
    """
    x = np.atleast_1d(np.asarray(x, dtype=float))
    g = np.empty_like(x)
    tiny = x < 1.0e-4
    if np.any(tiny):
        xt = x[tiny]
        g[tiny] = np.pi / 2 + xt * (0.5 * np.log(xt / 2) + 0.5 * np.euler_gamma - 0.75)
    if np.any(~tiny):
        xo = x[~tiny]
        ik0 = _ik0_struve(xo) if method == "struve" else np.array([_ik0_quad(v) for v in xo])
        g[~tiny] = 1.0 / xo - k1(xo) + ik0
    f1 = x * g
    return float(f1[0]) if f1.size == 1 else f1


# ================================================================== Layers 0-2: geometry -> column
def tau_at_rh(species: str, r_h_au):
    """Photodissociation lifetime [s] at heliocentric distance ``r_h_au`` (scales as r_h^2)."""
    return TAU_1AU[species] * np.asarray(r_h_au, dtype=float) ** 2


def column_per_unit_Q(species: str, r_h_au, params: ModelParams):
    """Molecules in the aperture per unit production rate, ``N_ap / Q = tau * f1(x)`` [s].

    Parameters
    ----------
    species : str
        One of ``config.SPECIES``.
    r_h_au : array_like
        Heliocentric distance of each measurement [au].
    params : ModelParams
        Supplies the aperture radius and the expansion-velocity law.
    """
    r_h = np.asarray(r_h_au, dtype=float)
    tau = tau_at_rh(species, r_h)
    v_g = params.v_g(r_h)
    x = params.rho_ap_km / (v_g * tau)
    return tau * np.atleast_1d(filling_factor(x))


def swings_factor(species: str, v_h_kms, params: ModelParams):
    """``g(v_h) / g(0)`` of the species (Layer 3): the CO table of the fluorescence database when
    ``params.co_swings`` is on, otherwise 1.  ``v_h_kms`` may be an array or None/NaN (-> 1)."""
    if not params.co_swings or v_h_kms is None or species not in _fl.SWINGS_SPECIES:
        return 1.0
    return _fl.swings_factor(species, v_h_kms, params.T_rot)


def amplitude_per_unit_Q(species: str, r_h_au, delta_au, params: ModelParams, v_h_kms=None):
    """Geometric amplitude multiplying a species' fixed spectral shape, per unit Q.

    Returns ``N_ap/Q * s(v_h) * r_h^-2 / (4 pi Delta^2)`` in SI-compatible units such that
    multiplying by :func:`species_shape_mjy` and by Q gives mJy.  ``s(v_h)`` is the Swings factor
    of :func:`swings_factor` (CO only); ``v_h_kms`` is the heliocentric radial velocity of each
    measurement, positive receding, and may be omitted (factor 1).
    """
    r_h = np.asarray(r_h_au, dtype=float)
    delta_m = np.asarray(delta_au, dtype=float) * AU_M
    return (column_per_unit_Q(species, r_h, params)
            * swings_factor(species, v_h_kms, params)
            * r_h ** -2.0
            / (4.0 * np.pi * delta_m ** 2))


# ================================================================================ Layer 3: opacity
def opacity_factor(species: str, Q: float, r_h_au: float, params: ModelParams) -> float:
    """Multiplicative correction to a species' emission for pump opacity.

    Mode 0 returns 1.0 exactly and is the only mode that leaves the model linear in Q, which is
    what the linear fitter requires. Modes 1 and 2 are provided so the systematic can be mapped
    after a fit, not used inside one.
    """
    if params.opacity_mode == 0 or Q <= 0:
        return 1.0

    r_h = float(r_h_au)
    tau = float(tau_at_rh(species, r_h))
    v_g = float(np.atleast_1d(params.v_g(r_h))[0])

    if params.opacity_mode == 1:                      # core excision
        rho_tau = RHO_TAU_REF_KM * (Q / 1.0e27)
        if rho_tau >= params.rho_ap_km:
            return 0.0
        x_ap = params.rho_ap_km / (v_g * tau)
        x_tau = rho_tau / (v_g * tau)
        f_ap = float(np.atleast_1d(filling_factor(x_ap))[0])
        f_tau = float(np.atleast_1d(filling_factor(x_tau))[0])
        return (f_ap - f_tau) / f_ap if f_ap > 0 else 0.0

    if params.opacity_mode == 2:                      # PSG first-order pump correction
        rr = np.linspace(1.0, params.rho_ap_km, 2000)
        # line-of-sight column [cm^-2] at impact parameter rr, thin-coma form with a decay factor
        n_per_km2 = Q / (4.0 * v_g * rr) * np.exp(-rr / (v_g * tau))
        n_col = n_per_km2 / 1.0e10
        kappa = np.mean([KAPPA_PUMP[b.key] for b in BANDS if b.species == species])
        tp = n_col * kappa
        corr = np.where(tp > 1e-8, (1.0 - np.exp(-tp)) / np.where(tp > 0, tp, 1.0), 1.0)
        return float(np.trapezoid(corr, rr) / (rr[-1] - rr[0]))

    raise ValueError(f"unknown opacity_mode {params.opacity_mode}")


# ================================================================= Layer 4: monochromatic emission
def band_profile(lam_um, band, T_rot: float = 70.0):
    """Unit-area Gaussian band profile :math:`\\Phi_b(\\lambda)` [um^-1] of the *legacy* model.

    Used only for ``ModelParams.profile_source = "gaussian"``; the main model takes the band
    shapes from the fluorescence database (:func:`species_shape_mjy`).  ``T_rot`` is accepted for
    signature compatibility and ignored.
    """
    sigma = band.fwhm_um / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    lam = np.asarray(lam_um, dtype=float)
    return np.exp(-0.5 * ((lam - band.lam_um) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def hires_grid(params: ModelParams) -> np.ndarray:
    """Log-spaced internal wavelength grid [um], uniform in resolving power."""
    n = int(np.log(params.lam_max_um / params.lam_min_um) * params.resolving_power) + 1
    return params.lam_min_um * np.exp(np.arange(n) / params.resolving_power)


def species_shape_mjy(species: str, lam_um, params: ModelParams) -> np.ndarray:
    """Geometry-independent spectral shape of one species, in mJy per unit amplitude.

    ``profile_source = "gfm"``: the database template ``g_lambda(lambda; T_rot)`` [photons s^-1
    molecule^-1 um^-1] rebinned onto ``lam_um`` (flux-conserving) times the photon energy --
    every band the species emits between 0.7 and 5.0 um, with its strength and shape, at 1 au and
    v_h = 0.  ``"gaussian"``: the eight bands of :data:`config.BANDS`, each ``g_b(1 au) *
    hc/lambda_b`` spread over a unit-area Gaussian.  Either way F_lambda is converted to F_nu so
    the result compares directly with a flux-density measurement.

    Multiply by :func:`amplitude_per_unit_Q` and by Q to get mJy.
    """
    lam = np.asarray(lam_um, dtype=float)
    if params.profile_source == "gfm":
        g_lam = _fl.g_lambda(species, lam, params.T_rot)
        e_photon = H_PLANCK * C_LIGHT / (lam * 1e-6)           # J per emitted photon, per cell
        out = g_lam * e_photon
    elif params.profile_source == "gaussian":
        out = np.zeros_like(lam)
        for b in BANDS:
            if b.species != species:
                continue
            e_photon = H_PLANCK * C_LIGHT / (b.lam_um * 1e-6)
            out += b.g1au * e_photon * band_profile(lam, b, params.T_rot)
    else:
        raise ValueError(f"unknown profile_source {params.profile_source!r} (gfm | gaussian)")
    return out * lam ** 2 * WM2UM_TO_MJY


def spectrum_mjy(lam_um, Q: dict, r_h_au: float, delta_au: float,
                 params: ModelParams, v_h_kms=None) -> dict:
    """Monochromatic gas spectrum [mJy] at a single geometry.

    Parameters
    ----------
    lam_um : array_like
        Wavelengths [um].
    Q : dict
        Production rates keyed by species [molecules s^-1]. Missing species are treated as zero.
    r_h_au, delta_au : float
        Heliocentric and observer distance of the measurement.
    params : ModelParams
    v_h_kms : float, optional
        Heliocentric radial velocity [km/s], positive receding; sets the CO Swings factor.

    Returns
    -------
    dict
        ``{"lam": ..., "total": ..., "H2O": ..., "CO2": ..., "CO": ...}`` -- the total and the
        per-species decomposition, all in mJy.
    """
    lam = np.asarray(lam_um, dtype=float)
    out = {"lam": lam, "total": np.zeros_like(lam)}
    for s in SPECIES:
        q = float(Q.get(s, 0.0) or 0.0)
        amp = float(np.atleast_1d(amplitude_per_unit_Q(s, r_h_au, delta_au, params, v_h_kms))[0])
        comp = q * amp * species_shape_mjy(s, lam, params) * opacity_factor(s, q, r_h_au, params)
        out[s] = comp
        out["total"] = out["total"] + comp
    return out
