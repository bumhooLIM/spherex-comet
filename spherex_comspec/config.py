"""
Every constant, window, threshold and placeholder of the ``spherex_comspec`` pipeline.

Nothing here computes anything.  The module exists so that a result can be
traced to the exact parameter set that produced it (``Variant.hash``), and so
that every value that is a *stand-in* rather than a measurement is declared as
such in one place (:data:`PLACEHOLDERS`).

Scientific lineage
------------------
The emission model is the AKARI/IRC method of Ootsubo et al. (2012, ApJ 752,
15): band-integrated fluorescence of H2O, CO2 and CO converted to production
rates with the Crovisier & Encrenaz (1983) band model, on a Haser coma with the
Yamamoto (1981) aperture filling factor.  The continuum subtraction and the
phase regrouping reproduce the methodology of the project notebooks
``continuum_subtraction.ipynb`` and ``phase_group_update.ipynb``.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple

__all__ = [
    "H_PLANCK", "C_LIGHT", "AU_M", "KM_M", "WM2UM_TO_MJY", "MJY_TO_WM2UM", "C_UM_S",
    "SPECIES", "TAU_1AU", "Band", "BANDS", "BAND_WINDOWS", "EMISSION_WINDOWS", "ALL_EM_WINDOWS",
    "KEY_RANGES", "H2O_HOT_RANGE", "BAND_CARRIER", "BAND_COLORS", "EMISSION_DTYPES", "SOURCEFLAG_PRIORITY",
    "RHO_TAU_REF_KM", "KAPPA_PUMP", "PLACEHOLDERS",
    "GroupingConfig", "ApertureConfig", "FlagPolicy", "ContinuumConfig", "ModelParams",
    "FitConfig", "Variant", "BASELINE_FLAGS", "DEFAULT_VARIANTS", "MAIN_VARIANT", "VARIANTS",
    "config_hash",
]

# ---------------------------------------------------------------------- physical constants
H_PLANCK = 6.62607015e-34          # J s
C_LIGHT = 2.99792458e8             # m / s
C_UM_S = 2.99792458e14             # um / s
AU_M = 1.495978707e11              # m
KM_M = 1.0e3

#: F_nu[mJy] = F_lambda[W m^-2 um^-1] * lam[um]^2 * this
WM2UM_TO_MJY = 1.0e23 / C_LIGHT
#: F_lambda[W m^-2 um^-1] = F_nu[mJy] * this / lam[um]^2
MJY_TO_WM2UM = 1.0e-29 * C_UM_S

#: dtypes to pin on every emission-CSV read: "2022E2" and "2024E1" are valid scientific
#: notation and pandas would otherwise turn them into floats.
EMISSION_DTYPES = {"target": str, "band": str, "role": str, "verdict": str,
                   "fail_reason": str, "notes": str, "sourceflag": str}

SPECIES = ("H2O", "CO2", "CO")

#: photodissociation lifetimes at 1 au [s]; scale as r_h^2 (Ootsubo et al. 2012, Table 2)
TAU_1AU = {"H2O": 8.3e4, "CO2": 5.0e5, "CO": 1.3e6}

#: Source-flag precedence written by ``spherex_apphot``: a > b > c > d > 0.
SOURCEFLAG_PRIORITY = ("a", "b", "c", "d", "0")


# ------------------------------------------------------------------------------- band table
#: The main model (``ModelParams.profile_source = "gfm"``) takes every band's g-factor *and*
#: shape from the reconstructed fluorescence database (``data/fluorescence``, see
#: ``doc/fluorescence_database.md``).  This table is the previous representation -- eight
#: Gaussian bands with the Ootsubo et al. (2012, Table 2) g-factors -- kept for the
#: ``profile_source = "gaussian"`` study variant and for ``KAPPA_PUMP``.  For reference, the
#: database gives at 70 K and 1 au (photons s^-1 molecule^-1): H2O nu3 3.14e-4, nu1 2.9e-5,
#: nu2+nu3-nu2 2.85e-5, nu1+nu3-nu1 1.9e-5, nu3-nu2 6.6e-6, nu1-nu2 4.4e-6 (plus ~1e-5 in the
#: 2.8-3.0 um hot bands); CO2 nu3 2.71e-3 (+1.0e-4 in hot bands); CO v(1-0) 1.92e-4 at v_h = 0
#: rising to 2.4-2.5e-4 for |v_h| > 10 km/s (the Swings effect, ``co_swings.csv``).
@dataclass(frozen=True)
class Band:
    """One vibrational emission band (Ootsubo et al. 2012, Table 2) of the *legacy* Gaussian model."""

    species: str
    label: str
    lam_um: float
    g1au: float          #: fluorescence efficiency at 1 au [photons s^-1 molecule^-1]
    fwhm_um: float       #: intrinsic envelope FWHM [um] -- PLACEHOLDER
    hot: bool = False

    @property
    def key(self) -> str:
        return f"{self.species}:{self.label}"


BANDS: Tuple[Band, ...] = (
    Band("H2O", "nu3",           2.66, 2.82e-4, 0.060),
    Band("H2O", "nu1",           2.73, 2.47e-5, 0.050),
    Band("H2O", "nu2+nu3-nu2",   2.66, 2.81e-5, 0.070, hot=True),
    Band("H2O", "nu1+nu3-nu1",   2.73, 2.19e-5, 0.060, hot=True),
    Band("H2O", "nu3-nu2",       4.63, 7.66e-6, 0.100, hot=True),
    Band("H2O", "nu1-nu2",       4.85, 7.35e-6, 0.100, hot=True),
    Band("CO2", "nu3",           4.26, 2.9e-3,  0.020),
    Band("CO",  "v(1-0)",        4.67, 2.6e-4,  0.045),
)

#: Emission and continuum windows per band [um], as in ``continuum_subtraction.ipynb``.
#: ``em`` is punched out of *every* continuum fit; ``cont`` bounds the sample for this band.
#: 2026-09-12: the 2.7 um emission window starts at 2.50 um (the H2O complex begins at 2.55 um;
#: 2.60 um lost 6.7 % of it) and the H2O and CO2 continuum windows are 0.1 um wider on both
#: sides; the 4.3 um continuum reaches into the CO emission window, which is punched out anyway.
BAND_WINDOWS: Dict[str, dict] = {
    "2.7um": dict(em=(2.50, 2.80), cont=(2.20, 3.10), lam_c=2.70),
    "4.3um": dict(em=(4.18, 4.35), cont=(3.90, 4.65), lam_c=4.26),
    "4.7um": dict(em=(4.55, 4.90), cont=(4.40, 5.00), lam_c=4.67),
}
EMISSION_WINDOWS = {b: w["em"] for b, w in BAND_WINDOWS.items()}
ALL_EM_WINDOWS = [w["em"] for w in BAND_WINDOWS.values()]

#: Diagnostic ranges that must be sampled before a species may be fitted at all.
#: Falling short is a *coverage* failure, categorically different from a non-detection.
KEY_RANGES = {
    "H2O": dict(lo=2.50, hi=2.80, min_points=3),
    "CO2": dict(lo=4.20, hi=4.30, min_points=2),
    "CO":  dict(lo=4.60, hi=4.70, min_points=2),
}
#: Fallback coverage for H2O when the 2.7 um main band is not covered (saturated or flagged
#: channels of bright, close comets such as 10P and 24P): the nu3-nu2 (4.63 um) and nu1-nu2
#: (4.85 um) hot bands.  ``min_points`` channels inside ``lo``-``hi`` and ``min_points_red``
#: of them beyond ``red_lo`` -- the 4.85 um band lies outside CO v(1-0), and that is what
#: keeps a hot-band Q(H2O) separable from Q(CO).  The hot-band g-factors come from the
#: reconstructed database (line-level validation to 10 %; the 4.85 um band is 40 % weaker than
#: the earlier harmonic estimate), but such a value rests on ~3 % of the water emission and
#: carries ``h2o_source = "hot"``.
#: ``max_rh_au`` (2026-09-12): beyond 3 au the hot-band values were all 1-2 sigma with errors
#: comparable to the value, so the fallback is not offered there (``FitConfig.h2o_hot_max_rh_au``).
H2O_HOT_RANGE = dict(lo=4.55, hi=4.90, min_points=3, red_lo=4.75, min_points_red=1, max_rh_au=3.0)

BAND_CARRIER = {
    "2.7um": r"H$_2$O $\nu_3+\nu_1$+hot",
    "4.3um": r"CO$_2$ $\nu_3$",
    "4.7um": r"CO v(1-0) $\oplus$ H$_2$O hot",
}
BAND_COLORS = {"2.7um": "tab:blue", "4.3um": "tab:orange", "4.7um": "tab:green"}

#: Opaque-core radius at Q = 1e27 s^-1 [km]; opacity mode 1 only.  PLACEHOLDER.
RHO_TAU_REF_KM = 1.0e3
#: S_p / v_p [cm^2] for the CO2 nu3 pump, mode 2 only; other bands scale with g.  PLACEHOLDER.
_KAPPA_CO2 = 3.0e-14
KAPPA_PUMP = {b.key: _KAPPA_CO2 * (b.g1au / 2.9e-3) for b in BANDS}


# ------------------------------------------------------------------------------ dataclasses
@dataclass(frozen=True)
class GroupingConfig:
    """
    Rules that subdivide a 28-day ``epoch`` into ``phase`` groups.

    Reproduces ``notebooks/phase_group_update.ipynb``.  In priority order:

    1. **r_h range (hard).**  ``(max - min) / mean < rh_tol`` inside every group;
       and, since 2026-09-12, ``(max - min) / mean < delta_tol`` for the observer
       distance (rule 1b), which also subdivides the manual r_h bins.
    2. **Perihelion (hard).**  Inbound and outbound epochs never share a group;
       applied only where the turn is resolved by more than ``arc_min_drh`` on
       both sides, so ephemeris round-off cannot trigger it.
    3. **Bands first.**  Epochs sampling one emission window that are mutually
       within tolerance stay together; an extra group is accepted to keep a
       band whole.
    4. **Fewest groups**, with cuts placed in real gaps: a cut between epochs
       closer than ``link_drh`` is penalised.
    """

    rh_tol: float = 0.10
    #: rule 1b (2026-09-12): ``(max - min) / mean`` of the observer distance inside every
    #: automatic group, and inside every manual r_h bin; None disables it.  Delta is not
    #: what the model integrates, but a group is one observing state only if it is bounded
    #: too (r_h^2 Delta^2 spread reached 183 % before).  0.20 splits 20 of 174 groups.
    delta_tol: Optional[float] = 0.20
    arc_min_drh: float = 0.001          # au
    arc_min_epoch: int = 2
    link_drh: float = 0.05              # au
    #: targets exempt from rule 2 (perihelion resolved but marginal and far out)
    manual_no_arc: Tuple[str, ...] = ("2023R1", "2023V1")
    #: hand-set descending r_h bin edges per arc, replacing rules 1-4 in those blocks
    manual_edges: Dict[str, Dict[str, List[float]]] = field(default_factory=lambda: {
        "24P": {"in": [1.85, 1.40, 1.20], "out": []},
        "2024E1": {"in": [3.40]},
    })


@dataclass(frozen=True)
class ApertureConfig:
    """
    One aperture per target, from its mean heliocentric distance.

    ``near_km`` inside ``rh_split_au``, ``far_km`` beyond -- the rule of
    ``continuum_subtraction.ipynb``.  The revised photometry additionally
    *refuses* an aperture that is smaller than the PSF or larger than the sky
    annulus, so a rule aperture can be absent for part or all of a distant
    target's exposures (2014 UN271 at 14.6 au has no 40 000 km measurement at
    all).  ``min_coverage`` therefore promotes the target to the smallest
    ``km`` aperture that exists for at least that fraction of its exposures.
    """

    rh_split_au: float = 3.0
    near_km: float = 20000.0
    far_km: float = 40000.0
    min_coverage: float = 0.95
    kind: str = "km"
    #: ``"snr"`` (2026-09-12): among the ``km`` apertures present for ``min_coverage`` of the
    #: exposures, at least ``min_psf_mult`` PSF FWHM wide and no larger than
    #: ``max_ap_frac_annulus`` of the sky annulus' inner radius (so the annulus stays outside
    #: the coma the aperture measures), take the smallest one whose median S/N over the
    #: star-free emission-window channels is within ``snr_tol`` of the best.  ``"rh"``: the previous
    #: rule, ``near_km`` inside ``rh_split_au`` and ``far_km`` beyond, promoted for coverage.
    rule: str = "snr"
    snr_tol: float = 0.10
    min_psf_mult: float = 2.0
    max_ap_frac_annulus: float = 1.0 / 3.0
    #: star-free (flag 0) emission-window channels an aperture needs before its S/N is trusted;
    #: with fewer at every aperture the smallest bounded aperture is taken (crowded field)
    min_clean_channels: int = 5


@dataclass(frozen=True)
class FlagPolicy:
    """
    Which photometry rows enter a spectrum.

    ``drop_badphot`` removes rows with a bad pixel inside the aperture -- the
    revised pipeline's definition, *not* the flux-sign cut the old ``badphot``
    encoded.  ``max_frac_badpix`` (when set) relaxes that to rows whose bad
    fraction exceeds the threshold, which recovers most of the sample for
    bright, close comets where any 50-pixel aperture contains one bad pixel.
    ``drop_flags`` lists source flags to exclude (``a`` bright blend, ``b``
    flux-ratio blend, ``c`` any source, ``d`` SNR < 1).
    """

    name: str = "all"
    drop_flags: Tuple[str, ...] = ()
    drop_badphot: bool = True
    max_frac_badpix: Optional[float] = None

    def describe(self) -> str:
        bp = ("badphot dropped" if self.drop_badphot and self.max_frac_badpix is None
              else f"frac_badpix > {self.max_frac_badpix} dropped" if self.drop_badphot
              else "badphot kept")
        fl = f"flags {'+'.join(self.drop_flags)} dropped" if self.drop_flags else "all flags kept"
        return f"{self.name}: {bp}, {fl}"


@dataclass(frozen=True)
class ContinuumConfig:
    """Local polynomial continuum subtraction (``continuum_subtraction.ipynb``)."""

    #: maximum (``order_mode="cv"``) or requested (``"fixed"``) polynomial order per band
    poly_orders: Dict[str, int] = field(default_factory=lambda: {"2.7um": 3, "4.3um": 3, "4.7um": 3})
    #: ``"cv"`` (2026-09-12): the order of every bracketed fit is the lowest one whose
    #: cross-validated RMSE is within ``cv_select_margin`` of the best among 1..max; the
    #: previous fixed 3 / 2 / 2 orders were contradicted by the CV in half of the fits.  The
    #: 10 % margin is the parsimony rule: a higher order must buy a real reduction of the
    #: leave-one-out error, not the ~5 % that noise alone produces with ~20 points.
    order_mode: str = "cv"
    cv_select_margin: float = 0.10
    #: when the continuum window has points on one side of the band only, extend the empty
    #: side to this distance from the emission edge (still punching out every emission
    #: window) before fitting; None keeps the window as configured (2026-09-12).
    one_sided_extend_um: Optional[float] = 1.0
    sigma: float = 3.0
    maxiters: int = 5
    n_min_poly: int = 6              # below this many continuum points -> 2-point fallback
    z_shape_max: float = 3.0
    pos_nsig_max: float = 1.0
    cv_tol: float = 1.25
    cv_loo_max: int = 50
    #: sufficiency gates on the group as a whole; a failing group is skipped, not FAILed
    min_points_phase: int = 5
    min_points_emission: int = 2
    plot_margin_um: float = 0.05
    max_order: int = 3


@dataclass(frozen=True)
class ModelParams:
    """Fixed parameters of the coma model (the three Q are what the fit solves for)."""

    rho_ap_km: float = 40000.0
    v_g_kms: Optional[float] = None      # None -> 0.8 r_h^-0.5
    T_rot: float = 70.0
    opacity_mode: int = 0
    lam_min_um: float = 2.0
    lam_max_um: float = 5.2
    resolving_power: int = 4000
    #: ``"gfm"``: species templates (g-factors and band shapes at ``T_rot``) from the reconstructed
    #: fluorescence database ``data/fluorescence`` (:mod:`fluorescence`); ``"gaussian"``: the
    #: previous eight Gaussian bands of :data:`BANDS`.
    profile_source: str = "gfm"
    #: scale g(CO) with the comet's heliocentric velocity (Swings effect; ``co_swings.csv``)
    co_swings: bool = True
    #: label of the database build the templates come from, so it is part of the variant hash
    fluorescence_db: str = "gfm-2026-09-11"

    def v_g(self, r_h_au):
        import numpy as np
        if self.v_g_kms is not None:
            return np.full_like(np.asarray(r_h_au, dtype=float), float(self.v_g_kms))
        return 0.8 * np.asarray(r_h_au, dtype=float) ** -0.5


@dataclass(frozen=True)
class FitConfig:
    """Options of the weighted linear least-squares solve."""

    scale_errors_by_chi2: bool = True
    #: detection tiers (2026-09-12): ``detected`` at >= detection_sigma, ``marginal`` between
    #: marginal_sigma and detection_sigma (value reported, limit quoted), ``upper_limit`` below;
    #: limits are ``Q_fit + upper_limit_sigma * err``
    detection_sigma: float = 3.0
    marginal_sigma: float = 1.0
    upper_limit_sigma: float = 3.0
    #: one-sided cut of channels below -k sigma; None (2026-09-12) keeps every channel, since the
    #: cut biased Q upward by ~0.3 sigma per channel and could empty a group
    drop_negative_sigma: Optional[float] = None
    require_key_coverage: bool = True
    #: use the 4.6-4.9 um hot bands for Q(H2O) when, and only when, 2.7 um is not covered
    h2o_hot_fallback: bool = True
    #: ... and only inside this heliocentric distance (None: no cap)
    h2o_hot_max_rh_au: Optional[float] = 3.0
    #: generalised least squares (2026-09-12): the continuum-model uncertainty is correlated
    #: across the channels of a band through the polynomial coefficients; their saved
    #: covariance builds the full data covariance.  False: diagonal errors as before.
    gls: bool = True
    bands: Optional[Tuple[str, ...]] = None
    accept_verdicts: Tuple[str, ...] = ("PASS", "WARN")
    clip_sigma: Optional[float] = None
    species: Tuple[str, ...] = SPECIES


@dataclass(frozen=True)
class Variant:
    """
    One complete pipeline run: a flag policy and a flux space.

    ``use_distcorr`` selects ``flux_distcorr_mjy`` (flux x r_h^2 x Delta^2, the
    flux the comet would show at 1 au / 1 au) for the continuum subtraction.
    The production rates are always retrieved in *physical* flux space -- the
    emission is divided back by each channel's own ``distcorr_factor`` before
    the design matrix is built -- so Q is the rate at the comet whichever space
    the continuum was fitted in.  See ``fitting.fit_production_rates``.
    """

    name: str
    flags: FlagPolicy = field(default_factory=FlagPolicy)
    use_distcorr: bool = True
    #: ``source_sum_err_empirical_mjy`` (placeholder 10, applied 2026-09-09): the formal error
    #: under-reports the annulus scatter (sky_excess_ratio ~ 1.1) and the fits ran at chi2_nu ~ 2.4
    error_column: str = "source_sum_err_empirical_mjy"     # or "source_sum_err_mjy"
    grouping: GroupingConfig = field(default_factory=GroupingConfig)
    aperture: ApertureConfig = field(default_factory=ApertureConfig)
    continuum: ContinuumConfig = field(default_factory=ContinuumConfig)
    fit: FitConfig = field(default_factory=FitConfig)
    model: ModelParams = field(default_factory=ModelParams)
    #: What the run is for.  ``main`` is the catalog result; ``flags``, ``distcorr``,
    #: ``badphot``, ``fluorescence``, ``errors`` and ``rules`` are its study partners;
    #: ``previous`` is an earlier baseline kept for a before/after comparison.  The driver selects study partners by
    #: role, so a new variant needs no change anywhere else.
    role: str = "study"

    @property
    def flux_column(self) -> str:
        return "flux_distcorr_mjy" if self.use_distcorr else "source_sum_mjy"

    @property
    def err_column(self) -> str:
        if not self.use_distcorr:
            return self.error_column
        return {"source_sum_err_mjy": "flux_distcorr_err_mjy",
                "source_sum_err_empirical_mjy": "flux_distcorr_err_empirical_mjy"}[self.error_column]

    def to_dict(self) -> dict:
        return json.loads(json.dumps(asdict(self), default=list))

    @property
    def hash(self) -> str:
        # the role is a label for the driver, not a parameter of the run
        return config_hash({k: v for k, v in self.to_dict().items() if k != "role"})


def config_hash(d: dict) -> str:
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------- the variants
#: The source-flag policy of the main result.  Placeholder 5 was applied on 2026-09-09: a row is
#: dropped only when more than 5 % of its aperture is bad (the strict rule discarded 52 % of
#: 24P at no gain), and flag ``a`` (a G < 13 star within r_ap + 2 FWHM) is dropped -- it costs
#: 89 of 3 815 channels and nothing in Q.
BASELINE_FLAGS = FlagPolicy("main", drop_flags=("a",), max_frac_badpix=0.05)

#: The runs the analysis compares.  ``dc_*`` fit the continuum in distance-corrected space;
#: ``raw_main`` is the baseline policy in physical flux space (distance-correction study);
#: ``dc_all`` is the 2026-09-08 baseline -- strict ``badphot``, every flag kept -- so the effect
#: of the placeholder switch is itself a product of the run.
DEFAULT_VARIANTS: Tuple[Variant, ...] = (
    Variant("dc_main",        BASELINE_FLAGS, use_distcorr=True, role="main"),
    Variant("dc_main_keep_a", FlagPolicy("keep_a", max_frac_badpix=0.05),
            use_distcorr=True, role="flags"),
    Variant("dc_main_no_b",   FlagPolicy("no_ab", drop_flags=("a", "b"), max_frac_badpix=0.05),
            use_distcorr=True, role="flags"),
    Variant("raw_main",       BASELINE_FLAGS, use_distcorr=False, role="distcorr"),
    Variant("dc_main_strict", FlagPolicy("strict_no_a", drop_flags=("a",)),
            use_distcorr=True, role="badphot"),
    Variant("dc_all",         FlagPolicy("all"), use_distcorr=True, role="previous"),
    # the emission model before 2026-09-11: Gaussian bands with the Ootsubo g-factors and a
    # velocity-independent g(CO); paired with the main run it isolates what the reconstructed
    # fluorescence database changes
    Variant("dc_main_gauss",  BASELINE_FLAGS, use_distcorr=True,
            model=ModelParams(profile_source="gaussian", co_swings=False), role="fluorescence"),
    # 2026-09-12: the fit with diagonal errors (no continuum covariance) -- what GLS changes
    Variant("dc_main_diag",   BASELINE_FLAGS, use_distcorr=True, fit=FitConfig(gls=False), role="errors"),
    # 2026-09-12: every rule of the 2026-09-11 run that a variant can carry -- r_h-based aperture,
    # fixed 3/2/2 orders without the one-sided extension, the 1 sigma negative cut and detection
    # tier, diagonal errors, no hot-band distance cap.  Windows and grouping are shared.
    Variant("dc_rules_previous", BASELINE_FLAGS, use_distcorr=True,
            aperture=ApertureConfig(rule="rh"),
            continuum=ContinuumConfig(poly_orders={"2.7um": 3, "4.3um": 2, "4.7um": 2},
                                      order_mode="fixed", one_sided_extend_um=None),
            fit=FitConfig(detection_sigma=1.0, upper_limit_sigma=1.0, drop_negative_sigma=1.0,
                          h2o_hot_max_rh_au=None, gls=False),
            role="rules"),
)
#: The variant whose products are the catalog's main result.
MAIN_VARIANT = "dc_main"
#: Name -> variant registry: the CLI and any script that drives one variant select it here,
#: so a variant added to ``DEFAULT_VARIANTS`` becomes addressable by name everywhere at once.
VARIANTS: Dict[str, Variant] = {v.name: v for v in DEFAULT_VARIANTS}


# --------------------------------------------------------------------- placeholder registry
#: Every value that is a stand-in, a convention, or an unattributed choice rather than a
#: measurement.  ``priority`` is the order in which they should be replaced.
PLACEHOLDERS: Tuple[dict, ...] = (
    dict(priority=1, quantity="SPHEREx line-spread function",
         value="Gaussian of FWHM = wlwidth (instrument.gaussian_bandpass)",
         role="channel bandpass used to band-average the model; the CO / H2O-hot-band "
              "separability at 4.7 um depends on its wings, not only its width",
         update="as-built LVF R(lambda) and LSF from the SPHEREx instrument model"),
    dict(priority=2, quantity="fluorescence g-factors and band profiles",
         value="reconstructed GSFC-style database at T_rot (data/fluorescence; applied 2026-09-11; "
               "was: 8 Gaussian bands with the Ootsubo et al. 2012 g-factors)",
         role="every band's strength and shape, including the 4.6-4.9 um H2O hot bands under CO; "
              "validated line by line to ~10 % against the published GSFC values",
         update="applied; residual: HITRAN hot-band completeness and the T_rot dependence of "
                "the hot-band ratios (doc/fluorescence_database.md section 5)"),
    dict(priority=2, quantity="g(CO) versus heliocentric velocity (Swings effect)",
         value="ratio g(v_h)/g(0) from data/fluorescence/co_swings.csv interpolated in T_rot, with "
               "v_h = d r_h/dt from the ephemeris r_hel(t) of each target's exposures "
               "(dataio.heliocentric_velocity)",
         role="g(CO) is 25 % lower at v_h = 0 than for |v_h| > 10 km/s; a constant g(CO) under-"
              "estimated Q(CO) near perihelion by that much",
         update="JPL Horizons radial rates per exposure if the finite-difference velocity is ever "
                "in doubt; the ratio depends on T_rot by up to 4 % at 2-10 km/s"),
    dict(priority=3, quantity="expansion-velocity law",
         value="v_g = 0.8 r_h^-0.5 km/s (ModelParams.v_g)",
         role="fixed, not fitted (F ~ Q/v_g is exactly degenerate); the largest single "
              "systematic in absolute Q, ~ +/-18 % per 20 % in v_g for CO2",
         update="species-specific velocity from a literature compilation, with its citation"),
    dict(priority=4, quantity="2.7 um emission window and continuum windows",
         value="emission 2.50-2.80 um, continuum 2.20-3.10 (H2O) and 3.90-4.65 um (CO2); an empty "
               "continuum side is extended to 1 um from the band edge (applied 2026-09-12; was: "
               "2.60 um edge, 2.30-3.00 and 4.00-4.55 um, no extension)",
         role="the H2O template starts at 2.55 um; wider windows and the extension cut the "
              "one-sided FAIL rate",
         update="applied; re-examine the extension's linear extrapolation over up to 1 um"),
    dict(priority=5, quantity="badphot policy",
         value="drop rows with frac_badpix_ap > 0.05 (FlagPolicy.max_frac_badpix; applied 2026-09-09, "
               "was: any aperture containing a bad pixel)",
         role="the revised photometry flags *any* bad pixel; for close, bright comets this removes "
              "~50 % of rows.  frac_badpix_ap is in the data for a threshold instead",
         update="applied; re-examine the 0.05 threshold once the flux loss per bad pixel is characterised"),
    dict(priority=6, quantity="aperture rule",
         value="S/N-driven per target: the smallest km aperture (>= 95 % coverage, >= 2 PSF FWHM, "
               "<= 1/3 of the sky-annulus inner radius) within 10 % of the best median emission-window "
               "S/N (ApertureConfig.rule = 'snr'; applied 2026-09-12; was: 20 000 / 40 000 km by r_h)",
         role="sets the coma column the model integrates; the annulus bound keeps the background "
              "outside the coma being measured",
         update="applied; the dc_rules_previous variant keeps the r_h rule"),
    dict(priority=7, quantity="continuum polynomial orders",
         value="chosen per fit by leave-one-out cross-validation up to order 3, lowest order within "
               "5 % of the best (ContinuumConfig.order_mode = 'cv'; applied 2026-09-12; was: fixed 3/2/2)",
         role="the CV preferred order 1 in half of the fits of the previous run",
         update="applied"),
    dict(priority=8, quantity="detection threshold",
         value="detected >= 3 sigma, marginal 1-3 sigma, limits at 3 sigma (FitConfig; applied "
               "2026-09-12; was: detected >= 1 sigma with 1 sigma limits)",
         role="gates entry into the census and the mixing ratios",
         update="applied"),
    dict(priority=9, quantity="negative-channel cut",
         value="none (FitConfig.drop_negative_sigma = None; applied 2026-09-12; was: channels more than "
               "1 sigma below zero dropped)",
         role="the one-sided cut biased Q upward by ~0.3 sigma per channel and could empty a group",
         update="applied"),
    dict(priority=10, quantity="photometric error column and error model",
         value="source_sum_err_empirical_mjy (applied 2026-09-09) with the continuum-coefficient "
               "covariance propagated in full (generalised least squares, FitConfig.gls; applied 2026-09-12)",
         role="the VARIANCE plane under-reports the scatter; the continuum error is correlated "
              "across a band and was propagated as diagonal before",
         update="applied; a per-target sky_excess_ratio rescaling remains an option"),
    dict(priority=11, quantity="grouping thresholds",
         value="rh_tol 10 %, arc_min_drh 0.001 au, link_drh 0.05 au, manual edges for 24P / 2024E1",
         role="define what counts as one physical state; the manual groups exceed 10 % (to 18.5 %)",
         update="re-examine against the activity light curves once Q(r_h) tracks exist"),
    dict(priority=12, quantity="sufficiency gates",
         value="min_points_phase 5, min_points_emission 2 (ContinuumConfig)",
         role="a group below either is skipped, not FAILed",
         update="tie to the KEY_RANGES minimum channel counts"),
    dict(priority=13, quantity="opacity calibration",
         value="RHO_TAU_REF_KM = 1000 km, KAPPA_PUMP CO2 = 3e-14 cm^2",
         role="modes 1-2 only, never inside the fit; map the optically-thin systematic",
         update="Debout et al. (2016) for rho_tau; PSG pump-line intensities for kappa"),
    dict(priority=14, quantity="rotational temperature",
         value="T_rot = 70 K (ModelParams)",
         role="selects the database template (line intensities within a band); the "
              "band-integrated g-factors change by < 3 % over 30-130 K, so Q is insensitive",
         update="a per-comet value from the literature, or leave fixed"),
    dict(priority=16, quantity="sky annulus (upstream)",
         value="inner radius 150 000 km at the comet, floored at 15 px and capped at 40 px, 5 px wide "
               "(spherex_apphot Config.annulus_r_in_km; applied 2026-09-12; was: fixed 15-20 px)",
         role="the fixed ring sat at 40 000 km for 24P, inside the coma, removing 3.5-9 % of the flux",
         update="applied; the physical radius is a convention -- test 100 000 and 200 000 km"),
    dict(priority=17, quantity="grouping: observer-distance rule",
         value="Delta spread < 20 % inside every group (GroupingConfig.delta_tol; applied 2026-09-12)",
         role="bounds the r_h^2 Delta^2 spread the distance-corrected continuum must absorb",
         update="applied; 20 % is a convention"),
    dict(priority=15, quantity="source-flag thresholds (upstream)",
         value="flag a: G_eff < 13 within r_ap + 2 FWHM; flag b: Gaia flux > 0.2 x comet flux "
               "within r_ap + FWHM (spherex_apphot.Config)",
         role="define the contamination tiers this package compares",
         update="from the survey-wide flag_effectiveness result in the apphot project"),
)
