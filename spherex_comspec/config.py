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
    "KEY_RANGES", "BAND_CARRIER", "BAND_COLORS", "EMISSION_DTYPES", "SOURCEFLAG_PRIORITY",
    "RHO_TAU_REF_KM", "KAPPA_PUMP", "PLACEHOLDERS",
    "GroupingConfig", "ApertureConfig", "FlagPolicy", "ContinuumConfig", "ModelParams",
    "FitConfig", "Variant", "DEFAULT_VARIANTS", "MAIN_VARIANT", "VARIANTS", "config_hash",
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
@dataclass(frozen=True)
class Band:
    """One vibrational emission band (Ootsubo et al. 2012, Table 2)."""

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
BAND_WINDOWS: Dict[str, dict] = {
    "2.7um": dict(em=(2.60, 2.80), cont=(2.30, 3.00), lam_c=2.70),
    "4.3um": dict(em=(4.18, 4.35), cont=(4.00, 4.55), lam_c=4.26),
    "4.7um": dict(em=(4.55, 4.90), cont=(4.40, 5.00), lam_c=4.67),
}
EMISSION_WINDOWS = {b: w["em"] for b, w in BAND_WINDOWS.items()}
ALL_EM_WINDOWS = [w["em"] for w in BAND_WINDOWS.values()]

#: Diagnostic ranges that must be sampled before a species may be fitted at all.
#: Falling short is a *coverage* failure, categorically different from a non-detection.
KEY_RANGES = {
    "H2O": dict(lo=2.60, hi=2.80, min_points=3),
    "CO2": dict(lo=4.20, hi=4.30, min_points=2),
    "CO":  dict(lo=4.60, hi=4.70, min_points=2),
}

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

    1. **r_h range (hard).**  ``(max - min) / mean < rh_tol`` inside every group.
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

    poly_orders: Dict[str, int] = field(default_factory=lambda: {"2.7um": 3, "4.3um": 2, "4.7um": 2})
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

    def v_g(self, r_h_au):
        import numpy as np
        if self.v_g_kms is not None:
            return np.full_like(np.asarray(r_h_au, dtype=float), float(self.v_g_kms))
        return 0.8 * np.asarray(r_h_au, dtype=float) ** -0.5


@dataclass(frozen=True)
class FitConfig:
    """Options of the weighted linear least-squares solve."""

    scale_errors_by_chi2: bool = True
    upper_limit_sigma: float = 1.0
    drop_negative_sigma: Optional[float] = 1.0
    require_key_coverage: bool = True
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
    error_column: str = "source_sum_err_mjy"     # or "source_sum_err_empirical_mjy"
    grouping: GroupingConfig = field(default_factory=GroupingConfig)
    aperture: ApertureConfig = field(default_factory=ApertureConfig)
    continuum: ContinuumConfig = field(default_factory=ContinuumConfig)
    fit: FitConfig = field(default_factory=FitConfig)
    model: ModelParams = field(default_factory=ModelParams)

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
        return config_hash(self.to_dict())


def config_hash(d: dict) -> str:
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()[:12]


# --------------------------------------------------------------------------- the variants
#: The runs the analysis compares.  ``dc_*`` fit the continuum in distance-corrected space;
#: ``raw_all`` is the same policy in physical flux space, for the distance-correction study.
DEFAULT_VARIANTS: Tuple[Variant, ...] = (
    Variant("dc_all",   FlagPolicy("all"), use_distcorr=True),
    Variant("dc_no_a",  FlagPolicy("no_a", drop_flags=("a",)), use_distcorr=True),
    Variant("dc_no_b",  FlagPolicy("no_b", drop_flags=("b",)), use_distcorr=True),
    Variant("dc_no_ab", FlagPolicy("no_ab", drop_flags=("a", "b")), use_distcorr=True),
    Variant("raw_all",  FlagPolicy("all"), use_distcorr=False),
    Variant("dc_all_lenient", FlagPolicy("all_lenient", max_frac_badpix=0.05), use_distcorr=True),
)
#: The variant whose products are the catalog's main result.
MAIN_VARIANT = "dc_all"
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
    dict(priority=2, quantity="band profiles Phi_b (8 bands)",
         value="Gaussian, FWHM 0.020-0.100 um (config.BANDS)",
         role="intrinsic band envelope; unit-normalised, so band-integrated flux and Q are "
              "unaffected at SPHEREx resolution, but any statement about band *shape* is not",
         update="PSG / GSFC Fluorescence Database templates at T_rot, unit area"),
    dict(priority=3, quantity="expansion-velocity law",
         value="v_g = 0.8 r_h^-0.5 km/s (ModelParams.v_g)",
         role="fixed, not fitted (F ~ Q/v_g is exactly degenerate); the largest single "
              "systematic in absolute Q, ~ +/-18 % per 20 % in v_g for CO2",
         update="species-specific velocity from a literature compilation, with its citation"),
    dict(priority=4, quantity="2.7 um emission window blue edge",
         value="2.60 um (config.BAND_WINDOWS)",
         role="loses 6.7 % of the convolved H2O complex and leaves a cont_used channel at 2.59 um "
              "carrying 23 % of the peak; the project handoff recommends 2.50 um",
         update="2.50 um after re-checking one-sided continuum rates (14 groups go one-sided)"),
    dict(priority=5, quantity="badphot policy",
         value="drop any aperture containing a bad pixel (FlagPolicy.drop_badphot)",
         role="the revised photometry flags *any* bad pixel; for close, bright comets this removes "
              "~50 % of rows.  frac_badpix_ap is in the data for a threshold instead",
         update="max_frac_badpix ~ 0.05 once the flux loss per bad pixel has been characterised"),
    dict(priority=6, quantity="aperture rule",
         value="20 000 km inside 3 au, 40 000 km beyond; promoted to the smallest valid "
               "km aperture with >= 95 % coverage (ApertureConfig)",
         role="sets the coma column the model integrates; a target's Q are aperture-consistent "
              "only because one radius is used per target",
         update="an S/N-driven choice per target, or a fixed angular aperture"),
    dict(priority=7, quantity="continuum polynomial orders",
         value="{2.7um: 3, 4.3um: 2, 4.7um: 2} (ContinuumConfig.poly_orders)",
         role="83 of 144 WARN verdicts in the previous run were only CV preferring another order",
         update="re-tune from the cv_best_order column of the emission summaries"),
    dict(priority=8, quantity="detection threshold",
         value="upper_limit_sigma = 1.0 (FitConfig)",
         role="'detected' currently means >= 1 sigma; it gates entry into the mixing ratios",
         update="3 sigma, or rename the tier -- see doc/fitting_methodology.md section 10.2b"),
    dict(priority=9, quantity="negative-channel cut",
         value="drop_negative_sigma = 1.0 (FitConfig)",
         role="one-sided cut biases Q upward by ~ +0.29 sigma per channel for pure noise",
         update="quantify jointly with the 1-sigma detection tier by noise injection"),
    dict(priority=10, quantity="photometric error column",
         value="source_sum_err_mjy (Variant.error_column)",
         role="the revised photometry reports sky_excess_ratio ~ 1.2, i.e. the VARIANCE plane "
              "under-reports the true scatter; this error is a lower bound",
         update="source_sum_err_empirical_mjy, or scale by sky_excess_ratio"),
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
         role="reaches only the band profile, which is unresolved at R <= 130",
         update="irrelevant until real band templates exist"),
    dict(priority=15, quantity="source-flag thresholds (upstream)",
         value="flag a: G_eff < 13 within r_ap + 2 FWHM; flag b: Gaia flux > 0.2 x comet flux "
               "within r_ap + FWHM (spherex_apphot.Config)",
         role="define the contamination tiers this package compares",
         update="from the survey-wide flag_effectiveness result in the apphot project"),
)
