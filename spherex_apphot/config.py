"""
Every tunable number in the pipeline, in one immutable object.

Rationale
---------
The primitive script hard-coded its aperture list, flag-bit list, magnitude
limit and annulus scaling at three different places in the file, and wrote none
of them into its output.  A CSV produced six months ago could therefore not be
reproduced or even interpreted.  :class:`Config` collects all of it, is frozen
(so a helper cannot quietly mutate a caller's settings), serialises to JSON, and
is written next to every result file together with a short hash that identifies
the parameter set (review items R9 and "no provenance in the outputs").
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

__all__ = ["Config", "FLAG_BIT_MEANINGS", "DEFAULT_BAD_FLAG_BITS"]


# --------------------------------------------------------------------------
# SPHEREx L2 FLAG plane
# --------------------------------------------------------------------------
#: Meaning of each bit of the ``FLAG`` extension.
#:
#: The SPHEREx cutouts carry no key to their own bit assignments, and the
#: primitive pipeline used a bare literal list with no comment.  The mapping
#: below was recovered empirically: for 120 cutouts of 2P the per-bit pixel
#: count was correlated against the frame-level ``L1/L2 N_*`` counters carried
#: in the Parquet index, scaled by the cutout/frame area ratio.  Entries marked
#: ``(inferred)`` reached both a high correlation and a count ratio near unity;
#: entries marked ``(unverified)`` were never populated in the sample and are
#: taken from the ordering of the ``N_*`` header keywords.
#:
#: This table is documentation only -- the pipeline masks whatever
#: :attr:`Config.bad_flag_bits` lists.  Confirm against the SPHEREx Explanatory
#: Supplement before publishing.
FLAG_BIT_MEANINGS: Dict[int, str] = {
    0:  "L1 TRANSIENT (inferred: r=0.70, count ratio 0.97)",
    1:  "L1 OVERFLOW (inferred: count ratio 1.10)",
    2:  "L1 SUR_ERROR (inferred: count ratio 0.91)",
    3:  "L1 PHANTOM (unverified)",
    4:  "L1 REFERENCE (unverified)",
    5:  "L1 NONFUNC (unverified)",
    6:  "L1 DICHROIC / NONFUNC (ambiguous in sample)",
    7:  "L1 MISSING (unverified)",
    8:  "L1 FULLSAMPLE (unverified)",
    9:  "L1 PHANMISS (unverified)",
    10: "reserved (unverified)",
    11: "L2 COLD (inferred: count ratio 0.97)",
    12: "reserved (unverified)",
    14: "reserved (unverified)",
    15: "reserved (unverified)",
    17: "L2 PERSIST (inferred: r=0.80, count ratio 1.04)",
    19: "L2 OUTLIER (inferred: count ratio 1.01)",
    21: "L2 SOURCE -- an astronomical source was detected here",
    22: "reserved (unverified)",
    24: "reserved (unverified)",
    26: "reserved (unverified)",
    27: "reserved (unverified)",
    28: "reserved (unverified)",
    29: "reserved (unverified)",
}

#: Bits treated as unusable pixels.  Inherited unchanged from the primitive
#: pipeline so results stay comparable.
#:
#: Bit 21 (``SOURCE``) and bit 19 (``OUTLIER``) are deliberately *absent*:
#: bit 21 fires on the comet itself, and masking it would delete the target.
DEFAULT_BAD_FLAG_BITS: Tuple[int, ...] = (
    2, 6, 7, 9, 10, 11, 12, 14, 15, 17, 22, 24, 26, 27, 28, 29,
)


def _default_km_apertures() -> List[float]:
    """10000-40000 km in 2000 km steps, then 60000, 80000, 100000 km."""
    return [float(r) for r in range(10_000, 40_001, 2_000)] + [60_000.0, 80_000.0, 100_000.0]


@dataclass(frozen=True)
class Config:
    """
    Immutable parameter set for one pipeline run.

    Notes
    -----
    Attribute defaults encode the decisions taken in response to the review in
    ``doc/apphot/code_review_primitive.md``; each group is annotated with the review
    item it answers.
    """

    # -- photometric apertures (review S2, S3) ------------------------------
    #: Fixed *physical* aperture radii [km].  A fixed cometocentric radius is
    #: what makes different epochs comparable.
    aperture_radii_km: List[float] = field(default_factory=_default_km_apertures)
    #: Fixed *angular* aperture radii [pixel], measured on every exposure
    #: regardless of geometry, for a distance-independent reference point.
    aperture_radii_pix: List[float] = field(default_factory=lambda: [2.0, 5.0])

    #: Sky annulus.  Since 2026-09-12 the inner radius is *physical*:
    #: ``r_in = max(annulus_r_in_km / kmpp, annulus_r_in_pix_min)`` pixels, capped
    #: at ``annulus_r_in_pix_max`` so the ring stays inside the 91-pixel cutout,
    #: with ``r_out = r_in + annulus_width_pix`` (capped at ``annulus_r_out_pix_max``).
    #: The fixed 15-20 px ring used before sat at 40 000 km for 24P and inside
    #: the coma of every close comet, removing 3.5-9 % of the clean-channel flux
    #: (``doc/comspec/apphot_comparison.md`` section 6); 150 000 km is the radius
    #: the earlier photometry used and where a Haser coma of the SPHEREx targets
    #: has fallen well below the sky.  Distant targets, for which 150 000 km is
    #: less than 15 px, keep the near ring, which there is already far outside the
    #: coma.  Set ``annulus_r_in_km = None`` to restore the fixed ``r_in_pix`` /
    #: ``r_out_pix`` ring for every exposure.
    annulus_r_in_km: Optional[float] = 150000.0
    annulus_r_in_pix_min: float = 15.0
    annulus_r_in_pix_max: float = 40.0
    annulus_width_pix: float = 5.0
    annulus_r_out_pix_max: float = 45.0
    #: The fixed ring, used when ``annulus_r_in_km`` is None (and as the floor
    #: of the physical rule, through ``annulus_r_in_pix_min``).
    r_in_pix: float = 15.0
    r_out_pix: float = 20.0

    #: Drop an aperture whose radius is below the PSF FWHM (unresolvable) or at
    #: or beyond ``r_in_pix`` (the aperture would overlap its own background).
    skip_ap_below_psf_fwhm: bool = True
    skip_ap_beyond_r_in: bool = True
    #: Unit in which ``skip_ap_below_psf_fwhm`` compares: "pixel" (default) uses
    #: ``PSF_FWHM / pix_scale``.  The apertures themselves are in pixels, so
    #: "pixel" is the dimensionally consistent choice.
    psf_fwhm_unit: str = "pixel"

    # -- background ---------------------------------------------------------
    sky_sigma: float = 3.0
    sky_maxiters: int = 5
    #: Minimum effective annulus area [pixel^2] for the sky estimate to be used.
    min_sky_area: float = 10.0

    # -- uncertainties (review S6) -----------------------------------------
    #: How the sky contributes to the flux uncertainty.
    #:
    #: ``"level"`` (default)
    #:     ``sigma^2 = sigma_pix^2 + A_eff^2 * ssky^2 / N_sky``.  The variance
    #:     plane already contains the background's photon noise, so adding a
    #:     further ``A_eff * ssky^2`` term (as DAOPHOT does when no variance
    #:     plane exists) would count it twice.  ``ssky`` still enters, through
    #:     the uncertainty on the *sky level* that was subtracted.
    #: ``"empirical"``
    #:     ``sigma^2 = A_eff * ssky^2 + A_eff^2 * ssky^2 / N_sky``.  Ignores the
    #:     variance plane and uses only the measured annulus scatter.
    #: ``"daophot"``
    #:     The primitive behaviour, retained so the change can be quantified.
    #:     Double-counts the background; do not use for science.
    #:
    #: Whichever is selected, all three are also reported column-by-column so a
    #: reader can recombine them.
    sky_noise_mode: str = "level"

    # -- Gaia source flagging (see "Notes for data masking and flagging") ----
    #: Only catalogue entries brighter than this are considered at all.
    gaia_gmag_limit: float = 18.5
    #: Radius of the sky cone loaded around each pointing [deg].
    gaia_search_radius_deg: float = 0.2
    #: Radius for the ``'a'`` test is ``r_ap + bright_psf_mult * PSF_FWHM``.
    sourceflag_bright_psf_mult: float = 2.0
    #: Radius for the ``'b'``/``'c'`` tests and for the reported ``gmag_*``
    #: columns is ``r_ap + near_psf_mult * PSF_FWHM``.
    sourceflag_near_psf_mult: float = 1.0
    #: ``'a'`` fires when the summed (effective) G magnitude inside the bright
    #: radius is brighter than this.
    sourceflag_bright_gmag: float = 13.0
    #: ``'b'`` fires when the blended Gaia flux inside the near radius exceeds
    #: this fraction of the comet's predicted flux::
    #:
    #:     10**(-0.4 * gmag_eff)  >  frac * 10**(-0.4 * vmag)
    #:
    #: i.e. the stars contribute more than ``frac`` of what the comet does.  A
    #: flux ratio, not a magnitude ratio -- the earlier ``vmag > 0.2*gmag_eff``
    #: form compared magnitudes on a scale where a factor of 0.2 has no
    #: photometric meaning, and fired for essentially every faint comet.
    sourceflag_b_flux_frac: float = 0.2
    #: ``'d'`` fires below this signal-to-noise ratio.
    sourceflag_snr_min: float = 1.0

    # -- bad pixels ---------------------------------------------------------
    bad_flag_bits: Tuple[int, ...] = DEFAULT_BAD_FLAG_BITS
    #: Treat non-finite science pixels as bad.  The sample data has ~16 % NaN
    #: science pixels that carry *no* flag bit at all, so this is not optional.
    mask_nonfinite_sci: bool = True
    #: Treat non-finite or non-positive variance as bad.  In the sample these
    #: coincide exactly with the NaN science pixels.
    mask_bad_variance: bool = True

    # -- epoch grouping (review "phase" naming) -----------------------------
    #: A gap larger than this starts a new epoch [day].
    epoch_gap_days: float = 28.0

    # -- photometric zero point --------------------------------------------
    #: Science pixels are mJy/pixel, so an aperture sum is already mJy and the
    #: AB magnitude is ``-2.5*log10(sum*1e-3) + 8.90``.  (Review item S1 argued
    #: for MJy/sr; the data provider confirms mJy/pixel, so no unit conversion
    #: is applied.  ``flux_unit`` is written into every output for the record.)
    flux_unit: str = "mJy/pixel"
    ab_zeropoint: float = 8.90

    # -- stacking (review S8) ----------------------------------------------
    stack_bands: Dict[str, Tuple[float, float]] = field(default_factory=lambda: {
        "Continuum": (1.3, 2.0),
        "H2O": (2.55, 2.8),
        "H2O_ICE": (2.8, 3.3),
        "CO2": (4.15, 4.35),
        "CO": (4.6, 4.8),
    })
    stack_radius_pix: int = 25
    stack_subtract_sky: bool = True
    stack_subpixel_shift: bool = True
    stack_align_north: bool = False
    stack_sigma: float = 3.0
    stack_maxiters: int = 5
    stack_min_frames: int = 3
    #: Statistic used to combine the sigma-clipped stack: ``"median"`` (default)
    #: or ``"mean"``.  The median is the more robust choice against the residual
    #: stars and cosmic rays that survive clipping in a crowded field; the mean
    #: is slightly more efficient on pure Gaussian noise.
    stack_combine: str = "median"

    # -- distance-corrected flux -------------------------------------------
    #: How ``flux_distcorr_mjy`` relates to ``source_sum_mjy``:
    #:
    #: ``"multiply"`` (default)
    #:     ``flux * r_hel**2 * r_obs**2`` -- the flux the comet would have shown
    #:     at ``r_hel = r_obs = 1 au``.  Reflected flux falls as
    #:     ``1 / (r_hel**2 * r_obs**2)``, so this exactly removes the geometric
    #:     dilution and leaves a quantity that is constant for an unchanging
    #:     coma.  It is what makes epochs at different distances comparable.
    #: ``"divide"``
    #:     ``flux / (r_hel**2 * r_obs**2)``.  Retained only for comparison; it
    #:     doubles the distance dependence rather than removing it.
    distcorr_mode: str = "multiply"

    # -- reflectance (review S11) ------------------------------------------
    #: Wavelength window used to normalise the reflectance [um].  Chosen to
    #: avoid the 2.7 um water band.
    refl_norm_window_um: Tuple[float, float] = (1.0, 2.5)
    #: Apply the ``r_hel^2 * r_obs^2`` geometric scaling.  Without it, epochs at
    #: different distances cannot be compared -- this was the S11 defect.
    refl_apply_distance: bool = True
    #: Optional linear phase correction ``10**(0.4 * beta * alpha)``.  Off by
    #: default: the phase function of a coma is target-dependent.
    refl_apply_phase: bool = False
    refl_phase_beta: float = 0.04

    # -- bookkeeping --------------------------------------------------------
    #: Emit no row for an aperture that fails the validity test, rather than a
    #: row of NaNs.  Counts are reported in the log and the status table.
    drop_invalid_apertures: bool = True
    float_format: Optional[str] = "%.8g"
    #: Julian dates need more than eight significant digits: "%.8g" rounds
    #: 2.46e6 to 0.1 day, which the catalog pipeline had to repair from
    #: ``date_obs``.  The ``jd_*`` columns are written with this format instead.
    jd_float_format: str = "%.9f"

    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        """Return a JSON-serialisable copy (tuples become lists)."""
        return json.loads(json.dumps(asdict(self), default=list))

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @property
    def hash(self) -> str:
        """Short stable digest of the parameter set, for provenance."""
        return hashlib.sha256(self.to_json(indent=0).encode()).hexdigest()[:12]

    @classmethod
    def from_dict(cls, d: dict) -> "Config":
        """Rebuild a Config from :meth:`to_dict` output (unknown keys ignored)."""
        known = {f for f in cls.__dataclass_fields__}
        kwargs = {k: v for k, v in d.items() if k in known}
        for key in ("bad_flag_bits", "refl_norm_window_um"):
            if key in kwargs and isinstance(kwargs[key], list):
                kwargs[key] = tuple(kwargs[key])
        if "stack_bands" in kwargs:
            kwargs["stack_bands"] = {k: tuple(v) for k, v in kwargs["stack_bands"].items()}
        return cls(**kwargs)

    def validate(self) -> None:
        """Raise ``ValueError`` on an internally inconsistent parameter set."""
        if self.r_in_pix >= self.r_out_pix:
            raise ValueError(f"r_in_pix ({self.r_in_pix}) must be < r_out_pix ({self.r_out_pix})")
        if self.annulus_r_in_km is not None:
            if self.annulus_r_in_km <= 0:
                raise ValueError("annulus_r_in_km must be positive or None")
            if not (0 < self.annulus_r_in_pix_min <= self.annulus_r_in_pix_max):
                raise ValueError("need 0 < annulus_r_in_pix_min <= annulus_r_in_pix_max")
            if self.annulus_width_pix <= 0 or self.annulus_r_in_pix_max >= self.annulus_r_out_pix_max:
                raise ValueError("need annulus_width_pix > 0 and annulus_r_in_pix_max < annulus_r_out_pix_max")
        if self.sky_noise_mode not in ("level", "empirical", "daophot"):
            raise ValueError(f"unknown sky_noise_mode {self.sky_noise_mode!r}")
        if self.psf_fwhm_unit not in ("pixel", "arcsec"):
            raise ValueError(f"unknown psf_fwhm_unit {self.psf_fwhm_unit!r}")
        if self.stack_combine not in ("median", "mean"):
            raise ValueError(f"unknown stack_combine {self.stack_combine!r}")
        if self.distcorr_mode not in ("divide", "multiply"):
            raise ValueError(f"unknown distcorr_mode {self.distcorr_mode!r}")
        if not self.aperture_radii_km and not self.aperture_radii_pix:
            raise ValueError("no apertures configured")
        for lo, hi in self.stack_bands.values():
            if lo >= hi:
                raise ValueError(f"stack band with lo >= hi: {(lo, hi)}")
