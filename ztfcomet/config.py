"""Configuration: observatory constants, survey calibration data, target records.

Every magic number that used to be pasted into a notebook cell lives here.  That
includes the JPL Horizons orbit record numbers, which are the single most
error-prone constant in the project: ``ztfquery_2P.ipynb`` shipped with
``target_id = 90000355  # 24P/Schaumasse`` in a notebook analysing 2P, so every
ephemeris annotation on its figures described the wrong comet.  Naming a target
once, here, makes that class of mistake impossible.
"""

from __future__ import annotations

import datetime
import logging

import numpy as np
from dataclasses import dataclass, field, replace

_log = logging.getLogger(__name__)

__all__ = [
    "ZTF_OBSCODE", "IRSA_SEARCH_URL", "IRSA_DATA_URL",
    "SOLAR_APPMAG_AB", "ZTF_FILTERS", "EPHEM_QUANTITIES_COARSE",
    "EPHEM_QUANTITIES_FULL", "QueryConfig", "PhotConfig", "Target",
    "TARGETS", "get_target",
]

#: Default start of the survey window.  ZTF cutouts before this are not part of
#: the current programme; override per target or with ``--start``.
DEFAULT_START_DATE = "2025-03-01"

#: Aperture radii at the comet, in km.  Af-rho is aperture-dependent, so a
#: survey reports several and quotes rho with every value.
RHO_KM_SET = (10_000.0, 15_000.0, 20_000.0, 30_000.0, 40_000.0)


def today() -> str:
    """Today's date as ``YYYY-MM-DD``.

    The default query window runs to *present*, so the end date is evaluated
    when a :class:`Target` is created rather than frozen into the source.
    """
    return datetime.date.today().isoformat()


#: MPC observatory code for Palomar / ZTF.
ZTF_OBSCODE = "I41"

#: IRSA IBE metadata search endpoint for ZTF science images.
IRSA_SEARCH_URL = "https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci"

#: IRSA IBE data endpoint; cutouts are this plus ``?center=&size=``.
IRSA_DATA_URL = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci"

#: Apparent AB magnitude of the Sun in the PS1 system (Willmer 2018, ApJS 236, 47).
#: ZTF photometry is calibrated onto PS1, so these are the correct reference
#: fluxes for the Af-rho flux ratio.
SOLAR_APPMAG_AB = {
    "ZTF_g": -26.54,
    "ZTF_r": -26.93,
    "ZTF_i": -27.05,
}

#: ``filtercode`` (IRSA metadata) to ``FILTER`` (FITS header) for each ZTF band.
ZTF_FILTERS = {"zg": "ZTF_g", "zr": "ZTF_r", "zi": "ZTF_i"}

#: Horizons quantities for the coarse survey pass: RA/Dec, rates, mag, helio dist.
EPHEM_QUANTITIES_COARSE = "1,3,9,19"

#: Horizons quantities at exact frame times, adding phase angle, elongation,
#: geocentric distance and the sun-target / velocity position angles.
#: Quantity 23 supplies ``elong``; the predecessor set ("1,3,4,8,9,16,18,19,20,24,27")
#: omitted it, and the old notebook only obtained elongation by calling
#: ``ephemerides()`` with no quantities at all and taking the service defaults.
EPHEM_QUANTITIES_FULL = "1,3,4,8,9,16,18,19,20,23,24,27"


@dataclass(frozen=True)
class QueryConfig:
    """Parameters controlling the ephemeris and image search.

    Attributes
    ----------
    interval_days : float
        Step of the coarse ephemeris, and the half-width of the ``obsjd``
        window searched around each step.  Also scales the search box, which is
        sized from the target's own apparent motion over this interval.
    rh_max, vmag_max : float
        Cuts applied to the coarse ephemeris before any image search.  The
        default ``vmag_max`` of 19 keeps epochs where the comet is plausibly
        detectable in a 30 s ZTF exposure (5-sigma limit ~20.5).
    cutout_size : str
        Fixed IRSA cutout size, used when ``adaptive_cutout`` is off.
    adaptive_cutout : bool
        Size each target's cutouts from its own ephemeris rather than using a
        fixed size.  A distant, faint comet needs only a small box; a bright or
        nearby one has a coma spanning far more sky.
    cutout_size_small, cutout_size_large : str
        The two sizes the adaptive choice picks between.
    bright_vmag : float
        A target reaching this magnitude or brighter at any epoch gets the large
        cutout.
    close_delta_au : float
        A target coming this close to the observer at any epoch gets the large
        cutout: angular coma size scales as 1/delta, so a nearby comet overflows
        a small box even when it is faint.
    is_cutout : bool
        When ``False``, full 3072x3080 quadrant images are requested instead of
        cutouts.  These are ~37 MB each; the SSD holds several such datasets.
    max_epochs_per_call : int
        astroquery sends Horizons epoch lists in the request URL, so the binding
        limit is URL length, not a documented epoch count.  Measured against the
        live service: 75 JDs succeed, 100 return HTTP 502.  Requests are chunked
        at this size, kept well below that.
    timeout, max_retries, backoff : float, int, float
        HTTP behaviour for IRSA and Horizons calls.
    """

    interval_days: float = 5.0
    rh_max: float = 10.0
    vmag_max: float = 19.0
    is_cutout: bool = True
    cutout_size: str = "10arcmin"
    adaptive_cutout: bool = True
    cutout_size_small: str = "5arcmin"
    cutout_size_large: str = "10arcmin"
    bright_vmag: float = 14.0
    close_delta_au: float = 1.0
    intersect: str = "overlaps"
    max_epochs_per_call: int = 50
    timeout: float = 120.0
    max_retries: int = 3
    backoff: float = 2.0


@dataclass(frozen=True)
class PhotConfig:
    """Parameters controlling aperture photometry and the Af-rho reduction.

    Attributes
    ----------
    rho_km : float
        Physical radius of the photometric aperture at the comet, in km.  Af-rho
        is aperture-dependent, so this must be reported with every result.
    rho_km_set : tuple of float
        Aperture radii for a multi-aperture reduction.  Each frame is measured
        at every radius that satisfies the scale test below.
    max_aperture_arcsec : float
        Upper bound on the on-sky aperture radius.  Together with the seeing
        FWHM this bounds a usable aperture:
        ``FWHM_arcsec < r_ap_arcsec < max_aperture_arcsec``.  Below the FWHM the
        aperture does not contain the PSF; above an arcminute it is dominated by
        sky and, in a 5 arcmin cutout, runs out of frame.
    sky_in_scale, sky_out_scale, sky_out_pad : float
        Sky annulus geometry: ``r_in = sky_in_scale * rho_pix`` and
        ``r_out = sky_out_scale * rho_pix + sky_out_pad`` (pixels).
    phase_beta, phase_beta_err : float
        Linear phase coefficient in mag/deg used to reduce to A(0 deg) f rho,
        and its uncertainty.  The default beta is a common dust-coma value, not
        a measurement for any particular comet; set it per target when better
        constrained.
    default_color_gr : float
        Comet ``g-r`` assumed when the colour term cannot be measured from a
        same-night pair.  Typical dust-dominated comae sit near 0.5.
    contam_flux_ratio : float
        Flag a frame when catalogued background sources inside the aperture
        carry at least this fraction of the comet's expected flux.  0.30 means
        ``G_eff <= Tmag + 1.31``.
    contam_radius_pad_fwhm : float
        The contamination search radius is ``rho_pix + contam_radius_pad_fwhm *
        FWHM`` pixels: a star just outside the aperture still spills flux into
        it through the PSF wings.
    check_contamination : bool
        Run the Gaia check at all.  Turned off automatically when no catalogue
        is available.
    min_rho_fwhm : float
        Quality threshold: apertures smaller than this many seeing FWHM do not
        contain the PSF, so the "rho_km aperture" is not measuring rho_km of
        coma.  Frames below it are flagged, never silently dropped.
    sigma_clip_sigma, sigma_clip_iters : float, int
        Sky background clipping.
    winpos_sig_scale : float
        ``sep.winpos`` Gaussian sigma, in units of the seeing FWHM.
    max_centroid_shift_fwhm : float
        Flag frames where the refined centroid moved further than this many FWHM
        from the ephemeris position — usually a field star capturing the centroid.
    anomaly_window_days, anomaly_neighbours, anomaly_min_dex, anomaly_sigma
        The single-frame anomaly test (:func:`ztfcomet.phot.flag_anomalies`):
        a frame is compared with up to *anomaly_neighbours* otherwise-clean
        frames on each side within *anomaly_window_days*, and flagged when it
        exceeds their median by more than *anomaly_min_dex* and more than
        *anomaly_sigma* times the robust scatter of such excesses, with
        neither adjacent frame sharing half the excess.
    """

    rho_km: float = 15_000.0
    rho_km_set: tuple = RHO_KM_SET
    max_aperture_arcsec: float = 60.0
    sky_in_scale: float = 3.0
    sky_out_scale: float = 4.0
    sky_out_pad: float = 20.0
    phase_beta: float = 0.03
    phase_beta_err: float = 0.0
    default_color_gr: float = 0.5
    min_rho_fwhm: float = 1.5
    contam_flux_ratio: float = 0.30
    contam_radius_pad_fwhm: float = 1.0
    check_contamination: bool = True
    sigma_clip_sigma: float = 3.0
    sigma_clip_iters: int = 5
    winpos_sig_scale: float = 3.0
    max_centroid_shift_fwhm: float = 3.0
    anomaly_window_days: float = 30.0
    anomaly_neighbours: int = 3
    anomaly_min_dex: float = 0.15
    anomaly_sigma: float = 5.0
    apply_color_term: bool = True
    apply_aperture_correction: bool = True


@dataclass(frozen=True)
class Target:
    """A comet, with everything needed to query and reduce it.

    Attributes
    ----------
    name : str
        Short designation used for directory names and plot labels, e.g. ``"24P"``.
    horizons_id : str or int
        Identifier passed to :class:`astroquery.jplhorizons.Horizons`.  For
        periodic comets a bare designation is ambiguous across apparitions, so
        prefer the numeric orbit record (e.g. ``90001204``).
    designation : str
        The name to resolve against Horizons, e.g. ``"240P"``.  Defaults to
        :attr:`name`.  Prefer this over a hardcoded record number: Horizons
        **renumbers** its small-body records, and the numbers the pre-merge
        notebooks carried for 240P now resolve to 233P and 234P — different
        comets.  See :mod:`ztfcomet.horizons`.
    allow_fragment : bool
        Permit Horizons to return a fragment (``240P-B``).  Default ``False``:
        asking for ``240P`` must give the parent body, and the predecessor's
        "take the last record" rule silently selected the fragment.
    orbit_records : dict
        Optional ``{jd_from: record_number}`` override for pinning specific
        apparitions.  **Deprecated** — record numbers go stale; leave it empty
        and let :mod:`ztfcomet.horizons` resolve per epoch.
    perihelion_jd : dict
        Optional ``{label: JD}`` times of perihelion, for ``T - T_p`` lightcurves.
    """

    name: str
    horizons_id: str | int | None = None
    designation: str | None = None
    allow_fragment: bool = False
    start_date: str = DEFAULT_START_DATE
    end_date: str = field(default_factory=today)
    orbit_records: dict[float, int] = field(default_factory=dict)
    perihelion_jd: dict[str, float] = field(default_factory=dict)
    query: QueryConfig = field(default_factory=QueryConfig)
    phot: PhotConfig = field(default_factory=PhotConfig)
    note: str = ""

    @property
    def query_designation(self) -> str:
        """Name to resolve against Horizons."""
        return str(self.designation or self.horizons_id or self.name)

    def resolve_orbit_record(self, jd: float):
        """Horizons id appropriate for an observation at *jd*.

        Resolves the designation dynamically (see :mod:`ztfcomet.horizons`),
        excluding fragments unless :attr:`allow_fragment` is set and preferring
        the orbit solution nearest the observation.

        An explicit :attr:`orbit_records` entry overrides the lookup, but that
        route is deprecated: Horizons renumbers records, so a number pinned
        today can name a different comet next year.
        """
        if self.orbit_records:
            applicable = [b for b in sorted(self.orbit_records) if jd >= b]
            record = self.orbit_records[applicable[-1] if applicable
                                        else min(self.orbit_records)]
            _log.warning("%s: using pinned Horizons record %s. Record numbers "
                         "are not stable — prefer designation lookup.",
                         self.name, record)
            return record

        from . import horizons as _horizons        # late: avoids a cycle
        return _horizons.resolve_target_id(
            self.query_designation, epoch_jd=jd, allow_fragment=self.allow_fragment)

    def with_dates(self, start_date: str, end_date: str) -> "Target":
        """Copy of this target over a different date range."""
        return replace(self, start_date=start_date, end_date=end_date)


#: Known targets.  Add entries here rather than pasting IDs into notebooks.
#:
#: Note the absence of hardcoded Horizons record numbers.  The pre-merge
#: notebooks carried ``90001203``/``90001204`` labelled "240P/NEAT"; those
#: records today resolve to 233P/La Sagra and 234P/LINEAR.  Designations are
#: stable, record numbers are not.
TARGETS: dict[str, Target] = {
    "24P": Target(
        name="24P",
        designation="24P",
        query=QueryConfig(interval_days=5, rh_max=9),
        note="24P/Schaumasse. Default test target for query.ipynb and afrho.ipynb.",
    ),
    "240P": Target(
        name="240P",
        designation="240P",
        start_date="2018-01-01",
        end_date="2026-06-30",
        # "240P" is ambiguous: Horizons lists two parent solutions (2014, 2024)
        # and the fragment 240P-B. allow_fragment=False keeps the parent, and
        # the solution nearest each observation is chosen per epoch.
        allow_fragment=False,
        perihelion_jd={"2018": 2458256.7471537665, "2025": 2461029.353067453},
        note="240P/NEAT (parent body; 240P-B is a fragment and is excluded).",
    ),
    "240P-B": Target(
        name="240P-B",
        designation="240P-B",
        allow_fragment=True,
        start_date="2024-01-01",
        end_date="2026-06-30",
        note="Fragment of 240P/NEAT. Requested explicitly, never by resolving 240P.",
    ),
    "2P": Target(
        name="2P",
        designation="2P",
        query=QueryConfig(interval_days=5, rh_max=9),
        note="2P/Encke.",
    ),
    "2014UN271": Target(
        name="2014UN271",
        designation="2014 UN271",
        # Bernardinelli-Bernstein: Tmag 16.1-16.6 at r_h 13.6-15.3 au through
        # the survey window.  The default rh_max of 10 au silently removed every
        # epoch before the magnitude cut ran; for this comet only the magnitude
        # cut applies.  With the cut lifted the IRSA search still returns
        # nothing: the comet is at declination -65 to -75 deg throughout
        # 2025-2026, and ZTF reaches only to about -30 deg.
        query=QueryConfig(interval_days=5, rh_max=100.0, vmag_max=20),
        note="C/2014 UN271 (Bernardinelli-Bernstein): bright at 14 au, no r_h cut; "
             "unobservable from Palomar (dec < -65 deg) in 2025-2026.",
    ),
    "2019Y3": Target(
        name="2019Y3",
        designation="2019 Y3",
        start_date="2025-01-01",
        end_date="2025-10-30",
        query=QueryConfig(interval_days=5, rh_max=10, vmag_max=20),
        note="C/2019 Y3 (ATLAS).",
    ),
    "2024E1": Target(
        name="2024E1",
        designation="2024 E1",
        start_date="2025-01-01",
        end_date="2025-10-30",
        note="C/2024 E1 (Wierzchos).",
    ),
}


def get_target(name: str) -> Target:
    """Look up a :class:`Target` by name, case- and whitespace-insensitively.

    Unknown names produce a bare :class:`Target` that passes *name* straight to
    Horizons, so ad-hoc objects still work without editing this file.
    """
    key = "".join(str(name).split())
    for candidate, target in TARGETS.items():
        if candidate.lower() == key.lower():
            return target
    return Target(name=key, designation=str(name),
                  note="Not in TARGETS; designation resolved against Horizons.")


# --------------------------------------------------------------------------- SPHEREx epochs
#: Per-(target, SPHEREx phase) adjustments of the Af-rho estimate at the SPHEREx epochs
#: (``activity.afrho_at_epoch`` through ``scripts/ztf/afrho_trends.py``), from the review memo
#: of 2026-09-15 (``doc/comspec/case_revisions.md``).  Keys are ``(target, phase)``; the values
#: are keyword arguments of :func:`ztfcomet.activity.afrho_at_epoch` plus a ``memo`` that is
#: appended to the estimate's ``note``:
#:
#: * ``leg``: evaluate the epoch with the fitted law of this leg whatever side of T_p it falls
#:   on -- 210P was observed by ZTF only after perihelion, and the memo asks for its four
#:   inbound SPHEREx phases to be read off the outbound law (a symmetric-activity assumption,
#:   said so in the note);
#: * ``extrap_grades``: laws that may be extrapolated -- 47P's rising law is grade D (a 0.09 au
#:   baseline) and S1 lies 0.005 dex beyond it; the memo asks for the extrapolation;
#: * ``exclude_outburst``: frames inside a detected outburst window do not enter the direct
#:   (in-window) mean -- 217P S1 had two of its five frames at the start of the 2.30 au
#:   outburst, which overestimated the quiescent value.
AFRHO_EPOCH_OVERRIDES = {
    ("47P", 1): dict(extrap_grades=("A", "B", "C", "D"),
                     memo="memo 2026-09-15: extrapolate the (grade D) rising law to S1"),
    ("210P", 1): dict(leg="outbound", memo="memo 2026-09-15: inbound epoch read off the outbound law"),
    ("210P", 2): dict(leg="outbound", memo="memo 2026-09-15: inbound epoch read off the outbound law"),
    ("210P", 3): dict(leg="outbound", memo="memo 2026-09-15: inbound epoch read off the outbound law"),
    ("210P", 4): dict(leg="outbound", memo="memo 2026-09-15: inbound epoch read off the outbound law"),
    ("217P", 1): dict(exclude_outburst=True,
                      memo="memo 2026-09-15: the in-window frames at the outburst start are excluded"),
}


def afrho_epoch_override(target, phase):
    """The :data:`AFRHO_EPOCH_OVERRIDES` entry of ``(target, phase)`` (a copy), or ``{}``."""
    key = ("".join(str(target).split()), int(phase))
    return dict(AFRHO_EPOCH_OVERRIDES.get(key, {}))

