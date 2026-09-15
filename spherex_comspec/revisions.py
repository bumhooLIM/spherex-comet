"""
Case revisions -- per-(target, phase) adjustments decided by hand.

The review memo of 2026-09-15 (``doc/notes_ver260915.xlsx``; the applied version is
``doc/comspec/case_revisions.md``) went through every fitted phase figure and asked, case
by case, for a different continuum window, a fixed polynomial order, the exclusion of
star-contaminated points, a waiver of the negative-continuum or one-sided-continuum
rejection, a fit on fewer channels than the coverage rule demands, or the rejection of a
detection the figure shows to be spurious.  Those decisions are science, so they live here
as data -- one :class:`CaseRevision` per (target, phase), each band's directives in a
:class:`BandRevision` -- and the pipeline applies them wherever the group is processed
(:func:`continuum.process_group`, :func:`fitting.fit_production_rates`, the figures), so
that the batch run, the study variants and the figures cannot drift apart.  A variant with
``revisions=False`` (``dc_main_norev``) runs without them, which is what measures their effect.

The phase regrouping of 240P ("regroup phases 2-3") is a grouping directive and lives in
:class:`config.GroupingConfig` (``manual_edges`` and ``delta_tol_exempt``); the ZTF Af-rho
directives (47P, 210P, 217P) live in ``ztfcomet.config.AFRHO_EPOCH_OVERRIDES``.  Both are
listed here in ``memo`` fields for the record only.

Semantics of a :class:`BandRevision`
------------------------------------
``cont_lo`` / ``cont_hi``     the continuum window edges of the band [um] (replace ``BAND_WINDOWS``)
``em_lo`` / ``em_hi``         the emission window edges [um]; also the punch-out every other band's
                              continuum applies, and the channels the band contributes to the fit
``order``                     a fixed polynomial order (no cross-validation for this band)
``exclude_top_blue`` / ``_red``  drop the N brightest continuum points on that side before the fit
                              (star-contaminated points the sigma clipping does not catch)
``ignore_negative``           the positivity check is waived: a continuum below zero under the band
                              is fitted and subtracted as it is
``one_sided_ok``              an unbracketed (one-sided) continuum is accepted, and the window is
                              used as given -- no extension to the empty side
``drop_flags``                rows with these source flags are removed from the band (continuum
                              and emission) before anything is fitted
``exclude_top_emission``      the N brightest emission channels of the band are excluded from the
                              production-rate fit (role ``excluded``)
``min_channels``              the species carried by the band is fitted when at least this many of
                              the band's emission channels are present, instead of the
                              ``KEY_RANGES`` rule (coverage, not significance)
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

from .config import BAND_WINDOWS

__all__ = ["BandRevision", "CaseRevision", "CASE_REVISIONS", "REVISIONS", "BAND_OF", "SPECIES_OF",
           "revision_for", "effective_windows", "MEMO_DATE"]

MEMO_DATE = "2026-09-15"
#: the band that carries each species' coverage, and back
BAND_OF = {"H2O": "2.7um", "CO2": "4.3um", "CO": "4.7um"}
SPECIES_OF = {v: k for k, v in BAND_OF.items()}


@dataclass(frozen=True)
class BandRevision:
    cont_lo: Optional[float] = None
    cont_hi: Optional[float] = None
    em_lo: Optional[float] = None
    em_hi: Optional[float] = None
    order: Optional[int] = None
    exclude_top_blue: int = 0
    exclude_top_red: int = 0
    ignore_negative: bool = False
    one_sided_ok: bool = False
    drop_flags: Tuple[str, ...] = ()
    exclude_top_emission: int = 0
    min_channels: Optional[int] = None
    memo: str = ""

    def describe(self) -> str:
        """The directives in words, for the summary tables and the figure notes."""
        p = []
        if self.cont_lo is not None or self.cont_hi is not None:
            p.append("continuum window " + ("%.2f" % self.cont_lo if self.cont_lo is not None else "..")
                     + "-" + ("%.2f" % self.cont_hi if self.cont_hi is not None else "..") + " um")
        if self.em_lo is not None or self.em_hi is not None:
            p.append("emission window " + ("%.2f" % self.em_lo if self.em_lo is not None else "..")
                     + "-" + ("%.2f" % self.em_hi if self.em_hi is not None else "..") + " um")
        if self.order is not None:
            p.append(f"order {self.order} fixed")
        if self.exclude_top_blue:
            p.append(f"{self.exclude_top_blue} brightest blue continuum point(s) excluded")
        if self.exclude_top_red:
            p.append(f"{self.exclude_top_red} brightest red continuum point(s) excluded")
        if self.ignore_negative:
            p.append("negative continuum accepted")
        if self.one_sided_ok:
            p.append("one-sided continuum accepted, no extension")
        if self.drop_flags:
            p.append("flag " + "/".join(self.drop_flags) + " rows dropped")
        if self.exclude_top_emission:
            p.append(f"{self.exclude_top_emission} brightest emission channel(s) excluded from the fit")
        if self.min_channels is not None:
            p.append(f"fitted with >= {self.min_channels} emission channel(s)")
        return "; ".join(p)


@dataclass(frozen=True)
class CaseRevision:
    target: str
    phase: int
    bands: Dict[str, BandRevision] = field(default_factory=dict)
    #: species whose detection the review rejected as spurious (status ``rejected``)
    reject: Tuple[str, ...] = ()
    #: the memo's text for directives applied elsewhere (grouping, ZTF Af-rho)
    memo: str = ""

    def describe(self) -> str:
        parts = [f"{b}: {r.describe()}" for b, r in self.bands.items() if r.describe()]
        if self.reject:
            parts.append("rejected: " + ", ".join(self.reject))
        return "; ".join(parts)

    def coverage(self) -> Dict[str, int]:
        """Species -> minimum number of emission channels that replaces the ``KEY_RANGES`` rule."""
        return {SPECIES_OF[b]: int(r.min_channels) for b, r in self.bands.items()
                if r.min_channels is not None and b in SPECIES_OF}

    def band(self, name: str) -> Optional[BandRevision]:
        return self.bands.get(name)


def effective_windows(rev: Optional[CaseRevision]) -> dict:
    """``BAND_WINDOWS`` with the case's window edges applied (a deep copy; the module constant
    is never touched)."""
    w = deepcopy(BAND_WINDOWS)
    if rev is None:
        return w
    for b, r in rev.bands.items():
        if b not in w:
            continue
        lo, hi = w[b]["cont"]
        w[b]["cont"] = (r.cont_lo if r.cont_lo is not None else lo, r.cont_hi if r.cont_hi is not None else hi)
        lo, hi = w[b]["em"]
        w[b]["em"] = (r.em_lo if r.em_lo is not None else lo, r.em_hi if r.em_hi is not None else hi)
    return w


def _B(**kw) -> BandRevision:
    return BandRevision(**kw)


NEG = "Ignore continuum to be negatived but fit it. and make PASS"

#: The memo of 2026-09-15, one entry per (target, phase); ``memo`` carries the directives that
#: are applied outside this module.
CASE_REVISIONS: Tuple[CaseRevision, ...] = (
    CaseRevision("2P", 1, {"2.7um": _B(cont_lo=2.40, memo="Move the left window to 2.4um"),
                           "4.7um": _B(cont_hi=4.95, memo="Move the right window to 4.95um")}),
    CaseRevision("2P", 2, {"4.3um": _B(ignore_negative=True, memo=NEG)}),
    CaseRevision("10P", 1, {"2.7um": _B(drop_flags=("b",), memo="Exclude flag b points")}),
    CaseRevision("24P", 7, {"4.7um": _B(order=2, memo="Increase the order 1 -> 2")}),
    CaseRevision("24P", 8, {"4.7um": _B(order=2, memo="Increase the order 1 -> 2")}),
    CaseRevision("43P", 1, {"4.7um": _B(order=1, memo="Reduce the order 2 -> 1")}),
    CaseRevision("47P", 1, memo="ztf_afrho: Extrapolate the Afrho trend to S1"),
    CaseRevision("47P", 2, {"2.7um": _B(cont_hi=2.95, em_hi=2.85,
                                        memo="Move the right window to 2.95 um (continuum), 2.85 um (emission)")}),
    CaseRevision("63P", 2, {"2.7um": _B(order=1, memo="Reduce the order 2 -> 1")}),
    CaseRevision("124P", 1, {"4.3um": _B(exclude_top_red=3, memo="Exclude top 3 data points in right window"),
                             "4.7um": _B(exclude_top_blue=3, memo="Exclude top 3 data points in left window")}),
    CaseRevision("124P", 3, {"4.7um": _B(order=1, cont_hi=4.95,
                                         memo="Reduce the order 2 -> 1; Move the right window to 4.95 um")}),
    CaseRevision("210P", 2, memo="ztf_afrho: Refit the Afrho trend to S1, S2, S3, and S4."),
    CaseRevision("217P", 1, memo="ztf_afrho: Overestimated. Refit the Afrho to the nearby point."),
    CaseRevision("217P", 3, {"4.3um": _B(ignore_negative=True, memo=NEG)}),
    CaseRevision("229P", 1, {"4.3um": _B(ignore_negative=True, memo=NEG)}),
    CaseRevision("240P", 2, memo="phase_note: Regroup phase 2-3 (GroupingConfig.manual_edges / delta_tol_exempt)"),
    CaseRevision("306P", 4, {"2.7um": _B(order=1, memo="Reduce the order 2 -> 1"),
                             "4.3um": _B(exclude_top_emission=3,
                                         memo="Exclude the top 3 data points before emission fitting.")}),
    CaseRevision("2019U5", 1, {"4.3um": _B(ignore_negative=True, memo=NEG)}),
    CaseRevision("2019U5", 2, {"2.7um": _B(exclude_top_blue=2, exclude_top_red=2, ignore_negative=True,
                                           memo="Exclude top 2 data points in right window and 2 data points "
                                                "in left window before continuum fit; " + NEG)}),
    CaseRevision("2022N2", 2, {"4.3um": _B(em_lo=4.13, min_channels=2,
                                           memo="Increase the left side emission window +0.05 um. Make a emission fit.")}),
    CaseRevision("2022QE78", 1, {"2.7um": _B(min_channels=2, memo="Make emission fit with 2 channels")}),
    CaseRevision("2022R6", 2, {"4.3um": _B(exclude_top_blue=4, memo="Exclude top 4 data points in left window")}),
    CaseRevision("2023A3", 2, {"4.3um": _B(exclude_top_red=1, memo="Exclude top 1 data point in right window")}),
    CaseRevision("2023R1", 1, {"2.7um": _B(cont_lo=2.40, memo="Move the left window to 2.4um")}),
    CaseRevision("2023R1", 4, {"2.7um": _B(cont_lo=2.00, cont_hi=2.50, order=1, one_sided_ok=True, min_channels=3,
                                           memo="Make 1D continuum fit only with left window from 2.0 to 2.5 um; "
                                                "Make emission fit with 3 channels")}),
    CaseRevision("2023RS61", 1, reject=("H2O",), memo="h2o_emission: This is fault detection. Reject the detection."),
    CaseRevision("2023U1", 2, {"4.7um": _B(em_lo=4.50, memo="Increase the emission window to 4.5 um.")}),
    CaseRevision("2023V1", 2, {"4.3um": _B(exclude_top_blue=3, memo="Exclude top 3 data points in left window")}),
    CaseRevision("2024A1", 1, {"4.3um": _B(exclude_top_blue=5, memo="Exclude top 5 data points in left window")}),
    CaseRevision("2024E1", 2, {"4.3um": _B(ignore_negative=True, order=1, memo=NEG + "; Reduce the order 2->1")}),
    CaseRevision("2024J2", 1, {"4.3um": _B(min_channels=3, memo="Make emission fit with 3 channels")}),
    CaseRevision("2024L5", 1, {"4.3um": _B(min_channels=5, memo="Make emission fit with 5 channels")}),
    CaseRevision("2024L5", 4, {"4.3um": _B(ignore_negative=True, memo=NEG)}),
    CaseRevision("2025L1", 1, {"4.3um": _B(min_channels=1,
                                           memo="Make emission fit with single channel + nearby 2 continuum point "
                                                "as a baseline (the two-point continuum through the nearest points)")}),
    CaseRevision("2025R1", 1, {"4.3um": _B(min_channels=2, memo="Make emission fit with 2 channels")}),
    CaseRevision("2025R1", 4, {"4.3um": _B(exclude_top_red=4, memo="Exclude top 4 data points in right window")},
                 reject=("H2O",), memo="h2o_emission: This is fault detection. Reject the detection."),
    CaseRevision("2025R2", 2, {"4.3um": _B(one_sided_ok=True, order=1,
                                           memo="Make the continuum fit only with right window and make PASS")}),
    CaseRevision("2025W2", 2, {"4.3um": _B(min_channels=2, memo="Make emission fit with 2 channels")}),
)

REVISIONS: Dict[Tuple[str, int], CaseRevision] = {(r.target, r.phase): r for r in CASE_REVISIONS}
assert len(REVISIONS) == len(CASE_REVISIONS), "duplicate (target, phase) in CASE_REVISIONS"


def revision_for(target: str, phase, enabled: bool = True) -> Optional[CaseRevision]:
    """The case revision of ``(target, phase)`` -- ``None`` when there is none, or when the
    variant runs without revisions."""
    if not enabled:
        return None
    key = (str(target).replace(" ", ""), int(phase))
    rev = REVISIONS.get(key)
    if rev is not None and not rev.bands and not rev.reject:
        return None            # a memo-only entry (grouping or ZTF directive): nothing to apply here
    return rev
