"""Resolve a comet designation to the right JPL Horizons orbit record.

Two problems make this harder than passing a designation straight to Horizons,
and both silently produce ephemerides for the *wrong object*:

**Fragments share the parent's designation.**  Asking Horizons for ``240P``
returns an ambiguity listing:

.. code-block:: text

    Record #  Epoch-yr  >MATCH DESIG<  Primary Desig  Name
    90001211    2014    240P           240P            NEAT
    90001212    2024    240P           240P            NEAT
    90001213    2025    240P-B         240P-B          NEAT     <- fragment

Taking the *last* record — which is what the predecessor's
``extract_lastrecnum`` did — selects ``240P-B``, a fragment of the parent
comet, not the parent itself.  Nothing in the returned table announces the
substitution; only ``targetname`` records it, and the old code never looked.

**Record numbers are not stable.**  Horizons renumbers its small-body records.
The numbers hardcoded in the pre-merge notebooks (``90001203``, ``90001204``,
labelled "240P/NEAT") today resolve to **233P/La Sagra** and **234P/LINEAR** —
different comets entirely.  So a record number is a cache, never an identity:
this module resolves by designation and verifies what came back.

Selection rule: drop fragments (unless explicitly requested), then take the
record whose orbit-solution epoch is nearest the observation.  For 240P the two
parent solutions differ by ~50 arcsec in 2025 — larger than the photometric
aperture — so the epoch choice is not cosmetic.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

from astroquery.jplhorizons import Horizons

from . import config as cfg

__all__ = [
    "HorizonsRecord", "parse_ambiguity_table", "is_fragment_designation",
    "split_designation", "resolve_record", "resolve_target_id",
    "verify_targetname", "clear_cache",
]

log = logging.getLogger(__name__)

#: A cometary fragment carries the parent designation plus a hyphen and one or
#: two capital letters: ``240P-B``, ``73P-C``, ``73P-BB``, ``C/2019 Y4-A``.
_FRAGMENT_RE = re.compile(r"-([A-Z]{1,2})$")

#: A row of the Horizons ambiguity listing.  Name may contain spaces, so the
#: first four fields are matched positionally and the rest is the name.
_RECORD_RE = re.compile(
    r"^\s*(?P<record>\d{6,9})\s+"
    r"(?P<epoch>\d{4}|n\.a\.|-+)\s+"
    r"(?P<match>\S+)\s+"
    r"(?P<primary>\S+)\s*"
    r"(?P<name>.*?)\s*$"
)

#: designation -> resolved record, so a multi-chunk run resolves once.
_CACHE: dict[tuple[str, int | None, bool], "HorizonsRecord"] = {}


@dataclass(frozen=True)
class HorizonsRecord:
    """One candidate orbit record from a Horizons ambiguity listing."""

    record: int
    epoch_yr: int | None
    match_desig: str
    primary_desig: str
    name: str = ""

    @property
    def is_fragment(self) -> bool:
        """True when this record is a fragment rather than the parent body."""
        return is_fragment_designation(self.primary_desig) or \
            is_fragment_designation(self.match_desig)

    @property
    def label(self) -> str:
        return f"{self.primary_desig}/{self.name}".rstrip("/")

    def __str__(self) -> str:
        epoch = self.epoch_yr if self.epoch_yr is not None else "----"
        tag = "  [FRAGMENT]" if self.is_fragment else ""
        return f"{self.record}  {epoch}  {self.label}{tag}"


def is_fragment_designation(designation: str) -> bool:
    """True if *designation* names a fragment of a parent comet.

    Examples
    --------
    >>> is_fragment_designation("240P-B")
    True
    >>> is_fragment_designation("73P-BB")
    True
    >>> is_fragment_designation("240P")
    False
    >>> is_fragment_designation("C/2024 E1")
    False
    """
    return bool(_FRAGMENT_RE.search(str(designation).strip()))


def split_designation(designation: str) -> tuple[str, str | None]:
    """Split a designation into its parent part and fragment letter.

    >>> split_designation("240P-B")
    ('240P', 'B')
    >>> split_designation("240P")
    ('240P', None)
    """
    text = str(designation).strip()
    match = _FRAGMENT_RE.search(text)
    if not match:
        return text, None
    return text[: match.start()], match.group(1)


def parse_ambiguity_table(message: str) -> list[HorizonsRecord]:
    """Parse the record listing out of a Horizons ambiguity error.

    Returns an empty list when *message* is not an ambiguity listing, which is
    how callers distinguish "several candidates" from "the service failed".
    """
    records: list[HorizonsRecord] = []
    for line in str(message).splitlines():
        if "Record" in line and "Desig" in line:      # the header row
            continue
        match = _RECORD_RE.match(line)
        if not match:
            continue
        epoch = match.group("epoch")
        records.append(HorizonsRecord(
            record=int(match.group("record")),
            epoch_yr=int(epoch) if epoch.isdigit() else None,
            match_desig=match.group("match"),
            primary_desig=match.group("primary"),
            name=match.group("name").strip(),
        ))
    return records


def _epoch_from_jd(jd):
    """Approximate calendar year of a Julian Date, for epoch matching."""
    if jd is None:
        return None
    try:
        return 2000.0 + (float(jd) - 2451545.0) / 365.25
    except (TypeError, ValueError):
        return None


def select_record(records, epoch_jd=None, allow_fragment=False,
                  designation=None):
    """Choose the best record from an ambiguity listing.

    Parameters
    ----------
    records : list of HorizonsRecord
    epoch_jd : float, optional
        Observation time.  The record whose orbit epoch is nearest this is
        preferred; without it the most recent epoch wins.
    allow_fragment : bool
        Permit fragment records.  Default ``False`` — the parent body is almost
        always what a comet study means by "240P".
    designation : str, optional
        When given, only records whose parent designation matches are
        considered, so ``240P`` never resolves to ``240P-B`` or to an unrelated
        object that happened to appear in the listing.

    Returns
    -------
    HorizonsRecord or None
    """
    candidates = list(records)
    if not candidates:
        return None

    if designation:
        parent, _ = split_designation(designation)
        wanted = parent.replace(" ", "").upper()
        matched = [
            r for r in candidates
            if split_designation(r.primary_desig)[0].replace(" ", "").upper() == wanted
            or split_designation(r.match_desig)[0].replace(" ", "").upper() == wanted
        ]
        if matched:
            candidates = matched

    if not allow_fragment:
        parents = [r for r in candidates if not r.is_fragment]
        if parents:
            dropped = [r for r in candidates if r.is_fragment]
            if dropped:
                log.info("Ignoring %d fragment record(s): %s",
                         len(dropped), ", ".join(r.primary_desig for r in dropped))
            candidates = parents
        else:
            log.warning("Every candidate for %r is a fragment; using it anyway",
                        designation)

    target_epoch = _epoch_from_jd(epoch_jd)
    dated = [r for r in candidates if r.epoch_yr is not None]
    if dated and target_epoch is not None:
        # Nearest orbit solution to the observation.  Comets carry
        # non-gravitational forces, so a contemporaneous solution beats the
        # newest one: for 240P the two parent records differ by ~50 arcsec in
        # 2025, more than the photometric aperture.
        return min(dated, key=lambda r: abs(r.epoch_yr - target_epoch))
    if dated:
        return max(dated, key=lambda r: r.epoch_yr)
    return candidates[-1]


def resolve_record(designation, epoch_jd=None, allow_fragment=False,
                   location=cfg.ZTF_OBSCODE, use_cache=True):
    """Resolve *designation* to a Horizons record, avoiding fragments.

    Probes Horizons once with the bare designation.  An unambiguous designation
    needs no record number and returns ``None``; an ambiguous one produces the
    candidate listing, which is filtered and ranked by :func:`select_record`.

    Parameters
    ----------
    designation : str or int
        Comet designation (``"240P"``) or an explicit record number, which is
        returned unchanged.
    epoch_jd : float, optional
        Observation time, used to pick the nearest orbit solution.
    allow_fragment : bool
        Set ``True`` only when the fragment really is the target.

    Returns
    -------
    HorizonsRecord or None
        ``None`` means the designation is already unambiguous — pass it straight
        to Horizons.

    Notes
    -----
    Results are cached per ``(designation, epoch decade, allow_fragment)``, so a
    run that queries in chunks resolves once rather than once per chunk.
    """
    text = str(designation).strip()
    if text.isdigit() and len(text) >= 6:
        return None                       # already an explicit record number

    epoch = _epoch_from_jd(epoch_jd)
    key = (text.upper(), int(epoch // 10) if epoch is not None else None, allow_fragment)
    if use_cache and key in _CACHE:
        return _CACHE[key]

    try:
        Horizons(id=text, location=location,
                 epochs=[epoch_jd or 2451545.0]).ephemerides(quantities="1")
        return None                       # unambiguous
    except ValueError as exc:
        records = parse_ambiguity_table(str(exc))
        if not records:
            log.debug("Horizons rejected %r without a candidate listing: %s", text, exc)
            return None
    except Exception as exc:              # noqa: BLE001 — transient; caller retries
        log.debug("Could not probe %r for ambiguity: %s", text, exc)
        return None

    chosen = select_record(records, epoch_jd=epoch_jd,
                           allow_fragment=allow_fragment, designation=text)
    if chosen is None:
        return None

    log.info("Resolved %r -> record %d (%s, epoch %s) from %d candidates",
             text, chosen.record, chosen.label, chosen.epoch_yr, len(records))
    if use_cache:
        _CACHE[key] = chosen
    return chosen


def resolve_target_id(designation, epoch_jd=None, allow_fragment=False,
                      location=cfg.ZTF_OBSCODE):
    """Horizons id to query for *designation* — a record number, or the name.

    Convenience wrapper around :func:`resolve_record` for callers that only need
    something to pass to :class:`~astroquery.jplhorizons.Horizons`.
    """
    record = resolve_record(designation, epoch_jd=epoch_jd,
                            allow_fragment=allow_fragment, location=location)
    return record.record if record is not None else designation


def verify_targetname(targetname, expected, allow_fragment=False):
    """Check that Horizons returned the object we asked for.

    The last line of defence: a record number can be stale, and a fragment can
    slip through if the listing was formatted unexpectedly.  ``targetname`` in
    the returned table is authoritative.

    Parameters
    ----------
    targetname : str
        The ``targetname`` column from a Horizons result, e.g. ``"240P/NEAT"``.
    expected : str
        The designation that was requested, e.g. ``"240P"``.

    Returns
    -------
    (ok, message) : tuple of bool and str
        ``ok`` is False when the result is a fragment (and fragments were not
        requested) or plainly a different object.
    """
    got = str(targetname).strip()
    want = str(expected).strip()
    if not got or not want:
        return True, ""

    # "240P/NEAT" or "240P-B/NEAT" -> designation part
    got_desig = got.split("(")[0].split("/")[0].strip() if "/" in got else got
    if got.startswith("C/") or got.startswith("P/") or got.startswith("X/"):
        got_desig = got.split("(")[0].strip()

    got_parent, got_frag = split_designation(got_desig)
    want_parent, want_frag = split_designation(want)

    if got_frag and not want_frag and not allow_fragment:
        return False, (f"Horizons returned fragment {got_desig!r} for {want!r}; "
                       f"pass allow_fragment=True if that is intended")

    def norm(text):
        return re.sub(r"[^A-Z0-9]", "", str(text).upper())

    if want_parent and norm(want_parent) and norm(got_parent):
        if norm(want_parent) not in norm(got) and norm(got_parent) not in norm(want_parent):
            return False, (f"Horizons returned {got!r} for {want!r} — "
                           f"a stale or wrong record number?")
    return True, ""


def clear_cache():
    """Forget resolved records.  Mostly for tests."""
    _CACHE.clear()
