"""JPL Horizons ephemerides and IRSA ZTF image-metadata search.

Ported from ``ztfssoquery.construct_fitsurl`` with four correctness fixes
(C10-C14 of ``doc/primitive_code_analysis.md``):

* The ephemeris block is now joined to the image metadata **on JD**, not by
  row position.  The old ``pd.concat(axis=1)`` was correct only because
  Horizons happens to return epochs sorted and the caller happened to sort the
  metadata first; one change in service behaviour would have silently attached
  every ephemeris to the wrong image.
* Colliding column names are namespaced.  ``airmass`` was delivered by both
  services, so the merged frame carried the label twice and ``df["airmass"]``
  returned a DataFrame.
* :func:`query_sso_ephemeris` no longer raises ``UnboundLocalError`` when the
  record-number retry fails.
* Failures are counted and reported rather than printed and forgotten, so a
  gap in coverage can be told apart from a dropped network call.

The valuable part of the original — recovering the right apparition record from
a Horizons ambiguity error — is preserved in :func:`extract_lastrecnum`.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass, field
from io import BytesIO

import numpy as np
import pandas as pd
import requests
from astropy.io.votable import parse as parse_votable
from astroquery.jplhorizons import Horizons
from tqdm.auto import tqdm

from . import config as cfg
from . import horizons as hz

__all__ = [
    "extract_lastrecnum", "query_sso_ephemeris", "query_ztf_metadata",
    "search_frames", "QueryReport",
]

log = logging.getLogger(__name__)

#: Horizons rounds ``datetime_jd`` in its output, so the JD join needs a
#: tolerance.  1e-5 d = 0.86 s, far below the ZTF cadence and well above the
#: rounding.
_JD_MERGE_TOL = 1e-5


@dataclass
class QueryReport:
    """Per-run accounting, so dropped epochs are visible instead of silent."""

    n_eph_steps: int = 0
    n_eph_kept: int = 0
    n_cut_rh: int = 0        # epochs removed by rh_max
    n_cut_vmag: int = 0      # epochs removed by vmag_max
    n_steps_queried: int = 0
    n_steps_no_coverage: int = 0
    n_steps_failed: int = 0
    n_frames_found: int = 0
    n_frames_containing_target: int = 0
    failures: list[str] = field(default_factory=list)

    def __str__(self) -> str:
        return "\n".join([
            f"ephemeris steps          : {self.n_eph_steps}",
            f"  surviving rh/Vmag cuts : {self.n_eph_kept}",
            f"steps queried at IRSA    : {self.n_steps_queried}",
            f"  no ZTF coverage        : {self.n_steps_no_coverage}",
            f"  FAILED (see .failures) : {self.n_steps_failed}",
            f"frames returned          : {self.n_frames_found}",
            f"  containing the target  : {self.n_frames_containing_target}",
        ])


def extract_lastrecnum(error_msg: str) -> str | None:
    """Recover the most recent orbit record number from a Horizons error.

    Horizons rejects an ambiguous comet designation with a listing of candidate
    apparition records.  The last row is the most recent solution, which is what
    we want for archival ZTF data.  Returns ``None`` when the message is not of
    that form.
    """
    data_lines = [
        line.strip() for line in str(error_msg).split("\n")
        if re.match(r"^\d{8}", line.strip())
    ]
    return data_lines[-1].split()[0] if data_lines else None


def query_sso_ephemeris(target_id, epochs, quantities=cfg.EPHEM_QUANTITIES_COARSE,
                        location=cfg.ZTF_OBSCODE, max_epochs_per_call=50,
                        max_retries=3, backoff=2.0, id_type=hz.SMALLBODY,
                        allow_fragment=False):
    """Query JPL Horizons for ephemerides of a small body.

    Parameters
    ----------
    target_id : str or int
        Designation or orbit record number.
    epochs : dict or array_like
        Either a ``{'start', 'stop', 'step'}`` dict or an explicit list of JDs.
        Long lists are chunked at *max_epochs_per_call*.
    quantities : str
        Comma-separated Horizons quantity codes.
    location : str
        MPC observatory code; defaults to Palomar/ZTF.
    id_type : str
        Passed to Horizons.  Defaults to ``"smallbody"``, which is **not**
        optional for comets: left to guess, Horizons resolves ``"2P"`` to Styx
        (905), a moon of Pluto, and returns a plausible ephemeris for it.
    max_retries, backoff : int, float
        Transient HTTP failures are retried with linear backoff.

    Returns
    -------
    pandas.DataFrame
        Ephemeris table, **one row per unique epoch** and sorted by JD.  Callers
        must join on JD rather than assume a row-for-row match with the input:
        duplicate epochs are collapsed before the request, and a partial table
        is returned when some chunks succeed and others do not.

    Notes
    -----
    astroquery passes the epoch list in the request URL, so a long list fails on
    URL length rather than any documented epoch limit: 75 JDs succeed against
    the live service and 100 return HTTP 502.  Hence the chunking.

    An ambiguous designation triggers one retry against the record number parsed
    from the error message (:func:`extract_lastrecnum`) — the mechanism that
    makes periodic comets work.  Unlike the original implementation, a failed
    retry returns an empty frame instead of falling through to an
    ``UnboundLocalError``.
    """
    if isinstance(epochs, dict):
        chunks = [epochs]
    else:
        jds = np.atleast_1d(np.asarray(epochs, dtype=float))
        jds = jds[np.isfinite(jds)]
        if jds.size == 0:
            return pd.DataFrame()
        # Deduplicate. Horizons rejects a TLIST whose entries are ALL identical
        # with "Bad dates -- start must be earlier than stop", which is what
        # happens when one IRSA step returns several frames from a single
        # exposure (the comet landing on two CCD quadrants). That killed 4 of 65
        # steps in a 24P run. Duplicates are safe to drop because every caller
        # joins the result back on JD, so both frames still get their ephemeris.
        jds = np.unique(jds).tolist()
        chunks = [jds[i:i + max_epochs_per_call]
                  for i in range(0, len(jds), max_epochs_per_call)]

    resolved_id = target_id
    frames, n_failed = [], 0

    for chunk in chunks:
        table = None
        for attempt in range(max_retries):
            try:
                table = Horizons(id=resolved_id, id_type=id_type, location=location,
                                 epochs=chunk).ephemerides(quantities=quantities)
                break
            except ValueError as exc:
                # astroquery raises ValueError both for an ambiguous designation
                # and for a plain service failure ("Query failed without known
                # error message"). Only the first carries a candidate listing;
                # the second is transient and must be retried, not treated as
                # "no such object" -- doing so silently dropped 4 of 65 epochs
                # in a 24P run, reported as "No valid Horizons record".
                candidates = hz.parse_ambiguity_table(str(exc))
                if not candidates:
                    if attempt == max_retries - 1:
                        log.warning("Horizons failed for %d epochs after %d attempts: %s",
                                    len(chunk), max_retries,
                                    str(exc).strip().splitlines()[0][:160])
                    else:
                        time.sleep(backoff * (attempt + 1))
                    continue

                # Ambiguous designation: pick the right apparition and retry.
                # This must go through select_record, not "take the last row" --
                # the last row is the fragment for 240P (see ztfcomet.horizons).
                chosen = hz.select_record(
                    candidates,
                    epoch_jd=(chunk[0] if not isinstance(chunk, dict) and len(chunk)
                              else None),
                    allow_fragment=allow_fragment,
                    designation=str(resolved_id))
                if chosen is None or str(chosen.record) == str(resolved_id):
                    log.warning("No usable Horizons record for %r among %d candidates",
                                resolved_id, len(candidates))
                    return pd.DataFrame()
                log.info("Ambiguous designation %r -> record %s (%s)",
                         resolved_id, chosen.record, chosen.label)
                resolved_id = chosen.record   # remember it for the later chunks
                id_type = None                # a record number needs no class hint
            except Exception as exc:          # noqa: BLE001 — transient HTTP
                if attempt == max_retries - 1:
                    log.warning("Horizons failed for %d epochs after %d attempts: %s",
                                len(chunk), max_retries, exc)
                else:
                    time.sleep(backoff * (attempt + 1))

        if table is None:
            n_failed += 1
            continue
        frames.append(table.to_pandas())

    if n_failed:
        log.warning("%d/%d Horizons chunks failed; the ephemeris table is partial",
                    n_failed, len(chunks))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def query_ztf_metadata(ra, dec, width, height, jd_start, jd_end,
                       intersect=cfg.QueryConfig.intersect,
                       timeout=120.0, max_retries=3, backoff=2.0):
    """Search IRSA for ZTF science images overlapping a box and time window.

    Parameters
    ----------
    ra, dec : float
        Box centre in degrees.
    width, height : float
        Box size in degrees.
    jd_start, jd_end : float
        Inclusive ``obsjd`` bounds.
    intersect : str
        IRSA intersection mode; ``"overlaps"`` is the permissive default and is
        narrowed afterwards by the footprint containment test in
        :func:`search_frames`.

    Returns
    -------
    pandas.DataFrame
        Image metadata; **empty** when the search succeeded but matched nothing.

    Raises
    ------
    requests.HTTPError
        On a non-200 response after *max_retries* attempts.  The original code
        returned ``None`` for both "no results" and "request failed", which made
        coverage gaps indistinguishable from network trouble.
    """
    params = {
        "POS": f"{ra},{dec}",
        "SIZE": f"{width},{height}",
        "INTERSECT": intersect,
        "WHERE": f"obsjd BETWEEN {jd_start} AND {jd_end}",
        "RESPONSEFORMAT": "VOTABLE",
    }

    last_exc: Exception | None = None
    for attempt in range(max_retries):
        try:
            response = requests.get(cfg.IRSA_SEARCH_URL, params=params, timeout=timeout)
            response.raise_for_status()
            table = parse_votable(BytesIO(response.content))
            df = table.get_first_table().to_table(use_names_over_ids=True).to_pandas()
            if df.empty:
                return df
            df["filefracday"] = df["filefracday"].astype(str)
            return df
        except Exception as exc:                                # noqa: BLE001
            last_exc = exc
            if attempt < max_retries - 1:
                time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"IRSA metadata query failed after {max_retries} attempts") from last_exc


def _merge_on_jd(ztf: pd.DataFrame, eph: pd.DataFrame) -> pd.DataFrame:
    """Attach an ephemeris block to image metadata by matching JD, not position.

    Colliding column names (``airmass`` is supplied by both services) get an
    ``eph_`` prefix.  ``obsjd`` and ``datetime_jd`` are matched with
    :data:`_JD_MERGE_TOL` because Horizons rounds its output JDs.
    """
    if eph.empty:
        return ztf

    eph = eph.copy()
    overlap = (set(eph.columns) & set(ztf.columns)) - {"datetime_jd"}
    eph = eph.rename(columns={c: f"eph_{c}" for c in overlap})

    merged = pd.merge_asof(
        ztf.sort_values("obsjd"),
        eph.sort_values("datetime_jd"),
        left_on="obsjd", right_on="datetime_jd",
        direction="nearest", tolerance=_JD_MERGE_TOL,
    )
    unmatched = int(merged["datetime_jd"].isna().sum())
    if unmatched:
        log.warning("%d/%d frames had no ephemeris within %g d", unmatched, len(merged), _JD_MERGE_TOL)
    return merged.reset_index(drop=True)


def _window_midpoint(start_date, end_date):
    """Mid-window Julian Date, used to pick the orbit solution for the run."""
    try:
        from astropy.time import Time
        return float(Time([str(start_date), str(end_date)]).jd.mean())
    except Exception:                                           # noqa: BLE001
        return None


def search_frames(target, progress=True):
    """Find every ZTF science frame containing *target*.

    Runs the two-stage search: a coarse ephemeris to decide where and when to
    look, then a per-step IRSA query whose box is sized from the target's own
    apparent motion, then an exact-time ephemeris for each frame returned, then
    a footprint containment test.

    Parameters
    ----------
    target : ztfcomet.config.Target
        Target and its :class:`~ztfcomet.config.QueryConfig`.
    progress : bool
        Show a tqdm bar.

    Returns
    -------
    eph : pandas.DataFrame
        The coarse ephemeris after the ``rh``/``Vmag`` cuts.
    frames : pandas.DataFrame
        Image metadata joined to exact-time ephemerides, one row per frame
        whose footprint contains the target.
    report : QueryReport
        Counts of what was queried, kept, and dropped.

    Notes
    -----
    The search box is ``interval_days`` times the apparent rate in each axis, so
    a fast-moving target near opposition is searched over a correspondingly
    wider region.  This is the right behaviour and is preserved unchanged from
    the original implementation.
    """
    qc = target.query
    report = QueryReport()

    # Resolve the designation to a concrete orbit record for this window before
    # anything else. Passing a bare designation would leave the object class to
    # Horizons ("2P" -> Styx, a moon of Pluto) and the apparition to chance.
    mid_jd = _window_midpoint(target.start_date, target.end_date)
    target_id = target.resolve_orbit_record(mid_jd) if mid_jd else target.query_designation
    log.info("%s: querying Horizons as %r", target.name, target_id)

    eph = query_sso_ephemeris(
        target_id,
        epochs={"start": target.start_date, "stop": target.end_date,
                "step": f"{qc.interval_days:g}d"},
        quantities=cfg.EPHEM_QUANTITIES_COARSE,
        max_retries=qc.max_retries, backoff=qc.backoff,
        allow_fragment=target.allow_fragment,
    )
    if eph.empty:
        log.warning("No ephemeris returned for %s", target.name)
        return eph, pd.DataFrame(), report

    if "targetname" in eph.columns and len(eph):
        returned = str(eph["targetname"].iloc[0])
        ok, why = hz.verify_targetname(returned, target.query_designation,
                                       allow_fragment=target.allow_fragment)
        if not ok:
            raise ValueError(f"{target.name}: {why}")
        log.info("%s: Horizons returned %s", target.name, returned)

    report.n_eph_steps = len(eph)
    n0 = len(eph)
    eph = eph[eph["r"] < qc.rh_max]
    report.n_cut_rh = n0 - len(eph)
    if "Tmag" in eph.columns:
        n1 = len(eph)
        eph = eph[eph["Tmag"] < qc.vmag_max]
        report.n_cut_vmag = n1 - len(eph)
    eph = eph.reset_index(drop=True)
    report.n_eph_kept = len(eph)

    collected = []
    rows = eph.iterrows()
    if progress:
        rows = tqdm(rows, total=len(eph), desc=f"{target.name}: IRSA search")

    for _, step in rows:
        report.n_steps_queried += 1
        try:
            # Apparent rates are arcsec/hr; convert to deg/day and span the step.
            ra_dot = step["RA_rate"] / 3600 * 24
            dec_dot = step["DEC_rate"] / 3600 * 24
            width = abs(qc.interval_days * ra_dot)
            height = abs(qc.interval_days * dec_dot)

            ztf = query_ztf_metadata(
                ra=step["RA"], dec=step["DEC"], width=width, height=height,
                jd_start=step["datetime_jd"] - 0.5 * qc.interval_days,
                jd_end=step["datetime_jd"] + 0.5 * qc.interval_days,
                intersect=qc.intersect, timeout=qc.timeout,
                max_retries=qc.max_retries, backoff=qc.backoff,
            )
        except Exception as exc:                                # noqa: BLE001
            report.n_steps_failed += 1
            report.failures.append(f"JD {step['datetime_jd']:.5f}: {exc}")
            log.warning("IRSA query failed at JD %.5f: %s", step["datetime_jd"], exc)
            continue

        if ztf.empty:
            report.n_steps_no_coverage += 1
            continue

        report.n_frames_found += len(ztf)
        ztf = ztf.sort_values("obsjd").reset_index(drop=True)

        eph_exact = query_sso_ephemeris(
            target_id, epochs=list(ztf["obsjd"]),
            quantities=cfg.EPHEM_QUANTITIES_FULL,
            max_epochs_per_call=qc.max_epochs_per_call,
            max_retries=qc.max_retries, backoff=qc.backoff,
        )
        merged = _merge_on_jd(ztf, eph_exact)
        if "RA" not in merged.columns:
            report.n_steps_failed += 1
            report.failures.append(f"JD {step['datetime_jd']:.5f}: no exact-time ephemeris")
            continue

        # Footprint containment: the IRSA "overlaps" mode returns frames whose
        # corners merely intersect the search box, so test the target explicitly.
        corners_ra = merged[["ra1", "ra2", "ra3", "ra4"]]
        corners_dec = merged[["dec1", "dec2", "dec3", "dec4"]]
        merged["ra_min"], merged["ra_max"] = corners_ra.min(axis=1), corners_ra.max(axis=1)
        merged["dec_min"], merged["dec_max"] = corners_dec.min(axis=1), corners_dec.max(axis=1)

        inside = merged[
            merged["RA"].between(merged["ra_min"], merged["ra_max"])
            & merged["DEC"].between(merged["dec_min"], merged["dec_max"])
        ]
        if not inside.empty:
            collected.append(inside)

    frames = (pd.concat(collected, ignore_index=True) if collected else pd.DataFrame())
    if not frames.empty:
        frames = (frames
                  .drop_duplicates(subset=["filefracday", "field", "ccdid", "qid", "filtercode"])
                  .sort_values("obsjd")
                  .reset_index(drop=True))
    report.n_frames_containing_target = len(frames)

    if report.n_steps_failed:
        log.warning("%d/%d steps FAILED — coverage is incomplete; see report.failures",
                    report.n_steps_failed, report.n_steps_queried)
    return eph, frames, report
