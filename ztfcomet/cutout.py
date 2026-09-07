"""Construct IRSA cutout URLs and download them with validation.

The original ``import_fitsurl`` wrote whatever came back whenever the status was
200.  IRSA answers some cutout requests with 200 and an HTML error body, so
245-byte ``<title>404 Not Found</title>`` documents ended up on disk named
``*_sciimg.fits`` — 2 of 113 in the ``2024E1`` sample.  Worse, the resume check
``if save_path.exists(): continue`` meant no later run could ever repair them:
they were permanently poisoned, and ``ccdproc`` silently dropped them from the
analysis without a word.

This module fixes that by (a) checking the payload really is FITS before
accepting it, (b) writing to a ``.part`` file and renaming only on success, so
an interrupted download can never masquerade as a complete one, and (c)
recording per-URL status in a manifest so failures are visible and retryable.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import requests
from tqdm.auto import tqdm

from . import config as cfg

__all__ = ["construct_fitsurl", "build_urls", "save_urls", "download_urls",
           "choose_cutout_size", "verify_downloads", "DownloadReport"]

log = logging.getLogger(__name__)

#: Every valid FITS file starts with this.  Cheapest possible payload check.
_FITS_MAGIC = b"SIMPLE"

#: Nothing smaller than one FITS header block can be a real image.
_MIN_FITS_BYTES = 2880


@dataclass
class DownloadReport:
    """Per-run download accounting."""

    n_total: int = 0
    n_downloaded: int = 0
    n_skipped_existing: int = 0
    n_failed: int = 0
    n_rejected_not_fits: int = 0
    failures: list[tuple[str, str]] = field(default_factory=list)

    def __str__(self) -> str:
        return "\n".join([
            f"URLs                 : {self.n_total}",
            f"  downloaded         : {self.n_downloaded}",
            f"  already on disk    : {self.n_skipped_existing}",
            f"  rejected (no FITS) : {self.n_rejected_not_fits}",
            f"  failed             : {self.n_failed}",
        ])


def construct_fitsurl(row, is_cutout=True, cutout_size="10arcmin",
                      ra_col="RA", dec_col="DEC"):
    """Build the IRSA URL for one ZTF science frame.

    Parameters
    ----------
    row : pandas.Series
        One row of IRSA image metadata; needs ``filefracday``, ``field``,
        ``filtercode``, ``ccdid``, ``imgtypecode``, ``qid``, and — for cutouts —
        the ephemeris position in *ra_col* / *dec_col*.
    is_cutout : bool
        Append the ``?center=&size=`` cutout query.  When ``False`` the full
        quadrant image (~37 MB) is requested.
    cutout_size : str
        IRSA size string, e.g. ``"10arcmin"``.

    Returns
    -------
    str
        Fully-qualified download URL.
    """
    ffd = str(row["filefracday"])
    year, monthday, fracday = ffd[:4], ffd[4:8], ffd[8:]

    name = (f"ztf_{ffd}_{int(row['field']):06d}_{row['filtercode']}"
            f"_c{int(row['ccdid']):02d}_{row['imgtypecode']}_q{int(row['qid'])}_sciimg.fits")
    url = f"{cfg.IRSA_DATA_URL}/{year}/{monthday}/{fracday}/{name}"

    if is_cutout:
        url += f"?center={row[ra_col]},{row[dec_col]}&size={cutout_size}&gzip=false"
    return url


def frame_filename(row) -> str:
    """Local filename for a frame, matching the IRSA basename."""
    return construct_fitsurl(row, is_cutout=False).rsplit("/", 1)[-1]


def build_urls(frames, is_cutout=True, cutout_size="10arcmin"):
    """Vectorised :func:`construct_fitsurl` over a metadata table.

    Returns
    -------
    pandas.DataFrame
        ``frames`` with ``fits_url`` and ``file`` columns added.

    Notes
    -----
    Cutouts of the same base image at different centres collapse to the same
    local filename.  Duplicates are dropped here rather than left to overwrite
    each other at download time.
    """
    if frames is None or frames.empty:
        return pd.DataFrame(columns=["fits_url", "file"])

    out = frames.copy()
    out["fits_url"] = out.apply(
        construct_fitsurl, axis=1, is_cutout=is_cutout, cutout_size=cutout_size)
    out["file"] = out.apply(frame_filename, axis=1)

    dupes = int(out["file"].duplicated().sum())
    if dupes:
        log.warning("%d frames share a local filename; keeping the first of each", dupes)
        out = out.drop_duplicates(subset="file").reset_index(drop=True)
    return out


def save_urls(frames, path):
    """Write the ``fits_url`` column to a newline-delimited text file."""
    path = Path(path)
    urls = list(frames["fits_url"]) if "fits_url" in frames else []
    path.write_text("\n".join(urls) + ("\n" if urls else ""))
    log.info("Wrote %d URLs to %s", len(urls), path)
    return path


def _looks_like_fits(path: Path) -> bool:
    """True if *path* is plausibly a FITS file rather than an error page."""
    try:
        if path.stat().st_size < _MIN_FITS_BYTES:
            return False
        with path.open("rb") as fh:
            return fh.read(len(_FITS_MAGIC)) == _FITS_MAGIC
    except OSError:
        return False


def download_urls(urls, outdir, overwrite=False, repair=True, progress=True,
                  timeout=120.0, max_retries=3, backoff=2.0, manifest=True):
    """Download FITS files, rejecting anything that is not really FITS.

    Parameters
    ----------
    urls : iterable of str or pandas.DataFrame
        URLs, or a table with a ``fits_url`` column.
    outdir : path-like
        Destination directory; created if absent.
    overwrite : bool
        Re-download files that already exist and are valid.
    repair : bool
        Re-download existing files that fail the FITS check.  This is what makes
        a corrupt earlier run recoverable; leave it on.
    manifest : bool
        Write ``download_manifest.csv`` recording per-URL status.

    Returns
    -------
    DownloadReport

    Notes
    -----
    The payload is written to ``<name>.part`` and renamed only after passing the
    FITS check, so an interrupted transfer leaves no file that a later run would
    mistake for complete.
    """
    if isinstance(urls, pd.DataFrame):
        urls = list(urls["fits_url"])
    urls = [u for u in (str(u).strip() for u in urls) if u]

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    report = DownloadReport(n_total=len(urls))
    records = []

    iterator = tqdm(urls, desc="Downloading FITS") if progress else urls
    for url in iterator:
        name = url.split("?")[0].rsplit("/", 1)[-1]
        dest = outdir / name
        part = dest.with_suffix(dest.suffix + ".part")

        if dest.exists() and not overwrite:
            if _looks_like_fits(dest):
                report.n_skipped_existing += 1
                records.append((name, url, "existing", ""))
                continue
            if not repair:
                report.n_rejected_not_fits += 1
                records.append((name, url, "corrupt-kept", "failed FITS check"))
                continue
            log.info("Repairing corrupt file %s", name)

        status, detail = "failed", ""
        for attempt in range(max_retries):
            try:
                with requests.get(url, stream=True, timeout=timeout) as response:
                    response.raise_for_status()
                    ctype = response.headers.get("Content-Type", "")
                    if "html" in ctype.lower() or "xml" in ctype.lower():
                        # IRSA answers some bad cutout requests 200 + HTML.
                        raise ValueError(f"server returned {ctype!r}, not FITS")
                    with part.open("wb") as fh:
                        for chunk in response.iter_content(chunk_size=1 << 16):
                            fh.write(chunk)

                if not _looks_like_fits(part):
                    raise ValueError(f"payload is not FITS ({part.stat().st_size} bytes)")

                part.replace(dest)
                status, detail = "downloaded", ""
                report.n_downloaded += 1
                break
            except Exception as exc:                            # noqa: BLE001
                detail = str(exc)
                part.unlink(missing_ok=True)
                if attempt < max_retries - 1:
                    time.sleep(backoff * (attempt + 1))

        if status == "failed":
            if "not FITS" in detail or "not, FITS" in detail or "SIMPLE" in detail:
                report.n_rejected_not_fits += 1
            else:
                report.n_failed += 1
            report.failures.append((url, detail))
            log.warning("Download failed: %s (%s)", name, detail)

        records.append((name, url, status, detail))

    if manifest and records:
        pd.DataFrame(records, columns=["file", "url", "status", "detail"]).to_csv(
            outdir / "download_manifest.csv", index=False)

    if report.n_failed or report.n_rejected_not_fits:
        log.warning("%d downloads did not yield valid FITS; rerun to retry",
                    report.n_failed + report.n_rejected_not_fits)
    return report


def choose_cutout_size(frames, query_config=None):
    """Pick a cutout size from the target's own ephemeris.

    Angular coma size scales as ``1/delta``, and a bright comet has a larger
    detectable coma, so a fixed box is either wasteful for the faint distant
    majority or too small for the few bright or nearby ones.

    The large size is used when **either** condition holds at any epoch:

    * the predicted total magnitude reaches ``bright_vmag`` or brighter, or
    * the observer distance ``delta`` falls below ``close_delta_au``.

    Parameters
    ----------
    frames : pandas.DataFrame
        Frame table carrying ``Tmag`` and ``delta`` (from the exact-time
        ephemeris).  Missing columns simply do not trigger their condition.
    query_config : ztfcomet.config.QueryConfig, optional

    Returns
    -------
    size : str
        The IRSA size string to request.
    reason : str
        Human-readable justification, for the run log.

    Notes
    -----
    ``delta`` is the *observer* distance, not the heliocentric distance: it is
    what sets the angular size of the coma on the sky.
    """
    qc = query_config or cfg.QueryConfig()
    if not qc.adaptive_cutout:
        return qc.cutout_size, "adaptive sizing off"
    if frames is None or len(frames) == 0:
        return qc.cutout_size_small, "no frames"

    def _min(column):
        if column not in frames:
            return None
        values = pd.to_numeric(frames[column], errors="coerce").dropna()
        return float(values.min()) if len(values) else None

    vmag_min = _min("Tmag")
    if vmag_min is None:
        vmag_min = _min("tmag")
    delta_min = _min("delta")

    triggers = []
    if vmag_min is not None and vmag_min < qc.bright_vmag:
        triggers.append(f"Vmag_min={vmag_min:.2f} < {qc.bright_vmag:g}")
    if delta_min is not None and delta_min < qc.close_delta_au:
        triggers.append(f"delta_min={delta_min:.3f} au < {qc.close_delta_au:g}")

    parts = []
    if vmag_min is not None:
        parts.append(f"Vmag_min={vmag_min:.2f}")
    if delta_min is not None:
        parts.append(f"delta_min={delta_min:.3f} au")
    summary = ", ".join(parts) if parts else "no Vmag/delta available"

    if triggers:
        return qc.cutout_size_large, f"{summary} -> large ({'; '.join(triggers)})"
    return qc.cutout_size_small, f"{summary} -> small"


def verify_downloads(urls, outdir):
    """Check that every requested URL produced a valid FITS file on disk.

    Run after :func:`download_urls`.  The download step already rejects
    non-FITS payloads, but this re-reads what is actually on disk, so a file
    truncated by a full volume or removed afterwards is still caught.

    Parameters
    ----------
    urls : iterable of str or pandas.DataFrame
        The URL list that was requested.
    outdir : path-like

    Returns
    -------
    dict
        ``n_expected``, ``n_present``, ``n_valid``, ``n_missing``,
        ``n_corrupt``, ``missing`` and ``corrupt`` (file-name lists), and
        ``complete`` (bool).
    """
    if isinstance(urls, pd.DataFrame):
        urls = list(urls["fits_url"]) if "fits_url" in urls else []
    outdir = Path(outdir)

    expected = [str(u).split("?")[0].rsplit("/", 1)[-1] for u in urls if str(u).strip()]
    missing, corrupt, valid = [], [], 0
    for name in expected:
        path = outdir / name
        if not path.exists():
            missing.append(name)
        elif not _looks_like_fits(path):
            corrupt.append(name)
        else:
            valid += 1

    return {
        "n_expected": len(expected),
        "n_present": len(expected) - len(missing),
        "n_valid": valid,
        "n_missing": len(missing),
        "n_corrupt": len(corrupt),
        "missing": missing,
        "corrupt": corrupt,
        "complete": not missing and not corrupt,
    }
