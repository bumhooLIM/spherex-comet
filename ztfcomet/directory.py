"""Single source of truth for every path used by :mod:`ztfcomet`.

Nothing else in the project may hardcode a path or derive one from
``Path.cwd()``.  Notebooks that did so (``WORKDIR = Path.cwd() / ".."``) broke
whenever they were executed from anywhere but ``notebooks/``; this module
locates the project by walking up for the ``pyproject.toml`` marker instead, so
imports behave identically from a notebook, a script, or a test.

Raw FITS live on an external SSD, not in the repository.  ``DATA_ROOT``
therefore resolves in this order:

1. ``$ZTFCOMET_DATA`` if set — explicit override, always wins.
2. The SSD path in :data:`SSD_DATA_ROOT` if that volume is mounted.
3. ``<project>/data`` — the fallback, so a fresh clone works with no SSD.

Layout under each root is ``<root>/<target_slug>/``, matching what is already
on the SSD.

Examples
--------
>>> from ztfcomet import directory as d
>>> d.data_dir("24P")            # doctest: +SKIP
PosixPath('/Volumes/T7/data/ztf-comet/24P')
>>> d.target_slug("2019 Y3")
'2019Y3'
"""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

log = logging.getLogger(__name__)

__all__ = [
    "using_fallback_root",
    "PROJECT_ROOT", "DATA_ROOT", "RESULT_ROOT", "FIG_ROOT", "DOC_ROOT",
    "SSD_DATA_ROOT", "GAIA_ROOT", "target_slug", "data_dir", "result_dir",
    "fig_dir", "describe",
]

#: Where raw FITS live when the external SSD is mounted.
SSD_DATA_ROOT = Path("/Volumes/T7/data/ztf-comet")

#: Default home of the local Gaia DR3 catalogue used for contamination
#: flagging.  Override with ``$ZTFCOMET_GAIA``.  Holds ``gaiadr3_all.npy`` and,
#: ideally, the ``gaiadr3_deccache/`` fast path.
DEFAULT_GAIA_ROOT = Path.home() / "Desktop" / "data" / "gaia_dr3"

_MARKER = "pyproject.toml"


def _find_project_root(start: Path | None = None) -> Path:
    """Walk up from *start* until the directory containing ``pyproject.toml``.

    Falls back to the package's own parent so that an un-installed, un-marked
    checkout still resolves to something sensible rather than raising.
    """
    start = (start or Path(__file__)).resolve()
    for candidate in (start, *start.parents):
        if (candidate / _MARKER).is_file():
            return candidate
    return Path(__file__).resolve().parent.parent


PROJECT_ROOT: Path = _find_project_root()


def _probe(path: Path, attempts: int = 3, delay: float = 1.0) -> bool:
    """Is *path* a readable directory, retrying a sleeping external volume?

    An external SSD that has spun down can fail its first ``is_dir()`` and
    answer normally a second later.  Taking the first answer as final sent a
    resumed survey run to the local fallback, where it found no status file and
    silently restarted the whole thing in the wrong place.
    """
    for attempt in range(attempts):
        try:
            if path.is_dir():
                return True
        except OSError:
            pass
        if attempt < attempts - 1:
            time.sleep(delay)
    return False


def _resolve_data_root() -> Path:
    env = os.environ.get("ZTFCOMET_DATA")
    if env:
        return Path(env).expanduser()
    if _probe(SSD_DATA_ROOT):
        return SSD_DATA_ROOT
    log.warning("SSD data root %s is not available; falling back to %s",
                SSD_DATA_ROOT, PROJECT_ROOT / "data")
    return PROJECT_ROOT / "data"


#: Raw FITS cutouts, ``eph.csv``, ``ztf.csv``, ``fits_urls.txt``.  Large; never committed.
DATA_ROOT: Path = _resolve_data_root()

#: Photometry tables and other clean, small outputs.  Never committed.
RESULT_ROOT: Path = PROJECT_ROOT / "results"

#: Figures.  Never committed.
FIG_ROOT: Path = PROJECT_ROOT / "fig"

#: Technical guidebooks and review documents.  Committed.
DOC_ROOT: Path = PROJECT_ROOT / "doc"


def _resolve_gaia_root() -> Path:
    env = os.environ.get("ZTFCOMET_GAIA")
    return Path(env).expanduser() if env else DEFAULT_GAIA_ROOT


#: Local Gaia DR3 catalogue.  Large and read-only; never inside the repository.
GAIA_ROOT: Path = _resolve_gaia_root()


def target_slug(targetname: str) -> str:
    """Filesystem-safe directory name for a target designation.

    Strips whitespace and replaces path separators, so ``"2019 Y3"`` and
    ``"C/2024 E1"`` become ``"2019Y3"`` and ``"C2024E1"``.  Matches the naming
    already used on the SSD.
    """
    slug = "".join(str(targetname).split())
    return slug.replace("/", "").replace("\\", "")


def _sub(root: Path, targetname: str | None, create: bool) -> Path:
    path = root if targetname is None else root / target_slug(targetname)
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def data_dir(targetname: str | None = None, create: bool = True) -> Path:
    """Directory holding raw FITS and query products for *targetname*."""
    return _sub(DATA_ROOT, targetname, create)


SUBJECTS = ("photometry", "profile", "afrho")


def _subject(subject: str) -> str:
    if subject not in SUBJECTS:
        raise ValueError(
            f"unknown output subject {subject!r}: results and figures are organised by "
            f"subject {SUBJECTS}, not by target.  Use result_path()/fig_path() with "
            "targetname= for a per-target file, or fig_dir(subject, targetname) for bulk images.")
    return subject


def result_dir(subject: str, create: bool = True) -> Path:
    """``results/<subject>/`` -- outputs are grouped by subject, not by target.

    A target name here is an error on purpose.  The tree used to be
    per-target, and a stale call site should fail rather than quietly rebuild
    that layout beside the new one.
    """
    return _sub(RESULT_ROOT / _subject(subject), None, create)


def fig_dir(subject: str, targetname: str | None = None, create: bool = True) -> Path:
    """``fig/<subject>/``, or ``fig/<subject>/<target>/`` for bulk per-frame images."""
    return _sub(FIG_ROOT / _subject(subject), targetname, create)


def result_path(subject: str, filename: str, targetname: str | None = None,
                create: bool = True) -> Path:
    """``results/<subject>/<target>_<filename>``, or ``<filename>`` for a survey-level table."""
    name = f"{target_slug(targetname)}_{filename}" if targetname else filename
    return result_dir(subject, create) / name


def fig_path(subject: str, filename: str, targetname: str | None = None,
             create: bool = True) -> Path:
    """``fig/<subject>/<target>_<filename>``, or ``<filename>`` for a survey-level figure."""
    name = f"{target_slug(targetname)}_{filename}" if targetname else filename
    return fig_dir(subject, None, create) / name


def photometry_path(targetname: str, create: bool = True) -> Path:
    """``results/photometry/<target>.csv`` -- the one table everything downstream reads."""
    return result_dir("photometry", create) / f"{target_slug(targetname)}.csv"


def profile_paths(targetname: str, create: bool = True) -> dict:
    """The three profile tables of a target, keyed ``profile``, ``summary``, ``stars``."""
    return {k: result_path("profile", f"{k}.csv", targetname, create)
            for k in ("profile", "summary", "stars")}


def describe() -> str:
    """Human-readable summary of the resolved roots, for notebook sanity checks."""
    src = (
        "$ZTFCOMET_DATA" if os.environ.get("ZTFCOMET_DATA")
        else "SSD" if DATA_ROOT == SSD_DATA_ROOT
        else "project fallback"
    )
    return "\n".join([
        f"PROJECT_ROOT : {PROJECT_ROOT}",
        f"DATA_ROOT    : {DATA_ROOT}   [{src}]{'' if DATA_ROOT.is_dir() else '  (MISSING)'}",
        f"RESULT_ROOT  : {RESULT_ROOT}",
        f"FIG_ROOT     : {FIG_ROOT}",
        f"DOC_ROOT     : {DOC_ROOT}",
        f"GAIA_ROOT    : {GAIA_ROOT}"
        f"{'' if GAIA_ROOT.is_dir() else '   (MISSING — contamination flagging disabled)'}",
    ])


def using_fallback_root() -> bool:
    """True when DATA_ROOT is the in-project fallback rather than the SSD.

    A long unattended run should refuse to start in that state: the SSD holds
    both the existing data and the resume checkpoint, so silently using the
    fallback means redoing everything into the wrong place.
    """
    return (not os.environ.get("ZTFCOMET_DATA")
            and DATA_ROOT == PROJECT_ROOT / "data")
