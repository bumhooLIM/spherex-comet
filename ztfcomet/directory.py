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

import os
from pathlib import Path

__all__ = [
    "PROJECT_ROOT", "DATA_ROOT", "RESULT_ROOT", "FIG_ROOT", "DOC_ROOT",
    "SSD_DATA_ROOT", "target_slug", "data_dir", "result_dir", "fig_dir",
    "describe",
]

#: Where raw FITS live when the external SSD is mounted.
SSD_DATA_ROOT = Path("/Volumes/T7/data/ztf-comet")

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


def _resolve_data_root() -> Path:
    env = os.environ.get("ZTFCOMET_DATA")
    if env:
        return Path(env).expanduser()
    if SSD_DATA_ROOT.is_dir():
        return SSD_DATA_ROOT
    return PROJECT_ROOT / "data"


#: Raw FITS cutouts, ``eph.csv``, ``ztf.csv``, ``fits_urls.txt``.  Large; never committed.
DATA_ROOT: Path = _resolve_data_root()

#: Photometry tables and other clean, small outputs.  Never committed.
RESULT_ROOT: Path = PROJECT_ROOT / "results"

#: Figures.  Never committed.
FIG_ROOT: Path = PROJECT_ROOT / "fig"

#: Technical guidebooks and review documents.  Committed.
DOC_ROOT: Path = PROJECT_ROOT / "doc"


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


def result_dir(targetname: str | None = None, create: bool = True) -> Path:
    """Directory holding photometry tables for *targetname*."""
    return _sub(RESULT_ROOT, targetname, create)


def fig_dir(targetname: str | None = None, create: bool = True) -> Path:
    """Directory holding figures for *targetname*."""
    return _sub(FIG_ROOT, targetname, create)


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
    ])
