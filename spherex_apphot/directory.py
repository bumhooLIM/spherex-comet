"""
Filesystem layout for the SPHEREx comet aperture-photometry pipeline.

Every path used anywhere in the package is resolved here, exactly once, so that
there is a single place to look when a run cannot find its inputs.  Each path
may be overridden with an environment variable, which is what makes the package
usable both against the production dataset on the external SSD and against the
small ``data/spherex_sample/`` tree checked into the repository.

Since 2026-09-15 the package is the second stage of the ``spherex-comet`` project
(``ztfcomet`` -> ``spherex_apphot`` -> ``spherex_comspec``); its outputs are namespaced
as ``results/apphot/`` and ``fig/apphot/``, and ``results/apphot/photometry/`` is what
``spherex_comspec`` reads.

Environment overrides
---------------------
``SPHEREX_T7_DIR``      root of the production dataset (default ``/Volumes/T7/data/spherex-comet``)
``SPHEREX_DB``          path to the cutout-index Parquet file
``SPHEREX_FITS_DIR``    directory (or ``os.pathsep``-separated list) holding the cutout FITS files
``SPHEREX_FITS_FILTERED_DIR``  per-target tree written by ``extract_targets.py``
``SPHEREX_GAIA_DIR``    directory holding the Gaia DR3 ``.npy`` catalogue and its cache
``SPHEREX_RESULT_DIR``  root of all pipeline outputs (default ``<project>/results/apphot``)
``SPHEREX_FIG_DIR``     root of all figures (default ``<project>/fig/apphot``)

Notes
-----
The primitive pipeline referenced ``directory.APPHOT_DIR`` and
``directory.COMBFITS_DIR``, neither of which existed, so it raised
``AttributeError`` before doing any work (review item B1).  Those names are
defined here, and :func:`ensure_dirs` creates every output directory so the
first ``savefig``/``to_csv`` cannot fail on a missing parent (review item B3).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List

__all__ = [
    "WORK_DIR", "T7_DIR", "DB_DIR", "DB_PATH", "FITS_DIR", "FITS_ROOTS",
    "FITS_FILTERED_DIR", "REFCAT_DIR", "GAIA_ALL_NPY", "GAIA_CACHE_DIR",
    "DATA_SAMPLE_DIR", "RESULT_DIR", "APPHOT_DIR", "COMBFITS_DIR", "LOG_DIR",
    "STATUS_CSV", "FIG_DIR", "DOC_DIR", "NOTEBOOK_DIR", "TARGET_LIST",
    "ensure_dirs", "use_sample_data", "describe",
]


def _env_path(key: str, default: Path) -> Path:
    val = os.environ.get(key)
    return Path(val).expanduser() if val else default


def _env_paths(key: str, default: List[Path]) -> List[Path]:
    val = os.environ.get(key)
    if not val:
        return list(default)
    return [Path(p).expanduser() for p in val.split(os.pathsep) if p]


# --- repository ------------------------------------------------------------
# directory.py lives at <project>/spherex_apphot/directory.py
WORK_DIR: Path = Path(__file__).resolve().parents[1]
#: Notes of this stage (the code review and the upgrade record).
DOC_DIR: Path = WORK_DIR / "doc" / "apphot"
#: Holds the shared ``rcparams.py`` (figure style) and the validation notebooks.
NOTEBOOK_DIR: Path = WORK_DIR / "notebooks"
#: Sample cutouts of 2P and 24P plus their ``db_filtered.parq`` slice, for the
#: end-to-end test when the SSD is not mounted.
DATA_SAMPLE_DIR: Path = WORK_DIR / "data" / "spherex_sample"
#: The 68-comet working list shared by the three stages.
TARGET_LIST: Path = WORK_DIR / "data" / "reference" / "sx_comet_list_ver2607.xlsx"

# --- production dataset (external SSD) -------------------------------------
T7_DIR: Path = _env_path("SPHEREX_T7_DIR", Path("/Volumes/T7/data/spherex-comet"))
DB_DIR: Path = T7_DIR / "comets_v5"
DB_PATH: Path = _env_path("SPHEREX_DB", DB_DIR / "db_filtered.parq")
FITS_DIR: Path = _env_path("SPHEREX_FITS_DIR", T7_DIR / "comets_v5_fits")

#: Per-target tree written by ``extract_targets.py``.  The flat ``comets_v5_fits``
#: holds ~139 000 files, and a directory that large is slow to stat on an
#: external drive; the filtered tree holds only the targets in the working list,
#: split one directory per comet.
FITS_FILTERED_DIR: Path = _env_path("SPHEREX_FITS_FILTERED_DIR",
                                    T7_DIR / "comets_v5_fits_filtered")

#: Roots searched for a cutout FITS file, in order.  Both the flat production
#: layout (``<root>/<filename>``) and the per-target layout
#: (``<root>/<objdesig>/<filename>``) are probed; see ``fitsio.FitsResolver``.
#: The filtered tree comes first so it is preferred when present.
FITS_ROOTS: List[Path] = _env_paths(
    "SPHEREX_FITS_DIR", [FITS_FILTERED_DIR, FITS_DIR, DATA_SAMPLE_DIR])

# --- reference catalogue ---------------------------------------------------
REFCAT_DIR: Path = _env_path("SPHEREX_GAIA_DIR", Path.home() / "Desktop" / "data" / "gaia_dr3")
GAIA_ALL_NPY: Path = REFCAT_DIR / "gaiadr3_all.npy"
#: Declination-sorted cache written by ``build_gaia_cache.py``; optional but
#: turns each per-target query from an 11 GB scan into a few-MB slice.
GAIA_CACHE_DIR: Path = REFCAT_DIR / "gaiadr3_deccache"

# --- outputs ---------------------------------------------------------------
# results/apphot/photometry/<slug>.csv is the product every later stage reads;
# spherex_comspec.directory.APPHOT_DIR points at the same directory.
RESULT_DIR: Path = _env_path("SPHEREX_RESULT_DIR", WORK_DIR / "results" / "apphot")
APPHOT_DIR: Path = RESULT_DIR / "photometry"
COMBFITS_DIR: Path = RESULT_DIR / "stacks"
LOG_DIR: Path = RESULT_DIR / "logs"
STATUS_CSV: Path = RESULT_DIR / "status.csv"
FIG_DIR: Path = _env_path("SPHEREX_FIG_DIR", WORK_DIR / "fig" / "apphot")

_OUTPUT_DIRS = (RESULT_DIR, APPHOT_DIR, COMBFITS_DIR, LOG_DIR, FIG_DIR)


def ensure_dirs() -> None:
    """Create every output directory.  Idempotent; safe to call repeatedly."""
    for d in _OUTPUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)


def use_sample_data() -> None:
    """
    Repoint the module at the in-repo ``data/spherex_sample/`` tree.

    Intended for the end-to-end test and for the notebook when the external SSD
    is not mounted.  Mutates the module-level globals in place so that modules
    which did ``from . import directory`` see the change.
    """
    global DB_PATH, FITS_ROOTS
    DB_PATH = DATA_SAMPLE_DIR / "db_filtered.parq"
    FITS_ROOTS = [DATA_SAMPLE_DIR]


def describe() -> str:
    """Return a human-readable summary of the resolved paths and their status."""
    rows = [
        ("WORK_DIR", WORK_DIR), ("DB_PATH", DB_PATH),
        ("FITS_FILTERED_DIR", FITS_FILTERED_DIR),
        ("FITS_ROOTS", os.pathsep.join(str(p) for p in FITS_ROOTS)),
        ("GAIA_ALL_NPY", GAIA_ALL_NPY), ("GAIA_CACHE_DIR", GAIA_CACHE_DIR),
        ("APPHOT_DIR", APPHOT_DIR), ("COMBFITS_DIR", COMBFITS_DIR),
        ("LOG_DIR", LOG_DIR), ("FIG_DIR", FIG_DIR),
    ]
    out = []
    for name, val in rows:
        if isinstance(val, Path):
            tag = "ok     " if val.exists() else "MISSING"
        else:
            tag = "       "
        out.append(f"  {name:<14s} [{tag}] {val}")
    return "\n".join(out)
