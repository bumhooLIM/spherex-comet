"""
Filesystem layout for ``spherex_comspec``.

Every path the package reads or writes is resolved here.  Inputs come from the
catalog's ``data/apphot_revised/`` (the output of the ``spherex_apphot``
photometry pipeline); outputs never touch that directory.  The primitive
notebooks wrote a ``phase_update`` column *into* the photometry CSVs, which the
project handoff warns is unsafe; this package writes its group assignment to a
separate file instead.

Environment overrides
---------------------
``COMSPEC_ROOT``       repository root (default: two levels above this file)
``COMSPEC_APPHOT_DIR`` directory of ``<target>.csv`` photometry tables
``COMSPEC_DATA_DIR``   root for intermediate products (default ``data/comspec``)
``COMSPEC_RESULT_DIR`` root for results (default ``results/comspec``)
``COMSPEC_FIG_DIR``    root for figures (default ``fig/comspec``)
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "ROOT", "APPHOT_DIR", "LEGACY_APPHOT_DIR", "DATA_DIR", "RESULT_DIR", "FIG_DIR",
    "LOG_DIR", "NOTEBOOK_DIR", "DOC_DIR", "ORBIT_CLASSES_CSV",
    "variant_dirs", "ensure_dirs", "describe",
]


def _env(key: str, default: Path) -> Path:
    v = os.environ.get(key)
    return Path(v).expanduser() if v else default


ROOT: Path = _env("COMSPEC_ROOT", Path(__file__).resolve().parents[2])

#: Revised aperture photometry -- one CSV per target, written by ``spherex_apphot``.
APPHOT_DIR: Path = _env("COMSPEC_APPHOT_DIR", ROOT / "data" / "apphot_revised")
#: The pre-revision tables, kept only so the two groupings can be compared.
LEGACY_APPHOT_DIR: Path = ROOT / "data" / "apphot"

DATA_DIR: Path = _env("COMSPEC_DATA_DIR", ROOT / "data" / "comspec")
RESULT_DIR: Path = _env("COMSPEC_RESULT_DIR", ROOT / "results" / "comspec")
FIG_DIR: Path = _env("COMSPEC_FIG_DIR", ROOT / "fig" / "comspec")
LOG_DIR: Path = RESULT_DIR / "logs"
NOTEBOOK_DIR: Path = ROOT / "notebooks"
DOC_DIR: Path = ROOT / "doc"
ORBIT_CLASSES_CSV: Path = ROOT / "data" / "comet_orbit_classes.csv"


def variant_dirs(name: str) -> dict:
    """
    Output directories for one pipeline variant.

    A variant is one (flag policy, flux space) combination; keeping each one in
    its own tree is what lets the contamination and distance-correction studies
    compare complete runs rather than overwritten ones.
    """
    return {
        "emission": DATA_DIR / name / "emission",
        "results": RESULT_DIR / name,
        "lines": RESULT_DIR / name / "gas_fit_lines",
        "fig": FIG_DIR / name,
        "fig_cont": FIG_DIR / name / "cont_subtract",
        "fig_fit": FIG_DIR / name / "emission_model",
    }


def ensure_dirs(variant: str | None = None) -> None:
    """Create the shared output directories, and a variant's if one is named."""
    for d in (DATA_DIR, RESULT_DIR, FIG_DIR, LOG_DIR, FIG_DIR / "phase_group"):
        d.mkdir(parents=True, exist_ok=True)
    if variant:
        for d in variant_dirs(variant).values():
            d.mkdir(parents=True, exist_ok=True)


def describe() -> str:
    rows = [("ROOT", ROOT), ("APPHOT_DIR", APPHOT_DIR), ("DATA_DIR", DATA_DIR),
            ("RESULT_DIR", RESULT_DIR), ("FIG_DIR", FIG_DIR), ("LOG_DIR", LOG_DIR)]
    return "\n".join(f"  {k:<12s} [{'ok     ' if v.exists() else 'MISSING'}] {v}" for k, v in rows)
