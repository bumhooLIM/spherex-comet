"""
Filesystem layout for ``spherex_comspec``.

Every path the package reads or writes is resolved here.  The layout is flat on
purpose -- the main result must be readable in one glance -- and every product
of a *study* variant (a flag-policy, badphot, distance-correction or previous
baseline run) is kept one level down under ``studies/`` so it can never be
mistaken for the catalog result::

    data/apphot/                 revised aperture photometry (input, never written to)
    data/phase_assignment.csv    per-exposure phase labels, shared by every variant
    data/emission/               continuum summaries + point spectra of the main variant
    data/studies/<variant>/emission/
    results/gas_fit.csv          THE production rates (main variant), gas_fit_lines/, run.meta.json,
                                 continuum_summary.csv, skipped_groups.csv, not_fitted.csv, apertures.csv,
                                 phase_map.csv, phase_cuts.csv, placeholders.csv, logs/
    results/afrho_ztf.csv        ZTF dust context per (target, phase) -- written by
                                 ../scripts/attach_afrho_ztf.py, attached to gas_fit.csv and phase_map.csv
    results/studies/<variant>/   the same tables for each study variant
    results/studies/             flag_policy_*.csv, distcorr_effect_*.csv, apphot_comparison/
    fig/emission_model/, fig/cont_subtract/, fig/phase_group/, fig/summary_*.png
    fig/studies/                 flag_policy_comparison.png, distcorr_effect.png, <variant>/, apphot_comparison/

Environment overrides: ``COMSPEC_ROOT``, ``COMSPEC_APPHOT_DIR``, ``COMSPEC_DATA_DIR``,
``COMSPEC_RESULT_DIR``, ``COMSPEC_FIG_DIR``.
"""

from __future__ import annotations

import os
from pathlib import Path

from .config import MAIN_VARIANT

__all__ = [
    "ROOT", "APPHOT_DIR", "DATA_DIR", "RESULT_DIR", "FIG_DIR", "LOG_DIR", "STUDY_DATA_DIR",
    "STUDY_RESULT_DIR", "STUDY_FIG_DIR", "NOTEBOOK_DIR", "DOC_DIR", "ORBIT_CLASSES_CSV",
    "AFRHO_ZTF_CSV", "variant_dirs", "ensure_dirs", "describe",
]


def _env(key: str, default: Path) -> Path:
    v = os.environ.get(key)
    return Path(v).expanduser() if v else default


ROOT: Path = _env("COMSPEC_ROOT", Path(__file__).resolve().parents[2])

#: Revised aperture photometry -- one CSV per target, written by ``spherex_apphot``.
APPHOT_DIR: Path = _env("COMSPEC_APPHOT_DIR", ROOT / "data" / "apphot")
DATA_DIR: Path = _env("COMSPEC_DATA_DIR", ROOT / "data")
RESULT_DIR: Path = _env("COMSPEC_RESULT_DIR", ROOT / "results")
FIG_DIR: Path = _env("COMSPEC_FIG_DIR", ROOT / "fig")
LOG_DIR: Path = RESULT_DIR / "logs"
STUDY_DATA_DIR: Path = DATA_DIR / "studies"
STUDY_RESULT_DIR: Path = RESULT_DIR / "studies"
STUDY_FIG_DIR: Path = FIG_DIR / "studies"
NOTEBOOK_DIR: Path = ROOT / "notebooks"
DOC_DIR: Path = ROOT / "doc"
ORBIT_CLASSES_CSV: Path = ROOT / "data" / "reference" / "comet_orbit_classes.csv"
#: ZTF dust context -- A(0°)fρ at each phase's mean r_h, from the ``ztf-comet`` project
#: (``../scripts/attach_afrho_ztf.py`` writes it).  Optional: attached to ``gas_fit.csv`` and
#: ``phase_map.csv`` by :func:`dataio.attach_afrho_ztf` whenever the file is present.
AFRHO_ZTF_CSV: Path = RESULT_DIR / "afrho_ztf.csv"


def variant_dirs(name: str) -> dict:
    """
    Output directories for one pipeline variant.

    The main variant writes at the top of ``data/``, ``results/`` and ``fig/``;
    every other variant writes under ``studies/<name>/``.  Keeping each study in
    its own tree is what lets the contamination, badphot and distance-correction
    studies compare complete runs rather than overwritten ones.
    """
    if name == MAIN_VARIANT:
        d, r, f = DATA_DIR, RESULT_DIR, FIG_DIR
    else:
        d, r, f = STUDY_DATA_DIR / name, STUDY_RESULT_DIR / name, STUDY_FIG_DIR / name
    return {
        "emission": d / "emission",
        "results": r,
        "lines": r / "gas_fit_lines",
        "fig": f,
        "fig_cont": f / "cont_subtract",
        "fig_fit": f / "emission_model",
    }


def ensure_dirs(variant: str | None = None) -> None:
    """Create the shared output directories, and a variant's if one is named."""
    for d in (DATA_DIR, RESULT_DIR, FIG_DIR, LOG_DIR, FIG_DIR / "phase_group",
              STUDY_DATA_DIR, STUDY_RESULT_DIR, STUDY_FIG_DIR):
        d.mkdir(parents=True, exist_ok=True)
    if variant:
        for d in variant_dirs(variant).values():
            d.mkdir(parents=True, exist_ok=True)


def describe() -> str:
    rows = [("ROOT", ROOT), ("APPHOT_DIR", APPHOT_DIR), ("DATA_DIR", DATA_DIR),
            ("RESULT_DIR", RESULT_DIR), ("FIG_DIR", FIG_DIR), ("STUDIES", STUDY_RESULT_DIR)]
    return "\n".join(f"  {k:<12s} [{'ok     ' if v.exists() else 'MISSING'}] {v}" for k, v in rows)
