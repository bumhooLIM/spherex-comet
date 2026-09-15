"""
Filesystem layout for ``spherex_comspec``.

Every path the package reads or writes is resolved here.  ``spherex_comspec`` is the
third stage of the ``spherex-comet`` project (``ztfcomet`` -> ``spherex_apphot`` ->
``spherex_comspec``), so its tree is namespaced under ``comspec/``.  The layout is flat
on purpose -- the main result must be readable in one glance -- and every product of a
*study* variant (a flag-policy, badphot, distance-correction or previous baseline run)
is kept one level down under ``studies/`` so it can never be mistaken for the catalog
result::

    results/apphot/photometry/          the SPHEREx aperture photometry of spherex_apphot (input, never written to)
    data/comspec/phase_assignment.csv   per-exposure phase labels, shared by every variant
    data/comspec/emission/              continuum summaries + point spectra of the main variant
    data/comspec/studies/<variant>/emission/
    data/reference/                     orbit classes, literature Q tables, the previous phase map, the working list
    data/fluorescence/                  the reconstructed g-factor database (species templates, CO Swings table)
    results/comspec/gas_fit.csv         THE production rates (main variant), gas_fit_lines/, run.meta.json,
                                        continuum_summary.csv, skipped_groups.csv, not_fitted.csv, apertures.csv,
                                        phase_map.csv, phase_cuts.csv, placeholders.csv, logs/
    results/comspec/afrho_ztf.csv       ZTF dust context per (target, phase) -- written by
                                        scripts/comspec/attach_afrho_ztf.py from results/ztf/afrho/spherex_afrho.csv,
                                        attached to gas_fit.csv and phase_map.csv
    results/comspec/studies/<variant>/  the same tables for each study variant
    results/comspec/studies/            flag_policy_*.csv, distcorr_effect_*.csv, apphot_comparison/, method_matrix/
    fig/comspec/emission_model/, cont_subtract/, phase_group/, phase_images/, summary_*.png
    fig/comspec/studies/                flag_policy_comparison.png, distcorr_effect.png, <variant>/, apphot_comparison/

Environment overrides: ``COMSPEC_ROOT``, ``COMSPEC_APPHOT_DIR``, ``COMSPEC_DATA_DIR``,
``COMSPEC_RESULT_DIR``, ``COMSPEC_FIG_DIR`` (an isolated run sets them before import).
"""

from __future__ import annotations

import os
from pathlib import Path

from .config import MAIN_VARIANT

__all__ = [
    "ROOT", "APPHOT_DIR", "DATA_DIR", "RESULT_DIR", "FIG_DIR", "LOG_DIR", "STUDY_DATA_DIR",
    "STUDY_RESULT_DIR", "STUDY_FIG_DIR", "NOTEBOOK_DIR", "DOC_DIR", "REFERENCE_DIR",
    "FLUORESCENCE_DIR", "ORBIT_CLASSES_CSV", "AFRHO_ZTF_CSV", "ZTF_AFRHO_SRC",
    "variant_dirs", "ensure_dirs", "describe",
]


def _env(key: str, default: Path) -> Path:
    v = os.environ.get(key)
    return Path(v).expanduser() if v else default


# directory.py lives at <project>/spherex_comspec/directory.py
ROOT: Path = _env("COMSPEC_ROOT", Path(__file__).resolve().parents[1])

#: The SPHEREx aperture photometry -- one CSV (+ .meta.json) per target, written by
#: ``spherex_apphot`` into ``results/apphot/photometry/``.  Input only: never written to here.
APPHOT_DIR: Path = _env("COMSPEC_APPHOT_DIR", ROOT / "results" / "apphot" / "photometry")
DATA_DIR: Path = _env("COMSPEC_DATA_DIR", ROOT / "data" / "comspec")
RESULT_DIR: Path = _env("COMSPEC_RESULT_DIR", ROOT / "results" / "comspec")
FIG_DIR: Path = _env("COMSPEC_FIG_DIR", ROOT / "fig" / "comspec")
LOG_DIR: Path = RESULT_DIR / "logs"
STUDY_DATA_DIR: Path = DATA_DIR / "studies"
STUDY_RESULT_DIR: Path = RESULT_DIR / "studies"
STUDY_FIG_DIR: Path = FIG_DIR / "studies"
#: Holds the shared ``rcparams.py`` (figure style); the study scripts are in ``notebooks/comspec/``.
NOTEBOOK_DIR: Path = ROOT / "notebooks"
DOC_DIR: Path = ROOT / "doc" / "comspec"
#: Shared reference tables (orbit classes, literature production rates, the previous
#: phase map, the 68-comet working list) and the reconstructed fluorescence database.
REFERENCE_DIR: Path = ROOT / "data" / "reference"
FLUORESCENCE_DIR: Path = _env("COMSPEC_FLUORESCENCE_DIR", ROOT / "data" / "fluorescence")
ORBIT_CLASSES_CSV: Path = REFERENCE_DIR / "comet_orbit_classes.csv"
#: ZTF dust context -- A(0°)fρ at each phase's mean r_h, from the ``ztfcomet`` stage
#: (``scripts/comspec/attach_afrho_ztf.py`` pivots ``ZTF_AFRHO_SRC`` into it).  Optional:
#: attached to ``gas_fit.csv`` and ``phase_map.csv`` by :func:`dataio.attach_afrho_ztf`
#: whenever the file is present.
AFRHO_ZTF_CSV: Path = RESULT_DIR / "afrho_ztf.csv"
ZTF_AFRHO_SRC: Path = _env("COMSPEC_ZTF_AFRHO", ROOT / "results" / "ztf" / "afrho" / "spherex_afrho.csv")


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
            ("RESULT_DIR", RESULT_DIR), ("FIG_DIR", FIG_DIR), ("STUDIES", STUDY_RESULT_DIR),
            ("REFERENCE", REFERENCE_DIR), ("FLUORESCENCE", FLUORESCENCE_DIR),
            ("ZTF_AFRHO", ZTF_AFRHO_SRC)]
    return "\n".join(f"  {k:<12s} [{'ok     ' if v.exists() else 'MISSING'}] {v}" for k, v in rows)
