"""
SPHEREx comet aperture photometry.

A rewrite of the primitive ``apphot.py`` prototype following the review in
``doc/apphot/code_review_primitive.md``.  What changed, and why, is summarised in
``doc/apphot/pipeline_upgrade_notes.md``.

Layout
------
================== =========================================================
:mod:`config`      every tunable number, frozen and serialisable
:mod:`directory`   filesystem layout, with environment overrides
:mod:`logging_utils` run logging
:mod:`status`      the ``results/apphot/status.csv`` run table and target slugs
:mod:`fitsio`      cutout FITS and Parquet index reading
:mod:`wcsutil`     WCS from a header, or rebuilt from an index row
:mod:`masking`     bad-pixel identification (flags, NaN, bad variance)
:mod:`catalog`     Gaia DR3 access with a declination-sorted cache
:mod:`apertures`   per-exposure aperture construction and validity
:mod:`sourceflag`  Gaia contamination flags (sources are recorded, not masked)
:mod:`phot`        the aperture photometry itself
:mod:`epochs`      grouping exposures into observing epochs
:mod:`stacking`    robust band stacking (sigma-clipped median)
:mod:`reflectance` reflectance spectra
:mod:`diagnostics` growth-curve metrics and source-flag effectiveness
:mod:`pipeline`    orchestration: :func:`pipeline.run_target`
:mod:`plotting`    figures -- imported by notebooks only, never by ``main.py``
================== =========================================================

Quick start
-----------
>>> from spherex_apphot import Config, run_target          # doctest: +SKIP
>>> res = run_target("24P", Config())                      # doctest: +SKIP
>>> res.phot.head()                                        # doctest: +SKIP
"""

from __future__ import annotations

__version__ = "2.2.0"

from .config import Config, FLAG_BIT_MEANINGS, DEFAULT_BAD_FLAG_BITS
from .pipeline import TargetResult, run_target, stack_target
from .status import StatusLog, slugify

__all__ = [
    "__version__",
    "Config", "FLAG_BIT_MEANINGS", "DEFAULT_BAD_FLAG_BITS",
    "TargetResult", "run_target", "stack_target",
    "StatusLog", "slugify",
]
