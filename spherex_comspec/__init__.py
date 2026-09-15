"""
spherex_comspec -- gas production rates from SPHEREx comet spectrophotometry.

End-to-end model from the revised aperture photometry (``results/apphot/photometry/``)
to production rates of H2O, CO2 and CO:

====================  =========================================================
:mod:`config`         every constant, threshold, window and placeholder
:mod:`directory`      filesystem layout
:mod:`dataio`         photometry reading, aperture choice, flag policy, all I/O
:mod:`grouping`       28-day epochs -> single-state ``phase`` groups
:mod:`continuum`      local polynomial continuum subtraction and validation
:mod:`fluorescence`   reconstructed GSFC-style g-factor templates and the CO Swings factor
:mod:`gasmodel`       Haser coma, Yamamoto filling factor, fluorescence (Layers 0-4)
:mod:`instrument`     SPHEREx channel bandpass (Layer 5)
:mod:`fitting`        design matrix, weighted linear solve, limits (Layer 6)
:mod:`pipeline`       orchestration: :func:`pipeline.run_variant`
:mod:`analysis`       flag-contamination and distance-correction studies
:mod:`plotting`       every figure; never imported by the batch path
====================  =========================================================

Quick start::

    from spherex_comspec import run_variant, DEFAULT_VARIANTS
    res = run_variant(DEFAULT_VARIANTS[0])       # dc_main
"""

from __future__ import annotations

__version__ = "1.3.0"

from .config import (BANDS, BAND_WINDOWS, DEFAULT_VARIANTS, EMISSION_WINDOWS, KEY_RANGES,
                     MAIN_VARIANT, PLACEHOLDERS, SPECIES, VARIANTS, ApertureConfig,
                     ContinuumConfig, FitConfig, FlagPolicy, GroupingConfig, ModelParams, Variant)
from .fitting import FitResult, fit_production_rates, model_curves
from .pipeline import run_grouping, run_variant

__all__ = [
    "__version__", "BANDS", "BAND_WINDOWS", "DEFAULT_VARIANTS", "EMISSION_WINDOWS", "KEY_RANGES",
    "MAIN_VARIANT", "PLACEHOLDERS", "SPECIES", "VARIANTS", "ApertureConfig", "ContinuumConfig",
    "FitConfig",
    "FlagPolicy", "GroupingConfig", "ModelParams", "Variant", "FitResult",
    "fit_production_rates", "model_curves", "run_grouping", "run_variant",
]
