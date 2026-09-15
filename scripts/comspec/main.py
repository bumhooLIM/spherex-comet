#!/usr/bin/env python
"""
End-to-end driver for ``spherex_comspec``.

Stages
------
``group``     regroup the 28-day epochs into single-state phases (once, shared by every variant)
``run``       continuum subtraction + production-rate fitting for one or more variants
``analyze``   flag-contamination and distance-correction studies across the variants
``figures``   every figure (grouping, continuum validation, emission-model fits, summaries)
``all``       the four in order

Examples
--------
Everything, all default variants::

    python scripts/comspec/main.py all

Only the main variant, then its figures::

    python scripts/comspec/main.py run --variants dc_main
    python scripts/comspec/main.py figures --variants dc_main

Two targets, for a quick look::

    python scripts/comspec/main.py run --variants dc_main --targets 24P 2P
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from spherex_comspec import __version__  # noqa: E402
from spherex_comspec import directory as _dir  # noqa: E402
from spherex_comspec.config import MAIN_VARIANT, VARIANTS  # noqa: E402
from spherex_comspec.logging_utils import get_logger, setup_logging  # noqa: E402


def _studies(names):
    """
    Which selected variants play which part in the cross-variant studies, by role.

    Returns ``(main, flag_variants, raw, extra)``: the flag census compares ``main`` with
    the ``flags`` runs, the distance-correction study pairs it with the ``distcorr`` run,
    and the ``badphot`` / ``previous`` / ``fluorescence`` / ``errors`` / ``rules`` runs join the
    census as extra comparisons.
    """
    main_name = MAIN_VARIANT if MAIN_VARIANT in names else names[0]

    def role(r):
        return [n for n in names if VARIANTS[n].role == r and n != main_name]

    raw = role("distcorr")
    return main_name, [main_name] + role("flags"), (raw[0] if raw else None), \
        role("badphot") + role("previous") + role("fluorescence") + role("errors") + role("rules")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stage", choices=["group", "run", "analyze", "figures", "all"])
    p.add_argument("--variants", nargs="+", default=None,
                   help=f"variant names (default: all of {list(VARIANTS)})")
    p.add_argument("--targets", nargs="+", default=None, help="restrict to these targets")
    p.add_argument("--max-groups", type=int, default=None,
                   help="figures: cap the per-group figures per variant")
    p.add_argument("--no-validation-figs", action="store_true")
    p.add_argument("--no-fit-figs", action="store_true")
    p.add_argument("--study-figs", action="store_true",
                   help="figures: also draw the per-group figures of the study variants "
                        "(by default only the main variant gets them; studies get their summaries)")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING"])
    p.add_argument("--version", action="version", version=f"spherex_comspec {__version__}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    _dir.ensure_dirs()
    setup_logging(_dir.LOG_DIR, run_name=f"comspec_{args.stage}",
                  level=getattr(logging, args.log_level))
    log = get_logger("main")
    log.info("spherex_comspec %s\n%s", __version__, _dir.describe())
    names = args.variants or list(VARIANTS)
    unknown = [n for n in names if n not in VARIANTS]
    if unknown:
        log.error("unknown variant(s) %s; choose from %s", unknown, list(VARIANTS))
        return 2
    t0 = time.time()

    if args.stage in ("group", "all"):
        from spherex_comspec.pipeline import run_grouping
        run_grouping(args.targets)

    if args.stage in ("run", "all"):
        from spherex_comspec.dataio import PhaseAssignment
        from spherex_comspec.pipeline import run_variant
        assignment = PhaseAssignment.load()
        for n in names:
            run_variant(VARIANTS[n], args.targets, assignment)

    if args.stage in ("analyze", "all"):
        from spherex_comspec.analysis import write_analysis
        main_name, flag_variants, raw, extra = _studies(names)
        write_analysis(main_name, flag_variants, raw, extra)

    if args.stage in ("figures", "all"):
        from spherex_comspec.dataio import PhaseAssignment
        from spherex_comspec.grouping import regroup_all
        from spherex_comspec.plotting import (apply_rcparams, plot_distcorr_effect,
                                              plot_flag_comparison, save_grouping_figures,
                                              save_variant_figures, savefig)
        from spherex_comspec.analysis import compare_variants, distcorr_effect
        apply_rcparams()
        assignment = PhaseAssignment.load()
        if args.targets:
            assignment = assignment[assignment.target.isin(args.targets)]
        main_name, flag_variants, raw, extra = _studies(names)
        main_v = VARIANTS[main_name]
        _, gmap, _, epochs = regroup_all(sorted(assignment.target.unique()), main_v.grouping)
        save_grouping_figures(epochs, gmap, assignment, main_v.aperture, main_v,
                              max_targets=None)
        for n in names:
            per_group = (n == main_name) or args.study_figs
            save_variant_figures(VARIANTS[n], assignment, max_groups=args.max_groups,
                                 validation=per_group and not args.no_validation_figs,
                                 fits=per_group and not args.no_fit_figs)
        partners = [n for n in flag_variants if n != main_name] + extra
        if partners:
            cmp = compare_variants(main_name, partners)
            fig = plot_flag_comparison(cmp["census"], cmp["pairs"], cmp["paired"], main_name)
            savefig(fig, _dir.STUDY_FIG_DIR / "flag_policy_comparison.png", dpi=200)
        if raw:
            dce = distcorr_effect(main_name, raw)
            fig = plot_distcorr_effect(dce["paired"], dce["band_rows"], main_name, raw)
            savefig(fig, _dir.STUDY_FIG_DIR / "distcorr_effect.png", dpi=200)

    log.info("%s finished in %.1f min", args.stage, (time.time() - t0) / 60)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
