#!/usr/bin/env python
"""
Figure generation for targets already processed by ``main.py``.

Kept separate from ``main.py`` on purpose: the batch photometry path draws
nothing, so a production run is never slowed down by rendering.  This script
reads the CSV that ``main.py`` wrote and produces the full figure set.

Per cutout (``fig/<target>/cutouts/``)
    ``<stem>_cutout.png``  raw science, FLAG plane, and masked science with every
                           Gaia source marked at a size proportional to its flux
    ``<stem>_growth.png``  enclosed flux and annulus surface brightness vs radius

Per target (``fig/<target>/``)
    spectrum and reflectance per epoch, band stacks, radial profiles, the
    contamination analysis, the flag census and the aperture-acceptance plot.

Examples
--------
Every figure for both sample targets::

    python scripts/apphot/make_figures.py --sample --objdesig 2P 24P

Target summaries only, skipping the ~2000 per-cutout files::

    python scripts/apphot/make_figures.py --objdesig 24P --no-cutouts
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd  # noqa: E402

from spherex_apphot import Config, __version__  # noqa: E402
from spherex_apphot import directory as _dir  # noqa: E402
from spherex_apphot.catalog import GaiaCatalog  # noqa: E402
from spherex_apphot.diagnostics import contamination_summary, growth_metrics  # noqa: E402
from spherex_apphot.fitsio import FitsResolver  # noqa: E402
from spherex_apphot.logging_utils import get_logger, setup_logging  # noqa: E402
from spherex_apphot.pipeline import stack_target  # noqa: E402
from spherex_apphot.status import slugify  # noqa: E402
from spherex_apphot.targetlist import read_target_list  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--objdesig", type=str, nargs="+",
                     help="targets to render, e.g. 2P 24P")
    sel.add_argument("--target-list", type=Path,
                     help="spreadsheet or text file of designations")
    p.add_argument("--workers", type=int, default=1,
                   help="render this many targets in parallel (default: %(default)s)")
    p.add_argument("--skip-existing", action="store_true",
                   help="skip figures that are already on disk (per exposure), "
                        "so an interrupted render resumes cheaply")
    p.add_argument("--sample", action="store_true",
                   help="use the in-repo data/spherex_sample/ tree")
    p.add_argument("--apphot-dir", type=Path, default=None,
                   help="where the photometry CSVs live (default: results/apphot)")
    p.add_argument("--fig-dir", type=Path, default=None,
                   help="figure root (default: fig/)")
    p.add_argument("--fits-root", type=Path, nargs="+", default=None)
    p.add_argument("--ap-label", type=str, default=None,
                   help="reference aperture (default: the smallest fixed-pixel one)")
    p.add_argument("--max-cutouts", type=int, default=None,
                   help="cap the per-cutout figures (default: every exposure)")
    p.add_argument("--no-cutouts", action="store_true",
                   help="skip the per-cutout figures entirely")
    p.add_argument("--no-target", action="store_true",
                   help="skip the per-target summary figures")
    p.add_argument("--cutout-dpi", type=int, default=50,
                   help="dpi for the batch per-cutout figures (default: %(default)s)")
    p.add_argument("--target-dpi", type=int, default=200,
                   help="dpi for the per-target figures (default: %(default)s)")
    p.add_argument("--version", action="version", version=f"spherex_apphot {__version__}")
    return p


def _render_one(target, args, cfg, apphot_dir, fig_root, resolver, gaia, log) -> int:
    """Render every figure for one target.  Returns the number written."""
    from spherex_apphot import plotting

    slug = slugify(target)
    csv = apphot_dir / f"{slug}.csv"
    if not csv.is_file():
        log.error("%s: no photometry at %s -- run main.py first", target, csv)
        return 0

    out_dir = fig_root / slug
    phot = pd.read_csv(csv, dtype={"sourceflag": str})
    log.info("%s: %d rows, %d exposures", target, len(phot), phot["filename"].nunique())

    subset = gaia.cone_search(phot["ra"], phot["dec"],
                              radius_deg=cfg.gaia_search_radius_deg,
                              gmag_limit=cfg.gaia_gmag_limit)
    n = 0
    if not args.no_cutouts:
        n += len(plotting.save_cutout_figures(
            phot, resolver, cfg, gaia_subset=subset, objdesig=target,
            out_dir=out_dir / "cutouts", ap_label=args.ap_label,
            max_cutouts=args.max_cutouts, dpi=args.cutout_dpi,
            skip_existing=args.skip_existing))

    if not args.no_target:
        stacks = stack_target(phot, cfg, objdesig=target, ap_label=args.ap_label,
                              fits_roots=args.fits_root, write=False)
        metrics = growth_metrics(phot, ap_kind="km")
        if not metrics.empty:
            log.info("%s contamination analysis:\n%s", target,
                     contamination_summary(metrics, threshold=3.0))
        n += len(plotting.save_target_figures(
            phot, cfg, objdesig=target, stacks=stacks, metrics=metrics,
            ap_label=args.ap_label, out_dir=out_dir, dpi=args.target_dpi))
    return n


def _render_worker(target: str, args) -> int:
    """Pool entry point: rebuild the per-process state, then render."""
    if args.sample:
        _dir.use_sample_data()
    setup_logging(_dir.LOG_DIR, run_name=f"figures_{slugify(target)}", console=False)
    log = get_logger("make_figures")
    import matplotlib
    matplotlib.use("Agg")
    from spherex_apphot import plotting
    plotting.apply_rcparams()

    cfg = Config()
    cfg.validate()
    return _render_one(
        target, args, cfg,
        args.apphot_dir or _dir.APPHOT_DIR,
        args.fig_dir or _dir.FIG_DIR,
        FitsResolver(args.fits_root if args.fits_root is not None else _dir.FITS_ROOTS),
        GaiaCatalog(cache_dir=_dir.GAIA_CACHE_DIR, src_npy=_dir.GAIA_ALL_NPY),
        log)


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.sample:
        _dir.use_sample_data()
    _dir.ensure_dirs()
    setup_logging(_dir.LOG_DIR, run_name="figures")
    log = get_logger("make_figures")

    # matplotlib is imported only here, so the photometry path never pays for it.
    import matplotlib
    matplotlib.use("Agg")
    from spherex_apphot import plotting
    plotting.apply_rcparams()

    cfg = Config()
    cfg.validate()
    apphot_dir = args.apphot_dir or _dir.APPHOT_DIR
    fig_root = args.fig_dir or _dir.FIG_DIR
    resolver = FitsResolver(args.fits_root if args.fits_root is not None else _dir.FITS_ROOTS)
    gaia = GaiaCatalog(cache_dir=_dir.GAIA_CACHE_DIR, src_npy=_dir.GAIA_ALL_NPY)

    targets = (read_target_list(args.target_list) if args.target_list
               else list(args.objdesig))
    log.info("%d target(s) to render", len(targets))

    t0 = time.time()
    n_written = 0
    if args.workers > 1 and len(targets) > 1:
        # Rendering is CPU-bound in matplotlib and completely independent per
        # target, so a process pool scales almost linearly.
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = {pool.submit(_render_worker, t, args): t for t in targets}
            for i, fut in enumerate(as_completed(futs), start=1):
                target = futs[fut]
                try:
                    n = fut.result()
                except Exception as exc:                # noqa: BLE001
                    log.error("%s: rendering failed: %s", target, exc)
                    n = 0
                n_written += n
                log.info("[%d/%d] %s -> %d figures", i, len(targets), target, n)
    else:
        for i, target in enumerate(targets, start=1):
            log.info("[%d/%d] %s", i, len(targets), target)
            n_written += _render_one(target, args, cfg, apphot_dir, fig_root,
                                     resolver, gaia, log)

    log.info("=" * 60)
    log.info("%d figures in %.1f min -> %s", n_written, (time.time() - t0) / 60.0, fig_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
