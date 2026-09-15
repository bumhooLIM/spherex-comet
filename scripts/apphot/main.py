#!/usr/bin/env python
"""
Batch driver for the SPHEREx comet aperture photometry pipeline.

Produces photometry only -- one CSV plus a JSON sidecar per target, and an
updated ``results/apphot/status.csv``.  It deliberately draws no figures: plotting
lives in ``notebooks/apphot/main.ipynb``, which calls the very same
``pipeline.run_target``.

Examples
--------
One target against the production dataset::

    python scripts/apphot/main.py --objdesig 24P

Every comet in the working list, four at a time, resuming a previous run::

    python scripts/apphot/main.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 4

A quick end-to-end check against the sample data in the repository::

    python scripts/apphot/main.py --sample --objdesig 2P --limit 40

Notes
-----
``--all`` loads the Gaia catalogue once for the whole run.  The primitive
orchestrator spawned a subprocess per target and so paid an 11.7 GB catalogue
read 445 times (review items P1 and P5).  With ``--workers`` each process opens
the declination cache as a memory map, which the operating system shares.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from spherex_apphot import Config, __version__  # noqa: E402
from spherex_apphot import directory as _dir  # noqa: E402
from spherex_apphot.catalog import GaiaCatalog  # noqa: E402
from spherex_apphot.fitsio import list_targets  # noqa: E402
from spherex_apphot.logging_utils import get_logger, log_config, setup_logging  # noqa: E402
from spherex_apphot.pipeline import run_target, stack_target  # noqa: E402
from spherex_apphot.status import StatusLog  # noqa: E402
from spherex_apphot.targetlist import read_target_list  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--objdesig", type=str, nargs="+",
                     help="one or more target designations, e.g. 24P '2021 G2'")
    sel.add_argument("--all", action="store_true", help="process every target in the index")
    sel.add_argument("--target-list", type=Path,
                     help="spreadsheet or text file of designations, e.g. "
                          "data/reference/sx_comet_list_ver2607.xlsx")

    p.add_argument("--sample", action="store_true",
                   help="use the in-repo data/spherex_sample/ tree instead of the external SSD")
    p.add_argument("--db", type=Path, default=None, help="override the cutout index path")
    p.add_argument("--fits-root", type=Path, nargs="+", default=None,
                   help="override the FITS search root(s)")
    p.add_argument("--out-dir", type=Path, default=None, help="override the output directory")

    p.add_argument("--limit", type=int, default=None,
                   help="process at most this many exposures per target")
    p.add_argument("--workers", type=int, default=1,
                   help="parallel target workers (default: %(default)s)")
    p.add_argument("--force", action="store_true",
                   help="reprocess targets already marked ok in the status table")
    p.add_argument("--stack", action="store_true",
                   help="also write the per-epoch band stacks (FITS)")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--version", action="version", version=f"spherex_apphot {__version__}")
    return p


def _resolve_targets(args, status: StatusLog, log) -> List[str]:
    db = args.db or _dir.DB_PATH
    if args.all:
        targets = list_targets(db)
    elif args.target_list:
        targets = read_target_list(args.target_list)
        known = set(list_targets(db))
        unknown = [t for t in targets if t not in known]
        if unknown:
            log.warning("%d designation(s) absent from the index and skipped: %s",
                        len(unknown), ", ".join(map(str, unknown)))
            targets = [t for t in targets if t in known]
    else:
        targets = list(args.objdesig)
    # Resuming applies to any multi-target run, not just --all.
    if (args.all or args.target_list) and not args.force:
        keep = [t for t in targets if not status.is_done(t)]
        skipped = len(targets) - len(keep)
        if skipped:
            log.info("skipping %d target(s) already marked ok (use --force to redo)", skipped)
        targets = keep
    return targets


def _run_one(objdesig: str, cfg: Config, args, gaia, status) -> str:
    log = get_logger("main")
    res = run_target(
        objdesig, cfg,
        db_path=args.db, fits_roots=args.fits_root,
        gaia=gaia, out_dir=args.out_dir, status=status,
        limit=args.limit,
    )
    if res.status == "ok" and args.stack and res.phot is not None:
        try:
            stacks = stack_target(res.phot, cfg, objdesig=objdesig,
                                  fits_roots=args.fits_root)
            log.info("%s: wrote %d stack(s)", objdesig, len(stacks))
        except Exception:                        # noqa: BLE001
            log.exception("%s: stacking failed (photometry was written)", objdesig)
    return res.status


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.sample:
        _dir.use_sample_data()
    if args.db is None:
        args.db = _dir.DB_PATH
    if args.fits_root is None:
        args.fits_root = _dir.FITS_ROOTS

    _dir.ensure_dirs()
    import logging
    level = getattr(logging, args.log_level)
    setup_logging(_dir.LOG_DIR, run_name="apphot", level=level, console_level=level)
    log = get_logger("main")
    log.info("spherex_apphot %s | pid %d", __version__, os.getpid())
    log.info("paths:\n%s", _dir.describe())

    cfg = Config()
    cfg.validate()
    log_config(cfg, log)

    status = StatusLog(_dir.STATUS_CSV)
    targets = _resolve_targets(args, status, log)
    if not targets:
        log.info("nothing to do")
        return 0
    log.info("%d target(s) to process", len(targets))

    t0 = time.time()
    outcomes: List[str] = []

    if args.workers > 1 and len(targets) > 1:
        # Each worker opens the declination cache as its own memory map; the
        # page cache means the bytes are shared rather than duplicated.
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_worker, t, args, cfg): t for t in targets}
            for i, fut in enumerate(as_completed(futures), start=1):
                target = futures[fut]
                try:
                    outcomes.append(fut.result())
                except Exception as exc:         # noqa: BLE001
                    log.error("%s: worker crashed: %s", target, exc)
                    outcomes.append("failed")
                log.info("[%d/%d] %s -> %s", i, len(targets), target, outcomes[-1])
    else:
        gaia = GaiaCatalog(cache_dir=_dir.GAIA_CACHE_DIR, src_npy=_dir.GAIA_ALL_NPY)
        for i, target in enumerate(targets, start=1):
            log.info("[%d/%d] %s", i, len(targets), target)
            outcomes.append(_run_one(target, cfg, args, gaia, status))

    dt = time.time() - t0
    ok = outcomes.count("ok")
    log.info("=" * 60)
    log.info("done: %d ok, %d not ok, in %.1f min", ok, len(outcomes) - ok, dt / 60.0)
    log.info("status table: %s -- %s", _dir.STATUS_CSV, status.summary())
    return 0 if ok == len(outcomes) else 1


def _worker(objdesig: str, args, cfg: Config) -> str:
    """Entry point for a pool worker: fresh logging, fresh catalogue handle."""
    if args.sample:
        _dir.use_sample_data()
    setup_logging(_dir.LOG_DIR, run_name=f"apphot_{objdesig}", console=False)
    gaia = GaiaCatalog(cache_dir=_dir.GAIA_CACHE_DIR, src_npy=_dir.GAIA_ALL_NPY)
    status = StatusLog(_dir.STATUS_CSV)
    return _run_one(objdesig, cfg, args, gaia, status)


if __name__ == "__main__":
    raise SystemExit(main())
