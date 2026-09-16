#!/usr/bin/env python
"""
Copy the cutout FITS for a set of targets into a per-target directory tree.

``comets_v5_fits`` is a single flat directory of ~139 000 files on an external
drive, which is slow to enumerate and slow to stat.  This script builds
``comets_v5_fits_filtered/<target>/`` holding only the exposures of the targets
you actually work on -- for the 68-comet working list that is ~28 000 files and
~3.6 GB, about a quarter of the archive.

The file list comes from the Parquet index, **not** from listing the source
directory, so the slow directory is never enumerated: each source path is
constructed directly from a known filename.

The copy is non-destructive (the archive is left untouched) and resumable
(existing files of the right size are skipped), so it is safe to interrupt and
re-run.

Examples
--------
Extract everything in the working list::

    python scripts/apphot/extract_targets.py --target-list data/reference/sx_comet_list_ver2607.xlsx

Just two comets, into a different tree::

    python scripts/apphot/extract_targets.py --objdesig 2P 24P --dest /tmp/subset
"""

from __future__ import annotations

import argparse
import os

# macOS writes an AppleDouble "._name" sidecar beside every file it copies to a
# filesystem without native extended-attribute support -- which exFAT, the T7
# format, is.  At 27 797 files and a 128 kB exFAT cluster each, those 4 kB
# sidecars cost ~3.5 GB of slack and double the entry count of every directory.
# Setting this before shutil is imported suppresses them.
os.environ.setdefault("COPYFILE_DISABLE", "1")

import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd  # noqa: E402

from spherex_apphot import __version__  # noqa: E402
from spherex_apphot import directory as _dir  # noqa: E402
from spherex_apphot.logging_utils import get_logger, setup_logging  # noqa: E402
from spherex_apphot.status import slugify  # noqa: E402
from spherex_apphot.targetlist import read_target_list  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--objdesig", type=str, nargs="+", help="target designations")
    sel.add_argument("--target-list", type=Path,
                     help="spreadsheet or text file of designations")
    sel.add_argument("--all", action="store_true", help="every target in the index")

    p.add_argument("--db", type=Path, default=None, help="cutout index (default: %(default)s)")
    p.add_argument("--src", type=Path, default=None,
                   help="flat source directory (default: comets_v5_fits)")
    p.add_argument("--dest", type=Path, default=None,
                   help="destination root (default: comets_v5_fits_filtered)")
    p.add_argument("--workers", type=int, default=8,
                   help="parallel copy threads (default: %(default)s)")
    p.add_argument("--link", action="store_true",
                   help="hard-link instead of copying (same filesystem only; "
                        "uses no extra space)")
    p.add_argument("--overwrite", action="store_true",
                   help="re-copy files that already exist at the destination")
    p.add_argument("--dry-run", action="store_true", help="report, copy nothing")
    p.add_argument("--clean-appledouble", action="store_true",
                   help="delete ._* sidecars under the destination and exit")
    p.add_argument("--version", action="version", version=f"spherex_apphot {__version__}")
    return p


def _copy_one(src: Path, dst: Path, link: bool, overwrite: bool) -> Tuple[str, int]:
    """
    Place one file at ``dst``.

    Returns
    -------
    (status, nbytes) : ('copied' | 'skipped' | 'missing' | 'error:...', int)

    Notes
    -----
    Copies to a ``.part`` sibling and renames, so an interrupted run never
    leaves a truncated FITS that a later resume would mistake for complete.
    """
    try:
        st = src.stat()
    except OSError:
        return "missing", 0
    if dst.exists() and not overwrite:
        try:
            if dst.stat().st_size == st.st_size:
                return "skipped", 0
        except OSError:
            pass
    try:
        tmp = dst.with_suffix(dst.suffix + ".part")
        if link:
            if tmp.exists():
                tmp.unlink()
            os.link(src, tmp)
        else:
            shutil.copyfile(src, tmp)
        os.replace(tmp, dst)
        return "copied", st.st_size
    except OSError as exc:
        return f"error:{exc}", 0


def _clean_appledouble(root: Path, log) -> int:
    """
    Delete AppleDouble ``._*`` sidecars under ``root``.

    Only files whose name starts with ``._`` **and** which sit beside a real
    file *or directory* of the same name are removed, so nothing that could be
    data is touched.  macOS writes a sidecar for directories too, which is easy
    to miss when the check only looks for a matching file.
    """
    if not root.is_dir():
        log.error("no such directory: %s", root)
        return 1
    n = freed = 0
    for side in root.rglob("._*"):
        if not side.is_file():
            continue
        real = side.with_name(side.name[2:])
        if not real.exists():
            log.warning("keeping %s: no matching file or directory", side.name)
            continue
        try:
            freed += side.stat().st_size
            side.unlink()
            n += 1
        except OSError as exc:
            log.warning("could not remove %s: %s", side.name, exc)
    log.info("removed %d AppleDouble sidecar(s), %.1f MB of payload", n, freed / 1e6)
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    _dir.ensure_dirs()
    setup_logging(_dir.LOG_DIR, run_name="extract_targets")
    log = get_logger("extract_targets")

    db_path = args.db or _dir.DB_PATH
    src_root = args.src or _dir.FITS_DIR
    dest_root = args.dest or _dir.FITS_FILTERED_DIR

    if args.target_list:
        targets = read_target_list(args.target_list)
    elif args.objdesig:
        targets = list(args.objdesig)
    else:
        targets = None

    log.info("index : %s", db_path)
    log.info("source: %s", src_root)
    log.info("dest  : %s", dest_root)

    if args.clean_appledouble:
        return _clean_appledouble(dest_root, log)

    db = pd.read_parquet(db_path, columns=["objdesig", "filename"])
    if targets is not None:
        known = set(db["objdesig"].unique())
        unknown = [t for t in targets if t not in known]
        if unknown:
            log.warning("%d designation(s) absent from the index and skipped: %s",
                        len(unknown), ", ".join(map(str, unknown)))
        db = db[db["objdesig"].isin(targets)]
    if db.empty:
        log.error("no rows matched; nothing to do")
        return 1

    per_target = db.groupby("objdesig").size()
    log.info("%d target(s), %d file(s) to place", len(per_target), len(db))

    # Build the work list from the index: the flat source directory, which holds
    # ~139 000 entries, is never enumerated.
    jobs: List[Tuple[Path, Path]] = []
    for objdesig, grp in db.groupby("objdesig"):
        out_dir = dest_root / slugify(str(objdesig))
        if not args.dry_run:
            out_dir.mkdir(parents=True, exist_ok=True)
        for fn in grp["filename"]:
            # AppleDouble sidecars are metadata, not data.
            if str(fn).startswith("._"):
                continue
            jobs.append((src_root / str(fn), out_dir / str(fn)))

    if args.dry_run:
        # Stat in parallel: a serial pass over ~28 000 files on the external
        # drive takes minutes, which defeats the purpose of a dry run.
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            present = sum(pool.map(lambda sd: sd[0].is_file(), jobs))
        log.info("dry run: %d of %d source file(s) present", present, len(jobs))
        for t, n in per_target.items():
            log.info("  %-14s %5d", t, n)
        return 0

    t0 = time.time()
    counts = {"copied": 0, "skipped": 0, "missing": 0, "error": 0}
    nbytes = 0
    missing_examples: List[str] = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_copy_one, s, d, args.link, args.overwrite): s
                   for s, d in jobs}
        for i, fut in enumerate(as_completed(futures), start=1):
            status, size = fut.result()
            nbytes += size
            if status.startswith("error"):
                counts["error"] += 1
                if counts["error"] <= 5:
                    log.warning("%s: %s", futures[fut].name, status)
            else:
                counts[status] += 1
                if status == "missing" and len(missing_examples) < 5:
                    missing_examples.append(futures[fut].name)
            if i % 2000 == 0:
                log.info("  %d/%d  (%.1f GB, %.0f files/s)", i, len(jobs),
                         nbytes / 1e9, i / max(time.time() - t0, 1e-9))

    dt = time.time() - t0
    log.info("=" * 60)
    log.info("copied %d, skipped %d, missing %d, errors %d  (%.2f GB in %.1f min)",
             counts["copied"], counts["skipped"], counts["missing"], counts["error"],
             nbytes / 1e9, dt / 60.0)
    if missing_examples:
        log.warning("missing source files, e.g.: %s", ", ".join(missing_examples))
    log.info("tree: %s", dest_root)
    return 0 if counts["error"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
