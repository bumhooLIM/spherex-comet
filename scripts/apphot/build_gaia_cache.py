#!/usr/bin/env python
"""
One-off builder for the declination-sorted Gaia cache.

The pipeline runs without it -- ``GaiaCatalog`` falls back to a chunked scan --
but each fallback query costs a full pass over the 11.7 GB catalogue.  Building
the cache once turns every later query into a few-MB slice.

Examples
--------
Build the full cache (~8 GB on disk, memory-bounded)::

    python scripts/apphot/build_gaia_cache.py

Keep only the sources the pipeline can use, for a smaller cache::

    python scripts/apphot/build_gaia_cache.py --gmag-limit 18.0
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from spherex_apphot import directory  # noqa: E402
from spherex_apphot.catalog import build_dec_cache  # noqa: E402
from spherex_apphot.logging_utils import setup_logging, get_logger  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", type=Path, default=directory.GAIA_ALL_NPY,
                    help="source Gaia .npy (default: %(default)s)")
    ap.add_argument("--out", type=Path, default=directory.GAIA_CACHE_DIR,
                    help="cache directory to write (default: %(default)s)")
    ap.add_argument("--gmag-limit", type=float, default=None,
                    help="drop sources fainter than this G magnitude")
    ap.add_argument("--chunk", type=int, default=20_000_000,
                    help="rows per streaming pass (default: %(default)s)")
    ap.add_argument("--bucket-deg", type=float, default=0.5,
                    help="declination bucket width for the external sort")
    args = ap.parse_args(argv)

    directory.ensure_dirs()
    setup_logging(directory.LOG_DIR, run_name="gaia_cache")
    log = get_logger("build_gaia_cache")

    t0 = time.time()
    try:
        out = build_dec_cache(args.src, args.out, gmag_limit=args.gmag_limit,
                              chunk=args.chunk, bucket_deg=args.bucket_deg)
    except Exception:
        log.exception("cache build failed")
        return 1
    log.info("done in %.1f min -> %s", (time.time() - t0) / 60.0, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
