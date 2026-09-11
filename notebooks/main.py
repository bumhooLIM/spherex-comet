#!/usr/bin/env python
"""End-to-end ZTF comet pipeline: query, download, photometry, figures.

Runs the whole chain for one or several targets, using :mod:`ztfcomet` for
everything.  This script holds no science logic of its own — it is a driver, so
that what runs in batch is the same code the notebooks validate interactively.

Examples
--------
Full run for the default test target::

    python notebooks/main.py 24P

Several targets, query and photometry only::

    python notebooks/main.py 24P 240P 2P --no-figures

Reduce data already on disk, skipping the network entirely::

    python notebooks/main.py 240P --steps phot figures

Re-fetch over a different window with a bigger aperture::

    python notebooks/main.py 2P --start 2025-07-01 --end 2026-05-31 --rho-km 20000

Stage names for ``--steps`` are ``query``, ``download``, ``phot``, ``profile``
and ``figures``; all five run by default. ``profile`` compares the comet's
radial surface-brightness profile with field stars on the same frame.

Quality checks include a Gaia DR3 background-source test: if catalogued stars
inside the aperture carry >=30% of the comet's predicted flux, the frame is
flagged as contaminated.

Comet designations are resolved against Horizons per epoch and **fragments are
excluded** — asking for 240P gives the parent body, not 240P-B. Use
``--allow-fragment`` (or the ``240P-B`` target) when the fragment is the object
of interest.

Notes
-----
Frames failing a quality check are **flagged, not dropped** (see
``doc/primitive_code_analysis.md``).  The summary prints the flag counts, and
``--only-good`` restricts the plots — never the saved table.
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
import sys
from pathlib import Path

# Allow running straight from a checkout without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
import pandas as pd

import ztfcomet as zc

log = logging.getLogger("ztfcomet.main")

STEPS = ("query", "download", "phot", "profile", "figures")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("targets", nargs="+",
                        help="target names, e.g. 24P 240P '2019 Y3'")
    parser.add_argument("--steps", nargs="+", choices=STEPS, default=list(STEPS),
                        help="pipeline stages to run (default: all)")
    parser.add_argument("--start", help="query start date, YYYY-MM-DD (overrides config)")
    parser.add_argument("--end", help="query end date, YYYY-MM-DD (overrides config)")
    parser.add_argument("--interval-days", type=float, help="coarse ephemeris step")
    parser.add_argument("--rh-max", type=float, help="max heliocentric distance, AU")
    parser.add_argument("--vmag-max", type=float, help="max total magnitude")
    parser.add_argument("--cutout-size", help="IRSA cutout size, e.g. 10arcmin")
    parser.add_argument("--full-frame", action="store_true",
                        help="download whole quadrant images instead of cutouts (~37 MB each)")
    parser.add_argument("--rho-km", type=float, help="Af-rho aperture radius at the comet, km")
    parser.add_argument("--no-color-term", action="store_true",
                        help="skip the CLRCOEFF colour correction")
    parser.add_argument("--no-aperture-correction", action="store_true",
                        help="skip the APCOR aperture correction")
    parser.add_argument("--no-contamination", action="store_true",
                        help="skip the Gaia background-source check")
    parser.add_argument("--contam-ratio", type=float,
                        help="flag when background flux reaches this fraction of "
                             "the comet's (default 0.30)")
    parser.add_argument("--allow-fragment", action="store_true",
                        help="permit Horizons to return a fragment (e.g. 240P-B) "
                             "instead of the parent comet")
    parser.add_argument("--only-good", action="store_true",
                        help="plot only unflagged frames (the saved table always keeps everything)")
    parser.add_argument("--no-figures", action="store_true", help="shorthand for dropping the figures step")
    parser.add_argument("--no-cutout-figures", action="store_true",
                        help="skip the per-frame PNGs; still make the lightcurve")
    parser.add_argument("--no-profile-plots", action="store_true",
                        help="skip per-frame radial-profile PNGs; the summary is still made")
    parser.add_argument("--overwrite", action="store_true", help="re-download existing files")
    parser.add_argument("--dpi", type=int, default=50, help="dpi for batch cutout PNGs (default 50)")
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args(argv)


def build_target(name, args):
    """Resolve a target from config and apply any command-line overrides."""
    target = zc.get_target(name)

    query_overrides = {
        k: v for k, v in {
            "interval_days": args.interval_days,
            "rh_max": args.rh_max,
            "vmag_max": args.vmag_max,
            "cutout_size": args.cutout_size,
            "is_cutout": False if args.full_frame else None,
        }.items() if v is not None
    }
    phot_overrides = {
        k: v for k, v in {
            "rho_km": args.rho_km,
            "apply_color_term": False if args.no_color_term else None,
            "apply_aperture_correction": False if args.no_aperture_correction else None,
            "check_contamination": False if args.no_contamination else None,
            "contam_flux_ratio": args.contam_ratio,
        }.items() if v is not None
    }

    changes = {}
    if query_overrides:
        changes["query"] = dataclasses.replace(target.query, **query_overrides)
    if phot_overrides:
        changes["phot"] = dataclasses.replace(target.phot, **phot_overrides)
    if args.start:
        changes["start_date"] = args.start
    if args.end:
        changes["end_date"] = args.end
    if args.allow_fragment:
        changes["allow_fragment"] = True
    return dataclasses.replace(target, **changes) if changes else target


def run_target(target, steps, args):
    """Run the requested stages for one target; return its photometry table."""
    datadir = zc.data_dir(target.name)
    figdir = zc.fig_dir("photometry", target.name)

    print(f"\n{'=' * 72}\n{target.name}   {target.start_date} .. {target.end_date}"
          f"\n  designation : {target.query_designation}"
          f"{'  [fragments allowed]' if target.allow_fragment else ''}"
          f"\n  data        : {datadir}\n{'=' * 72}")

    urls = None
    if "query" in steps:
        eph, frames, report = zc.search_frames(target)
        print(report)
        if frames.empty:
            log.warning("%s: no frames contain the target; nothing to reduce", target.name)
            return None
        eph.to_csv(datadir / "eph.csv", index=False)
        urls = zc.build_urls(frames, is_cutout=target.query.is_cutout,
                             cutout_size=target.query.cutout_size)
        urls.to_csv(datadir / "ztf.csv", index=False)
        zc.save_urls(urls, datadir / "fits_urls.txt")

    if "download" in steps:
        source = urls
        if source is None:
            url_file = datadir / "fits_urls.txt"
            if not url_file.exists():
                log.warning("%s: no fits_urls.txt; run the query step first", target.name)
                return None
            source = url_file.read_text().split()
        report = zc.download_urls(source, datadir, overwrite=args.overwrite, repair=True)
        print(report)

    table = None
    if "phot" in steps:
        table = zc.run_photometry(target, datadir=datadir, phot_config=target.phot)
        if table is None or table.empty:
            log.warning("%s: no frames to reduce in %s", target.name, datadir)
            return None
        summarise(target, table)
    elif "figures" in steps or "profile" in steps:
        path = zc.photometry_path(target.name)
        if path.exists():
            table = pd.read_csv(path, dtype={"target": str})
        else:
            log.warning("%s: %s not found; run the phot step first", target.name, path)

    if "profile" in steps and table is not None and not table.empty:
        from ztfcomet import profile as profile_mod
        _, psum = profile_mod.run_profiles(target, table, datadir,
                                           plots=not args.no_profile_plots)
        if psum is not None and not psum.empty:
            good = psum[psum["quality_ok"].astype(bool)] if "quality_ok" in psum else psum
            print(f"  coma profile   : slope median {good['slope_comet'].median():+.2f} "
                  f"(stars {psum['slope_star'].median():+.2f}; steady-state coma = -1), "
                  f"{len(psum)} frames, excess @3 FWHM median {good['excess_at_3fwhm'].median():.0f}x")

    if "figures" in steps and table is not None and not table.empty:
        if not args.no_cutout_figures:
            zc.save_all_cutouts(table, datadir, figdir / "cutout", target=target,
                                phot_config=target.phot, dpi=args.dpi)
        make_lightcurve(target, table, figdir, only_good=args.only_good)

    return table


def summarise(target, table):
    """Print per-filter counts and the quality-flag tally."""
    print(f"\n  frames reduced : {len(table)}")
    counts = table["filter"].value_counts().to_dict()
    print(f"  by filter      : {counts}")

    flagged = {c.replace('flag_', ''): int(table[c].sum())
               for c in zc.FLAG_COLUMNS if c in table and table[c].any()}
    n_ok = int(table["quality_ok"].sum()) if "quality_ok" in table else len(table)
    print(f"  unflagged      : {n_ok}/{len(table)}")
    if flagged:
        print(f"  flags          : {flagged}")
        print("                   (flagged rows are kept in the table, not dropped)")

    if "contam_n_sources" in table and table["flag_contaminated"].any():
        hit = table[table["flag_contaminated"]]
        print(f"  contamination  : {len(hit)} frame(s) with Gaia sources in the aperture; "
              f"flux ratio {hit.contam_ratio.min():.2f}-{hit.contam_ratio.max():.2f}")

    good = table[table["quality_ok"]] if "quality_ok" in table else table
    if not good.empty and good["afrho0_cm"].notna().any():
        print(f"  A(0)frho (cm)  : {good['afrho0_cm'].min():.1f} .. {good['afrho0_cm'].max():.1f}"
              f"  (rho = {target.phot.rho_km:.0f} km)")


def make_lightcurve(target, table, figdir, only_good=False):
    """Save the Af-rho lightcurve, versus T-Tp where perihelion is known."""
    import matplotlib.pyplot as plt

    bands = [b for b in ("ZTF_g", "ZTF_r", "ZTF_i") if (table["filter"] == b).any()]
    if not bands:
        return

    x, tp = "rh", None
    if target.perihelion_jd:
        tp = list(target.perihelion_jd.values())[-1]
        x = "tp"

    ax = zc.plot_afrho({target.name: table}, filters=bands, x=x,
                       perihelion_jd={target.name: tp} if tp else None,
                       only_good=only_good)
    ax.set_title(f"{target.name}   " + r"$\rho$ = " + f"{target.phot.rho_km:.0f} km")
    outpath = zc.fig_kind_path("afrho", "lightcurve", target.name)
    ax.figure.savefig(outpath, dpi=200)
    plt.close(ax.figure)
    print(f"  lightcurve     : {outpath}")


def main(argv=None):
    args = parse_args(argv)
    zc.setup_logging(logging.WARNING if args.quiet else logging.INFO)

    steps = [s for s in args.steps if not (args.no_figures and s == "figures")]
    matplotlib.use("Agg")   # batch: never try to open a window

    print(zc.directory.describe())

    tables = {}
    for name in args.targets:
        target = build_target(name, args)
        try:
            table = run_target(target, steps, args)
        except Exception as exc:                                # noqa: BLE001
            log.exception("%s failed: %s", name, exc)
            continue
        if table is not None and not table.empty:
            tables[target.name] = table

    if len(tables) > 1 and "figures" in steps:
        import matplotlib.pyplot as plt
        ax = zc.plot_afrho(tables, filters=["ZTF_r"], x="rh", only_good=args.only_good)
        ax.set_title(r"$A(0\degree)f\rho$ — all targets (ZTF_r)")
        outpath = zc.fig_path("afrho", "all_targets.png")
        ax.figure.savefig(outpath, dpi=200)
        plt.close(ax.figure)
        print(f"\ncombined lightcurve: {outpath}")

    print(f"\nDone: {len(tables)}/{len(args.targets)} target(s) reduced.")
    return 0 if tables else 1


if __name__ == "__main__":
    raise SystemExit(main())
