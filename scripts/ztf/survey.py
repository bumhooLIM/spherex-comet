#!/usr/bin/env python
"""Batch ZTF survey: query, download, verify and reduce a list of comets.

Runs the full chain for every comet in a spreadsheet, unattended.  Written to
survive a long unsupervised run: each target is independent, failures are
recorded rather than fatal, and the run is resumable.

Per target
----------
1. Resolve the designation against Horizons (fragment-aware) and search IRSA.
2. Choose the cutout size from the target's own ephemeris -- the small box
   unless the comet is bright or nearby (see :func:`ztfcomet.choose_cutout_size`).
3. Download with FITS validation, then **re-verify every requested URL against
   what is on disk** before anything downstream runs.
4. Multi-aperture Af-rho photometry, apertures filtered by the scale test.
5. Radial surface-brightness profile of the comet against field stars on the
   same frame -- is the coma extended, and does it fall as 1/rho?
6. Figures: Af-rho vs r_h - q (clean frames only), aperture comparison,
   profile summary.

A target whose ephemeris never reaches ``vmag_max`` gets an empty directory and
a ``no_epochs`` status -- that is a result, not a failure.

Outputs, all under the data root
--------------------------------
``survey_log_<stamp>.log``     full procedure log
``survey_status.csv``          one row per target, machine-readable
``survey_summary.md``          human-readable summary written at the end

Examples
--------
    python notebooks/survey.py                        # all 68, resumable
    python notebooks/survey.py --targets 24P 2P       # a subset
    python notebooks/survey.py --no-resume            # redo everything
    python notebooks/survey.py --steps phot figures    # reduce what is on disk
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import logging
import os
import socket
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import numpy as np
import ztfcomet as zc
from ztfcomet import orbit
from ztfcomet import profile as profile_mod

log = logging.getLogger("ztfcomet.survey")

# Last-resort guard for an unattended run.  astroquery sets its own 30 s
# timeout and the downloader sets 120 s, but a socket that goes away silently
# can still leave a read blocked with no timeout of its own: the run stalled
# 54 minutes on one Horizons call at 0% CPU while the service answered other
# clients in under a second.  A default timeout bounds every socket that does
# not set one, so a dead peer costs minutes instead of the rest of the night.
socket.setdefaulttimeout(300)

STEPS = ("query", "download", "phot", "profile", "figures", "cutouts")
DEFAULT_LIST = "doc/sx_comet_list_ver2607.xlsx"


# --------------------------------------------------------------------- helpers
def read_designations(path):
    """Comet designations from the survey spreadsheet.

    Ignores blank rows and any stray numeric cell (the sheet carries the row
    count in the last populated cell of the designation column).
    """
    df = pd.read_excel(path)
    column = "desig" if "desig" in df.columns else df.columns[1]
    out = []
    for value in df[column].dropna():
        text = str(value).strip()
        if text and not text.isdigit():
            out.append(text)
    return out


def setup_logging(logfile, verbose=True):
    """Log to both the run file and the console."""
    root = logging.getLogger("ztfcomet")
    root.handlers[:] = []
    root.setLevel(logging.INFO if verbose else logging.WARNING)
    root.propagate = False

    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s",
                            "%Y-%m-%d %H:%M:%S")
    handler = logging.FileHandler(logfile, encoding="utf-8")
    handler.setFormatter(fmt)
    root.addHandler(handler)

    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter("%(levelname)-7s %(message)s"))
    root.addHandler(console)
    return root


def status_path(root):
    return Path(root) / "survey_status.csv"


def load_status(root):
    path = status_path(root)
    if not path.exists():
        return {}
    try:
        # "2024E1" parses as the float 20240.0 unless the column is pinned to
        # str -- which silently breaks resume for any \d{4}E\d designation.
        df = pd.read_csv(path, dtype={"target": str, "designation": str})
        return {str(r["target"]): dict(r) for _, r in df.iterrows()}
    except Exception:                                           # noqa: BLE001
        log.warning("Could not read %s; starting a fresh status table", path)
        return {}


def save_status(root, status):
    if not status:
        return
    df = pd.DataFrame(list(status.values()))
    lead = [c for c in ("target", "status", "n_frames", "n_downloaded",
                        "n_missing", "n_corrupt", "cutout_size") if c in df.columns]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df.to_csv(status_path(root), index=False)


def is_done(entry):
    """A target needs no rework when it completed or has nothing to fetch."""
    return str(entry.get("status", "")) in {"ok", "no_epochs", "no_frames"}


# ------------------------------------------------------------------ one target
def process(designation, args, root, prior=None):
    """Run the pipeline for one comet.  Returns a status dict.

    *prior* is the target's existing status entry, if any; a post-hoc pass
    (``--steps profile figures``) starts from it so the query and download
    counts recorded by the full run are preserved rather than overwritten.
    """
    started = time.time()
    target = zc.get_target(designation)
    if args.vmag_max is not None or args.interval_days is not None:
        overrides = {}
        if args.vmag_max is not None:
            overrides["vmag_max"] = args.vmag_max
        if args.interval_days is not None:
            overrides["interval_days"] = args.interval_days
        target = dataclasses.replace(
            target, query=dataclasses.replace(target.query, **overrides))
    if args.start:
        target = dataclasses.replace(target, start_date=args.start)
    # --end defaults to today, as its help text says.  Applying it only when
    # given let the per-target end_date in config.py silently win: 2024E1
    # carried 2025-10-30, so a run asking for "2025-03-01 .. present" queried
    # ten months less than that -- including the comet's 2026-01-20 perihelion
    # and every post-perihelion epoch.
    target = dataclasses.replace(target, end_date=args.end or zc.config.today())

    entry = dict(prior) if prior else {}
    entry.update(target=target.name, designation=designation, status="running",
                 started=datetime.now(timezone.utc).isoformat(timespec="seconds"))

    datadir = zc.data_dir(target.name)          # created even when empty
    figdir = zc.fig_dir("afrho")
    entry["datadir"] = str(datadir)

    log.info("=" * 70)
    log.info("%s  (designation %r)  %s .. %s  Vmag<%g",
             target.name, target.query_designation, target.start_date,
             target.end_date, target.query.vmag_max)

    urls = None
    # ---------------------------------------------------------------- query
    if "query" in args.steps:
        eph, frames, report = zc.search_frames(target)
        entry.update(n_eph_steps=report.n_eph_steps, n_eph_kept=report.n_eph_kept,
                     n_steps_failed=report.n_steps_failed,
                     n_frames_found=report.n_frames_found,
                     n_frames=report.n_frames_containing_target)
        log.info("%s: %d/%d epochs pass the cuts, %d frames contain the target",
                 target.name, report.n_eph_kept, report.n_eph_steps,
                 report.n_frames_containing_target)
        if report.n_steps_failed:
            log.warning("%s: %d query steps FAILED -- coverage incomplete",
                        target.name, report.n_steps_failed)

        if eph is not None and not eph.empty:
            eph.to_csv(datadir / "eph.csv", index=False)
        if report.n_eph_kept == 0:
            # Say which cut did it: C/2014 UN271 at Tmag 16 was reported as
            # "never brighter than Vmag 20" when rh_max = 10 au had removed it.
            why = (f"r_h >= {target.query.rh_max:g} au on every epoch" if report.n_cut_rh == report.n_eph_steps
                   else f"never brighter than Vmag {target.query.vmag_max:g}" if report.n_cut_vmag and not report.n_cut_rh
                   else f"{report.n_cut_rh} epochs beyond r_h {target.query.rh_max:g} au, "
                        f"{report.n_cut_vmag} fainter than Vmag {target.query.vmag_max:g}")
            entry.update(status="no_epochs", note=why)
            log.info("%s: no epoch passes the cuts (%s) -- leaving the directory empty", target.name, why)
            return _finish(entry, started)
        if frames.empty:
            entry.update(status="no_frames", note="epochs pass the cuts but ZTF has no coverage")
            log.info("%s: no ZTF frames contain the target", target.name)
            return _finish(entry, started)

        size, reason = zc.choose_cutout_size(frames, target.query)
        entry.update(cutout_size=size, cutout_reason=reason)
        log.info("%s: cutout %s (%s)", target.name, size, reason)

        urls = zc.build_urls(frames, is_cutout=target.query.is_cutout, cutout_size=size)
        urls.to_csv(datadir / "ztf.csv", index=False)
        zc.save_urls(urls, datadir / "fits_urls.txt")

    # ------------------------------------------------------------- download
    if "download" in args.steps:
        if urls is None:
            url_file = datadir / "fits_urls.txt"
            if not url_file.exists():
                entry.update(status="no_urls", note="run the query step first")
                return _finish(entry, started)
            urls = [u for u in url_file.read_text().split() if u]
        report = zc.download_urls(urls, datadir, overwrite=args.overwrite, repair=True,
                                  progress=False)
        entry.update(n_downloaded=report.n_downloaded, n_existing=report.n_skipped_existing,
                     n_dl_failed=report.n_failed + report.n_rejected_not_fits)
        log.info("%s: downloaded %d, existing %d, failed %d",
                 target.name, report.n_downloaded, report.n_skipped_existing,
                 report.n_failed + report.n_rejected_not_fits)

        # Completeness is re-checked against the disk, not inferred from the
        # download report: a file can be truncated or removed afterwards.
        check = zc.verify_downloads(urls, datadir)
        entry.update(n_expected=check["n_expected"], n_valid=check["n_valid"],
                     n_missing=check["n_missing"], n_corrupt=check["n_corrupt"],
                     complete=check["complete"])
        if not check["complete"]:
            log.warning("%s: INCOMPLETE -- %d missing, %d corrupt of %d expected",
                        target.name, check["n_missing"], check["n_corrupt"],
                        check["n_expected"])
            for name in (check["missing"] + check["corrupt"])[:5]:
                log.warning("    %s", name)
        else:
            log.info("%s: all %d files verified as valid FITS",
                     target.name, check["n_valid"])
        if check["n_valid"] == 0:
            entry.update(status="no_data", note="nothing valid on disk")
            return _finish(entry, started)

    # ----------------------------------------------------------- photometry
    table = None
    if "phot" in args.steps:
        table, skipped = zc.run_multi_aperture(target, datadir=datadir, progress=False)
        if table is None or table.empty:
            entry.update(status="no_photometry", note="no aperture passed the scale test")
            return _finish(entry, started)
        entry["aperture_skips"] = json.dumps({int(k): int(v) for k, v in skipped.items()})
        entry["n_rows"] = len(table)
        entry["n_apertures"] = int(table["rho_km"].nunique())
        entry["n_clean"] = int(table["quality_ok"].sum())
        for flag in zc.FLAG_COLUMNS:
            if flag in table:
                entry[flag] = int(table[flag].sum())
        log.info("%s: %d rows over %d apertures, %d clean",
                 target.name, len(table), entry["n_apertures"], entry["n_clean"])
    elif any(k in args.steps for k in ("figures", "profile", "cutouts")):
        path = zc.photometry_path(target.name, create=False)
        if path.exists():
            table = pd.read_csv(path, dtype={"target": str})

    elements = None
    if table is not None and not table.empty and ("profile" in args.steps or "figures" in args.steps):
        elements = fetch_elements_for(target, table)

    # -------------------------------------------------------- radial profile
    if "profile" in args.steps and table is not None and not table.empty:
        try:
            _, psum = profile_mod.run_profiles(
                target, table, datadir, progress=False,
                plots=not args.no_profile_plots, elements=elements)
            if psum is not None and not psum.empty:
                good = psum[psum["quality_ok"].astype(bool)] if "quality_ok" in psum else psum
                entry.update(
                    profile_frames=len(psum),
                    profile_stars_median=float(psum["n_stars"].median()),
                    profile_slope_median=float(good["slope_comet"].median()) if len(good) else np.nan,
                    profile_slope_p16=float(good["slope_comet"].quantile(0.16)) if len(good) else np.nan,
                    profile_slope_p84=float(good["slope_comet"].quantile(0.84)) if len(good) else np.nan,
                    profile_star_slope_median=float(psum["slope_star"].median()),
                    profile_excess_median=float(good["excess_at_3fwhm"].median()) if len(good) else np.nan,
                )
                log.info("%s: coma slope median %+.2f (stars %+.2f) over %d clean frames",
                         target.name, entry["profile_slope_median"],
                         entry["profile_star_slope_median"], len(good))
        except Exception as exc:                                # noqa: BLE001
            log.warning("%s: radial profiles failed: %s", target.name, exc)
            entry["profile"] = f"failed: {exc}"

    # -------------------------------------------------------------- figures
    if "cutouts" in args.steps and table is not None and not table.empty:
        # One PNG per frame, drawn with the smallest aperture that passed the
        # scale test on that target; the annotated cutouts live with the
        # photometry, under fig/photometry/<target>/cutout/.
        ref = table[np.isclose(table["rho_km"], table["rho_km"].min())]
        written = zc.save_all_cutouts(ref, datadir, zc.fig_dir("photometry", target.name) / "cutout",
                                      target=target, dpi=50, progress=False)
        entry["n_cutout_figures"] = len(written)
        log.info("%s: %d cutout figures", target.name, len(written))
    if "figures" in args.steps and table is not None and not table.empty:
        try:
            make_figures(target, table, figdir, elements=elements)
            entry["figures"] = "ok"
        except Exception as exc:                                # noqa: BLE001
            log.warning("%s: figures failed: %s", target.name, exc)
            entry["figures"] = f"failed: {exc}"

    entry["status"] = "ok"
    return _finish(entry, started)


def _finish(entry, started):
    entry["elapsed_s"] = round(time.time() - started, 1)
    entry["finished"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return entry


def fetch_elements_for(target, table):
    """Perihelion elements osculating near the observations, or None."""
    try:
        record = target.resolve_orbit_record(float(table["obsjd"].median()))
        return orbit.fetch_elements(record, epoch_jd=float(table["obsjd"].median()))
    except Exception as exc:                                    # noqa: BLE001
        log.warning("%s: no orbital elements (%s); plotting plain r_h", target.name, exc)
        return None


def make_figures(target, table, figdir, elements=None):
    """Af-rho vs r_h - q (clean only), aperture comparison, diagnostics."""
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    slug = zc.target_slug(target.name)
    bands = [b for b in ("ZTF_g", "ZTF_r", "ZTF_i") if (table["filter"] == b).any()]
    if not bands:
        return

    clean = table[table["quality_ok"]] if "quality_ok" in table else table
    if clean.empty:
        log.warning("%s: every frame is flagged; no science figure", target.name)
        return

    # Reference aperture: the one with the most clean measurements.
    counts = clean.groupby("rho_km").size()
    rho_ref = float(counts.idxmax())

    ax = zc.plot_afrho_vs_rh({target.name: table}, rho_km=rho_ref, filters=bands,
                             only_good=True, elements=elements)
    ax.set_title(f"{target.name}   " + r"$\rho$ = " + f"{rho_ref:.0f} km   (clean frames only)",
                 pad=34 if elements is not None else 12)
    ax.figure.tight_layout()
    ax.figure.savefig(zc.fig_kind_path("afrho", "rh", target.name), dpi=200)
    plt.close(ax.figure)

    if table["rho_km"].nunique() > 1:
        band = "ZTF_r" if "ZTF_r" in bands else bands[0]
        ax = zc.plot_afrho_apertures(table, filters=[band], only_good=True,
                                     elements=elements)
        ax.set_title(f"{target.name} — apertures ({band}, clean)",
                     pad=34 if elements is not None else 12)
        ax.figure.tight_layout()
        ax.figure.savefig(zc.fig_kind_path("afrho", "apertures", target.name), dpi=200)
        plt.close(ax.figure)


# --------------------------------------------------------------------- driver
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--list", default=DEFAULT_LIST, help="spreadsheet of designations")
    p.add_argument("--targets", nargs="+", help="override the list with these designations")
    p.add_argument("--steps", nargs="+", choices=STEPS, default=list(STEPS))
    p.add_argument("--vmag-max", type=float, default=20.0)
    p.add_argument("--interval-days", type=float)
    p.add_argument("--start", default="2025-03-01")
    p.add_argument("--end", help="default: today")
    p.add_argument("--no-resume", action="store_true", help="reprocess completed targets")
    p.add_argument("--overwrite", action="store_true", help="re-download existing files")
    p.add_argument("--limit", type=int, help="process at most this many targets")
    p.add_argument("--no-profile-plots", action="store_true",
                   help="skip the per-frame radial-profile PNGs (summary still made)")
    p.add_argument("--allow-fallback-root", action="store_true",
                   help="proceed even if the SSD data root is unavailable")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    root = zc.directory.DATA_ROOT

    # Refuse to start on the in-project fallback. The SSD holds both the data
    # and survey_status.csv, so running without it silently restarts the whole
    # survey in the wrong place -- which is exactly what happened once when the
    # volume was briefly unresolvable at import.
    # A pinned root that is not there is not a fallback, so the check below
    # would let the run start -- and an external SSD that has dropped off the
    # bus (twice, with I/O errors, during long runs) leaves exactly that state.
    pinned = os.environ.get("ZTFCOMET_DATA")
    if pinned and not Path(pinned).is_dir():
        sys.stderr.write(f"REFUSING TO START: ZTFCOMET_DATA={pinned} is not a directory. "
                         "Mount the drive, or unset the variable to use the local data/ root.\n")
        return 2
    if zc.directory.using_fallback_root() and not args.allow_fallback_root:
        sys.stderr.write(
            f"REFUSING TO START: data root resolved to {root}, the in-project\n"
            f"fallback, not {zc.directory.SSD_DATA_ROOT}.\n"
            f"The SSD holds the existing data and the resume checkpoint.\n"
            f"Mount it, or set ZTFCOMET_DATA, or pass --allow-fallback-root.\n")
        return 2

    root.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    logfile = root / f"survey_log_{stamp}.log"
    setup_logging(logfile)

    designations = args.targets or read_designations(args.list)
    if args.limit:
        designations = designations[: args.limit]

    log.info("ZTF comet survey — %d targets", len(designations))
    log.info("window %s .. %s, Vmag < %g", args.start, args.end or "today", args.vmag_max)
    log.info("data root %s", root)
    log.info("log file  %s", logfile)

    status = {} if args.no_resume else load_status(root)
    n_skipped = 0

    for i, designation in enumerate(designations, 1):
        name = zc.target_slug(designation)
        if not args.no_resume and name in status and is_done(status[name]):
            n_skipped += 1
            log.info("[%d/%d] %s — already done (%s), skipping",
                     i, len(designations), name, status[name].get("status"))
            continue

        log.info("[%d/%d] %s", i, len(designations), designation)
        try:
            entry = process(designation, args, root, prior=status.get(name))
        except KeyboardInterrupt:
            log.warning("interrupted by user")
            save_status(root, status)
            raise
        except Exception as exc:                                # noqa: BLE001
            log.error("%s FAILED: %s", designation, exc)
            log.debug(traceback.format_exc())
            entry = {"target": name, "designation": designation, "status": "error",
                     "note": str(exc)[:300],
                     "finished": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        status[entry["target"]] = entry
        save_status(root, status)          # checkpoint after every target
        log.info("[%d/%d] %s -> %s (%.0f s)", i, len(designations), entry["target"],
                 entry.get("status"), entry.get("elapsed_s", 0))

    write_summary(root, status, n_skipped, args)
    log.info("Survey complete. Log: %s", logfile)
    return 0


def write_summary(root, status, n_skipped, args):
    """Write survey_summary.md: process counts and flag statistics."""
    if not status:
        return
    df = pd.DataFrame(list(status.values()))
    lines = [
        "# ZTF comet survey — summary", "",
        f"Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"Window {args.start} .. {args.end or 'today'}, Vmag < {args.vmag_max:g}",
        f"Data root `{root}`", "",
        "## Targets", "",
        f"- processed: **{len(df)}**" + (f" ({n_skipped} skipped as already done)" if n_skipped else ""),
        "",
        "| status | targets | meaning |",
        "|---|---|---|",
    ]
    meaning = {
        "ok": "reduced end to end",
        "no_epochs": "never brighter than the Vmag cut — directory left empty",
        "no_frames": "bright enough, but ZTF has no coverage",
        "no_data": "nothing valid on disk after download",
        "no_photometry": "no aperture passed the scale test",
        "error": "failed — see the log",
    }
    for name, n in df["status"].value_counts().items():
        lines.append(f"| `{name}` | {n} | {meaning.get(name, '')} |")

    ok = df[df["status"] == "ok"] if "status" in df else df.iloc[:0]
    if len(ok):
        def total(column):
            """Column sum, tolerant of a column that no run produced."""
            if column not in ok.columns:
                return 0
            return int(pd.to_numeric(ok[column], errors="coerce").fillna(0).sum())

        lines += ["", "## Data", "",
                  f"- frames containing a target: **{total('n_frames')}**",
                  f"- files downloaded this run: **{total('n_downloaded')}**",
                  f"- verified valid on disk: **{total('n_valid')}**",
                  f"- missing: {total('n_missing')}   corrupt: {total('n_corrupt')}",
                  f"- photometry rows (frame x aperture): **{total('n_rows')}**",
                  f"- clean rows: **{total('n_clean')}**"]
        if "cutout_size" in ok:
            lines += ["", "| cutout | targets |", "|---|---|"]
            for size, n in ok["cutout_size"].value_counts().items():
                lines.append(f"| {size} | {n} |")

        flags = [f for f in zc.FLAG_COLUMNS if f in ok.columns]
        if flags:
            n_rows = total("n_rows") or 1
            lines += ["", "## Flag statistics", "",
                      "Rows are (frame x aperture). Flagged rows are kept in the",
                      "tables; only the science figures exclude them.", "",
                      "| flag | rows | % | kind |", "|---|---|---|---|"]
            for f in sorted(flags, key=lambda c: -total(c)):
                kind = "critical" if f in zc.CRITICAL_FLAGS else "advisory"
                lines.append(f"| `{f.replace('flag_', '')}` | {total(f)} | "
                             f"{100 * total(f) / n_rows:.1f}% | {kind} |")

    if "profile_slope_median" in df.columns:
        prof = df[np.isfinite(pd.to_numeric(df["profile_slope_median"], errors="coerce"))]
        if len(prof):
            coma = prof[prof["profile_slope_median"].between(-1.3, -0.7)]
            lines += ["", "## Coma radial profiles", "",
                      "Power-law slope of surface brightness vs radius outside the PSF core,",
                      "median over clean frames. A steady-state coma gives -1; field stars on",
                      "the same frames give about -4.", "",
                      f"- targets profiled: **{len(prof)}**",
                      f"- median slope within [-1.3, -0.7] (steady-state-like): **{len(coma)}**", "",
                      "| target | frames | slope (16-84%) | stars | excess @3 FWHM |",
                      "|---|---|---|---|---|"]
            for _, r in prof.sort_values("profile_slope_median").iterrows():
                lines.append(
                    f"| {r['target']} | {int(r.get('profile_frames', 0))} | "
                    f"{r['profile_slope_median']:+.2f} ({r.get('profile_slope_p16', np.nan):+.2f} .. "
                    f"{r.get('profile_slope_p84', np.nan):+.2f}) | "
                    f"{r.get('profile_star_slope_median', np.nan):+.1f} | "
                    f"{r.get('profile_excess_median', np.nan):.0f}x |")

    incomplete = df[df.get("complete") == False] if "complete" in df else df.iloc[:0]  # noqa: E712
    if len(incomplete):
        lines += ["", "## Incomplete downloads", "",
                  "| target | expected | valid | missing | corrupt |", "|---|---|---|---|---|"]
        for _, r in incomplete.iterrows():
            lines.append(f"| {r['target']} | {r.get('n_expected')} | {r.get('n_valid')} | "
                         f"{r.get('n_missing')} | {r.get('n_corrupt')} |")

    errors = df[df["status"] == "error"] if "status" in df else df.iloc[:0]
    if len(errors):
        lines += ["", "## Errors", ""]
        for _, r in errors.iterrows():
            lines.append(f"- **{r['target']}** — {r.get('note', '')}")

    (Path(root) / "survey_summary.md").write_text("\n".join(lines) + "\n")
    log.info("Wrote %s", Path(root) / "survey_summary.md")


if __name__ == "__main__":
    raise SystemExit(main())
