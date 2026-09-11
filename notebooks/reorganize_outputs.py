#!/usr/bin/env python
"""Move ``results/`` and ``fig/`` from the per-target tree to the per-subject tree.

The layout used to be ``results/<target>/photometry_<target>.csv`` and
``fig/<target>/afrho_rh_<target>.png``; it is now ``results/photometry/<target>.csv``
and ``fig/afrho/<target>_rh.png`` -- see ``ztfcomet.directory``.  Idempotent:
a file already in place is left alone, and empty per-target directories are
removed afterwards.  Both trees are gitignored, so there is nothing for git to
track here.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import ztfcomet as zc

# (regex on the old file name, subject, new file-name template); <T> is the target slug
RESULT_RULES = [
    (r"photometry_(?P<t>[^/]+)\.csv$", "photometry", "{t}.csv"),
    (r"profile_summary_(?P<t>[^/]+)\.csv$", "profile", "{t}_summary.csv"),
    (r"profile_stars_(?P<t>[^/]+)\.csv$", "profile", "{t}_stars.csv"),
    (r"profile_resolution_(?P<t>[^/]+)\.csv$", "profile", "{t}_resolution.csv"),
    (r"profile_(?P<t>[^/]+)\.csv$", "profile", "{t}_profile.csv"),
    (r"afrho_(?P<t>[^/]+)\.csv$", "afrho", "{t}_afrho.csv"),
]
FIG_RULES = [
    (r"afrho_rh_(?P<t>[^/]+)\.png$", "afrho", "{t}_rh.png"),
    (r"afrho_apertures_(?P<t>[^/]+)\.png$", "afrho", "{t}_apertures.png"),
    (r"afrho_trend_(?P<t>[^/]+)\.png$", "afrho", "{t}_trend.png"),
    (r"afrho_color_(?P<t>[^/]+)\.png$", "afrho", "{t}_colour.png"),
    (r"afrho_(?P<t>[^/]+)\.png$", "afrho", "{t}_lightcurve.png"),
    (r"profile_summary_(?P<t>[^/]+)\.png$", "profile", "{t}_summary.png"),
    (r"profile_resolution_(?P<t>[^/]+)\.png$", "profile", "{t}_resolution.png"),
    (r"profile_validity_(?P<t>[^/]+)\.png$", "profile", "{t}_validity.png"),
    (r"profile_correction_(?P<t>[^/]+)\.png$", "profile", "{t}_correction.png"),
    (r"cutout_grid_(?P<t>[^/]+)\.png$", "photometry", "{t}_cutout_grid.png"),
    (r"diagnostics_(?P<t>[^/]+)\.png$", "photometry", "{t}_diagnostics.png"),
    (r"flagged_(?P<t>[^/]+)\.png$", "photometry", "{t}_flagged.png"),
]
SURVEY_RESULTS = {"profile_survey_fits.csv": ("profile", "survey_fits.csv"),
                  "profile_survey_ratio.csv": ("profile", "survey_ratio.csv"),
                  "profile_survey_targets.csv": ("profile", "survey_targets.csv"),
                  "profile_feasibility_survey.csv": ("profile", "feasibility.csv")}
ACTIVITY = {"afrho_trends.csv": "trends.csv", "afrho_peaks.csv": "peaks.csv", "elements.csv": "elements.csv"}
SURVEY_FIGS = {"profile_survey_slope.png": ("profile", "survey_slope.png"),
               "profile_survey_ratio.png": ("profile", "survey_ratio.png"),
               "afrho_trends_overview.png": ("afrho", "survey_overview.png"),
               "afrho_all_targets.png": ("afrho", "all_targets.png")}


def move(src: Path, dst: Path, log: list):
    if not src.exists() or dst.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    log.append((src, dst))


def apply_rules(path: Path, rules, root_fn, log):
    for pat, subject, template in rules:
        m = re.search(pat, path.name)
        if m:
            move(path, root_fn(subject) / template.format(t=m.group("t")), log)
            return True
    return False


def main():
    moved = []
    skip = set(zc.SUBJECTS) | {"survey", "activity"}
    # results/<target>/...
    for d in sorted(p for p in zc.RESULT_ROOT.iterdir() if p.is_dir() and p.name not in skip):
        for f in sorted(d.glob("*")):
            if f.is_file() and not f.name.startswith("."):
                apply_rules(f, RESULT_RULES, lambda s: zc.result_dir(s), moved)
    for name, (subject, new) in SURVEY_RESULTS.items():
        move(zc.RESULT_ROOT / name, zc.result_dir(subject) / new, moved)
    for name, new in ACTIVITY.items():
        move(zc.RESULT_ROOT / "activity" / name, zc.result_dir("afrho") / new, moved)
    # fig/<target>/...
    for d in sorted(p for p in zc.FIG_ROOT.iterdir() if p.is_dir() and p.name not in skip):
        t = d.name
        for f in sorted(d.glob("*.png")):
            if not apply_rules(f, FIG_RULES, lambda s: zc.fig_dir(s), moved) and f.name.startswith("cutout_ztf_"):
                move(f, zc.fig_dir("photometry", t) / f.name, moved)
        for f in sorted((d / "profile").glob("*.png")):
            move(f, zc.fig_dir("profile", t) / f.name, moved)
        for f in sorted((d / "cutout").glob("*.png")):
            move(f, zc.fig_dir("photometry", t) / "cutout" / f.name, moved)
    for name, (subject, new) in SURVEY_FIGS.items():
        move(zc.FIG_ROOT / "survey" / name, zc.fig_dir(subject) / new, moved)
        move(zc.FIG_ROOT / name, zc.fig_dir(subject) / new, moved)
    # fig/afrho/<T>_<kind>.png -> fig/afrho/<kind>/<T>.png
    for f in sorted(zc.fig_dir("afrho", create=False).glob("*_*.png")):
        m = re.match(r"(?P<t>.+)_(?P<k>rh|apertures|trend|lightcurve|colour)\.png$", f.name)
        if m:
            move(f, zc.fig_kind_path("afrho", m.group("k"), m.group("t")), moved)
    # prune what is now empty
    pruned = 0
    for root in (zc.RESULT_ROOT, zc.FIG_ROOT):
        for d in sorted((p for p in root.rglob("*") if p.is_dir()), key=lambda p: -len(p.parts)):
            if d.name in zc.SUBJECTS and d.parent == root:
                continue
            leftovers = [f for f in d.iterdir() if not f.name.startswith(".")]
            if not leftovers:
                for junk in d.iterdir():
                    junk.unlink()
                d.rmdir(); pruned += 1
    print(f"moved {len(moved)} files, removed {pruned} empty directories")
    for root in (zc.RESULT_ROOT, zc.FIG_ROOT):
        stray = [p for p in root.iterdir() if p.is_dir() and p.name not in zc.SUBJECTS]
        if stray:
            print(f"  left behind under {root.name}/: " + ", ".join(p.name for p in stray))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
