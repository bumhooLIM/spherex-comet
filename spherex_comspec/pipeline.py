"""
Orchestration: grouping once, then continuum subtraction and fitting per variant.

The grouping is computed **once** and shared by every variant.  It depends on
r_h, time and which bands each exposure sampled -- none of which the flag policy
or the flux space touches -- so holding it fixed is what makes the variants a
controlled comparison: the same groups, measured under different selections.

No figure is drawn here.  :mod:`plotting` is imported only by ``main.py --figures``
and by the notebooks.
"""

from __future__ import annotations

import json
import time
import warnings
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import directory as _dir
from .config import (EMISSION_DTYPES, MAIN_VARIANT, ApertureConfig, FitConfig, GroupingConfig,
                     ModelParams, Variant)
from .continuum import insufficient, n_in_emission, process_group
from .dataio import (PhaseAssignment, aperture_table, attach_afrho_ztf, list_catalog, list_targets, load_apphot,
                     load_fit_input, save_emission, save_fit_lines, save_fit_table,
                     select_spectrum)
from .fitting import fit_production_rates, model_curves
from .revisions import revision_for
from .grouping import regroup_all
from .logging_utils import get_logger, utcnow_iso

__all__ = ["run_grouping", "run_variant", "VariantResult", "choose_apertures"]

log = get_logger("pipeline")


# ------------------------------------------------------------------------------ grouping
def run_grouping(targets: Optional[Sequence[str]] = None,
                 cfg: Optional[GroupingConfig] = None, write: bool = True):
    """
    Regroup the catalog and persist the assignment.

    Returns
    -------
    assignment, group_map, cuts, epochs
        See :func:`grouping.regroup_all`.
    """
    cfg = cfg or GroupingConfig()
    assignment, group_map, cuts, epochs = regroup_all(targets, cfg)
    group_map = attach_afrho_ztf(group_map)                  # ZTF dust context, when the table exists
    if write:
        p1, p2 = PhaseAssignment.save(assignment, group_map)
        _dir.ensure_dirs()
        cuts.to_csv(_dir.RESULT_DIR / "phase_cuts.csv", index=False)
        log.info("wrote %s (%d exposures) and %s (%d groups)", p1.name, len(assignment),
                 p2.name, len(group_map))
    return assignment, group_map, cuts, epochs


def choose_apertures(targets: Sequence[str], cfg: ApertureConfig,
                     assignment: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """
    The aperture of every (target, phase), with the evidence that justified it
    (``dataio.aperture_table``).  Without an assignment each target is one phase (``-1``).
    """
    tables = [aperture_table(t, cfg, assignment) for t in targets]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


# --------------------------------------------------------------------------- one variant
@dataclass
class VariantResult:
    variant: Variant
    apertures: pd.DataFrame
    skipped: pd.DataFrame
    cont_summary: pd.DataFrame
    fits: pd.DataFrame
    not_fitted: pd.DataFrame
    elapsed_s: float = 0.0
    paths: Dict[str, Path] = field(default_factory=dict)


def _clear(variant: str) -> int:
    """Remove a variant's previous products: file names encode (target, aperture, phase), so a
    stale file from an earlier grouping would otherwise survive beside the new one."""
    d = _dir.variant_dirs(variant)
    n = 0
    for pat, key in (("*.csv", "emission"), ("*.csv", "lines"), ("gas_fit*.csv", "results"),
                     ("*.json", "results")):
        for f in d[key].glob(pat):
            f.unlink()
            n += 1
    return n


def run_variant(variant: Variant, targets: Optional[Sequence[str]] = None,
                assignment: Optional[pd.DataFrame] = None, *, clear: bool = True,
                write: bool = True, progress_every: int = 10) -> VariantResult:
    """
    Continuum subtraction and production-rate fitting for every target under one variant.

    Parameters
    ----------
    variant : Variant
        Flag policy, flux space and all parameters.
    targets : sequence of str, optional
        Defaults to every table in ``APPHOT_DIR``.
    assignment : DataFrame, optional
        The phase assignment; loaded from disk when omitted.
    clear : bool
        Remove the variant's previous outputs first.
    write : bool
        Persist emission files, fit tables and model curves.

    Returns
    -------
    VariantResult
    """
    t0 = time.time()
    targets = list(targets) if targets else list_targets()
    assignment = assignment if assignment is not None else PhaseAssignment.load()
    _dir.ensure_dirs(variant.name)
    if clear and write:
        n = _clear(variant.name)
        if n:
            log.info("[%s] cleared %d stale file(s)", variant.name, n)
    log.info("[%s] %s | flux space %s | %d target(s)", variant.name, variant.flags.describe(),
             "distcorr" if variant.use_distcorr else "physical", len(targets))

    # ---- apertures ------------------------------------------------------------
    apertures = choose_apertures(targets, variant.aperture, assignment)
    space = "distcorr" if variant.use_distcorr else "physical"

    # ---- continuum subtraction ---------------------------------------------------
    # One aperture per phase (2026-09-14): the rows of a target are read once, and each
    # distinct aperture label selects its own spectrum for the phases assigned to it.
    skipped, summaries = [], []
    by_target = list(apertures.groupby("target", sort=False))
    for k, (t, ap_t) in enumerate(by_target, start=1):
        df = load_apphot(t)
        a_t = assignment[assignment.target == t]
        target_sums, target_pts = [], []
        for lab, ap_rows in ap_t.groupby("ap_label", sort=False):
            r_ap_km = float(ap_rows.r_ap_km.iloc[0])
            phases = set(int(x) for x in ap_rows.phase)
            spec = select_spectrum(df, lab, variant, a_t)
            rej = spec.attrs["rejection"]
            spec = spec[spec.phase.isin(phases)]
            # A phase that the flag policy (or the aperture's coverage) emptied entirely never
            # reaches the groupby below; it is still a skipped group and must be counted as one,
            # or the census under-reports the cut.
            for ph in sorted(phases - set(spec.phase.unique())):
                skipped.append(dict(target=t, r_ap_km=r_ap_km, phase=int(ph), n_points=0,
                                    n_emission=0, reason="0 usable points (all rows rejected)"))
            for ph, raw in spec.groupby("phase", sort=True):
                raw = raw.reset_index(drop=True)
                why = insufficient(raw, variant.continuum)
                if why is not None:
                    skipped.append(dict(target=t, r_ap_km=r_ap_km, phase=int(ph), n_points=len(raw),
                                        n_emission=n_in_emission(raw), reason=why))
                    continue
                rev = revision_for(t, ph, variant.revisions)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    out = process_group(raw, t, r_ap_km, int(ph), int(raw.epoch.iloc[0]),
                                        variant.continuum, space, rev=rev)
                summ = out["summary"]
                summ["n_rejected_badphot"] = rej["n_badphot"]
                summ["n_rejected_flag"] = rej["n_flag"]
                target_sums.append(summ)
                target_pts.append(out["points"])
        if target_sums:
            summaries.append(pd.concat(target_sums, ignore_index=True))
            if write:
                save_emission(variant.name, summaries[-1], pd.concat(target_pts, ignore_index=True))
        if progress_every and (k % progress_every == 0 or k == len(by_target)):
            log.info("[%s] continuum %d/%d  (%d band-rows, %d groups skipped)", variant.name, k,
                     len(by_target), sum(len(s) for s in summaries), len(skipped))

    cont = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    skipped = pd.DataFrame(skipped)

    # ---- fitting ------------------------------------------------------------------
    rows, not_fitted = [], []
    if len(cont):
        catalog = _catalog_from(cont)
        for i, r in enumerate(catalog.itertuples(index=False), start=1):
            params = replace(variant.model, rho_ap_km=float(r.r_ap_km))
            try:
                pts = load_fit_input(variant.name, r.target, r.r_ap_km, int(r.phase), variant.fit) \
                    if write else _fit_input_from(cont, r, variant.fit)
            except (ValueError, FileNotFoundError) as exc:
                not_fitted.append(dict(target=r.target, phase=int(r.phase),
                                       reason="no band survived continuum validation"
                                       if "no usable band" in str(exc) else str(exc)[:100]))
                continue
            try:
                rev = revision_for(r.target, r.phase, variant.revisions)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = fit_production_rates(pts, params, variant.fit, space="physical", rev=rev)
                    curves = model_curves(fit, pts, params)
                if write:
                    save_fit_lines(variant.name, fit, pts, curves)
                row = fit.to_row()
                row["epoch"] = int(pts.epoch.iloc[0]) if "epoch" in pts else -1
                row["n_flag_a"] = int((pts.sourceflag.astype(str) == "a").sum())
                row["n_flag_b"] = int((pts.sourceflag.astype(str) == "b").sum())
                rows.append(row)
            except Exception as exc:                        # noqa: BLE001 - recorded, not hidden
                not_fitted.append(dict(target=r.target, phase=int(r.phase), reason=repr(exc)[:120]))
    fits = attach_afrho_ztf(pd.DataFrame(rows), jd_col="jd_utc_mean")   # ZTF dust context, when the table exists
    not_fitted = pd.DataFrame(not_fitted)

    paths = {}
    if write:
        paths["gas_fit"] = save_fit_table(variant.name, fits)
        paths["continuum"] = save_fit_table(variant.name, cont, "continuum_summary.csv")
        paths["skipped"] = save_fit_table(variant.name, skipped, "skipped_groups.csv")
        paths["not_fitted"] = save_fit_table(variant.name, not_fitted, "not_fitted.csv")
        paths["apertures"] = save_fit_table(variant.name, apertures, "apertures.csv")
        meta = dict(variant=variant.to_dict(), hash=variant.hash, written_utc=utcnow_iso(),
                    n_targets=len(targets), n_groups_analysed=int(cont.groupby(["target", "phase"]).ngroups) if len(cont) else 0,
                    n_groups_skipped=len(skipped), n_band_rows=len(cont), n_fitted=len(fits),
                    n_not_fitted=len(not_fitted), elapsed_s=round(time.time() - t0, 1))
        p = _dir.variant_dirs(variant.name)["results"] / "run.meta.json"
        p.write_text(json.dumps(meta, indent=2, default=str))
        paths["meta"] = p

    el = time.time() - t0
    log.info("[%s] done: %d groups analysed, %d skipped, %d band-rows, %d fitted, %d not "
             "fittable, %.1f s", variant.name,
             int(cont.groupby(["target", "phase"]).ngroups) if len(cont) else 0, len(skipped),
             len(cont), len(fits), len(not_fitted), el)
    return VariantResult(variant, apertures, skipped, cont, fits, not_fitted, el, paths)


def _catalog_from(cont: pd.DataFrame) -> pd.DataFrame:
    g = cont.groupby(["target", "r_ap_km", "phase"], sort=True)
    return g.size().reset_index()[["target", "r_ap_km", "phase"]]


def _fit_input_from(cont, r, cfg):      # pragma: no cover - only for write=False
    raise ValueError("in-memory fit input requires write=True in this version")
