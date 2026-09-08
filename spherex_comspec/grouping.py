"""
Subdividing the 28-day ``epoch`` groups into ``phase`` groups of one physical state.

Reproduces ``notebooks/phase_group_update.ipynb``.  The 28-day rule that defines
``epoch`` is purely temporal and has no notion of where the comet was: 2024 E1
epoch 1 spans r_h = 4.00-2.74 au inside one "spectrum".  Each epoch is cut
under four rules in priority order (see :class:`config.GroupingConfig`).

**A greedy sweep cannot build this.**  The rules are not hereditary and need
lookahead, so :func:`split_groups` is a dynamic programme over
``(band breaks, n_groups, cut penalty, spread)`` compared lexicographically.

Provenance
----------
The output never modifies the photometry tables.  Every exposure's label is
written to ``data/comspec/phase_assignment.csv``; the per-group table goes to
``results/comspec/phase_map.csv``.
"""

from __future__ import annotations

import time
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import EMISSION_WINDOWS, GroupingConfig
from .dataio import exposure_table, list_targets, slug
from .logging_utils import get_logger

__all__ = ["orbital_arc", "band_runs", "cut_costs", "split_groups", "manual_segments",
           "subdivide", "group_summary", "regroup_all"]

log = get_logger("grouping")


def orbital_arc(r_hel, cfg: GroupingConfig) -> Tuple[np.ndarray, "int | None"]:
    """
    Label each epoch inbound (0) or outbound (1) about the minimum of ``r_hel``.

    A perihelion passage is accepted only when ``r_hel`` both falls and rises
    by more than ``cfg.arc_min_drh`` around its minimum, with at least
    ``cfg.arc_min_epoch`` epochs on each side.  Without that guard the last
    point of a monotone series -- a fraction of a micro-au above its neighbour
    through ephemeris round-off -- would be labelled a separate arc.
    """
    r = np.asarray(r_hel, float)
    k = int(np.argmin(r))
    resolved = (r[0] - r[k] > cfg.arc_min_drh and r[-1] - r[k] > cfg.arc_min_drh
                and k >= cfg.arc_min_epoch and len(r) - k >= cfg.arc_min_epoch)
    if not resolved:
        return np.zeros(len(r), int), None
    arc = np.zeros(len(r), int)
    arc[k:] = 1
    return arc, k


def band_runs(r_hel, sampled, tol: float) -> List[Tuple[int, int]]:
    """
    Maximal runs of band epochs that can be held in one group, as index ranges.

    A run grows while the *whole interval it spans* still satisfies the
    tolerance, so a run is always feasible by itself.  Single-epoch runs are
    dropped: a lone point cannot be split and constrains nothing.
    """
    hits = np.flatnonzero(sampled)
    runs, start = [], None
    for pos, i in enumerate(hits):
        if start is None:
            start = i
            continue
        seg = r_hel[start:i + 1]
        if (seg.max() - seg.min()) / seg.mean() < tol:
            continue
        runs.append((start, hits[pos - 1]))
        start = i
    if start is not None:
        runs.append((start, hits[-1]))
    return [(a, b) for a, b in runs if b > a]


def cut_costs(r_hel, sampled_by_band: Dict[str, np.ndarray], cfg: GroupingConfig):
    """
    Per-boundary cost of cutting: (band runs broken, link penalty).

    Boundary ``b`` sits between epochs ``b-1`` and ``b``.  Counting incidences
    rather than distinct broken runs keeps the cost additive over cuts and
    penalises shattering one run across several groups more than splitting it
    once.
    """
    n = len(r_hel)
    n_break = np.zeros(n, int)
    for sampled in sampled_by_band.values():
        for s, e in band_runs(r_hel, sampled, cfg.rh_tol):
            n_break[s + 1:e + 1] += 1
    link = np.zeros(n)
    d = np.abs(np.diff(r_hel))
    link[1:] = np.maximum(0.0, 1.0 - d / cfg.link_drh)
    return n_break, link


def split_groups(r_hel, weight, n_break, link, tol: float) -> List[Tuple[int, int]]:
    """
    Cut a time-ordered epoch block into groups, minimising the four-part cost.

    ``dp[i]`` is the best partition of the first ``i`` epochs as the tuple
    ``(band breaks, n_groups, link penalty, total spread)``, compared
    lexicographically -- band integrity outranks group count, group count
    outranks where the cuts land, and the spread only breaks ties.  The inner
    loop stops on the monotone bound ``1 - r_min/r_max >= tol``; feasibility is
    then tested with the true (measurement-weighted) criterion.

    Returns
    -------
    list of (start, stop)
        Half-open index ranges in time order, covering the input exactly.
    """
    n = len(r_hel)
    unreachable = (n + 1, n + 1, np.inf, np.inf)
    dp = [unreachable] * (n + 1)
    dp[0] = (0, 0, 0.0, 0.0)
    back = [-1] * (n + 1)
    for i in range(1, n + 1):
        r_max, r_min, w_sum, rw_sum = -np.inf, np.inf, 0.0, 0.0
        for j in range(i - 1, -1, -1):
            r = r_hel[j]
            r_max = max(r_max, r)
            r_min = min(r_min, r)
            w_sum += weight[j]
            rw_sum += r * weight[j]
            if (r_max - r_min) >= tol * r_max:
                break
            spread = (r_max - r_min) / (rw_sum / w_sum)
            if spread >= tol:
                continue
            b, l = (n_break[j], link[j]) if j > 0 else (0, 0.0)
            cand = (dp[j][0] + b, dp[j][1] + 1, dp[j][2] + l, dp[j][3] + spread)
            if cand < dp[i]:
                dp[i], back[i] = cand, j
    out, i = [], n
    while i > 0:
        j = back[i]
        out.append((j, i))
        i = j
    return out[::-1]


def manual_segments(r_hel, edges: Sequence[float]) -> List[Tuple[int, int]]:
    """
    Hand-specified descending r_h bins as index ranges in time order.

    ``edges = [a, b, c]`` gives ``(a, inf), (b, a], (c, b], (-inf, c]``.  Inside
    one arc r_h is monotone in time, so each bin is a contiguous stretch; that
    is asserted rather than assumed.
    """
    lab = np.zeros(len(r_hel), int)
    for i, cut in enumerate(sorted(edges, reverse=True)):
        lab[np.asarray(r_hel) <= cut] = i + 1
    bounds = [0, *(np.flatnonzero(lab[1:] != lab[:-1]) + 1).tolist(), len(lab)]
    segs = list(zip(bounds[:-1], bounds[1:]))
    assert len({lab[a] for a, _ in segs}) == len(segs), "an r_h bin is not contiguous in time"
    return segs


def subdivide(target: str, ep: pd.DataFrame, cfg: GroupingConfig):
    """
    Assign a ``phase`` to every epoch of one target.

    Blocks are ``(epoch, arc)`` pairs -- both mandatory cuts -- and each block
    is partitioned by :func:`split_groups` or, for a target in
    ``cfg.manual_edges``, by :func:`manual_segments`.  Labels run 1..N over the
    whole target in time order.
    """
    key = slug(target)
    ep = ep.copy()
    arc, k = orbital_arc(ep.r_hel.to_numpy(), cfg)
    ep["arc_seen"] = arc
    ep["arc"] = 0 if key in cfg.manual_no_arc else arc
    ep.attrs["perihelion_seen"] = k
    ep.attrs["perihelion_index"] = None if key in cfg.manual_no_arc else k

    label = np.zeros(len(ep), int)
    is_manual = np.zeros(len(ep), bool)
    cuts = []
    n = 0
    for _, g in ep.groupby(["epoch", "arc"], sort=True):
        idx = g.index.to_numpy()
        r_hel = g.r_hel.to_numpy(float)
        n_break, link = cut_costs(r_hel, {b: g[b].to_numpy() for b in EMISSION_WINDOWS}, cfg)
        edges = cfg.manual_edges.get(key, {}).get("out" if g.arc.iloc[0] else "in")
        manual = edges is not None
        segs = (manual_segments(r_hel, edges) if manual
                else split_groups(r_hel, g.weight.to_numpy(float), n_break, link, cfg.rh_tol))
        for a, b in segs:
            n += 1
            label[idx[a:b]] = n
            is_manual[idx[a:b]] = manual
            if a > 0:
                cuts.append(dict(target=key, epoch=int(g.epoch.iloc[0]), arc=int(g.arc.iloc[0]),
                                 new_group=n, mode="manual" if manual else "auto",
                                 drh_at_cut=abs(r_hel[a] - r_hel[a - 1]),
                                 bands_broken=int(n_break[a]), link_penalty=float(link[a])))
    ep["phase"] = label
    ep["manual"] = is_manual
    return ep, cuts


def group_summary(target: str, ep: pd.DataFrame) -> List[dict]:
    """Per-group r_h, timing and band-sampling statistics for the epoch and phase groupings."""
    rows = []
    for keys, name in ((["epoch"], "epoch"), (["epoch", "arc", "phase"], "phase")):
        for k, g in ep.groupby(keys, sort=True):
            k = k if isinstance(k, tuple) else (k,)
            w = g.weight.to_numpy(float)
            mean = float(np.average(g.r_hel, weights=w))
            row = dict(target=slug(target), grouping=name, epoch=int(k[0]), group=int(k[-1]),
                       arc=("out" if k[1] else "in") if len(k) > 2 else "",
                       n_epoch=len(g), n_meas=int(w.sum()),
                       r_hel_mean=mean, r_hel_min=float(g.r_hel.min()),
                       r_hel_max=float(g.r_hel.max()),
                       spread=(float(g.r_hel.max()) - float(g.r_hel.min())) / mean,
                       r_obs_mean=float(np.average(g.r_obs, weights=w)),
                       r_obs_min=float(g.r_obs.min()), r_obs_max=float(g.r_obs.max()),
                       jd_start=float(g.jd_utc.min()), jd_end=float(g.jd_utc.max()),
                       tspan_day=float(g.jd_utc.max() - g.jd_utc.min()),
                       manual=bool(g.manual.any()) if "manual" in g else False)
            row.update({f"n_{b}": int(g[b].sum()) for b in EMISSION_WINDOWS})
            rows.append(row)
    return rows


def regroup_all(targets: Sequence[str] | None = None, cfg: GroupingConfig | None = None):
    """
    Run the regrouping over the catalog.

    Returns
    -------
    assignment : DataFrame
        One row per exposure: ``target, filename, jd_utc, epoch, arc, phase, manual``.
    group_map : DataFrame
        One row per ``phase`` group with its r_h / r_obs / timing / band-sampling statistics,
        ``parent_spread``, ``was_split`` and ``manual``.
    cuts : DataFrame
        Every cut made and what it cost.
    epochs : dict of DataFrame
        The per-target epoch tables, for figures.
    """
    cfg = cfg or GroupingConfig()
    targets = list(targets) if targets else list_targets()
    t0 = time.time()
    epochs, cuts, rows, assign = {}, [], [], []
    for t in targets:
        ep, c = subdivide(t, exposure_table(t), cfg)
        epochs[t] = ep
        cuts.extend(c)
        rows.extend(group_summary(t, ep))
        df = ep[["obsid", "jd_utc", "epoch", "arc", "phase", "manual"]].copy()
        df.insert(0, "target", slug(t))
        assign.append(df)

    groups = pd.DataFrame(rows)
    old = groups[groups.grouping == "epoch"].set_index(["target", "epoch"])
    new = groups[groups.grouping == "phase"].drop(columns="grouping")
    new = new.rename(columns={"group": "phase"})
    n_sub = new.groupby(["target", "epoch"]).size()
    new["parent_spread"] = new.set_index(["target", "epoch"]).index.map(old.spread)
    new["was_split"] = new.set_index(["target", "epoch"]).index.isin(n_sub[n_sub > 1].index)
    new = new.sort_values(["target", "phase"]).reset_index(drop=True)

    # Because the per-exposure table carries one row per exposure *time* and a spectrum is
    # measured once per exposure, the assignment must reach every filename of the target.
    assignment = pd.concat(assign, ignore_index=True)
    full = []
    for t in targets:
        from .dataio import load_apphot
        fn = load_apphot(t).drop_duplicates("filename")[["filename", "obsid", "jd_utc_exact"]]
        lab = assignment[assignment.target == slug(t)].drop(columns=["jd_utc"]).set_index("obsid")
        m = fn.merge(lab, left_on="obsid", right_index=True, how="left") \
              .rename(columns={"jd_utc_exact": "jd_utc"})
        assert m.phase.notna().all(), f"{t}: exposures without a phase label"
        full.append(m.assign(target=slug(t)))
    assignment = pd.concat(full, ignore_index=True)[
        ["target", "filename", "obsid", "jd_utc", "epoch", "arc", "phase", "manual"]]
    assignment["phase"] = assignment.phase.astype(int)

    auto = new[~new.manual]
    assert (auto.spread < cfg.rh_tol).all(), "an automatic group violates rule 1"
    log.info("regrouped %d targets: %d epochs -> %d phases (%d epochs subdivided, %d manual "
             "groups) in %.1f s", len(targets), len(old), len(new), int((n_sub > 1).sum()),
             int(new.manual.sum()), time.time() - t0)
    return assignment, new, pd.DataFrame(cuts), epochs
