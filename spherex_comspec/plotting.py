"""
Every figure of the pipeline.  Imported by ``main.py --figures`` and the notebooks only.

Style follows ``notebooks/rcparams.py`` (base font 20 pt, imported for its side
effects); dense multi-panel figures are enlarged, never shrunk.  Batch output
is written at 50 dpi per the project convention, one-off figures at 200.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import warnings

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from . import directory as _dir
from .config import (BAND_CARRIER, BAND_COLORS, BAND_WINDOWS, EMISSION_WINDOWS, KEY_RANGES,
                     SPECIES, Variant)
from .logging_utils import get_logger

__all__ = ["apply_rcparams", "savefig", "plot_phase_comparison", "plot_phase_summary",
           "plot_raw_spectrum", "plot_validation_grid", "plot_fit", "plot_summary_Q",
           "plot_mixing_ratios", "plot_flag_comparison", "plot_distcorr_effect",
           "save_variant_figures", "save_grouping_figures"]

log = get_logger("plotting")
DPI_BATCH = 50
DPI_ONE = 200
VERDICT_COLORS = {"PASS": "tab:green", "WARN": "tab:orange", "FAIL": "tab:red"}
SP_COLORS = {"H2O": "tab:blue", "CO2": "tab:orange", "CO": "tab:green"}
PRETTY = {"H2O": r"H$_2$O", "CO2": r"CO$_2$", "CO": "CO"}
_PALETTE = list(plt.get_cmap("tab10").colors) + list(plt.get_cmap("Set2").colors)


def apply_rcparams() -> bool:
    nb = str(_dir.NOTEBOOK_DIR)
    if nb not in sys.path:
        sys.path.insert(0, nb)
    try:
        import rcparams  # noqa: F401
        plt.rcParams.update({"figure.dpi": 100})
        return True
    except ImportError:
        log.warning("notebooks/rcparams.py not importable; matplotlib defaults in use")
        return False


def savefig(fig, path: Path, dpi: int = DPI_BATCH) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


# ------------------------------------------------------------------------ phase grouping
def _colors(values):
    return {v: _PALETTE[i % len(_PALETTE)] for i, v in enumerate(sorted(pd.unique(values)))}


def _robust_ylim(y, frac=0.02, pad=0.12):
    lo, hi = np.percentile(y, [100 * frac, 100 * (1 - frac)])
    if hi <= lo:
        lo, hi = float(np.min(y)), float(np.max(y)) + 1e-9
    m = pad * (hi - lo)
    return min(lo - m, 0.0), hi + m


def _fmt_r_ap(r_ap_km) -> str:
    """``20000 -> "20,000 km"``; a set of per-phase apertures -> ``"20,000 / 60,000 km"``."""
    r = np.atleast_1d(np.asarray(r_ap_km, float))
    return " / ".join(f"{x:,.0f}" for x in np.unique(r)) + " km"


def plot_phase_comparison(target: str, ep: pd.DataFrame, spec: pd.DataFrame,
                          r_ap_km, flux_col: str = "flux_raw"):
    """
    Before/after grouping figure: (r_h vs time, raw spectrum) x (epoch, phase).

    ``ep`` is the per-exposure table from the grouping (with ``epoch``, ``phase``
    and the perihelion index in ``.attrs``); ``spec`` the selected spectrum with
    a ``phase`` column.
    """
    k_peri = ep.attrs.get("perihelion_index")
    t0 = ep.jd_utc.min()
    y = spec[flux_col].to_numpy(float)
    ylo, yhi = _robust_ylim(y)
    n_clip = int(((y < ylo) | (y > yhi)).sum())
    fig, axes = plt.subplots(2, 2, figsize=(24, 15), gridspec_kw=dict(width_ratios=[1, 1.4]))
    for row, key in enumerate(["epoch", "phase"]):
        colors = _colors(ep[key])
        ax_t, ax_s = axes[row]
        for v, g in ep.groupby(key):
            sp = (g.r_hel.max() - g.r_hel.min()) / g.r_hel.mean()
            ax_t.scatter(g.jd_utc - t0, g.r_hel, s=80, color=colors[v], edgecolors="0.2", lw=0.7,
                         zorder=3, label=f"{v}:  {int(g.weight.sum())} meas,  {100 * sp:.1f} %")
        if key == "phase":
            o = ep.sort_values("jd_utc")
            lab, jd = o.phase.to_numpy(), o.jd_utc.to_numpy()
            for i in np.flatnonzero(lab[1:] != lab[:-1]):
                ax_t.axvline(0.5 * (jd[i] + jd[i + 1]) - t0, color="0.35", ls="--", lw=2, zorder=1)
        if k_peri is not None:
            ax_t.axvline(ep.jd_utc.iloc[k_peri] - t0, color="tab:red", lw=2.5, alpha=0.7, zorder=2)
            ax_t.text(ep.jd_utc.iloc[k_peri] - t0, ax_t.get_ylim()[1], " perihelion",
                      color="tab:red", fontsize=15, va="top", ha="left", fontweight="bold")
        ax_t.set(xlabel="days since first observation", ylabel=r"$r_h$ [au]")
        ax_t.set_title(f"{key}:  {ep[key].nunique()} group{'s' if ep[key].nunique() > 1 else ''}",
                       pad=12)
        ax_t.legend(fontsize=13, ncol=2 if ep[key].nunique() > 6 else 1, framealpha=0.9,
                    title=r"group:  measurements,  $r_h$ spread", title_fontsize=13)
        for b, (lo, hi) in EMISSION_WINDOWS.items():
            ax_s.axvspan(lo, hi, color=BAND_COLORS[b], alpha=0.16, zorder=0)
            if row == 0:
                ax_s.text(0.5 * (lo + hi), yhi, f"\n{b}", ha="center", va="top", fontsize=14,
                          color=BAND_COLORS[b], fontweight="bold")
        for v, g in spec.groupby(key):
            ax_s.errorbar(g.wl, g[flux_col], yerr=g.err_raw if flux_col == "flux_raw" else g.err,
                          fmt="none", ecolor="0.8", lw=1.0, zorder=1)
            ax_s.scatter(g.wl, g[flux_col], s=50, color=colors.get(v, "k"), edgecolors="0.25",
                         lw=0.5, zorder=3)
        ax_s.axhline(0, color="k", lw=1.0, ls=":")
        ax_s.set(xlabel=r"wavelength [$\mu$m]", ylabel="source sum [mJy]", xlim=(0.7, 5.05),
                 ylim=(ylo, yhi))
        ax_s.set_title(f"raw spectrum,  $r_{{ap}}$ = {_fmt_r_ap(r_ap_km)},  N = {len(spec)}"
                       + (f"  ({n_clip} off scale)" if n_clip else ""), pad=12)
    n_old, n_new = ep.epoch.nunique(), ep.phase.nunique()
    arc = "  split at perihelion" if k_peri is not None else ""
    fig.suptitle(f"{target}    phase grouping  {n_old} $\\rightarrow$ {n_new}"
                 f"    $r_h$ = {ep.r_hel.min():.3f}$-${ep.r_hel.max():.3f} au{arc}", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    return fig


def plot_phase_summary(group_map: pd.DataFrame, epoch_map: pd.DataFrame, rh_tol: float = 0.10,
                       min_points: int = 5):
    """Catalog-level view: groups per epoch, r_h spread before/after, group size vs gate."""
    NEW, OLD = group_map, epoch_map
    fig, axes = plt.subplots(1, 3, figsize=(24, 7.5))
    ax = axes[0]
    vc = NEW.groupby(["target", "epoch"]).size().value_counts().sort_index()
    bars = ax.bar(vc.index.astype(str), vc.values, 0.65, color="tab:purple", alpha=0.85,
                  edgecolor="0.3", lw=1.5)
    for b, v in zip(bars, vc.values):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v}", ha="center", va="bottom", fontsize=16)
    ax.set(xlabel="groups the epoch became", ylabel="parent epochs")
    ax.set_title(f"{len(OLD)} epochs $\\rightarrow$ {len(NEW)} phases", pad=14)
    ax = axes[1]
    bins = np.linspace(0, 36, 37)
    ax.hist(100 * OLD.spread, bins=bins, color="0.75", edgecolor="0.4", lw=1.2, label="epoch")
    ax.hist(100 * NEW[~NEW.manual].spread, bins=bins, color="tab:blue", alpha=0.8,
            label="phase (automatic)")
    ax.hist(100 * NEW[NEW.manual].spread, bins=bins, color="tab:red", alpha=0.9,
            label="phase (manual)")
    ax.axvline(100 * rh_tol, color="k", lw=3, ls="--")
    ax.set(xlabel=r"$(r_h^{max}-r_h^{min})\,/\,\langle r_h\rangle$  [%]", ylabel="groups",
           yscale="log")
    ax.set_title("$r_h$ spread within a group", pad=14)
    ax.legend(fontsize=14)
    ax = axes[2]
    bins = np.logspace(0, np.log10(max(NEW.n_meas.max(), OLD.n_meas.max())) + 0.05, 34)
    ax.hist(OLD.n_meas, bins=bins, color="0.75", edgecolor="0.4", lw=1.2, label="epoch")
    ax.hist(NEW.n_meas, bins=bins, color="tab:blue", alpha=0.8, label="phase")
    ax.axvline(min_points, color="tab:red", lw=3, ls="--")
    ax.set(xlabel="measurements in the group", ylabel="groups", xscale="log")
    ax.set_title(f"group size   ({int((NEW.n_meas < min_points).sum())} below the gate)", pad=14)
    ax.legend(fontsize=14)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------- continuum figures
def plot_raw_spectrum(raw: pd.DataFrame, target: str, r_ap_km: float, phase, flux_space: str,
                      ax=None):
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 7))
    else:
        fig = ax.get_figure()
    ax.errorbar(raw.wl, raw.flux, yerr=raw.err, fmt="o", ms=5, lw=1.2, color="0.35",
                ecolor="0.75", label=f"{flux_space} flux")
    for name, b in BAND_WINDOWS.items():
        ax.axvspan(*b["em"], color=BAND_COLORS[name], alpha=0.20)
        ax.axvspan(b["cont"][0], b["cont"][1], color=BAND_COLORS[name], alpha=0.06)
        ax.text(b["lam_c"], ax.get_ylim()[1] * 0.92, name, ha="center", fontsize=18,
                color=BAND_COLORS[name])
    ax.axhline(0, color="k", lw=0.6, ls=":")
    unit = r"$F\,r_h^2\Delta^2$ [mJy at 1 au]" if flux_space == "distcorr" else "source sum [mJy]"
    ax.set(xlabel=r"wavelength [$\mu$m]", ylabel=unit, xlim=(0.7, 5.05),
           title=f"{target}  |  $r_{{ap}}$ = {r_ap_km:,.0f} km  |  phase {phase}  |  N = {len(raw)}")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return fig, ax


def _band_panel(raw, fit, val, sub, cfg, axes, flux_space):
    b, col = fit["band"], BAND_COLORS[fit["band"]]
    ax0, ax1, ax2 = axes
    lo, hi = fit["cont"][0] - cfg.plot_margin_um, fit["cont"][1] + cfg.plot_margin_um
    win = raw[(raw.wl >= lo) & (raw.wl <= hi)]
    lam_g = np.linspace(lo, hi, 400)
    from .continuum import continuum_at
    ax0.errorbar(win.wl, win.flux, yerr=win.err, fmt="o", ms=6, lw=1.2, color="0.55", ecolor="0.8",
                 zorder=1, label="all points")
    if fit["ok"]:
        keep = fit["keep"]
        ax0.errorbar(fit["xc"][keep], fit["yc"][keep], yerr=fit["ec"][keep], fmt="o", ms=9, lw=1.6,
                     color="k", ecolor="0.4", zorder=3, label="continuum sample")
        if (~keep).any():
            lbl = (f"not used by 2-point fit ({(~keep).sum()})" if fit["method"] == "2point"
                   else f"sigma-clipped ({(~keep).sum()})")
            ax0.plot(fit["xc"][~keep], fit["yc"][~keep], "x", ms=15, mew=3, color="tab:red",
                     zorder=4, label=lbl)
        cg, sg = continuum_at(fit, lam_g)
        ax0.plot(lam_g, cg, color=col, lw=2, zorder=5,
                 label=f"continuum ({fit['method']}, order {fit['order_used']:.0f})")
        ax0.fill_between(lam_g, cg - sg, cg + sg, color=col, alpha=0.25, lw=0, zorder=2)
    # source-flagged points are ringed so contamination is visible where it lands
    if "sourceflag" in win.columns:
        for fl, mk in (("a", "^"), ("b", "v")):
            w = win[win.sourceflag.astype(str) == fl]
            if len(w):
                ax0.plot(w.wl, w.flux, mk, ms=13, mfc="none", mec="crimson", mew=2, zorder=6,
                         label=f"flag {fl} ({len(w)})")
    ax0.axvspan(*fit["em"], color=col, alpha=0.13, zorder=0)
    if len(win):
        span = np.nanpercentile(win.flux, [2, 98])
        pad = 0.25 * max(np.ptp(span), 1e-3)
        ax0.set_ylim(span[0] - pad, span[1] + pad)
    ax0.set(xlim=(lo, hi), xlabel=r"$\lambda$ [$\mu$m]",
            ylabel="flux [mJy]" if flux_space == "physical" else r"$F r_h^2\Delta^2$ [mJy]",
            title=f"{b}  {BAND_CARRIER[b]}")
    ax0.legend(loc="best", framealpha=0.9, fontsize=13)

    sw = sub[(sub.wl >= lo) & (sub.wl <= hi)]
    ax1.axhline(0, color="k", lw=1, ls="--")
    ax1.axvspan(*fit["em"], color=col, alpha=0.13)
    out, ins = sw[~sw.in_emission], sw[sw.in_emission]
    ax1.errorbar(out.wl, out.emis_mjy, yerr=out.emis_err_mjy, fmt="o", ms=6, lw=1.2, color="0.55",
                 ecolor="0.8", label="continuum region")
    ax1.errorbar(ins.wl, ins.emis_mjy, yerr=ins.emis_err_mjy, fmt="o", ms=10, lw=1.6, color=col,
                 ecolor=col, label=f"emission window (N = {len(ins)})")
    if len(ins) and np.isfinite(ins.emis_mjy).any():
        k = ins.emis_mjy.idxmax()
        ax1.plot(ins.wl[k], ins.emis_mjy[k], "*", ms=30, color="crimson", zorder=6,
                 label=f"peak {ins.emis_mjy[k]:.2f} mJy @ {ins.wl[k]:.3f} µm")
    ax1.set(xlim=(lo, hi), xlabel=r"$\lambda$ [$\mu$m]", ylabel="flux $-$ continuum [mJy]",
            title="continuum-subtracted")
    ax1.legend(loc="best", framealpha=0.9, fontsize=13)

    if fit["ok"]:
        keep = fit["keep"]
        chi = (fit["yc"] - np.polyval(fit["coeffs"], fit["xc"] - fit["lam_c"])) / fit["ec"]
        ax2.axhline(0, color="k", lw=1, ls="--")
        for s in (-cfg.sigma, cfg.sigma):
            ax2.axhline(s, color="tab:red", lw=0.9, ls=":")
        ax2.plot(fit["xc"][keep], chi[keep], "o", ms=8, color="k")
        if (~keep).any():
            ax2.plot(fit["xc"][~keep], chi[~keep], "x", ms=15, mew=3, color="tab:red")
        ax2.set_ylim(np.array([-1, 1]) * max(4.0, 1.15 * np.nanmax(np.abs(chi[keep]))))
    c2 = "n/a" if not np.isfinite(val["chi2_red"]) else f"{val['chi2_red']:.1f}"
    es = "n/a" if not np.isfinite(val["err_scale"]) else f"{val['err_scale']:.1f}"
    ax2.set(xlabel=r"$\lambda$ [$\mu$m]", ylabel=r"$(y-\hat{y})/\sigma$")
    ax2.set_title(r"$\chi^2_\nu$ = %s    $\sigma$-scale = %s" % (c2, es), pad=14)
    ax0.text(0.0, 1.14, val["verdict"], transform=ax0.transAxes, va="bottom", ha="left",
             fontsize=20, fontweight="bold", color=VERDICT_COLORS[val["verdict"]],
             bbox=dict(fc="w", ec=VERDICT_COLORS[val["verdict"]], lw=1.5, alpha=0.95))
    if val["notes"]:
        ax1.text(0.5, -0.30, textwrap.fill(val["notes"], 95), transform=ax1.transAxes,
                 ha="center", va="top", fontsize=15, color="tab:red", style="italic")


def plot_validation_grid(out: dict, target: str, r_ap_km: float, phase, cfg,
                         flux_space: str = "physical"):
    """The three-band continuum validation grid for one group (``out`` from ``process_group``)."""
    fits, vals, subs, raw = out["fits"], out["valid"], out["subs"], out["raw"]
    fig, axes = plt.subplots(len(fits), 3, figsize=(24, 7.4 * len(fits)),
                             gridspec_kw=dict(width_ratios=[2.1, 2.1, 1.0], wspace=0.28, hspace=0.62))
    axes = np.atleast_2d(axes)
    for k, name in enumerate(fits):
        _band_panel(raw, fits[name], vals[name], subs[name], cfg, axes[k], flux_space)
    fig.suptitle(f"{target}   $r_{{ap}}$ = {r_ap_km:,.0f} km   phase {phase}   |   "
                 f"continuum in {flux_space} space", y=1.0, fontsize=26)
    return fig


# ------------------------------------------------------------------------ fit figures
def plot_fit(fit, points: pd.DataFrame, curves: dict, fit_cfg):
    """Intrinsic and instrument-convolved models over the continuum-subtracted data."""
    from .fitting import channel_masks
    windows = EMISSION_WINDOWS
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), gridspec_kw=dict(wspace=0.24))
    lam, intr, conv = curves["lam"], curves["intrinsic"], curves["convolved"]
    DASHES = {"H2O": (0, (6, 3)), "CO2": (0, (2, 2)), "CO": (0, (1, 1.6))}
    for ax, (bname, (lo, hi)) in zip(axes, windows.items()):
        pad = 0.06 * (hi - lo)
        m = (lam >= lo - 4 * pad) & (lam <= hi + 4 * pad)
        d = points[points.band == bname]
        ax.axvspan(lo, hi, color=BAND_COLORS[bname], alpha=0.12, zorder=0)
        ax.axhline(0, color="k", lw=1.2, ls=":", zorder=1)
        ax.plot(lam[m], intr[m], color="0.35", lw=2.0, zorder=3)
        for s, v in curves["per_species"].items():
            if np.nanmax(np.abs(v[m])) > 0:
                ax.plot(lam[m], v[m], ls=DASHES[s], lw=2.6, color=SP_COLORS[s], alpha=0.9, zorder=4)
        ax.plot(lam[m], conv[m], color="crimson", lw=3.5, zorder=5)
        n_used = n_excl = 0
        if len(d):
            mk = channel_masks(d, fit_cfg, "physical")
            used, excl = d[mk["used"]], d[mk["negative"]]
            n_used, n_excl = len(used), len(excl)
            ax.errorbar(used.wl, used.emis_raw_mjy, yerr=used.emis_raw_err_mjy, fmt="o", ms=10,
                        lw=2, color="k", ecolor="0.4", zorder=6)
            if "sourceflag" in used:
                for fl, mkr in (("a", "^"), ("b", "v")):
                    w = used[used.sourceflag.astype(str) == fl]
                    if len(w):
                        ax.plot(w.wl, w.emis_raw_mjy, mkr, ms=15, mfc="none", mec="crimson",
                                mew=2, zorder=7)
            if len(excl):
                ax.errorbar(excl.wl, excl.emis_raw_mjy, yerr=excl.emis_raw_err_mjy, fmt="none",
                            ecolor="tab:red", lw=1.6, alpha=0.65, zorder=6)
                ax.plot(excl.wl, excl.emis_raw_mjy, "x", ms=15, mew=3.5, color="tab:red", zorder=7)
        ax.set(xlim=(lo - 4 * pad, hi + 4 * pad), xlabel=r"$\lambda$ [$\mu$m]",
               ylabel="continuum-subtracted flux [mJy]" if bname == "2.7um" else "")
        in_fit = bname in fit.bands_used
        ax.set_title(f"{bname}   {BAND_CARRIER[bname]}", pad=38)
        sub = (f"{n_used} channels fitted" + (f", {n_excl} excluded" if n_excl else "")) if in_fit \
            else "not fitted"
        ax.text(0.5, 1.02, sub, transform=ax.transAxes, ha="center", va="bottom", fontsize=14,
                color="k" if in_fit else "tab:red")
    qtxt = []
    for k, s in enumerate(fit.species):
        st = fit.status[k]
        if st == "not_covered":
            qtxt.append(f"Q({PRETTY[s]}) not covered")
        elif st in ("upper_limit", "negative_fit"):
            qtxt.append(f"Q({PRETTY[s]}) $<$ {fit.Q_limit[k]:.2e}")
        elif st == "marginal":
            qtxt.append(f"Q({PRETTY[s]}) = {fit.Q[k]:.2e} $\\pm$ {fit.Q_err[k]:.1e} (marginal)")
        elif st == "rejected":
            qtxt.append(f"Q({PRETTY[s]}) rejected (case revision)")
        else:
            qtxt.append(f"Q({PRETTY[s]}) = {fit.Q[k]:.2e} $\\pm$ {fit.Q_err[k]:.1e}")
    fig.suptitle(f"{fit.target}   phase {fit.phase}   $r_{{ap}}$ = {fit.r_ap_km:,.0f} km   "
                 f"$\\langle r_h \\rangle$ = {fit.geometry['r_hel_mean']:.3f} au", y=1.15)
    anchor = ("" if fit.h2o_source == "main" else
              "   |   H$_2$O from 4.6–4.9 µm hot bands (2.7 µm not covered)" if fit.h2o_source == "hot"
              else "   |   H$_2$O not covered")
    dropped = f"   |   {fit.n_dropped_negative} channel(s) excluded" if fit.n_dropped_negative else ""
    fig.text(0.5, 1.045, "      ".join(qtxt) + r"   [s$^{-1}$]", ha="center", fontsize=16)
    fig.text(0.5, 1.005, f"$\\chi^2_\\nu$ = {fit.chi2_red:.1f}{anchor}{dropped}", ha="center",
             fontsize=15, color="k" if fit.h2o_source == "main" else "tab:orange" if fit.h2o_source == "hot" else "tab:red")
    handles = [
        plt.Line2D([], [], color="k", marker="o", ls="none", ms=10, label="fitted channels"),
        plt.Line2D([], [], color="tab:red", marker="x", ls="none", ms=13, mew=3.5,
                   label=r"excluded ($>1\sigma$ below zero)"),
        plt.Line2D([], [], color="crimson", marker="^", ls="none", ms=12, mfc="none", mew=2,
                   label="flag a / b (ring)"),
        plt.Line2D([], [], color="crimson", lw=3.5, label="model convolved with SPHEREx channel"),
        plt.Line2D([], [], color="0.35", lw=2.0, label="model intrinsic"),
        plt.Line2D([], [], color="tab:blue", lw=2.0, ls="--", label=r"H$_2$O"),
        plt.Line2D([], [], color="tab:orange", lw=2.0, ls="--", label=r"CO$_2$"),
        plt.Line2D([], [], color="tab:green", lw=2.0, ls="--", label="CO"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=14,
               bbox_to_anchor=(0.5, -0.14))
    return fig


# --------------------------------------------------------------------- summary figures
def _aperture_class(r_ap_km) -> np.ndarray:
    """Three marker classes: the per-aperture legend of the S/N rule (19 radii) is unreadable."""
    r = np.asarray(r_ap_km, float)
    return np.where(r <= 20000, 0, np.where(r <= 40000, 1, 2))


_AP_MARKERS = ("o", "s", "D")


def _aperture_class_labels(r_ap_km) -> list:
    """One label per class from the radii actually present: ``"20 000 km"`` when the class holds a
    single radius (the fixed rule), ``"22 000–40 000 km"`` when it holds several (the S/N rule)."""
    r = np.asarray(r_ap_km, float)
    cls = _aperture_class(r)
    labels = []
    for k in range(3):
        rk = np.unique(r[cls == k])
        if len(rk) == 0:
            labels.append(("≤ 20 000 km", "22 000–40 000 km", "≥ 60 000 km")[k])
        elif len(rk) == 1:
            labels.append(f"{rk[0]:,.0f} km")
        else:
            labels.append(f"{rk.min():,.0f}–{rk.max():,.0f} km")
    return labels


def plot_summary_Q(fits: pd.DataFrame, title: str = ""):
    """Q vs mean r_h per species: ≥ 3σ detections filled, marginal (1–3σ) hollow, 3σ limits grey;
    the marker shape is the aperture class."""
    fig, axes = plt.subplots(1, 3, figsize=(24, 9))
    _AP_CLASS = list(zip(_AP_MARKERS, _aperture_class_labels(fits.r_ap_km)))
    for ax, s in zip(axes, SPECIES):
        det = fits[(fits[f"Q_{s}_status"] == "detected") & (fits[f"Q_{s}_n_eff"] >= 2)]
        mar = fits[fits[f"Q_{s}_status"] == "marginal"]
        lim = fits[fits[f"Q_{s}_status"].isin(["upper_limit", "negative_fit"])]
        hot = (det.h2o_source == "hot") if (s == "H2O" and "h2o_source" in det) else pd.Series(False, index=det.index)
        for k, (mk, lab) in enumerate(_AP_CLASS):
            dd = det[(_aperture_class(det.r_ap_km) == k) & ~hot.to_numpy()]
            if len(dd):
                ax.errorbar(dd.r_hel_mean, dd[f"Q_{s}"], yerr=dd[f"Q_{s}_err"], fmt=mk, ms=11, lw=1.5,
                            capsize=4, color=SP_COLORS[s], label=f"≥ 3σ, n_eff ≥ 2 ({lab})")
            hh = det[(_aperture_class(det.r_ap_km) == k) & hot.to_numpy()]
            if len(hh):
                ax.errorbar(hh.r_hel_mean, hh[f"Q_{s}"], yerr=hh[f"Q_{s}_err"], fmt=mk, ms=11, lw=1.5,
                            capsize=4, color=SP_COLORS[s], mfc="none", mew=2.5,
                            label=f"≥ 3σ from hot bands ({lab})")
            mm = mar[_aperture_class(mar.r_ap_km) == k]
            if len(mm):
                ax.errorbar(mm.r_hel_mean, mm[f"Q_{s}"], yerr=mm[f"Q_{s}_err"], fmt=mk, ms=8, lw=1.0,
                            capsize=3, color=SP_COLORS[s], mfc="none", mew=1.2, alpha=0.55,
                            label=f"marginal 1–3σ ({lab})")
            uu = lim[_aperture_class(lim.r_ap_km) == k]
            if len(uu):
                ax.errorbar(uu.r_hel_mean, uu[f"Q_{s}_upper_limit"], yerr=0.35 * uu[f"Q_{s}_upper_limit"],
                            fmt=mk, ms=7, lw=1.0, uplims=True, color="0.6", alpha=0.8,
                            label=f"3σ limit ({lab})")
        ax.set(yscale="log", xscale="log", xlabel=r"$\langle r_h \rangle$ [au]",
               ylabel=r"$Q$ [molecules s$^{-1}$]" if s == "H2O" else "")
        ax.set_title(f"Q({PRETTY[s]}):  {len(det)} robust ≥ 3σ,  {len(mar)} marginal", pad=14)
        ax.grid(alpha=0.3, which="both")
        # one entry per (status, aperture class): at most 12, below the axes
        h, l = ax.get_legend_handles_labels()
        ax.legend(h, l, fontsize=11, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=False)
    fig.suptitle(title or "Gas production rates vs heliocentric distance", y=1.02)
    fig.tight_layout()
    return fig


def plot_mixing_ratios(fits: pd.DataFrame, title: str = ""):
    fig, axes = plt.subplots(1, 2, figsize=(22, 8))
    for ax, num in zip(axes, ("CO2", "CO")):
        g = fits[(fits[f"Q_{num}_status"] == "detected") & (fits.Q_H2O_status == "detected")
                 & (fits[f"{num}_H2O"] > 0)
                 & (fits[f"Q_{num}_n_eff"] >= 2) & (fits.Q_H2O_n_eff >= 2)]
        hot = (g.h2o_source == "hot") if "h2o_source" in g else pd.Series(False, index=g.index)
        gm, gh = g[~hot], g[hot]
        ax.errorbar(gm.r_hel_mean, 100 * gm[f"{num}_H2O"], yerr=100 * gm[f"{num}_H2O_err"], fmt="o",
                    ms=11, lw=1.5, capsize=4, color="tab:purple", label="H$_2$O from the 2.7 µm band")
        if len(gh):
            ax.errorbar(gh.r_hel_mean, 100 * gh[f"{num}_H2O"], yerr=100 * gh[f"{num}_H2O_err"], fmt="o",
                        ms=11, lw=1.5, capsize=4, color="tab:purple", mfc="none", mew=2,
                        label="H$_2$O from hot bands (2.7 µm not covered)")
        ax.legend(fontsize=13, loc="lower right")
        ax.axhspan(2, 20, color="0.85", zorder=0)
        ax.text(0.02, 0.95, "2-20 %: the range other comets occupy", transform=ax.transAxes,
                fontsize=14, color="0.35", va="top")
        ax.set(xscale="log", yscale="log", xlabel=r"$\langle r_h \rangle$ [au]",
               ylabel=f"{PRETTY[num]}/H$_2$O  [%]")
        ax.set_title(f"{PRETTY[num]}/H$_2$O   ({len(g)} joint robust detections)", pad=14)
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(title or "Mixing ratios", y=1.02)
    fig.tight_layout()
    return fig


def plot_flag_comparison(census: pd.DataFrame, pairs: pd.DataFrame, paired: pd.DataFrame,
                         ref: str):
    """Flag-policy study: sample size, verdicts, and Q ratios of the groups kept by both."""
    fig, axes = plt.subplots(1, 3, figsize=(26, 8))
    ax = axes[0]
    x = np.arange(len(census))
    for i, (v, c) in enumerate((("PASS", "tab:green"), ("WARN", "tab:orange"), ("FAIL", "tab:red"))):
        bottom = census[[f"n_{w}" for w in ("PASS", "WARN", "FAIL")[:i]]].sum(axis=1) if i else 0
        ax.bar(x, census[f"n_{v}"], bottom=bottom, color=c, label=v)
    ax.set_xticks(x)
    ax.set_xticklabels(census.variant, rotation=20)
    ax.set(ylabel="band-rows")
    ax.set_title("continuum verdicts per policy", pad=14)
    ax.legend(fontsize=13)
    ax = axes[1]
    w = 0.25
    for i, s in enumerate(SPECIES):
        ax.bar(x + (i - 1) * w, census[f"n_robust_{s}"], w, color=SP_COLORS[s], label=PRETTY[s])
    ax.set_xticks(x)
    ax.set_xticklabels(census.variant, rotation=20)
    ax.set(ylabel="robust detections ($n_{eff} \\geq 2$)")
    ax.set_title("detections per policy", pad=14)
    ax.legend(fontsize=13)
    ax = axes[2]
    if len(paired):
        d = paired[(paired.status_ref == "detected") & (paired.status_oth == "detected")]
        rng = np.random.default_rng(0)
        ticks, labels = [], []
        # one slot per (comparison, species), a gap slot between comparisons; the tick
        # positions are the data positions, so labels cannot drift against the points
        for i, (cmp_, g) in enumerate(d.groupby("comparison", sort=False)):
            for j, s in enumerate(SPECIES):
                x = i * 4 + j
                gg = g[g.species == s]
                if len(gg):
                    ax.scatter(np.full(len(gg), x) + rng.normal(0, 0.08, len(gg)), gg.ratio,
                               s=30, color=SP_COLORS[s], alpha=0.7,
                               label=PRETTY[s] if i == 0 else None)
                ticks.append(x)
                labels.append(f"{cmp_.split(' vs ')[0]}\n{PRETTY[s]}")
        ax.axhline(1, color="k", ls="--")
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, rotation=90, fontsize=11)
        ax.set(yscale="log", ylabel=f"Q / Q({ref})")
        ax.set_title("Q of groups detected under both policies", pad=14)
        ax.legend(fontsize=12)
    fig.suptitle(f"Source-flag contamination study (reference: {ref})", y=1.02)
    fig.tight_layout()
    return fig


def plot_distcorr_effect(paired: pd.DataFrame, band_rows: pd.DataFrame, dc: str, raw: str):
    """Q(dc-continuum)/Q(physical-continuum) against the in-group spread of r_h^2 Delta^2."""
    fig, axes = plt.subplots(1, 3, figsize=(26, 8))
    d = paired[(paired.status_ref == "detected") & (paired.status_oth == "detected")]
    ax = axes[0]
    for s in SPECIES:
        g = d[d.species == s]
        if len(g):
            ax.scatter(100 * g.dc_spread, g.ratio, s=40, color=SP_COLORS[s], alpha=0.75,
                       label=f"{PRETTY[s]} (n={len(g)})")
    ax.axhline(1, color="k", ls="--")
    ax.set(xlabel=r"spread of $r_h^2\Delta^2$ inside the group [%]", ylabel=f"Q({dc}) / Q({raw})",
           yscale="log", xscale="log")
    ax.set_title("Q change vs geometric spread", pad=14)
    ax.legend(fontsize=13)
    ax = axes[1]
    bins = np.linspace(-6, 6, 37)
    for s in SPECIES:
        g = d[d.species == s]
        if len(g):
            ax.hist(g.nsig.clip(-6, 6), bins=bins, histtype="step", lw=2.2, color=SP_COLORS[s],
                    label=f"{PRETTY[s]}: median {g.nsig.median():+.2f}")
    ax.axvline(0, color="k", ls="--")
    ax.set(xlabel=r"$(Q_{dc} - Q_{raw}) / \sigma$", ylabel="groups")
    ax.set_title("significance of the change", pad=14)
    ax.legend(fontsize=13)
    ax = axes[2]
    if len(band_rows):
        for b in BAND_WINDOWS:
            g = band_rows[(band_rows.band == b) & np.isfinite(band_rows.snr_ratio)]
            if len(g):
                ax.hist(g.snr_ratio.clip(0, 3), bins=np.linspace(0, 3, 31), histtype="step", lw=2.2,
                        color=BAND_COLORS[b], label=f"{b}: median {g.snr_ratio.median():.2f}")
    ax.axvline(1, color="k", ls="--")
    ax.set(xlabel="peak S/N (dc) / peak S/N (physical)", ylabel="band-rows")
    ax.set_title("continuum-subtracted peak S/N", pad=14)
    ax.legend(fontsize=13)
    fig.suptitle(f"Effect of fitting the continuum in distance-corrected space ({dc} vs {raw})",
                 y=1.02)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------ batch writers
def save_grouping_figures(epochs: Dict[str, pd.DataFrame], group_map: pd.DataFrame,
                          assignment: pd.DataFrame, aperture_cfg, variant: Variant,
                          max_targets: Optional[int] = None) -> List[Path]:
    """Per-target before/after grouping figures plus the catalog summary."""
    from .dataio import load_apphot, phase_spectra
    from .grouping import group_summary
    apply_rcparams()
    out_dir = _dir.FIG_DIR / "phase_group"
    written = []
    for i, (t, ep) in enumerate(epochs.items()):
        if max_targets and i >= max_targets:
            break
        apt, parts = phase_spectra(t, aperture_cfg, variant, assignment)
        spec = pd.concat([s for _, s in parts], ignore_index=True).sort_values("wl") if parts \
            else pd.DataFrame()
        if len(spec) == 0:
            continue
        fig = plot_phase_comparison(t, ep, spec, sorted(apt.r_ap_km.unique()))
        written.append(savefig(fig, out_dir / f"{t}_phase_group.png", dpi=150))
    rows = [r for t, ep in epochs.items() for r in group_summary(t, ep)]
    g = pd.DataFrame(rows)
    fig = plot_phase_summary(g[g.grouping == "phase"], g[g.grouping == "epoch"])
    written.append(savefig(fig, out_dir / "summary_phase_group.png", dpi=DPI_ONE))
    log.info("%d grouping figures -> %s", len(written), out_dir)
    return written


def save_variant_figures(variant: Variant, assignment: pd.DataFrame,
                         max_groups: Optional[int] = None, validation: bool = True,
                         fits: bool = True) -> List[Path]:
    """Continuum validation grids and emission-model figures for one variant."""
    from .analysis import load_fits
    from .continuum import insufficient, process_group
    from .dataio import aperture_label, load_fit_input, phase_spectra
    from .fitting import fit_production_rates, model_curves
    from .revisions import revision_for
    from .config import ModelParams
    apply_rcparams()
    dirs = _dir.variant_dirs(variant.name)
    space = "distcorr" if variant.use_distcorr else "physical"
    written = []
    n = 0
    targets = sorted(assignment.target.unique())
    for t in targets:
        _, parts = phase_spectra(t, variant.aperture, variant, assignment)
        for r_ap, raw in parts:
            ph = int(raw.phase.iloc[0])
            raw = raw.reset_index(drop=True)
            if insufficient(raw, variant.continuum) is not None:
                continue
            if max_groups and n >= max_groups:
                break
            n += 1
            stem = f"{t}_{aperture_label(r_ap)}km_ph{ph}"
            rev = revision_for(t, ph, variant.revisions)
            if validation:
                out = process_group(raw, t, r_ap, int(ph), int(raw.epoch.iloc[0]),
                                    variant.continuum, space, rev=rev)
                fig = plot_validation_grid(out, t, r_ap, ph, variant.continuum, space)
                written.append(savefig(fig, dirs["fig_cont"] / f"{stem}_validation.png"))
                fig, _ = plot_raw_spectrum(raw, t, r_ap, ph, space)
                written.append(savefig(fig, dirs["fig_cont"] / f"{stem}_raw.png"))
            if fits:
                try:
                    pts = load_fit_input(variant.name, t, r_ap, int(ph), variant.fit)
                    params = ModelParams(rho_ap_km=r_ap)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        f = fit_production_rates(pts, params, variant.fit, rev=rev)
                        c = model_curves(f, pts, params)
                except (ValueError, FileNotFoundError) as exc:
                    # the pipeline lists these groups in not_fitted.csv: no fit, no fit figure
                    log.info("[%s] %s: no fit figure (%s)", variant.name, stem, str(exc)[:90])
                    continue
                fig = plot_fit(f, pts, c, variant.fit)
                written.append(savefig(fig, dirs["fig_fit"] / f"{stem}.png"))
    F = load_fits(variant.name)
    if len(F):
        fig = plot_summary_Q(F, f"{variant.name}: {variant.flags.describe()}, continuum in {space} space")
        written.append(savefig(fig, dirs["fig"] / "summary_Q_vs_rhel.png", dpi=DPI_ONE))
        fig = plot_mixing_ratios(F, f"{variant.name}")
        written.append(savefig(fig, dirs["fig"] / "summary_mixing_ratios.png", dpi=DPI_ONE))
    log.info("[%s] %d figures", variant.name, len(written))
    return written
