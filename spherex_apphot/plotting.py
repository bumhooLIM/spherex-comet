"""
Figures.

**This module is imported by notebooks only.**  ``main.py`` never touches it, so
a batch run writes photometry and nothing else -- the primitive pipeline spent a
large fraction of its wall time rendering 200-inch-wide PNGs inside the
production loop.

Keeping the drawing code here rather than pasted into the notebook means the
figures are version-controlled, reviewable, and identical between sessions;
the notebook stays a narrative that calls them.

Styling follows ``notebooks/rcparams.py``, which is imported for its side
effects as the project convention requires.  A batch of hundreds of figures
should pass ``dpi=50``; the project default of 200 is for one-off plots.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from matplotlib.gridspec import GridSpec

from . import directory as _dir
from .logging_utils import get_logger

__all__ = [
    "apply_rcparams", "savefig", "FLAG_MARKERS", "VOLATILE_BANDS", "bands_from_config",
    "plot_coverage", "plot_spectrum", "plot_reflectance", "plot_spectrum_by_flag",
    "plot_cutout_grid", "plot_cutout_with_gaia", "plot_growth_curve",
    "plot_stacks", "plot_radial_profile", "plot_aperture_validity",
    "plot_sourceflag_summary", "plot_contamination_analysis",
    "save_cutout_figures", "save_target_figures", "plot_flag_survey",
]

log = get_logger("plotting")

#: Marker per source flag, for the spectrum plots.  ``badphot`` takes precedence
#: over any contamination flag: it is a data-quality failure rather than a
#: contamination warning, so it is drawn as a cross whatever else is set.
FLAG_MARKERS = {
    "0": ("o", "tab:green",  "flag 0 (clean)"),
    "a": ("^", "tab:red",    "flag a (bright blend)"),
    "b": ("v", "tab:orange", "flag b (flux-ratio blend)"),
    "c": ("s", "gold",       "flag c (any source)"),
    "d": ("D", "tab:blue",   "flag d (SNR < 1)"),
}
#: Marker for rows with a bad pixel inside the aperture.
BADPHOT_MARKER = ("x", "black", "badphot (bad pixel in aperture)")

#: Fallback band shading, used only when no ``Config`` is available.
#: Prefer :func:`bands_from_config`, which reads ``Config.stack_bands`` so the
#: shaded regions cannot drift from the bands actually being stacked.
VOLATILE_BANDS: Tuple[Tuple[str, float, float], ...] = (
    (r"H$_2$O (2.7 $\mu$m)", 2.55, 2.8),
    (r"CO$_2$ (4.3 $\mu$m)", 4.15, 4.35),
    (r"CO (4.7 $\mu$m)", 4.6, 4.8),
)

#: Band names that are continuum rather than a volatile species.
_CONTINUUM_BANDS = {"continuum", "cont", "dust_cont_1", "dust_cont_2", "dust"}


def bands_from_config(cfg=None) -> Tuple[Tuple[str, float, float], ...]:
    """
    Band shading taken from ``cfg.stack_bands``.

    Parameters
    ----------
    cfg : Config, optional
        ``None`` returns :data:`VOLATILE_BANDS`.

    Returns
    -------
    tuple of (label, lo, hi)
        Continuum bands are dropped -- shading them would cover most of the
        spectrum and obscure the features the plot exists to show.

    Notes
    -----
    Reading the bands from the config is what stops the shaded regions drifting
    away from the bands actually being stacked when someone edits
    ``Config.stack_bands``.
    """
    if cfg is None or not getattr(cfg, "stack_bands", None):
        return VOLATILE_BANDS
    pretty = {"h2o": r"H$_2$O", "co2": r"CO$_2$", "co": "CO",
              "h2o_ice": r"H$_2$O ice", "ch4": r"CH$_4$", "oh": "OH"}
    out = []
    for name, (lo, hi) in cfg.stack_bands.items():
        key = str(name).strip().lower()
        if key in _CONTINUUM_BANDS:
            continue
        label = pretty.get(key, str(name))
        out.append((rf"{label} ({lo:g}-{hi:g} $\mu$m)", float(lo), float(hi)))
    return tuple(out) or VOLATILE_BANDS


def apply_rcparams() -> bool:
    """
    Apply the project matplotlib style from ``notebooks/rcparams.py``.

    Returns
    -------
    bool
        ``True`` if the module was found and imported.

    Notes
    -----
    ``rcparams.py`` lives in ``notebooks/`` and is imported for its side
    effects, per the project convention of not restating the settings.
    """
    nb = str(_dir.NOTEBOOK_DIR)
    if nb not in sys.path:
        sys.path.insert(0, nb)
    try:
        import rcparams  # noqa: F401
        return True
    except ImportError:
        log.warning("notebooks/rcparams.py not importable; using matplotlib defaults")
        return False


def savefig(fig, path: Path, dpi: Optional[int] = None) -> Path:
    """Save a figure, creating the parent directory, and close it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", **({"dpi": dpi} if dpi else {}))
    plt.close(fig)
    log.info("figure written: %s", path)
    return path


# ---------------------------------------------------------------------------
def plot_coverage(work: pd.DataFrame, gaia=None, objdesig: str = "",
                  max_footprints: int = 400, figsize_per_epoch: float = 7.0):
    """
    Sky coverage per epoch: cutout footprints, target track and Gaia field.

    Parameters
    ----------
    work : pandas.DataFrame
        Epoch-grouped index rows; needs ``epoch``, ``ra``, ``dec`` and the WCS
        keywords used by :func:`wcsutil.reconstruct_cutout_wcs`.
    gaia : GaiaSubset, optional
        Reference stars to draw behind the footprints.
    objdesig : str
        Title.
    max_footprints : int
        Footprints drawn per epoch.  A comet can have thousands of exposures at
        essentially the same pointing; drawing every one is slow and adds
        nothing, so they are evenly subsampled.
    figsize_per_epoch : float
        Width allotted to each epoch panel [inch].

    Returns
    -------
    matplotlib.figure.Figure
    """
    from .wcsutil import footprint_radec, reconstruct_cutout_wcs

    epochs = sorted(work["epoch"].unique())
    n = len(epochs)
    fig = plt.figure(figsize=(figsize_per_epoch * n, figsize_per_epoch))
    gs = GridSpec(1, n, figure=fig)

    for i, epoch in enumerate(epochs):
        ax = fig.add_subplot(gs[i])
        sub = work[work["epoch"] == epoch]

        if gaia is not None and len(gaia):
            ax.scatter(gaia.ra, gaia.dec, s=8, alpha=0.4, color="tab:green",
                       label="Gaia" if i == 0 else None)

        step = max(1, len(sub) // max_footprints)
        for j, row in enumerate(sub.iloc[::step].itertuples(index=False)):
            try:
                wcs = reconstruct_cutout_wcs(row)
                size = int(getattr(row, "cutout_size"))
                ra, dec = footprint_radec(wcs, size, size)
            except Exception:                   # noqa: BLE001 - a bad row must not kill the plot
                continue
            ax.plot(ra, dec, color="tab:red", lw=1.0, alpha=0.25,
                    label="cutout FoV" if (i == 0 and j == 0) else None)

        ax.scatter(sub["ra"], sub["dec"], s=25, marker="*", color="gold",
                   edgecolor="black", linewidth=0.4, zorder=5,
                   label="target" if i == 0 else None)

        if "date_obs" in sub.columns:
            d = pd.to_datetime(sub["date_obs"], errors="coerce", utc=True)
            span = f"\n{d.min():%Y-%m-%d} to {d.max():%Y-%m-%d}" if d.notna().any() else ""
        else:
            span = ""
        ax.set_title(f"epoch {epoch} (N={len(sub)}){span}")
        ax.set_xlabel(r"$\alpha$ (deg)")
        ax.set_ylabel(r"$\delta$ (deg)")
        ax.grid(True, ls="--", alpha=0.3)
        if i == 0:
            ax.legend(loc="best")

    fig.suptitle(objdesig, y=1.02)
    fig.tight_layout()
    return fig


def plot_spectrum(df: pd.DataFrame, objdesig: str = "", epoch=None,
                  value: str = "source_sum_mjy", error: str = "source_sum_err_mjy",
                  color_by: str = "r_hel", binned: Optional[pd.DataFrame] = None,
                  ax=None, cfg=None):
    """
    Flux (or any column) against wavelength, with the volatile bands marked.

    Parameters
    ----------
    df : pandas.DataFrame
        Photometry for a single aperture and epoch.
    objdesig, epoch : str, int
        Title parts.
    value, error : str
        Columns to plot.
    color_by : str
        Column mapped to the colour bar.
    binned : pandas.DataFrame, optional
        Output of :func:`reflectance.bin_spectrum`, overplotted as the
        inverse-variance weighted estimate of the underlying spectrum.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    Negative fluxes are plotted, not hidden.  They are real measurements of a
    faint source and suppressing them biases the eye exactly as dropping them
    biases the mean (review item S5).
    """
    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(15, 6))
    else:
        fig = ax.figure

    d = df.sort_values("wl")
    ax.errorbar(d["wl"], d[value], yerr=d.get(error), fmt="none",
                ecolor="gray", alpha=0.5, zorder=2)
    sc = ax.scatter(d["wl"], d[value], c=d[color_by] if color_by in d else None,
                    cmap="RdBu", s=28, zorder=3)
    if color_by in d:
        cbar = fig.colorbar(sc, ax=ax, pad=0.02)
        cbar.set_label({"r_hel": r"$r_{\rm hel}$ (au)"}.get(color_by, color_by))

    if binned is not None and not binned.empty:
        ax.errorbar(binned["wl"], binned["value"], yerr=binned["error"],
                    fmt="o-", color="black", ms=4, lw=1.4, zorder=5,
                    label="binned (inv-var weighted)")

    for label, lo, hi in bands_from_config(cfg):
        ax.axvspan(lo, hi, color="skyblue", alpha=0.25, zorder=0, label=label)

    ax.axhline(0.0, color="k", lw=0.8, ls=":", zorder=1)
    ax.set_xlabel(r"wavelength ($\mu$m)")
    ax.set_ylabel({"source_sum_mjy": "aperture flux (mJy)"}.get(value, value))
    title = objdesig if epoch is None else f"{objdesig} | epoch {epoch}"
    if "ap_label" in d.columns and d["ap_label"].nunique() == 1:
        title += f" | {d['ap_label'].iloc[0]}"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=10)
    if created:
        fig.tight_layout()
    return fig


def plot_reflectance(df: pd.DataFrame, objdesig: str = "", epoch=None,
                     binned: Optional[pd.DataFrame] = None, ylim=(0.0, 2.0),
                     cfg=None):
    """Normalised reflectance against wavelength.  See :func:`plot_spectrum`."""
    fig, ax = plt.subplots(figsize=(15, 6))
    plot_spectrum(df, objdesig=objdesig, epoch=epoch, value="refl_norm",
                  error="refl_norm_err", binned=binned, ax=ax, cfg=cfg)
    ax.set_ylabel("normalised reflectance")
    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    if ylim:
        ax.set_ylim(*ylim)
    fig.tight_layout()
    return fig


def plot_cutout_grid(df: pd.DataFrame, resolver, cfg, objdesig: str = "",
                     ncols: int = 5, max_panels: int = 25, zoom_pix: int = 20):
    """
    A grid of individual cutouts with the aperture and annulus drawn on.

    Parameters
    ----------
    df : pandas.DataFrame
        Photometry rows, one per panel (already reduced to a single aperture).
    resolver : FitsResolver
    cfg : Config
    objdesig : str
    ncols : int
        Panels per row.
    max_panels : int
        Hard cap.  The primitive version drew one panel per exposure, producing
        a 20x160 inch canvas at 200 dpi for a busy epoch.
    zoom_pix : int
        Half-width of the displayed region [pixel].

    Returns
    -------
    matplotlib.figure.Figure
    """
    from .fitsio import read_cutout
    from .masking import build_badpix_mask

    d = df.head(int(max_panels))
    n = len(d)
    if n == 0:
        log.warning("nothing to draw")
        return plt.figure()

    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for i, row in enumerate(d.itertuples(index=False)):
        ax = axes[i]
        path = resolver.find(getattr(row, "filename"), objdesig)
        if path is None:
            ax.set_title("file missing", fontsize=8, color="red")
            ax.axis("off")
            continue
        cut = read_cutout(path)
        badpix, _ = build_badpix_mask(cut.sci, cut.var, cut.flag, cfg.bad_flag_bits,
                                      mask_nonfinite_sci=cfg.mask_nonfinite_sci,
                                      mask_bad_variance=cfg.mask_bad_variance)
        img = np.where(badpix, np.nan, cut.sci)

        xc, yc = float(getattr(row, "xcen")), float(getattr(row, "ycen"))
        try:
            vmin, vmax = ZScaleInterval().get_limits(img[np.isfinite(img)])
        except (ValueError, IndexError):
            vmin, vmax = np.nanpercentile(img, [5, 95])

        # extent uses pixel *edges*, so a pixel centre at integer i spans
        # [i-0.5, i+0.5]; without the half-pixel the overlaid circles sit
        # half a pixel off the data they describe.
        ny, nx = img.shape
        ax.imshow(img, origin="lower", cmap="magma", vmin=vmin, vmax=vmax,
                  extent=(-0.5, nx - 0.5, -0.5, ny - 0.5))
        ax.add_patch(plt.Circle((xc, yc), getattr(row, "r_ap_pix"),
                                fill=False, color="red", lw=1.4))
        ax.add_patch(plt.Circle((xc, yc), cfg.r_in_pix, fill=False,
                                color="cyan", lw=1.0, ls="--"))
        ax.add_patch(plt.Circle((xc, yc), cfg.r_out_pix, fill=False,
                                color="cyan", lw=1.0, ls="--"))
        ax.set_xlim(xc - zoom_pix, xc + zoom_pix)
        ax.set_ylim(yc - zoom_pix, yc + zoom_pix)
        ax.set_title(f"{getattr(row, 'wl'):.2f} $\\mu$m | flag {getattr(row, 'sourceflag')}",
                     fontsize=9)
        ax.tick_params(labelsize=7)

    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(f"{objdesig} cutouts (red: aperture, cyan: sky annulus)", y=1.0)
    fig.tight_layout()
    return fig


def plot_stacks(stacks: Dict[tuple, object], objdesig: str = "", epoch=None,
                cfg=None):
    """
    Band stacks for one epoch, with the count map beneath each panel.

    Parameters
    ----------
    stacks : dict
        ``{(epoch, band): StackResult}`` from :func:`pipeline.stack_target`.
    objdesig : str
    epoch : int, optional
        Restrict to one epoch; defaults to the first present.
    cfg : Config, optional
        Used to draw the reference aperture radius.

    Returns
    -------
    matplotlib.figure.Figure
    """
    keys = sorted(k for k in stacks if epoch is None or k[0] == epoch)
    if not keys:
        log.warning("no stacks to draw")
        return plt.figure()
    epoch = keys[0][0]
    keys = [k for k in keys if k[0] == epoch]

    n = len(keys)
    fig, axes = plt.subplots(2, n, figsize=(3.4 * n, 7.0),
                             gridspec_kw={"height_ratios": [3, 2]}, squeeze=False)
    for j, key in enumerate(keys):
        res = stacks[key]
        r = res.radius
        ext = (-r - 0.5, r + 0.5, -r - 0.5, r + 0.5)

        ax = axes[0][j]
        finite = res.image[np.isfinite(res.image)]
        try:
            vmin, vmax = ZScaleInterval().get_limits(finite)
        except (ValueError, IndexError):
            vmin, vmax = (np.nanpercentile(finite, [5, 95]) if finite.size else (0, 1))
        ax.imshow(res.image, origin="lower", cmap="magma", vmin=vmin, vmax=vmax, extent=ext)
        if cfg is not None:
            ax.add_patch(plt.Circle((0, 0), cfg.aperture_radii_pix[0],
                                    fill=False, color="red", lw=1.3))
        ax.set_title(f"{key[1]}\nN={res.n_frames}", fontsize=11)
        ax.set_xlabel(r"$\Delta x$ (pix)")
        if j == 0:
            ax.set_ylabel(r"$\Delta y$ (pix)")

        axc = axes[1][j]
        im = axc.imshow(res.count, origin="lower", cmap="viridis", extent=ext)
        fig.colorbar(im, ax=axc, pad=0.02, fraction=0.046)
        axc.set_title("frames per pixel", fontsize=9)
        axc.tick_params(labelsize=7)

    fig.suptitle(f"{objdesig} | epoch {epoch} | sigma-clipped mean, "
                 f"sky-subtracted, sub-pixel registered", y=1.01)
    fig.tight_layout()
    return fig


def plot_aperture_validity(df: pd.DataFrame, cfg, objdesig: str = ""):
    """
    Which apertures survived, as a function of observer distance.

    Makes the ``PSF_FWHM <= r_ap < r_in`` acceptance window visible, and shows
    directly why a nearby target loses its large apertures and a distant one
    loses its small ones.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    km = df[df["ap_kind"] == "km"]
    sc = ax.scatter(km["r_obs"], km["r_ap_pix"], c=km["r_ap_km"] / 1e3,
                    cmap="viridis", s=10, alpha=0.7)
    fig.colorbar(sc, ax=ax, pad=0.02).set_label(r"$r_{\rm ap}$ ($10^3$ km)")

    ax.axhline(cfg.r_in_pix, color="red", ls="--", lw=2,
               label=fr"$r_{{\rm in}} = {cfg.r_in_pix:g}$ pix (upper limit)")
    if "psf_fwhm_pix" in df.columns:
        ax.axhline(float(df["psf_fwhm_pix"].median()), color="orange", ls="--", lw=2,
                   label="median PSF FWHM (lower limit)")
    for r in cfg.aperture_radii_pix:
        ax.axhline(r, color="gray", ls=":", lw=1)

    ax.set_xlabel(r"$r_{\rm obs}$ (au)")
    ax.set_ylabel(r"$r_{\rm ap}$ (pix)")
    ax.set_yscale("log")
    ax.set_title(f"{objdesig}: aperture radii retained")
    ax.legend(fontsize=11)
    fig.tight_layout()
    return fig


def plot_sourceflag_summary(df: pd.DataFrame, objdesig: str = ""):
    """
    Contamination-flag census and where the flags land in wavelength.

    Returns
    -------
    matplotlib.figure.Figure
    """
    order = ["a", "b", "c", "d", "0"]
    colors = {"a": "tab:red", "b": "tab:orange", "c": "gold",
              "d": "tab:blue", "0": "tab:green"}
    present = [f for f in order if (df["sourceflag"].astype(str) == f).any()]

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    counts = [int((df["sourceflag"].astype(str) == f).sum()) for f in present]
    axes[0].bar(present, counts, color=[colors[f] for f in present])
    axes[0].set_xlabel("sourceflag")
    axes[0].set_ylabel("rows")
    axes[0].set_title("flag census")
    for x, c in zip(present, counts):
        axes[0].text(x, c, f"{c}", ha="center", va="bottom", fontsize=11)

    for f in present:
        sel = df["sourceflag"].astype(str) == f
        axes[1].hist(df.loc[sel, "wl"], bins=40, histtype="step", lw=2,
                     color=colors[f], label=f"'{f}' (n={int(sel.sum())})")
    axes[1].set_xlabel(r"wavelength ($\mu$m)")
    axes[1].set_ylabel("rows")
    axes[1].set_title("flags vs wavelength")
    axes[1].legend(fontsize=11)

    fig.suptitle(f"{objdesig}: Gaia contamination flags "
                 "(a: bright blend, b: vmag test, c: any source, d: SNR<1)", y=1.02)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Per-cutout figures
# ---------------------------------------------------------------------------
def _gaia_marker_sizes(gmag, ref: float = 12.0, base: float = 40.0,
                       lo: float = 12.0, hi: float = 900.0) -> np.ndarray:
    """
    Marker area proportional to Gaia flux, clipped to a drawable range.

    Sizing by flux rather than by magnitude is what makes the brightness
    *visible*: a G = 10 star should look overwhelmingly larger than a G = 17
    one, which is exactly the 10**(-0.4 dG) ratio.
    """
    g = np.asarray(gmag, dtype=float)
    size = base * np.power(10.0, -0.4 * (g - ref))
    return np.clip(np.nan_to_num(size, nan=lo), lo, hi)


def plot_cutout_with_gaia(cut, badpix, row, gaia_field, cfg, *,
                          objdesig: str = "", r_ap_pix: Optional[float] = None,
                          show_ids: bool = True):
    """
    Three-panel view of one cutout: raw science, flag plane, masked + Gaia.

    Parameters
    ----------
    cut : Cutout
        From :func:`fitsio.read_cutout`.
    badpix : ndarray of bool
        Bad-pixel mask for this cutout.
    row : pandas.Series or namedtuple
        The photometry/index row supplying ``xcen``, ``ycen``, ``wl`` and the
        source-flag columns.
    gaia_field : GaiaPixelField
        Projected Gaia sources, from :func:`sourceflag.project_gaia`.
    cfg : Config
    objdesig : str
    r_ap_pix : float, optional
        Aperture to draw.  Defaults to ``row.r_ap_pix`` when present.
    show_ids : bool
        Annotate each Gaia source with its G magnitude.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    The right-hand panel is the one that matters: it shows what the photometry
    actually integrated, with every catalogued star marked at the size its flux
    warrants.  Because the pipeline does not mask stars, this is the only place a
    blend becomes visible before it turns into an anomalous flux.
    """
    def _get(name, default=np.nan):
        if hasattr(row, name):
            return getattr(row, name)
        try:
            return row[name]
        except Exception:                       # noqa: BLE001
            return default

    xc, yc = float(_get("xcen")), float(_get("ycen"))
    r_ap = float(r_ap_pix if r_ap_pix is not None else _get("r_ap_pix", 2.0))
    ny, nx = cut.sci.shape
    ext = (-0.5, nx - 0.5, -0.5, ny - 0.5)

    masked = np.where(badpix, np.nan, cut.sci)
    finite = masked[np.isfinite(masked)]
    try:
        vmin, vmax = ZScaleInterval().get_limits(finite)
    except (ValueError, IndexError):
        vmin, vmax = (np.nanpercentile(finite, [5, 95]) if finite.size else (0.0, 1.0))

    fig, axes = plt.subplots(1, 3, figsize=(19, 6.4))

    axes[0].imshow(cut.sci, origin="lower", cmap="magma", vmin=vmin, vmax=vmax, extent=ext)
    axes[0].set_title("science (raw)")

    # log2 of the flag value keeps the many single-bit values distinguishable.
    with np.errstate(divide="ignore"):
        flagshow = np.where(cut.flag > 0, np.log2(cut.flag.astype(float)), np.nan)
    im = axes[1].imshow(flagshow, origin="lower", cmap="viridis", extent=ext)
    fig.colorbar(im, ax=axes[1], pad=0.02, fraction=0.046).set_label(r"$\log_2$(FLAG)")
    axes[1].set_title(f"FLAG plane ({int(badpix.sum())} px bad)")

    axes[2].imshow(masked, origin="lower", cmap="magma", vmin=vmin, vmax=vmax, extent=ext)
    axes[2].set_title("masked + Gaia sources")

    if gaia_field is not None and len(gaia_field):
        sizes = _gaia_marker_sizes(gaia_field.gmag)
        sc = axes[2].scatter(gaia_field.x, gaia_field.y, s=sizes, facecolors="none",
                             edgecolors="cyan", linewidths=1.4, zorder=6)
        if show_ids:
            for gx, gy, gm in zip(gaia_field.x, gaia_field.y, gaia_field.gmag):
                if -0.5 <= gx <= nx and -0.5 <= gy <= ny:
                    axes[2].annotate(f"{gm:.1f}", (gx, gy), xytext=(4, 4),
                                     textcoords="offset points", fontsize=8,
                                     color="cyan")

    for ax in axes:
        ax.add_patch(plt.Circle((xc, yc), r_ap, fill=False, color="red", lw=1.8))
        ax.add_patch(plt.Circle((xc, yc), cfg.r_in_pix, fill=False, color="deepskyblue",
                                lw=1.2, ls="--"))
        ax.add_patch(plt.Circle((xc, yc), cfg.r_out_pix, fill=False, color="deepskyblue",
                                lw=1.2, ls="--"))
        ax.set_xlim(-0.5, nx - 0.5); ax.set_ylim(-0.5, ny - 0.5)
        ax.set_xlabel("x (pix)")
    axes[0].set_ylabel("y (pix)")

    flag = _get("sourceflag", "?")
    ngaia = _get("n_gaia", np.nan)
    geff = _get("gmag_eff", np.nan)
    fig.suptitle(
        rf"{objdesig}  {_get('wl', float('nan')):.3f} $\mu$m  |  "
        rf"flag '{flag}'  n_gaia={ngaia:.0f}  $G_{{\rm eff}}$={geff:.2f}  |  "
        rf"red: $r_{{\rm ap}}$={r_ap:.2f} px, blue: sky annulus, "
        rf"cyan: Gaia (size $\propto$ flux)", y=1.02, fontsize=14)
    fig.tight_layout()
    return fig


def plot_growth_curve(rows: pd.DataFrame, objdesig: str = "", ax=None,
                      annotate_flag: bool = True):
    """
    Enclosed flux and annulus surface brightness against aperture radius.

    Parameters
    ----------
    rows : pandas.DataFrame
        All aperture rows of a single exposure.
    objdesig : str
    ax : matplotlib.axes.Axes, optional
        Ignored when the two-panel default layout is used.
    annotate_flag : bool
        Mark the radius of the nearest Gaia source.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    The lower panel is the diagnostic one.  For a steady-state coma the annulus
    surface brightness falls monotonically outward; a field star entering the
    aperture makes it *rise*, which is the signature
    :func:`diagnostics.growth_metrics` quantifies.
    """
    g = rows[rows["ap_kind"] == "km"].sort_values("r_ap_pix")
    if g.empty:
        g = rows.sort_values("r_ap_pix")
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2]})

    axes[0].errorbar(g["r_ap_pix"], g["source_sum_mjy"], yerr=g["source_sum_err_mjy"],
                     marker="o", ms=5, lw=1.4, capsize=3, color="tab:blue")
    bad = g[g["badphot"].astype(bool)]
    if len(bad):
        axes[0].scatter(bad["r_ap_pix"], bad["source_sum_mjy"], marker="x", s=70,
                        color="black", zorder=6, label="badphot")
        axes[0].legend(fontsize=11)
    axes[0].axhline(0.0, color="k", lw=0.8, ls=":")
    axes[0].set_ylabel("enclosed flux (mJy)")
    axes[0].grid(alpha=0.3)

    r = g["r_ap_pix"].to_numpy(float)
    F = g["source_sum_mjy"].to_numpy(float)
    A = g["aperture_area_eff_pix2"].to_numpy(float)
    if r.size > 1:
        dA = np.diff(A)
        ok = dA > 0
        sb = np.where(ok, np.diff(F) / np.where(ok, dA, np.nan), np.nan)
        axes[1].plot(0.5 * (r[1:] + r[:-1]), sb, marker="s", ms=4, lw=1.3,
                     color="tab:purple")
    axes[1].axhline(0.0, color="k", lw=0.8, ls=":")
    axes[1].set_xlabel(r"$r_{\rm ap}$ (pix)")
    axes[1].set_ylabel("annulus SB\n(mJy/pix)")
    axes[1].grid(alpha=0.3)

    if annotate_flag and "dist_gmag_nearest" in g.columns:
        d = float(g["dist_gmag_nearest"].iloc[-1])
        gm = float(g["gmag_nearest"].iloc[-1])
        if np.isfinite(d) and d <= float(np.nanmax(r)) * 1.5:
            for ax in axes:
                ax.axvline(d, color="tab:red", ls="--", lw=1.5, alpha=0.8)
            axes[0].annotate(f"nearest Gaia\nG={gm:.1f} at {d:.1f} px",
                             xy=(d, axes[0].get_ylim()[1]), xytext=(4, -18),
                             textcoords="offset points", color="tab:red", fontsize=11,
                             va="top")

    row0 = g.iloc[0]
    fig.suptitle(rf"{objdesig}  {row0['wl']:.3f} $\mu$m  |  flag '{row0['sourceflag']}'  "
                 rf"|  growth curve", y=1.0)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Spectra by flag
# ---------------------------------------------------------------------------
def plot_spectrum_by_flag(df: pd.DataFrame, objdesig: str = "", epoch=None,
                          value: str = "flux_distcorr_mjy",
                          error: str = "flux_distcorr_err_mjy",
                          ax=None, logy: bool = False, cfg=None):
    """
    Spectrum with a distinct marker per source flag.

    Markers follow :data:`FLAG_MARKERS`: circle for flag ``0``, upper triangle
    for ``a``, lower triangle for ``b``, square for ``c``, diamond for ``d``,
    and a cross for ``badphot`` — which takes precedence, being a data-quality
    failure rather than a contamination warning.

    Parameters
    ----------
    df : pandas.DataFrame
        Photometry for one aperture and epoch.
    objdesig, epoch : str, int
    value, error : str
        Columns to plot; defaults to the distance-corrected flux.
    ax : matplotlib.axes.Axes, optional
    logy : bool
        Symmetric-log y axis, useful when a few contaminated points dominate.

    Returns
    -------
    matplotlib.figure.Figure
    """
    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(15, 6.5))
    else:
        fig = ax.figure

    d = df.sort_values("wl")
    flags = d["sourceflag"].astype(str)
    badphot = d["badphot"].astype(bool)

    for label, sel, (marker, color, legend) in (
            [(f, (flags == f) & ~badphot, FLAG_MARKERS[f]) for f in FLAG_MARKERS]
            + [("badphot", badphot, BADPHOT_MARKER)]):
        if not sel.any():
            continue
        sub = d[sel]
        ax.errorbar(sub["wl"], sub[value],
                    yerr=sub[error] if error in sub else None,
                    fmt="none", ecolor=color, alpha=0.35, zorder=2)
        # 'x' is an unfilled marker: matplotlib ignores edgecolor on it, so it
        # takes `color` while the closed markers stay hollow for legibility
        # where points overlap.
        style = (dict(color=color) if marker == "x"
                 else dict(facecolors="none", edgecolors=color))
        ax.scatter(sub["wl"], sub[value], marker=marker, s=42, zorder=3,
                   linewidths=1.5, label=f"{legend}  (n={int(sel.sum())})", **style)

    for label, lo, hi in bands_from_config(cfg):
        ax.axvspan(lo, hi, color="skyblue", alpha=0.22, zorder=0, label=label)
    ax.axhline(0.0, color="k", lw=0.8, ls=":", zorder=1)
    if logy:
        ax.set_yscale("symlog", linthresh=np.nanpercentile(np.abs(d[value]), 20) or 1e-3)

    ax.set_xlabel(r"wavelength ($\mu$m)")
    ax.set_ylabel({
        # flux normalised to r_hel = r_obs = 1 au, so the unit is still mJy
        "flux_distcorr_mjy": r"$F \times r_{\rm hel}^2 r_{\rm obs}^2$  (mJy at 1 au)",
        "source_sum_mjy": "aperture flux (mJy)",
    }.get(value, value))
    title = objdesig if epoch is None else f"{objdesig} | epoch {epoch}"
    if "ap_label" in d.columns and d["ap_label"].nunique() == 1:
        title += f" | {d['ap_label'].iloc[0]}"
    ax.set_title(title)
    ax.legend(fontsize=10, ncol=2, loc="best")
    if created:
        fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Stacks
# ---------------------------------------------------------------------------
def radial_profile(image: np.ndarray, center=None, nbins: int = 20):
    """
    Azimuthally averaged radial profile of an image.

    Parameters
    ----------
    image : ndarray
        2-D image; NaN pixels are ignored.
    center : (x, y), optional
        Defaults to the geometric centre.
    nbins : int
        Number of radial bins out to the half-width.

    Returns
    -------
    r, value, error, n : ndarray
        Bin centre [pixel], median, standard error of the median, and the number
        of contributing pixels.
    """
    ny, nx = image.shape
    cx, cy = center if center is not None else ((nx - 1) / 2.0, (ny - 1) / 2.0)
    y, x = np.mgrid[:ny, :nx]
    rr = np.hypot(x - cx, y - cy)
    rmax = min(nx, ny) / 2.0
    edges = np.linspace(0.0, rmax, int(nbins) + 1)
    r_c, val, err, cnt = [], [], [], []
    for i in range(len(edges) - 1):
        sel = (rr >= edges[i]) & (rr < edges[i + 1]) & np.isfinite(image)
        n = int(sel.sum())
        if n == 0:
            continue
        v = image[sel]
        r_c.append(0.5 * (edges[i] + edges[i + 1]))
        val.append(float(np.median(v)))
        # sqrt(pi/2) converts the standard error of the mean to that of the median.
        err.append(float(np.std(v, ddof=1) / np.sqrt(n) * np.sqrt(np.pi / 2)) if n > 1 else np.nan)
        cnt.append(n)
    return (np.array(r_c), np.array(val), np.array(err), np.array(cnt))


def plot_radial_profile(stacks: Dict[tuple, object], objdesig: str = "", epoch=None,
                        pix_scale_arcsec: Optional[float] = None):
    r"""
    Radial profiles of the band stacks for one epoch, on one axis.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    Plotted on log-log axes, where a canonical steady-state coma
    (:math:`\Sigma \propto 1/\rho`) is a straight line of slope -1; a
    reference line of that slope is drawn for comparison.
    """
    keys = sorted(k for k in stacks if epoch is None or k[0] == epoch)
    if not keys:
        log.warning("no stacks to profile")
        return plt.figure()
    epoch = keys[0][0]
    keys = [k for k in keys if k[0] == epoch]

    fig, ax = plt.subplots(figsize=(11, 7))
    for key in keys:
        res = stacks[key]
        r, v, e, n = radial_profile(res.image, nbins=res.radius)
        ok = r > 0
        ax.errorbar(r[ok], v[ok], yerr=e[ok], marker="o", ms=4, lw=1.4, capsize=2,
                    label=f"{key[1]} (N={res.n_frames})")

    pos = [ln for ln in ax.get_lines() if len(ln.get_ydata())]
    if pos:
        rr = np.array([1.0, 20.0])
        ydata = np.concatenate([np.asarray(ln.get_ydata(), float) for ln in pos])
        scale = np.nanmax(ydata[ydata > 0]) if np.any(ydata > 0) else 1.0
        ax.plot(rr, scale / rr, color="gray", ls="--", lw=1.5,
                label=r"$\Sigma \propto \rho^{-1}$ (steady-state coma)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\rho$ (pix)")
    ax.set_ylabel("surface brightness (mJy/pix)")
    ax.set_title(f"{objdesig} | epoch {epoch} | stacked radial profiles")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=11)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Contamination analysis
# ---------------------------------------------------------------------------
def plot_contamination_analysis(metrics: pd.DataFrame, threshold: float = 3.0,
                                objdesig: str = ""):
    """
    How the source flags relate to growth-curve contamination.

    Four panels: the distribution of the contamination indicator per flag; the
    indicator against the nearest-star distance; recall and precision of each
    cumulative cut; and the surviving-sample purity.

    Parameters
    ----------
    metrics : pandas.DataFrame
        Output of :func:`diagnostics.growth_metrics`.
    threshold : float
        ``sb_rise_sigma`` above which an exposure counts as contaminated.
    objdesig : str

    Returns
    -------
    matplotlib.figure.Figure
    """
    from .diagnostics import FLAG_PRIORITY, flag_effectiveness

    eff = flag_effectiveness(metrics, threshold=threshold)
    base = eff.attrs["base_rate"]
    colors = {f: FLAG_MARKERS[f][1] for f in FLAG_MARKERS}

    fig, axes = plt.subplots(2, 2, figsize=(16, 11))

    ax = axes[0][0]
    present = [f for f in FLAG_PRIORITY if (metrics["flag_worst"] == f).any()]
    data = [metrics.loc[metrics["flag_worst"] == f, "sb_rise_sigma"].to_numpy() for f in present]
    lo = float(np.nanpercentile(metrics["sb_rise_sigma"], 1))
    hi = float(np.nanpercentile(metrics["sb_rise_sigma"], 99))
    bins = np.linspace(min(lo, -5), max(hi, 10), 45)
    for f, d in zip(present, data):
        ax.hist(d, bins=bins, histtype="step", lw=2.2, color=colors[f],
                label=f"'{f}' (n={len(d)}, median {np.nanmedian(d):+.1f})")
    ax.axvline(threshold, color="k", ls="--", lw=1.8,
               label=rf"contaminated  > {threshold:g}$\sigma$")
    ax.set_xlabel(r"$\Sigma$-rise significance ($\sigma$)")
    ax.set_ylabel("exposures")
    ax.set_title("contamination indicator by flag")
    ax.legend(fontsize=10)

    ax = axes[0][1]
    m = metrics[np.isfinite(metrics["dist_gmag_nearest"])]
    sc = ax.scatter(m["dist_gmag_nearest"], m["sb_rise_sigma"],
                    c=m["gmag_brightest_outer"], cmap="viridis_r", s=22, alpha=0.8)
    fig.colorbar(sc, ax=ax, pad=0.02).set_label("brightest G in aperture")
    ax.axhline(threshold, color="k", ls="--", lw=1.5)
    ax.set_xlabel("distance to nearest Gaia source (pix)")
    ax.set_ylabel(r"$\Sigma$-rise significance ($\sigma$)")
    ax.set_title("contamination vs. star proximity")
    ax.set_yscale("symlog", linthresh=5)

    cuts = eff[eff["cut"].str.startswith("cut")]
    ax = axes[1][0]
    xs = np.arange(len(cuts))
    ax.bar(xs - 0.2, cuts["recall"], width=0.4, label="recall", color="tab:blue")
    ax.bar(xs + 0.2, cuts["precision"], width=0.4, label="precision", color="tab:orange")
    ax.axhline(base, color="k", ls=":", lw=2, label=f"base rate ({base:.2f})")
    ax.set_xticks(xs); ax.set_xticklabels(cuts["cut"], rotation=15)
    ax.set_ylim(0, 1.05); ax.set_ylabel("fraction")
    ax.set_title("does the cut find contamination?")
    ax.legend(fontsize=11)

    ax = axes[1][1]
    ax.bar(xs, cuts["purity_kept"], width=0.6, color="tab:green")
    ax.axhline(1.0 - base, color="k", ls=":", lw=2,
               label=f"no cut ({1 - base:.3f})")
    ax.set_xticks(xs); ax.set_xticklabels(cuts["cut"], rotation=15)
    ax.set_ylabel("clean fraction of what survives")
    ax.set_ylim(min(0.9 * (1 - base), 0.5), 1.005)
    ax.set_title("purity of the retained sample")
    for x, (v, n) in enumerate(zip(cuts["purity_kept"], cuts["n_kept"])):
        ax.text(x, v, f"n={n}", ha="center", va="bottom", fontsize=10)
    ax.legend(fontsize=11)

    fig.suptitle(f"{objdesig}: source-flag effectiveness against growth-curve "
                 f"contamination ({eff.attrs['n_contaminated']}/{eff.attrs['n_total']} "
                 f"exposures contaminated)", y=1.01)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Batch figure writers
# ---------------------------------------------------------------------------
def save_cutout_figures(phot: pd.DataFrame, resolver, cfg, gaia_subset=None, *,
                        objdesig: str = "", out_dir: Optional[Path] = None,
                        ap_label: Optional[str] = None, max_cutouts: Optional[int] = None,
                        dpi: int = 50, skip_existing: bool = False) -> List[Path]:
    """
    Write the two per-cutout figures for every exposure: the three-panel
    cutout view and the growth curve.

    Parameters
    ----------
    phot : pandas.DataFrame
        Photometry table for one target.
    resolver : FitsResolver
    cfg : Config
    gaia_subset : GaiaSubset, optional
        Sources for the target region; without it no Gaia markers are drawn.
    objdesig : str
    out_dir : Path, optional
        Defaults to ``FIG_DIR/<objdesig>/cutouts``.
    ap_label : str, optional
        Aperture whose radius is drawn on the cutout.  Defaults to the smallest
        fixed-pixel aperture.
    max_cutouts : int, optional
        Cap, for interactive use.  ``None`` writes every exposure.
    dpi : int
        50 by project convention for batch output; 200 is for one-off figures.
    skip_existing : bool
        Skip an exposure whose two figures are already on disk.  Rendering the
        full working list takes hours, so a resumed run should not redo it --
        and the skip is checked before the FITS is read, which is where the time
        actually goes on an external drive.

    Returns
    -------
    list of Path
    """
    from .catalog import NeighborIndex
    from .fitsio import read_cutout
    from .masking import build_badpix_mask
    from .sourceflag import project_gaia
    from .wcsutil import wcs_from_header

    out_dir = Path(out_dir) if out_dir else (_dir.FIG_DIR / (objdesig or "target") / "cutouts")
    out_dir.mkdir(parents=True, exist_ok=True)

    if ap_label is None:
        pix = sorted(phot.loc[phot["ap_kind"] == "pix", "ap_label"].unique())
        ap_label = pix[0] if pix else str(phot["ap_label"].iloc[0])

    files = list(dict.fromkeys(phot["filename"]))
    if max_cutouts:
        files = files[:int(max_cutouts)]
    neighbors = NeighborIndex(gaia_subset) if gaia_subset is not None else None

    written: List[Path] = []
    for i, fname in enumerate(files, start=1):
        stem = Path(fname).stem
        cut_png = out_dir / f"{stem}_cutout.png"
        growth_png = out_dir / f"{stem}_growth.png"
        if skip_existing and cut_png.is_file() and growth_png.is_file():
            continue

        rows = phot[phot["filename"] == fname]
        ref = rows[rows["ap_label"] == ap_label]
        ref = ref.iloc[0] if len(ref) else rows.iloc[0]
        path = resolver.find(fname, objdesig)
        if path is None:
            log.warning("skipping figures for missing %s", fname)
            continue
        cut = read_cutout(path)
        badpix, _ = build_badpix_mask(cut.sci, cut.var, cut.flag, cfg.bad_flag_bits,
                                      mask_nonfinite_sci=cfg.mask_nonfinite_sci,
                                      mask_bad_variance=cfg.mask_bad_variance)
        field = None
        if neighbors is not None:
            ny, nx = cut.sci.shape
            rad_deg = (0.5 * np.hypot(nx, ny) + 30.0) * float(ref["pix_scale"]) / 3600.0
            near = neighbors.query(float(ref["ra"]), float(ref["dec"]), rad_deg)
            field = project_gaia(near, wcs_from_header(cut.header),
                                 float(ref["xcen"]), float(ref["ycen"]), cut.sci.shape)

        fig = plot_cutout_with_gaia(cut, badpix, ref, field, cfg, objdesig=objdesig,
                                    r_ap_pix=float(ref["r_ap_pix"]))
        written.append(savefig(fig, cut_png, dpi=dpi))

        fig = plot_growth_curve(rows, objdesig=objdesig)
        written.append(savefig(fig, growth_png, dpi=dpi))

        if i % 50 == 0:
            log.info("  per-cutout figures: %d/%d exposures", i, len(files))

    log.info("wrote %d per-cutout figures to %s", len(written), out_dir)
    return written


def save_target_figures(phot: pd.DataFrame, cfg, *, objdesig: str = "",
                        stacks: Optional[Dict[tuple, object]] = None,
                        metrics: Optional[pd.DataFrame] = None,
                        ap_label: Optional[str] = None,
                        out_dir: Optional[Path] = None,
                        dpi: int = 200) -> List[Path]:
    """
    Write the per-target summary figures: spectrum, reflectance, stacks,
    radial profiles, and — when ``metrics`` is supplied — the contamination
    analysis.

    Parameters
    ----------
    phot : pandas.DataFrame
    cfg : Config
    objdesig : str
    stacks : dict, optional
        ``{(epoch, band): StackResult}``.
    metrics : pandas.DataFrame, optional
        Output of :func:`diagnostics.growth_metrics`.
    ap_label : str, optional
    out_dir : Path, optional
        Defaults to ``FIG_DIR/<objdesig>``.
    dpi : int
        200 — these are one-off figures, not a batch.

    Returns
    -------
    list of Path
    """
    from .reflectance import add_reflectance, bin_spectrum, normalise_reflectance

    out_dir = Path(out_dir) if out_dir else (_dir.FIG_DIR / (objdesig or "target"))
    out_dir.mkdir(parents=True, exist_ok=True)

    if ap_label is None:
        pix = sorted(phot.loc[phot["ap_kind"] == "pix", "ap_label"].unique())
        ap_label = pix[0] if pix else str(phot["ap_label"].iloc[0])

    refl = normalise_reflectance(add_reflectance(phot, cfg), cfg)
    written: List[Path] = []

    for epoch in sorted(phot["epoch"].unique()):
        sel = (phot["ap_label"] == ap_label) & (phot["epoch"] == epoch)
        spec = phot[sel]
        if spec.empty:
            continue
        fig = plot_spectrum_by_flag(spec, objdesig=objdesig, epoch=epoch, cfg=cfg)
        written.append(savefig(fig, out_dir / f"spectrum_{objdesig}_epoch{epoch}.png", dpi=dpi))

        rsel = (refl["ap_label"] == ap_label) & (refl["epoch"] == epoch)
        rsub = refl[rsel]
        rbin = bin_spectrum(rsub, value_col="refl_norm", err_col="refl_norm_err", n_bins=50)
        fig = plot_reflectance(rsub, objdesig=objdesig, epoch=epoch, binned=rbin, cfg=cfg)
        written.append(savefig(fig, out_dir / f"reflectance_{objdesig}_epoch{epoch}.png", dpi=dpi))

        if stacks:
            if any(k[0] == epoch for k in stacks):
                fig = plot_stacks(stacks, objdesig=objdesig, epoch=epoch, cfg=cfg)
                written.append(savefig(fig, out_dir / f"stacks_{objdesig}_epoch{epoch}.png",
                                       dpi=dpi))
                fig = plot_radial_profile(stacks, objdesig=objdesig, epoch=epoch)
                written.append(savefig(fig, out_dir / f"radialprofile_{objdesig}_epoch{epoch}.png",
                                       dpi=dpi))

    if metrics is not None and not metrics.empty:
        fig = plot_contamination_analysis(metrics, objdesig=objdesig)
        written.append(savefig(fig, out_dir / f"contamination_{objdesig}.png", dpi=dpi))

    fig = plot_sourceflag_summary(phot, objdesig=objdesig)
    written.append(savefig(fig, out_dir / f"sourceflags_{objdesig}.png", dpi=dpi))
    fig = plot_aperture_validity(phot, cfg, objdesig=objdesig)
    written.append(savefig(fig, out_dir / f"apertures_{objdesig}.png", dpi=dpi))

    log.info("wrote %d target figures to %s", len(written), out_dir)
    return written


def plot_flag_survey(survey: pd.DataFrame, threshold: float = 3.0):
    """
    Flag effectiveness across every target in the survey.

    Parameters
    ----------
    survey : pandas.DataFrame
        Output of :func:`diagnostics.survey_flag_effectiveness`.
    threshold : float
        The ``sb_rise_sigma`` cut used, for the title.

    Returns
    -------
    matplotlib.figure.Figure

    Notes
    -----
    The decisive panel is lift against retained fraction: a cut worth applying
    sits high (it finds contamination) *and* to the right (it keeps data).
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    ax = axes[0]
    ax.hist(survey["base_rate"], bins=20, color="tab:purple", alpha=0.8)
    ax.axvline(survey["base_rate"].median(), color="k", ls="--", lw=2,
               label=f"median {survey['base_rate'].median():.3f}")
    ax.set_xlabel("contaminated fraction of exposures")
    ax.set_ylabel("targets")
    ax.set_title("how contaminated is each target?")
    ax.legend(fontsize=11)

    ax = axes[1]
    for tag, color in (("a", "tab:red"), ("a+b", "tab:orange")):
        if f"lift_{tag}" not in survey:
            continue
        ax.scatter(survey[f"kept_{tag}"], survey[f"lift_{tag}"], s=30, alpha=0.75,
                   color=color, label=f"cut {tag}")
    ax.axhline(1.0, color="k", ls=":", lw=2, label="uninformative (lift = 1)")
    ax.set_xlabel("fraction of data kept")
    ax.set_ylabel("lift  (precision / base rate)")
    ax.set_yscale("log")
    ax.set_title("informative, and at what cost")
    ax.legend(fontsize=11)

    ax = axes[2]
    for tag, color in (("a", "tab:red"), ("a+b", "tab:orange")):
        if f"recall_{tag}" not in survey:
            continue
        ax.scatter(survey[f"recall_{tag}"], survey[f"lift_{tag}"], s=30, alpha=0.75,
                   color=color, label=f"cut {tag}")
    ax.axhline(1.0, color="k", ls=":", lw=2)
    ax.set_xlabel("recall (fraction of contamination removed)")
    ax.set_ylabel("lift")
    ax.set_yscale("log")
    ax.set_title("recall vs. precision trade")
    ax.legend(fontsize=11)

    fig.suptitle(f"source-flag effectiveness across {len(survey)} targets "
                 f"(contaminated = growth-curve rise > {threshold:g}$\\sigma$)", y=1.02)
    fig.tight_layout()
    return fig
