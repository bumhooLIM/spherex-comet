"""Figures: annotated cutout images and Af-rho lightcurves.

The annotated-cutout block was pasted four times across the old notebooks and
the copies had already drifted apart (one gained an ``Elong`` line, another lost
a blank line).  It is one function here, :func:`plot_cutout`.

All figures import :mod:`ztfcomet.rcparams` for the project style.  Batch jobs
should pass ``dpi=50`` — writing hundreds of 200-dpi PNGs is slow and large.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS
from photutils.aperture import CircularAnnulus, CircularAperture
from tqdm.auto import tqdm

from . import rcparams  # noqa: F401  — imported for its rcParams side effects
from . import config as cfg

__all__ = ["plot_cutout", "plot_cutout_grid", "plot_afrho", "plot_afrho_vs_rh",
           "plot_afrho_apertures", "signed_rh", "save_all_cutouts"]

log = logging.getLogger(__name__)

_ARROW_STYLE = dict(facecolor="white", edgecolor="white", width=3,
                    head_width=10, head_length=5)
_ARROW_TEXT = dict(color="white", ha="center", va="center", fontsize=15, weight="bold")


def _sky_axes(wcs, shape):
    """Unit pixel vectors pointing north and east, measured from the WCS.

    The predecessor hardcoded ``N = (0, -1)`` and ``E = (-1, 0)``.  That is only
    right for a perfectly aligned field: measured across the ZTF sample, frames
    carry up to ~5 deg of rotation, so the drawn compass — and every position-angle
    vector referred to it — was off by that much.  Deriving the axes from the WCS
    costs nothing and is correct for any orientation or instrument.

    Returns
    -------
    (north, east) : tuple of ndarray
        Unit vectors in pixel space, or the conventional fallback if the WCS is
        unusable.
    """
    ny, nx = shape
    try:
        centre = wcs.pixel_to_world(nx / 2, ny / 2)
        step = 30 * u.arcsec
        origin = np.array([nx / 2, ny / 2])

        north = np.array(wcs.world_to_pixel(
            SkyCoord(centre.ra, centre.dec + step))) - origin
        east = np.array(wcs.world_to_pixel(
            SkyCoord(centre.ra + step / np.cos(centre.dec.radian), centre.dec))) - origin

        north = north / np.linalg.norm(north)
        east = east / np.linalg.norm(east)
        if not (np.all(np.isfinite(north)) and np.all(np.isfinite(east))):
            raise ValueError("degenerate WCS")
    except Exception:                                           # noqa: BLE001
        log.debug("Could not derive sky axes from the WCS; using N down / E left")
        north, east = np.array([0.0, -1.0]), np.array([-1.0, 0.0])
    return north, east


def _draw_compass(ax, wcs, shape, sun_pa=None, vel_pa=None,
                  base=(0.2, 0.8), length_frac=0.13):
    """Draw the N/E axes plus the anti-solar and anti-velocity vectors.

    Position angles run east of north, so a vector at PA theta is
    ``N cos(theta) + E sin(theta)`` in the WCS-derived frame.  Horizons'
    ``sunTargetPA`` and ``velocityPA`` are already the *extended* Sun-target
    radius vector and the *negative* heliocentric velocity, which is why they are
    labelled -(sun) and -V.
    """
    ny, nx = shape
    origin = np.array([base[0] * nx, base[1] * ny])
    length = length_frac * min(nx, ny)
    north, east = _sky_axes(wcs, shape)

    for vector, label in [(north, "N"), (east, "E")]:
        delta = vector * length
        ax.arrow(origin[0], origin[1], delta[0], delta[1], **_ARROW_STYLE)
        ax.text(*(origin + delta * 1.4), label, **_ARROW_TEXT)

    for pa, label in [(vel_pa, r"$-V$"), (sun_pa, r"$-\odot$")]:
        if pa is None or not np.isfinite(pa):
            continue
        rad = np.deg2rad(pa)
        delta = (north * np.cos(rad) + east * np.sin(rad)) * length
        ax.arrow(origin[0], origin[1], delta[0], delta[1], **_ARROW_STYLE)
        ax.text(*(origin + delta * 1.4), label, **_ARROW_TEXT)


def plot_cutout(path, row=None, target=None, ax=None, show_apertures=True,
                phot_config=None, cmap="gray", annotate=True, title=None):
    """Plot one ZTF cutout with WCS axes, an ephemeris panel and a compass.

    Parameters
    ----------
    path : path-like
        FITS cutout to display.
    row : pandas.Series, optional
        Matching row of a photometry table.  Supplies the ephemeris values, the
        measured centroid and the aperture geometry; without it only header
        information is shown.
    target : ztfcomet.config.Target or str, optional
        Used for the ``OBJECT`` label.
    show_apertures : bool
        Overlay the photometric aperture and sky annulus when *row* has them.
    phot_config : ztfcomet.config.PhotConfig, optional
        Only needed when *row* lacks the aperture columns.

    Returns
    -------
    matplotlib.axes.Axes
    """
    path = Path(path)
    with fits.open(path) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(projection=wcs)

    vmin, vmax = ZScaleInterval().get_limits(data)
    ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin="lower")

    obstime = Time(header["OBSJD"], format="jd")
    name = getattr(target, "name", target) or header.get("OBJECT", "—")

    sun_pa = vel_pa = None
    if row is not None:
        sun_pa = row.get("sunTargetPA")
        vel_pa = row.get("velocityPA")
    _draw_compass(ax, wcs, data.shape, sun_pa=sun_pa, vel_pa=vel_pa)

    if show_apertures and row is not None and np.isfinite(row.get("x_center", np.nan)):
        pc = phot_config or cfg.PhotConfig()
        centre = (row["x_center"], row["y_center"])
        rho = row.get("rho_pix", np.nan)
        if np.isfinite(rho):
            CircularAperture(centre, r=rho).plot(ax=ax, color="cyan", lw=1.5)
            CircularAnnulus(centre,
                            r_in=row.get("sky_in_pix", pc.sky_in_scale * rho),
                            r_out=row.get("sky_out_pix", pc.sky_out_scale * rho + pc.sky_out_pad)
                            ).plot(ax=ax, color="yellow", lw=1.0, ls="--")

    if annotate:
        lines = [
            f"OBSERVAT = {header.get('ORIGIN', '—')}",
            f"DATE-OBS = {obstime.isot}",
            f"FILTER   = {header.get('FILTER', '—')}",
            f"EXPTIME  = {header.get('EXPTIME', '—')} sec",
            "",
            f"OBJECT   = {name}",
        ]
        if row is not None:
            for key, label, fmt in [("r", "Rh", "{:.3f} AU"), ("delta", "Delta", "{:.3f} AU"),
                                    ("alpha", "Phase", "{:.1f} deg"), ("elong", "Elong", "{:.1f} deg"),
                                    ("sunTargetPA", "S-T PA", "{:.1f} deg"),
                                    ("velocityPA", "Vel PA", "{:.1f} deg")]:
                value = row.get(key)
                if value is not None and np.isfinite(value):
                    lines.append(f"{label:8s} = {fmt.format(value)}")
            if np.isfinite(row.get("filter_mag", np.nan)):
                lines += ["", f"{'mag':8s} = {row['filter_mag']:.2f} +/- {row.get('filter_mag_err', np.nan):.2f}"]
            if np.isfinite(row.get("afrho0_cm", np.nan)):
                lines.append(f"{'A(0)frho':8s} = {row['afrho0_cm']:.1f} cm")
            flags = [c for c in row.index if c.startswith("flag_") and bool(row[c])]
            if flags:
                lines += ["", "FLAGS    = " + ", ".join(f.replace("flag_", "") for f in flags)]
        ax.annotate("\n".join(lines), xy=(1.05, 0.95), xycoords="axes fraction",
                    fontsize=13, color="k", va="top", ha="left", family="monospace")

    # Only claim sky coordinates when the axes actually carry the WCS
    # projection; a bare subplot shows pixels, and labelling those "RA"/"Dec"
    # misrepresents the figure.
    if hasattr(ax, "coords"):
        ax.set_xlabel("RA")
        ax.set_ylabel("Dec")
    else:
        ax.set_xlabel("x (pix)")
        ax.set_ylabel("y (pix)")
    if title:
        ax.set_title(title)
    return ax


def save_all_cutouts(table, datadir, figdir, target=None, phot_config=None,
                     dpi=50, progress=True):
    """Write one annotated PNG per frame.

    Uses ``dpi=50`` by default: the project style sets 200 dpi for one-off
    figures, which is wasteful for hundreds of batch files.

    Returns
    -------
    list of pathlib.Path
    """
    datadir, figdir = Path(datadir), Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    written = []
    rows = table.iterrows()
    if progress:
        rows = tqdm(rows, total=len(table), desc="Saving cutouts")

    for _, row in rows:
        path = datadir / row["file"]
        if not path.exists():
            continue
        try:
            ax = plot_cutout(path, row=row, target=target, phot_config=phot_config)
            outpath = figdir / f"{path.stem}.png"
            ax.figure.savefig(outpath, dpi=dpi)
            written.append(outpath)
        except Exception as exc:                                # noqa: BLE001
            log.warning("Could not plot %s: %s", row["file"], exc)
        finally:
            plt.close("all")

    log.info("Wrote %d cutout figures to %s", len(written), figdir)
    return written


def plot_cutout_grid(table, datadir, target=None, ncols=4, nmax=12,
                     phot_config=None, figsize_per=(3.4, 3.4)):
    """Contact sheet of up to *nmax* cutouts, for a quick look at a whole run."""
    datadir = Path(datadir)
    subset = table.head(nmax)
    n = len(subset)
    if n == 0:
        raise ValueError("nothing to plot")

    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    fig = plt.figure(figsize=(figsize_per[0] * ncols, figsize_per[1] * nrows))

    for i, (_, row) in enumerate(subset.iterrows(), start=1):
        path = datadir / row["file"]
        if not path.exists():
            continue
        with fits.open(path) as hdul:
            wcs = WCS(hdul[0].header)
        ax = fig.add_subplot(nrows, ncols, i, projection=wcs)
        plot_cutout(path, row=row, target=target, ax=ax, annotate=False,
                    phot_config=phot_config,
                    title=f"{row['isot'][:10]} {row['filter']}")
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.tight_layout()
    return fig


def signed_rh(table):
    """Heliocentric distance signed by orbital leg: negative inbound.

    An Af-rho vs r_h curve folds the pre- and post-perihelion legs on top of
    each other, hiding the hysteresis that is often the point of the plot.
    Signing r_h by the heliocentric range rate separates them on one axis:
    ``r_rate < 0`` (approaching, pre-perihelion) becomes negative, ``r_rate > 0``
    (receding, post-perihelion) stays positive, so the comet traces left to
    right through perihelion at x = 0.

    Falls back to unsigned r_h when ``r_rate`` is unavailable.
    """
    rh = pd.to_numeric(table["r"], errors="coerce")
    if "r_rate" not in table:
        return rh
    rate = pd.to_numeric(table["r_rate"], errors="coerce")
    sign = np.where(rate < 0, -1.0, 1.0)
    return rh * np.where(np.isfinite(rate), sign, 1.0)


def plot_afrho_vs_rh(tables, rho_km=None, filters=("ZTF_r",), only_good=True,
                     split_perihelion=True, ax=None, colors=None, markers=None,
                     legend=True, annotate_legs=True):
    """Af-rho against heliocentric distance, clean points only.

    Parameters
    ----------
    tables : DataFrame or mapping of label -> DataFrame
    rho_km : float, optional
        Select one aperture from a multi-aperture table.  Required when the
        table holds more than one, since mixing apertures on one axis is
        meaningless.
    only_good : bool
        Default **True** here: this figure is the science result, so flagged
        frames are excluded rather than drawn as open symbols.
    split_perihelion : bool
        Sign r_h by orbital leg (see :func:`signed_rh`) so the inbound and
        outbound branches do not overlap.  A marker at x = 0 is perihelion.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if not isinstance(tables, dict):
        tables = {"target": tables}

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    palette = colors or plt.rcParams["axes.prop_cycle"].by_key().get("color", ["k"])
    marker_cycle = markers or ["o", "s", "^", "D", "v", "P"]

    def _selected(table):
        """The rows this call will actually draw, for the pre-scan."""
        sub = table
        if rho_km is not None and "rho_km" in sub:
            sub = sub[np.isclose(sub["rho_km"], rho_km)]
        if only_good and "quality_ok" in sub:
            sub = sub[sub["quality_ok"]]
        return sub

    # Decide about signing BEFORE plotting. Signing exists to separate two
    # orbital legs; with only one leg it would render every r_h negative under
    # an axis labelled r_h. Most targets in a 18-month window are single-leg.
    saw_inbound = saw_outbound = False
    for _t in tables.values():
        if _t is None or len(_t) == 0 or "r_rate" not in _t:
            continue
        rate = pd.to_numeric(_selected(_t)["r_rate"], errors="coerce")
        saw_inbound |= bool((rate < 0).any())
        saw_outbound |= bool((rate >= 0).any())
    use_signed = bool(split_perihelion and saw_inbound and saw_outbound)

    for i, (label, table) in enumerate(tables.items()):
        if table is None or len(table) == 0:
            continue
        sub = table
        if rho_km is not None and "rho_km" in sub:
            sub = sub[np.isclose(sub["rho_km"], rho_km)]
        elif "rho_km" in sub and sub["rho_km"].nunique() > 1:
            raise ValueError("table holds several apertures; pass rho_km=")

        if only_good and "quality_ok" in sub:
            sub = sub[sub["quality_ok"]]
        sub = sub[np.isfinite(pd.to_numeric(sub.get("afrho0_cm"), errors="coerce"))]
        if sub.empty:
            continue

        x = signed_rh(sub) if use_signed else pd.to_numeric(sub["r"], errors="coerce")

        for j, band in enumerate(filters):
            keep = sub["filter"] == band
            if not keep.any():
                continue
            name = f"{label} ({band})" if len(tables) > 1 or len(filters) > 1 else band
            ax.errorbar(x[keep], sub.loc[keep, "afrho0_cm"],
                        yerr=sub.loc[keep, "afrho0_cm_err"],
                        fmt=marker_cycle[j % len(marker_cycle)], ms=6,
                        color=palette[(j if len(tables) == 1 else i) % len(palette)],
                        ecolor=palette[(j if len(tables) == 1 else i) % len(palette)],
                        elinewidth=1, capsize=2, ls="none", label=name)

    if use_signed:
        ax.axvline(0, color="0.5", ls=":", lw=1.5)
        if annotate_legs:
            ax.annotate("pre-perihelion", xy=(0.02, 0.02), xycoords="axes fraction",
                        fontsize=12, ha="left", va="bottom", color="0.35")
            ax.annotate("post-perihelion", xy=(0.98, 0.02), xycoords="axes fraction",
                        fontsize=12, ha="right", va="bottom", color="0.35")
        ax.set_xlabel(r"$-r_\mathrm{h}$  |  $+r_\mathrm{h}$  (AU)")
        # Show distance, not the sign, on the tick labels.
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _pos: f"{abs(v):g}"))
    else:
        # Single leg: plain r_h, labelled with which leg, so the figure is not
        # silently ambiguous about where the comet was in its orbit.
        leg = ("inbound, pre-perihelion" if saw_inbound and not saw_outbound
               else "outbound, post-perihelion" if saw_outbound and not saw_inbound
               else None)
        ax.set_xlabel(r"$r_\mathrm{h}$ (AU)" + (f"   [{leg}]" if leg else ""))

    ax.set_ylabel(r"$A(0\degree)f\rho$ (cm)")
    if legend:
        ax.legend(fontsize=11, frameon=False)
    return ax


def plot_afrho_apertures(table, filters=("ZTF_r",), only_good=True,
                         split_perihelion=True, ax=None, legend=True):
    """Af-rho vs r_h for every aperture in a multi-aperture table.

    One colour per ``rho_km``.  Because Af-rho for a steady-state ``1/rho`` coma
    is independent of aperture, the curves separating tells you the coma is not
    in steady state -- which is why plotting them together is worth doing.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))
    if table is None or len(table) == 0 or "rho_km" not in table:
        raise ValueError("need a table with an rho_km column")

    radii = sorted(table["rho_km"].dropna().unique())
    cmap = plt.get_cmap("viridis")
    for i, rho in enumerate(radii):
        colour = cmap(i / max(len(radii) - 1, 1))
        plot_afrho_vs_rh({f"{rho / 1000:g}k km": table}, rho_km=rho, filters=filters,
                         only_good=only_good, split_perihelion=split_perihelion,
                         ax=ax, colors=[colour], legend=False,
                         annotate_legs=(i == 0))
    if legend:
        ax.legend(fontsize=10, frameon=False, title=r"$\rho$")
    return ax


def plot_afrho(tables, filters=("ZTF_r",), x="jd", perihelion_jd=None,
               only_good=False, ax=None, colors=None, markers=None,
               show_flagged=True, legend=True):
    """Af-rho lightcurve for one or several targets.

    Parameters
    ----------
    tables : DataFrame or mapping of label -> DataFrame
        Photometry tables from :func:`ztfcomet.phot.run_photometry`.
    filters : sequence of str
        Which ZTF bands to draw.
    x : {"jd", "tp", "rh"}
        Abscissa: Julian date, ``T - T_p`` in days (needs *perihelion_jd*), or
        heliocentric distance in AU.
    perihelion_jd : float or mapping, optional
        Time of perihelion, per label when a mapping.
    only_good : bool
        Plot only rows with ``quality_ok``.  Default ``False``: flagged points
        are drawn as open symbols so nothing disappears without being seen.
    show_flagged : bool
        Draw flagged points (open symbols) when *only_good* is ``False``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if not isinstance(tables, dict):
        tables = {getattr(tables, "attrs", {}).get("target", "target"): tables}

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    palette = colors or plt.rcParams["axes.prop_cycle"].by_key().get("color", ["k"])
    marker_cycle = markers or ["o", "s", "^", "D", "v", "P"]

    for i, (label, table) in enumerate(tables.items()):
        if table is None or table.empty:
            continue
        for j, band in enumerate(filters):
            # With a single target, colour encodes the band; with several it
            # encodes the target and the marker encodes the band.
            color = (palette[j % len(palette)] if len(tables) == 1
                     else palette[i % len(palette)])
            sub = table[table["filter"] == band]
            if sub.empty:
                continue

            if x == "tp":
                tp = perihelion_jd.get(label) if isinstance(perihelion_jd, dict) else perihelion_jd
                if tp is None:
                    raise ValueError(f"x='tp' needs perihelion_jd for {label!r}")
                xv = sub["obsjd"] - tp
            elif x == "rh":
                xv = sub["r"]
            else:
                xv = sub["obsjd"]

            good = sub["quality_ok"] if "quality_ok" in sub else np.ones(len(sub), bool)
            marker = marker_cycle[j % len(marker_cycle)]
            name = f"{label} ({band})" if len(filters) > 1 or len(tables) > 1 else band

            ax.errorbar(xv[good], sub.loc[good, "afrho0_cm"],
                        yerr=sub.loc[good, "afrho0_cm_err"],
                        fmt=marker, ms=6, color=color, ecolor=color,
                        elinewidth=1, capsize=2, ls="none", label=name)

            if show_flagged and not only_good and (~good).any():
                ax.errorbar(xv[~good], sub.loc[~good, "afrho0_cm"],
                            yerr=sub.loc[~good, "afrho0_cm_err"],
                            fmt=marker, ms=6, mfc="none", color=color, ecolor=color,
                            elinewidth=1, capsize=2, ls="none", alpha=0.6,
                            label=f"{name} [flagged]")

    xlabels = {"tp": r"$T-T_\mathrm{p}$ (d)", "rh": r"$r_\mathrm{h}$ (AU)", "jd": "JD"}
    ax.set_xlabel(xlabels.get(x, x))
    ax.set_ylabel(r"$A(0\degree)f\rho$ (cm)")
    if legend:
        ax.legend(fontsize=11, frameon=False)
    return ax
