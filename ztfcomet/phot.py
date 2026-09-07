"""Aperture photometry and Af-rho for ZTF comet cutouts.

This replaces the pipeline that lived in ``notebooks/afrho_240P.ipynb``.  Four
corrections separate it from that version (C1, C4, C5, C8 of
``doc/primitive_code_analysis.md``):

**C1 — per-frame zeropoint.**  The notebook wrote
``phot_target["inst_mag"] + row.zpmag`` where ``row`` was the loop variable left
over from the previous cell, so one frame's ``MAGZP`` was broadcast over the
whole table.  Measured spread in the ``2024E1`` sample was 3.32 mag in ZTF_g —
a factor of 21 in flux — and it silently mixed g and r.  Every Af-rho value
produced before this module was wrong by up to that factor.

**C4 — colour term.**  ZTF calibration is
``mag = -2.5 log10(DN) + MAGZP + CLRCOEFF * colour``.  ``CLRCOEFF`` ranges over
[-0.20, +0.13] in the sample and was never applied.  Comets are redder than the
stellar calibrators, so this is a systematic that varies with epoch and filter,
distorting lightcurve *shape*, not just normalisation.  Colour is measured from
same-night g/r pairs where they exist and falls back to a configured default,
flagged in ``flag_color_default``.

**C5 — aperture correction.**  ``MAGZP`` is the zeropoint for *PSF* photometry
(the header says so).  ``APCOR1..6`` give ``PSF_mag - AP_mag`` at aperture
diameters 2, 3, 4, 6, 10 and 14 px; this module interpolates them in
log-diameter.  Note the caveat: aperture corrections are derived from point
sources, so applying them to an extended coma is approximate.  The correction is
small (|APCOR| < 0.01 mag at 14 px diameter) for typical apertures and can be
switched off with ``PhotConfig.apply_aperture_correction``.

**C8 — uncertainties.**  ``afrho_err`` propagates photon noise, sky noise,
``MAGZPRMS``, the colour-term uncertainty (including its covariance with the
zeropoint, ``ZPCLRCOV``) and optionally the phase-coefficient uncertainty.  The
notebook computed ``source_sum_err`` and then dropped it.

Quality problems are reported as ``flag_*`` columns and never cause a row to be
dropped, so a frame that fails a cut stays in the table with the reason attached.
:data:`CRITICAL_FLAGS` mark measurements that are untrustworthy and drive
``quality_ok``; :data:`ADVISORY_FLAGS` merely qualify a good measurement.
"""

from __future__ import annotations

import dataclasses
import logging
import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import sep
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.time import Time
from astropy.wcs import WCS
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from tqdm.auto import tqdm

from . import config as cfg
from . import horizons as hz
from . import query as qry

__all__ = [
    "build_frame_table", "attach_ephemerides", "measure_photometry",
    "calibrate", "compute_afrho", "run_photometry", "flag_contamination",
    "aperture_scale_ok", "run_multi_aperture",
    "FLAG_COLUMNS", "CRITICAL_FLAGS", "ADVISORY_FLAGS",
]

log = logging.getLogger(__name__)

#: Aperture *diameters* (px) tabulated by the ZTF ``APCOR1..6`` keywords.
_APCOR_DIAMETERS = np.array([2.0, 3.0, 4.0, 6.0, 10.0, 14.0])
_APCOR_KEYS = [f"APCOR{i}" for i in range(1, 7)]

#: Flags that make the *measurement itself* untrustworthy.  ``quality_ok`` is
#: the negation of these, so they are what excludes a point from a default plot.
CRITICAL_FLAGS = [
    "flag_outside",        # target falls outside the cutout array
    "flag_aperture_edge",  # photometric aperture is clipped by the array edge
    "flag_sky_edge",       # sky annulus is clipped by the array edge
    "flag_undersampled",   # aperture smaller than min_rho_fwhm seeing FWHM
    "flag_centroid",       # winpos moved too far from the ephemeris position
    "flag_negative_flux",  # background-subtracted sum <= 0 (mag undefined)
    "flag_lowsnr",         # SNR < 3
    "flag_contaminated",   # catalogued background source(s) inside the aperture
]

#: Flags that qualify a measurement without invalidating it.  Both add a small
#: systematic or inflate the error, so they belong in the table and in any
#: uncertainty budget — but excluding a whole single-band night because its
#: colour had to be assumed (typically 0.015 mag) would throw away good data.
ADVISORY_FLAGS = [
    "flag_color_default",  # colour term used the assumed, not a measured, colour
    "flag_nan_pixels",     # negative pixels were clipped in the variance model
]

#: Every flag, critical and advisory.
FLAG_COLUMNS = CRITICAL_FLAGS + ADVISORY_FLAGS

_HEADER_KEYS = {
    "exptime": "EXPTIME", "egain": "GAIN", "readnoise": "READNOI",
    "pixscale": "PIXSCALE", "fwhm_pix": "SEEING", "filter": "FILTER",
    "zpmag": "MAGZP", "zpmagunc": "MAGZPUNC", "zpmagrms": "MAGZPRMS",
    "clrcoeff": "CLRCOEFF", "clrcounc": "CLRCOUNC", "zpclrcov": "ZPCLRCOV",
    "clrmed": "CLRMED", "maglim": "MAGLIM", "obsjd": "OBSJD", "origin": "ORIGIN",
}


def build_frame_table(datadir, pattern="*.fits", progress=True):
    """Inventory the FITS cutouts in *datadir* and pull their header metadata.

    Replaces ``ccdproc.ImageFileCollection``, which silently skipped unreadable
    files.  Files that cannot be opened are reported and counted here.

    Parameters
    ----------
    datadir : path-like
        Directory of ZTF cutouts.
    pattern : str
        Glob for the FITS files.

    Returns
    -------
    pandas.DataFrame
        One row per readable frame, sorted by ``obsjd``, with the header
        keywords in :data:`_HEADER_KEYS` plus ``naxis1``/``naxis2``.
    """
    datadir = Path(datadir)
    paths = sorted(datadir.glob(pattern))
    if not paths:
        log.warning("No files matching %r in %s", pattern, datadir)
        return pd.DataFrame()

    rows, unreadable = [], []
    iterator = tqdm(paths, desc="Reading headers") if progress else paths
    for path in iterator:
        try:
            header = fits.getheader(path)
            if "SIMPLE" not in header and "XTENSION" not in header:
                raise ValueError("no FITS signature")
        except Exception as exc:                                # noqa: BLE001
            unreadable.append((path.name, str(exc)))
            continue
        row = {"file": path.name}
        row.update({col: header.get(key) for col, key in _HEADER_KEYS.items()})
        row["naxis1"], row["naxis2"] = header.get("NAXIS1"), header.get("NAXIS2")
        for i, key in enumerate(_APCOR_KEYS):
            row[f"apcor{i + 1}"] = header.get(key)
        rows.append(row)

    if unreadable:
        log.warning("%d/%d files were not readable FITS (e.g. %s) — rerun the "
                    "download with repair=True", len(unreadable), len(paths),
                    unreadable[0][0])

    table = pd.DataFrame(rows)
    if not table.empty:
        table["isot"] = Time(table["obsjd"].to_numpy(), format="jd").isot
        table = table.sort_values("obsjd").reset_index(drop=True)
    return table


def attach_ephemerides(table, target, progress=True):
    """Add exact-time Horizons ephemerides to a frame table.

    Uses :meth:`~ztfcomet.config.Target.resolve_orbit_record`, so comets whose
    orbit solution changes between apparitions (240P/NEAT) get the right record
    for each epoch without any per-notebook special-casing.

    Returns
    -------
    pandas.DataFrame
        *table* with ``ra``, ``dec``, ``r``, ``delta``, ``alpha``, ``elong``,
        ``sunTargetPA``, ``velocityPA`` and ``orb_id`` columns added.
    """
    if table.empty:
        return table

    out = table.copy()
    out["orb_id"] = [target.resolve_orbit_record(jd) for jd in out["obsjd"]]

    # Horizons column -> our column.  Requested quantities do not all come back
    # for every object, so missing ones are filled with NaN rather than raising.
    column_map = {
        "RA": "ra", "DEC": "dec", "r": "r", "delta": "delta", "alpha": "alpha",
        "elong": "elong", "sunTargetPA": "sunTargetPA", "velocityPA": "velocityPA",
        # Predicted total magnitude, the reference brightness the Gaia
        # contamination check compares against.
        "Tmag": "tmag", "Nmag": "nmag",
        # Heliocentric range rate: negative inbound, positive outbound. Used to
        # separate the pre- and post-perihelion legs of an Af-rho vs r_h curve,
        # which otherwise fold on top of each other.
        "r_rate": "r_rate",
    }
    for col in column_map.values():
        out[col] = np.nan

    # Group by orbit record so each apparition is one batched Horizons call.
    groups = out.groupby("orb_id").groups
    iterator = (tqdm(groups.items(), total=len(groups), desc="Ephemerides")
                if progress else groups.items())
    for orb_id, idx in iterator:
        eph = qry.query_sso_ephemeris(
            orb_id, epochs=list(out.loc[idx, "obsjd"]),
            quantities=cfg.EPHEM_QUANTITIES_FULL,
            max_epochs_per_call=target.query.max_epochs_per_call,
            max_retries=target.query.max_retries, backoff=target.query.backoff)
        if eph.empty:
            log.warning("No ephemeris for orbit record %s (%d epochs)", orb_id, len(idx))
            continue

        # Confirm Horizons returned the object we asked for.  A record number
        # can be stale, and an ambiguous designation can resolve to a fragment
        # (240P -> 240P-B); "targetname" is the only place that shows up.
        if "targetname" in eph.columns and len(eph):
            returned = str(eph["targetname"].iloc[0])
            ok, why = hz.verify_targetname(returned, target.query_designation,
                                           allow_fragment=target.allow_fragment)
            out.loc[idx, "targetname"] = returned
            if not ok:
                log.error("%s: %s", target.name, why)
                raise ValueError(f"{target.name}: {why}")
            log.info("%s: Horizons record %s -> %s", target.name, orb_id, returned)

        available = {src: dst for src, dst in column_map.items() if src in eph.columns}
        missing = set(column_map) - set(available)
        if missing:
            log.warning("Horizons did not return %s for record %s",
                        sorted(missing), orb_id)

        # Always join on JD.  Even when the row counts match, matching by
        # position would rely on the service returning epochs in the order they
        # were sent -- the assumption that made the predecessor fragile (C10).
        merged = pd.merge_asof(
            out.loc[idx, ["obsjd"]].sort_values("obsjd"),
            eph.sort_values("datetime_jd")[["datetime_jd", *available]],
            left_on="obsjd", right_on="datetime_jd",
            direction="nearest", tolerance=1e-5)
        merged.index = out.loc[idx].sort_values("obsjd").index

        unmatched = int(merged["datetime_jd"].isna().sum())
        if unmatched:
            log.warning("%d/%d frames for record %s got no ephemeris within 0.9 s",
                        unmatched, len(merged), orb_id)

        for src, dst in available.items():
            out.loc[merged.index, dst] = merged[src].to_numpy()

    return out


def _variance_map(data, gain, readnoise):
    """Poisson + read-noise variance in DN, with negatives clipped.

    ZTF science images are not background-subtracted, so pixels are usually
    positive — but not always.  The notebook's ``np.sqrt(data/gain + ...)``
    produced NaN across 97.5% of one frame in the sample and propagated it
    silently into the SNR.  Clipping at zero keeps the variance finite; the
    number of clipped pixels inside the aperture is reported so the frame can be
    flagged rather than quietly trusted.
    """
    negative = data < 0
    poisson = np.clip(data, 0.0, None) / gain
    return poisson + (readnoise / gain) ** 2, negative


def _apcor_for_radius(row, radius_pix):
    """Interpolate ``APCOR1..6`` (PSF_mag - AP_mag) to an aperture *radius*.

    Interpolation is linear in log-diameter, where the curve is smoothest.
    Beyond the largest tabulated diameter (14 px) the correction is clamped to
    ``APCOR6``, which is already within 0.01 mag of zero.

    Returns ``0.0`` when the keywords are absent.
    """
    values = np.array([row.get(f"apcor{i}", np.nan) for i in range(1, 7)], dtype=float)
    good = np.isfinite(values)
    if good.sum() < 2:
        return 0.0
    diameter = np.clip(2.0 * radius_pix, _APCOR_DIAMETERS[good][0], _APCOR_DIAMETERS[good][-1])
    return float(np.interp(np.log(diameter),
                           np.log(_APCOR_DIAMETERS[good]), values[good]))


def aperture_scale_ok(rho_km, delta_au, pixscale, fwhm_pix, phot_config=None):
    """Is an aperture of *rho_km* usable on this frame?

    Requires ``FWHM < r_ap < max_aperture_arcsec`` on the sky.  Below the seeing
    FWHM the aperture does not contain the PSF, so it is not measuring the coma
    it claims to; above an arcminute it is dominated by sky and, in a 5 arcmin
    cutout, leaves no room for the sky annulus.

    Parameters
    ----------
    rho_km : float
        Aperture radius at the comet, km.
    delta_au : float
        Observer distance, au.
    pixscale : float
        Arcsec per pixel.
    fwhm_pix : float
        Seeing FWHM in pixels.

    Returns
    -------
    (ok, r_arcsec, fwhm_arcsec) : tuple of bool and float
    """
    pc = phot_config or cfg.PhotConfig()
    if not all(np.isfinite([rho_km, delta_au, pixscale, fwhm_pix])):
        return False, np.nan, np.nan
    delta_km = float(delta_au) * (1 * u.au).to_value(u.km)
    if delta_km <= 0:
        return False, np.nan, np.nan
    r_arcsec = float(rho_km) / delta_km * (1 * u.rad).to_value(u.arcsec)
    fwhm_arcsec = float(fwhm_pix) * float(pixscale)
    return bool(fwhm_arcsec < r_arcsec < pc.max_aperture_arcsec), r_arcsec, fwhm_arcsec


def measure_photometry(table, datadir, phot_config=None, progress=True):
    """Aperture photometry on every frame, with quality flags.

    For each frame: project the ephemeris position through the WCS, refine it
    with ``sep.winpos``, sum a circular aperture of physical radius
    ``rho_km`` at the comet, and subtract a sigma-clipped median sky measured in
    a surrounding annulus.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`attach_ephemerides`.
    datadir : path-like
        Directory holding the FITS named in ``table["file"]``.
    phot_config : ztfcomet.config.PhotConfig, optional

    Returns
    -------
    pandas.DataFrame
        *table* with geometry, instrumental photometry and :data:`FLAG_COLUMNS`.

    Notes
    -----
    Nothing is dropped.  A frame whose aperture falls off the array, or whose
    aperture is narrower than the seeing disc, keeps its measurement and gains a
    flag explaining why it should not be trusted.
    """
    pc = phot_config or cfg.PhotConfig()
    if table.empty:
        return table

    out = table.copy()
    datadir = Path(datadir)

    # Aperture geometry.  Small-angle: 1 arcsec at delta AU subtends
    # (1 arcsec in rad) * delta_km at the comet.
    km_per_pixel = (out["pixscale"] * (1 * u.arcsec).to_value(u.rad)
                    * out["delta"] * (1 * u.au).to_value(u.km))
    out["rho_km"] = pc.rho_km
    out["km_per_pixel"] = km_per_pixel
    out["rho_pix"] = out["rho_km"] / km_per_pixel
    out["sky_in_pix"] = pc.sky_in_scale * out["rho_pix"]
    out["sky_out_pix"] = pc.sky_out_scale * out["rho_pix"] + pc.sky_out_pad
    out["rho_fwhm"] = out["rho_pix"] / out["fwhm_pix"]

    for col in FLAG_COLUMNS:
        out[col] = False
    out["flag_undersampled"] = out["rho_fwhm"] < pc.min_rho_fwhm

    n_winpos_failed = 0
    last_winpos_error = None

    rows = out.iterrows()
    if progress:
        rows = tqdm(rows, total=len(out), desc="Aperture photometry")

    for idx, row in rows:
        path = datadir / row["file"]
        try:
            with fits.open(path) as hdul:
                data = hdul[0].data.astype(np.float64)
                header = hdul[0].header
                wcs = WCS(header)
        except Exception as exc:                                # noqa: BLE001
            log.warning("Cannot read %s: %s", row["file"], exc)
            continue

        variance, negative = _variance_map(data, row["egain"], row["readnoise"])
        error = np.sqrt(variance)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            x0, y0 = wcs.world_to_pixel(
                SkyCoord(ra=row["ra"] * u.deg, dec=row["dec"] * u.deg, frame="icrs"))
        x0, y0 = float(x0), float(y0)

        ny, nx = data.shape
        inside_array = (0 <= x0 < nx) and (0 <= y0 < ny)
        out.at[idx, "flag_outside"] = not inside_array

        if inside_array:
            # sep needs a C-contiguous native-endian float32 view, and winpos
            # returns THREE values (x, y, flag) -- unpacking only two raises,
            # and silently falling back would disable centroiding everywhere.
            try:
                xw, yw, winflag = sep.winpos(
                    np.ascontiguousarray(data, dtype=np.float32),
                    xinit=x0, yinit=y0,
                    sig=pc.winpos_sig_scale * row["fwhm_pix"])
                xw, yw = float(xw), float(yw)
                out.at[idx, "winpos_flag"] = int(np.atleast_1d(winflag)[0])
            except Exception as exc:                            # noqa: BLE001
                n_winpos_failed += 1
                last_winpos_error = exc
                xw, yw = x0, y0
        else:
            xw, yw = x0, y0

        shift = float(np.hypot(xw - x0, yw - y0))
        out.at[idx, "x_center"], out.at[idx, "y_center"] = xw, yw
        out.at[idx, "x_ephem"], out.at[idx, "y_ephem"] = x0, y0
        out.at[idx, "centroid_shift_pix"] = shift
        out.at[idx, "flag_centroid"] = shift > pc.max_centroid_shift_fwhm * row["fwhm_pix"]

        edge_dist = float(min(xw, yw, nx - 1 - xw, ny - 1 - yw))
        out.at[idx, "edge_dist_pix"] = edge_dist
        out.at[idx, "flag_aperture_edge"] = edge_dist < row["rho_pix"]
        out.at[idx, "flag_sky_edge"] = edge_dist < row["sky_out_pix"]

        aperture = CircularAperture((xw, yw), r=row["rho_pix"])
        annulus = CircularAnnulus((xw, yw), r_in=row["sky_in_pix"], r_out=row["sky_out_pix"])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            phot = aperture_photometry(data, aperture, error=error)

        # Count clipped (negative) pixels inside the aperture so the variance
        # approximation is visible rather than assumed harmless.
        ap_mask = aperture.to_mask(method="center")
        ap_neg = ap_mask.multiply(negative.astype(float))
        out.at[idx, "n_negative_pix"] = float(np.nansum(ap_neg)) if ap_neg is not None else np.nan
        out.at[idx, "flag_nan_pixels"] = bool(out.at[idx, "n_negative_pix"])

        sky_mask = annulus.to_mask(method="center")
        sky_cut = sky_mask.multiply(data)
        if sky_cut is None:
            msky = ssky = np.nan
            nsky = 0
        else:
            sky_pixels = sky_cut[sky_mask.data == 1]
            clipped = sigma_clip(sky_pixels, sigma=pc.sigma_clip_sigma,
                                 maxiters=pc.sigma_clip_iters, masked=True)
            msky = float(np.ma.median(clipped))
            ssky = float(np.ma.std(clipped))
            nsky = int(clipped.count())

        area = float(aperture.area)
        source_sum = float(phot["aperture_sum"][0]) - msky * area
        # Sky term: aperture photon noise plus the uncertainty on the sky level
        # itself, scaled by the aperture area.
        sky_var = (area * ssky ** 2) * (1.0 + area / nsky) if nsky else np.nan
        source_err = float(np.sqrt(float(phot["aperture_sum_err"][0]) ** 2 + sky_var))

        out.at[idx, "aparea"] = area
        out.at[idx, "msky"] = msky
        out.at[idx, "ssky"] = ssky
        out.at[idx, "nsky"] = nsky
        out.at[idx, "source_sum"] = source_sum
        out.at[idx, "source_sum_err"] = source_err
        out.at[idx, "snr"] = source_sum / source_err if source_err > 0 else np.nan

    if n_winpos_failed:
        # A systematic failure means every aperture sits on the raw ephemeris
        # position. That must never pass unnoticed.
        level = log.error if n_winpos_failed > 0.5 * len(out) else log.warning
        level("sep.winpos failed on %d/%d frames (%s); those apertures use the "
              "unrefined ephemeris position", n_winpos_failed, len(out), last_winpos_error)

    out["flag_negative_flux"] = ~(out["source_sum"] > 0)
    out["flag_lowsnr"] = ~(out["snr"] >= 3)

    with np.errstate(divide="ignore", invalid="ignore"):
        out["inst_mag"] = -2.5 * np.log10(out["source_sum"].where(out["source_sum"] > 0))
        out["inst_mag_err"] = 2.5 / np.log(10) / out["snr"]
    return out


def flag_contamination(table, phot_config=None, catalogue=None, progress=True):
    """Flag apertures containing catalogued background sources.

    A comet drifts across the star field, so on some frames a background star
    falls inside the photometric aperture.  The extra flux is indistinguishable
    from cometary activity in the image alone, and reads out as a spurious
    Af-rho spike.

    For each frame this sums the Gaia DR3 flux inside the aperture into a single
    effective magnitude

    .. math:: G_\\mathrm{eff} = -2.5\\log_{10}\\sum_i 10^{-0.4 G_i}

    and compares it with the comet's predicted total magnitude ``Tmag`` from the
    JPL ephemeris.  The frame is flagged when the contaminating flux reaches
    ``PhotConfig.contam_flux_ratio`` of the comet's — 30% by default, i.e.
    ``G_eff <= Tmag + 1.31``.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`measure_photometry`; needs ``ra``/``dec``, ``rho_pix``,
        ``fwhm_pix``, ``pixscale`` and ``tmag``.
    phot_config : ztfcomet.config.PhotConfig, optional
    catalogue : ztfcomet.gaia.GaiaCatalog, optional
        Defaults to the catalogue at :data:`ztfcomet.directory.GAIA_ROOT`.

    Returns
    -------
    pandas.DataFrame
        *table* with ``contam_radius_arcsec``, ``contam_n_sources``,
        ``contam_g_eff``, ``contam_g_brightest``, ``contam_ratio`` and
        ``flag_contaminated``.

    Notes
    -----
    The search radius is the aperture **plus** ``contam_radius_pad_fwhm`` seeing
    FWHM, because a star just outside the aperture still spills flux into it
    through the PSF wings.

    Two limits are worth stating with any result.  The catalogue is complete
    only to ``G = 18.5``, so "uncontaminated" means "no catalogued source" —
    fainter stars, and galaxies, are invisible to this test.  And Gaia ``G`` is
    compared directly with a visual ``Tmag``; the passbands differ by of order
    0.1-0.2 mag for typical stellar colours, which is small against a 30%
    (0.28 mag) threshold but not zero.
    """
    from . import gaia as gaia_module

    pc = phot_config or cfg.PhotConfig()
    if table.empty:
        return table

    out = table.copy()
    for col, value in [("contam_radius_arcsec", np.nan), ("contam_n_sources", 0),
                       ("contam_g_eff", np.inf), ("contam_g_brightest", np.inf),
                       ("contam_ratio", 0.0)]:
        out[col] = value
    if "flag_contaminated" not in out:
        out["flag_contaminated"] = False

    if not pc.check_contamination:
        log.info("Contamination check disabled by configuration")
        return out

    catalogue = catalogue or gaia_module.GaiaCatalog()
    if not catalogue.available:
        log.warning("%s — every frame left unflagged for contamination",
                    catalogue.describe())
        return out
    log.info(catalogue.describe())

    if "tmag" not in out or out["tmag"].isna().all():
        log.warning("No predicted magnitude (Tmag) available; "
                    "contamination cannot be assessed")
        return out

    rows = out.iterrows()
    if progress:
        rows = tqdm(rows, total=len(out), desc="Gaia contamination")

    for idx, row in rows:
        if not (np.isfinite(row.get("ra", np.nan)) and np.isfinite(row.get("dec", np.nan))):
            continue
        rho_pix = row.get("rho_pix", np.nan)
        fwhm_pix = row.get("fwhm_pix", np.nan)
        pixscale = row.get("pixscale", np.nan)
        if not np.isfinite(rho_pix) or not np.isfinite(pixscale):
            continue
        if not np.isfinite(fwhm_pix):
            fwhm_pix = 0.0

        radius = (rho_pix + pc.contam_radius_pad_fwhm * fwhm_pix) * pixscale
        result = catalogue.check_aperture(
            ra=row["ra"], dec=row["dec"], radius_arcsec=radius,
            comet_mag=row.get("tmag", np.nan),
            flux_ratio_threshold=pc.contam_flux_ratio,
        )

        out.at[idx, "contam_radius_arcsec"] = result.radius_arcsec
        out.at[idx, "contam_n_sources"] = result.n_sources
        out.at[idx, "contam_g_eff"] = result.g_eff
        out.at[idx, "contam_g_brightest"] = result.g_brightest
        out.at[idx, "contam_ratio"] = result.ratio
        out.at[idx, "flag_contaminated"] = result.contaminated

    n = int(out["flag_contaminated"].sum())
    if n:
        log.info("%d/%d frames flagged for background contamination "
                 "(>= %.0f%% of the comet flux)", n, len(out), 100 * pc.contam_flux_ratio)
    return out


def _measure_night_colors(table):
    """Median ``g - r`` per night from provisional (colour-free) magnitudes.

    Colour calibration is circular — the colour term needs a colour, which needs
    calibrated magnitudes.  Resolving it in two passes (calibrate without the
    term, measure the colour, then re-calibrate) converges immediately because
    the term is small.

    Returns a Series of ``g - r`` indexed by integer night, or an empty Series.
    """
    provisional = table["inst_mag"] + table["zpmag"]
    night = np.floor(table["obsjd"] - 0.5).astype("Int64")
    frame = pd.DataFrame({"night": night, "filter": table["filter"], "mag": provisional})
    frame = frame.dropna(subset=["mag", "night"])
    if frame.empty:
        return pd.Series(dtype=float)

    per_night = frame.groupby(["night", "filter"])["mag"].median().unstack()
    if "ZTF_g" not in per_night or "ZTF_r" not in per_night:
        return pd.Series(dtype=float)
    return (per_night["ZTF_g"] - per_night["ZTF_r"]).dropna()


def calibrate(table, phot_config=None):
    """Convert instrumental magnitudes to calibrated PS1 AB magnitudes.

    Applies, in order: the **per-frame** zeropoint ``MAGZP``, the colour term
    ``CLRCOEFF * colour``, and the interpolated aperture correction.  Propagates
    ``MAGZPRMS``, ``CLRCOUNC`` and the zeropoint-colour covariance ``ZPCLRCOV``
    into ``mag_err``.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`measure_photometry`.
    phot_config : ztfcomet.config.PhotConfig, optional

    Returns
    -------
    pandas.DataFrame
        *table* with ``color_gr``, ``apcor``, ``filter_mag``, ``filter_mag_err``
        and the ``flag_color_default`` flag resolved.

    Notes
    -----
    ``colour`` is ``g - r`` for the g and r bands and ``r - i`` for i, following
    the ZTF calibration convention.  Where a night has both g and r frames the
    colour is measured from their median magnitudes; otherwise
    ``PhotConfig.default_color_gr`` is used and the row is flagged.
    """
    pc = phot_config or cfg.PhotConfig()
    if table.empty:
        return table

    out = table.copy()

    night = np.floor(out["obsjd"] - 0.5).astype("Int64")
    measured = _measure_night_colors(out) if pc.apply_color_term else pd.Series(dtype=float)
    colors = night.map(measured) if len(measured) else pd.Series(np.nan, index=out.index)
    out["flag_color_default"] = colors.isna() & pc.apply_color_term
    out["color_gr"] = colors.fillna(pc.default_color_gr)
    # Colour uncertainty: the scatter of measured colours where we have them,
    # otherwise a generous 0.2 mag standing in for the unknown comet colour.
    color_err = float(np.nanstd(measured)) if len(measured) > 1 else 0.2
    out["color_gr_err"] = np.where(out["flag_color_default"], 0.2, max(color_err, 0.02))

    if pc.apply_color_term:
        color_term = out["clrcoeff"].fillna(0.0) * out["color_gr"]
    else:
        color_term = pd.Series(0.0, index=out.index)
    out["color_term"] = color_term

    if pc.apply_aperture_correction:
        out["apcor"] = [
            _apcor_for_radius(row, row["rho_pix"]) if np.isfinite(row.get("rho_pix", np.nan)) else 0.0
            for _, row in out.iterrows()
        ]
    else:
        out["apcor"] = 0.0

    out["filter_mag"] = out["inst_mag"] + out["zpmag"] + out["color_term"] + out["apcor"]

    # var = photon + ZP scatter + (colour * d_clrcoeff)^2 + (clrcoeff * d_colour)^2
    #       + 2 * colour * cov(ZP, clrcoeff)
    variance = out["inst_mag_err"].fillna(0.0) ** 2
    variance = variance + out["zpmagrms"].fillna(0.0) ** 2
    if pc.apply_color_term:
        variance = (variance
                    + (out["color_gr"] * out["clrcounc"].fillna(0.0)) ** 2
                    + (out["clrcoeff"].fillna(0.0) * out["color_gr_err"]) ** 2
                    + 2.0 * out["color_gr"] * out["zpclrcov"].fillna(0.0))
    out["filter_mag_err"] = np.sqrt(variance.clip(lower=0.0))
    return out


def compute_afrho(table, phot_config=None):
    """Compute Af-rho and its phase-corrected form A(0 deg) f rho.

    Uses the A'Hearn et al. (1984) definition

    .. math::  Af\\rho = \\frac{4 \\Delta^2 r_h^2}{\\rho} \\,
               10^{-0.4 (m_\\mathrm{comet} - m_\\odot)}

    with :math:`\\Delta` and :math:`\\rho` as lengths and :math:`r_h` in AU as a
    bare number, then reduces to zero phase angle with a linear phase law,
    :math:`\\Phi = 10^{0.4 \\beta \\alpha}`.

    Parameters
    ----------
    table : pandas.DataFrame
        Output of :func:`calibrate`.
    phot_config : ztfcomet.config.PhotConfig, optional

    Returns
    -------
    pandas.DataFrame
        *table* with ``solmag``, ``flux_ratio_sol``, ``afrho_cm``,
        ``afrho_cm_err``, ``phase_corr``, ``afrho0_cm`` and ``afrho0_cm_err``.

    Notes
    -----
    Af-rho is aperture-dependent by construction, so ``rho_km`` must be quoted
    alongside any value from this function.  The fractional error follows from
    :math:`\\mathrm{d}(Af\\rho)/Af\\rho = 0.4 \\ln 10 \; \\sigma_m`.
    """
    pc = phot_config or cfg.PhotConfig()
    if table.empty:
        return table

    out = table.copy()
    out["solmag"] = out["filter"].map(cfg.SOLAR_APPMAG_AB)
    if out["solmag"].isna().any():
        log.warning("No solar magnitude for filters: %s",
                    sorted(out.loc[out["solmag"].isna(), "filter"].dropna().unique()))

    out["magdiff_sol"] = out["filter_mag"] - out["solmag"]
    out["flux_ratio_sol"] = 10 ** (-0.4 * out["magdiff_sol"])

    delta = out["delta"].to_numpy() * u.au
    rh = out["r"].to_numpy()                    # AU, dimensionless by definition
    rho = out["rho_km"].to_numpy() * u.km
    afrho = (4 * delta ** 2 * rh ** 2 / rho * out["flux_ratio_sol"].to_numpy()).to_value(u.cm)
    out["afrho_cm"] = afrho

    frac_err = 0.4 * np.log(10) * out["filter_mag_err"].to_numpy()
    out["afrho_cm_err"] = np.abs(afrho) * frac_err

    alpha = out["alpha"].to_numpy()
    out["phase_corr"] = 10 ** (0.4 * pc.phase_beta * alpha)
    out["afrho0_cm"] = out["afrho_cm"] * out["phase_corr"]

    # Adding the phase-coefficient uncertainty in quadrature; zero by default.
    phase_frac = 0.4 * np.log(10) * pc.phase_beta_err * alpha
    out["afrho0_cm_err"] = np.abs(out["afrho0_cm"]) * np.sqrt(frac_err ** 2 + phase_frac ** 2)

    out["quality_ok"] = ~out[CRITICAL_FLAGS].any(axis=1)
    out["flags"] = out[FLAG_COLUMNS].apply(
        lambda r: ",".join(c.replace("flag_", "") for c in FLAG_COLUMNS if r[c]), axis=1)
    return out


def run_multi_aperture(target, datadir=None, phot_config=None, rho_km_set=None,
                       progress=True, save=True):
    """Reduce one target at several aperture radii.

    Af-rho depends on the aperture, so a survey measures a set of radii and
    quotes rho with every value.  The expensive shared work -- reading headers
    and fetching exact-time ephemerides -- is done once; only the aperture-
    dependent steps repeat.

    An aperture is measured on a frame only where
    ``FWHM < r_ap < max_aperture_arcsec`` holds (see :func:`aperture_scale_ok`).
    Frames failing that test for a given radius produce **no row** for it, which
    is deliberate: the measurement would not mean what its label says.  The
    counts are returned so nothing is lost silently.

    Parameters
    ----------
    target : ztfcomet.config.Target
    datadir : path-like, optional
    phot_config : ztfcomet.config.PhotConfig, optional
    rho_km_set : sequence of float, optional
        Defaults to ``phot_config.rho_km_set``.

    Returns
    -------
    table : pandas.DataFrame
        Long format: one row per (frame, aperture), with an ``rho_km`` column.
    skipped : dict
        ``{rho_km: n_frames_skipped}`` from the scale test.
    """
    from . import directory as d

    pc = phot_config or target.phot
    rho_set = list(rho_km_set if rho_km_set is not None else pc.rho_km_set)
    datadir = Path(datadir) if datadir else d.data_dir(target.name)

    base = build_frame_table(datadir, progress=progress)
    if base.empty:
        log.warning("No frames for %s in %s", target.name, datadir)
        return base, {}
    base.insert(1, "target", target.name)
    base = attach_ephemerides(base, target, progress=progress)

    parts, skipped = [], {}
    for rho_km in rho_set:
        ok = np.array([
            aperture_scale_ok(rho_km, row.get("delta"), row.get("pixscale"),
                              row.get("fwhm_pix"), pc)[0]
            for _, row in base.iterrows()])
        skipped[rho_km] = int((~ok).sum())
        if not ok.any():
            log.info("%s: rho=%.0f km unusable on every frame (scale test)",
                     target.name, rho_km)
            continue

        sub = base[ok].reset_index(drop=True)
        cfg_rho = dataclasses.replace(pc, rho_km=float(rho_km))
        sub = measure_photometry(sub, datadir, cfg_rho, progress=progress)
        sub = flag_contamination(sub, cfg_rho, progress=progress)
        sub = calibrate(sub, cfg_rho)
        sub = compute_afrho(sub, cfg_rho)
        sub["rho_km"] = float(rho_km)
        parts.append(sub)
        log.info("%s: rho=%5.0f km -> %d frames (%d skipped by scale test)",
                 target.name, rho_km, len(sub), skipped[rho_km])

    if not parts:
        log.warning("%s: no aperture passed the scale test on any frame", target.name)
        return pd.DataFrame(), skipped

    table = pd.concat(parts, ignore_index=True).sort_values(["rho_km", "obsjd"])
    table = table.reset_index(drop=True)

    if save:
        outpath = d.result_dir(target.name) / f"photometry_{d.target_slug(target.name)}.csv"
        table.to_csv(outpath, index=False)
        log.info("Wrote %s (%d rows, %d apertures)", outpath, len(table),
                 table["rho_km"].nunique())
    return table, skipped


def run_photometry(target, datadir=None, phot_config=None, progress=True, save=True):
    """End-to-end photometry for one target: inventory, ephemerides, Af-rho.

    Parameters
    ----------
    target : ztfcomet.config.Target
    datadir : path-like, optional
        Defaults to :func:`ztfcomet.directory.data_dir` for the target.
    phot_config : ztfcomet.config.PhotConfig, optional
        Defaults to ``target.phot``.
    save : bool
        Write ``photometry_<target>.csv`` into the target's results directory.

    Returns
    -------
    pandas.DataFrame
        The full photometry table, flagged but not filtered.
    """
    from . import directory as d

    pc = phot_config or target.phot
    datadir = Path(datadir) if datadir else d.data_dir(target.name)

    table = build_frame_table(datadir, progress=progress)
    if table.empty:
        log.warning("No frames for %s in %s", target.name, datadir)
        return table

    table.insert(1, "target", target.name)
    table = attach_ephemerides(table, target, progress=progress)
    table = measure_photometry(table, datadir, pc, progress=progress)
    table = flag_contamination(table, pc, progress=progress)
    table = calibrate(table, pc)
    table = compute_afrho(table, pc)

    n_flagged = int((~table["quality_ok"]).sum())
    if n_flagged:
        counts = {c: int(table[c].sum()) for c in FLAG_COLUMNS if table[c].any()}
        log.info("%s: %d/%d frames flagged %s", target.name, n_flagged, len(table), counts)

    if save:
        outpath = d.result_dir(target.name) / f"photometry_{d.target_slug(target.name)}.csv"
        table.to_csv(outpath, index=False)
        log.info("Wrote %s", outpath)
    return table
