"""Radial surface-brightness profiles: is the coma active, and is it 1/rho?

Aperture photometry measures how much light is in the coma; it says nothing
about how that light is distributed.  A steady-state coma -- dust leaving the
nucleus at constant rate and speed -- has surface brightness falling as
``1/rho``, while a point source falls off with the PSF, far more steeply.
Comparing the comet's radial profile with the profile of field stars on the
*same frame* therefore answers two questions at once: is there extended
emission at all, and does it have the slope a simple coma should?

Procedure, per frame
--------------------
1. Comet profile in annuli at radii 0.5 .. 10 px (0.5 px step), each the
   sigma-clipped **mean** of the sky-subtracted pixels, so a field star that
   drifts into one annulus is rejected rather than averaged in.  The plain
   mean is kept alongside as the check on that rejection.
2. The comet cutout is **oversampled** first (bilinear, flux-conserving) so
   that 0.5 px annuli on 1"/px ZTF pixels contain enough samples to be
   well-defined.  Bilinear rather than cubic: cubic splines overshoot at a
   sharp core and can push an annulus mean negative.
3. Up to 20 field stars -- unsaturated, S/N > 10, isolated, point-like -- are
   profiled at native resolution, each normalised to its central surface
   brightness, and stacked with a sigma-clipped **median** at every radius.
4. A power law is fitted to each profile outside the PSF core.  ``m = -1`` is
   the steady-state coma; a star gives ``m`` of -3 or steeper.

Everything is written to ``results/`` (long-format CSV) and ``fig/``.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import sep
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

from . import config as cfg

__all__ = [
    "DEFAULT_RADII", "KM_PER_ARCSEC_AU", "radial_profile", "select_field_stars",
    "stack_star_profiles", "fit_powerlaw", "profile_frame", "run_profiles",
    "ProfileConfig", "physical_scales", "fit_coma_model", "psf_kernel_from_stack",
]

#: Kilometres subtended by one arcsecond at one au.
KM_PER_ARCSEC_AU = 1.495978707e8 * np.deg2rad(1.0 / 3600.0)


def physical_scales(frame, fit_rmin_pix=None, rmax_pix=10.0):
    """Pixel-to-kilometre scales for a frame (row with ``pixscale``, ``delta``, ``fwhm_pix``).

    Returns a dict: ``km_per_pix``, ``fwhm_km``, ``rho_fit_rmin_km`` (inner fit
    limit), ``rho_rmax_km`` (outer annulus), and ``n_pix_per_1e4km``.  These are
    what decide whether a 10 px profile is probing 4 000 km of coma or 15 000.
    """
    pixscale = float(frame.get("pixscale", np.nan))
    delta = float(frame.get("delta", np.nan))
    fwhm = float(frame.get("fwhm_pix", np.nan))
    km_per_pix = pixscale * delta * KM_PER_ARCSEC_AU
    rmin = float(fit_rmin_pix) if fit_rmin_pix is not None else np.nan
    return dict(km_per_pix=km_per_pix, fwhm_km=fwhm * km_per_pix,
                rho_fit_rmin_km=rmin * km_per_pix, rho_rmax_km=rmax_pix * km_per_pix,
                n_pix_per_1e4km=1.0e4 / km_per_pix if km_per_pix > 0 else np.nan)

log = logging.getLogger(__name__)

#: Annulus centres, px: 0.5, 1.0, ..., 10.0.  Half-width 0.25 makes them tile.
DEFAULT_RADII = np.round(np.arange(0.5, 10.0 + 1e-9, 0.5), 2)


class ProfileConfig:
    """Tunables for the radial-profile step.  Plain attributes, no magic."""

    def __init__(self, radii=None, half_width=0.25, oversample=4, sigma=3.0,
                 maxiters=5, max_stars=20, snr_min=10.0, saturation_fraction=0.8,
                 detect_sigma=5.0, exclude_comet_pix=30.0, min_pixels=3,
                 fit_rmin_fwhm=1.5, fit_rmax_pix=10.0, star_isolation_pix=None,
                 oversample_stars=True):
        self.radii = np.asarray(DEFAULT_RADII if radii is None else radii, float)
        self.half_width = float(half_width)
        self.oversample = int(oversample)
        self.sigma = float(sigma)
        self.maxiters = int(maxiters)
        self.max_stars = int(max_stars)
        self.snr_min = float(snr_min)
        self.saturation_fraction = float(saturation_fraction)
        self.detect_sigma = float(detect_sigma)
        self.exclude_comet_pix = float(exclude_comet_pix)
        self.min_pixels = int(min_pixels)
        self.fit_rmin_fwhm = float(fit_rmin_fwhm)
        self.fit_rmax_pix = float(fit_rmax_pix)
        # Profile the field stars on the same oversampled grid as the comet, so
        # the two are strictly like-for-like. Off, stars use native pixels.
        self.oversample_stars = bool(oversample_stars)
        # Another source inside this radius disqualifies a star; default is the
        # outer annulus edge plus a margin, so the profile is never blended.
        self.star_isolation_pix = (float(star_isolation_pix) if star_isolation_pix
                                   else float(self.radii.max() + self.half_width + 3.0))

    @property
    def rmax(self):
        return float(self.radii.max() + self.half_width)


# --------------------------------------------------------------------- profile
def _cutout(data, x, y, pad):
    """Square box around (x, y) with the slice used, clipped to the array."""
    ny, nx = data.shape
    x0, y0 = int(round(x)), int(round(y))
    ys = slice(max(y0 - pad, 0), min(y0 + pad + 1, ny))
    xs = slice(max(x0 - pad, 0), min(x0 + pad + 1, nx))
    return data[ys, xs], ys, xs


def radial_profile(data, x, y, radii=None, half_width=0.25, sky=0.0, oversample=1,
                   sigma=3.0, maxiters=5, min_pixels=3):
    """Azimuthally averaged surface brightness in annuli around (x, y).

    Parameters
    ----------
    data : ndarray
        Image, in DN.  NaN pixels are ignored.
    x, y : float
        Centre in pixel coordinates (0-based, pixel centres at integers).
    radii : array_like, optional
        Annulus centres in px; default :data:`DEFAULT_RADII`.
    half_width : float
        Annulus half-width in px.
    sky : float
        Sky level to subtract before averaging.
    oversample : int
        Bilinear, flux-conserving zoom factor applied to the cutout first.  1
        means native pixels.
    sigma, maxiters : float, int
        Sigma-clipping parameters for the robust mean.
    min_pixels : int
        Annuli with fewer samples than this are reported as NaN.

    Returns
    -------
    pandas.DataFrame
        One row per radius: ``r_pix``, ``mean_clip``, ``std_clip``, ``n_clip``,
        ``mean``, ``std``, ``n``.  Surface brightness is DN per *original*
        pixel regardless of *oversample* -- the zoom conserves the mean.

    Notes
    -----
    The oversampled grid uses ``grid_mode=True``, so sub-pixel ``j`` sits at
    original coordinate ``(j + 0.5)/oversample - 0.5``; pixel centres stay at
    integers and the radius of every sample is exact.
    """
    radii = np.asarray(DEFAULT_RADII if radii is None else radii, float)
    pad = int(np.ceil(radii.max() + half_width)) + 2
    box, ys, xs = _cutout(np.asarray(data, float), x, y, pad)
    if box.size == 0:
        return _empty_profile(radii)

    box = box - float(sky)
    valid = np.isfinite(box)
    box = np.where(valid, box, 0.0)

    if oversample > 1:
        box = ndimage.zoom(box, oversample, order=1, grid_mode=True, mode="nearest")
        valid = ndimage.zoom(valid.astype(float), oversample, order=1,
                             grid_mode=True, mode="nearest") > 0.99
        step = 1.0 / oversample
        yy, xx = np.mgrid[0:box.shape[0], 0:box.shape[1]]
        # sub-pixel index -> original pixel coordinate
        px = xs.start + (xx + 0.5) * step - 0.5
        py = ys.start + (yy + 0.5) * step - 0.5
    else:
        yy, xx = np.mgrid[ys, xs]
        px, py = xx.astype(float), yy.astype(float)

    rr = np.hypot(px - x, py - y)

    rows = []
    for r in radii:
        sel = valid & (np.abs(rr - r) < half_width)
        vals = box[sel]
        if vals.size < min_pixels:
            rows.append(dict(r_pix=r, mean_clip=np.nan, std_clip=np.nan, n_clip=int(vals.size),
                             mean=np.nan, std=np.nan, n=int(vals.size)))
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mean_c, _, std_c = sigma_clipped_stats(vals, sigma=sigma, maxiters=maxiters)
            n_c = int(np.sum(np.abs(vals - mean_c) <= sigma * max(std_c, 1e-12))) \
                if std_c > 0 else int(vals.size)
        rows.append(dict(r_pix=r, mean_clip=float(mean_c), std_clip=float(std_c), n_clip=n_c,
                         mean=float(vals.mean()), std=float(vals.std()), n=int(vals.size)))
    return pd.DataFrame(rows)


def _empty_profile(radii):
    return pd.DataFrame({"r_pix": radii, "mean_clip": np.nan, "std_clip": np.nan,
                         "n_clip": 0, "mean": np.nan, "std": np.nan, "n": 0})


def central_sb(profile, r_core=1.5, column="mean_clip"):
    """Central surface brightness: pixel-weighted mean of the bins with r <= r_core.

    1.5 px spans three annuli (~12 native pixels), which keeps the normalisation
    stable for stars profiled without oversampling, where the innermost 0.5 px
    annulus often holds fewer than three pixels and is reported as NaN.
    """
    core = profile[(profile["r_pix"] <= r_core) & np.isfinite(profile[column])]
    if core.empty or core["n"].sum() == 0:
        return np.nan
    return float(np.average(core[column], weights=core["n"]))


def fit_powerlaw(r, sb, rmin, rmax):
    """Slope of log(SB) vs log(r) over [rmin, rmax]; NaN with < 3 points.

    Returns ``(slope, intercept, n_points)``.  For a steady-state coma the
    slope is -1; a stellar PSF is much steeper.
    """
    r = np.asarray(r, float)
    sb = np.asarray(sb, float)
    keep = np.isfinite(r) & np.isfinite(sb) & (sb > 0) & (r >= rmin) & (r <= rmax)
    if keep.sum() < 3:
        return np.nan, np.nan, int(keep.sum())
    slope, intercept = np.polyfit(np.log10(r[keep]), np.log10(sb[keep]), 1)
    return float(slope), float(intercept), int(keep.sum())


# ----------------------------------------------------------------- field stars
def select_field_stars(data, gain, readnoise, saturate, exclude_xy=None,
                       profile_config=None):
    """Isolated, unsaturated, point-like stars with S/N > threshold.

    Detection runs on a sep background-subtracted image.  The cuts, in order:

    * S/N > ``snr_min`` in a 3 px aperture;
    * peak below ``saturation_fraction`` of ``SATURATE``;
    * sep ``flag == 0`` -- no blending, no truncation;
    * elongation ``a/b < 1.5`` and a size within a factor of the median
      candidate (rejects galaxies and cosmic rays);
    * farther than ``exclude_comet_pix`` from the comet and than ``rmax + 2``
      from any edge;
    * no other detection within ``star_isolation_pix``.

    Returns
    -------
    stars : pandas.DataFrame
        Up to ``max_stars`` rows, brightest first: ``x``, ``y``, ``flux``,
        ``snr``, ``peak``, ``a``, ``b``.
    background : sep.Background
        So the caller can profile on the same subtracted image.
    """
    pc = profile_config or ProfileConfig()
    img = np.ascontiguousarray(np.where(np.isfinite(data), data, 0.0), dtype=np.float32)
    bkg = sep.Background(img)
    sub = img - bkg
    err = np.sqrt(np.clip(img, 0, None) / gain + (readnoise / gain) ** 2).astype(np.float32)

    try:
        objs = sep.extract(sub, pc.detect_sigma, err=bkg.globalrms, minarea=5)
    except Exception as exc:                                    # noqa: BLE001
        log.debug("sep.extract failed: %s", exc)
        return pd.DataFrame(), bkg
    if len(objs) == 0:
        return pd.DataFrame(), bkg

    flux, fluxerr, _ = sep.sum_circle(sub, objs["x"], objs["y"], 3.0, err=err)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = flux / fluxerr
    cand = pd.DataFrame({
        "x": objs["x"], "y": objs["y"], "flux": flux, "snr": snr,
        "peak": objs["peak"] + bkg.globalback, "a": objs["a"], "b": objs["b"],
        "flag": objs["flag"],
    })

    ny, nx = data.shape
    margin = pc.rmax + 2.0
    keep = ((cand["snr"] > pc.snr_min)
            & (cand["peak"] < pc.saturation_fraction * float(saturate))
            & (cand["flag"] == 0)
            & (cand["a"] / cand["b"].clip(lower=1e-3) < 1.5)
            & (cand["x"] > margin) & (cand["x"] < nx - 1 - margin)
            & (cand["y"] > margin) & (cand["y"] < ny - 1 - margin))
    if exclude_xy is not None and np.all(np.isfinite(exclude_xy)):
        d_comet = np.hypot(cand["x"] - exclude_xy[0], cand["y"] - exclude_xy[1])
        keep &= d_comet > pc.exclude_comet_pix
    stars = cand[keep].copy()
    if stars.empty:
        return stars, bkg

    # Size cut relative to the population: galaxies are larger, cosmic rays smaller.
    size = np.sqrt(stars["a"] * stars["b"])
    med = float(np.median(size))
    stars = stars[(size > 0.6 * med) & (size < 1.6 * med)]

    # Isolation against EVERY detection, not just the survivors.
    ax_, ay_ = cand["x"].to_numpy(), cand["y"].to_numpy()
    iso = []
    for _, s in stars.iterrows():
        d = np.hypot(ax_ - s["x"], ay_ - s["y"])
        iso.append(int(np.sum(d < pc.star_isolation_pix)) == 1)   # itself only
    stars = stars[np.array(iso, bool)]

    stars = stars.sort_values("flux", ascending=False).head(pc.max_stars)
    return stars.reset_index(drop=True), bkg


def stack_star_profiles(sub, stars, profile_config=None, oversample=1):
    """Normalised, sigma-clipped-median radial profile of the field stars.

    Each star is profiled at native resolution on the background-subtracted
    image, normalised by its central surface brightness, and the stack is the
    sigma-clipped median at every radius across stars.

    Returns
    -------
    stack : pandas.DataFrame
        ``r_pix``, ``median``, ``std``, ``n_stars``.
    per_star : pandas.DataFrame
        Long format, ``star`` index plus the normalised profile, for the CSV.
    """
    pc = profile_config or ProfileConfig()
    profiles, rows = [], []
    for i, s in stars.iterrows():
        prof = radial_profile(sub, s["x"], s["y"], pc.radii, pc.half_width, sky=0.0,
                              oversample=oversample, sigma=pc.sigma, maxiters=pc.maxiters,
                              min_pixels=pc.min_pixels)
        c0 = central_sb(prof)
        if not np.isfinite(c0) or c0 <= 0:
            continue
        norm = prof["mean_clip"].to_numpy() / c0
        profiles.append(norm)
        for r, v in zip(prof["r_pix"], norm):
            rows.append(dict(star=i, x=s["x"], y=s["y"], snr=s["snr"], r_pix=r, sb_norm=v))

    if not profiles:
        return pd.DataFrame({"r_pix": pc.radii, "median": np.nan, "std": np.nan, "n_stars": 0}), \
            pd.DataFrame(rows)

    arr = np.vstack(profiles)
    med, std, n = [], [], []
    for col in arr.T:
        good = col[np.isfinite(col)]
        if good.size == 0:
            med.append(np.nan); std.append(np.nan); n.append(0); continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, m, s = sigma_clipped_stats(good, sigma=pc.sigma, maxiters=pc.maxiters)
        med.append(float(m)); std.append(float(s)); n.append(int(good.size))
    stack = pd.DataFrame({"r_pix": pc.radii, "median": med, "std": std, "n_stars": n})
    return stack, pd.DataFrame(rows)


# ------------------------------------------------------------------ per frame
def profile_frame(path, x, y, sky, gain, readnoise, saturate, fwhm_pix,
                  profile_config=None):
    """Comet and field-star profiles for one frame.

    Returns
    -------
    comet : pandas.DataFrame
        Oversampled clipped/unclipped profile plus a native clipped one, all
        normalised to the comet's central SB (``*_norm`` columns).
    stars : pandas.DataFrame
        The stacked star profile.
    per_star : pandas.DataFrame
    summary : dict
        Slopes, central SB, star count, and the comet/star excess at 3 FWHM.
    """
    pc = profile_config or ProfileConfig()
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)

    stars, bkg = select_field_stars(data, gain, readnoise, saturate, exclude_xy=(x, y),
                                    profile_config=pc)
    sub = np.ascontiguousarray(np.where(np.isfinite(data), data, 0.0), np.float32) - bkg
    stack, per_star = stack_star_profiles(
        sub, stars, pc, oversample=pc.oversample if pc.oversample_stars else 1)

    sky_level = float(sky) if np.isfinite(sky) else float(bkg.globalback)
    comet = radial_profile(data, x, y, pc.radii, pc.half_width, sky=sky_level,
                           oversample=pc.oversample, sigma=pc.sigma, maxiters=pc.maxiters,
                           min_pixels=pc.min_pixels)
    native = radial_profile(data, x, y, pc.radii, pc.half_width, sky=sky_level,
                            oversample=1, sigma=pc.sigma, maxiters=pc.maxiters,
                            min_pixels=pc.min_pixels)
    comet["mean_clip_native"] = native["mean_clip"].to_numpy()

    c0 = central_sb(comet)
    for col in ("mean_clip", "mean", "mean_clip_native"):
        comet[f"{col}_norm"] = comet[col] / c0 if np.isfinite(c0) and c0 > 0 else np.nan
    comet["star_median_norm"] = stack["median"].to_numpy()
    comet["star_std_norm"] = stack["std"].to_numpy()
    comet["n_stars"] = stack["n_stars"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        comet["excess"] = comet["mean_clip_norm"] / comet["star_median_norm"]

    rmin = max(pc.fit_rmin_fwhm * float(fwhm_pix), 1.5) if np.isfinite(fwhm_pix) else 2.0
    slope_c, _, n_c = fit_powerlaw(comet["r_pix"], comet["mean_clip_norm"], rmin, pc.fit_rmax_pix)
    slope_u, _, _ = fit_powerlaw(comet["r_pix"], comet["mean_norm"], rmin, pc.fit_rmax_pix)
    slope_s, _, n_s = fit_powerlaw(stack["r_pix"], stack["median"], rmin, pc.fit_rmax_pix)

    r_ex = 3.0 * float(fwhm_pix) if np.isfinite(fwhm_pix) else 6.0
    idx = int(np.argmin(np.abs(comet["r_pix"].to_numpy() - r_ex)))
    summary = dict(
        central_sb=c0, sky=sky_level,
        n_stars=int(stack["n_stars"].max()) if len(stack) else 0,
        n_stars_selected=int(len(stars)),
        slope_comet=slope_c, slope_comet_unclipped=slope_u, slope_star=slope_s,
        n_fit_comet=n_c, n_fit_star=n_s, fit_rmin_pix=rmin, fit_rmax_pix=pc.fit_rmax_pix,
        excess_r_pix=float(comet["r_pix"].iloc[idx]),
        excess_at_3fwhm=float(comet["excess"].iloc[idx]),
    )
    return comet, stack, per_star, summary


# ------------------------------------------------------- PSF-convolved model
def _annulus_means(grid, rr, valid, radii, half_width):
    """Plain annulus means on a grid whose sample radii are already known."""
    out = np.full(len(radii), np.nan)
    for i, r in enumerate(radii):
        sel = valid & (np.abs(rr - r) < half_width)
        if sel.any():
            out[i] = grid[sel].mean()
    return out


def psf_kernel_from_stack(stack, oversample=4, half_size=14.0, slope_fallback=-4.0):
    """Azimuthally symmetric PSF on an oversampled grid, from the star stack.

    The stacked star profile *is* the PSF measurement for the frame, so use it
    directly rather than assume a Gaussian or Moffat.  Beyond the last measured
    radius the kernel is extrapolated with the star profile's own power law.

    Returns
    -------
    kernel : ndarray
        Normalised so that ``kernel.sum() == 1`` (unit total flux).
    """
    r_k = stack["r_pix"].to_numpy(float)
    k = stack["median"].to_numpy(float)
    good = np.isfinite(k) & (k > 0)
    r_k, k = r_k[good], k[good]
    if r_k.size < 4:
        raise ValueError("star stack too sparse for a kernel")
    slope, intercept, _ = fit_powerlaw(r_k, k, max(r_k.min(), 3.0), r_k.max())
    if not np.isfinite(slope):
        slope, intercept = slope_fallback, np.log10(k[-1]) - slope_fallback * np.log10(r_k[-1])

    n = int(2 * half_size * oversample) + 1
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(xx - c, yy - c) / oversample
    kern = np.interp(rr, r_k, k, left=k[0])
    tail = rr > r_k.max()
    kern[tail] = 10 ** (intercept + slope * np.log10(rr[tail]))
    return kern / kern.sum()


def fit_coma_model(comet, stack, oversample=4, half_size=14.0, radii=None,
                   half_width=0.25, fit_rmax=10.0, m_fixed=None, soft_pix=0.25,
                   fit_sky=False):
    """Forward-model the comet profile as nucleus + PSF-convolved power-law coma.

    A power law fitted *outside* the PSF core is biased wherever the core is a
    large fraction of the usable range -- which is exactly the situation at
    large observer distance, where 10 px is many thousands of km and the fit
    window shrinks toward the PSF.  Modelling the PSF explicitly removes that
    bias: the intrinsic profile is

    .. math:: I(\\rho) = F_n\\,\\delta(\\rho) + C\\,\\rho^{-m}

    convolved with the empirical PSF from the field-star stack, then averaged
    in the same annuli as the data.  ``m`` is the corrected coma slope, and
    ``F_n`` the nucleus (point-source) flux.

    Parameters
    ----------
    comet : pandas.DataFrame
        From :func:`radial_profile` (oversampled), sky-subtracted DN/px.
    stack : pandas.DataFrame
        Normalised star stack from :func:`stack_star_profiles`.
    m_fixed : float, optional
        Fix the coma slope (e.g. 1.0 for a steady-state test) and fit only the
        two amplitudes; the returned chi-square can be compared with the free fit.
    fit_sky : bool
        Add a free constant ``s`` to the model.  On a faint comet the outer
        annuli sit within a few sigma of the sky level, so an error of a few DN
        in the sky estimate steepens (or flattens) the whole outer profile; a
        free offset absorbs it, and the inner annuli -- far above the sky --
        still pin the slope.

    Returns
    -------
    dict
        ``m``, ``F_nuc``, ``C``, ``sky_offset``, ``nucleus_fraction`` (of model
        flux inside ``fit_rmax``), ``chi2_red``, ``n_points``, ``model`` (model
        annulus means), ``success``.
    """
    from scipy.optimize import least_squares

    radii = np.asarray(DEFAULT_RADII if radii is None else radii, float)
    obs = comet["mean_clip"].to_numpy(float)
    err = (comet["std_clip"].to_numpy(float) / np.sqrt(np.maximum(comet["n_clip"].to_numpy(float), 1)))
    keep = np.isfinite(obs) & np.isfinite(err) & (radii <= fit_rmax)
    if keep.sum() < 5:
        return dict(m=np.nan, F_nuc=np.nan, C=np.nan, sky_offset=np.nan, nucleus_fraction=np.nan,
                    chi2_red=np.nan, n_points=int(keep.sum()), model=None, success=False)
    err = np.where(err > 0, err, np.nanmedian(err[err > 0]) if np.any(err > 0) else 1.0)
    err = np.maximum(err, 0.01 * np.abs(obs) + 1e-3)

    kernel = psf_kernel_from_stack(stack, oversample, half_size)
    n = kernel.shape[0]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(xx - c, yy - c) / oversample          # px
    valid = np.ones_like(rr, bool)
    rho = np.maximum(rr, soft_pix)

    from scipy.signal import fftconvolve

    # Nucleus: a point source of unit flux, in DN/px on the sub-pixel grid.
    # A unit-flux kernel already IS the PSF image of one unit of flux spread
    # over sub-pixels; per original pixel the SB is os^2 times that.
    nuc_sb = kernel * oversample ** 2
    # Coma: SB map (DN/px) convolved with the unit-flux kernel.
    def coma_sb(m):
        return fftconvolve(rho ** (-m), kernel, mode="same")

    ann_nuc = _annulus_means(nuc_sb, rr, valid, radii, half_width)

    # No caching keyed on m: least_squares probes derivatives with ~1e-8 steps,
    # and a rounded cache key returned the same profile for every probe -- a
    # zero gradient, so the slope never moved from its starting value.
    def model_profile(F, C, m, s=0.0):
        return F * ann_nuc + C * _annulus_means(coma_sb(m), rr, valid, radii, half_width) + s

    # Starting values: coma amplitude from r ~ 5 px, nucleus from the excess
    # of the innermost point over the coma there.
    i5 = int(np.argmin(np.abs(radii - 5.0)))
    m0 = 1.0 if m_fixed is None else float(m_fixed)
    prof_c0 = _annulus_means(coma_sb(m0), rr, valid, radii, half_width)
    C0 = max(obs[i5] / prof_c0[i5], 1e-6) if np.isfinite(obs[i5]) and prof_c0[i5] > 0 else 1.0
    F0 = max((obs[keep][0] - C0 * prof_c0[keep][0]) / max(ann_nuc[keep][0], 1e-9), 0.0)

    # Sky offset bounded to a few times the outermost observed SB: enough to
    # absorb a mis-estimated sky, not enough to swallow the coma.
    s_bound = 3.0 * float(np.nanmax(np.abs(obs[keep][-3:]))) + 1e-3
    x0, lo, hi, scale = [np.log10(F0 + 1e-3), np.log10(C0)], [-6, -6], [8, 8], [1.0, 1.0]
    if m_fixed is None:
        x0.append(m0); lo.append(0.0); hi.append(3.0); scale.append(0.3)
    if fit_sky:
        x0.append(0.0); lo.append(-s_bound); hi.append(s_bound); scale.append(max(s_bound / 3, 1e-3))

    def unpack(pv):
        F, C = 10 ** pv[0], 10 ** pv[1]
        i = 2
        m = m0
        if m_fixed is None:
            m = pv[i]; i += 1
        s = pv[i] if fit_sky else 0.0
        return F, C, m, s

    def resid(pv):
        F, C, m, s = unpack(pv)
        return (model_profile(F, C, m, s)[keep] - obs[keep]) / err[keep]

    try:
        sol = least_squares(resid, x0, bounds=(lo, hi), max_nfev=400,
                            diff_step=1e-4, x_scale=scale)
    except Exception as exc:                                    # noqa: BLE001
        log.debug("coma model fit failed: %s", exc)
        return dict(m=np.nan, F_nuc=np.nan, C=np.nan, sky_offset=np.nan, nucleus_fraction=np.nan,
                    chi2_red=np.nan, n_points=int(keep.sum()), model=None, success=False)

    F, C, m, s = unpack(sol.x)
    m = float(m)
    model = model_profile(F, C, m, s)
    dof = max(int(keep.sum()) - len(sol.x), 1)
    chi2 = float(np.sum(sol.fun ** 2))

    # Flux fractions inside fit_rmax, from the model components.
    inside = rr <= fit_rmax
    f_nuc = F * nuc_sb[inside].sum() / oversample ** 2
    f_com = C * coma_sb(m)[inside].sum() / oversample ** 2
    return dict(m=m, F_nuc=float(F), C=float(C), sky_offset=float(s),
                nucleus_fraction=float(f_nuc / (f_nuc + f_com)) if (f_nuc + f_com) > 0 else np.nan,
                chi2_red=chi2 / dof, n_points=int(keep.sum()), model=model,
                success=bool(sol.success))


# --------------------------------------------------------------------- driver
def run_profiles(target, table, datadir, profile_config=None, progress=True,
                 save=True, plots=True, dpi=60, elements=None):
    """Profiles for every frame of a target, from its photometry table.

    Parameters
    ----------
    target : ztfcomet.config.Target
    table : pandas.DataFrame
        Photometry table.  Only one row per frame is used (the smallest
        aperture), for the centroid and sky level.
    datadir : path-like
    plots : bool
        Per-frame figures into ``fig/<target>/profile/`` plus one summary.
    elements : ztfcomet.orbit.PerihelionInfo, optional
        Lets the summary plot use ``r_h - q`` on its abscissa, matching the
        Af-rho figures.

    Returns
    -------
    profiles : pandas.DataFrame
        Long format, one row per (frame, radius).
    summary : pandas.DataFrame
        One row per frame with the fitted slopes and star counts.
    """
    from . import directory as d
    from tqdm.auto import tqdm

    pc = profile_config or ProfileConfig()
    datadir = Path(datadir)
    if table is None or table.empty:
        return pd.DataFrame(), pd.DataFrame()

    frames = table.copy()
    if "rho_km" in frames:
        frames = frames[np.isclose(frames["rho_km"], frames["rho_km"].min())]
    frames = frames.drop_duplicates("file").reset_index(drop=True)

    carry = [c for c in ("obsjd", "isot", "filter", "r", "delta", "r_rate", "alpha",
                         "fwhm_pix", "pixscale", "quality_ok", "flags", "afrho0_cm")
             if c in frames.columns]

    long_rows, sum_rows, star_rows = [], [], []
    rows = frames.iterrows()
    if progress:
        rows = tqdm(rows, total=len(frames), desc=f"{target.name}: profiles")
    for _, row in rows:
        path = datadir / row["file"]
        x, y = row.get("x_center", np.nan), row.get("y_center", np.nan)
        if not path.exists() or not (np.isfinite(x) and np.isfinite(y)):
            continue
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
            comet, stack, per_star, summ = profile_frame(
                path, x, y, row.get("msky", np.nan), hdr.get("GAIN", row.get("egain", 6.2)),
                hdr.get("READNOI", row.get("readnoise", 9.7)), hdr.get("SATURATE", 6e4),
                row.get("fwhm_pix", np.nan), pc)
        except Exception as exc:                                # noqa: BLE001
            log.warning("%s: profile failed on %s: %s", target.name, row["file"], exc)
            continue

        base = {"target": target.name, "file": row["file"], **{c: row[c] for c in carry}}
        scales = physical_scales(row, fit_rmin_pix=summ["fit_rmin_pix"], rmax_pix=pc.rmax)
        for _, pr in comet.iterrows():
            long_rows.append({**base, **pr.to_dict(), "rho_km": pr["r_pix"] * scales["km_per_pix"]})
        sum_rows.append({**base, **summ, **scales})
        if not per_star.empty:
            per_star = per_star.assign(file=row["file"])
            star_rows.append(per_star)

        if plots:
            try:
                figdir = d.fig_dir(target.name) / "profile"
                figdir.mkdir(parents=True, exist_ok=True)
                _plot_frame(comet, stack, summ, base, figdir / f"{Path(row['file']).stem}.png", dpi)
            except Exception as exc:                            # noqa: BLE001
                log.debug("profile plot failed for %s: %s", row["file"], exc)

    profiles = pd.DataFrame(long_rows)
    summary = pd.DataFrame(sum_rows)
    stars = pd.concat(star_rows, ignore_index=True) if star_rows else pd.DataFrame()

    if save and not profiles.empty:
        slug = d.target_slug(target.name)
        out = d.result_dir(target.name)
        profiles.to_csv(out / f"profile_{slug}.csv", index=False)
        summary.to_csv(out / f"profile_summary_{slug}.csv", index=False)
        if not stars.empty:
            stars.to_csv(out / f"profile_stars_{slug}.csv", index=False)
        log.info("%s: profiles for %d frames -> %s", target.name, len(summary), out)
    if plots and not summary.empty:
        try:
            _plot_summary(target, profiles, summary, d.fig_dir(target.name)
                          / f"profile_summary_{d.target_slug(target.name)}.png",
                          elements=elements)
        except Exception as exc:                                # noqa: BLE001
            log.warning("%s: profile summary plot failed: %s", target.name, exc)
    return profiles, summary


# -------------------------------------------------------------------- figures
def _plot_frame(comet, stack, summ, meta, outpath, dpi=60):
    import matplotlib.pyplot as plt
    from . import rcparams  # noqa: F401

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    r = comet["r_pix"].to_numpy()
    ax.plot(r, comet["mean_clip_norm"], "o-", color="k", ms=4, label="comet, clipped mean (oversampled)")
    ax.plot(r, comet["mean_norm"], ":", color="0.4", lw=1.5, label="comet, plain mean")
    ax.plot(r, comet["mean_clip_native_norm"], "--", color="tab:blue", lw=1.2, label="comet, native pixels")
    if stack["n_stars"].max() > 0:
        m, s = stack["median"].to_numpy(), stack["std"].to_numpy()
        ax.plot(r, m, "s-", color="tab:red", ms=4,
                label=f"field stars, median (N={int(stack['n_stars'].max())}, same grid)")
        ax.fill_between(r, np.clip(m - s, 1e-4, None), m + s, color="tab:red", alpha=0.15, lw=0)
    # 1/rho reference anchored on the comet at the inner fit radius.
    rmin = summ.get("fit_rmin_pix", 2.0)
    anchor = comet.loc[(comet["r_pix"] - rmin).abs().idxmin()]
    if np.isfinite(anchor["mean_clip_norm"]) and anchor["mean_clip_norm"] > 0:
        rr = np.linspace(max(rmin, 0.5), r.max(), 50)
        ax.plot(rr, anchor["mean_clip_norm"] * anchor["r_pix"] / rr, color="tab:green", lw=1.2,
                ls="-.", label=r"$\rho^{-1}$ (steady-state coma)")
    ax.set(xscale="log", yscale="log", xlabel="radius (pix)", ylabel="SB / central SB")
    ax.set_ylim(bottom=max(1e-4, np.nanmin(comet["mean_clip_norm"].clip(lower=1e-4)) / 3))
    lines = [f"{meta.get('isot', '')[:10]}  {meta.get('filter', '')}",
             f"r_h = {meta.get('r', np.nan):.2f} AU   FWHM = {meta.get('fwhm_pix', np.nan):.2f} px",
             f"slope: comet {summ['slope_comet']:+.2f}  stars {summ['slope_star']:+.2f}  (coma = -1)",
             f"excess at {summ['excess_r_pix']:.1f} px: {summ['excess_at_3fwhm']:.1f}x"]
    ax.annotate("\n".join(lines), xy=(0.02, 0.02), xycoords="axes fraction", fontsize=9,
                ha="left", va="bottom", family="monospace",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.set_title(f"{meta.get('target', '')}  {Path(meta.get('file', '')).stem[:24]}", fontsize=11)
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def _plot_summary(target, profiles, summary, outpath, elements=None):
    import matplotlib.pyplot as plt
    from . import rcparams  # noqa: F401

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.8))
    good = summary["quality_ok"].astype(bool) if "quality_ok" in summary else np.ones(len(summary), bool)

    # (a) every frame's normalised comet profile, the median, and the star stack
    piv = profiles.pivot_table(index="file", columns="r_pix", values="mean_clip_norm")
    r = piv.columns.to_numpy(float)
    for _, prof in piv.iterrows():
        ax1.plot(r, prof.to_numpy(), color="0.75", lw=0.6, alpha=0.6)
    ax1.plot(r, np.nanmedian(piv.to_numpy(), axis=0), "o-", color="k", ms=4,
             label=f"comet, median of {len(piv)} frames")
    spiv = profiles.pivot_table(index="file", columns="r_pix", values="star_median_norm")
    ax1.plot(r, np.nanmedian(spiv.to_numpy(), axis=0), "s-", color="tab:red", ms=4,
             label="field stars, median stack")
    ref = np.nanmedian(piv.to_numpy(), axis=0)
    i0 = int(np.argmin(np.abs(r - 3.0)))
    if np.isfinite(ref[i0]) and ref[i0] > 0:
        ax1.plot(r[r >= 2], ref[i0] * 3.0 / r[r >= 2], "-.", color="tab:green", label=r"$\rho^{-1}$")
    ax1.set(xscale="log", yscale="log", xlabel="radius (pix)", ylabel="SB / central SB",
            title=f"{target.name} — radial profiles")
    ax1.set_ylim(bottom=1e-3)
    ax1.legend(fontsize=9, frameon=False)

    # (b) fitted slope through the apparition, on the same abscissa as the
    #     Af-rho figures: r_h - q signed by leg, perihelion at zero.
    have_q = elements is not None and np.isfinite(elements.get("q", np.nan))
    x = summary["r"] if "r" in summary else np.arange(len(summary))
    if "r_rate" in summary and have_q:
        x = (summary["r"] - elements["q"]) * np.where(summary["r_rate"] < 0, -1, 1)
        ax2.axvline(0, color="crimson", ls="--", lw=1.2, zorder=0)
        ax2.set_xlabel(r"$r_\mathrm{h} - q$ (AU)   $\leftarrow$ pre    post $\rightarrow$")
    elif "r_rate" in summary:
        x = summary["r"] * np.where(summary["r_rate"] < 0, -1, 1)
        ax2.set_xlabel(r"$\pm r_\mathrm{h}$ (AU; negative = inbound)")
    else:
        ax2.set_xlabel(r"$r_\mathrm{h}$ (AU)")
    ax2.scatter(x[good], summary.loc[good, "slope_comet"], c="k", s=22, label="comet (clean)")
    ax2.scatter(x[~good], summary.loc[~good, "slope_comet"], facecolors="none", edgecolors="k",
                s=22, label="comet (flagged)")
    ax2.scatter(x, summary["slope_star"], c="tab:red", s=10, alpha=0.6, label="field stars")
    ax2.axhline(-1, color="tab:green", ls="-.", label="steady-state coma")
    ax2.set(ylabel="power-law slope of SB(r)", title="profile slope through the apparition")
    ax2.legend(fontsize=9, frameon=False)

    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
