"""
Aperture photometry for one cutout.

Design decisions and the review items they answer
-------------------------------------------------
**Units (S1).**  Science pixels are mJy/pixel, so an aperture sum is already in
mJy and no solid-angle conversion is applied.  ``Config.flux_unit`` records this
in every output file so the assumption is never implicit again.

**Fixed annulus (S2).**  The background is measured in a fixed 15-20 pixel
annulus for every aperture and every exposure.  Scaling the annulus to the
largest aperture, as the primitive code did, pushed it outside the 91-pixel
cutout for any target closer than ~1.5 au.  A fixed annulus also keeps the
background sample close to the target, which matters because SPHEREx's central
wavelength varies across the detector: a distant annulus samples a different
bandpass than the aperture it is meant to correct.

**No centroiding (S4).**  ``xcen``/``ycen`` from the index are used as given.

**Negative fluxes are kept (S5).**  A negative aperture sum is a legitimate
measurement of a faint source.  Discarding them, as the primitive ``badphot``
did, truncates the noise distribution and biases any subsequent average upward.
``abmag`` is ``nan`` where the flux is not positive -- the magnitude is
undefined there -- but the flux and its uncertainty are always reported.

**Uncertainties (S6).**  The variance plane already contains the background's
photon noise, so the DAOPHOT ``A * ssky^2`` term would count it a second time.
The default combines the propagated pixel variance with the uncertainty on the
*sky level* only.  All three components are reported separately, plus a
``sky_excess_ratio`` diagnostic that says whether the variance plane is
under-reporting the real scatter.

**Bad pixels.**  Excluded, not zeroed: the aperture sums the good pixels and
its effective area shrinks accordingly.  No flux is invented to fill the gap;
``badphot`` marks any aperture containing one, while bad pixels falling only in
the annulus do not raise it.
"""

from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
from astropy.stats import SigmaClip
from photutils.aperture import (
    ApertureStats,
    CircularAnnulus,
    CircularAperture,
    aperture_photometry,
)

from .apertures import ApertureSpec, annulus_radii_pix
from .logging_utils import get_logger
from .sourceflag import SOURCEFLAG_COLUMNS, GaiaPixelField, evaluate_sourceflag

__all__ = ["PHOT_COLUMNS", "DISTCORR_COLUMNS", "SkyEstimate", "estimate_sky",
           "measure_exposure", "measure_exposure_records",
           "add_distance_corrected_flux"]

log = get_logger("phot")

#: Columns produced per (exposure, aperture).  Metadata columns are appended by
#: :mod:`pipeline`; this list is the photometric payload only.
PHOT_COLUMNS = [
    "ap_label", "ap_kind", "r_ap_pix", "r_ap_km", "r_ap_arcsec",
    "r_in_pix", "r_out_pix", "pixel_scale_km",
    "aperture_area_pix2", "aperture_area_eff_pix2",
    "n_badpix_ap", "badpix_area_ap_pix2", "frac_badpix_ap", "badphot",
    "sky_median_mjy_per_pix", "sky_std_mjy_per_pix", "sky_area_pix2",
    "n_badpix_sky", "sky_excess_ratio",
    "aperture_sum_mjy", "source_sum_mjy",
    "err_pix_mjy", "err_skylevel_mjy", "err_skyscatter_mjy",
    "source_sum_err_mjy", "source_sum_err_empirical_mjy",
    "snr", "abmag", "abmag_err",
] + SOURCEFLAG_COLUMNS

#: Added after the index metadata is merged on, because they need ``r_hel`` and
#: ``r_obs``, which live in the index rather than in the cutout.
DISTCORR_COLUMNS = ["flux_distcorr_mjy", "flux_distcorr_err_mjy", "distcorr_factor"]


class SkyEstimate:
    """Sigma-clipped background statistics from the annulus."""

    __slots__ = ("median", "std", "area", "n_badpix", "var_median", "ok", "r_in", "r_out")

    def __init__(self, median: float, std: float, area: float,
                 n_badpix: int, var_median: float, ok: bool,
                 r_in: float = float("nan"), r_out: float = float("nan")):
        self.median = median
        self.std = std
        self.area = area
        self.n_badpix = n_badpix
        self.var_median = var_median
        self.ok = ok
        self.r_in = r_in          #: annulus radii actually used [pixel]
        self.r_out = r_out

    @property
    def excess_ratio(self) -> float:
        """
        Measured annulus scatter divided by the scatter the variance plane predicts.

        A value near 1 means the variance plane is a complete noise model; a
        value well above 1 means it is optimistic (confusion, flat-field or
        zodiacal structure that the plane does not describe) and that the
        reported ``source_sum_err_mjy`` should be treated as a lower bound.
        """
        if not (np.isfinite(self.std) and np.isfinite(self.var_median)) or self.var_median <= 0:
            return float("nan")
        return float(self.std / np.sqrt(self.var_median))

    def __repr__(self) -> str:
        return (f"SkyEstimate(median={self.median:.4g}, std={self.std:.4g}, "
                f"area={self.area:.1f}, ok={self.ok})")


def _scalar(x) -> float:
    """Unwrap a length-1 photutils result, which may carry an astropy unit."""
    v = x[0] if np.ndim(x) else x
    return float(getattr(v, "value", v))


def estimate_sky(
    sci: np.ndarray,
    var: np.ndarray,
    badpix: np.ndarray,
    xcen: float,
    ycen: float,
    cfg,
    r_in: Optional[float] = None,
    r_out: Optional[float] = None,
) -> SkyEstimate:
    """
    Measure the local background in the sky annulus.

    Parameters
    ----------
    sci, var, badpix : ndarray
        Science plane, variance plane and bad-pixel mask.
    xcen, ycen : float
        Target centroid, 0-indexed.
    cfg : Config
        Supplies ``sky_sigma``, ``sky_maxiters`` and ``min_sky_area``, and the
        fixed ``r_in_pix`` / ``r_out_pix`` ring when no radii are given.
    r_in, r_out : float, optional
        Annulus radii [pixel] for this exposure, normally from
        :func:`apertures.annulus_radii_pix` (the physical 150 000 km rule).

    Returns
    -------
    SkyEstimate
        ``ok`` is ``False`` when the clipped annulus area falls below
        ``cfg.min_sky_area``; the median and scatter are then ``nan`` so a bad
        background cannot silently propagate into every aperture.
    """
    if r_in is None or r_out is None:
        r_in, r_out = float(cfg.r_in_pix), float(cfg.r_out_pix)
    ann = CircularAnnulus([(float(xcen), float(ycen))], r_in=float(r_in), r_out=float(r_out))
    sigclip = SigmaClip(sigma=float(cfg.sky_sigma), maxiters=int(cfg.sky_maxiters))

    stats = ApertureStats(sci, ann, mask=badpix, sigma_clip=sigclip)
    median = _scalar(stats.median)
    std = _scalar(stats.std)
    area = _scalar(stats.sum_aper_area)

    # Bad pixels are reported for diagnostics only: by the project's rule they
    # never raise `badphot`, which concerns the photometric aperture alone.
    # Whole-pixel counting, to match `n_badpix_ap`.
    n_bad = _scalar(aperture_photometry(badpix.astype(np.float64), ann,
                                        method="center")["aperture_sum"])

    var_stats = ApertureStats(np.where(badpix, np.nan, var), ann, mask=badpix)
    var_median = _scalar(var_stats.median)

    ok = np.isfinite(median) and np.isfinite(area) and area >= float(cfg.min_sky_area)
    if not ok:
        median = std = float("nan")
    return SkyEstimate(median=median, std=std, area=area,
                       n_badpix=int(round(n_bad)), var_median=var_median, ok=ok,
                       r_in=float(r_in), r_out=float(r_out))


def measure_exposure_records(
    sci: np.ndarray,
    var: np.ndarray,
    badpix: np.ndarray,
    xcen: float,
    ycen: float,
    specs: Sequence[ApertureSpec],
    cfg,
    *,
    pixel_scale_km: float = float("nan"),
    psf_fwhm_pix: float = float("nan"),
    vmag: float = float("nan"),
    gaia_field: Optional[GaiaPixelField] = None,
    sky: Optional[SkyEstimate] = None,
) -> List[dict]:
    """
    Measure every valid aperture on one cutout, returning plain dicts.

    Parameters
    ----------
    sci, var, badpix : ndarray
        Science plane [mJy/pixel], variance plane [(mJy/pixel)^2] and bad-pixel
        mask, all the same shape.
    xcen, ycen : float
        Target centroid, 0-indexed pixel coordinates, used as given (S4).
    specs : sequence of ApertureSpec
        Apertures resolved for this exposure.  Entries with ``valid=False`` are
        not measured; whether they appear as NaN rows is set by
        ``cfg.drop_invalid_apertures``.
    cfg : Config
    pixel_scale_km, psf_fwhm_pix, vmag : float
        Per-exposure geometry, copied into the output and used by the source
        flags.
    gaia_field : GaiaPixelField, optional
        Projected Gaia sources.  ``None`` yields empty contamination columns.
    sky : SkyEstimate, optional
        Pre-computed background; recomputed here when omitted.

    Returns
    -------
    list of dict
        One entry per measured aperture, keyed by :data:`PHOT_COLUMNS`.  Plain
        dicts rather than a DataFrame because the pipeline measures ~100k
        exposures: building one small DataFrame each and concatenating them is
        far slower than accumulating records and building one frame at the end.
        :func:`measure_exposure` wraps this for interactive use.

    Notes
    -----
    Cost is three ``aperture_photometry`` calls and three ``ApertureStats``
    calls per *exposure*, not per aperture: photutils accepts a list of
    apertures and returns ``aperture_sum_0 ... aperture_sum_n``.  The primitive
    code built four ``ApertureStats`` objects per aperture, which for 21
    apertures over the full catalogue is millions of redundant objects.
    """
    valid = [s for s in specs if s.valid]
    if not valid:
        return []

    if not np.isfinite(xcen) or not np.isfinite(ycen):
        log.warning("non-finite centroid (%r, %r); skipping exposure", xcen, ycen)
        return []

    sci = np.asarray(sci, dtype=np.float64)
    var = np.asarray(var, dtype=np.float64)
    badpix = np.asarray(badpix, dtype=bool)

    if sky is None:
        r_in, r_out = annulus_radii_pix(cfg, pixel_scale_km)
        sky = estimate_sky(sci, var, badpix, xcen, ycen, cfg, r_in, r_out)

    # Sanitise before photutils sees the arrays: masked pixels are set to zero
    # explicitly so that a NaN can never leak into a sum even if photutils'
    # own mask handling changes.
    sci_clean = np.where(badpix, 0.0, np.nan_to_num(sci, nan=0.0, posinf=0.0, neginf=0.0))
    var_clean = np.where(badpix, 0.0, np.nan_to_num(var, nan=0.0, posinf=0.0, neginf=0.0))
    err_clean = np.sqrt(var_clean)
    badpix_f = badpix.astype(np.float64)

    pos = [(float(xcen), float(ycen))]
    aps = [CircularAperture(pos, r=float(s.r_ap_pix)) for s in valid]

    phot = aperture_photometry(sci_clean, aps, error=err_clean, mask=badpix)
    # Exact-overlap and whole-pixel bad-pixel accounting.  `A_geom - exact
    # overlap` reproduces the masked sum's effective area to machine precision,
    # which is why the area does not need a separate ApertureStats call.
    bad_exact = aperture_photometry(badpix_f, aps)
    bad_count = aperture_photometry(badpix_f, aps, method="center")

    ny, nx = sci.shape
    records: List[dict] = []
    for i, spec in enumerate(valid):
        # photutils appends `_i` whenever a *list* of apertures is passed --
        # including a list of one, which is why the suffix is unconditional.
        sfx = f"_{i}"
        ap_sum = float(phot[f"aperture_sum{sfx}"][0])
        err_pix = float(phot[f"aperture_sum_err{sfx}"][0])
        area_geom = float(aps[i].area)
        badpix_area = float(bad_exact[f"aperture_sum{sfx}"][0])
        n_badpix = int(round(float(bad_count[f"aperture_sum{sfx}"][0])))
        area_eff = max(area_geom - badpix_area, 0.0)

        # An aperture that reaches outside the array would have an effective
        # area smaller than pi r^2 for a reason photutils does not report; the
        # fixed r_in = 15 px annulus makes this essentially impossible, so it is
        # a guard rather than an expected path.
        if (xcen - spec.r_ap_pix < 0 or ycen - spec.r_ap_pix < 0
                or xcen + spec.r_ap_pix > nx or ycen + spec.r_ap_pix > ny):
            log.warning("aperture %s (r=%.2f px) extends past the array edge at "
                        "(%.2f, %.2f) in a %dx%d cutout",
                        spec.label, spec.r_ap_pix, xcen, ycen, nx, ny)

        source_sum = ap_sum - area_eff * sky.median if sky.ok else float("nan")

        # --- uncertainty budget (S6) ---------------------------------------
        if sky.ok and np.isfinite(sky.std):
            err_skylevel = (area_eff * sky.std / np.sqrt(sky.area)) if sky.area > 0 else float("nan")
            err_skyscatter = np.sqrt(area_eff) * sky.std
        else:
            err_skylevel = err_skyscatter = float("nan")

        if cfg.sky_noise_mode == "level":
            total_sq = err_pix ** 2 + err_skylevel ** 2
        elif cfg.sky_noise_mode == "empirical":
            total_sq = err_skyscatter ** 2 + err_skylevel ** 2
        else:  # "daophot" -- retained only to quantify the change; double-counts
            total_sq = err_pix ** 2 + err_skyscatter ** 2 + err_skylevel ** 2
        source_err = float(np.sqrt(total_sq)) if np.isfinite(total_sq) else float("nan")
        source_err_emp = float(np.sqrt(err_skyscatter ** 2 + err_skylevel ** 2)) \
            if np.isfinite(err_skyscatter) and np.isfinite(err_skylevel) else float("nan")

        snr = source_sum / source_err if (np.isfinite(source_err) and source_err > 0) else float("nan")

        # AB magnitude is undefined for a non-positive flux; the flux itself is
        # still reported (S5).
        if np.isfinite(source_sum) and source_sum > 0:
            abmag = -2.5 * np.log10(source_sum * 1e-3) + cfg.ab_zeropoint
            abmag_err = (2.5 / np.log(10.0)) * (source_err / source_sum) \
                if np.isfinite(source_err) else float("nan")
        else:
            abmag = abmag_err = float("nan")

        sflag = evaluate_sourceflag(
            gaia_field if gaia_field is not None else GaiaPixelField(
                np.empty(0), np.empty(0), np.empty(0), np.empty(0)),
            r_ap_pix=spec.r_ap_pix, psf_fwhm_pix=psf_fwhm_pix,
            vmag=vmag, snr=snr, cfg=cfg,
        )

        rec = {
            "ap_label": spec.label, "ap_kind": spec.kind,
            "r_ap_pix": spec.r_ap_pix, "r_ap_km": spec.r_ap_km,
            "r_ap_arcsec": spec.r_ap_arcsec,
            "r_in_pix": float(sky.r_in), "r_out_pix": float(sky.r_out),
            "pixel_scale_km": float(pixel_scale_km),
            "aperture_area_pix2": area_geom,
            "aperture_area_eff_pix2": area_eff,
            "n_badpix_ap": n_badpix,
            "badpix_area_ap_pix2": badpix_area,
            "frac_badpix_ap": badpix_area / area_geom if area_geom > 0 else float("nan"),
            # By the project's rule, only bad pixels inside the photometric
            # aperture raise badphot; those in the annulus do not.
            "badphot": bool(badpix_area > 0.0) or area_eff <= 0.0,
            "sky_median_mjy_per_pix": sky.median,
            "sky_std_mjy_per_pix": sky.std,
            "sky_area_pix2": sky.area,
            "n_badpix_sky": sky.n_badpix,
            "sky_excess_ratio": sky.excess_ratio,
            "aperture_sum_mjy": ap_sum,
            "source_sum_mjy": source_sum,
            "err_pix_mjy": err_pix,
            "err_skylevel_mjy": err_skylevel,
            "err_skyscatter_mjy": err_skyscatter,
            "source_sum_err_mjy": source_err,
            "source_sum_err_empirical_mjy": source_err_emp,
            "snr": snr, "abmag": abmag, "abmag_err": abmag_err,
        }
        rec.update(sflag.as_dict())
        records.append(rec)

    if not cfg.drop_invalid_apertures:
        for spec in (s for s in specs if not s.valid):
            blank = {c: np.nan for c in PHOT_COLUMNS}
            blank.update(ap_label=spec.label, ap_kind=spec.kind,
                         r_ap_pix=spec.r_ap_pix, r_ap_km=spec.r_ap_km,
                         r_ap_arcsec=spec.r_ap_arcsec, badphot=True,
                         sourceflag="0", n_gaia=0)
            records.append(blank)

    return records


def measure_exposure(*args, **kwargs) -> pd.DataFrame:
    """
    DataFrame form of :func:`measure_exposure_records`.

    Convenient for the notebook and for tests; the batch pipeline uses the
    record form directly.
    """
    return pd.DataFrame(measure_exposure_records(*args, **kwargs), columns=PHOT_COLUMNS)


def add_distance_corrected_flux(df: "pd.DataFrame", cfg) -> "pd.DataFrame":
    """
    Add the distance-corrected flux columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Photometry rows carrying ``source_sum_mjy``, ``source_sum_err_mjy``,
        ``r_hel`` and ``r_obs``.
    cfg : Config
        ``distcorr_mode`` selects the sense of the correction.

    Returns
    -------
    pandas.DataFrame
        A copy with :data:`DISTCORR_COLUMNS` added.  ``distcorr_factor`` is
        recorded alongside so the correction can be undone or re-derived.

    Notes
    -----
    ``distcorr_mode="divide"`` (the configured default) computes

    .. math:: F_{\rm corr} = F / (r_{\rm hel}^2 r_{\rm obs}^2)

    exactly as specified for this project.  Be aware of what that does to the
    geometry: a comet's reflected flux already scales as
    :math:`1/(r_{\rm hel}^2 r_{\rm obs}^2)`, so dividing by the same factor
    *doubles* the distance dependence rather than removing it, and the result
    varies as :math:`1/(r_{\rm hel}^4 r_{\rm obs}^4)` between epochs.  The
    quantity that is constant for an unchanging coma -- and the one
    :mod:`reflectance` uses -- is the ``"multiply"`` form,
    :math:`F \times r_{\rm hel}^2 r_{\rm obs}^2`.  Switching is a one-field
    change to the config; nothing else in the pipeline depends on the choice.
    """
    out = df.copy()
    need = {"source_sum_mjy", "r_hel", "r_obs"} - set(out.columns)
    if need:
        raise KeyError(f"distance correction needs column(s) {sorted(need)}")

    r_hel = out["r_hel"].to_numpy(dtype=float)
    r_obs = out["r_obs"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        d2 = (r_hel ** 2) * (r_obs ** 2)
        factor = np.where(np.isfinite(d2) & (d2 > 0),
                          1.0 / d2 if cfg.distcorr_mode == "divide" else d2,
                          np.nan)

    out["distcorr_factor"] = factor
    out["flux_distcorr_mjy"] = out["source_sum_mjy"].to_numpy(float) * factor
    if "source_sum_err_mjy" in out.columns:
        # Linear in flux, so the uncertainty scales by the same factor.
        out["flux_distcorr_err_mjy"] = np.abs(
            out["source_sum_err_mjy"].to_numpy(float) * factor)
    else:
        out["flux_distcorr_err_mjy"] = np.nan
    return out
