"""
Stacked images of every fitted SPHEREx phase group, for the summary slides.

For each (target, phase) of ``results/comspec/gas_fit.csv`` (every phase of
``results/comspec/phase_map.csv`` with ``--all-phases``) five comet-centred stacks are built and
drawn on one figure, ``fig/comspec/phase_images/<target>_<ap>km_ph<phase>.png``:

  ZTF r       the clean ZTF r-band frames within ``--pad`` days of the SPHEREx window
              (``ztfcomet`` stage, ``results/ztf/photometry/<T>.csv`` and the FITS on its
              data root); with none there, the nearest clean frames within 30 d (up to 5),
              then within 120 d (up to 3), and the time offset is stated on the panel
  SPHEREx     the phase's own exposures whose channel wavelength falls in a window:
              the dust continuum 1.2-2.5 um, and the pipeline's H2O, CO2 and CO emission
              windows (``spherex_comspec.config.EMISSION_WINDOWS``: 2.55-2.80, 4.18-4.35 and
              4.55-4.90 um -- the last one also holds the H2O hot bands)

Every frame is background-subtracted (the photometry's own annulus median for SPHEREx, a
sigma-clipped ring median of the stamp for ZTF, and for a SPHEREx frame without one),
reprojected north-up onto a tangent plane centred on the comet -- the ephemeris position
at *that* frame's time, so the stack is in the comet's rest frame and field stars scatter
-- at the instrument's native pixel scale (ZTF 1.0", SPHEREx 6.2") over the same 310"
field, and combined as a sigma-clipped median with a count map.  Fluxes stay per pixel:
mJy/pixel for SPHEREx, uJy/pixel for ZTF (from each frame's MAGZP).  Each stack is also
written to ``results/comspec/phase_stacks/<stem>_<band>.fits`` (IMAGE, COUNT, a north-up WCS).

The figure carries the fit aperture (red circle), on the ZTF panel also the 10 000 km
aperture of the Afrho series (cyan, dashed), a scale bar in km, the N/E compass, and the
anti-solar and anti-velocity directions of the ZTF ephemeris when ZTF frames exist.

Usage
-----
    python scripts/comspec/phase_images.py                       # every group of gas_fit.csv
    python scripts/comspec/phase_images.py 24P:3 2023RS61 161P   # a group, or every phase of a target
    python scripts/comspec/phase_images.py --workers 4 --pad 5
"""
from __future__ import annotations

import argparse
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from astropy.io import fits                           # noqa: E402
from astropy.stats import sigma_clip, sigma_clipped_stats  # noqa: E402
from astropy.time import Time                         # noqa: E402
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval  # noqa: E402
from astropy.wcs import WCS, FITSFixedWarning         # noqa: E402
from scipy import ndimage                             # noqa: E402

# scripts/comspec/phase_images.py -> the project root holds the three packages;
# notebooks/ holds the shared rcparams.py
ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "notebooks"):
    sys.path.insert(0, str(p))

import rcparams                                       # noqa: E402,F401  (the project's figure style)
from spherex_comspec import directory as cd           # noqa: E402
from spherex_comspec.config import EMISSION_WINDOWS   # noqa: E402
from spherex_apphot import directory as ad            # noqa: E402
from spherex_apphot.config import Config as ApConfig  # noqa: E402
from spherex_apphot.fitsio import FitsResolver, read_cutout  # noqa: E402
from spherex_apphot.masking import build_badpix_mask  # noqa: E402
from spherex_apphot.wcsutil import wcs_from_header    # noqa: E402

try:
    import ztfcomet as zc                             # noqa: E402
    HAVE_ZTF = True
except Exception as exc:                              # noqa: BLE001
    HAVE_ZTF = False
    _ZTF_ERR = str(exc)

warnings.simplefilter("ignore", FITSFixedWarning)

FIG_DIR = cd.FIG_DIR / "phase_images"
STACK_DIR = cd.RESULT_DIR / "phase_stacks"
S = {"target": str}

FOV_RADIUS_ARCSEC = 155.0           # half-width of every stack: 25 SPHEREx pixels
SPHEREX_PIXSCALE = 6.2              # arcsec / pixel of the output SPHEREx grid
ZTF_PIXSCALE = 1.0
KM_PER_ARCSEC_AU = 725.27           # km subtended by 1" at 1 au
ZTF_AB_ZP = 23.9                    # AB magnitude of 1 uJy
WINDOWS = {
    "cont": ("continuum", (1.2, 2.5)),
    "H2O": ("H$_2$O", tuple(EMISSION_WINDOWS["2.7um"])),
    "CO2": ("CO$_2$", tuple(EMISSION_WINDOWS["4.3um"])),
    "CO": ("CO (+ H$_2$O hot bands)", tuple(EMISSION_WINDOWS["4.7um"])),
}
BAR_CHOICES_KM = (1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6)


@dataclass
class Stack:
    image: np.ndarray
    count: np.ndarray
    pixscale: float
    n: int
    unit: str
    files: List[str] = field(default_factory=list)
    ra0: float = np.nan
    dec0: float = np.nan
    note: str = ""


# ------------------------------------------------------------------ geometry
def tan_wcs(ra0: float, dec0: float, pixscale: float, radius: int) -> WCS:
    """North-up, east-left tangent plane centred on (ra0, dec0), *radius* pixels each way."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [float(ra0), float(dec0)]
    w.wcs.crpix = [radius + 1.0, radius + 1.0]
    w.wcs.cdelt = [-pixscale / 3600.0, pixscale / 3600.0]
    w.wcs.cunit = ["deg", "deg"]
    return w


def reproject(img: np.ndarray, wcs_in: WCS, wcs_out: WCS, size: int) -> np.ndarray:
    """
    Resample *img* onto the *wcs_out* grid, bilinear and NaN-safe.

    A zero-filled copy and a weight image are interpolated separately and divided, so
    a bad pixel neither spreads into its neighbours nor biases the result low; where
    the weight drops below one half the output is NaN.  The input WCS may carry SIP
    (SPHEREx) or TPV (ZTF) distortion -- ``all_world2pix`` handles both.
    """
    yy, xx = np.mgrid[0:size, 0:size]
    ra, dec = wcs_out.wcs_pix2world(xx, yy, 0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x, y = wcs_in.all_world2pix(ra, dec, 0, quiet=True)
    good = np.isfinite(img)
    filled = np.where(good, img, 0.0)
    kw = dict(order=1, mode="constant", cval=0.0, prefilter=False)
    num = ndimage.map_coordinates(filled, [y, x], **kw)
    den = ndimage.map_coordinates(good.astype(float), [y, x], **kw)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.5, num / den, np.nan)


def ring_median(stamp: np.ndarray, inner_frac: float = 0.8) -> float:
    """Sigma-clipped median of the stamp's outer ring, the display background."""
    n = stamp.shape[0]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.hypot(xx - c, yy - c)
    ring = stamp[(r >= inner_frac * c) & (r <= c) & np.isfinite(stamp)]
    if ring.size < 20:
        return 0.0
    _, med, _ = sigma_clipped_stats(ring, sigma=3.0, maxiters=5)
    return float(med)


def combine(stamps: List[np.ndarray], sigma: float = 3.0, maxiters: int = 5):
    """Sigma-clipped median along the stack and the number of frames per pixel."""
    cube = np.array(stamps, dtype=np.float64)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*invalid value.*")
        warnings.filterwarnings("ignore", message=".*All-NaN.*")
        clipped = sigma_clip(cube, sigma=sigma, maxiters=maxiters, axis=0, masked=True, copy=False)
        count = (~clipped.mask & np.isfinite(cube)).sum(axis=0).astype(np.int32)
        image = np.ma.median(clipped, axis=0).filled(np.nan)
    image[count == 0] = np.nan
    return image, count


def write_stack(path: Path, st: Stack, target: str, phase: int, band: str, wl: Optional[tuple]):
    path.parent.mkdir(parents=True, exist_ok=True)
    radius = (st.image.shape[0] - 1) // 2
    hdr = tan_wcs(st.ra0, st.dec0, st.pixscale, radius).to_header()
    hdr["TARGET"] = target
    hdr["PHASE"] = int(phase)
    hdr["BAND"] = band
    if wl:
        hdr["WL_MIN"], hdr["WL_MAX"] = (float(wl[0]), "[um]"), (float(wl[1]), "[um]")
    hdr["NCOMBINE"] = st.n
    hdr["BUNIT"] = st.unit
    hdr["PIXSCALE"] = (st.pixscale, "[arcsec/pixel]")
    hdr["COMBINE"] = "sigma-clipped median, north up, comet-centred"
    hdr["COMMENT"] = "CRVAL is the mean ephemeris position of the frames; each frame was centred on its own"
    fits.HDUList([fits.PrimaryHDU(header=hdr), fits.ImageHDU(st.image.astype(np.float32), header=hdr, name="IMAGE"),
                  fits.ImageHDU(st.count.astype(np.int16), name="COUNT")]).writeto(path, overwrite=True)


# ------------------------------------------------------------------ SPHEREx
_APPHOT_CACHE: Dict[str, pd.DataFrame] = {}
APPHOT_COLS = ["filename", "wl", "xcen", "ycen", "sky_median_mjy_per_pix", "pix_scale", "r_obs", "r_hel",
               "alpha", "jd_utc", "ra", "dec", "psf_fwhm_pix", "pixel_scale_km"]


def apphot_rows(target: str) -> pd.DataFrame:
    if target not in _APPHOT_CACHE:
        p = cd.APPHOT_DIR / f"{target}.csv"
        _APPHOT_CACHE[target] = (pd.read_csv(p, usecols=APPHOT_COLS).drop_duplicates("filename")
                                 if p.exists() else pd.DataFrame(columns=APPHOT_COLS))
    return _APPHOT_CACHE[target]


def spherex_stacks(target: str, phase: int, assignment: pd.DataFrame, resolver: FitsResolver,
                   cfg: ApConfig, radius: int) -> Dict[str, Optional[Stack]]:
    files = assignment[(assignment["target"] == target) & (assignment["phase"] == phase)]["filename"]
    rows = apphot_rows(target)
    rows = rows[rows["filename"].isin(set(files))]
    out: Dict[str, Optional[Stack]] = {}
    size = 2 * radius + 1
    for key, (_, (lo, hi)) in WINDOWS.items():
        sel = rows[(rows["wl"] >= lo) & (rows["wl"] <= hi)]
        stamps, used, cen = [], [], []
        for r in sel.itertuples(index=False):
            path = resolver.find(r.filename, target)
            if path is None:
                continue
            try:
                cut = read_cutout(path)
            except Exception:                         # noqa: BLE001
                continue
            mask, _ = build_badpix_mask(cut.sci, cut.var, cut.flag, cfg.bad_flag_bits)
            img = np.where(mask, np.nan, cut.sci)
            w_in = wcs_from_header(cut.header)
            xc, yc = cut.header.get("XCEN", r.xcen), cut.header.get("YCEN", r.ycen)
            ra0, dec0 = (float(v) for v in w_in.all_pix2world(xc, yc, 0))
            stamp = reproject(img, w_in, tan_wcs(ra0, dec0, SPHEREX_PIXSCALE, radius), size)
            sky = float(r.sky_median_mjy_per_pix) if np.isfinite(r.sky_median_mjy_per_pix) else ring_median(stamp)
            stamps.append(stamp - sky)
            used.append(r.filename)
            cen.append((ra0, dec0))
        if not stamps:
            out[key] = None
            continue
        image, count = combine(stamps)
        out[key] = Stack(image, count, SPHEREX_PIXSCALE, len(stamps), "mJy/pixel", used,
                         float(np.mean([c[0] for c in cen])), float(np.mean([c[1] for c in cen])))
    return out


# ------------------------------------------------------------------ ZTF
def ztf_select(ph: pd.DataFrame, jd_lo: float, jd_hi: float, pad: float):
    """Which r-band frames represent the window, and how to describe them."""
    r = ph[ph["filter"] == "ZTF_r"].copy()
    if r.empty:
        return r, "no ZTF r-band frames"
    r["dt"] = np.where(r["obsjd"] < jd_lo, r["obsjd"] - jd_lo, np.where(r["obsjd"] > jd_hi, r["obsjd"] - jd_hi, 0.0))
    clean = r[r["quality_ok"].astype(bool)]
    near = clean[clean["dt"].abs() <= pad]
    if len(near):
        return near.reindex(near["dt"].abs().sort_values().index).head(12), f"within {pad:g} d of the window"
    usable = r[~r["flag_outside"].astype(bool)] if "flag_outside" in r else r
    near = usable[usable["dt"].abs() <= pad]
    if len(near):
        return near.reindex(near["dt"].abs().sort_values().index).head(12), f"within {pad:g} d, flagged frames"
    for limit, nmax in ((30.0, 5), (120.0, 3)):
        near = clean[clean["dt"].abs() <= limit]
        if len(near):
            near = near.reindex(near["dt"].abs().sort_values().index).head(nmax)
            return near, f"nearest clean frames, {near['dt'].median():+.0f} d from the window"
    return r.iloc[0:0], "no clean ZTF r-band frame within 120 d"


def ztf_stack(target: str, jd_lo: float, jd_hi: float, pad: float, radius: int):
    """The ZTF stack, the frames used, the ephemeris PAs and a note; (None, ...) without data."""
    if not HAVE_ZTF:
        return None, pd.DataFrame(), {}, f"ztfcomet unavailable ({_ZTF_ERR})"
    p = zc.photometry_path(target, create=False)
    if not p.exists():
        return None, pd.DataFrame(), {}, "not in the ZTF survey"
    ph = pd.read_csv(p, dtype=S)
    # a frame is usable when any aperture of it is clean: the 10 000 km rows of a distant
    # comet fail the scale test while the 20 000 km rows of the same frame are fine
    ph["quality_ok"] = ph.groupby("file")["quality_ok"].transform("any")
    ph = ph.drop_duplicates("file")
    sel, note = ztf_select(ph, jd_lo, jd_hi, pad)
    if sel.empty:
        return None, sel, {}, note
    ddir = zc.data_dir(target, create=False)
    size = 2 * radius + 1
    stamps, used, cen = [], [], []
    for r in sel.itertuples(index=False):
        f = ddir / r.file
        if not f.exists():
            continue
        try:
            with fits.open(f, memmap=False) as h:
                data = np.asarray(h[0].data, dtype=np.float64)
                hdr = h[0].header
        except Exception:                             # noqa: BLE001
            continue
        zp = float(hdr.get("MAGZP", np.nan))
        if not np.isfinite(zp):
            continue
        data = data * 10 ** (0.4 * (ZTF_AB_ZP - zp))    # DN -> uJy/pixel
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w_in = WCS(hdr)
        if bool(r.quality_ok) and np.isfinite(r.x_center) and np.isfinite(r.y_center):
            ra0, dec0 = (float(v) for v in w_in.all_pix2world(r.x_center, r.y_center, 0))
        else:
            ra0, dec0 = float(r.ra), float(r.dec)
        stamp = reproject(data, w_in, tan_wcs(ra0, dec0, ZTF_PIXSCALE, radius), size)
        stamps.append(stamp - ring_median(stamp))
        used.append(r.file)
        cen.append((ra0, dec0))
    if not stamps:
        return None, sel, {}, "ZTF frames not on disk"
    image, count = combine(stamps)
    st = Stack(image, count, ZTF_PIXSCALE, len(stamps), "uJy/pixel", used,
               float(np.mean([c[0] for c in cen])), float(np.mean([c[1] for c in cen])), note)
    pas = dict(sun_pa=float(sel["sunTargetPA"].median()) if "sunTargetPA" in sel else np.nan,
               vel_pa=float(sel["velocityPA"].median()) if "velocityPA" in sel else np.nan)
    return st, sel, pas, note


# ------------------------------------------------------------------ drawing
def _bar_km(fov_km: float) -> float:
    ok = [b for b in BAR_CHOICES_KM if b <= fov_km / 3.5]
    return ok[-1] if ok else BAR_CHOICES_KM[0]


def _show(ax, st: Optional[Stack], title: str, delta_au: float, ap_arcsec: float,
          extra_circle: Optional[tuple] = None, pas: Optional[dict] = None, empty_text: str = "",
          smooth_px: float = 0.0):
    ax.set_title(title, fontsize=15)
    if st is None:
        ax.text(0.5, 0.5, empty_text or "no exposure in this window", ha="center", va="center",
                transform=ax.transAxes, fontsize=16, color="0.4", wrap=True)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("0.7")
        return
    R = (st.image.shape[0] - 1) // 2
    half = (R + 0.5) * st.pixscale
    ext = (-half, half, -half, half)
    shown = st.image
    if smooth_px > 0:
        # display only: a light Gaussian brings a faint coma above the pixel noise
        good = np.isfinite(shown)
        num = ndimage.gaussian_filter(np.where(good, shown, 0.0), smooth_px)
        den = ndimage.gaussian_filter(good.astype(float), smooth_px)
        with np.errstate(invalid="ignore", divide="ignore"):
            shown = np.where(den > 0.5, num / den, np.nan)
    finite = shown[np.isfinite(shown)]
    if finite.size > 20:
        # the stretch follows the noise, not the brightest star: black at -1.5 sigma below the
        # sky, white at 12 sigma or the 99.9th percentile, whichever is brighter, asinh in between
        _, med, sig = sigma_clipped_stats(finite, sigma=3.0, maxiters=5)
        vmax = max(med + 12 * sig, float(np.percentile(finite, 99.9)))
        norm = ImageNormalize(vmin=med - 1.5 * sig, vmax=vmax, stretch=AsinhStretch(0.08))
    else:
        norm = None
    im = ax.imshow(shown, origin="lower", extent=ext, cmap="magma", norm=norm, interpolation="nearest")
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=11)
    cb.set_label(st.unit, fontsize=12)
    if np.isfinite(ap_arcsec):
        ax.add_patch(plt.Circle((0, 0), ap_arcsec, fill=False, ec="red", lw=1.8))
    if extra_circle is not None:
        rad, col, lab = extra_circle
        ax.add_patch(plt.Circle((0, 0), rad, fill=False, ec=col, lw=1.4, ls="--"))
    # scale bar, bottom left
    fov_km = 2 * half * KM_PER_ARCSEC_AU * delta_au
    bar = _bar_km(fov_km)
    bar_as = bar / (KM_PER_ARCSEC_AU * delta_au)
    x0, y0 = -half * 0.88, -half * 0.86
    ax.plot([x0, x0 + bar_as], [y0, y0], color="w", lw=3, solid_capstyle="butt")
    ax.text(x0 + bar_as / 2, y0 + half * 0.05, f"{bar:,.0f} km", color="w", ha="center", va="bottom", fontsize=11)
    # compass, top left: north up, east left after reprojection
    o = np.array([-half * 0.78, half * 0.66])
    L = half * 0.18
    style = dict(color="w", width=half * 0.008, head_width=half * 0.05, length_includes_head=True)
    ax.arrow(o[0], o[1], 0, L, **style); ax.text(o[0], o[1] + L * 1.25, "N", color="w", ha="center", va="bottom", fontsize=11)
    ax.arrow(o[0], o[1], -L, 0, **style); ax.text(o[0] - L * 1.25, o[1], "E", color="w", ha="right", va="center", fontsize=11)
    if pas:
        for key, lab, col in (("vel_pa", r"$-V$", "cyan"), ("sun_pa", r"$-\odot$", "yellow")):
            pa = pas.get(key, np.nan)
            if np.isfinite(pa):
                d = np.array([-np.sin(np.radians(pa)), np.cos(np.radians(pa))]) * L * 1.1
                ax.arrow(o[0], o[1], d[0], d[1], **{**style, "color": col})
                ax.text(*(o + d * 1.3), lab, color=col, ha="center", va="center", fontsize=11)
    ax.text(0.98, 0.03, f"N = {st.n}", transform=ax.transAxes, ha="right", va="bottom", color="w", fontsize=12)
    ax.set_xlabel(r"$\leftarrow$ E   offset ($''$)   W $\rightarrow$", fontsize=13)
    ax.set_ylabel(r"offset ($''$)   N $\uparrow$", fontsize=13)
    ax.tick_params(labelsize=11)


def _q_line(row: pd.Series, key: str, pretty: str) -> str:
    st = str(row.get(f"Q_{key}_status", ""))
    if st in ("detected", "marginal"):
        return f"Q({pretty}) = {row[f'Q_{key}']:.2e} ± {row[f'Q_{key}_err']:.1e} s$^{{-1}}$" + (" (marginal)" if st == "marginal" else "")
    lim = row.get(f"Q_{key}_upper_limit", np.nan)
    return f"Q({pretty}) < {lim:.1e} s$^{{-1}}$" if np.isfinite(lim) else f"Q({pretty}): {st.replace('_', ' ')}"


def _afrho_line(row: pd.Series, tag: str) -> str:
    v, e, m = row.get(f"afrho_{tag}_cm", np.nan), row.get(f"afrho_{tag}_err_cm", np.nan), str(row.get(f"afrho_{tag}_method", ""))
    if np.isfinite(v):
        return f"{tag}: {v:.3g} ± {e:.2g} cm ({m})"
    return f"{tag}: —"


def draw_figure(stem: str, row: pd.Series, info: dict, ztf: Optional[Stack], ztf_sel: pd.DataFrame, pas: dict,
                ztf_note: str, sx: Dict[str, Optional[Stack]], out: Path, dpi: int = 80):
    fig, axes = plt.subplots(2, 3, figsize=(19, 11.6))
    delta = info["delta"]
    ap_as = float(row["r_ap_km"]) / (KM_PER_ARCSEC_AU * delta)
    ztf_title = f"ZTF r-band   N = {ztf.n}   (shown smoothed 1$''$)\n{ztf_note}" if ztf is not None else "ZTF r-band"
    _show(axes[0, 0], ztf, ztf_title, delta, ap_as, extra_circle=(1e4 / (KM_PER_ARCSEC_AU * delta), "cyan", "10,000 km"),
          pas=pas, empty_text=ztf_note, smooth_px=1.0)
    order = [("cont", axes[0, 1]), ("H2O", axes[0, 2]), ("CO2", axes[1, 0]), ("CO", axes[1, 1])]
    for key, ax in order:
        label, (lo, hi) = WINDOWS[key]
        st = sx.get(key)
        n_win = info["n_window"].get(key, 0)
        _show(ax, st, f"SPHEREx {label}  {lo:.2f}–{hi:.2f} µm" + (f"   N = {st.n}" if st else ""), delta, ap_as, pas=pas,
              empty_text=("no exposure of this phase in the window" if n_win == 0 else "exposures not readable"))
    # info panel
    ax = axes[1, 2]; ax.axis("off")
    t0, t1 = Time(info["jd_lo"], format="jd").iso[:10], Time(info["jd_hi"], format="jd").iso[:10]
    lines = [f"{row['target']}  phase {int(row['phase'])}:  {t0} → {t1}",
             f"⟨r_h⟩ = {row['r_hel_mean']:.2f} au,  Δ = {delta:.2f} au,  α = {info['alpha']:.1f}°",
             f"fit aperture {float(row['r_ap_km']):,.0f} km = {ap_as:.0f}\" = {ap_as / SPHEREX_PIXSCALE:.1f} SPHEREx px (red)",
             f"SPHEREx exposures in the phase: {info['n_phase']}",
             "  in the windows: " + ", ".join(f"{k} {info['n_window'].get(k, 0)}" for k in WINDOWS),
             ""]
    if ztf is not None:
        d = ztf_sel.sort_values("obsjd")
        lines += [f"ZTF r frames: {ztf.n}, {Time(d['obsjd'].min(), format='jd').iso[:10]} → "
                  f"{Time(d['obsjd'].max(), format='jd').iso[:10]}",
                  f"  {ztf_note}", "  cyan dashed: the 10,000 km aperture of the Afρ series"]
    else:
        lines += [f"ZTF: {ztf_note}"]
    lines += ["", "A(0°)fρ at ⟨r_h⟩ from ZTF:", "  " + _afrho_line(row, "10k"), "  " + _afrho_line(row, "20k"),
              "", _q_line(row, "H2O", "H$_2$O"), _q_line(row, "CO2", "CO$_2$"), _q_line(row, "CO", "CO"), "",
              "stacks: sigma-clipped median, background-subtracted,",
              "north up / east left, comet-centred on each frame's",
              "ephemeris; 310\" field at the native pixel scale"]
    ax.text(0.0, 0.98, "\n".join(lines), transform=ax.transAxes, va="top", ha="left", fontsize=12.5,
            linespacing=1.4)
    fig.suptitle(f"{row['target']}  phase {int(row['phase'])}  —  ZTF and SPHEREx stacks at the phase", fontsize=20, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------ driver
def aperture_label(r_km: float) -> str:
    e = int(np.floor(np.log10(float(r_km))))
    return f"{float(r_km) / 10 ** e:g}e{e}"


def process_group(args) -> str:
    row, pad, radius_sx, write_fits = args
    row = pd.Series(row)
    target, phase = str(row["target"]), int(row["phase"])
    stem = f"{target}_{aperture_label(row['r_ap_km'])}km_ph{phase}"
    assignment = _assignment()
    grp = assignment[(assignment["target"] == target) & (assignment["phase"] == phase)]
    ap = apphot_rows(target)
    apg = ap[ap["filename"].isin(set(grp["filename"]))]
    jd_lo, jd_hi = float(grp["jd_utc"].min()), float(grp["jd_utc"].max())
    info = dict(jd_lo=jd_lo, jd_hi=jd_hi, n_phase=len(grp),
                delta=float(apg["r_obs"].median()) if len(apg) else float(row.get("r_obs_mean", np.nan)),
                alpha=float(apg["alpha"].median()) if len(apg) else np.nan,
                n_window={k: int(((apg["wl"] >= lo) & (apg["wl"] <= hi)).sum()) for k, (_, (lo, hi)) in WINDOWS.items()})
    resolver = FitsResolver(ad.FITS_ROOTS)
    sx = spherex_stacks(target, phase, assignment, resolver, ApConfig(), radius_sx)
    radius_ztf = int(round(FOV_RADIUS_ARCSEC / ZTF_PIXSCALE))
    ztf, sel, pas, note = ztf_stack(target, jd_lo, jd_hi, pad, radius_ztf)
    if write_fits:
        for key, st in sx.items():
            if st is not None:
                write_stack(STACK_DIR / f"{stem}_{key}.fits", st, target, phase, key, WINDOWS[key][1])
        if ztf is not None:
            write_stack(STACK_DIR / f"{stem}_ZTFr.fits", ztf, target, phase, "ZTF_r", None)
    draw_figure(stem, row, info, ztf, sel, pas, note, sx, FIG_DIR / f"{stem}.png")
    made = [k for k, v in sx.items() if v is not None]
    return f"{stem}: ZTF {ztf.n if ztf else 0} ({note}); SPHEREx " + ", ".join(f"{k}={sx[k].n}" for k in made)


_ASSIGN: Optional[pd.DataFrame] = None


def _assignment() -> pd.DataFrame:
    global _ASSIGN
    if _ASSIGN is None:
        _ASSIGN = pd.read_csv(cd.DATA_DIR / "phase_assignment.csv", dtype=S)
    return _ASSIGN


def select(df: pd.DataFrame, args: List[str]) -> pd.DataFrame:
    if not args:
        return df
    keep = pd.Series(False, index=df.index)
    for a in args:
        if ":" in a:
            t, ph = a.split(":")
            keep |= (df["target"] == t) & (df["phase"].astype(int) == int(ph))
        else:
            keep |= df["target"] == a
    return df[keep]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("groups", nargs="*", help="targets or target:phase (default: every group of gas_fit.csv)")
    ap.add_argument("--all-phases", action="store_true", help="every phase of phase_map.csv, not only the fitted ones")
    ap.add_argument("--pad", type=float, default=5.0, help="days around the SPHEREx window for the ZTF frames")
    ap.add_argument("--radius", type=int, default=25, help="SPHEREx stamp half-width [pixel]; the ZTF stamp matches it on the sky")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--no-fits", action="store_true", help="figures only, no stack FITS")
    args = ap.parse_args(argv)
    global FOV_RADIUS_ARCSEC
    FOV_RADIUS_ARCSEC = args.radius * SPHEREX_PIXSCALE
    if args.all_phases:
        pm = pd.read_csv(cd.RESULT_DIR / "phase_map.csv", dtype=S)
        apt = pd.read_csv(cd.RESULT_DIR / "apertures.csv", dtype=S)[["target", "phase", "r_ap_km"]]
        df = pm.merge(apt, on=["target", "phase"], how="left")
        df["r_ap_km"] = df["r_ap_km"].fillna(20000.0)
    else:
        df = pd.read_csv(cd.RESULT_DIR / "gas_fit.csv", dtype=S)
    sel = select(df, args.groups)
    if sel.empty:
        print(f"no group matched {args.groups}")
        return 1
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    jobs = [(r.to_dict(), args.pad, args.radius, not args.no_fits) for _, r in sel.iterrows()]
    if args.workers > 1:
        with ProcessPoolExecutor(args.workers) as ex:
            for line in ex.map(process_group, jobs):
                print(line, flush=True)
    else:
        for j in jobs:
            print(process_group(j), flush=True)
    print(f"{len(jobs)} figure(s) -> {FIG_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
