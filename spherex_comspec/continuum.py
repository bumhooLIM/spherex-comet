"""
Local polynomial continuum subtraction under each emission band.

Reproduces ``notebooks/continuum_subtraction.ipynb``.  For every band a weighted
polynomial in ``lam - lam_c`` is fitted to the continuum window with *all three*
emission windows punched out, sigma-clipped, and validated on four checks
(bracketing and positivity are hard; residual shape and cross-validation are
soft).  Its order is capped at 1 whenever the continuum sample does not bracket
the band: curvature you have to extrapolate is curvature you cannot justify.

Flux space
----------
The spectrum handed in carries ``flux``/``err`` in the variant's flux space --
physical mJy, or distance-corrected ``F x r_h^2 Delta^2`` -- and every fitted
and subtracted quantity is in that space.  ``emis_raw_mjy`` is then the
emission divided back by each channel's own ``distcorr_factor``: the physical
flux the fitter needs.  When the variant fits in physical space the two are
identical.

Why fit the continuum in distance-corrected space at all: SPHEREx scans a
moving target non-simultaneously, so the channels of one group were taken at
different geometries.  Within a 28-day epoch of 24P the factor r_h^2 Delta^2
varies by 3x; the raw spectrum therefore carries a geometric gradient across
wavelength that a polynomial continuum would partly absorb into the band.  The
corrected spectrum has every channel on the same footing.
"""

from __future__ import annotations

import warnings
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from astropy.stats import mad_std, sigma_clip

from .config import ALL_EM_WINDOWS, BAND_WINDOWS, MJY_TO_WM2UM, ContinuumConfig
from .dataio import slug
from .logging_utils import get_logger

__all__ = ["continuum_mask", "fit_continuum", "continuum_at", "cv_rmse", "select_order", "validate_fit",
           "subtract_continuum", "aggregate_band", "insufficient", "n_in_emission",
           "process_group", "continuum_model_columns", "rebuild_continuum"]

log = get_logger("continuum")


# ------------------------------------------------------------------------------- fitting
def continuum_mask(wl, cont_win, exclude=ALL_EM_WINDOWS) -> np.ndarray:
    """Points inside the continuum window and outside every emission window."""
    m = (wl >= cont_win[0]) & (wl <= cont_win[1])
    for lo, hi in exclude:
        m &= ~((wl >= lo) & (wl <= hi))
    return m


def _wpolyfit(x, y, err, order):
    """Weighted polynomial fit -> (coeffs, covariance), scaled covariance when n > order + 2."""
    cov_mode = True if len(x) > order + 2 else "unscaled"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c, cov = np.polyfit(x, y, order, w=1.0 / err, cov=cov_mode)
    return c, cov


def _poly_with_sigma(c, cov, x):
    model = np.polyval(c, x)
    V = np.vander(np.atleast_1d(x), len(c))
    sig = np.sqrt(np.einsum("ij,jk,ik->i", V, cov, V).clip(min=0))
    return model, sig


def fit_continuum(raw: pd.DataFrame, band: str, cfg: ContinuumConfig, windows=None,
                  rev=None) -> dict:
    """
    Window, sigma-clip and fit the local continuum for one band.

    Returns a dict with the fitted model and every piece of bookkeeping the
    validation, the subtraction and the figures need.

    ``windows`` replaces ``BAND_WINDOWS`` (a case revision's effective windows, see
    :func:`revisions.effective_windows`); ``rev`` is the band's :class:`revisions.BandRevision`:
    a fixed order, the exclusion of the brightest continuum points on a side, and the
    acceptance of a one-sided window as given (no extension).
    """
    W = windows if windows is not None else BAND_WINDOWS
    b = W[band]
    exclude = [w["em"] for w in W.values()]
    lam_c = b["lam_c"]
    fixed_order = rev is not None and rev.order is not None
    order_req = int(rev.order) if fixed_order else int(cfg.poly_orders.get(band, 2))
    wl, y, e = raw.wl.to_numpy(float), raw.flux.to_numpy(float), raw.err.to_numpy(float)
    cont = tuple(b["cont"])
    cmask = continuum_mask(wl, cont, exclude)
    extended = ""
    one_sided_ok = rev is not None and rev.one_sided_ok
    if cfg.one_sided_extend_um is not None and cmask.sum() > 0 and not one_sided_ok:
        # a window with points on one side only: reach out to `one_sided_extend_um` from the
        # band edge on the empty side (every emission window is still punched out)
        n_blue0 = int((wl[cmask] < b["em"][0]).sum())
        n_red0 = int((wl[cmask] > b["em"][1]).sum())
        lo, hi = cont
        if n_blue0 == 0 and n_red0 > 0:
            lo, extended = b["em"][0] - cfg.one_sided_extend_um, "blue"
        elif n_red0 == 0 and n_blue0 > 0:
            hi, extended = b["em"][1] + cfg.one_sided_extend_um, "red"
        if extended:
            cont = (lo, hi)
            cmask = continuum_mask(wl, cont, exclude)
    # case revision: the N brightest points of a side of the continuum window are
    # star-contaminated and are taken out before the clipping, which would not catch them
    # against a sparse baseline.  The reviewer counted them on the validation figure, where
    # a side shows every point of the window -- those inside a neighbouring band's emission
    # window included (the red side of the 4.3 um window overlaps the CO band) -- so the N
    # brightest are ranked over all of them, and the continuum candidates among them are dropped.
    n_excl = {"blue": 0, "red": 0}
    if rev is not None and (rev.exclude_top_blue or rev.exclude_top_red):
        in_win = (wl >= cont[0]) & (wl <= cont[1])
        for side, n, sel in (("blue", rev.exclude_top_blue, wl < b["em"][0]),
                             ("red", rev.exclude_top_red, wl > b["em"][1])):
            if n:
                idx = np.where(in_win & sel)[0]
                top = idx[np.argsort(y[idx])[::-1][:int(n)]]
                drop = top[cmask[top]]
                cmask[drop] = False
                n_excl[side] = int(len(drop))
    xc, yc, ec = wl[cmask], y[cmask], e[cmask]

    res = dict(band=band, lam_c=lam_c, em=b["em"], cont=cont, order_req=order_req,
               n_cont_avail=int(cmask.sum()), extended=extended,
               n_blue=int((xc < b["em"][0]).sum()), n_red=int((xc > b["em"][1]).sum()),
               cmask=cmask, xc=xc, yc=yc, ec=ec, order_fixed=fixed_order,
               n_excl_blue=n_excl["blue"], n_excl_red=n_excl["red"],
               revision=rev.describe() if rev is not None else "")

    if len(xc) < 2:
        res.update(method="none", ok=False, reason="fewer than 2 continuum points",
                   keep=np.zeros(len(xc), bool), order_used=np.nan, coeffs=None, cov=None,
                   order_capped=False, bracketed=False, n_cont_used=0, n_clipped=0)
        return res

    # a fixed order (case revision) is fitted on as few as order + 2 points; the default keeps
    # the two-point fallback below ``n_min_poly``
    n_min_poly = max(order_req + 2, 3) if fixed_order else cfg.n_min_poly
    if len(xc) < n_min_poly:
        left, right = xc < b["em"][0], xc > b["em"][1]
        if left.any() and right.any():
            i = np.where(left)[0][np.argmax(xc[left])]
            j = np.where(right)[0][np.argmin(xc[right])]
            sel, bracketed = np.array([i, j]), True
        else:
            sel, bracketed = np.argsort(np.abs(xc - lam_c))[:2], False
        keep = np.zeros(len(xc), bool)
        keep[sel] = True
        c, cov = _wpolyfit(xc[sel] - lam_c, yc[sel], ec[sel], 1)
        res.update(method="2point", ok=True, reason="", keep=keep, order_used=1, coeffs=c,
                   cov=cov, bracketed=bracketed, order_capped=False, n_cont_used=2, n_clipped=0)
        return res

    has_left = bool((xc < b["em"][0]).any())
    has_right = bool((xc > b["em"][1]).any())
    order_cap = order_req if (has_left and has_right) else 1
    res["order_capped"] = order_cap != order_req
    order_cap = int(np.clip(order_cap, 1, max(1, len(xc) - (2 if fixed_order else 3))))
    if cfg.order_mode == "cv" and order_cap > 1 and not fixed_order:
        order_cap = select_order(xc - lam_c, yc, ec, order_cap, cfg)
    res["order_selected"] = order_cap
    order = order_cap
    keep = np.ones(len(xc), bool)
    for _ in range(cfg.maxiters):
        c, cov = _wpolyfit(xc[keep] - lam_c, yc[keep], ec[keep], order)
        chi = (yc - np.polyval(c, xc - lam_c)) / ec
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            new = ~np.asarray(sigma_clip(chi, sigma=cfg.sigma, maxiters=1, stdfunc=mad_std,
                                         masked=True).mask, bool)
        if new.sum() < order + (2 if fixed_order else 3) or (new == keep).all():
            break
        keep = new
        order = int(np.clip(order_cap, 1, max(1, keep.sum() - (2 if fixed_order else 3))))
    c, cov = _wpolyfit(xc[keep] - lam_c, yc[keep], ec[keep], order)
    res.update(method="poly", ok=True, reason="", keep=keep, order_used=order, coeffs=c, cov=cov,
               bracketed=bool((xc[keep] < b["em"][0]).any() and (xc[keep] > b["em"][1]).any()),
               n_cont_used=int(keep.sum()), n_clipped=int((~keep).sum()))
    return res


def continuum_at(fit: dict, lam):
    """Continuum model and its 1-sigma uncertainty at ``lam``."""
    lam = np.atleast_1d(np.asarray(lam, float))
    if fit["coeffs"] is None:
        return np.full(lam.shape, np.nan), np.full(lam.shape, np.nan)
    return _poly_with_sigma(fit["coeffs"], fit["cov"], lam - fit["lam_c"])


# ---------------------------------------------------------------------------- validation
def select_order(x, y, err, max_order: int, cfg: ContinuumConfig) -> int:
    """
    Polynomial order by leave-one-out cross-validation: the lowest order in 1..max_order whose
    CV RMSE is within ``cfg.cv_select_margin`` of the best (parsimony breaks near-ties).  Falls
    back to 1 when no order can be cross-validated.
    """
    cand = {o: cv_rmse(x, y, err, o, cfg) for o in range(1, int(max_order) + 1)}
    cand = {o: r for o, r in cand.items() if np.isfinite(r)}
    if not cand:
        return 1
    best = min(cand.values())
    for o in sorted(cand):
        if cand[o] <= (1.0 + cfg.cv_select_margin) * best:
            return int(o)
    return int(min(cand, key=cand.get))


def cv_rmse(x, y, err, order, cfg: ContinuumConfig, seed: int = 0) -> float:
    """Cross-validated RMSE (leave-one-out, or 10-fold above ``cfg.cv_loo_max`` points)."""
    n = len(x)
    if n < order + 3:
        return np.nan
    if n <= cfg.cv_loo_max:
        folds = [np.array([i]) for i in range(n)]
    else:
        folds = np.array_split(np.random.default_rng(seed).permutation(n), 10)
    sq = []
    for f in folds:
        m = np.ones(n, bool)
        m[f] = False
        if m.sum() < order + 2:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c = np.polyfit(x[m], y[m], order, w=1.0 / err[m])
        sq.extend((y[f] - np.polyval(c, x[f])) ** 2)
    return float(np.sqrt(np.mean(sq))) if sq else np.nan


def validate_fit(fit: dict, cfg: ContinuumConfig, rev=None) -> dict:
    """Automatic quality checks -> PASS / WARN / FAIL with the reasons.

    A case revision (``rev``) can waive the two hard checks: ``one_sided_ok`` accepts an
    unbracketed continuum, ``ignore_negative`` a continuum below zero under the band; both
    are recorded in ``notes``."""
    v = dict(band=fit["band"], method=fit["method"])
    blank = dict(chi2_red=np.nan, rms_mjy=np.nan, err_scale=np.nan, z_shape=np.nan,
                 cv_rmse=np.nan, cv_best_order=np.nan, cont_min=np.nan, cont_min_nsig=np.nan,
                 pass_shape=False, pass_cv=False, pass_bracket=False, pass_positive=False)
    if not fit["ok"]:
        note = fit["reason"] + ("; case revision: " + fit["revision"] if fit.get("revision") else "")
        return {**v, **blank, "verdict": "FAIL", "notes": note}

    keep = fit["keep"]
    x = fit["xc"][keep] - fit["lam_c"]
    y, e = fit["yc"][keep], fit["ec"][keep]
    resid = y - np.polyval(fit["coeffs"], x)
    dof = len(x) - (fit["order_used"] + 1)
    v["chi2_red"] = float(np.sum((resid / e) ** 2) / dof) if dof > 0 else np.nan
    v["rms_mjy"] = float(np.sqrt(np.mean(resid ** 2)))
    v["err_scale"] = float(v["rms_mjy"] / np.median(e)) if dof > 0 else np.nan

    if dof > 0 and len(x) >= 4:
        half = x < np.median(x)
        zs = [abs(resid[m].mean()) / (v["rms_mjy"] / np.sqrt(m.sum()) + 1e-12)
              for m in (half, ~half) if m.sum() >= 2]
        v["z_shape"] = float(max(zs)) if zs else np.nan
    else:
        v["z_shape"] = np.nan

    cand = {o: cv_rmse(x, y, e, o, cfg) for o in (1, 2, 3)}
    cand = {o: r for o, r in cand.items() if np.isfinite(r)}
    v["cv_rmse"] = cand.get(fit["order_used"], np.nan)
    v["cv_best_order"] = min(cand, key=cand.get) if cand else np.nan
    best = min(cand.values()) if cand else np.nan

    lam_em = np.linspace(*fit["em"], 60)
    cont_em, cont_sig = continuum_at(fit, lam_em)
    imin = int(np.argmin(cont_em))
    v["cont_min"] = float(cont_em[imin])
    if cont_em[imin] >= 0:
        v["cont_min_nsig"] = 0.0
    elif np.isfinite(cont_sig[imin]) and cont_sig[imin] > 0:
        v["cont_min_nsig"] = float(-cont_em[imin] / cont_sig[imin])
    else:
        v["cont_min_nsig"] = np.inf

    v["pass_shape"] = bool(not np.isfinite(v["z_shape"]) or v["z_shape"] <= cfg.z_shape_max)
    v["pass_cv"] = bool(not np.isfinite(v["cv_rmse"]) or not np.isfinite(best)
                        or v["cv_rmse"] <= cfg.cv_tol * best)
    v["pass_bracket"] = bool(fit.get("bracketed", False))
    v["pass_positive"] = bool(v["cont_min_nsig"] <= cfg.pos_nsig_max)

    notes = []
    if not v["pass_bracket"]:
        if rev is not None and rev.one_sided_ok:
            v["pass_bracket"] = True
            notes.append("one-sided continuum accepted as fitted (case revision)")
        else:
            notes.append("EXTRAPOLATED (no continuum on one side)")
    if not v["pass_positive"]:
        if rev is not None and rev.ignore_negative:
            v["pass_positive"] = True
            notes.append(f"continuum negative by {v['cont_min_nsig']:.1f} sigma under the band "
                         "-- accepted (case revision)")
        else:
            notes.append(f"continuum negative by {v['cont_min_nsig']:.1f} sigma under the band")
    if not v["pass_shape"]:
        notes.append(f"half-window residual offset {v['z_shape']:.1f} sigma")
    if not v["pass_cv"]:
        notes.append(f"CV favours order {v['cv_best_order']:.0f}")
    if fit["method"] == "2point":
        notes.append("2-point fallback (sparse continuum, 0 dof -- unverifiable)")
    if fit.get("extended"):
        notes.append(f"continuum window extended on the {fit['extended']} side to "
                     f"{fit['cont'][0]:.2f}-{fit['cont'][1]:.2f} um")
    if fit.get("order_capped"):
        notes.append(f"order capped {fit['order_req']}->1 (one-sided continuum)")
    elif fit.get("order_fixed") and fit["method"] == "poly":
        notes.append(f"order {fit['order_used']} fixed (case revision)")
    elif fit["order_used"] != fit["order_req"] and fit["method"] == "poly":
        notes.append(f"order {fit['order_used']} of max {fit['order_req']} by CV"
                     if cfg.order_mode == "cv" else
                     f"order reduced {fit['order_req']}->{fit['order_used']}")
    if np.isfinite(v["err_scale"]) and v["err_scale"] > 3:
        notes.append(f"formal errors understate scatter by x{v['err_scale']:.0f}")
    if fit.get("revision"):
        notes.append("case revision: " + fit["revision"])

    hard = v["pass_bracket"] and v["pass_positive"]
    soft = v["pass_shape"] and v["pass_cv"]
    v["verdict"] = "PASS" if (hard and soft) else ("WARN" if hard else "FAIL")
    if v["verdict"] == "PASS" and fit["method"] == "2point":
        v["verdict"] = "WARN"
    v["notes"] = "; ".join(notes)
    return v


# --------------------------------------------------------------------------- subtraction
def subtract_continuum(raw: pd.DataFrame, fit: dict) -> pd.DataFrame:
    """
    Continuum-subtract every point of the spectrum and label its ``role``.

    Adds, in the fit's flux space, ``cont_mjy``, ``cont_err_mjy``, ``emis_mjy``,
    ``emis_err_mjy`` and ``chi``; and in physical space ``cont_raw_mjy``,
    ``emis_raw_mjy``, ``emis_raw_err_mjy`` -- the same quantities divided by
    each channel's ``distcorr_factor``, which is 1 when the space is physical.
    """
    lam = raw.wl.to_numpy(float)
    cont, cont_err = continuum_at(fit, lam)
    out = raw.copy()
    out["cont_mjy"] = cont
    out["cont_err_mjy"] = cont_err
    out["emis_mjy"] = out.flux - cont
    out["emis_err_mjy"] = np.sqrt(out.err ** 2 + cont_err ** 2)
    out["chi"] = (out.flux - cont) / out.err
    f = out["distcorr_factor"].to_numpy(float) if "distcorr_factor" in out else np.ones(len(out))
    space_factor = f if _is_distcorr(raw) else np.ones(len(out))
    out["cont_raw_mjy"] = out.cont_mjy / space_factor
    out["emis_raw_mjy"] = out.emis_mjy / space_factor
    out["emis_raw_err_mjy"] = out.emis_err_mjy / space_factor
    out["in_emission"] = (lam >= fit["em"][0]) & (lam <= fit["em"][1])
    role = np.full(len(out), "window", dtype=object)
    role[out.in_emission.to_numpy()] = "emission"
    idx = np.where(fit["cmask"])[0]
    if len(idx) and fit["keep"].size == len(idx):
        role[idx[fit["keep"]]] = "cont_used"
        role[idx[~fit["keep"]]] = "cont_clipped"
    out["role"] = role
    out["band"] = fit["band"]
    return out


def _is_distcorr(raw: pd.DataFrame) -> bool:
    return bool(raw.attrs.get("flux_space", "physical") == "distcorr")


def _trapz_band_flux(lam, f_lam, f_lam_err, cont_err_lam):
    """Trapezoid integral of F_lambda with point noise in quadrature and the continuum term linear."""
    if len(lam) < 2:
        return np.nan, np.nan
    o = np.argsort(lam)
    lam, f_lam, f_lam_err, cont_err_lam = lam[o], f_lam[o], f_lam_err[o], cont_err_lam[o]
    w = np.zeros_like(lam)
    dl = np.diff(lam)
    w[:-1] += dl / 2
    w[1:] += dl / 2
    total = float(np.sum(w * f_lam))
    err_pts = float(np.sqrt(np.sum((w * f_lam_err) ** 2)))
    err_cont = float(np.sum(w * cont_err_lam))
    return total, float(np.hypot(err_pts, err_cont))


def fail_reason(row: dict) -> str:
    if row.get("verdict") != "FAIL":
        return ""
    r = []
    if row.get("fit_method") == "none":
        r.append("no_continuum")
    if row.get("n_emission", 0) == 0:
        r.append("no_emission_data")
    if not row.get("pass_bracket", True):
        r.append("extrapolated")
    if not row.get("pass_positive", True):
        r.append("negative_continuum")
    return "+".join(r) if r else "other"


def continuum_model_columns(fit: dict, max_order: int = 3) -> dict:
    """Coefficients and covariance in an increasing-power basis (numpy's is decreasing)."""
    out = {"cont_lam_ref_um": fit["lam_c"]}
    n = max_order + 1
    c = np.full(n, np.nan)
    cov = np.full((n, n), np.nan)
    if fit.get("coeffs") is not None:
        k = len(fit["coeffs"])
        c[:k] = fit["coeffs"][::-1]
        cov[:k, :k] = np.asarray(fit["cov"])[::-1, ::-1]
    out.update({f"cont_c{i}": c[i] for i in range(n)})
    out.update({f"cont_cov_{i}{j}": cov[i, j] for i in range(n) for j in range(i, n)})
    return out


def rebuild_continuum(row, lam, max_order: int = 3):
    """Evaluate a saved continuum model and its 1-sigma band at ``lam`` from CSV columns alone."""
    dl = np.asarray(lam, float) - row["cont_lam_ref_um"]
    c = np.nan_to_num(np.array([row[f"cont_c{i}"] for i in range(max_order + 1)]))
    C = np.zeros((max_order + 1,) * 2)
    for i in range(max_order + 1):
        for j in range(i, max_order + 1):
            v = row[f"cont_cov_{i}{j}"]
            C[i, j] = C[j, i] = 0.0 if not np.isfinite(v) else v
    V = np.stack([dl ** i for i in range(max_order + 1)], axis=1)
    return V @ c, np.sqrt(np.clip(np.einsum("ij,jk,ik->i", V, C, V), 0, None))


def aggregate_band(sub: pd.DataFrame, fit: dict, val: dict, target: str, r_ap_km: float,
                   phase: int, epoch: int, flux_space: str, cfg: ContinuumConfig) -> dict:
    """Per-band summary row: fit provenance, validation, the emission detection, band flux."""
    em = sub[sub.in_emission]
    row = dict(
        target=slug(target), r_ap_km=float(r_ap_km), phase=int(phase), epoch=int(epoch),
        band=fit["band"], flux_space=flux_space,
        lam_center_um=fit["lam_c"], em_lo_um=fit["em"][0], em_hi_um=fit["em"][1],
        cont_lo_um=fit["cont"][0], cont_hi_um=fit["cont"][1],
        fit_method=fit["method"], poly_order_req=fit["order_req"],
        poly_order_used=fit["order_used"],
        n_cont_avail=fit["n_cont_avail"], n_cont_blue=fit["n_blue"], n_cont_red=fit["n_red"],
        n_cont_used=fit["n_cont_used"], n_cont_clipped=fit["n_clipped"],
        bracketed=bool(fit.get("bracketed", False)), order_capped=bool(fit.get("order_capped", False)),
        cont_extended=str(fit.get("extended", "")),
        chi2_red=val["chi2_red"], rms_cont_mjy=val["rms_mjy"], err_scale=val["err_scale"],
        z_shape=val["z_shape"], cv_rmse_mjy=val["cv_rmse"], cv_best_order=val["cv_best_order"],
        pass_bracket=val["pass_bracket"], pass_positive=val["pass_positive"],
        pass_shape=val["pass_shape"], pass_cv=val["pass_cv"],
        cont_min_em_mjy=val["cont_min"], cont_min_nsig=val["cont_min_nsig"],
        **continuum_model_columns(fit, cfg.max_order),
        verdict=val["verdict"], notes=val["notes"], n_emission=len(em),
        n_flag_a=int((em.sourceflag.astype(str) == "a").sum()) if "sourceflag" in em else 0,
        n_flag_b=int((em.sourceflag.astype(str) == "b").sum()) if "sourceflag" in em else 0,
    )
    usable = len(em) > 0 and bool(np.isfinite(em.emis_mjy.to_numpy()).any())
    keys = ("wl_peak_um", "f_peak_mjy", "f_peak_err_mjy", "snr_peak", "r_hel_min", "r_hel_mean",
            "r_hel_max", "r_hel_peak", "r_obs_mean", "jd_utc_mean", "band_flux_W_m2",
            "band_flux_err_W_m2", "cont_at_center_mjy", "cont_err_at_center_mjy",
            "distcorr_factor_mean")
    if not usable:
        row.update({k: np.nan for k in keys})
        why = "no data in emission window" if len(em) == 0 else "no continuum fit, emission undefined"
        row["notes"] = (row["notes"] + "; " if row["notes"] else "") + why
        row["verdict"] = "FAIL"
        row["fail_reason"] = fail_reason(row)
        return row

    ipk = int(np.nanargmax(em.emis_mjy.to_numpy()))
    pk = em.iloc[ipk]
    # the band flux is a physical quantity, so it is integrated over the physical emission
    lam = em.wl.to_numpy(float)
    f_lam = em.emis_raw_mjy.to_numpy(float) * MJY_TO_WM2UM / lam ** 2
    f_lam_e = em.emis_raw_err_mjy.to_numpy(float) * MJY_TO_WM2UM / lam ** 2
    fac = em.distcorr_factor.to_numpy(float) if flux_space == "distcorr" else np.ones(len(em))
    c_lam_e = (em.cont_err_mjy.to_numpy(float) / fac) * MJY_TO_WM2UM / lam ** 2
    bf, bfe = _trapz_band_flux(lam, f_lam, f_lam_e, c_lam_e)
    c0, c0e = continuum_at(fit, fit["lam_c"])
    row.update(wl_peak_um=float(pk.wl), f_peak_mjy=float(pk.emis_mjy),
               f_peak_err_mjy=float(pk.emis_err_mjy),
               snr_peak=float(pk.emis_mjy / pk.emis_err_mjy),
               r_hel_min=float(em.r_hel.min()), r_hel_mean=float(em.r_hel.mean()),
               r_hel_max=float(em.r_hel.max()), r_hel_peak=float(pk.r_hel),
               r_obs_mean=float(em.r_obs.mean()), jd_utc_mean=float(em.jd_utc.mean()),
               band_flux_W_m2=bf, band_flux_err_W_m2=bfe,
               cont_at_center_mjy=float(c0[0]), cont_err_at_center_mjy=float(c0e[0]),
               distcorr_factor_mean=float(em.distcorr_factor.mean()))
    row["fail_reason"] = fail_reason(row)
    return row


# ---------------------------------------------------------------------------- one group
def n_in_emission(raw: pd.DataFrame) -> int:
    m = np.zeros(len(raw), bool)
    for lo, hi in ALL_EM_WINDOWS:
        m |= (raw.wl >= lo) & (raw.wl <= hi)
    return int(m.sum())


def insufficient(raw: pd.DataFrame, cfg: ContinuumConfig) -> Optional[str]:
    """Why a group cannot be analysed at all, or ``None``.  A skipped group is not a FAIL row."""
    if len(raw) < cfg.min_points_phase:
        return f"{len(raw)} usable points (< {cfg.min_points_phase})"
    n_em = n_in_emission(raw)
    if n_em < cfg.min_points_emission:
        return f"{n_em} points inside the emission windows (< {cfg.min_points_emission})"
    return None


def process_group(raw: pd.DataFrame, target: str, r_ap_km: float, phase: int, epoch: int,
                  cfg: ContinuumConfig, flux_space: str = "physical", rev=None) -> Dict[str, object]:
    """
    Steps 2-6 for one (target, aperture, phase): fit, validate and subtract every band.

    ``rev`` is the group's :class:`revisions.CaseRevision` (``None`` for none): its windows
    replace the module constants for every band of the group, its per-band directives are
    passed to :func:`fit_continuum` / :func:`validate_fit`, rows carrying a dropped source
    flag leave the band before the fit, and the brightest emission channels it excludes are
    labelled ``role = "excluded"`` (out of the production-rate fit and the band flux).

    Returns
    -------
    dict
        ``fits``, ``valid``, ``subs`` keyed by band, ``summary`` (one row per band)
        and ``points`` (the plotted neighbourhood of every band, for the file).
    """
    from .revisions import effective_windows
    raw = raw.copy()
    raw.attrs["flux_space"] = flux_space
    W = effective_windows(rev)
    fits, vals, subs, revs = {}, {}, {}, {}
    for b in BAND_WINDOWS:
        br = rev.band(b) if rev is not None else None
        revs[b] = br
        raw_b = raw
        if br is not None and br.drop_flags and "sourceflag" in raw:
            raw_b = raw[~raw.sourceflag.astype(str).isin(br.drop_flags)].reset_index(drop=True)
            raw_b.attrs["flux_space"] = flux_space
        fits[b] = fit_continuum(raw_b, b, cfg, W, br)
        vals[b] = validate_fit(fits[b], cfg, br)
        sub = subtract_continuum(raw_b, fits[b])
        if br is not None and br.exclude_top_emission:
            em = sub.in_emission.to_numpy() & np.isfinite(sub.emis_mjy.to_numpy(float))
            idx = np.where(em)[0]
            drop = idx[np.argsort(sub.emis_mjy.to_numpy(float)[idx])[::-1][:int(br.exclude_top_emission)]]
            sub.loc[sub.index[drop], "role"] = "excluded"
            sub.loc[sub.index[drop], "in_emission"] = False
        subs[b] = sub
    rows = []
    for b in BAND_WINDOWS:
        row = aggregate_band(subs[b], fits[b], vals[b], target, r_ap_km, phase, epoch, flux_space, cfg)
        row["revision"] = revs[b].describe() if revs[b] is not None else ""
        rows.append(row)
    summary = pd.DataFrame(rows)
    chunks = []
    for b, d in subs.items():
        # the saved neighbourhood covers the continuum window *and* the emission window: a
        # revised continuum window need not bracket the band (2023 R1 phase 4)
        lo = min(W[b]["cont"][0], W[b]["em"][0])
        hi = max(W[b]["cont"][1], W[b]["em"][1])
        m = (d.wl >= lo - cfg.plot_margin_um) & (d.wl <= hi + cfg.plot_margin_um)
        if m.any():
            chunks.append(d[m])
    if chunks:
        pts = pd.concat(chunks, ignore_index=True)
    else:
        pts = next(iter(subs.values())).iloc[:0].copy()
    pts.insert(0, "target", slug(target))
    pts.insert(1, "r_ap_km", float(r_ap_km))
    pts["phase"] = int(phase)
    pts["epoch"] = int(epoch)
    pts["flux_space"] = flux_space
    cols = ["target", "r_ap_km", "phase", "epoch", "band", "role", "flux_space", "wl", "wlwidth",
            "flux", "err", "flux_raw", "err_raw", "distcorr_factor",
            "cont_mjy", "cont_err_mjy", "emis_mjy", "emis_err_mjy", "chi",
            "cont_raw_mjy", "emis_raw_mjy", "emis_raw_err_mjy",
            "r_hel", "r_obs", "v_hel_kms", "jd_utc", "detector", "sourceflag", "badphot",
            "frac_badpix_ap", "n_gaia", "gmag_eff", "filename"]
    pts = pts[[c for c in cols if c in pts.columns]]
    return dict(fits=fits, valid=vals, subs=subs, summary=summary, points=pts, raw=raw)
