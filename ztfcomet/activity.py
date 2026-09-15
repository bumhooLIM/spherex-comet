"""Heliocentric dependence of Af-rho: power-law trends per orbital leg, and the
time of peak activity relative to perihelion.

The trend quantity is the phase-corrected A(0 deg)f-rho, because the phase
angle changes along a leg and would otherwise masquerade as a heliocentric
trend.  Fits are done in log space, ``log Afrho = a - x log r_h``, weighted by
the propagated per-point errors, and every slope carries three uncertainties
because each answers a different question:

* *formal* -- from the weighted normal equations; right only if the errors
  are right and there is no intrinsic scatter;
* *scaled* -- formal times ``sqrt(chi2_red)`` when the fit is worse than the
  errors allow; an active comet varies, and the errors carry no jitter term;
* *bootstrap* -- the 16-84% interval from resampling the points, which
  assumes neither Gaussian residuals nor correct errors.

A power-law slope needs a lever arm: ``sigma_x ~ sigma_y / (sqrt(N) *
std(log r_h))``.  Several survey comets span < 0.02 dex in r_h, where any
slope is arithmetic rather than a measurement, so the reliability grade leads
with the baseline, not the point count.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

BANDS = {"r": "ZTF_r", "g": "ZTF_g"}
LN10 = np.log(10.0)


# ------------------------------------------------------------- selection
def select_points(phot, band, rho_km, max_centroid_frac=0.5):
    """Clean rows for one band and aperture, with the quality bookkeeping.

    Rows whose photometry centre sits more than *max_centroid_frac* of the
    aperture radius from the ephemeris are dropped and counted.  The
    ``flag_centroid`` gate fires at 3 FWHM -- a star-capture rule -- which at
    the 10,000 km aperture leaves ~5% of clean rows with the nucleus more than
    half an aperture radius off centre and a 10-40% low bias on a 1/rho coma.

    Returns
    -------
    kept : pandas.DataFrame
        With ``log_afrho``, ``log_afrho_err``, ``log_rh`` and ``leg``
        (``inbound`` where ``r_rate < 0``, else ``outbound``).
    counts : dict
        ``n_all``, ``n_clean``, ``n_contaminated``, ``n_offcentre``.
    """
    d = phot[np.isclose(phot["rho_km"], rho_km) & (phot["filter"] == BANDS[band])]
    n_contam = int(d["flag_contaminated"].fillna(False).astype(bool).sum()) if "flag_contaminated" in d else 0
    ok = (d["quality_ok"].astype(bool) & np.isfinite(d["afrho0_cm"]) & (d["afrho0_cm"] > 0)
          & np.isfinite(d["afrho0_cm_err"]) & (d["afrho0_cm_err"] > 0) & np.isfinite(d["r"]))
    clean = d[ok]
    if "centroid_shift_pix" in clean and "rho_pix" in clean:
        off = (clean["centroid_shift_pix"] / clean["rho_pix"]) > max_centroid_frac
    else:
        off = pd.Series(False, index=clean.index)
    kept = clean[~off].copy()
    kept["log_afrho"] = np.log10(kept["afrho0_cm"])
    kept["log_afrho_err"] = kept["afrho0_cm_err"] / (kept["afrho0_cm"] * LN10)
    kept["log_rh"] = np.log10(kept["r"])
    kept["leg"] = np.where(kept["r_rate"] < 0, "inbound", "outbound")
    return kept, dict(n_all=int(len(d)), n_clean=int(len(clean)),
                      n_contaminated=n_contam, n_offcentre=int(off.sum()))


# ------------------------------------------------------------------ fit
def _wls(x, y, w):
    W = w.sum()
    xm, ym = (w * x).sum() / W, (w * y).sum() / W
    sxx = (w * (x - xm) ** 2).sum()
    if sxx <= 0:
        return np.nan, np.nan, np.nan, np.nan
    b = (w * (x - xm) * (y - ym)).sum() / sxx
    a = ym - b * xm
    return a, b, np.sqrt(1.0 / W + xm ** 2 / sxx), np.sqrt(1.0 / sxx)


def fit_powerlaw(log_rh, log_afrho, log_err, n_boot=1000, seed=0):
    """Weighted fit of ``log Afrho = a - x log r_h`` with three error estimates.

    Returns a dict with ``n``, ``x`` (the power-law index, positive for
    activity that falls with distance), ``x_err`` (formal),
    ``x_err_scaled``, ``x_boot_lo``/``x_boot_hi`` (16-84%), ``a``
    (log10 Afrho at 1 au) and ``a_err``, ``chi2_red``, ``rms_dex`` (unweighted
    residual scatter), ``dlog_rh`` (baseline), ``median_log_err``,
    ``log_rh_mean`` (the weighted-mean abscissa) and ``level_err`` (the scaled
    uncertainty of the law there).  All fit fields are NaN with fewer than
    three points or no baseline.
    """
    x = np.asarray(log_rh, float)
    y = np.asarray(log_afrho, float)
    e = np.asarray(log_err, float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    x, y, e = x[keep], y[keep], e[keep]
    n = int(keep.sum())
    out = dict(n=n, x=np.nan, x_err=np.nan, x_err_scaled=np.nan, x_boot_lo=np.nan,
               x_boot_hi=np.nan, a=np.nan, a_err=np.nan, chi2_red=np.nan, rms_dex=np.nan,
               dlog_rh=float(np.ptp(x)) if n else np.nan,
               median_log_err=float(np.median(e)) if n else np.nan,
               log_rh_mean=np.nan, level_err=np.nan)
    if n < 3 or out["dlog_rh"] <= 0:
        return out
    w = 1.0 / e ** 2
    a, b, a_err, b_err = _wls(x, y, w)
    if not np.isfinite(b):
        return out
    resid = y - (a + b * x)
    dof = n - 2
    chi2_red = float((w * resid ** 2).sum() / dof) if dof > 0 else np.nan
    scale = np.sqrt(max(chi2_red, 1.0)) if np.isfinite(chi2_red) else 1.0
    # the level and the slope are uncorrelated at the weighted-mean abscissa;
    # from there the law can be evaluated anywhere with a two-term variance
    # (see evaluate_trend), which a_err -- the level at 1 au -- cannot give.
    xm = float((w * x).sum() / w.sum())
    level_err = float(scale / np.sqrt(w.sum()))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 2:
            continue
        _, bb, _, _ = _wls(x[idx], y[idx], w[idx])
        if np.isfinite(bb):
            boots.append(bb)
    b_lo, b_hi = (np.percentile(boots, [16, 84]) if len(boots) >= 50 else (np.nan, np.nan))
    out.update(x=-b, x_err=b_err, x_err_scaled=b_err * scale, x_boot_lo=-b_hi, x_boot_hi=-b_lo,
               a=a, a_err=a_err * scale, chi2_red=chi2_red,
               rms_dex=float(np.sqrt(np.mean(resid ** 2))), log_rh_mean=xm, level_err=level_err)
    return out


# ----------------------------------------------------------------- peak
def find_peak(t_days, log_afrho, log_err, kernel_days=None, n_boot=500, seed=0,
              edge_frac=0.1, min_drop_sigma=2.0, n_grid=400):
    """Time of peak activity, as the maximum of a kernel-smoothed log Afrho.

    The peak is *bracketed* -- and only then worth reporting -- when it lies
    inside the central ``1 - 2*edge_frac`` of the sampled span and the
    smoothed curve sits at least ``min_drop_sigma`` times the residual RMS
    below it at both ends.  A maximum at the edge of the data is a censoring
    statement, not a peak.  The kernel is a Gaussian of width
    ``max(10 d, 2 * median sampling gap)``; weights are ``1/err^2`` capped at
    ten times their median so one precise point cannot pull the curve.

    Returns a dict with ``t_peak`` (days from the reference time; negative
    is earlier), ``t_lo``/``t_hi`` (bootstrap 16-84%), ``bracketed``,
    ``kernel_days``, ``n``, ``peak_log_afrho``, ``rms_dex``, ``drop_lo``,
    ``drop_hi`` (how far the curve falls to the first and last sample).
    """
    t = np.asarray(t_days, float)
    y = np.asarray(log_afrho, float)
    e = np.asarray(log_err, float)
    keep = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    t, y, e = t[keep], y[keep], e[keep]
    n = int(keep.sum())
    out = dict(t_peak=np.nan, t_lo=np.nan, t_hi=np.nan, bracketed=False, interior=False,
               kernel_days=np.nan, n=n, peak_log_afrho=np.nan, rms_dex=np.nan, drop_lo=np.nan,
               drop_hi=np.nan)
    if n < 6 or np.ptp(t) <= 0:
        return out
    order = np.argsort(t)
    gap = float(np.median(np.diff(t[order]))) if n > 1 else 0.0
    h = float(kernel_days) if kernel_days else max(10.0, 2.0 * gap)
    grid = np.linspace(t.min(), t.max(), n_grid)
    w = 1.0 / e ** 2
    w = np.minimum(w, 10.0 * np.median(w))

    def smooth(tt, yy, ww):
        k = np.exp(-0.5 * ((grid[:, None] - tt[None, :]) / h) ** 2) * ww[None, :]
        return (k * yy[None, :]).sum(1) / np.maximum(k.sum(1), 1e-300)

    s = smooth(t, y, w)
    i = int(np.argmax(s))
    rms = float(np.std(y - np.interp(t, grid, s)))
    span = float(np.ptp(t))
    interior = bool((t.min() + edge_frac * span) < grid[i] < (t.max() - edge_frac * span))
    drop_lo, drop_hi = float(s[i] - s[0]), float(s[i] - s[-1])
    # bracketed: the curve falls on both sides.  interior: the maximum is not
    # at the edge of the coverage.  A maximum on a plateau (240P: flat before,
    # fading after) is interior but not bracketed, and still the right place
    # to split the phases.
    bracketed = bool(interior and drop_lo > min_drop_sigma * rms and drop_hi > min_drop_sigma * rms)
    rng = np.random.default_rng(seed)
    peaks = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.ptp(t[idx]) <= 0:
            continue
        peaks.append(grid[int(np.argmax(smooth(t[idx], y[idx], w[idx])))])
    lo, hi = (np.percentile(peaks, [16, 84]) if len(peaks) >= 50 else (np.nan, np.nan))
    out.update(t_peak=float(grid[i]), t_lo=float(lo), t_hi=float(hi), bracketed=bracketed,
               interior=interior, kernel_days=h, peak_log_afrho=float(s[i]), rms_dex=rms,
               drop_lo=drop_lo, drop_hi=drop_hi)
    return out


# ---------------------------------------------------------------- grade
def grade_fit(fit, counts=None):
    """Reliability grade for one fitted leg, with the reasons.

    A: N >= 10, baseline >= 0.2 dex, scaled slope error <= 0.5, chi2_red <= 3.
    B: N >= 5, baseline >= 0.1 dex, scaled slope error <= 1.0.
    C: N >= 5, baseline >= 0.05 dex, scaled slope error <= 2.0.
    D: a fit exists but constrains nothing.
    -: no fit.
    A grade is lowered one step when more than half of the band's frames
    were rejected for contamination or off-centre apertures: the surviving
    points are then a selection, not a sample.
    """
    n, dl, xe, c2 = fit.get("n", 0), fit.get("dlog_rh", np.nan), fit.get("x_err_scaled", np.nan), fit.get("chi2_red", np.nan)
    reasons = []
    if not np.isfinite(fit.get("x", np.nan)):
        return "-", ["no fit"]
    if n >= 10 and dl >= 0.2 and xe <= 0.5 and (not np.isfinite(c2) or c2 <= 3):
        g = "A"
    elif n >= 5 and dl >= 0.1 and xe <= 1.0:
        g = "B"
    elif n >= 5 and dl >= 0.05 and xe <= 2.0:
        g = "C"
    else:
        g = "D"
    if n < 5:
        reasons.append(f"only {n} points")
    if dl < 0.1:
        reasons.append(f"r_h baseline {dl:.3f} dex")
    if xe > 1.0:
        reasons.append(f"slope error {xe:.1f}")
    if np.isfinite(c2) and c2 > 3:
        reasons.append(f"chi2_red {c2:.1f}: scatter beyond the errors")
    if counts:
        rej = counts.get("n_contaminated", 0) + counts.get("n_offcentre", 0)
        tot = max(counts.get("n_all", 0), 1)
        if rej / tot > 0.5:
            reasons.append(f"{100 * rej / tot:.0f}% of frames rejected (contamination/off-centre)")
            g = {"A": "B", "B": "C", "C": "D", "D": "D"}[g]
    return g, reasons


# ------------------------------------------------------- broken power law
def _wlstsq(A, y, w):
    """Weighted linear least squares; returns coefficients, chi-square, covariance."""
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    resid = y - A @ coef
    cov = np.linalg.pinv((A * w[:, None]).T @ A)
    return coef, float((w * resid ** 2).sum()), cov


def _broken_design(x, xb):
    """Continuous broken line: y = a + b1 (x - xb) + (b2 - b1) max(x - xb, 0)."""
    return np.column_stack([np.ones_like(x), x - xb, np.maximum(x - xb, 0.0)])


def fit_broken(x, y, e, min_side=4, fixed_xb=None, n_boot=300, seed=0):
    """Two-slope continuous fit in log space, against a single slope.

    With *fixed_xb* the break is imposed (the r_h = 3 au hypothesis); otherwise
    it is searched over every abscissa that leaves at least *min_side* points
    on each side, and the best is compared with the single slope by
    ``dbic = BIC_single - BIC_broken`` (parameters 2 against 4), with the
    errors first inflated so the single slope has ``chi2_red = 1`` -- on raw
    chi-square, intrinsic scatter makes any kink look decisive.  ``preferred``
    is ``dbic > 6`` -- the break has to earn its two extra parameters.  The
    break's 16-84% interval comes from bootstrapping the search.

    Returns slopes ``b1`` (below the break) and ``b2`` (above), their errors
    scaled by ``sqrt(chi2_red)`` when the fit is worse than the errors, the
    break in the fitted abscissa (``xb``) and in linear units (``r_break``),
    both chi-squares, ``dbic``, ``preferred`` and ``testable``.
    """
    x, y, e = (np.asarray(v, float) for v in (x, y, e))
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    x, y, e = x[keep], y[keep], e[keep]
    n = int(keep.sum())
    out = dict(n=n, testable=False, preferred=False, xb=np.nan, r_break=np.nan,
               r_break_lo=np.nan, r_break_hi=np.nan, b1=np.nan, b1_err=np.nan, b2=np.nan,
               b2_err=np.nan, chi2_broken=np.nan, chi2_single=np.nan, dbic=np.nan)
    if n < 2 * min_side + 1:
        return out
    w = 1.0 / e ** 2
    order = np.sort(x)
    if fixed_xb is not None:
        if (x < fixed_xb).sum() < min_side or (x > fixed_xb).sum() < min_side:
            return out
        grid = np.array([float(fixed_xb)])
    else:
        grid = np.unique(order[min_side - 1:n - min_side + 1])
        if grid.size < 2:
            return out
    # single slope
    _, chi2_s, _ = _wlstsq(np.column_stack([np.ones_like(x), x]), y, w)

    def best(xv, yv, wv, g):
        chi = np.array([_wlstsq(_broken_design(xv, xb), yv, wv)[1] for xb in g])
        return g[int(np.argmin(chi))]

    xb = best(x, y, w, grid)
    coef, chi2_b, cov = _wlstsq(_broken_design(x, xb), y, w)
    dof = max(n - (3 if fixed_xb is not None else 4), 1)
    scale = np.sqrt(max(chi2_b / dof, 1.0))
    b1, d = coef[1], coef[2]
    b1_err = np.sqrt(max(cov[1, 1], 0)) * scale
    b2_err = np.sqrt(max(cov[1, 1] + cov[2, 2] + 2 * cov[1, 2], 0)) * scale
    k_b = 3 if fixed_xb is not None else 4
    # BIC on raw chi-square is meaningless when the errors carry no intrinsic
    # scatter: with chi2_red ~ 10 every kink in the last four points is
    # "decisive".  Inflate the errors so the *simpler* model has chi2_red = 1,
    # then ask whether the break still earns its parameters.
    inflate = max(chi2_s / max(n - 2, 1), 1.0)
    dbic = (chi2_s / inflate + 2 * np.log(n)) - (chi2_b / inflate + k_b * np.log(n))
    lo = hi = np.nan
    if fixed_xb is None and n_boot:
        rng = np.random.default_rng(seed)
        bs = []
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            xs = np.sort(x[idx])
            g = np.unique(xs[min_side - 1:n - min_side + 1])
            if g.size >= 2:
                bs.append(best(x[idx], y[idx], w[idx], g))
        if len(bs) >= 30:
            lo, hi = np.percentile(bs, [16, 84])
    out.update(testable=True, preferred=bool(dbic > 6), xb=float(xb), r_break=float(10 ** xb),
               r_break_lo=float(10 ** lo) if np.isfinite(lo) else np.nan,
               r_break_hi=float(10 ** hi) if np.isfinite(hi) else np.nan,
               b1=float(b1), b1_err=float(b1_err), b2=float(b1 + d), b2_err=float(b2_err),
               chi2_broken=chi2_b, chi2_single=chi2_s, dbic=float(dbic))
    return out


def fit_broken_powerlaw(log_rh, log_afrho, log_err, **kw):
    """:func:`fit_broken` for Afrho: indices ``x_inner`` (r_h below the break)
    and ``x_outer`` are the negated slopes."""
    f = fit_broken(log_rh, log_afrho, log_err, **kw)
    return dict(n=f["n"], testable=f["testable"], preferred=f["preferred"],
                r_break=f["r_break"], r_break_lo=f["r_break_lo"], r_break_hi=f["r_break_hi"],
                x_inner=-f["b1"], x_inner_err=f["b1_err"], x_outer=-f["b2"], x_outer_err=f["b2_err"],
                chi2_broken=f["chi2_broken"], chi2_single=f["chi2_single"], dbic=f["dbic"])


# ---------------------------------------------------------------- colour
def colour_pairs(phot, rho_km, max_dt_days=1.0, max_centroid_frac=0.5):
    """g and r measurements of the same aperture within *max_dt_days*.

    The colour is reported two ways.  ``excess`` is
    ``-2.5 log10(Afrho_g / Afrho_r)``: the dust colour relative to the Sun,
    which needs no solar magnitudes and in which the phase correction
    cancels.  ``gr`` is the observed g - r from the calibrated aperture
    magnitudes, where those exist.  Each g frame is paired with the nearest r
    frame; an r frame is used at most once.
    """
    g, _ = select_points(phot, "g", rho_km, max_centroid_frac)
    r, _ = select_points(phot, "r", rho_km, max_centroid_frac)
    if g.empty or r.empty:
        return pd.DataFrame()
    r = r.sort_values("obsjd").reset_index(drop=True)
    used, rows = set(), []
    for _, gg in g.sort_values("obsjd").iterrows():
        d = np.abs(r["obsjd"].to_numpy() - gg["obsjd"])
        for j in np.argsort(d):
            if d[j] > max_dt_days:
                break
            if j in used:
                continue
            used.add(j)
            rr = r.iloc[j]
            ex = -2.5 * np.log10(gg["afrho0_cm"] / rr["afrho0_cm"])
            ex_err = 2.5 / LN10 * np.hypot(gg["afrho0_cm_err"] / gg["afrho0_cm"],
                                            rr["afrho0_cm_err"] / rr["afrho0_cm"])
            gr = gr_err = np.nan
            if "filter_mag" in gg and "filter_mag" in rr and np.isfinite(gg.get("filter_mag", np.nan)) \
                    and np.isfinite(rr.get("filter_mag", np.nan)):
                gr = gg["filter_mag"] - rr["filter_mag"]
                gr_err = float(np.hypot(gg.get("filter_mag_err", np.nan), rr.get("filter_mag_err", np.nan)))
            rows.append(dict(obsjd=0.5 * (gg["obsjd"] + rr["obsjd"]), dt_days=float(d[j]),
                             r=0.5 * (gg["r"] + rr["r"]), log_rh=0.5 * (gg["log_rh"] + rr["log_rh"]),
                             r_rate=gg["r_rate"], excess=float(ex), excess_err=float(ex_err),
                             gr=gr, gr_err=gr_err, file_g=gg["file"], file_r=rr["file"]))
            break
    return pd.DataFrame(rows)


def fit_step(x, y, e, min_side=4, n_boot=300, seed=0):
    """Two-level step in *y* at a change point in *x*, against one level.

    The right shape for a colour: a dust colour that changes with distance
    does so between two plateaus, and a continuous kink cannot hold two
    plateaus at once.  The split is searched over every abscissa leaving
    *min_side* points on each side; weighted means on each side; compared
    with a single weighted mean by ``dbic = BIC_one - BIC_two`` (parameters 1
    against 3) on errors inflated so the one-level model has ``chi2_red = 1``,
    ``preferred`` when ``dbic > 6``.  Errors on the two levels are
    scaled by ``sqrt(chi2_red)``; the change point's 16-84% interval comes
    from bootstrapping the search.
    """
    x, y, e = (np.asarray(v, float) for v in (x, y, e))
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    x, y, e = x[keep], y[keep], e[keep]
    n = int(keep.sum())
    out = dict(n=n, testable=False, preferred=False, x_change=np.nan, x_change_lo=np.nan,
               x_change_hi=np.nan, level_below=np.nan, level_below_err=np.nan, level_above=np.nan,
               level_above_err=np.nan, delta=np.nan, delta_err=np.nan, chi2_two=np.nan, chi2_one=np.nan,
               dbic=np.nan)
    if n < 2 * min_side + 1:
        return out
    w = 1.0 / e ** 2
    mean_one = (w * y).sum() / w.sum()
    chi2_one = float((w * (y - mean_one) ** 2).sum())

    def search(xv, yv, wv):
        xs = np.sort(xv)
        cands = np.unique(xs[min_side - 1:len(xs) - min_side])   # split between cand and the next
        best, best_chi = np.nan, np.inf
        for c in cands:
            lo, hi = xv <= c, xv > c
            if hi.sum() < min_side:
                continue
            ml = (wv[lo] * yv[lo]).sum() / wv[lo].sum(); mh = (wv[hi] * yv[hi]).sum() / wv[hi].sum()
            chi = (wv[lo] * (yv[lo] - ml) ** 2).sum() + (wv[hi] * (yv[hi] - mh) ** 2).sum()
            if chi < best_chi:
                best, best_chi = c, chi
        return best, best_chi

    c, chi2_two = search(x, y, w)
    if not np.isfinite(c):
        return out
    lo, hi = x <= c, x > c
    ml, mh = (w[lo] * y[lo]).sum() / w[lo].sum(), (w[hi] * y[hi]).sum() / w[hi].sum()
    scale = np.sqrt(max(chi2_two / max(n - 3, 1), 1.0))
    el, eh = scale / np.sqrt(w[lo].sum()), scale / np.sqrt(w[hi].sum())
    # the change point sits between the last point below and the first above
    xc = 0.5 * (c + x[hi].min())
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        cb, _ = search(x[idx], y[idx], w[idx])
        if np.isfinite(cb):
            bs.append(cb)
    lo_q, hi_q = (np.percentile(bs, [16, 84]) if len(bs) >= 30 else (np.nan, np.nan))
    inflate = max(chi2_one / max(n - 1, 1), 1.0)      # see fit_broken: BIC on scaled errors
    dbic = (chi2_one / inflate + np.log(n)) - (chi2_two / inflate + 3 * np.log(n))
    out.update(testable=True, preferred=bool(dbic > 6),
               x_change=float(xc), x_change_lo=float(lo_q), x_change_hi=float(hi_q),
               level_below=float(ml), level_below_err=float(el), level_above=float(mh),
               level_above_err=float(eh), delta=float(mh - ml), delta_err=float(np.hypot(el, eh)),
               chi2_two=float(chi2_two), chi2_one=chi2_one, dbic=float(dbic))
    return out


def fit_colour_trend(log_rh, colour, err, min_side=4, n_boot=300):
    """Linear colour trend with log r_h, and a step change-point test.

    Returns the linear slope (mag per dex, scaled error) and RMS, and the
    :func:`fit_step` result in r_h: ``r_change`` with its interval, the colour
    level on either side, ``delta_colour`` (above minus below, with error),
    and whether the change is ``change_preferred`` by BIC.
    """
    lin = fit_powerlaw(log_rh, colour, err, n_boot=n_boot)
    st = fit_step(log_rh, colour, err, min_side=min_side, n_boot=n_boot)
    tolin = lambda v: float(10 ** v) if np.isfinite(v) else np.nan
    return dict(n=lin["n"], slope=-lin["x"], slope_err=lin["x_err_scaled"], rms=lin["rms_dex"],
                dlog_rh=lin["dlog_rh"], r_change=tolin(st["x_change"]), r_change_lo=tolin(st["x_change_lo"]),
                r_change_hi=tolin(st["x_change_hi"]), colour_below=st["level_below"],
                colour_below_err=st["level_below_err"], colour_above=st["level_above"],
                colour_above_err=st["level_above_err"], delta_colour=st["delta"],
                delta_colour_err=st["delta_err"], dbic=st["dbic"], change_preferred=st["preferred"])


# ------------------------------------------------------- tails and outbursts
def sparse_tail(log_rh, max_n=4, max_frac=0.10, min_gap=0.10, min_core=5):
    """Mark a small group of points isolated at either end of the r_h range.

    A power-law fit is a lever: a handful of points far from the rest carry
    most of the leverage and pull the slope away from what the bulk of the
    data shows -- 10P has four frames at 3.6 au beyond a 0.135 dex gap from
    the other 56, and they turn a steep onset into a moderate single law.
    A group of at most ``max(max_n, max_frac * n)`` points separated from
    the rest by a gap of at least *min_gap* dex is a tail; it is reported
    and fitted separately, not silently dropped.

    Returns a boolean mask, True on tail points.
    """
    x = np.asarray(log_rh, float)
    n = len(x)
    tail = np.zeros(n, bool)
    if n < min_core + 1:
        return tail
    order = np.argsort(x)
    xs = x[order]
    kmax = int(max(max_n, np.floor(max_frac * n)))
    for end in ("high", "low"):
        for k in range(1, kmax + 1):
            if n - k < min_core:
                break
            gap = (xs[n - k] - xs[n - k - 1]) if end == "high" else (xs[k] - xs[k - 1])
            if gap >= min_gap:
                idx = order[n - k:] if end == "high" else order[:k]
                tail[idx] = True
                break
    return tail


def _extrapolate(t_prev, y_prev, t_at, max_gap):
    """Recent level extended to *t_at*: a weighted line through the previous
    points when there are enough of them and the gap is short, otherwise
    their median.  A steep but smooth trend must extrapolate to itself."""
    t_prev, y_prev = np.asarray(t_prev, float), np.asarray(y_prev, float)
    if len(t_prev) >= 3 and np.ptp(t_prev) >= 2.0 and (t_at - t_prev.max()) <= max_gap:
        A = np.column_stack([np.ones_like(t_prev), t_prev - t_prev.max()])
        c, *_ = np.linalg.lstsq(A, y_prev, rcond=None)
        return float(c[0] + c[1] * (t_at - t_prev.max()))
    return float(np.median(y_prev))


def detect_outbursts(t_days, log_afrho, min_rise_dex=0.30, baseline_n=5, baseline_days=60.0,
                     end_tol_dex=0.10, confirm=True, min_points=3):
    """Find brightenings that are too sudden and too sustained to be the trend.

    An outburst is a *jump* above the recent trend that the next frame
    shares and that then decays: 217P at 3.10 au went from ~29 to 162 cm in
    one step and took forty days to come most of the way back.  The
    single-frame anomaly flag (:func:`ztfcomet.phot.flag_anomalies`) has
    already removed lone spikes, so what is left that jumps *and stays up*
    is real.

    The comparison is with the recent trend *extrapolated*, not the recent
    median: a comet brightening steeply but smoothly toward its peak sits
    above the median of its previous frames at every step, and a median
    test opened a window at the start of such a series that never closed.
    Walking forward in time, the previous *baseline_n* frames within
    *baseline_days* (at least three) give a line; a frame more than
    *min_rise_dex* above that line at its own time, confirmed by the next
    frame keeping at least half the excess above the same line, opens a
    window.  The window closes at the first frame whose three-frame running
    median is within *end_tol_dex* of the pre-outburst line extrapolated to
    it, or at the end of the data (``ended`` False).  An outburst *decays*:
    a window that has not closed and whose second half is not fainter than
    its first is a rapid onset of activity, not an outburst -- 261P jumps
    0.4 dex and keeps rising to perihelion -- and is left to the trend.
    Windows shorter than *min_points* frames are ignored: the single-frame
    anomaly flag already handles spikes, and two frames cannot show a decay.

    Returns a list of dicts (``i_start``, ``i_end`` -- inclusive indices in
    time order -- ``t_start``, ``t_end``, ``rise_dex``, ``ended``,
    ``n_points``) and the time-sorted index array, so callers can map back.
    """
    t = np.asarray(t_days, float)
    y = np.asarray(log_afrho, float)
    order = np.argsort(t)
    t, y = t[order], y[order]
    n = len(t)
    windows = []
    i = 1
    while i < n:
        prev = [j for j in range(i - 1, -1, -1) if t[i] - t[j] <= baseline_days][:baseline_n]
        if len(prev) < 3:
            i += 1
            continue
        tp, yp = t[prev], y[prev]
        level = lambda tt: _extrapolate(tp, yp, tt, baseline_days / 2)   # noqa: E731
        rise = y[i] - level(t[i])
        if rise < min_rise_dex or (confirm and (i + 1 >= n or y[i + 1] - level(t[i + 1]) < 0.5 * rise)):
            i += 1
            continue
        j = i + 1
        ended = False
        while j < n:
            run = float(np.median(y[max(j - 1, i):min(j + 2, n)]))
            if run - level(t[j]) <= end_tol_dex:
                ended = True
                break
            j += 1
        i_end = j - 1 if ended else n - 1
        seg = y[i:i_end + 1]
        decays = ended or (len(seg) >= 4 and np.median(seg[len(seg) // 2:]) < np.median(seg[:len(seg) // 2]) - 0.05)
        if len(seg) >= min_points and decays:
            windows.append(dict(i_start=int(i), i_end=int(i_end), t_start=float(t[i]), t_end=float(t[i_end]),
                                rise_dex=float(rise), baseline_log=float(level(t[i])), ended=ended,
                                n_points=int(i_end - i + 1)))
            i = i_end + 1
        else:
            i += 1
    return windows, order


def outburst_mask(t_days, log_afrho, **kw):
    """Boolean mask of points inside outburst windows, in the input order."""
    windows, order = detect_outbursts(t_days, log_afrho, **kw)
    mask = np.zeros(len(order), bool)
    for w in windows:
        mask[order[w["i_start"]:w["i_end"] + 1]] = True
    return mask, windows


def smooth_trend(log_rh, log_afrho, log_err, n_grid=120, bandwidth=None, min_weight=0.5):
    """Kernel-smoothed log Afrho against log r_h over every point, for the eye.

    A Gaussian kernel in log r_h, weights ``1/err^2`` capped at ten times
    their median, bandwidth ``max(0.03 dex, 2 x median spacing)``.  Points
    from both phases at the same r_h are averaged -- this is a guide to the
    overall shape, not a model.  Returned only where the kernel-weighted
    count exceeds *min_weight* points, so the curve stops with the data.
    """
    x = np.asarray(log_rh, float); y = np.asarray(log_afrho, float); e = np.asarray(log_err, float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    x, y, e = x[keep], y[keep], e[keep]
    if len(x) < 4 or np.ptp(x) <= 0:
        return np.array([]), np.array([])
    xs = np.unique(x)
    h = float(bandwidth) if bandwidth else max(0.03, 2.0 * float(np.median(np.diff(xs)))) if len(xs) > 1 else 0.03
    w = 1.0 / e ** 2
    w = np.minimum(w, 10.0 * np.median(w)) / np.median(w)
    grid = np.linspace(x.min(), x.max(), n_grid)
    k = np.exp(-0.5 * ((grid[:, None] - x[None, :]) / h) ** 2)
    kw = k * w[None, :]
    count = k.sum(1)
    sm = (kw * y[None, :]).sum(1) / np.maximum(kw.sum(1), 1e-300)
    ok = count >= min_weight
    return grid[ok], sm[ok]


# --------------------------------------------------------------- per target
def _fit_row(base, split, leg, sub, counts, n_boot, extra=None):
    fit = fit_powerlaw(sub["log_rh"], sub["log_afrho"], sub["log_afrho_err"], n_boot=n_boot)
    g, why = grade_fit(fit, counts)
    row = dict(**base, split=split, leg=leg, rh_min=float(sub["r"].min()), rh_max=float(sub["r"].max()),
               grade=g, reasons="; ".join(why), **fit)
    if extra:
        row.update(extra)
    return row, fit


def analyse(phot, target, tp_jd=None, bands=("r", "g"), rhos=(10000, 20000),
            min_side=3, n_boot=1000, break_min_points=10, break_min_dlog=0.10):
    """Trends, peak, breaks, colour and outbursts for one comet.

    Returns ``(trends, peaks, breaks, colours, outbursts)``.

    Phases are defined by the activity peak: a comet sampled on both sides
    of T_p gets a smoothed peak, and when that maximum is *interior* to the
    coverage the series splits into ``rising`` and ``fading`` there (a
    plateau before the peak still splits -- ``peak_bracketed`` records
    whether the curve falls on both sides).  When the maximum sits at an
    edge the whole series is one phase.  A one-sided comet is one leg named
    by its orbital direction.  Perihelion-split rows are kept under
    ``split == "perihelion"`` for reference.

    Before any primary fit, outburst windows (:func:`detect_outbursts`) and
    sparse tails (:func:`sparse_tail`) are removed from the leg and counted
    (``n_outburst``, ``n_tail``); the tail and the pre-/post-outburst
    stretches are fitted separately under ``split == "tail"`` and
    ``split == "outburst"``.  Where a free-break law is preferred on a
    primary leg, the two sides are also fitted as plain power laws under
    ``split == "segment"`` -- the "divided" trend.
    """
    trends, peaks, breaks, colours, outbursts = [], [], [], [], []
    for band in bands:
        for rho in rhos:
            pts, counts = select_points(phot, band, rho)
            base = dict(target=target, band=band, rho_km=rho, **counts)
            if len(pts) == 0:
                continue
            pts = pts.sort_values("obsjd").reset_index(drop=True)
            # outbursts on the whole clean series, in time
            omask, windows = outburst_mask(pts["obsjd"].to_numpy(), pts["log_afrho"].to_numpy())
            for w in windows:
                sel = pts[(pts["obsjd"] >= w["t_start"]) & (pts["obsjd"] <= w["t_end"])]
                if sel.empty:
                    continue
                outbursts.append(dict(target=target, band=band, rho_km=rho, jd_start=w["t_start"],
                                      jd_end=w["t_end"], t_start=(w["t_start"] - tp_jd) if tp_jd else np.nan,
                                      t_end=(w["t_end"] - tp_jd) if tp_jd else np.nan,
                                      rh_start=float(sel["r"].iloc[0]), rh_end=float(sel["r"].iloc[-1]),
                                      rise_dex=w["rise_dex"], n_points=w["n_points"], ended=w["ended"]))
            pts["outburst"] = omask
            quiet = pts[~pts["outburst"]]
            inb, outb = quiet[quiet["leg"] == "inbound"], quiet[quiet["leg"] == "outbound"]
            two_sided = len(inb) >= min_side and len(outb) >= min_side
            legs = [("all", "all", quiet)]
            for name, sub in (("inbound", inb), ("outbound", outb)):
                if len(sub) >= 3:
                    legs.append(("perihelion", name, sub))
            pk, t = None, None
            if tp_jd is not None and np.isfinite(tp_jd) and len(quiet) >= 6:
                t = quiet["obsjd"].to_numpy() - tp_jd
                pk = find_peak(t, quiet["log_afrho"], quiet["log_afrho_err"], n_boot=max(n_boot // 2, 200))
            if two_sided and pk is not None:
                peaks.append(dict(target=target, band=band, rho_km=rho, n_inbound=len(inb),
                                  n_outbound=len(outb), t_min=float(t.min()), t_max=float(t.max()), **pk))
                before, after = quiet[t < pk["t_peak"]], quiet[t >= pk["t_peak"]]
                if pk["interior"] and len(before) >= min_side and len(after) >= min_side:
                    legs += [("peak", "rising", before), ("peak", "fading", after)]
                else:
                    legs.append(("peak", "rising" if pk["t_peak"] > np.median(t) else "fading", quiet))
            for split, leg, sub in legs:
                primary = (split == "peak") if two_sided else (split == "perihelion")
                tail = sparse_tail(sub["log_rh"].to_numpy()) if primary else np.zeros(len(sub), bool)
                core = sub[~tail]
                n_ob = int(pts["outburst"].sum()) if primary else 0
                row, fit = _fit_row(base, split, leg, core, counts, n_boot, dict(
                    primary=primary, two_sided=two_sided, t_peak=(pk["t_peak"] if pk else np.nan),
                    peak_bracketed=(bool(pk["bracketed"]) if pk else False),
                    peak_interior=(bool(pk["interior"]) if pk else False),
                    n_tail=int(tail.sum()), n_outburst=n_ob, segment_of="", has_segments=False))
                trends.append(row)
                if tail.sum() >= 3:
                    trow, _ = _fit_row(base, "tail", f"{leg} tail", sub[tail], counts, n_boot,
                                       dict(primary=False, two_sided=two_sided, segment_of=leg))
                    trends.append(trow)
                if primary and n_ob:
                    lo = pts[pts["obsjd"] < min(w["t_start"] for w in windows)]
                    hi = pts[pts["obsjd"] > max(w["t_end"] for w in windows)]
                    for name, seg in (("pre-outburst", lo), ("post-outburst", hi)):
                        seg = seg[seg["leg"].isin(sub["leg"].unique())] if split == "perihelion" else seg
                        if len(seg) >= 5:
                            orow, _ = _fit_row(base, "outburst", name, seg, counts, n_boot,
                                               dict(primary=False, two_sided=two_sided, segment_of=leg))
                            trends.append(orow)
                if fit["n"] >= break_min_points and fit["dlog_rh"] >= break_min_dlog:
                    free = fit_broken_powerlaw(core["log_rh"], core["log_afrho"], core["log_afrho_err"],
                                               n_boot=max(n_boot // 3, 100))
                    at3 = fit_broken_powerlaw(core["log_rh"], core["log_afrho"], core["log_afrho_err"],
                                              fixed_xb=np.log10(3.0), n_boot=0)
                    breaks.append(dict(target=target, band=band, rho_km=rho, split=split, leg=leg,
                                       primary=primary, grade=row["grade"], rh_min=float(core["r"].min()),
                                       rh_max=float(core["r"].max()), x_single=fit["x"],
                                       x_single_err=fit["x_err_scaled"], **free,
                                       at3_testable=at3["testable"], x_lt3=at3["x_inner"],
                                       x_lt3_err=at3["x_inner_err"], x_gt3=at3["x_outer"],
                                       x_gt3_err=at3["x_outer_err"], at3_dbic=at3["dbic"]))
                    if primary and free["preferred"]:
                        row["has_segments"] = True
                        rb = free["r_break"]
                        for name, seg in ((f"{leg} r_h<{rb:.2f}", core[core["r"] < rb]),
                                          (f"{leg} r_h>{rb:.2f}", core[core["r"] >= rb])):
                            if len(seg) >= 3:
                                srow, _ = _fit_row(base, "segment", name, seg, counts, n_boot,
                                                   dict(primary=True, two_sided=two_sided, segment_of=leg,
                                                        r_break=rb))
                                trends.append(srow)
        if band != "r":
            continue
        for rho in rhos:
            pairs = colour_pairs(phot, rho)
            if len(pairs) < 8:
                continue
            series = [("all", pairs)]
            row = [r for r in trends if r["band"] == "r" and r["rho_km"] == rho and r["split"] == "peak"]
            tpk = row[0]["t_peak"] if row and row[0]["peak_interior"] and tp_jd is not None else None
            if tpk is not None:
                tt = pairs["obsjd"].to_numpy() - tp_jd
                series += [("rising", pairs[tt < tpk]), ("fading", pairs[tt >= tpk])]
            for leg, sub in series:
                if len(sub) < 8 or np.ptp(sub["log_rh"]) < 0.15:
                    continue
                ct = fit_colour_trend(sub["log_rh"], sub["excess"], sub["excess_err"], n_boot=max(n_boot // 3, 100))
                colours.append(dict(target=target, rho_km=rho, leg=leg, n_pairs=len(sub),
                                    rh_min=float(sub["r"].min()), rh_max=float(sub["r"].max()),
                                    excess_median=float(sub["excess"].median()),
                                    gr_median=float(sub["gr"].median()) if sub["gr"].notna().any() else np.nan, **ct))
    return (pd.DataFrame(trends), pd.DataFrame(peaks), pd.DataFrame(breaks),
            pd.DataFrame(colours), pd.DataFrame(outbursts))


# ------------------------------------------------------- value at an epoch
EPOCH_METHODS = ("direct", "trend", "trend_extrap", "none")
_PRIMARY_SPLITS = ("peak", "perihelion")


def evaluate_trend(row, log_rh0, include_scatter=True):
    """log Afrho of one fitted power law at *log_rh0*, with its uncertainty.

    ``log Afrho = a - x log r_h`` has its level and slope uncorrelated at the
    weighted-mean abscissa of the fit, so the variance anywhere else is
    ``level_err**2 + (log_rh0 - log_rh_mean)**2 * x_err_scaled**2`` -- the
    scaled errors, since the fits are usually worse than their errors.  With
    *include_scatter* the RMS of the points about the law is added in
    quadrature: a comet at one epoch sits that far from its own trend, and
    the question asked is its brightness then, not the mean law.  Rows fitted
    before ``log_rh_mean`` existed fall back to ``a_err`` (the level at 1 au)
    and the slope error with no correlation term, an overestimate.

    Returns ``(log_afrho, log_afrho_err)``, NaN without a fit.
    """
    a, x = float(row["a"]), float(row["x"])
    if not (np.isfinite(a) and np.isfinite(x)):
        return np.nan, np.nan
    xm, le = float(row.get("log_rh_mean", np.nan)), float(row.get("level_err", np.nan))
    xe = float(row.get("x_err_scaled", np.nan))
    xe = xe if np.isfinite(xe) else 0.0
    if np.isfinite(xm) and np.isfinite(le):
        var = le ** 2 + (log_rh0 - xm) ** 2 * xe ** 2
    else:
        var = float(row.get("a_err", 0.0)) ** 2 + log_rh0 ** 2 * xe ** 2
    rms = float(row.get("rms_dex", np.nan))
    if include_scatter and np.isfinite(rms):
        var += rms ** 2
    return float(a - x * log_rh0), float(np.sqrt(var))


def _leg_at_epoch(trends, jd0, tp_jd=None, arc=None):
    """Which fitted phase of the comet the epoch *jd0* belongs to.

    Peak-split comets divide at ``t_peak``; a single-phase comet (one leg
    fitted, or a maximum at the edge of its coverage) is one phase; a
    perihelion-split comet divides at T_p, or by the caller's *arc*
    (``"in"`` / ``"out"``) when T_p is unknown.  Returns the leg name and the
    primary rows.
    """
    prim = trends[trends["primary"].astype(bool) & trends["split"].isin(_PRIMARY_SPLITS)]
    if prim.empty:
        return "", prim
    has_tp = tp_jd is not None and np.isfinite(tp_jd)
    if (prim["split"] == "peak").any():
        rows = prim[prim["split"] == "peak"]
        if set(rows["leg"]) >= {"rising", "fading"} and has_tp:
            leg = "rising" if (jd0 - tp_jd) < float(rows["t_peak"].iloc[0]) else "fading"
        else:
            leg = str(rows["leg"].iloc[0])
        return leg, rows
    if has_tp:
        leg = "inbound" if jd0 < tp_jd else "outbound"
    elif arc:
        leg = "inbound" if str(arc).lower().startswith("in") else "outbound"
    else:
        leg = str(prim["leg"].iloc[0])
    return leg, prim


def _epoch_rows(trends, leg):
    """The fitted laws of one leg: its primary row, its segments and its tail."""
    if not leg or trends is None or len(trends) == 0:
        return pd.DataFrame()
    seg_of = trends["segment_of"].fillna("").astype(str) if "segment_of" in trends else pd.Series("", index=trends.index)
    m = ((trends["primary"].astype(bool) & (trends["leg"] == leg) & trends["split"].isin(_PRIMARY_SPLITS))
         | (trends["split"].isin(["segment", "tail"]) & (seg_of == leg)))
    rows = trends[m]
    return rows[np.isfinite(rows["a"]) & np.isfinite(rows["x"])]


def _with_linear(out):
    y, e = out["log_afrho"], out["log_afrho_err"]
    if np.isfinite(y):
        out.update(afrho_cm=float(10 ** y), afrho_err_cm=float(10 ** y * LN10 * e) if np.isfinite(e) else np.nan,
                   afrho_lo_cm=float(10 ** (y - e)) if np.isfinite(e) else np.nan,
                   afrho_hi_cm=float(10 ** (y + e)) if np.isfinite(e) else np.nan)
    else:
        out.update(afrho_cm=np.nan, afrho_err_cm=np.nan, afrho_lo_cm=np.nan, afrho_hi_cm=np.nan)
    return out


def afrho_at_epoch(pts, trends, rh0, jd0, jd_lo=None, jd_hi=None, tp_jd=None, arc=None,
                   outbursts=None, pad_days=5.0, min_direct=2, max_extrap_dex=0.10,
                   extrap_grades=("A", "B", "C"), tol=0.005):
    """Af-rho of one comet at one observing epoch, from its ZTF series.

    The epoch is a window ``[jd_lo, jd_hi]`` around *jd0* (a SPHEREx phase
    group, typically days to a month) at mean heliocentric distance *rh0*.
    Three estimates are tried in order of how much they assume:

    1. ``direct`` -- at least *min_direct* clean frames fall within
       *pad_days* of the window.  Each frame is moved to *rh0* along the
       leg's fitted slope (the most local law of grade A-C; otherwise left
       where it is), and the weighted mean is the value: this is the
       measurement itself, and a trend is only needed to correct for the
       r_h the comet covered inside the window.  The error is the scaled
       error of the mean plus the slope's contribution; ``rms_dex`` is the
       scatter of the frames.
    2. ``trend`` -- no frames then, but *rh0* lies inside a fitted law of the
       phase the epoch falls in (:func:`_leg_at_epoch`): the most local one
       among the primary law, its segments and its tail, evaluated with
       :func:`evaluate_trend`, scatter included.  Inside the data any grade
       serves -- the level of a flat, well-sampled series is measured even
       where its slope is not.
    3. ``trend_extrap`` -- *rh0* is beyond the primary law's r_h range by at
       most *max_extrap_dex* and the law is at least grade C (an
       unconstrained slope cannot be extended).  ``dlog_extrap`` records
       how far.

    Otherwise ``none``, and ``note`` says why: no ZTF frames at this
    aperture, ZTF fitted only the other phase of the orbit, no trend at
    all, or *rh0* too far outside the fitted range.  A single frame inside
    the window is used only when no trend can be evaluated, with an error
    floor of 0.1 dex.

    Parameters
    ----------
    pts : pandas.DataFrame
        Clean frames of one band and aperture from :func:`select_points`.
    trends : pandas.DataFrame
        The trend rows of the same comet, band and aperture.
    outbursts : pandas.DataFrame, optional
        Outburst windows (``jd_start``, ``jd_end``) of the same series, to
        count the frames of the window that lie in one.

    Returns
    -------
    dict
        ``method``, ``leg``, ``log_afrho``, ``log_afrho_err``, ``afrho_cm``,
        ``afrho_err_cm``, ``afrho_lo_cm``, ``afrho_hi_cm`` (the 1-sigma range
        in linear units), ``n`` (frames or fit points used), ``n_window``,
        ``n_outburst_window``, ``grade`` (of the law used), ``rms_dex``,
        ``dlog_extrap``, ``note``.
    """
    log_rh0 = float(np.log10(rh0))
    out = dict(method="none", leg="", log_afrho=np.nan, log_afrho_err=np.nan, n=0, n_window=0,
               n_outburst_window=0, grade="", rms_dex=np.nan, dlog_extrap=0.0, note="")
    if pts is None or len(pts) == 0:
        out["note"] = "no clean ZTF r-band frames at this aperture"
        return _with_linear(out)
    jd_lo = jd0 if jd_lo is None else jd_lo
    jd_hi = jd0 if jd_hi is None else jd_hi
    win = pts[(pts["obsjd"] >= jd_lo - pad_days) & (pts["obsjd"] <= jd_hi + pad_days)]
    out["n_window"] = int(len(win))
    if outbursts is not None and len(outbursts) and len(win):
        m = np.zeros(len(win), bool)
        for _, w in outbursts.iterrows():
            m |= ((win["obsjd"] >= w["jd_start"]) & (win["obsjd"] <= w["jd_end"])).to_numpy()
        out["n_outburst_window"] = int(m.sum())
    have_trends = trends is not None and len(trends) > 0
    leg, prim = _leg_at_epoch(trends, jd0, tp_jd, arc) if have_trends else ("", pd.DataFrame())
    out["leg"] = leg
    rows = _epoch_rows(trends, leg) if have_trends else pd.DataFrame()
    primary = rows[rows["split"].isin(_PRIMARY_SPLITS)] if len(rows) else rows
    primary = primary.iloc[0] if len(primary) else None
    inside = rows[(rows["rh_min"] * (1 - tol) <= rh0) & (rh0 <= rows["rh_max"] * (1 + tol))] if len(rows) else rows
    local = inside.iloc[int(np.argmin(np.log10(inside["rh_max"] / inside["rh_min"])))] if len(inside) else None

    # 1. the frames inside the window
    if len(win) >= min_direct:
        src = local if local is not None and local["grade"] in extrap_grades else \
            (primary if primary is not None and primary["grade"] in extrap_grades else None)
        slope = float(src["x"]) if src is not None else 0.0
        slope_err = float(src["x_err_scaled"]) if src is not None else 0.0
        d = win["log_rh"].to_numpy(float) - log_rh0
        y = win["log_afrho"].to_numpy(float) + slope * d
        e = win["log_afrho_err"].to_numpy(float)
        w = 1.0 / e ** 2
        W = w.sum()
        ym = float((w * y).sum() / W)
        chi2_red = float((w * (y - ym) ** 2).sum() / (len(y) - 1)) if len(y) > 1 else np.nan
        scale = float(np.sqrt(max(chi2_red, 1.0))) if np.isfinite(chi2_red) else 1.0
        dm = float((w * d).sum() / W)
        err = float(np.sqrt(scale ** 2 / W + (dm * slope_err) ** 2))
        note = f"{len(y)} ZTF frames within {pad_days:g} d of the window"
        if src is not None:
            note += f", moved to <r_h> along x = {slope:.2f} ({src['leg']}, grade {src['grade']})"
        if out["n_outburst_window"]:
            note += f"; {out['n_outburst_window']} of them in an outburst"
        out.update(method="direct", log_afrho=ym, log_afrho_err=err, n=int(len(y)),
                   rms_dex=float(np.std(y - ym)) if len(y) > 1 else np.nan,
                   grade=str(src["grade"]) if src is not None else "", note=note)
        return _with_linear(out)

    def _single_or(note):
        if len(win) == 1:
            f = win.iloc[0]
            out.update(method="direct", log_afrho=float(f["log_afrho"]), n=1,
                       log_afrho_err=float(np.hypot(f["log_afrho_err"], 0.1)),
                       note=f"single ZTF frame within {pad_days:g} d of the window ({note})")
        else:
            out["note"] = note
        return _with_linear(out)

    # 2. the fitted law the epoch falls in
    if primary is None:
        if not have_trends or prim.empty:
            return _single_or("no fitted ZTF trend at this aperture")
        fitted = "/".join(sorted(set(prim["leg"])))
        return _single_or(f"ZTF fits only the {fitted} phase; the SPHEREx epoch is {leg}")
    if local is not None:
        y, e = evaluate_trend(local, log_rh0)
        out.update(method="trend", log_afrho=y, log_afrho_err=e, n=int(local["n"]), grade=str(local["grade"]),
                   rms_dex=float(local["rms_dex"]),
                   note=f"{local['leg']} law (n = {int(local['n'])}, grade {local['grade']}) at <r_h>")
        return _with_linear(out)
    # 3. a bounded extrapolation of the primary law
    lo, hi = np.log10(float(primary["rh_min"])), np.log10(float(primary["rh_max"]))
    below = log_rh0 < lo
    dlog = float(lo - log_rh0) if below else float(log_rh0 - hi)
    out["dlog_extrap"] = dlog
    rng = f"{primary['rh_min']:.2f}-{primary['rh_max']:.2f} au"
    if primary["grade"] in extrap_grades and dlog <= max_extrap_dex:
        y, e = evaluate_trend(primary, log_rh0)
        out.update(method="trend_extrap", log_afrho=y, log_afrho_err=e, n=int(primary["n"]),
                   grade=str(primary["grade"]), rms_dex=float(primary["rms_dex"]),
                   note=f"{leg} law extrapolated {dlog:.2f} dex {'inward' if below else 'outward'} of its range {rng}")
        return _with_linear(out)
    if primary["grade"] not in extrap_grades:
        return _single_or(f"{leg} law is unconstrained (grade {primary['grade']}) and <r_h> = {rh0:.2f} au is outside its range {rng}")
    return _single_or(f"<r_h> = {rh0:.2f} au is {dlog:.2f} dex beyond the ZTF range {rng} (limit {max_extrap_dex:g})")
