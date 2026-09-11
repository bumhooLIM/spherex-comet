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
    residual scatter), ``dlog_rh`` (baseline), ``median_log_err``.  All fit
    fields are NaN with fewer than three points or no baseline.
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
               median_log_err=float(np.median(e)) if n else np.nan)
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
               rms_dex=float(np.sqrt(np.mean(resid ** 2))))
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
    out = dict(t_peak=np.nan, t_lo=np.nan, t_hi=np.nan, bracketed=False, kernel_days=np.nan,
               n=n, peak_log_afrho=np.nan, rms_dex=np.nan, drop_lo=np.nan, drop_hi=np.nan)
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
    interior = (t.min() + edge_frac * span) < grid[i] < (t.max() - edge_frac * span)
    drop_lo, drop_hi = float(s[i] - s[0]), float(s[i] - s[-1])
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
               kernel_days=h, peak_log_afrho=float(s[i]), rms_dex=rms, drop_lo=drop_lo, drop_hi=drop_hi)
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


# --------------------------------------------------------------- per target
def analyse(phot, target, tp_jd=None, bands=("r", "g"), rhos=(10000, 20000),
            min_side=3, n_boot=1000):
    """Trends and peak for one comet.  Returns ``(trends, peaks)`` DataFrames.

    For every band and aperture: the whole-data fit, the perihelion-split
    fits (inbound / outbound by the sign of ``r_rate``), and -- when both
    sides have at least *min_side* points and ``tp_jd`` is known -- the
    activity peak and the fits split at the peak instead.  One-sided comets
    get the single-leg fit only, as the peak is then not measurable.
    """
    trends, peaks = [], []
    for band in bands:
        for rho in rhos:
            pts, counts = select_points(phot, band, rho)
            base = dict(target=target, band=band, rho_km=rho, **counts)
            if len(pts) == 0:
                continue
            legs = [("all", "all", pts)]
            inb, outb = pts[pts["leg"] == "inbound"], pts[pts["leg"] == "outbound"]
            two_sided = len(inb) >= min_side and len(outb) >= min_side
            if len(inb) >= 3:
                legs.append(("perihelion", "inbound", inb))
            if len(outb) >= 3:
                legs.append(("perihelion", "outbound", outb))
            pk = None
            if two_sided and tp_jd is not None and np.isfinite(tp_jd):
                t = pts["obsjd"].to_numpy() - tp_jd
                pk = find_peak(t, pts["log_afrho"], pts["log_afrho_err"], n_boot=max(n_boot // 2, 200))
                peaks.append(dict(target=target, band=band, rho_km=rho, n_inbound=len(inb),
                                  n_outbound=len(outb), t_min=float(t.min()), t_max=float(t.max()), **pk))
                if pk["bracketed"]:
                    before, after = pts[t < pk["t_peak"]], pts[t >= pk["t_peak"]]
                    if len(before) >= 3:
                        legs.append(("peak", "rising", before))
                    if len(after) >= 3:
                        legs.append(("peak", "fading", after))
            for split, leg, sub in legs:
                fit = fit_powerlaw(sub["log_rh"], sub["log_afrho"], sub["log_afrho_err"], n_boot=n_boot)
                g, why = grade_fit(fit, counts)
                trends.append(dict(**base, split=split, leg=leg, rh_min=float(sub["r"].min()),
                                   rh_max=float(sub["r"].max()), two_sided=two_sided,
                                   t_peak=(pk["t_peak"] if pk else np.nan),
                                   grade=g, reasons="; ".join(why), **fit))
    return pd.DataFrame(trends), pd.DataFrame(peaks)
