"""Heliocentric Af-rho trends: fit, peak and grade on synthetic data."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ztfcomet import activity as ac


def _series(x_true=2.5, a=3.0, n=40, err=0.03, seed=0, rh=(1.2, 3.0)):
    rng = np.random.default_rng(seed)
    log_rh = np.log10(rng.uniform(*rh, n))
    y = a - x_true * log_rh + rng.normal(0, err, n)
    return log_rh, y, np.full(n, err)


def test_powerlaw_recovers_slope_within_errors():
    lr, y, e = _series()
    f = ac.fit_powerlaw(lr, y, e)
    assert abs(f["x"] - 2.5) < 3 * f["x_err"]
    assert f["x_boot_lo"] < 2.5 < f["x_boot_hi"]
    assert 0.5 < f["chi2_red"] < 2.0 and f["n"] == 40


def test_scaled_error_grows_with_excess_scatter():
    lr, y, e = _series(err=0.15)
    f = ac.fit_powerlaw(lr, y, e / 5)                    # errors understated 5x
    assert f["chi2_red"] > 10
    assert f["x_err_scaled"] > 3 * f["x_err"]


def test_no_baseline_means_no_slope():
    y = np.array([3.0, 3.1, 2.9, 3.0])
    f = ac.fit_powerlaw(np.zeros(4), y, np.full(4, 0.05))
    assert np.isnan(f["x"]) and f["n"] == 4


def test_peak_is_found_and_bracketed():
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(-120, 150, 60))
    y = 2.5 - 0.5 * ((t - 25.0) / 60.0) ** 2 + rng.normal(0, 0.03, 60)   # peak at +25 d
    pk = ac.find_peak(t, y, np.full(60, 0.03))
    assert pk["bracketed"]
    assert abs(pk["t_peak"] - 25.0) < 10
    assert pk["t_lo"] <= pk["t_peak"] <= pk["t_hi"]


def test_monotonic_series_has_no_bracketed_peak():
    t = np.linspace(-100, 40, 30)
    y = 2.0 + 0.004 * t
    pk = ac.find_peak(t, y, np.full(30, 0.03))
    assert not pk["bracketed"]


def test_grade_leads_with_the_baseline():
    fit_wide = dict(x=2.0, n=30, dlog_rh=0.3, x_err_scaled=0.2, chi2_red=1.2)
    fit_narrow = dict(x=2.0, n=30, dlog_rh=0.02, x_err_scaled=3.0, chi2_red=1.2)
    assert ac.grade_fit(fit_wide)[0] == "A"
    assert ac.grade_fit(fit_narrow)[0] == "D"
    g, why = ac.grade_fit(fit_wide, dict(n_all=100, n_contaminated=60, n_offcentre=0))
    assert g == "B" and any("rejected" in r for r in why)


def test_select_points_drops_offcentre_and_counts_contamination():
    n = 10
    d = pd.DataFrame(dict(rho_km=10000.0, filter="ZTF_r", quality_ok=True, afrho0_cm=100.0,
                          afrho0_cm_err=5.0, r=2.0, r_rate=-1.0, rho_pix=6.0,
                          centroid_shift_pix=np.r_[np.full(8, 0.5), 4.0, 5.0],
                          flag_contaminated=np.r_[np.full(7, False), True, False, False],
                          obsjd=np.arange(n) + 2.46e6))
    kept, counts = ac.select_points(d, "r", 10000)
    assert counts == dict(n_all=10, n_clean=10, n_contaminated=1, n_offcentre=2)
    assert len(kept) == 8 and (kept["leg"] == "inbound").all()


def _broken_series(x1=2.0, x2=5.0, rb=3.0, n=60, err=0.03, seed=4, rh=(1.5, 6.0)):
    rng = np.random.default_rng(seed)
    lr = np.log10(rng.uniform(*rh, n)); xb = np.log10(rb)
    y = 3.0 - x1 * (lr - xb) - (x2 - x1) * np.maximum(lr - xb, 0) + rng.normal(0, err, n)
    return lr, y, np.full(n, err)


def test_broken_powerlaw_finds_the_break_and_both_indices():
    lr, y, e = _broken_series()
    f = ac.fit_broken_powerlaw(lr, y, e)
    assert f["testable"] and f["preferred"]
    assert abs(f["r_break"] - 3.0) < 0.6 and f["r_break_lo"] <= f["r_break"] <= f["r_break_hi"]
    assert abs(f["x_inner"] - 2.0) < 0.4 and abs(f["x_outer"] - 5.0) < 0.5


def test_single_law_is_not_broken():
    lr, y, e = _series(n=60)
    f = ac.fit_broken_powerlaw(lr, y, e)
    assert f["testable"] and not f["preferred"] and f["dbic"] < 6


def test_fixed_break_at_3_au_gives_the_two_indices():
    lr, y, e = _broken_series()
    f = ac.fit_broken_powerlaw(lr, y, e, fixed_xb=np.log10(3.0), n_boot=0)
    assert f["testable"] and abs(f["x_inner"] - 2.0) < 0.4 and abs(f["x_outer"] - 5.0) < 0.5


def test_colour_pairs_and_change_point():
    rng = np.random.default_rng(5)
    n = 40
    jd = 2.46e6 + np.sort(rng.uniform(0, 300, n)); rh = np.linspace(1.5, 5.0, n)
    afr = 500 * rh ** -2.0
    excess_true = np.where(rh < 3.0, 0.10, 0.35)                  # redder beyond 3 au
    afg = afr * 10 ** (-0.4 * excess_true)
    rows = []
    for i in range(n):
        for band, a in (("ZTF_r", afr[i]), ("ZTF_g", afg[i])):
            rows.append(dict(obsjd=jd[i] + (0.01 if band == "ZTF_g" else 0.0), filter=band, rho_km=10000.0,
                             quality_ok=True, afrho0_cm=a * (1 + rng.normal(0, 0.02)), afrho0_cm_err=0.02 * a,
                             r=rh[i], r_rate=1.0, rho_pix=6.0, centroid_shift_pix=0.3, flag_contaminated=False,
                             file=f"{band}_{i}", filter_mag=np.nan, filter_mag_err=np.nan))
    phot = pd.DataFrame(rows)
    pairs = ac.colour_pairs(phot, 10000)
    assert len(pairs) == n and (pairs["dt_days"] < 0.02).all()
    ct = ac.fit_colour_trend(pairs["log_rh"], pairs["excess"], pairs["excess_err"])
    assert ct["change_preferred"] and abs(ct["r_change"] - 3.0) < 0.5
    assert abs(ct["delta_colour"] - 0.25) < 0.05 and ct["delta_colour"] > 5 * ct["delta_colour_err"]


def test_analyse_uses_the_peak_to_define_legs():
    rng = np.random.default_rng(6)
    n = 80
    t = np.sort(rng.uniform(-150, 150, n)); tp = 2.46e6
    rh = 1.2 + (t / 100.0) ** 2 * 0.8                              # symmetric orbit, q = 1.2 au
    logafr = 2.5 - 0.5 * ((t - 30.0) / 70.0) ** 2                 # peak 30 d after perihelion
    rows = [dict(obsjd=tp + t[i], filter="ZTF_r", rho_km=10000.0, quality_ok=True,
                 afrho0_cm=10 ** logafr[i] * (1 + rng.normal(0, 0.03)), afrho0_cm_err=0.03 * 10 ** logafr[i],
                 r=rh[i], r_rate=np.sign(t[i]) or 1.0, rho_pix=6.0, centroid_shift_pix=0.3,
                 flag_contaminated=False, file=f"f{i}") for i in range(n)]
    trends, peaks, breaks, colours, outbursts = ac.analyse(pd.DataFrame(rows), "X", tp_jd=tp, bands=("r",), rhos=(10000,), n_boot=200)
    assert peaks.iloc[0]["bracketed"] and abs(peaks.iloc[0]["t_peak"] - 30.0) < 12
    prim = trends[trends["primary"]]
    assert set(prim[prim["split"] == "peak"]["leg"]) == {"rising", "fading"}
    assert set(prim["split"]) <= {"peak", "segment"}          # segments of a preferred break are primary too


def test_understated_errors_do_not_manufacture_a_break():
    lr, y, e = _series(n=60, err=0.15)
    f = ac.fit_broken_powerlaw(lr, y, e / 10)                   # errors 10x too small
    assert f["testable"] and not f["preferred"]
    assert f["dbic"] < 6


def test_sparse_tail_marks_a_small_isolated_group():
    x = np.log10(np.r_[np.linspace(1.4, 2.6, 56), [3.58, 3.60, 3.62, 3.63]])   # 10P-like
    tail = ac.sparse_tail(x)
    assert tail.sum() == 4 and tail[-4:].all()
    assert ac.sparse_tail(np.log10(np.linspace(1.4, 3.6, 60))).sum() == 0    # no gap, no tail


def test_outburst_is_detected_and_windowed():
    rng = np.random.default_rng(7)
    t = np.arange(0, 120, 2.0); y = 1.5 - 0.002 * t + rng.normal(0, 0.03, len(t))
    on = 30                                                    # jump x5 at t=60, decay over ~40 d
    y[on:] += 0.7 * np.exp(-(t[on:] - t[on]) / 20.0)
    win, order = ac.detect_outbursts(t, y)
    assert len(win) == 1 and abs(win[0]["t_start"] - 60.0) < 2.1 and win[0]["ended"]
    mask, _ = ac.outburst_mask(t, y)
    assert mask[on] and not mask[on - 1] and 10 < mask.sum() < 40
    assert ac.detect_outbursts(t, 1.5 - 0.002 * t + rng.normal(0, 0.03, len(t)))[0] == []
    # a steep but smooth rise (0.5 dex over 100 d) is a trend, not an outburst
    tt = np.linspace(-150, 30, 60); smooth = 2.5 - 0.5 * ((tt - 30.0) / 70.0) ** 2 + rng.normal(0, 0.02, 60)
    assert ac.detect_outbursts(tt, smooth)[0] == []
    # a jump that keeps rising is an onset, not an outburst
    onset = np.r_[1.5 + rng.normal(0, 0.02, 20), 2.1 + 0.01 * np.arange(20) + rng.normal(0, 0.02, 20)]
    assert ac.detect_outbursts(np.arange(40) * 2.0, onset)[0] == []


def test_smooth_trend_follows_the_data_and_stops_with_it():
    lr, y, e = _series(n=80)
    g, s = ac.smooth_trend(lr, y, e)
    assert len(g) > 50 and np.all(np.isfinite(s))
    assert g.min() >= lr.min() - 1e-9 and g.max() <= lr.max() + 1e-9
    assert abs(np.polyfit(g, s, 1)[0] + 2.5) < 0.3                 # tracks the underlying slope


def test_plateau_peak_still_splits_the_phases():
    rng = np.random.default_rng(8); n = 70
    t = np.sort(rng.uniform(-150, 150, n)); tp = 2.46e6
    rh = 1.2 + (t / 100.0) ** 2 * 0.8
    # a plateau tilted 0.03 dex toward the peak (a truly flat one has no defined maximum), then fading
    logafr = np.where(t < -20, 2.47 + 0.03 * (t + 150) / 130.0, 2.5 - 0.6 * ((t + 20) / 120.0)) + rng.normal(0, 0.02, n)
    rows = [dict(obsjd=tp + t[i], filter="ZTF_r", rho_km=10000.0, quality_ok=True,
                 afrho0_cm=10 ** logafr[i], afrho0_cm_err=0.02 * 10 ** logafr[i], r=rh[i],
                 r_rate=np.sign(t[i]) or 1.0, rho_pix=6.0, centroid_shift_pix=0.3,
                 flag_contaminated=False, file=f"f{i}") for i in range(n)]
    trends, peaks, *_ = ac.analyse(pd.DataFrame(rows), "P", tp_jd=tp, bands=("r",), rhos=(10000,), n_boot=200)
    pk = peaks.iloc[0]
    assert pk["interior"] and not pk["bracketed"]
    assert set(trends[trends["primary"] & (trends["split"] == "peak")]["leg"]) == {"rising", "fading"}


# ------------------------------------------------------------ epoch value
def _epoch_setup(x_true=2.5, a=3.0, tp=2461000.0):
    """A one-leg inbound series over 1.2-3 au and its trend rows, as
    ``analyse`` would write them."""
    lr, y, e = _series(x_true, a, n=40, err=0.03)
    rh = 10 ** lr
    # inbound: r_h falls with time, 150 days before perihelion at 1.2 au
    jd = tp - 150.0 - 300.0 * (lr - lr.min()) / np.ptp(lr)
    pts = pd.DataFrame(dict(obsjd=jd, r=rh, log_rh=lr, log_afrho=y, log_afrho_err=e, leg="inbound"))
    fit = ac.fit_powerlaw(lr, y, e, n_boot=100)
    row = dict(target="T", band="r", rho_km=10000, split="perihelion", leg="inbound", primary=True,
               two_sided=False, grade="A", segment_of="", has_segments=False,
               rh_min=float(rh.min()), rh_max=float(rh.max()), t_peak=np.nan, **fit)
    return pts, pd.DataFrame([row]), tp


def test_epoch_value_from_the_trend_inside_the_range():
    pts, tr, tp = _epoch_setup()
    rh0 = 2.0
    far = pts["obsjd"].min() - 400.0                 # no frames anywhere near this epoch
    r = ac.afrho_at_epoch(pts, tr, rh0, far, far, far + 10, tp_jd=tp)
    assert r["method"] == "trend" and r["leg"] == "inbound"
    assert abs(r["log_afrho"] - (3.0 - 2.5 * np.log10(rh0))) < 3 * r["log_afrho_err"]
    assert r["afrho_lo_cm"] < r["afrho_cm"] < r["afrho_hi_cm"]


def test_epoch_prefers_the_frames_inside_the_window():
    pts, tr, tp = _epoch_setup()
    jd0 = float(pts["obsjd"].iloc[20])
    rh0 = float(pts["r"].iloc[20])
    # the comet is 0.3 dex above its own law during this window (an outburst);
    # the estimator pads the window by 5 d, so lift every frame it will see
    up = pts.copy()
    win = (up["obsjd"] >= jd0 - 8) & (up["obsjd"] <= jd0 + 8)
    assert win.sum() >= 2
    up.loc[win, "log_afrho"] += 0.3
    r = ac.afrho_at_epoch(up, tr, rh0, jd0, jd0 - 3, jd0 + 3, tp_jd=tp, pad_days=5.0)
    assert r["method"] == "direct" and r["n"] == int(win.sum())
    assert abs(r["log_afrho"] - (3.0 - 2.5 * np.log10(rh0) + 0.3)) < 0.05
    assert "moved to <r_h>" in r["note"]


def test_epoch_extrapolation_is_bounded_and_widens_the_error():
    pts, tr, tp = _epoch_setup()
    far = pts["obsjd"].min() - 400.0
    inside = ac.afrho_at_epoch(pts, tr, 2.0, far, tp_jd=tp)
    rh_max = float(tr["rh_max"].iloc[0])
    near = ac.afrho_at_epoch(pts, tr, rh_max * 10 ** 0.05, far, tp_jd=tp)   # 0.05 dex beyond
    assert near["method"] == "trend_extrap" and abs(near["dlog_extrap"] - 0.05) < 1e-6
    assert near["log_afrho_err"] > inside["log_afrho_err"]
    gone = ac.afrho_at_epoch(pts, tr, rh_max * 10 ** 0.3, far, tp_jd=tp)    # 0.3 dex beyond
    assert gone["method"] == "none" and "beyond the ZTF range" in gone["note"]
    assert np.isnan(gone["afrho_cm"])


def test_epoch_on_the_unfitted_leg_gives_nothing():
    pts, tr, tp = _epoch_setup()
    r = ac.afrho_at_epoch(pts, tr, 2.0, tp + 100.0, tp_jd=tp)              # outbound, ZTF fitted inbound only
    assert r["method"] == "none" and r["leg"] == "outbound"
    assert "fits only the inbound phase" in r["note"]


def test_epoch_uses_a_grade_d_law_only_inside_its_data():
    pts, tr, tp = _epoch_setup()
    tr = tr.assign(grade="D")
    far = pts["obsjd"].min() - 400.0
    assert ac.afrho_at_epoch(pts, tr, 2.0, far, tp_jd=tp)["method"] == "trend"
    out = ac.afrho_at_epoch(pts, tr, float(tr["rh_max"].iloc[0]) * 10 ** 0.05, far, tp_jd=tp)
    assert out["method"] == "none" and "unconstrained" in out["note"]


def test_epoch_without_frames_or_trends():
    r = ac.afrho_at_epoch(pd.DataFrame(), pd.DataFrame(), 2.0, 2461000.0)
    assert r["method"] == "none" and "no clean ZTF" in r["note"]


def test_evaluate_trend_error_is_smallest_at_the_data():
    lr, y, e = _series()
    f = ac.fit_powerlaw(lr, y, e, n_boot=100)
    _, e_mid = ac.evaluate_trend(f, f["log_rh_mean"], include_scatter=False)
    _, e_edge = ac.evaluate_trend(f, lr.max(), include_scatter=False)
    _, e_sc = ac.evaluate_trend(f, f["log_rh_mean"])
    assert e_mid < e_edge and e_sc > e_mid
