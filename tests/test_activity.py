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
