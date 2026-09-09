"""
Unit and reproduction tests for ``spherex_comspec``.

Run with::

    python spherex-comspec/tests/test_comspec.py

The reproduction tests need ``data/apphot_revised`` and the old
``results/phase_update_map.csv``; they are skipped when either is absent.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spherex_comspec import directory as _dir  # noqa: E402
from spherex_comspec.config import (ApertureConfig, FitConfig, FlagPolicy, GroupingConfig,  # noqa: E402
                                    ModelParams, Variant, DEFAULT_VARIANTS)
from spherex_comspec.gasmodel import filling_factor, spectrum_mjy  # noqa: E402
from spherex_comspec.grouping import split_groups, band_runs  # noqa: E402
from spherex_comspec.fitting import fit_production_rates, build_design_matrix  # noqa: E402
from spherex_comspec.continuum import process_group, rebuild_continuum  # noqa: E402

HAVE_DATA = (_dir.APPHOT_DIR / "24P.csv").is_file()
OLD_MAP = _dir.ROOT / "data" / "reference" / "phase_update_map_previous.csv"
HAVE_OLD = OLD_MAP.is_file()


# --------------------------------------------------------------------- physics limits
def test_filling_factor_limits():
    assert abs(filling_factor(1e-6) / (np.pi * 1e-6 / 2) - 1) < 1e-4
    assert abs(filling_factor(1e5) - 1.0) < 1e-6
    # "quad" is scipy's adaptive integrator at its default ~1e-8 absolute tolerance, which at
    # x ~ 0.05 is a few 1e-6 relative on int_x^inf K0; it bounds the *check*, not the model.
    x = np.logspace(-3, 1.5, 40)
    assert np.allclose(filling_factor(x, "struve"), filling_factor(x, "quad"), rtol=1e-5)


def test_model_is_linear_in_Q():
    p = ModelParams(rho_ap_km=20000)
    lam = np.linspace(2.5, 5.0, 500)
    a = spectrum_mjy(lam, {"H2O": 1e28, "CO2": 1e27, "CO": 1e26}, 1.5, 1.0, p)["total"]
    b = spectrum_mjy(lam, {"H2O": 3e28, "CO2": 3e27, "CO": 3e26}, 1.5, 1.0, p)["total"]
    assert np.allclose(b, 3 * a)


# --------------------------------------------------------------------------- grouping
def test_split_groups_respects_tolerance_and_prefers_fewest_groups():
    r = np.array([2.00, 1.99, 1.98, 1.70, 1.69, 1.68])
    w = np.ones(6)
    segs = split_groups(r, w, np.zeros(6, int), np.zeros(6), tol=0.10)
    assert segs == [(0, 3), (3, 6)]
    segs = split_groups(r, w, np.zeros(6, int), np.zeros(6), tol=0.30)
    assert segs == [(0, 6)]


def test_band_runs_keep_a_band_whole():
    r = np.array([2.0, 1.98, 1.96, 1.94])
    assert band_runs(r, np.array([True, False, False, True]), 0.10) == [(0, 3)]
    assert band_runs(r, np.array([True, True, False, False]), 0.10) == [(0, 1)]


# ------------------------------------------------------------------- flux-space equivalence
def _synthetic_points(n=12, seed=0):
    rng = np.random.default_rng(seed)
    # every species must clear its KEY_RANGES minimum: CO needs >= 2 channels in 4.60-4.70 um
    wl = np.concatenate([np.linspace(2.62, 2.78, 4), np.linspace(4.20, 4.32, 4),
                         np.array([4.60, 4.64, 4.68, 4.84])])
    band = ["2.7um"] * 4 + ["4.3um"] * 4 + ["4.7um"] * 4
    r_hel = 1.5 + 0.05 * rng.standard_normal(n)
    r_obs = 1.0 + 0.10 * rng.standard_normal(n)
    pts = pd.DataFrame(dict(target="X", r_ap_km=20000.0, phase=1, band=band, wl=wl,
                            wlwidth=wl / 60, r_hel=r_hel, r_obs=r_obs, jd_utc=2460000.0,
                            distcorr_factor=r_hel ** 2 * r_obs ** 2))
    p = ModelParams(rho_ap_km=20000)
    A, _ = build_design_matrix(pts, p, ("H2O", "CO2", "CO"), "physical")
    truth = np.array([1e28, 1e27, 2e26])
    y = A @ truth
    err = 0.05 * np.abs(y).max() * np.ones(n)
    y = y + err * rng.standard_normal(n)
    pts["emis_raw_mjy"], pts["emis_raw_err_mjy"] = y, err
    pts["emis_mjy"] = y * pts.distcorr_factor
    pts["emis_err_mjy"] = err * pts.distcorr_factor
    return pts, p, truth


def test_Q_identical_in_physical_and_distcorr_space():
    pts, p, truth = _synthetic_points()
    cfg = FitConfig(drop_negative_sigma=None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = fit_production_rates(pts, p, cfg, space="physical")
        b = fit_production_rates(pts, p, cfg, space="distcorr")
    assert np.allclose(a.Q_fit, b.Q_fit, rtol=1e-12)
    assert np.allclose(a.Q_err_formal, b.Q_err_formal, rtol=1e-12)
    # Recovery: the two strong species within 25 %; CO -- blended with the H2O hot bands and
    # carried by ~1.6 effective channels in this synthetic -- within its own 3 sigma.
    assert np.allclose(a.Q_fit[:2], truth[:2], rtol=0.25)
    assert np.all(np.abs(a.Q_fit - truth) < 3 * a.Q_err)


def test_fit_reports_a_group_with_no_surviving_channel():
    """An accepted band whose continuum sits above every emission channel must not fit: the
    negative-channel cut empties the point list and the caller records the group as not fitted
    (499P phase 3 of the lenient variant is the real case)."""
    pts, p, _ = _synthetic_points()
    pts["emis_raw_mjy"] = -5.0 * pts["emis_raw_err_mjy"]
    pts["emis_mjy"] = pts["emis_raw_mjy"] * pts.distcorr_factor
    with pytest.raises(ValueError, match="no channel survived"):
        fit_production_rates(pts, p, FitConfig(drop_negative_sigma=1.0))
    assert filling_factor(np.array([])).size == 0


def test_h2o_falls_back_to_hot_bands_only_without_the_main_band():
    """With the 2.7 um channels absent the 4.6-4.9 um hot bands carry Q(H2O), labelled "hot";
    with them present the main band anchors it ("main"); with the fallback off it is not covered."""
    pts, p, _ = _synthetic_points()
    cfg = FitConfig(drop_negative_sigma=None)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = fit_production_rates(pts, p, cfg)
        no27 = fit_production_rates(pts[pts.band != "2.7um"].reset_index(drop=True), p, cfg)
        off = fit_production_rates(pts[pts.band != "2.7um"].reset_index(drop=True), p,
                                   FitConfig(drop_negative_sigma=None, h2o_hot_fallback=False))
    assert full.h2o_source == "main" and full.h2o_anchored
    assert no27.h2o_source == "hot" and not no27.h2o_anchored and no27.covered[0]
    assert np.isfinite(no27.Q_fit[0]) and "hot bands" in "; ".join(no27.caveats())
    assert off.h2o_source == "none" and not off.covered[0]


def test_default_variants_are_well_formed():
    """The driver selects study partners by role, so the registry must have unique names,
    exactly one main variant, and that one must be MAIN_VARIANT."""
    from spherex_comspec.config import DEFAULT_VARIANTS, MAIN_VARIANT, VARIANTS
    names = [v.name for v in DEFAULT_VARIANTS]
    assert len(set(names)) == len(names) and set(VARIANTS) == set(names)
    mains = [v.name for v in DEFAULT_VARIANTS if v.role == "main"]
    assert mains == [MAIN_VARIANT]
    assert len({v.hash for v in DEFAULT_VARIANTS}) == len(names)


# ------------------------------------------------------------------ data-dependent tests
def test_regrouping_reproduces_the_old_map():
    if not (HAVE_DATA and HAVE_OLD):
        return
    from spherex_comspec.grouping import regroup_all
    _, gmap, _, _ = regroup_all(["24P", "2024E1", "10P", "172P"], GroupingConfig())
    old = pd.read_csv(OLD_MAP)
    for t in ("24P", "2024E1", "10P", "172P"):
        o, n = old[old.target == t], gmap[gmap.target == t]
        assert len(o) == len(n), t
        assert np.allclose(o.r_hel_mean.values, n.r_hel_mean.values, rtol=1e-5), t
        assert (o.n_meas.values == n.n_meas.values).all(), t


def test_aperture_promotion_for_a_distant_target():
    if not (_dir.APPHOT_DIR / "2014UN271.csv").is_file():
        return
    from spherex_comspec.dataio import aperture_for
    r, lab, cov = aperture_for("2014UN271", ApertureConfig())
    assert r == 80000.0 and cov >= 0.95
    r, lab, cov = aperture_for("24P", ApertureConfig())
    assert r == 20000.0 and cov >= 0.99


def test_flag_policies_select_nested_samples():
    if not HAVE_DATA:
        return
    from spherex_comspec.dataio import load_apphot, select_spectrum
    df = load_apphot("24P")
    n = {}
    for name, fl in (("all", ()), ("no_a", ("a",)), ("no_ab", ("a", "b"))):
        v = Variant(name, FlagPolicy(name, drop_flags=fl))
        n[name] = len(select_spectrum(df, "km20000", v))
    assert n["all"] >= n["no_a"] >= n["no_ab"]
    lenient = Variant("l", FlagPolicy("l", max_frac_badpix=0.05))
    assert len(select_spectrum(df, "km20000", lenient)) > n["all"]


def test_continuum_round_trips_from_its_saved_columns():
    if not HAVE_DATA:
        return
    from spherex_comspec.config import ContinuumConfig
    from spherex_comspec.dataio import load_apphot, select_spectrum
    from spherex_comspec.grouping import regroup_all
    assign, _, _, _ = regroup_all(["24P"])
    v = DEFAULT_VARIANTS[0]
    spec = select_spectrum(load_apphot("24P"), "km20000", v, assign)
    raw = spec[spec.phase == 2].reset_index(drop=True)
    out = process_group(raw, "24P", 20000, 2, int(raw.epoch.iloc[0]), ContinuumConfig(), "distcorr")
    for _, row in out["summary"].iterrows():
        if not np.isfinite(row.cont_c0):
            continue
        d = out["points"][(out["points"].band == row.band)]
        model, sig = rebuild_continuum(row, d.wl.to_numpy())
        assert np.allclose(model, d.cont_mjy, rtol=1e-8, atol=1e-12)
        assert np.allclose(sig, d.cont_err_mjy, rtol=1e-8, atol=1e-12)
        # physical-space emission is the corrected one divided by the per-channel factor
        assert np.allclose(d.emis_raw_mjy * d.distcorr_factor, d.emis_mjy, rtol=1e-10)


if __name__ == "__main__":
    import traceback
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:
            failed += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)
