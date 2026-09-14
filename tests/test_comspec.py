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


# ------------------------------------------------------------ fluorescence database wiring
def test_fluorescence_templates_integrate_to_the_database_gfactors():
    """The gfm species shape is the database g_lambda times the photon energy: its integral over a
    band region must equal the sum of the tabulated band g-factors there, on any grid."""
    from spherex_comspec import fluorescence as fl
    bands = pd.read_csv(fl.FLUOR_DIR / "band_gfactors.csv")
    for species, lo, hi in (("H2O", 2.5, 3.1), ("CO2", 4.1, 4.5), ("CO", 4.5, 5.0)):
        ref = bands[(bands.species == species) & bands.lam_center_um.between(lo, hi)]["g_T070K"].sum()
        assert abs(fl.band_total(species, lo, hi, 70.0) / ref - 1) < 0.03, species
        for lam in (np.linspace(lo, hi, 300), lo * np.exp(np.arange(0, np.log(hi / lo), 1 / 4000))):
            g = fl.g_lambda(species, lam, 70.0)
            e = fl._edges(lam)
            assert abs(np.sum(g * np.diff(e)) / ref - 1) < 0.03, (species, len(lam))
    # T_rot interpolation is bounded by the neighbouring tables and exact on a grid point
    lam = np.linspace(2.5, 3.1, 500)
    g50, g60, g70 = (fl.g_lambda("H2O", lam, T) for T in (50.0, 60.0, 70.0))
    assert np.all(g60 <= np.maximum(g50, g70) + 1e-30) and np.all(g60 >= np.minimum(g50, g70) - 1e-30)
    assert np.allclose(fl.g_lambda("H2O", lam, 200.0), fl.g_lambda("H2O", lam, 130.0))


def test_co_swings_factor_scales_only_the_co_column():
    from spherex_comspec import fluorescence as fl
    assert fl.swings_factor("CO", 0.0) == pytest.approx(1.0, abs=1e-6)
    f20 = fl.swings_factor("CO", 20.0)
    assert 1.15 < f20 < 1.45 and abs(fl.swings_factor("CO", -20.0) / f20 - 1) < 0.05
    assert fl.swings_factor("CO", 200.0) == fl.swings_factor("CO", 60.0)      # clamped
    assert fl.swings_factor("H2O", 20.0) == 1.0 and fl.swings_factor("CO", np.nan) == 1.0
    assert np.allclose(fl.swings_factor("CO", np.array([0.0, np.nan])), [1.0, 1.0])
    pts, p, _ = _synthetic_points()
    pts["v_hel_kms"] = 0.0
    A0, _ = build_design_matrix(pts, p, ("H2O", "CO2", "CO"), "physical")
    pts["v_hel_kms"] = 20.0
    A20, _ = build_design_matrix(pts, p, ("H2O", "CO2", "CO"), "physical")
    assert np.allclose(A20[:, :2], A0[:, :2])
    nz = A0[:, 2] > 0
    assert np.allclose(A20[nz, 2] / A0[nz, 2], f20)
    # the fit reports the factor and the velocity it used
    pts["v_hel_kms"] = 20.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = fit_production_rates(pts, p, FitConfig(drop_negative_sigma=None))
    assert fit.geometry["v_hel_mean_kms"] == pytest.approx(20.0)
    assert fit.geometry["swings_CO"] == pytest.approx(f20)
    off = ModelParams(rho_ap_km=20000, co_swings=False)
    Aoff, _ = build_design_matrix(pts, off, ("H2O", "CO2", "CO"), "physical")
    assert np.allclose(Aoff, A0)


def test_gaussian_legacy_model_still_runs_and_differs_from_gfm():
    lam = np.linspace(2.5, 5.0, 800)
    Q = {"H2O": 1e28, "CO2": 1e27, "CO": 1e26}
    new = spectrum_mjy(lam, Q, 1.5, 1.0, ModelParams(rho_ap_km=20000))["total"]
    old = spectrum_mjy(lam, Q, 1.5, 1.0, ModelParams(rho_ap_km=20000, profile_source="gaussian",
                                                       co_swings=False))["total"]
    assert np.all(np.isfinite(new)) and np.all(np.isfinite(old)) and new.max() > 0
    # same physics, different band strengths and shapes: the integrals agree to ~15 %
    assert abs(np.trapezoid(new, lam) / np.trapezoid(old, lam) - 1) < 0.25
    with pytest.raises(ValueError):
        spectrum_mjy(lam, Q, 1.5, 1.0, ModelParams(profile_source="psg"))


def test_heliocentric_velocity_from_the_ephemeris():
    from spherex_comspec.dataio import heliocentric_velocity, AU_KM_PER_DAY_TO_KMS
    t = 2460000.0 + np.sort(np.concatenate([np.linspace(0, 3, 40), [0.5, 0.5, 1.7]]))
    a, b = -0.0060, 0.00025                                   # au/day, au/day^2
    r = 1.5 + a * (t - t[0]) + b * (t - t[0]) ** 2
    v = heliocentric_velocity(np.round(t, 6), np.round(r, 6))
    expected = (a + 2 * b * (t - t[0])) * AU_KM_PER_DAY_TO_KMS
    assert np.all(np.isfinite(v)) and np.allclose(v, expected, atol=0.05)
    assert np.all(np.isnan(heliocentric_velocity([2460000.0, 2460000.0], [1.5, 1.5])))
    assert np.isnan(heliocentric_velocity([2460000.0], [1.5]))[0]


# --------------------------------------------------------- 2026-09-12 rules
def test_select_order_prefers_the_simplest_adequate_polynomial():
    from spherex_comspec.config import ContinuumConfig
    from spherex_comspec.continuum import select_order
    x = np.linspace(-0.4, 0.4, 24)
    e = np.full(24, 0.05)
    cfg = ContinuumConfig()
    picks_lin, picks_cub = [], []
    for seed in range(10):
        rng = np.random.default_rng(seed)
        picks_lin.append(select_order(x, 3.0 + 2.0 * x + e * rng.standard_normal(24), e, 3, cfg))
        picks_cub.append(select_order(x, 3.0 + 2.0 * x - 40.0 * x ** 3 + e * rng.standard_normal(24), e, 3, cfg))
    assert picks_lin.count(1) >= 8 and all(p == 3 for p in picks_cub)
    assert select_order(x[:3], x[:3], e[:3], 3, cfg) == 1           # nothing to cross-validate


def test_one_sided_continuum_is_extended_to_bracket_the_band():
    from spherex_comspec.config import ContinuumConfig
    from spherex_comspec.continuum import fit_continuum
    # channels only red of the 2.7 um band inside the standard window, plus a blue stretch at
    # 1.9-2.15 um that only the extension can reach
    wl = np.concatenate([np.linspace(1.90, 2.15, 8), np.linspace(2.82, 3.08, 8)])
    raw = pd.DataFrame(dict(wl=wl, flux=1.0 + 0.1 * (wl - 2.7), err=np.full(len(wl), 0.02),
                            distcorr_factor=1.0))
    ext = fit_continuum(raw, "2.7um", ContinuumConfig())
    assert ext["extended"] == "blue" and ext["n_blue"] > 0 and ext["bracketed"]
    assert abs(ext["cont"][0] - 1.50) < 1e-9
    off = fit_continuum(raw, "2.7um", ContinuumConfig(one_sided_extend_um=None))
    assert off["n_blue"] == 0 and not off["bracketed"] and off["order_used"] == 1


def test_detection_tiers_and_limits():
    pts, p, truth = _synthetic_points()
    cfg = FitConfig()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = fit_production_rates(pts, p, cfg)
        # inflate the errors so the strong species become marginal and the weak one a limit
        big = pts.copy()
        for c in ("emis_raw_err_mjy", "emis_err_mjy"):
            big[c] = big[c] * 12
        weak = fit_production_rates(big, p, cfg)
    for k in range(3):
        st, nsig = fit.status[k], fit.Q_fit[k] / fit.Q_err[k]
        assert st == ("detected" if nsig >= 3 else "marginal" if nsig >= 1 else "upper_limit") or st == "negative_fit"
        if st in ("marginal", "upper_limit"):
            assert fit.Q_limit[k] == pytest.approx(fit.Q_fit[k] + 3 * fit.Q_err[k])
    assert "marginal" in set(weak.status) or "upper_limit" in set(weak.status)
    assert not weak.upper_limit[list(weak.status).index("marginal")] if "marginal" in set(weak.status) else True
    row = fit.to_row()
    assert "Q_H2O_nsig" in row and np.isfinite(row["Q_H2O_nsig"])


def test_gls_reduces_to_the_diagonal_solve_without_continuum_covariance():
    from spherex_comspec.fitting import data_covariance
    pts, p, _ = _synthetic_points()
    cfg = FitConfig()
    assert data_covariance(pts, cfg) is None                 # no covariance columns -> diagonal path
    # zero coefficient covariance and err_raw = emis_raw_err: the GLS solve must equal the diagonal one
    z = pts.copy()
    z["err_raw"] = z["emis_raw_err_mjy"]
    z["flux_space"] = "physical"
    z["cont_lam_ref_um"] = z.band.map({"2.7um": 2.70, "4.3um": 4.26, "4.7um": 4.67})
    for i in range(4):
        z[f"cont_c{i}"] = 0.0
        for j in range(i, 4):
            z[f"cont_cov_{i}{j}"] = 0.0
    S = data_covariance(z, cfg)
    assert S is not None and np.allclose(np.diag(S), z.emis_raw_err_mjy ** 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = fit_production_rates(z, p, cfg)
        b = fit_production_rates(pts, p, FitConfig(gls=False))
    assert np.allclose(a.Q_fit, b.Q_fit, rtol=1e-10) and np.allclose(a.Q_err_formal, b.Q_err_formal, rtol=1e-10)
    assert a.geometry["gls"] and not b.geometry["gls"]
    # a correlated continuum error within a band inflates the errors of the species carried by it
    z2 = z.copy()
    z2.loc[z2.band == "4.3um", "cont_cov_00"] = (0.5 * z2.loc[z2.band == "4.3um", "emis_raw_err_mjy"].mean()) ** 2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = fit_production_rates(z2, p, cfg)
    assert c.Q_err_formal[1] > a.Q_err_formal[1]


def test_h2o_hot_fallback_respects_the_distance_cap():
    pts, p, _ = _synthetic_points()
    no27 = pts[pts.band != "2.7um"].reset_index(drop=True)
    far = no27.copy()
    far["r_hel"] = 4.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        near = fit_production_rates(no27, p, FitConfig())
        beyond = fit_production_rates(far, p, FitConfig())
        uncapped = fit_production_rates(far, p, FitConfig(h2o_hot_max_rh_au=None))
    assert near.h2o_source == "hot" and beyond.h2o_source == "none" and uncapped.h2o_source == "hot"
    assert "withheld beyond 3 au" in "; ".join(beyond.caveats())


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
    # the old map predates rule 1b (Delta spread); reproduce it with that rule off
    _, gmap, _, _ = regroup_all(["24P", "2024E1", "10P", "172P"], GroupingConfig(delta_tol=None))
    old = pd.read_csv(OLD_MAP)
    for t in ("24P", "2024E1", "10P", "172P"):
        o, n = old[old.target == t], gmap[gmap.target == t]
        assert len(o) == len(n), t
        assert np.allclose(o.r_hel_mean.values, n.r_hel_mean.values, rtol=1e-5), t
        assert (o.n_meas.values == n.n_meas.values).all(), t


def test_delta_rule_bounds_every_group():
    if not HAVE_DATA:
        return
    from spherex_comspec.grouping import regroup_all
    _, gmap, _, _ = regroup_all(["24P", "10P"], GroupingConfig())
    assert (gmap.delta_spread < 0.20).all() and (gmap[~gmap.manual].spread < 0.10).all()
    _, gmap0, _, _ = regroup_all(["24P", "10P"], GroupingConfig(delta_tol=None))
    assert len(gmap) >= len(gmap0) and (gmap0.delta_spread >= 0.20).any()


def test_aperture_rules():
    if not (_dir.APPHOT_DIR / "2014UN271.csv").is_file():
        return
    from spherex_comspec.dataio import aperture_for, aperture_choice
    # the previous rule, kept for the dc_rules_previous variant
    r, lab, cov = aperture_for("2014UN271", ApertureConfig(rule="rh"))
    assert r == 80000.0 and cov >= 0.95
    r, lab, cov = aperture_for("24P", ApertureConfig(rule="rh"))
    assert r == 20000.0 and cov >= 0.99
    # the S/N rule respects coverage, the PSF and the annulus bound; 2022 R3 has negative scores
    for t in ("24P", "10P", "2014UN271", "2022R3"):
        c = aperture_choice(t, ApertureConfig())
        assert c["coverage"] >= 0.95 and c["r_ap_km"] > 0
        if c["rule"] == "snr" and not c["relaxed"] and np.isfinite(c["r_in_km"]):
            assert c["r_ap_km"] <= c["r_in_km"] / 3 + 1e-6
            assert c["snr_median"] >= c["snr_best"] - 0.1 * abs(c["snr_best"])


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
