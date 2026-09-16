"""
Unit tests for the pieces that are cheap to test and easy to get wrong.

The review noted that the prototype had no tests at all, while carrying exactly
the kind of arithmetic -- bit masks, unit conversions, error budgets -- where a
silent sign or factor error is both easy to introduce and invisible in the
output.  These run in a second or two and need no data files.

Run with (no pytest needed -- the module has its own runner)::

    python tests/test_apphot_pipeline.py

or, if pytest is installed::

    python -m pytest tests/test_apphot_pipeline.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from spherex_apphot.apertures import build_apertures, pixel_scale_km  # noqa: E402
from spherex_apphot.config import Config  # noqa: E402
from spherex_apphot.epochs import group_epochs  # noqa: E402
from spherex_apphot.masking import build_badpix_mask, flag_to_mask  # noqa: E402
from spherex_apphot.phot import estimate_sky, measure_exposure  # noqa: E402
from spherex_apphot.reflectance import add_reflectance, bin_spectrum  # noqa: E402
from spherex_apphot.diagnostics import (  # noqa: E402
    flag_effectiveness, growth_metrics, worst_flag,
)
from spherex_apphot.phot import add_distance_corrected_flux  # noqa: E402
from spherex_apphot.sourceflag import (  # noqa: E402
    GaiaPixelField, effective_gmag, evaluate_sourceflag, mag_to_flux,
)
from spherex_apphot.stacking import StackInput, extract_stamp, shift_subpixel, stack_frames  # noqa: E402
from spherex_apphot.status import slugify  # noqa: E402


# --------------------------------------------------------------------- masking
def test_flag_to_mask_selects_only_requested_bits():
    flag = np.zeros((3, 3), dtype=np.int32)
    flag[0, 0] = 1 << 2
    flag[0, 1] = 1 << 21          # SOURCE: never masked, or the comet vanishes
    flag[0, 2] = (1 << 2) | (1 << 21)
    m = flag_to_mask(flag, [2, 6, 7])
    assert m[0, 0] and not m[0, 1] and m[0, 2]
    assert m.sum() == 2


def test_flag_to_mask_handles_high_bit_without_sign_extension():
    flag = np.array([[-1]], dtype=np.int32)          # every bit set
    assert flag_to_mask(flag, [31])[0, 0]
    assert not flag_to_mask(np.array([[0]], dtype=np.int32), [31])[0, 0]


def test_badpix_mask_catches_unflagged_nan_and_zero_variance():
    """The sample data has ~16 % NaN science pixels carrying no flag bit."""
    sci = np.ones((4, 4)); sci[0, 0] = np.nan
    var = np.ones((4, 4)); var[1, 1] = 0.0; var[2, 2] = -1.0
    flag = np.zeros((4, 4), dtype=np.int32)
    mask, rep = build_badpix_mask(sci, var, flag, bad_bits=[2])
    assert mask[0, 0] and mask[1, 1] and mask[2, 2]
    assert rep["n_badpix"] == 3 and rep["n_flag"] == 0


# ------------------------------------------------------------------ apertures
def test_pixel_scale_km_matches_small_angle_formula():
    kmpp = pixel_scale_km(6.2, 1.0)
    assert abs(kmpp - 4496.1) < 1.0                  # 6.2" at 1 au
    assert np.isnan(pixel_scale_km(6.2, -1.0))


def test_aperture_validity_window():
    cfg = Config(aperture_radii_km=[10_000.0, 100_000.0],
                 aperture_radii_pix=[2.0], r_in_pix=15.0, annulus_r_in_km=None)
    # r_obs = 1 au -> 4496 km/pix: 10^4 km = 2.2 px (ok), 10^5 km = 22 px (too big)
    specs = {s.label: s for s in build_apertures(
        cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)}
    assert specs["km10000"].valid
    assert not specs["km100000"].valid
    assert specs["km100000"].reject_reason == "beyond_r_in"
    assert specs["pix02.0"].valid

    # A distant target loses its small apertures to the PSF instead.
    far = {s.label: s for s in build_apertures(
        cfg, pix_scale_arcsec=6.2, r_obs_au=4.0, psf_fwhm_pix=1.0)}
    assert not far["km10000"].valid
    assert far["km10000"].reject_reason == "below_psf_fwhm"


# ----------------------------------------------------------------- sourceflag
def test_effective_gmag_sums_in_flux():
    assert abs(effective_gmag(np.array([15.0, 15.0])) - 14.247) < 1e-3
    assert np.isnan(effective_gmag(np.array([])))


def _field(gmag, dist):
    return GaiaPixelField(np.array([0.0]), np.array([0.0]),
                          np.array([float(gmag)]), np.array([float(dist)]))


def test_sourceflag_priority_is_a_over_b_over_c_over_d():
    cfg = Config()
    # 'a': G = 12 < 13 inside r_ap + 2*FWHM = 4 px
    assert evaluate_sourceflag(_field(12.0, 3.0), 2.0, 1.0,
                               vmag=18.0, snr=50, cfg=cfg).sourceflag == "a"
    # 'b': too faint for 'a', but contributes > 20 % of a faint comet's flux
    assert evaluate_sourceflag(_field(19.0, 2.0), 2.0, 1.0,
                               vmag=21.0, snr=50, cfg=cfg).sourceflag == "b"
    # 'c': present, but negligible against a much brighter comet
    assert evaluate_sourceflag(_field(21.0, 2.0), 2.0, 1.0,
                               vmag=18.0, snr=50, cfg=cfg).sourceflag == "c"
    empty = GaiaPixelField(np.empty(0), np.empty(0), np.empty(0), np.empty(0))
    assert evaluate_sourceflag(empty, 2.0, 1.0, vmag=20.0, snr=0.5, cfg=cfg).sourceflag == "d"
    assert evaluate_sourceflag(empty, 2.0, 1.0, vmag=20.0, snr=50.0, cfg=cfg).sourceflag == "0"


def test_flag_a_uses_a_two_fwhm_radius_and_a_g13_threshold():
    cfg = Config()
    # Inside r_ap + 2*FWHM = 4 px and brighter than G = 13 -> 'a'
    assert evaluate_sourceflag(_field(12.9, 3.9), 2.0, 1.0,
                               vmag=18.0, snr=50, cfg=cfg).sourceflag == "a"
    # Same star just outside the 2*FWHM radius -> not 'a'
    assert evaluate_sourceflag(_field(12.9, 4.5), 2.0, 1.0,
                               vmag=18.0, snr=50, cfg=cfg).sourceflag != "a"
    # Inside the radius but fainter than the threshold -> not 'a'
    assert evaluate_sourceflag(_field(13.5, 3.0), 2.0, 1.0,
                               vmag=8.0, snr=50, cfg=cfg).sourceflag != "a"


def test_flag_b_is_a_flux_ratio_not_a_magnitude_ratio():
    """G_eff must be brighter than vmag - 2.5*log10(0.2) = vmag + 1.747."""
    cfg = Config()
    vmag = 20.0
    edge = vmag - 2.5 * np.log10(cfg.sourceflag_b_flux_frac)
    assert evaluate_sourceflag(_field(edge - 0.2, 2.0), 2.0, 1.0,
                               vmag=vmag, snr=50, cfg=cfg).sourceflag == "b"
    assert evaluate_sourceflag(_field(edge + 0.2, 2.0), 2.0, 1.0,
                               vmag=vmag, snr=50, cfg=cfg).sourceflag == "c"
    # The flux relation itself
    assert abs(mag_to_flux(edge) / mag_to_flux(vmag) - cfg.sourceflag_b_flux_frac) < 1e-9


def test_sourceflag_only_counts_sources_inside_the_test_radius():
    cfg = Config()
    far = GaiaPixelField(np.array([0.0]), np.array([0.0]),
                         np.array([12.0]), np.array([40.0]))
    res = evaluate_sourceflag(far, 2.0, 1.0, vmag=18.0, snr=50, cfg=cfg)
    assert res.n_gaia == 0 and res.sourceflag == "0"
    assert res.dist_gmag_nearest == 40.0          # still reported


# ---------------------------------------------------------------- photometry
def _synthetic(sky=10.0, flux=100.0, sigma=1.2, size=91, noise=0.0, seed=0):
    """
    A Gaussian source of known total flux on a sky of known level.

    ``noise`` adds Gaussian pixel noise; the variance plane is set to match it,
    which is what the error-budget tests need -- on a perfectly flat sky the
    annulus scatter is zero and every ``sky_noise_mode`` collapses to the same
    number.
    """
    y, x = np.mgrid[:size, :size]
    c = size // 2
    src = flux * np.exp(-((x - c) ** 2 + (y - c) ** 2) / (2 * sigma ** 2)) / (2 * np.pi * sigma ** 2)
    img = src + sky
    if noise:
        img = img + np.random.default_rng(seed).normal(0.0, noise, img.shape)
    var = np.full((size, size), max(noise, 0.1) ** 2)
    return img, var, float(c)


def test_photometry_recovers_a_known_flux_and_subtracts_the_sky():
    cfg = Config(aperture_radii_km=[], aperture_radii_pix=[6.0])
    sci, var, c = _synthetic(sky=10.0, flux=100.0, sigma=1.2)
    bad = np.zeros_like(sci, dtype=bool)
    specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
    df = measure_exposure(sci, var, bad, c, c, specs, cfg)
    assert len(df) == 1
    row = df.iloc[0]
    # A 6-pixel aperture on a sigma = 1.2 px Gaussian encloses essentially all of it.
    assert abs(row["source_sum_mjy"] - 100.0) < 1.0
    assert abs(row["sky_median_mjy_per_pix"] - 10.0) < 0.05
    assert row["n_badpix_ap"] == 0 and not row["badphot"]


def test_bad_pixels_shrink_the_effective_area_and_raise_badphot():
    cfg = Config(aperture_radii_km=[], aperture_radii_pix=[6.0])
    sci, var, c = _synthetic()
    bad = np.zeros_like(sci, dtype=bool)
    ci = int(c)
    bad[ci + 3, ci + 3] = True                    # one bad pixel inside the aperture
    specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
    row = measure_exposure(sci, var, bad, c, c, specs, cfg).iloc[0]
    assert row["n_badpix_ap"] == 1
    assert row["badphot"]
    assert abs(row["aperture_area_pix2"] - row["aperture_area_eff_pix2"] - 1.0) < 1e-6


def test_bad_pixels_in_the_annulus_alone_do_not_raise_badphot():
    cfg = Config(aperture_radii_km=[], aperture_radii_pix=[3.0])
    sci, var, c = _synthetic()
    bad = np.zeros_like(sci, dtype=bool)
    bad[int(c) + 17, int(c)] = True               # inside 15-20 px, outside r_ap = 3
    specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
    row = measure_exposure(sci, var, bad, c, c, specs, cfg).iloc[0]
    assert row["n_badpix_ap"] == 0
    assert not row["badphot"]
    assert row["n_badpix_sky"] >= 1


def test_negative_flux_is_kept_and_magnitude_is_undefined():
    cfg = Config(aperture_radii_km=[], aperture_radii_pix=[4.0])
    sci, var, c = _synthetic(sky=10.0, flux=-50.0)
    bad = np.zeros_like(sci, dtype=bool)
    specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
    row = measure_exposure(sci, var, bad, c, c, specs, cfg).iloc[0]
    assert row["source_sum_mjy"] < 0              # kept, not dropped (S5)
    assert np.isnan(row["abmag"])
    assert np.isfinite(row["source_sum_err_mjy"])


def test_sky_noise_mode_level_is_smaller_than_the_double_counted_form():
    sci, var, c = _synthetic(noise=0.3)
    bad = np.zeros_like(sci, dtype=bool)
    kw = dict(aperture_radii_km=[], aperture_radii_pix=[5.0])
    out = {}
    for mode in ("level", "empirical", "daophot"):
        cfg = Config(sky_noise_mode=mode, **kw)
        specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
        out[mode] = measure_exposure(sci, var, bad, c, c, specs, cfg).iloc[0]["source_sum_err_mjy"]
    assert out["level"] < out["daophot"]          # the S6 double count is gone
    assert out["empirical"] < out["daophot"]


def test_sky_estimate_reports_not_ok_when_the_annulus_is_unusable():
    cfg = Config(min_sky_area=1e6)                # impossible to satisfy
    sci, var, c = _synthetic()
    sky = estimate_sky(sci, var, np.zeros_like(sci, bool), c, c, cfg)
    assert not sky.ok and np.isnan(sky.median)


def test_physical_annulus_rule_is_bounded_by_the_cutout():
    from spherex_apphot.apertures import annulus_radii_pix
    cfg = Config()                                     # 150 000 km, 15-40 px, 5 px wide
    r_in, r_out = annulus_radii_pix(cfg, 4496.1)
    assert abs(r_in - 33.36) < 0.01 and abs(r_out - 38.36) < 0.01
    assert annulus_radii_pix(cfg, 20_000.0) == (15.0, 20.0)          # distant: the near ring
    assert annulus_radii_pix(cfg, 2_763.0) == (40.0, 45.0)           # 24P at 0.4 au: capped
    assert annulus_radii_pix(cfg, float("nan")) == (15.0, 20.0)
    assert annulus_radii_pix(Config(annulus_r_in_km=None), 2_763.0) == (15.0, 20.0)
    # the aperture validity follows the exposure's own ring: 10^5 km = 22 px is now inside 33 px
    specs = {s.label: s for s in build_apertures(
        Config(aperture_radii_km=[100_000.0], aperture_radii_pix=[]),
        pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)}
    assert specs["km100000"].valid
    # and the records report the ring that was used
    cfg = Config(aperture_radii_km=[], aperture_radii_pix=[3.0])
    sci, var, c = _synthetic()
    bad = np.zeros_like(sci, dtype=bool)
    specs = build_apertures(cfg, pix_scale_arcsec=6.2, r_obs_au=1.0, psf_fwhm_pix=1.0)
    row = measure_exposure(sci, var, bad, c, c, specs, cfg, pixel_scale_km=4496.1).iloc[0]
    assert abs(row["r_in_pix"] - 33.36) < 0.01 and abs(row["r_out_pix"] - 38.36) < 0.01
    try:
        Config(annulus_r_in_pix_max=50.0, annulus_r_out_pix_max=45.0).validate()
        raise AssertionError("an inner cap beyond the outer cap must be rejected")
    except ValueError:
        pass


def test_jd_columns_are_written_at_full_precision():
    import pandas as pd
    from spherex_apphot.pipeline import format_jd_columns
    df = pd.DataFrame({"jd_utc": [2460797.659195463, np.nan], "jd_tdb": [2460797.66, 1.0], "wl": [1.23456789012, 2.0]})
    out = format_jd_columns(df, Config())
    assert out.jd_utc.iloc[0] == "2460797.659195463" and out.jd_utc.iloc[1] == ""
    assert out.wl.dtype == float                       # other columns untouched


# -------------------------------------------------------------------- epochs
def test_group_epochs_splits_on_the_configured_gap():
    import pandas as pd
    df = pd.DataFrame({"DATE-OBS": ["2025-01-01", "2025-01-02", "2025-06-01"],
                       "wl": [1.0, 2.0, 3.0]})
    out = group_epochs(df, gap_days=28.0)
    assert list(out["epoch"]) == [1, 1, 2]
    assert list(out["epoch_n"]) == [2, 2, 1]


def test_group_epochs_sorts_chronologically():
    import pandas as pd
    df = pd.DataFrame({"DATE-OBS": ["2025-06-01", "2025-01-01"], "wl": [3.0, 1.0]})
    out = group_epochs(df, gap_days=28.0)
    assert list(out["wl"]) == [1.0, 3.0]


# --------------------------------------------------------------------- slugs
def test_slugify_is_the_single_definition_shared_by_writer_and_resume_check():
    """The prototype wrote '2021G2.csv' but looked for '2021 G2.csv' (R1)."""
    assert slugify("2021 G2") == "2021G2"
    assert slugify("24P") == "24P"
    assert slugify(" 2014 UN271 ") == "2014UN271"


# ------------------------------------------------------------------ stacking
def test_shift_subpixel_conserves_flux_and_does_not_spread_nan():
    img = np.zeros((21, 21)); img[10, 10] = 1.0; img[0, 0] = np.nan
    out = shift_subpixel(img, 0.5, 0.5)
    assert abs(np.nansum(out) - 1.0) < 1e-9
    assert np.isfinite(out[10, 10])               # NaN stayed in its corner


def test_extract_stamp_pads_past_the_edge_with_nan():
    img = np.ones((21, 21))
    stamp = extract_stamp(img, 1.0, 1.0, radius=5, subpixel=False)
    assert stamp.shape == (11, 11)
    assert np.isnan(stamp[0, 0]) and stamp[6, 6] == 1.0


def test_stack_removes_per_frame_sky_and_counts_contributions():
    rng = np.random.default_rng(1)
    items = [StackInput(sci=rng.normal(0, 0.01, (41, 41)) + 100.0 * i,
                        badpix=np.zeros((41, 41), bool),
                        xcen=20.0, ycen=20.0, sky_median=100.0 * i, label=f"f{i}")
             for i in range(5)]
    res = stack_frames(items, radius=6, subtract_sky=True, subpixel=False)
    assert res.n_frames == 5
    assert abs(float(np.nanmean(res.image))) < 0.05      # wildly different skies removed
    assert int(res.count.max()) == 5


# --------------------------------------------------------------- reflectance
def test_distance_corrected_flux_normalises_to_one_au():
    """flux_distcorr is the flux the comet would show at r_hel = r_obs = 1 au."""
    import pandas as pd
    assert Config().distcorr_mode == "multiply"
    df = pd.DataFrame({"source_sum_mjy": [8.0], "source_sum_err_mjy": [0.8],
                       "r_hel": [2.0], "r_obs": [3.0]})
    mul = add_distance_corrected_flux(df, Config())
    assert abs(mul["flux_distcorr_mjy"].iloc[0] - 8.0 * 4.0 * 9.0) < 1e-12
    assert abs(mul["flux_distcorr_err_mjy"].iloc[0] - 0.8 * 4.0 * 9.0) < 1e-12
    # At unit distances the correction is the identity.
    unit = add_distance_corrected_flux(
        pd.DataFrame({"source_sum_mjy": [8.0], "source_sum_err_mjy": [0.8],
                      "r_hel": [1.0], "r_obs": [1.0]}), Config())
    assert abs(unit["flux_distcorr_mjy"].iloc[0] - 8.0) < 1e-12
    div = add_distance_corrected_flux(df, Config(distcorr_mode="divide"))
    assert abs(div["flux_distcorr_mjy"].iloc[0] - 8.0 / 36.0) < 1e-12


def test_distance_correction_removes_the_geometric_trend():
    """
    The same intrinsic coma seen at two distances must give the same corrected
    flux -- that is the whole point of the correction.
    """
    import pandas as pd
    intrinsic = 5.0
    geom = pd.DataFrame({"r_hel": [1.0, 2.5, 4.0], "r_obs": [1.2, 2.0, 3.5]})
    geom["source_sum_mjy"] = intrinsic / (geom.r_hel ** 2 * geom.r_obs ** 2)
    geom["source_sum_err_mjy"] = 0.01
    out = add_distance_corrected_flux(geom, Config())
    assert np.allclose(out["flux_distcorr_mjy"], intrinsic, rtol=1e-12)
    # The uncorrected fluxes span more than an order of magnitude.
    assert out["source_sum_mjy"].max() / out["source_sum_mjy"].min() > 10


def test_stack_median_rejects_a_minority_of_contaminated_frames():
    """Sigma clipping alone does not remove 2 bright frames out of 10."""
    rng = np.random.default_rng(3)
    items = []
    for i in range(10):
        img = rng.normal(0.0, 1.0, (41, 41)) + 5.0
        if i >= 8:
            img[24, 24] += 500.0
        items.append(StackInput(sci=img, badpix=np.zeros((41, 41), bool),
                                xcen=20.0, ycen=20.0, sky_median=5.0, label=f"f{i}"))
    med = stack_frames(items, radius=8, subpixel=False, combine="median")
    mean = stack_frames(items, radius=8, subpixel=False, combine="mean")
    assert abs(med.image[12, 12]) < 5.0
    assert mean.image[12, 12] > 50.0
    assert med.combine == "median"


def _growth_frame(fname, star=False, flag="0"):
    import pandas as pd
    r = np.arange(1.0, 13.0)
    F = 1.0 * r + (5.0 * (r >= 6) if star else 0.0)   # 1/rho coma, plus a star
    return pd.DataFrame(dict(filename=fname, ap_kind="km", r_ap_pix=r,
                             source_sum_mjy=F, source_sum_err_mjy=0.02 * np.ones_like(r),
                             aperture_area_eff_pix2=np.pi * r ** 2, badphot=False,
                             sourceflag=flag, snr=F / 0.02, epoch=1, wl=3.0,
                             gmag_brightest=np.nan, dist_gmag_nearest=np.nan, n_gaia=0))


def test_growth_metrics_separates_a_star_from_a_smooth_coma():
    import pandas as pd
    df = pd.concat([_growth_frame(f"clean{i}") for i in range(8)]
                   + [_growth_frame(f"star{i}", star=True, flag="a") for i in range(2)],
                   ignore_index=True)
    m = growth_metrics(df)
    assert len(m) == 10
    clean = m[m.flag_worst == "0"].sb_rise_sigma
    dirty = m[m.flag_worst == "a"].sb_rise_sigma
    assert clean.max() < 3.0 < dirty.min()      # a smooth coma never rises


def test_flag_effectiveness_scores_a_perfect_flag_perfectly():
    import pandas as pd
    df = pd.concat([_growth_frame(f"clean{i}") for i in range(8)]
                   + [_growth_frame(f"star{i}", star=True, flag="a") for i in range(2)],
                   ignore_index=True)
    eff = flag_effectiveness(growth_metrics(df), threshold=3.0)
    cut_a = eff[eff.cut == "cut a"].iloc[0]
    assert cut_a["recall"] == 1.0 and cut_a["precision"] == 1.0
    assert cut_a["purity_kept"] == 1.0


def test_status_log_replaces_rows_and_emits_no_warnings():
    """A re-run must update a target's row, not append a second one."""
    import tempfile, warnings
    from spherex_apphot.status import StatusLog
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # any warning fails the test
        with tempfile.TemporaryDirectory() as d:
            log = StatusLog(Path(d) / "status.csv")
            log.update("2P", "running")
            log.update("2021 G2", "ok", n_rows_out=42)
            log.update("2P", "ok", n_rows_out=7)
            df = log.read()
            assert len(df) == 2
            assert int((df.slug == "2P").sum()) == 1
            assert df.loc[df.slug == "2P", "status"].iloc[0] == "ok"
            assert log.is_done("2021 G2") and not log.is_done("10P")


def test_target_list_drops_the_summary_row():
    """The v2607 sheet ends with a Total row whose desig cell holds the count."""
    import tempfile
    import pandas as pd
    from spherex_apphot.targetlist import read_target_list
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / "list.csv"
        pd.DataFrame({"Unnamed: 0": [None, None, "Total"],
                      "desig": ["2P", "2021 G2", 2]}).to_csv(f, index=False)
        assert read_target_list(f) == ["2P", "2021 G2"]


def test_worst_flag_respects_priority():
    assert worst_flag(["0", "c", "a", "d"]) == "a"
    assert worst_flag(["0", "d"]) == "d"
    assert worst_flag(["0", "0"]) == "0"


def test_reflectance_applies_the_distance_scaling():
    import pandas as pd
    cfg = Config(refl_apply_distance=True)
    df = pd.DataFrame({"source_sum_mjy": [1.0, 1.0], "source_sum_err_mjy": [0.1, 0.1],
                       "sun_jy": [1e14, 1e14], "r_hel": [1.0, 2.0], "r_obs": [1.0, 1.0]})
    out = add_reflectance(df, cfg)
    assert abs(out["refl"].iloc[1] / out["refl"].iloc[0] - 4.0) < 1e-9   # r_hel^2
    # Errors scale with the same factor -- the relation is linear in flux.
    assert abs(out["refl_err"].iloc[1] / out["refl_err"].iloc[0] - 4.0) < 1e-9


def test_bin_spectrum_weights_by_inverse_variance():
    import pandas as pd
    df = pd.DataFrame({"wl": [1.0, 1.01, 1.02],
                       "refl_norm": [1.0, 1.0, 10.0],
                       "refl_norm_err": [0.01, 0.01, 10.0]})
    out = bin_spectrum(df, wl_bins=np.array([0.9, 1.1]))
    assert len(out) == 1
    assert abs(out["value"].iloc[0] - 1.0) < 0.05   # the noisy outlier carries no weight


if __name__ == "__main__":
    import traceback
    tests = [(n, f) for n, f in sorted(globals().items())
             if n.startswith("test_") and callable(f)]
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
