"""Tests for the radial-profile diagnostic, on synthetic images.

The point of the module is to tell a 1/rho coma from a stellar PSF, so the
tests build both and check the fitted slopes come out where they should.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ztfcomet import profile as pf


def _grid(n=61):
    y, x = np.mgrid[0:n, 0:n]
    return x.astype(float), y.astype(float), (n - 1) / 2.0


def gaussian(amp=1000.0, fwhm=2.5, n=61, x0=None, y0=None):
    x, y, c = _grid(n)
    x0 = c if x0 is None else x0
    y0 = c if y0 is None else y0
    s = fwhm / 2.3548
    return amp * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * s * s))


def coma(amp=1000.0, n=61):
    """Surface brightness proportional to 1/rho, softened inside 0.5 px."""
    x, y, c = _grid(n)
    r = np.hypot(x - c, y - c)
    return amp / np.maximum(r, 0.5)


def test_star_profile_is_steep_and_coma_is_minus_one():
    n = 61
    c = (n - 1) / 2.0
    star = pf.radial_profile(gaussian(fwhm=2.5), c, c, oversample=4)
    com = pf.radial_profile(coma(), c, c, oversample=4)

    m_star, _, _ = pf.fit_powerlaw(star["r_pix"], star["mean_clip"], 3.75, 10.0)
    m_coma, _, _ = pf.fit_powerlaw(com["r_pix"], com["mean_clip"], 3.75, 10.0)

    assert m_coma == pytest.approx(-1.0, abs=0.05), f"coma slope {m_coma}"
    assert m_star < -3.0, f"star slope {m_star} not steep"


def test_oversampling_preserves_surface_brightness():
    """Bilinear zoom must not change the annulus mean of a smooth profile.

    Checked on the 1/rho coma, where SB varies gently across a 0.5 px annulus.
    It is deliberately NOT checked on the Gaussian: native pixels sample only
    the discrete radii that happen to exist (eight pixels at r = 3.61 for the
    r = 3.5 annulus), so for a steep PSF the native "mean" is the value at one
    radius, not the area-weighted mean over the annulus -- the oversampled
    value differs from it by design.  That is the bias oversampling removes.
    """
    n = 61
    c = (n - 1) / 2.0
    img = coma()
    native = pf.radial_profile(img, c, c, oversample=1)
    over = pf.radial_profile(img, c, c, oversample=4)

    both = (native["r_pix"] >= 3.5) & np.isfinite(native["mean_clip"]) & np.isfinite(over["mean_clip"])
    assert both.sum() >= 8
    ratio = (over.loc[both, "mean_clip"] / native.loc[both, "mean_clip"]).to_numpy()
    assert np.allclose(ratio, 1.0, atol=0.05), ratio

    # Native sampling leaves some annuli with too few pixels to define a mean;
    # oversampling fills every one of them.
    assert native["mean_clip"].isna().any()
    assert over["mean_clip"].notna().all()
    assert (over["n"] > native["n"]).all()


def test_zoom_is_flux_conserving():
    """The oversampled cutout must hold the same total flux per original pixel."""
    from scipy import ndimage
    img = gaussian(amp=1000.0, fwhm=2.5, n=31)
    z = ndimage.zoom(img, 4, order=1, grid_mode=True, mode="nearest")
    assert z.sum() / 16 == pytest.approx(img.sum(), rel=1e-6)


def test_sigma_clipping_rejects_an_intruding_star():
    """A field star inside one annulus must be clipped out of the mean."""
    n = 61
    c = (n - 1) / 2.0
    img = coma() + gaussian(amp=5000.0, fwhm=2.0, x0=c + 5.0, y0=c)   # star at r = 5
    prof = pf.radial_profile(img, c, c, oversample=4)
    clean = pf.radial_profile(coma(), c, c, oversample=4)
    at5 = prof["r_pix"] == 5.0

    plain = float(prof.loc[at5, "mean"].iloc[0])
    clipped = float(prof.loc[at5, "mean_clip"].iloc[0])
    truth = float(clean.loc[at5, "mean_clip"].iloc[0])

    assert plain > 1.5 * truth, "test set-up: the star did not perturb the plain mean"
    assert abs(clipped - truth) / truth < 0.15, f"clipped {clipped} vs truth {truth}"
    assert clipped < plain


def test_fit_powerlaw_needs_three_points():
    slope, _, n = pf.fit_powerlaw([1, 2], [1, 0.5], 0.5, 10)
    assert np.isnan(slope) and n == 2


def test_central_sb_is_pixel_weighted():
    prof = pd.DataFrame({"r_pix": [0.5, 1.0, 1.5, 2.0], "mean_clip": [10.0, 8.0, 6.0, 4.0],
                         "n": [4, 12, 20, 28]})
    c0 = pf.central_sb(prof, r_core=1.5)
    assert c0 == pytest.approx((10 * 4 + 8 * 12 + 6 * 20) / 36)


def _field(n=400, seed=1, saturate=60000.0):
    """Synthetic sky with stars of a range of brightness, plus special cases."""
    rng = np.random.default_rng(seed)
    x, y, _ = _grid(n)
    img = rng.normal(200.0, 3.0, (n, n))
    stars = []
    for _ in range(40):
        sx, sy = rng.uniform(30, n - 30, 2)
        amp = 10 ** rng.uniform(2.3, 3.8)
        img += gaussian(amp=amp, fwhm=2.5, n=n, x0=sx, y0=sy)
        stars.append((sx, sy, amp))
    img += gaussian(amp=saturate * 2, fwhm=2.5, n=n, x0=200, y0=300)    # saturated
    img += gaussian(amp=3000, fwhm=2.5, n=n, x0=5, y0=200)              # at the edge
    img += gaussian(amp=3000, fwhm=2.5, n=n, x0=120, y0=120)            # "comet" position
    return img, stars


def test_field_star_selection_respects_the_cuts():
    img, _ = _field()
    stars, _ = pf.select_field_stars(img, gain=6.0, readnoise=9.0, saturate=60000.0,
                                     exclude_xy=(120.0, 120.0))
    assert 0 < len(stars) <= 20
    assert (stars["snr"] > 10).all()
    assert (stars["peak"] < 0.8 * 60000).all(), "saturated star kept"
    assert not ((abs(stars["x"] - 200) < 3) & (abs(stars["y"] - 300) < 3)).any()
    assert (stars["x"] > 12).all(), "edge star kept"
    assert (np.hypot(stars["x"] - 120, stars["y"] - 120) > 30).all(), "star next to the comet kept"


def test_star_stack_is_normalised_and_steep():
    img, _ = _field()
    stars, bkg = pf.select_field_stars(img, 6.0, 9.0, 60000.0, exclude_xy=(120.0, 120.0))
    sub = np.ascontiguousarray(img, np.float32) - bkg
    stack, per_star = pf.stack_star_profiles(sub, stars)

    assert stack["n_stars"].max() == len(per_star["star"].unique())
    core = stack[stack["r_pix"] <= 1.5]
    assert core["median"].mean() == pytest.approx(1.0, abs=0.35)
    slope, _, _ = pf.fit_powerlaw(stack["r_pix"], stack["median"], 3.75, 10.0)
    assert slope < -2.5


# --------------------------------------------------------- PSF-convolved model
def _synthetic_scene(m_true=1.0, nuc_flux=0.0, coma_amp=300.0, fwhm=2.8, n=81, seed=0):
    """A coma rho^-m plus optional nucleus, convolved with a Gaussian PSF, plus noise."""
    from scipy.signal import fftconvolve
    x, y, c = _grid(n)
    r = np.hypot(x - c, y - c)
    intrinsic = coma_amp * np.maximum(r, 0.25) ** (-m_true)
    if nuc_flux > 0:
        intrinsic = intrinsic.copy()
        intrinsic[int(c), int(c)] += nuc_flux
    psf = gaussian(amp=1.0, fwhm=fwhm, n=n)
    psf /= psf.sum()
    img = fftconvolve(intrinsic, psf, mode="same")
    rng = np.random.default_rng(seed)
    img = img + rng.normal(0, 0.5, img.shape) + 100.0        # sky 100, noise 0.5
    return img, psf, c


def _star_stack_from_psf(psf, c, n_stars=12):
    """Pretend the PSF image is a field star and build a 'stack' from it."""
    prof = pf.radial_profile(psf * 1e4, c, c, oversample=4)
    c0 = pf.central_sb(prof)
    return pd.DataFrame({"r_pix": prof["r_pix"], "median": prof["mean_clip"] / c0,
                         "std": 0.02 * prof["mean_clip"] / c0, "n_stars": n_stars})


# Coma flux inside 10 px is 2*pi*300*10 ~ 19 000 DN, so 6000 DN of nucleus is
# ~24% -- comfortably above the 15% threshold, not on top of it.
@pytest.mark.parametrize("m_true,nuc", [(1.0, 0.0), (1.0, 6000.0), (1.5, 0.0), (0.7, 0.0)])
def test_coma_model_recovers_intrinsic_slope(m_true, nuc):
    """A power law fitted outside the core is biased by the PSF and nucleus;
    the forward model must recover the intrinsic slope regardless."""
    img, psf, c = _synthetic_scene(m_true=m_true, nuc_flux=nuc, fwhm=2.8)
    comet = pf.radial_profile(img, c, c, oversample=4, sky=100.0)
    stack = _star_stack_from_psf(psf, c)

    fit = pf.fit_coma_model(comet, stack, oversample=4)
    assert fit["success"]
    assert fit["m"] == pytest.approx(m_true, abs=0.12), f"m={fit['m']}"
    if nuc > 0:
        assert fit["nucleus_fraction"] > 0.15, fit
    else:
        assert fit["nucleus_fraction"] < 0.10, fit


def test_naive_slope_is_biased_when_the_core_is_large_but_model_is_not():
    """A strong nucleus under a wide PSF biases the naive slope steep; the
    forward model removes it.

    Calibrated numerically: with the nucleus carrying ~70% of the flux inside
    10 px and FWHM 4.5 px, the naive fit from 1.5 FWHM outward gives ~-1.23 and
    the model ~-0.92.  Note the size of the effect -- even this extreme case
    cannot push the naive slope below about -1.25, which is itself a finding:
    nucleus + PSF leakage is not how a profile gets to -2.
    """
    img, psf, c = _synthetic_scene(m_true=1.0, nuc_flux=40000.0, fwhm=4.5)
    comet = pf.radial_profile(img, c, c, oversample=4, sky=100.0)
    stack = _star_stack_from_psf(psf, c)

    naive, _, _ = pf.fit_powerlaw(comet["r_pix"], comet["mean_clip"], 1.5 * 4.5, 10.0)
    fit = pf.fit_coma_model(comet, stack, oversample=4)
    assert naive < -1.15, f"expected the naive slope to be biased steep, got {naive}"
    assert fit["m"] == pytest.approx(1.0, abs=0.15), f"model m={fit['m']}"
    assert -fit["m"] > naive + 0.15, "model did not remove the bias"
    assert fit["nucleus_fraction"] > 0.5


def test_fixed_slope_fit_is_worse_when_slope_is_wrong():
    img, psf, c = _synthetic_scene(m_true=1.5, nuc_flux=0.0)
    comet = pf.radial_profile(img, c, c, oversample=4, sky=100.0)
    stack = _star_stack_from_psf(psf, c)
    free = pf.fit_coma_model(comet, stack, oversample=4)
    fixed = pf.fit_coma_model(comet, stack, oversample=4, m_fixed=1.0)
    assert fixed["chi2_red"] > 2 * free["chi2_red"]



def test_sky_over_subtraction_steepens_naive_slope_and_fit_sky_recovers_it():
    """The mechanism behind steep faint-comet profiles: the outer annuli sit
    within a few sigma of the sky, so a small sky error tilts the whole tail.
    Naive fit: biased.  Model with a free sky term: recovers m and the error."""
    img, psf, c = _synthetic_scene(m_true=1.0, nuc_flux=0.0, coma_amp=40.0, fwhm=2.8)
    # True sky is 100; subtract 102 (over-subtraction of 2 DN, ~ the outer SB).
    comet = pf.radial_profile(img, c, c, oversample=4, sky=102.0)
    stack = _star_stack_from_psf(psf, c)

    naive, _, _ = pf.fit_powerlaw(comet["r_pix"], comet["mean_clip"], 1.5 * 2.8, 10.0)
    plain = pf.fit_coma_model(comet, stack, oversample=4)
    with_sky = pf.fit_coma_model(comet, stack, oversample=4, fit_sky=True)

    assert naive < -1.3, f"sky error did not steepen the naive slope: {naive}"
    assert with_sky["m"] == pytest.approx(1.0, abs=0.15), with_sky
    assert with_sky["sky_offset"] == pytest.approx(-2.0, abs=0.7), with_sky
    assert abs(with_sky["m"] - 1.0) < abs(plain["m"] - 1.0)


# ---------------------------------------------------------------- centring
def test_refine_centre_recovers_a_sub_pixel_peak_from_an_offset_start():
    n = 61
    x0, y0 = 30.3, 29.8
    img = gaussian(amp=2000.0, fwhm=2.5, n=n, x0=x0, y0=y0) + 100.0
    img += np.random.default_rng(1).normal(0, 0.5, img.shape)
    cen = pf.refine_centre(img, x0 + 1.5, y0 - 1.0)
    assert cen["refined"]
    assert abs(cen["x"] - x0) < 0.15 and abs(cen["y"] - y0) < 0.15
    assert 1.5 < cen["shift_pix"] < 2.0


def test_refine_centre_finds_the_nucleus_of_a_coma_started_off_peak():
    img, _, c = _synthetic_scene(m_true=1.0, nuc_flux=6000.0)
    cen = pf.refine_centre(img, c + 1.5, c + 0.5)
    assert cen["refined"]
    assert np.hypot(cen["x"] - c, cen["y"] - c) < 0.25


def test_refine_centre_keeps_the_input_on_pure_noise():
    rng = np.random.default_rng(2)
    img = 100.0 + rng.normal(0, 1.0, (61, 61))
    cen = pf.refine_centre(img, 30.0, 30.0)
    assert not cen["refined"]
    assert cen["shift_pix"] == 0.0 and cen["x"] == 30.0


def test_refine_centre_leaves_an_on_peak_comet_alone():
    img, _, c = _synthetic_scene(m_true=1.0, nuc_flux=6000.0)
    cen = pf.refine_centre(img, c, c)
    assert cen["shift_pix"] < 0.2


def test_refine_centre_rarely_moves_on_pure_noise():
    # The contrast statistic is a max over patches minus one sample, so it is
    # biased high on noise; the guard must hold well below the rate that a
    # ring-based noise estimate allowed (17.5% at 3 sigma).
    rng = np.random.default_rng(3)
    moved = 0
    for _ in range(200):
        img = 100.0 + rng.normal(0, 1.0, (41, 41))
        moved += pf.refine_centre(img, 20.0 + rng.uniform(-.4, .4), 20.0 + rng.uniform(-.4, .4))["refined"]
    assert moved <= 4


def test_refine_centre_reaches_a_far_nucleus_only_with_the_ephemeris():
    # winpos dragged 3 px down a tail; the ephemeris sits 0.3 px from the nucleus
    img, _, c = _synthetic_scene(m_true=1.0, nuc_flux=6000.0)
    alone = pf.refine_centre(img, c + 3.0, c)
    with_ref = pf.refine_centre(img, c + 3.0, c, x_ref=c - 0.3, y_ref=c + 0.2)
    assert np.hypot(alone["x"] - c, alone["y"] - c) > 0.8          # 2 px disc cannot reach it
    # separation 3.31 px -> half 1.65 + 2.0 = 3.65, rounded up to the half pixel
    assert with_ref["refined"] and with_ref["search_pix"] == 4.0
    assert np.hypot(with_ref["x"] - c, with_ref["y"] - c) < 0.25
