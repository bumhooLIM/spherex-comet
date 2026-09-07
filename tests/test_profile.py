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
