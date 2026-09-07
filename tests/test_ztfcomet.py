"""Offline tests for :mod:`ztfcomet`.

No network: every test here runs against synthetic frames or fixed inputs, so
the suite is fast and deterministic.  The tests named ``test_regression_c*``
pin the specific defects catalogued in ``doc/primitive_code_analysis.md`` — they
exist to make sure those bugs cannot come back.
"""

from __future__ import annotations

import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import ztfcomet as zc
from ztfcomet import config as cfg
from ztfcomet import cutout, directory, phot, query


# --------------------------------------------------------------------- directory
def test_target_slug_strips_whitespace_and_separators():
    assert directory.target_slug("2019 Y3") == "2019Y3"
    assert directory.target_slug("C/2024 E1") == "C2024E1"
    assert directory.target_slug("24P") == "24P"


def test_project_root_is_the_repository(tmp_path):
    assert (directory.PROJECT_ROOT / "pyproject.toml").is_file()
    assert (directory.PROJECT_ROOT / "ztfcomet").is_dir()


def test_data_root_honours_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("ZTFCOMET_DATA", str(tmp_path))
    assert directory._resolve_data_root() == tmp_path


def test_data_dir_creates_per_target_subdirectory(monkeypatch, tmp_path):
    monkeypatch.setattr(directory, "DATA_ROOT", tmp_path)
    path = directory.data_dir("C/2024 E1")
    assert path == tmp_path / "C2024E1" and path.is_dir()


# ------------------------------------------------------------------------ config
def test_get_target_is_case_insensitive():
    assert zc.get_target("24p").horizons_id == zc.get_target("24P").horizons_id


def test_unknown_target_falls_through_to_horizons():
    target = zc.get_target("C/2099 Z9")
    assert target.query_designation == "C/2099 Z9"
    assert target.name == "C/2099Z9"


def test_orbit_record_delegates_to_the_designation_resolver(monkeypatch):
    """With no pinned records, resolution is dynamic and fragment-aware.

    Record numbers are deliberately *not* hardcoded any more: the ones the
    pre-merge notebooks carried for 240P now resolve to different comets.
    """
    from ztfcomet import horizons as hz

    seen = {}

    def fake(designation, epoch_jd=None, allow_fragment=False, location=None):
        seen.update(designation=designation, epoch_jd=epoch_jd,
                    allow_fragment=allow_fragment)
        return 90001212

    monkeypatch.setattr(hz, "resolve_target_id", fake)
    target = zc.get_target("240P")
    assert target.resolve_orbit_record(2460900.0) == 90001212
    assert seen["designation"] == "240P"
    assert seen["allow_fragment"] is False


def test_pinned_orbit_records_still_override(monkeypatch):
    """The escape hatch works, but it is deprecated and must not hit the network."""
    from ztfcomet import horizons as hz

    def boom(*a, **k):
        raise AssertionError("must not resolve when records are pinned")

    monkeypatch.setattr(hz, "resolve_target_id", boom)
    target = cfg.Target(name="X", designation="X",
                        orbit_records={0.0: 111111, 2459580.5: 222222})
    assert target.resolve_orbit_record(2458000.0) == 111111
    assert target.resolve_orbit_record(2460900.0) == 222222


def test_regression_c9_each_target_has_its_own_designation():
    """``ztfquery_2P.ipynb`` shipped 24P's orbit ID inside the 2P notebook."""
    designations = [t.query_designation for t in cfg.TARGETS.values()]
    assert len(designations) == len(set(designations)), "two targets share a designation"
    assert cfg.TARGETS["2P"].query_designation != cfg.TARGETS["24P"].query_designation


def test_no_target_pins_a_horizons_record_number():
    """Record numbers are not stable; designations are."""
    pinned = {name: t.orbit_records for name, t in cfg.TARGETS.items() if t.orbit_records}
    assert not pinned, f"pinned record numbers will go stale: {pinned}"


# ------------------------------------------------------------------------- query
def test_extract_lastrecnum_takes_the_most_recent_record():
    message = (
        "Ambiguous target name; provide unique id:\n"
        "  Record #  Epoch-yr  Primary Desig\n"
        "  90001203    2018     240P\n"
        "  90001204    2025     240P\n"
    )
    assert query.extract_lastrecnum(message) == "90001204"


def test_extract_lastrecnum_returns_none_for_other_errors():
    assert query.extract_lastrecnum("Unable to connect") is None


def _frames_and_eph():
    ztf = pd.DataFrame({
        "obsjd": [2460002.5, 2460000.5, 2460001.5],
        "airmass": [1.1, 1.2, 1.3],                   # also supplied by Horizons
        "field": [1, 2, 3],
    })
    # Horizons returns epochs sorted, i.e. NOT in the caller's order.
    eph = pd.DataFrame({
        "datetime_jd": [2460000.5, 2460001.5, 2460002.5],
        "RA": [10.0, 20.0, 30.0],
        "airmass": [9.1, 9.2, 9.3],
    })
    return ztf, eph


def test_regression_c10_ephemeris_joins_on_jd_not_row_position():
    """The predecessor's ``pd.concat(axis=1)`` matched by position.

    With an unsorted metadata table that silently attaches every ephemeris to
    the wrong frame; the JD join must get it right regardless of order.
    """
    ztf, eph = _frames_and_eph()
    merged = query._merge_on_jd(ztf, eph)

    for _, row in merged.iterrows():
        assert row["obsjd"] == pytest.approx(row["datetime_jd"])
    assert dict(zip(merged["obsjd"], merged["RA"])) == {
        2460000.5: 10.0, 2460001.5: 20.0, 2460002.5: 30.0}


def test_regression_c11_colliding_columns_are_namespaced():
    """Both services supply ``airmass``; a duplicate label breaks scalar access."""
    ztf, eph = _frames_and_eph()
    merged = query._merge_on_jd(ztf, eph)

    assert not merged.columns.duplicated().any()
    assert "eph_airmass" in merged.columns
    assert isinstance(merged["airmass"], pd.Series)


def test_regression_c12_failed_record_retry_returns_empty_not_unboundlocal(monkeypatch):
    """A failing retry used to fall through to ``return eph`` with eph unbound."""
    class Boom:
        def __init__(self, *a, **k): pass
        def ephemerides(self, **k):
            raise ValueError("Ambiguous target name\n90001204  2025  240P\n")

    monkeypatch.setattr(query, "Horizons", Boom)
    result = query.query_sso_ephemeris("240P", epochs=[2460000.5])
    assert isinstance(result, pd.DataFrame) and result.empty


def test_ephemeris_chunking_splits_long_epoch_lists(monkeypatch):
    """astroquery puts epochs in the URL, so long lists must be chunked."""
    seen = []

    class Fake:
        def __init__(self, id=None, location=None, epochs=None):
            seen.append(len(epochs))
            self._n = len(epochs)
        def ephemerides(self, **k):
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self): return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(query, "Horizons", Fake)
    out = query.query_sso_ephemeris("x", epochs=list(range(120)), max_epochs_per_call=50)
    assert seen == [50, 50, 20]
    assert len(out) == 120


# ------------------------------------------------------------------------ cutout
def _metadata_row(**kw):
    row = dict(filefracday="20250116517338", field=724, filtercode="zg", ccdid=9,
               imgtypecode="o", qid=3, RA=263.3, DEC=40.8)
    row.update(kw)
    return pd.Series(row)


def test_construct_fitsurl_matches_the_irsa_layout():
    url = cutout.construct_fitsurl(_metadata_row(), is_cutout=False)
    assert url == (
        "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"
        "2025/0116/517338/ztf_20250116517338_000724_zg_c09_o_q3_sciimg.fits")


def test_construct_fitsurl_appends_the_cutout_query():
    url = cutout.construct_fitsurl(_metadata_row(), is_cutout=True, cutout_size="10arcmin")
    assert url.endswith("?center=263.3,40.8&size=10arcmin&gzip=false")


def test_build_urls_drops_filename_collisions():
    """Two cutout centres on one base image collapse to the same local name."""
    frames = pd.DataFrame([_metadata_row(RA=263.3), _metadata_row(RA=263.9)])
    urls = cutout.build_urls(frames)
    assert len(urls) == 1


@pytest.mark.parametrize("payload,expected", [
    (b"SIMPLE  =                    T" + b" " * 4000, True),
    (b"<html><title>404 Not Found</title></html>", False),   # the real failure mode
    (b"", False),
    (b"SIMPLE", False),                                       # too short to be real
])
def test_regression_c2_fits_validation_rejects_error_pages(tmp_path, payload, expected):
    path = tmp_path / "x.fits"
    path.write_bytes(payload)
    assert cutout._looks_like_fits(path) is expected


def test_regression_c2_corrupt_file_is_repaired_not_skipped(tmp_path, monkeypatch):
    """``if save_path.exists(): continue`` made bad downloads permanent."""
    dest = tmp_path / "ztf_1_sciimg.fits"
    dest.write_bytes(b"<html>404</html>")          # poisoned by an earlier run

    good = b"SIMPLE  =                    T" + b" " * 4000

    class Response:
        headers = {"Content-Type": "application/fits"}
        def raise_for_status(self): pass
        def iter_content(self, chunk_size=None): yield good
        def __enter__(self): return self
        def __exit__(self, *a): return False

    monkeypatch.setattr(cutout.requests, "get", lambda *a, **k: Response())
    report = cutout.download_urls(["http://x/ztf_1_sciimg.fits"], tmp_path,
                                  repair=True, progress=False)

    assert report.n_downloaded == 1
    assert dest.read_bytes() == good
    assert not list(tmp_path.glob("*.part")), "left a partial file behind"


def test_download_rejects_html_served_with_status_200(tmp_path, monkeypatch):
    class Response:
        headers = {"Content-Type": "text/html"}
        def raise_for_status(self): pass
        def iter_content(self, chunk_size=None): yield b"<html>404</html>"
        def __enter__(self): return self
        def __exit__(self, *a): return False

    monkeypatch.setattr(cutout.requests, "get", lambda *a, **k: Response())
    report = cutout.download_urls(["http://x/a_sciimg.fits"], tmp_path,
                                  progress=False, max_retries=1)

    assert report.n_downloaded == 0
    assert not (tmp_path / "a_sciimg.fits").exists()


# -------------------------------------------------------------------------- phot
def test_variance_map_clips_negatives_instead_of_producing_nan():
    """``sqrt(data/gain + ...)`` gave 97.5% NaN on one real frame."""
    data = np.array([[-100.0, 0.0], [100.0, 400.0]])
    variance, negative = phot._variance_map(data, gain=6.2, readnoise=9.7)

    assert np.isfinite(variance).all()
    assert (variance > 0).all()
    assert negative.sum() == 1


def test_apcor_interpolates_between_tabulated_diameters():
    row = {f"apcor{i}": v for i, v in enumerate(
        [-1.6, -0.93, -0.55, -0.21, -0.042, -0.0095], start=1)}

    # radius 1 px  -> diameter 2 px, the first tabulated point
    assert phot._apcor_for_radius(row, 1.0) == pytest.approx(-1.6)
    # radius 7 px  -> diameter 14 px, the last one
    assert phot._apcor_for_radius(row, 7.0) == pytest.approx(-0.0095)
    # beyond the table it is clamped, and small
    assert phot._apcor_for_radius(row, 50.0) == pytest.approx(-0.0095)
    # in between, monotonic and bracketed
    mid = phot._apcor_for_radius(row, 2.5)
    assert -0.55 < mid < -0.21


def test_apcor_is_zero_without_the_keywords():
    assert phot._apcor_for_radius({}, 5.0) == 0.0


def _synthetic_table(n=4):
    """Frames with deliberately *different* zeropoints and filters."""
    return pd.DataFrame({
        "obsjd": 2460000.5 + np.arange(n),
        "file": [f"f{i}.fits" for i in range(n)],
        "filter": ["ZTF_g", "ZTF_r"] * (n // 2),
        "inst_mag": np.full(n, -10.0),
        "inst_mag_err": np.full(n, 0.01),
        "zpmag": 26.0 + np.arange(n) * 0.5,        # 26.0, 26.5, 27.0, 27.5
        "zpmagrms": np.full(n, 0.05),
        "clrcoeff": np.full(n, -0.03),
        "clrcounc": np.full(n, 1e-4),
        "zpclrcov": np.full(n, 0.0),
        "rho_pix": np.full(n, 5.0),
        "rho_km": np.full(n, 15000.0),
        "r": np.full(n, 2.0),
        "delta": np.full(n, 1.5),
        "alpha": np.full(n, 10.0),
        **{c: np.zeros(n, dtype=bool) for c in phot.FLAG_COLUMNS},
    })


def test_sep_winpos_returns_three_values():
    """``sep.winpos`` returns (x, y, flag), not (x, y).

    Unpacking two silently raised, was swallowed by the surrounding
    ``except``, and left every aperture on the unrefined ephemeris position —
    a mis-centring of ~25% of the aperture radius on 79 of 111 real frames,
    visible only as a centroid-shift histogram that was exactly zero everywhere.
    """
    import sep

    # A single Gaussian source offset from the seed position.
    y, x = np.mgrid[0:41, 0:41]
    image = np.exp(-((x - 22.0) ** 2 + (y - 19.0) ** 2) / (2 * 2.0 ** 2)).astype(np.float32)

    result = sep.winpos(np.ascontiguousarray(image), xinit=20.0, yinit=20.0, sig=3.0)
    assert len(result) == 3, "sep.winpos no longer returns (x, y, flag)"

    xw, yw, _flag = result
    assert float(xw) == pytest.approx(22.0, abs=0.5)
    assert float(yw) == pytest.approx(19.0, abs=0.5)


def test_measure_photometry_actually_refines_the_centroid(tmp_path):
    """End-to-end guard: the centroid must move off the ephemeris position."""
    from astropy.io import fits
    from astropy.wcs import WCS

    ny = nx = 81
    y, x = np.mgrid[0:ny, 0:nx]
    # Source deliberately offset by ~3 px from the frame centre.
    image = 1000.0 * np.exp(-((x - 43.0) ** 2 + (y - 38.0) ** 2) / (2 * 2.0 ** 2)) + 100.0

    header = fits.Header({
        "NAXIS": 2, "NAXIS1": nx, "NAXIS2": ny,
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CRPIX1": 41.0, "CRPIX2": 41.0, "CRVAL1": 100.0, "CRVAL2": 10.0,
        "CD1_1": -2.8e-4, "CD1_2": 0.0, "CD2_1": 0.0, "CD2_2": 2.8e-4,
    })
    path = tmp_path / "frame.fits"
    fits.writeto(path, image.astype(np.float32), header, overwrite=True)

    centre = WCS(header).pixel_to_world(40.0, 40.0)
    table = pd.DataFrame([{
        "file": "frame.fits", "obsjd": 2460000.5, "filter": "ZTF_r",
        "ra": centre.ra.deg, "dec": centre.dec.deg,
        "egain": 6.2, "readnoise": 9.7, "pixscale": 1.0, "fwhm_pix": 4.7,
        "delta": 1.0, "r": 2.0,
    }])

    out = phot.measure_photometry(table, tmp_path, cfg.PhotConfig(rho_km=1.0e6),
                                  progress=False)
    shift = float(out.loc[0, "centroid_shift_pix"])
    assert shift > 0.5, f"centroid was not refined (shift={shift})"
    assert out.loc[0, "x_center"] == pytest.approx(43.0, abs=1.0)
    assert out.loc[0, "y_center"] == pytest.approx(38.0, abs=1.0)


def test_regression_c1_zeropoint_is_applied_per_frame():
    """THE bug: ``inst_mag + row.zpmag`` broadcast one frame's MAGZP to all.

    Four frames with zeropoints 26.0/26.5/27.0/27.5 must produce four distinct
    magnitudes differing by exactly those offsets.
    """
    table = _synthetic_table()
    out = phot.calibrate(table, cfg.PhotConfig(apply_color_term=False,
                                               apply_aperture_correction=False))

    assert out["filter_mag"].nunique() == 4, "one zeropoint was used for every frame"
    np.testing.assert_allclose(out["filter_mag"].to_numpy(),
                               table["inst_mag"] + table["zpmag"])
    # And specifically not the last frame's zeropoint everywhere (the old bug).
    buggy = table["inst_mag"] + table["zpmag"].iloc[-1]
    assert not np.allclose(out["filter_mag"], buggy)


def test_regression_c4_colour_term_is_applied():
    table = _synthetic_table()
    without = phot.calibrate(table, cfg.PhotConfig(apply_color_term=False,
                                                   apply_aperture_correction=False))
    with_it = phot.calibrate(table, cfg.PhotConfig(apply_color_term=True,
                                                   apply_aperture_correction=False))

    shift = (with_it["filter_mag"] - without["filter_mag"]).abs()
    assert (shift > 0).all()
    # clrcoeff = -0.03 and the default colour 0.5 -> 0.015 mag
    assert shift.iloc[0] == pytest.approx(0.015, abs=1e-6)


def test_regression_c8_uncertainties_reach_afrho():
    table = phot.calibrate(_synthetic_table(), cfg.PhotConfig())
    out = phot.compute_afrho(table, cfg.PhotConfig())

    assert "afrho0_cm_err" in out
    assert (out["afrho0_cm_err"] > 0).all()
    # d(Afrho)/Afrho = 0.4 ln10 sigma_m
    expected = 0.4 * np.log(10) * out["filter_mag_err"]
    np.testing.assert_allclose(out["afrho_cm_err"] / out["afrho_cm"], expected, rtol=1e-9)
    # The zeropoint scatter must be in there, so the error exceeds photon noise alone.
    assert (out["filter_mag_err"] > table["inst_mag_err"]).all()


def test_afrho_matches_the_ahearn_definition():
    table = phot.calibrate(_synthetic_table(), cfg.PhotConfig())
    out = phot.compute_afrho(table, cfg.PhotConfig())
    row = out.iloc[0]

    expected = (4 * (row.delta * u.au) ** 2 * row.r ** 2 / (row.rho_km * u.km)
                * 10 ** (-0.4 * (row.filter_mag - row.solmag))).to_value(u.cm)
    assert row.afrho_cm == pytest.approx(expected, rel=1e-9)


def test_phase_correction_brightens_towards_zero_phase():
    table = phot.calibrate(_synthetic_table(), cfg.PhotConfig())
    out = phot.compute_afrho(table, cfg.PhotConfig(phase_beta=0.03))
    row = out.iloc[0]

    assert row.phase_corr == pytest.approx(10 ** (0.4 * 0.03 * row.alpha))
    assert row.afrho0_cm > row.afrho_cm       # alpha > 0, so A(0) exceeds A(alpha)


def test_quality_flags_do_not_remove_rows():
    """Flag, never drop — the user's explicit choice for this pipeline."""
    table = _synthetic_table()
    table.loc[0, "flag_sky_edge"] = True
    table.loc[1, "flag_undersampled"] = True

    out = phot.compute_afrho(phot.calibrate(table, cfg.PhotConfig()), cfg.PhotConfig())

    assert len(out) == len(table), "a flagged row was dropped"
    assert out["quality_ok"].tolist() == [False, False, True, True]
    assert out.loc[0, "afrho0_cm"] > 0, "flagged rows must keep their measurement"
    assert "sky_edge" in out.loc[0, "flags"]


def test_advisory_flags_do_not_reject_a_measurement():
    """An assumed colour costs ~0.015 mag; it must not exclude the whole night."""
    table = _synthetic_table(n=2)
    table["filter"] = ["ZTF_r", "ZTF_r"]          # single band -> colour assumed

    out = phot.compute_afrho(phot.calibrate(table, cfg.PhotConfig()), cfg.PhotConfig())

    assert out["flag_color_default"].all()
    assert out["quality_ok"].all(), "an advisory flag rejected a good measurement"
    assert set(phot.CRITICAL_FLAGS) & set(phot.ADVISORY_FLAGS) == set()
    assert set(phot.FLAG_COLUMNS) == set(phot.CRITICAL_FLAGS) | set(phot.ADVISORY_FLAGS)


def test_night_colour_measured_from_same_night_g_r_pair():
    table = _synthetic_table(n=2)
    table["obsjd"] = [2460000.9, 2460000.95]      # same night
    table["inst_mag"] = [-10.0, -10.4]
    table["zpmag"] = [26.0, 26.0]

    out = phot.calibrate(table, cfg.PhotConfig())
    assert not out["flag_color_default"].any()
    assert out["color_gr"].iloc[0] == pytest.approx(0.4, abs=1e-9)


def test_colour_falls_back_and_flags_when_only_one_band():
    table = _synthetic_table(n=2)
    table["filter"] = ["ZTF_r", "ZTF_r"]
    out = phot.calibrate(table, cfg.PhotConfig(default_color_gr=0.5))

    assert out["flag_color_default"].all()
    assert (out["color_gr"] == 0.5).all()


def test_solar_magnitudes_cover_every_ztf_band():
    assert set(cfg.SOLAR_APPMAG_AB) == set(cfg.ZTF_FILTERS.values())


def test_empty_input_survives_the_whole_chain():
    empty = pd.DataFrame()
    assert phot.calibrate(empty, cfg.PhotConfig()).empty
    assert phot.compute_afrho(empty, cfg.PhotConfig()).empty
    assert cutout.build_urls(empty).empty
