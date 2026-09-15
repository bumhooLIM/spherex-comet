"""Tests for Gaia contamination flagging and Horizons fragment handling.

Offline: the Gaia tests build a small synthetic catalogue on disk in the same
layout as the real declination-sorted cache, and the Horizons tests parse a
captured ambiguity listing rather than calling the service.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ztfcomet import config as cfg
from ztfcomet import gaia, horizons, phot

# The real listing Horizons returns for "240P", captured 2026-09.
AMBIGUITY_240P = """Ambiguous target name; provide unique id:
    Record #  Epoch-yr  >MATCH DESIG<  Primary Desig  Name
    --------  --------  -------------  -------------  -------------------------
    90001211    2014    240P           240P            NEAT
    90001212    2024    240P           240P            NEAT
    90001213    2025    240P-B         240P-B          NEAT"""


# ====================================================================== Gaia
def test_effective_magnitude_matches_the_specified_example():
    """Two G=15 sources combine to G_eff ~ 14.25."""
    assert gaia.effective_magnitude([15.0, 15.0]) == pytest.approx(14.2474, abs=1e-4)


def test_effective_magnitude_of_a_single_source_is_itself():
    assert gaia.effective_magnitude([16.3]) == pytest.approx(16.3)


def test_effective_magnitude_of_nothing_is_infinitely_faint():
    assert gaia.effective_magnitude([]) == np.inf
    assert gaia.effective_magnitude([np.nan]) == np.inf


def test_effective_magnitude_is_dominated_by_the_brightest_source():
    bright_only = gaia.effective_magnitude([12.0])
    with_faint = gaia.effective_magnitude([12.0, 18.0, 18.4])
    assert with_faint < bright_only                     # brighter (smaller mag)
    assert with_faint == pytest.approx(bright_only, abs=0.02)


@pytest.mark.parametrize("contaminant,comet,expected", [
    (15.0, 15.0, 1.0),          # equal flux
    (16.31, 15.0, 0.30),        # the 30% threshold, +1.31 mag
    (np.inf, 15.0, 0.0),        # nothing there
])
def test_flux_ratio(contaminant, comet, expected):
    assert gaia.flux_ratio(contaminant, comet) == pytest.approx(expected, abs=0.005)


def test_thirty_percent_threshold_is_1_31_magnitudes():
    """The stated criterion in magnitude terms, so the constant is pinned."""
    assert -2.5 * np.log10(0.30) == pytest.approx(1.3072, abs=1e-4)


def _write_catalogue(tmp_path, ra, dec, gmag):
    """Build a synthetic dec-sorted cache in the real on-disk layout."""
    cache = tmp_path / "gaiadr3_deccache"
    cache.mkdir(parents=True, exist_ok=True)
    order = np.argsort(dec)
    np.save(cache / "ra.npy", np.asarray(ra, float)[order])
    np.save(cache / "dec.npy", np.asarray(dec, float)[order])
    np.save(cache / "gmag.npy", np.asarray(gmag, np.float32)[order])
    (cache / "meta.json").write_text(json.dumps(
        {"sorted_by": "dec", "n_rows": len(ra)}))
    return gaia.GaiaCatalog(path=tmp_path)


def test_cone_search_finds_only_sources_inside_the_radius(tmp_path):
    # 0.005 deg = 18", 0.02 deg = 72"
    cat = _write_catalogue(tmp_path,
                           ra=[100.0, 100.0, 100.0, 200.0],
                           dec=[10.0, 10.005, 10.02, -30.0],
                           gmag=[15.0, 16.0, 12.0, 11.0])
    assert cat.has_cache and cat.available

    gmag, sep = cat.cone_search(100.0, 10.0, radius_arcsec=30.0)
    assert sorted(np.round(gmag, 1)) == [15.0, 16.0]     # the 0.02 deg one is out
    assert sep.max() < 30.0


def test_cone_search_returns_sources_brightest_first(tmp_path):
    cat = _write_catalogue(tmp_path, ra=[50.0] * 3, dec=[5.0, 5.001, 5.002],
                           gmag=[17.0, 13.0, 15.0])
    gmag, _ = cat.cone_search(50.0, 5.0, radius_arcsec=60.0)
    assert list(np.round(gmag, 1)) == [13.0, 15.0, 17.0]


def test_cone_search_handles_the_ra_wrap_at_zero(tmp_path):
    """A cone at RA=0 must find sources at RA=359.99."""
    cat = _write_catalogue(tmp_path, ra=[359.995, 0.005], dec=[0.0, 0.0],
                           gmag=[14.0, 15.0])
    gmag, _ = cat.cone_search(0.0, 0.0, radius_arcsec=60.0)
    assert len(gmag) == 2


def test_cone_search_near_the_pole_does_not_explode(tmp_path):
    """cos(dec) -> 0 at the pole; the flat-sky RA window must be guarded."""
    cat = _write_catalogue(tmp_path, ra=[0.0, 180.0], dec=[89.999, 89.999],
                           gmag=[14.0, 15.0])
    gmag, sep = cat.cone_search(0.0, 89.999, radius_arcsec=30.0)
    assert np.all(np.isfinite(sep))
    assert len(gmag) <= 2                                # never raises


def test_cone_search_respects_a_magnitude_limit(tmp_path):
    cat = _write_catalogue(tmp_path, ra=[50.0] * 2, dec=[5.0, 5.001],
                           gmag=[14.0, 18.2])
    gmag, _ = cat.cone_search(50.0, 5.0, radius_arcsec=60.0, mag_limit=17.0)
    assert list(np.round(gmag, 1)) == [14.0]


def test_missing_catalogue_is_reported_not_raised(tmp_path):
    cat = gaia.GaiaCatalog(path=tmp_path / "nothing_here")
    assert not cat.available
    assert "NOT FOUND" in cat.describe()
    gmag, sep = cat.cone_search(10.0, 10.0, 30.0)
    assert len(gmag) == 0 and len(sep) == 0


def test_check_aperture_flags_a_bright_intruder(tmp_path):
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[15.0])
    result = cat.check_aperture(100.0, 10.0, radius_arcsec=30.0, comet_mag=17.0)

    assert result.n_sources == 1
    assert result.g_eff == pytest.approx(15.0)
    assert result.ratio == pytest.approx(10 ** (0.4 * 2.0), rel=1e-6)   # ~6.3x
    assert result.contaminated


def test_check_aperture_passes_a_clean_field(tmp_path):
    """A source 100x fainter than the comet is 1%, well under the threshold."""
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[20.0])
    result = cat.check_aperture(100.0, 10.0, radius_arcsec=30.0, comet_mag=15.0)
    assert not result.contaminated
    assert result.ratio < 0.30


def test_check_aperture_sits_exactly_on_the_threshold(tmp_path):
    """G_eff = comet + 1.31 is the 30% boundary; just brighter must flag."""
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[16.25])
    assert cat.check_aperture(100.0, 10.0, 30.0, comet_mag=15.0).contaminated

    cat2 = _write_catalogue(tmp_path / "b", ra=[100.0], dec=[10.0], gmag=[16.40])
    assert not cat2.check_aperture(100.0, 10.0, 30.0, comet_mag=15.0).contaminated


def test_two_faint_sources_can_together_exceed_the_threshold(tmp_path):
    """The criterion is on combined flux, not the brightest source alone."""
    cat = _write_catalogue(tmp_path, ra=[100.0, 100.0], dec=[10.0, 10.001],
                           gmag=[16.9, 16.9])
    result = cat.check_aperture(100.0, 10.0, radius_arcsec=30.0, comet_mag=15.0)

    assert gaia.flux_ratio(16.9, 15.0) < 0.30            # neither alone
    assert result.g_eff == pytest.approx(16.9 - 0.7526, abs=1e-3)
    assert result.contaminated                            # but together they do


# ------------------------------------------------------- integration with phot
def _phot_table(n=3):
    return pd.DataFrame({
        "obsjd": 2460000.5 + np.arange(n),
        "file": [f"f{i}.fits" for i in range(n)],
        "filter": ["ZTF_r"] * n,
        "ra": np.full(n, 100.0), "dec": [10.0, 40.0, 70.0],
        "rho_pix": np.full(n, 5.0), "fwhm_pix": np.full(n, 3.0),
        "pixscale": np.full(n, 1.0), "tmag": np.full(n, 17.0),
    })


def test_flag_contamination_marks_only_the_affected_frame(tmp_path):
    # A bright star sits at the first frame's position only.
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[14.0])
    out = phot.flag_contamination(_phot_table(), cfg.PhotConfig(), catalogue=cat,
                                  progress=False)

    assert out["flag_contaminated"].tolist() == [True, False, False]
    assert out.loc[0, "contam_n_sources"] == 1
    assert out.loc[1, "contam_ratio"] == 0.0
    assert len(out) == 3, "contamination flagging must not drop rows"


def test_contamination_radius_includes_the_seeing_pad(tmp_path):
    """radius = (rho_pix + pad*FWHM) * pixscale, so PSF wings are covered."""
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[14.0])
    out = phot.flag_contamination(_phot_table(), cfg.PhotConfig(contam_radius_pad_fwhm=1.0),
                                  catalogue=cat, progress=False)
    assert out.loc[0, "contam_radius_arcsec"] == pytest.approx(8.0)   # (5+1*3)*1.0

    out2 = phot.flag_contamination(_phot_table(), cfg.PhotConfig(contam_radius_pad_fwhm=0.0),
                                   catalogue=cat, progress=False)
    assert out2.loc[0, "contam_radius_arcsec"] == pytest.approx(5.0)


def test_contamination_can_be_switched_off(tmp_path):
    cat = _write_catalogue(tmp_path, ra=[100.0], dec=[10.0], gmag=[10.0])
    out = phot.flag_contamination(_phot_table(), cfg.PhotConfig(check_contamination=False),
                                  catalogue=cat, progress=False)
    assert not out["flag_contaminated"].any()


def test_contamination_without_a_catalogue_leaves_frames_unflagged(tmp_path):
    cat = gaia.GaiaCatalog(path=tmp_path / "absent")
    out = phot.flag_contamination(_phot_table(), cfg.PhotConfig(), catalogue=cat,
                                  progress=False)
    assert not out["flag_contaminated"].any()
    assert len(out) == 3


def test_contamination_is_a_critical_flag():
    """A star in the aperture invalidates the measurement, so it must gate quality_ok."""
    assert "flag_contaminated" in phot.CRITICAL_FLAGS
    assert "flag_contaminated" not in phot.ADVISORY_FLAGS


# ================================================================== Horizons
def test_parse_the_real_240p_ambiguity_listing():
    records = horizons.parse_ambiguity_table(AMBIGUITY_240P)
    assert [r.record for r in records] == [90001211, 90001212, 90001213]
    assert [r.epoch_yr for r in records] == [2014, 2024, 2025]
    assert [r.is_fragment for r in records] == [False, False, True]
    assert all(r.name == "NEAT" for r in records)


def test_parse_returns_nothing_for_a_non_ambiguity_error():
    assert horizons.parse_ambiguity_table("Unable to connect to Horizons") == []


@pytest.mark.parametrize("designation,expected", [
    ("240P-B", True), ("73P-C", True), ("73P-BB", True), ("C/2019 Y4-A", True),
    ("240P", False), ("73P", False), ("C/2024 E1", False), ("2P", False),
])
def test_fragment_designation_detection(designation, expected):
    assert horizons.is_fragment_designation(designation) is expected


@pytest.mark.parametrize("designation,parent,letter", [
    ("240P-B", "240P", "B"), ("73P-BB", "73P", "BB"), ("240P", "240P", None),
])
def test_split_designation(designation, parent, letter):
    assert horizons.split_designation(designation) == (parent, letter)


def test_regression_240p_resolves_to_the_parent_not_the_fragment():
    """THE bug: 'take the last record' selected 240P-B, a fragment.

    240P-B is ~2.9 mag fainter than the parent, so the whole reduction would
    have been of the wrong object.
    """
    records = horizons.parse_ambiguity_table(AMBIGUITY_240P)
    chosen = horizons.select_record(records, epoch_jd=2460900.5, designation="240P")

    assert chosen.record != 90001213, "selected the fragment 240P-B"
    assert not chosen.is_fragment
    assert chosen.primary_desig == "240P"


def test_orbit_solution_nearest_the_observation_is_preferred():
    """The two parent solutions differ by ~50 arcsec in 2025, so epoch matters."""
    records = horizons.parse_ambiguity_table(AMBIGUITY_240P)

    # 2018-06 -> the 2014 solution; 2025-09 -> the 2024 solution
    assert horizons.select_record(records, epoch_jd=2458300.5,
                                  designation="240P").record == 90001211
    assert horizons.select_record(records, epoch_jd=2460900.5,
                                  designation="240P").record == 90001212


def test_fragment_is_selectable_when_explicitly_requested():
    records = horizons.parse_ambiguity_table(AMBIGUITY_240P)
    chosen = horizons.select_record(records, epoch_jd=2460900.5,
                                    designation="240P-B", allow_fragment=True)
    assert chosen.record == 90001213 and chosen.is_fragment


def test_designation_filter_excludes_unrelated_records():
    records = horizons.parse_ambiguity_table(AMBIGUITY_240P) + [
        horizons.HorizonsRecord(90009999, 2024, "241P", "241P", "OTHER")]
    chosen = horizons.select_record(records, epoch_jd=2460900.5, designation="240P")
    assert chosen.primary_desig == "240P"


@pytest.mark.parametrize("returned,requested,ok", [
    ("240P/NEAT", "240P", True),
    ("240P-B/NEAT", "240P", False),      # the fragment substitution
    ("234P/LINEAR", "240P", False),      # a stale record number
    ("233P/La Sagra", "240P", False),
])
def test_verify_targetname_catches_wrong_objects(returned, requested, ok):
    assert horizons.verify_targetname(returned, requested)[0] is ok


def test_verify_targetname_accepts_a_deliberately_requested_fragment():
    assert horizons.verify_targetname("240P-B/NEAT", "240P-B")[0] is True
    assert horizons.verify_targetname("240P-B/NEAT", "240P", allow_fragment=True)[0] is True


def test_regression_designations_are_queried_as_small_bodies(monkeypatch):
    """A bare designation must not be matched against major bodies.

    Horizons guesses the object class. Left to guess, ``"2P"`` resolves to
    **Styx (905)**, a moon of Pluto, and returns a full, plausible-looking
    ephemeris for it — no error, no ambiguity listing, just the wrong object.
    Every lookup must force ``id_type="smallbody"``.
    """
    seen = {}

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            seen["id"], seen["id_type"] = id, id_type
        def ephemerides(self, **k):
            raise ValueError("Ambiguous target name\n"
                             "  90000091    2022    2P             2P    Encke\n")

    monkeypatch.setattr(horizons, "Horizons", Fake)
    horizons.clear_cache()
    horizons.resolve_record("2P", epoch_jd=2461000.0)

    assert seen["id_type"] == horizons.SMALLBODY, \
        "designation probed without id_type='smallbody'; 2P would resolve to Styx"


def test_ephemeris_queries_also_force_smallbody(monkeypatch):
    from ztfcomet import query as qmod

    seen = {}

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            seen["id_type"] = id_type
            self._n = len(epochs)
        def ephemerides(self, **k):
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Fake)
    qmod.query_sso_ephemeris("2P", epochs=[2461000.0])
    assert seen["id_type"] == horizons.SMALLBODY


def test_2p_apparition_listing_selects_a_recent_solution():
    """2P/Encke has ~60 apparition records going back to 1786."""
    listing = "Ambiguous target name; provide unique id:\n" + "\n".join(
        f"    {90000031 + i}    {1786 + i * 4}    2P             2P              Encke"
        for i in range(60))
    records = horizons.parse_ambiguity_table(listing)
    assert len(records) == 60

    chosen = horizons.select_record(records, epoch_jd=2461000.0, designation="2P")
    assert chosen.epoch_yr >= 2000, "picked a 19th-century orbit solution"
    assert not chosen.is_fragment


def test_explicit_record_numbers_pass_through_unchanged():
    assert horizons.resolve_record(90001212) is None
    assert horizons.resolve_target_id(90001212) == 90001212


def test_config_240p_uses_a_designation_not_a_stale_record():
    """90001203/90001204 now resolve to 233P and 234P — different comets."""
    target = cfg.get_target("240P")
    assert target.query_designation == "240P"
    assert not target.allow_fragment
    assert not target.orbit_records, "pinned record numbers go stale"


def test_config_exposes_the_fragment_as_its_own_target():
    fragment = cfg.get_target("240P-B")
    assert fragment.allow_fragment
    assert horizons.is_fragment_designation(fragment.query_designation)


def test_search_frames_uses_the_designation_not_the_empty_horizons_id(monkeypatch):
    """Targets carry a designation; ``horizons_id`` is None for all of them.

    ``search_frames`` read ``target.horizons_id`` directly and passed None to
    Horizons, which aborted every query with "'id' parameter not set" — a whole
    run producing zero frames and only a warning.
    """
    from ztfcomet import query as qmod

    seen = {}

    def fake_eph(target_id, epochs, **kw):
        seen["target_id"] = target_id
        return pd.DataFrame()

    monkeypatch.setattr(qmod, "query_sso_ephemeris", fake_eph)
    target = cfg.get_target("24P")
    qmod.search_frames(target, progress=False)

    assert seen["target_id"] is not None, "passed None as the Horizons id"
    assert str(seen["target_id"]) not in ("", "None")


def test_ambiguity_fallback_inside_query_is_fragment_aware(monkeypatch):
    """The retry path inside query_sso_ephemeris must not take the last record."""
    from ztfcomet import query as qmod

    tried = []

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            tried.append(id)
            self._id = id
        def ephemerides(self, **k):
            if str(self._id) == "240P":
                raise ValueError(AMBIGUITY_240P)
            class T:
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": [2460900.5]})
            return T()

    monkeypatch.setattr(qmod, "Horizons", Fake)
    qmod.query_sso_ephemeris("240P", epochs=[2460900.5])

    assert 90001213 not in tried, "retried with the fragment 240P-B"
    assert 90001212 in tried or 90001211 in tried


def test_transient_horizons_valueerror_is_retried_not_fatal(monkeypatch):
    """astroquery raises ValueError for a service failure as well as ambiguity.

    Only the ambiguous case carries a candidate listing. Treating the other as
    "no such object" silently dropped 4 of 65 epochs in a real 24P run, logged
    as "No valid Horizons record for 90000356".
    """
    from ztfcomet import query as qmod

    calls = {"n": 0}

    class Flaky:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            self._n = len(epochs)
        def ephemerides(self, **k):
            calls["n"] += 1
            if calls["n"] < 3:
                raise ValueError("Query failed without known error message; "
                                 "received the following response:\n<html>")
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Flaky)
    out = qmod.query_sso_ephemeris(90000356, epochs=[2460900.5], backoff=0.0)

    assert calls["n"] == 3, "did not retry a transient failure"
    assert not out.empty, "gave up on a recoverable error"


def test_unparsable_valueerror_does_not_abort_other_chunks(monkeypatch):
    """One bad chunk must not discard the chunks that succeeded."""
    from ztfcomet import query as qmod

    seen = {"n": 0}

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            self._n = len(epochs)
        def ephemerides(self, **k):
            seen["n"] += 1
            if seen["n"] <= 3:          # first chunk fails all 3 attempts
                raise ValueError("Query failed without known error message")
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Fake)
    out = qmod.query_sso_ephemeris(90000356, epochs=list(range(100)),
                                   max_epochs_per_call=50, backoff=0.0)
    assert len(out) == 50, "surviving chunk was discarded with the failing one"


def test_regression_all_identical_epochs_are_deduplicated(monkeypatch):
    """Horizons rejects a TLIST whose entries are ALL identical.

    It answers "Bad dates -- start must be earlier than stop", which is what
    happens when one IRSA step returns several frames from a single exposure
    (the comet landing on two CCD quadrants of the same image). That killed
    4 of 65 steps in a real 24P run.
    """
    from ztfcomet import query as qmod

    sent = []

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            sent.append(list(epochs))
            self._n = len(epochs)
        def ephemerides(self, **k):
            if len(set(sent[-1])) == 1 and len(sent[-1]) > 1:
                raise ValueError("Bad dates -- start must be earlier than stop")
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Fake)
    out = qmod.query_sso_ephemeris(90000356, epochs=[2460900.5] * 3, backoff=0.0)

    assert sent[0] == [2460900.5], "duplicate epochs were sent to Horizons"
    assert len(out) == 1


def test_epoch_deduplication_preserves_distinct_values(monkeypatch):
    from ztfcomet import query as qmod

    sent = []

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            sent.append(list(epochs)); self._n = len(epochs)
        def ephemerides(self, **k):
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Fake)
    qmod.query_sso_ephemeris(1, epochs=[2460901.5, 2460900.5, 2460900.5])
    assert sent[0] == [2460900.5, 2460901.5], "dedup lost a distinct epoch or the sort"


def test_non_finite_epochs_are_dropped(monkeypatch):
    """A NaN in the TLIST makes Horizons fail with an opaque loader error."""
    from ztfcomet import query as qmod

    sent = []

    class Fake:
        def __init__(self, id=None, id_type=None, location=None, epochs=None):
            sent.append(list(epochs)); self._n = len(epochs)
        def ephemerides(self, **k):
            class T:
                def __init__(self, n): self._n = n
                def to_pandas(self):
                    return pd.DataFrame({"datetime_jd": np.arange(self._n)})
            return T(self._n)

    monkeypatch.setattr(qmod, "Horizons", Fake)
    qmod.query_sso_ephemeris(1, epochs=[2460900.5, float("nan")])
    assert sent[0] == [2460900.5]


@pytest.mark.parametrize("returned,requested,ok", [
    # Horizons clips targetname at 31 chars, losing the designation entirely.
    ("Bernardinelli-Bernstein (C/2014", "2014 UN271", True),
    ("Tsuchinshan-ATLAS (C/2023 A3)", "2023 A3", True),
    ("Ye (P/2025 UX109)", "2025 UX109", True),
    ("ATLAS (C/2022 QE78)", "2022 QE78", True),
    ("29P/Schwassmann-Wachmann 1", "29P", True),      # hyphen in the NAME, not a fragment
    ("PANSTARRS (C/2017 K2)", "2023 A3", False),      # genuinely the wrong comet
    ("Styx (905)", "2P", False),
    ("240P-B/NEAT", "240P", False),
])
def test_verify_targetname_handles_truncation_and_long_names(returned, requested, ok):
    """A clipped name must not fail verification, but a wrong one still must."""
    assert horizons.verify_targetname(returned, requested)[0] is ok


def test_aperture_scale_test_bounds():
    """FWHM < r_ap < 1 arcmin. Both bounds must actually bind."""
    from ztfcomet import phot as ph

    # delta = 1.5 au, 1.012"/px, FWHM 2.85 px -> 2.9"
    ok, r, fwhm = ph.aperture_scale_ok(10_000, 1.5, 1.012, 2.85)
    assert ok and fwhm < r < 60

    # Too large: 40000 km at 0.3 au subtends 184"
    ok, r, _ = ph.aperture_scale_ok(40_000, 0.3, 1.012, 2.85)
    assert not ok and r > 60

    # Too small: 1000 km at 5 au is well under the seeing disc
    ok, r, fwhm = ph.aperture_scale_ok(1_000, 5.0, 1.012, 2.85)
    assert not ok and r < fwhm


def test_aperture_scale_test_rejects_bad_input():
    from ztfcomet import phot as ph
    assert ph.aperture_scale_ok(np.nan, 1.0, 1.0, 3.0)[0] is False
    assert ph.aperture_scale_ok(10_000, 0.0, 1.0, 3.0)[0] is False


def test_signed_rh_separates_the_orbital_legs():
    """Pre-perihelion goes negative so the two legs do not overlap."""
    from ztfcomet import plotting as plotting_mod

    table = pd.DataFrame({"r": [2.0, 1.5, 1.5, 2.0],
                          "r_rate": [-5.0, -1.0, 1.0, 5.0]})
    x = plotting_mod.signed_rh(table)
    assert list(np.sign(x)) == [-1, -1, 1, 1]
    assert list(np.abs(x)) == [2.0, 1.5, 1.5, 2.0]


def test_signed_rh_falls_back_without_r_rate():
    from ztfcomet import plotting as plotting_mod
    table = pd.DataFrame({"r": [2.0, 1.5]})
    assert list(plotting_mod.signed_rh(table)) == [2.0, 1.5]


def _leg_table(rates, rho=15000.0):
    n = len(rates)
    return pd.DataFrame({
        "r": np.linspace(3.0, 1.5, n), "r_rate": rates,
        "filter": ["ZTF_r"] * n, "rho_km": [rho] * n,
        "afrho0_cm": np.linspace(50, 300, n), "afrho0_cm_err": [5.0] * n,
        "quality_ok": [True] * n,
    })


def test_single_leg_plot_keeps_rh_positive():
    """Signing r_h with only one orbital leg would make every x negative.

    10P is entirely inbound; the first survey figure showed r_h running
    -3.9 to -1.4 under an axis labelled r_h.
    """
    import matplotlib
    matplotlib.use("Agg")
    from ztfcomet import plotting as pl

    # No elements: falls back to plain r_h, which must stay positive.
    ax = pl.plot_afrho_vs_rh({"10P": _leg_table([-5.0] * 6)}, rho_km=15000.0)
    xs = _plotted_x(ax)
    assert (xs > 0).all(), "single-leg plot rendered negative r_h"
    assert "pre-perihelion" in ax.get_xlabel()
    plt_close(ax)


def test_two_leg_plot_signs_the_inbound_branch():
    import matplotlib
    matplotlib.use("Agg")
    from ztfcomet import plotting as pl

    from ztfcomet import orbit

    el = orbit.PerihelionInfo(q=1.2, e=0.7, Tp_jd=2461048.8)
    ax = pl.plot_afrho_vs_rh({"24P": _leg_table([-5.0, -5.0, -5.0, 5.0, 5.0, 5.0])},
                             rho_km=15000.0, elements=el)
    xs = _plotted_x(ax)
    assert (xs < 0).any() and (xs > 0).any(), "two-leg plot did not separate the legs"
    plt_close(ax)


def plt_close(ax):
    import matplotlib.pyplot as plt
    plt.close(ax.figure)


def _plotted_x(ax):
    """x of the drawn data points, excluding the perihelion marker line."""
    xs = [ln.get_xdata() for ln in ax.lines
          if len(ln.get_xdata()) and ln.get_marker() not in ("", "None", None)]
    return np.concatenate(xs) if xs else np.array([])


# ---------------------------------------------------------------- orbit / axis
def test_kepler_time_from_perihelion_matches_elliptic_geometry():
    """t(r_h) from Kepler, checked against the definition at known points."""
    from ztfcomet import orbit

    q, e = 1.1839, 0.7083                       # 24P/Schaumasse
    # At perihelion the time offset is zero.
    assert orbit.time_from_perihelion([q], q, e)[0] == pytest.approx(0.0, abs=1e-6)
    # Monotonic: further out is longer from perihelion.
    t = orbit.time_from_perihelion([q + 0.1, q + 0.5, q + 1.0], q, e)
    assert list(t) == sorted(t)
    # Aphelion of this orbit is a(1+e); half the period from perihelion.
    a = q / (1 - e)
    period = 2 * np.pi / (orbit.GAUSS_K / a ** 1.5)
    assert orbit.time_from_perihelion([a * (1 + e)], q, e)[0] == pytest.approx(
        period / 2, rel=1e-6)


def test_kepler_handles_all_three_conic_cases():
    from ztfcomet import orbit

    for e in (0.9, 0.9999, 1.0, 1.5, 3.0):      # elliptic, near-parabolic, hyperbolic
        t = orbit.time_from_perihelion([2.0, 5.0], q=1.0, e=e)
        assert np.all(np.isfinite(t)), f"e={e} produced non-finite times"
        assert t[1] > t[0] > 0, f"e={e} not monotonic in r_h"


def test_beyond_aphelion_is_unreachable():
    """A closed orbit never exceeds a(1+e); such radii must not saturate."""
    from ztfcomet import orbit
    q, e = 1.0, 0.3                                  # aphelion = 1.857 au
    t = orbit.time_from_perihelion([1.5, 5.0], q, e)
    assert np.isfinite(t[0])
    assert np.isnan(t[1]), "beyond aphelion silently clamped to half a period"


def test_time_from_perihelion_signs_by_leg():
    from ztfcomet import orbit
    t = orbit.time_from_perihelion([2.0, 2.0], q=1.0, e=0.5, signed_by=[-1.0, 1.0])
    assert t[0] < 0 < t[1] and abs(t[0]) == pytest.approx(abs(t[1]))


def test_inside_perihelion_is_unreachable_not_nan():
    from ztfcomet import orbit
    t = orbit.time_from_perihelion([0.5], q=1.0, e=0.5)
    assert np.isfinite(t[0]) and t[0] == 0.0


def test_afrho_axis_is_offset_from_perihelion_not_signed_rh():
    """x must be r_h - q, so no point sits in the unreachable band below q."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ztfcomet import orbit, plotting as pl

    q = 1.2
    table = _leg_table([-5.0, -5.0, 5.0, 5.0])
    table["r"] = [2.0, 1.5, 1.5, 2.0]
    el = orbit.PerihelionInfo(q=q, e=0.7, Tp_jd=2461048.8)

    ax = pl.plot_afrho_vs_rh({"t": table}, rho_km=15000.0, elements=el)
    xs = _plotted_x(ax)
    assert np.isclose(np.abs(xs).min(), 0.3)        # 1.5 - 1.2
    assert (xs < 0).any() and (xs > 0).any()
    assert "q$" in ax.get_xlabel() or "q" in ax.get_xlabel()
    plt.close(ax.figure)



def test_target_names_survive_a_csv_round_trip(tmp_path):
    """"2024E1" is valid scientific notation; pandas turns it into 20240.0.

    The trap bites only when the column is homogeneous -- one target's own
    photometry or profile table -- because a mixed column stays object dtype.
    survey_status.csv is the resume checkpoint, so a float-parsed name means
    the target is not recognised as done and gets reprocessed, silently.
    """
    single = tmp_path / "photometry_2024E1.csv"
    pd.DataFrame({"target": ["2024E1"] * 3, "obsjd": [1.0, 2.0, 3.0]}).to_csv(single, index=False)
    naive = pd.read_csv(single)
    assert naive["target"].iloc[0] != "2024E1", "expected pandas to mangle a homogeneous column"
    assert pd.read_csv(single, dtype={"target": str})["target"].iloc[0] == "2024E1"

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts" / "ztf"))
    import survey
    for names in (["2024E1"], ["2024E1", "2022E2", "24P"]):
        pd.DataFrame({"target": names, "status": ["ok"] * len(names)}).to_csv(
            tmp_path / "survey_status.csv", index=False)
        assert set(survey.load_status(tmp_path)) == set(names)
