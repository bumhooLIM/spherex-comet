# `tests/` — the offline test suites

```bash
pytest tests/ -q                    # everything (226 tests, ~20 s); pythonpath is set in pyproject.toml
pytest tests/test_comspec.py -q     # one stage
python tests/test_apphot_pipeline.py   # the apphot suite also runs without pytest
```

| file | package | what it pins |
|---|---|---|
| `test_ztfcomet.py` | `ztfcomet` | paths, query, cutout validation, photometry, calibration, Afρ; the `test_regression_c*` cases pin the defects of `doc/ztf/primitive_code_analysis.md` (per-frame zeropoint, JD joins, FITS-signature downloads, `sep.winpos` three values, …) |
| `test_contamination.py` | `ztfcomet` | Gaia contamination flagging, Horizons resolution (fragments, renumbered records), the survey status reader |
| `test_profile.py` | `ztfcomet` | radial profiles, star stacking, power-law fits, centring |
| `test_activity.py` | `ztfcomet` | trend fits, peaks, breaks, outbursts, the Afρ at a SPHEREx epoch |
| `test_apphot_pipeline.py` | `spherex_apphot` | bit masks, aperture sets, sky and error budget, source flags, stacking (median vs clipped mean), reflectance, slugs — no data files needed |
| `test_comspec.py` | `spherex_comspec` | grouping, continuum, design matrix and fits, fluorescence-database wiring; the reproduction tests use `results/apphot/photometry/` and `data/reference/phase_update_map_previous.csv` and are skipped when absent |

All suites are deterministic and offline (no network, no SSD).  One ZTF test
(`test_measure_photometry_actually_refines_the_centroid`) needs more than 3 GB in a bare
Linux VM because of the sep/photutils build there; it passes in the `spherex` environment.
