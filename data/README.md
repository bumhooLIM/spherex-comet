# `data/`

Inputs and intermediates.  Raw data live on the external SSD (see the root `README.md`);
only `reference/` and the READMEs are tracked in git.

| folder | content | written by |
|---|---|---|
| `reference/` | `sx_comet_list_ver2607.xlsx` (the 68-comet working list shared by all three stages), `comet_orbit_classes.csv` (dynamical classes), `gicquel2023_neowise_production_rates.csv` and `harrington_pinto2022_gas_production_rates.csv` (literature Q values, built by `notebooks/comspec/build_*_table.py`), `phase_update_map_previous.csv` (the previous phase grouping, for the regrouping regression test) | by hand / the builders |
| `fluorescence/` | the reconstructed GSFC-style fluorescence database, 0.7–5.0 µm: band g-factors, species summary, spectral profiles per T_rot (`profiles/`), line lists (`lines/`), the CO Swings table (`co_swings.csv`), the solar pump, validation, `build_meta.json`; `inputs/` holds the HITRAN and solar spectra it is built from.  Formats in its own `README.md`, physics in `doc/comspec/fluorescence_database.md` | `notebooks/comspec/fluorescence_gfm/build_fluorescence_db.py` (25 s) |
| `comspec/` | `phase_assignment.csv` (the phase of every SPHEREx exposure, shared by every variant), `emission/` (continuum summaries `<comet>_<aperture>.csv` and continuum-subtracted point spectra `*_points.csv` of the main variant), `studies/<variant>/emission/` (the same per study variant) | `scripts/comspec/main.py` |
| `spherex_sample/` | 2P (226) and 24P (765) SPHEREx cutout FITS plus their slice of `db_filtered.parq`, for the end-to-end test when the SSD is not mounted (`scripts/apphot/main.py --sample`, `notebooks/apphot/main.ipynb`) | copied from the archive |
| `ztf/` | the local fallback for ZTF raw FITS, `eph.csv`, `ztf.csv` when neither `ZTFCOMET_DATA` nor the SSD is available; empty by design (`survey.py` refuses to run here) | `ztfcomet` |

The SPHEREx aperture photometry that the catalog stage reads is **not** here: it is the
product of stage 2 and lives once, in `results/apphot/photometry/` (`spherex_comspec.directory.APPHOT_DIR`).
Its columns are documented in `results/README.md`.
