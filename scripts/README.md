# `scripts/` — batch drivers and standalone processing scripts

Run everything from the project root (`spherex-comet/`); each script puts the root on
`sys.path` itself.  Science logic lives in the packages — a script here only orchestrates,
so what runs in batch is what the notebooks validate.

## `ztf/` — stage 1 (`ztfcomet`)

| script | purpose |
|---|---|
| `main.py 24P [240P …] [--steps query download phot figures]` | end-to-end driver for one or several targets |
| `survey.py [--targets …] [--steps …] [--no-resume]` | unattended 68-comet batch, resumable from the status file on the data root; refuses to start on the local fallback |
| `afrho_trends.py [--figures-only] [--targets …]` | heliocentric Afρ trends per comet (peaks, breaks, outbursts, colour), the SPHEREx windows, the Afρ at every SPHEREx phase → `results/ztf/afrho/`, `fig/ztf/afrho/trend*/` |
| `afrho_trends_report.py` | rebuilds the generated tables of `doc/ztf/afrho_heliocentric_trends.md` |
| `spherex_windows.py` | JD / r_h range of every SPHEREx phase group (from the comspec phase assignment and the apphot tables) → `results/ztf/afrho/spherex_windows.csv`; called by `afrho_trends.py` |
| `profile_survey.py`, `profile_resolution.py` | the 1/ρ coma-law studies (`doc/ztf/profile_survey_1rho.md`, `profile_resolution_24P.md`) |
| `aperture_systematics.py` | the C/2024 E1 aperture study (`doc/ztf/afrho_aperture_systematics.md`) |
| `reflag_anomalies.py` | re-applies the single-frame anomaly flag to the photometry tables |

Needs the T7 SSD (`/Volumes/T7/data/ztf-comet`) or `ZTFCOMET_DATA`, and Gaia DR3 at
`ZTFCOMET_GAIA` for contamination flags.  Working list: `data/reference/sx_comet_list_ver2607.xlsx`.

## `apphot/` — stage 2 (`spherex_apphot`)

| script | purpose |
|---|---|
| `main.py --objdesig 24P` / `--target-list LIST --workers 4` / `--sample` | photometry only, never plots; resumes from `results/apphot/status.csv`; `--force` reruns |
| `make_figures.py … [--no-cutouts] [--max-cutouts N]` | per-target summary figures (and, on request, two PNGs per exposure: ~8 GB for the survey) → `fig/apphot/<slug>/` |
| `extract_targets.py --target-list LIST` | copies the working list's cutouts out of the flat 139k-file archive into `comets_v5_fits_filtered/<target>/` on the SSD (once per list; `--clean-appledouble`) |
| `build_gaia_cache.py` | one-off declination-sorted Gaia cache (~5 min, ~8 GB) — 100× faster queries |
| `h5cut2fits.py` | data preparation: expands the HDF5 cutout archive into individual FITS with full headers (needs `h5py`, `click` and the private `pqfilt`) |

Needs the T7 SSD (`/Volumes/T7/data/spherex-comet`) or `SPHEREX_T7_DIR` / `SPHEREX_DB` /
`SPHEREX_FITS_DIR`; `--sample` runs on `data/spherex_sample/`.  Outputs go to
`results/apphot/` (`SPHEREX_RESULT_DIR`) and `fig/apphot/` (`SPHEREX_FIG_DIR`).

## `comspec/` — stage 3 (`spherex_comspec`)

| script | purpose |
|---|---|
| `main.py all` / `group` / `run --variants dc_main [--targets …]` / `analyze` / `figures` | the catalog pipeline: phases → continuum → fits → cross-variant tables → figures (≈ 10 min for `all`) |
| `attach_afrho_ztf.py [--ztf-afrho PATH]` | pivots `results/ztf/afrho/spherex_afrho.csv` to `results/comspec/afrho_ztf.csv` and attaches the `afrho_*` columns to `gas_fit.csv` and `phase_map.csv` |
| `phase_images.py [--workers 4] [--all-phases] [24P:3 …]` | comet-centred north-up median stacks (ZTF r, SPHEREx continuum / H₂O / CO₂ / CO windows) per fitted phase → `fig/comspec/phase_images/`, `results/comspec/phase_stacks/` |
| `make_figure_slides.py [2P:1 …]` | the review deck `doc/figures.pptx`, three slides per fitted phase (needs `python-pptx`, `Pillow`) |

Reads `results/apphot/photometry/` (`COMSPEC_APPHOT_DIR`), writes `data/comspec/`,
`results/comspec/`, `fig/comspec/` (`COMSPEC_DATA_DIR`, `COMSPEC_RESULT_DIR`, `COMSPEC_FIG_DIR`
— set all three before import for an isolated run).  **After a regroup, rerun
`ztf/afrho_trends.py` first, then `attach_afrho_ztf.py`.**
