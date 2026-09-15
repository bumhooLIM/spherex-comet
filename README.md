# spherex-comet

Production of the **SPHEREx comet spectrum catalog**: gas production rates of H₂O, CO₂ and
CO for 68 comets from SPHEREx spectrophotometry, with the dust activity of the same epochs
from ZTF imaging.  Three stages, three packages, one tree:

| stage | package | does | main product |
|---|---|---|---|
| 1 | `ztfcomet/` | ZTF cutout query + download (IRSA, Horizons), Afρ photometry, coma radial profiles, heliocentric Afρ trends, Afρ at every SPHEREx phase | `results/ztf/afrho/spherex_afrho.csv` |
| 2 | `spherex_apphot/` | aperture photometry of the SPHEREx cutouts: bad-pixel masks, physical sky annulus, 21 apertures, Gaia contamination flags, band stacks | `results/apphot/photometry/<comet>.csv` |
| 3 | `spherex_comspec/` | phase grouping, continuum subtraction, joint H₂O/CO₂/CO fluorescence fit, study variants, figures, review deck | **`results/comspec/gas_fit.csv`** |

`doc/summary.md` is the scientific overview (background, concepts, pipeline, results, open
items); `handoff.md` is the working state for the next session; `CLAUDE.md` holds the
conventions.  The tree was merged on 2026-09-15 from the former `ztf-comet`,
`spherex-comet-apphot` and `spherex-comet-catalog` projects.

---

## Install

Python ≥ 3.10 in the `spherex` conda environment.

```bash
conda activate spherex
pip install -e .            # the three packages; optional: pip install -e ".[dev,deck]"
pytest tests/ -q            # 226 offline tests
```

Nothing needs to be installed to run from the project root: every driver and notebook puts
the root on `sys.path`.  Raw data are **not** in the tree (see *Data* below).

## Quick start

```bash
# stage 1 — ZTF (needs the T7 SSD, or ZTFCOMET_DATA, and Gaia DR3 for contamination flags)
python scripts/ztf/main.py 24P                                   # query, download, reduce, plot one comet
python scripts/ztf/survey.py                                     # all 68, resumable, unattended
python scripts/ztf/afrho_trends.py                               # trends + Afρ at the SPHEREx phases

# stage 2 — SPHEREx photometry (needs the T7 SSD; --sample uses data/spherex_sample/)
python scripts/apphot/main.py --sample --objdesig 2P 24P --stack # end-to-end check, ~1 min
python scripts/apphot/main.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 4
python scripts/apphot/make_figures.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 6 --no-cutouts

# stage 3 — production rates (runs on results/apphot/photometry/, ~10 min)
python scripts/comspec/main.py all                               # group -> variants -> analyze -> figures
python scripts/comspec/attach_afrho_ztf.py                       # ZTF dust context into gas_fit.csv
python scripts/comspec/phase_images.py --workers 4               # stacked images per fitted phase
python scripts/comspec/make_figure_slides.py                     # doc/figures.pptx
```

Each stage's `README.md` (`ztfcomet/`, `spherex_apphot/`, `spherex_comspec/`) is its tool
guide: options, configuration, output columns.

## Layout

```
ztfcomet/  spherex_apphot/  spherex_comspec/    the packages (README.md in each)
scripts/{ztf,apphot,comspec}/                   batch drivers and standalone scripts
notebooks/rcparams.py, {ztf,apphot,comspec}/     figure style; validation notebooks; study scripts
tests/                                          the three offline suites
data/                                           reference/ · fluorescence/ · comspec/ · spherex_sample/ · ztf/
results/{ztf,apphot,comspec}/                   main clean products of each stage
fig/{ztf,apphot,comspec}/                       all figures
doc/                                            summary.md · README.md (index) · {ztf,apphot,comspec}/ · literature/ · figures.pptx
_archive/                                       superseded prototype code (for the code reviews only)
```

`results/`, `fig/` and most of `data/` are gitignored (regenerable or large); the reference
tables and every README are tracked.

## Where the results are

| what | where |
|---|---|
| **production rates per (comet, phase)** | `results/comspec/gas_fit.csv` (`Q_H2O`, `Q_CO2`, `Q_CO` with errors, tiers, `n_eff`, `h2o_source`, caveats; the ZTF dust context as `afrho_*`) |
| continuum fits, skipped / not-fitted groups, apertures, phases | `results/comspec/continuum_summary.csv`, `skipped_groups.csv`, `not_fitted.csv`, `apertures.csv`, `phase_map.csv`, `phase_cuts.csv` |
| the run that produced them | `results/comspec/run.meta.json` (hash `a6463b4c610b`), `logs/`; open stand-ins `placeholders.csv` |
| study variants and the method matrix | `results/comspec/studies/`, `fig/comspec/studies/` |
| SPHEREx photometry (one row per exposure × aperture) | `results/apphot/photometry/<comet>.csv` + `.meta.json`; `results/apphot/status.csv`; stacks `results/apphot/stacks/` |
| ZTF photometry, coma profiles, Afρ trends | `results/ztf/photometry/<comet>.csv`, `results/ztf/profile/`, `results/ztf/afrho/` |
| figures | `fig/comspec/{emission_model,cont_subtract,phase_group,phase_images}/`, `fig/comspec/summary_*.png`; `fig/apphot/<comet>/`; `fig/ztf/{afrho,photometry,profile}/` |
| the review deck | `doc/figures.pptx` — three slides per fitted phase |

## Data

Raw data live on the external SSD and are found through each package's `directory.py`
(environment overrides in every case):

- SPHEREx cutouts and index: `/Volumes/T7/data/spherex-comet/` (`SPHEREX_T7_DIR`, `SPHEREX_DB`,
  `SPHEREX_FITS_DIR`); the in-tree sample `data/spherex_sample/` for tests.
- ZTF frames: `/Volumes/T7/data/ztf-comet/<comet>/` (`ZTFCOMET_DATA`; falls back to `data/ztf/`).
- Gaia DR3: `~/Desktop/data/gaia_dr3/` (`SPHEREX_GAIA_DIR`, `ZTFCOMET_GAIA`); build the
  declination cache once with `scripts/apphot/build_gaia_cache.py`.
- Reference tables (`data/reference/`) and the fluorescence database (`data/fluorescence/`)
  are in the tree.

## Documents

`doc/README.md` indexes everything.  Start with `doc/summary.md`, then
`doc/comspec/pipeline_decisions.md` (what the catalog pipeline does and why),
`doc/ztf/afrho_heliocentric_trends.md` (the dust side) and
`doc/apphot/pipeline_upgrade_notes.md` (the photometry).

## Git

One repository.  The `ztfcomet` history is the trunk; the `spherex-comspec` history was
merged in (browse it with `git log spherex-comspec-import -- spherex_comspec/<file>`);
`spherex_apphot` was added at the merge.  Remotes of the former repositories:
`github.com/bumhooLIM/ztfcomet`, `github.com/bumhooLIM/spherex-comspec`.
