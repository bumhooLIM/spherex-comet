# spherex-comet-catalog

Gas production rates of comets from SPHEREx spectrophotometry: H₂O, CO₂ and CO
for 68 comets, from the revised aperture photometry (`data/apphot/`) through one
pipeline, `spherex-comspec/`.

## Where the results are

| what | where |
|---|---|
| **production rates per (comet, phase)** | `results/gas_fit.csv` (`Q_H2O`, `Q_CO2`, `Q_CO` with errors, status, `h2o_source`; the ZTF dust context as `afrho_*`) |
| dust context per (comet, phase) | `results/afrho_ztf.csv`: ZTF r-band A(0°)fρ at each phase's ⟨r_h⟩ (ρ = 10 000 and 20 000 km) with how it was obtained, from the `ztf-comet` trends; attached to `gas_fit.csv` and `phase_map.csv` as `afrho_rh_au`, `afrho_10k_cm`, `afrho_10k_err_cm`, `afrho_10k_method`, `afrho_20k_*`, `afrho_note` |
| model curves and residuals per fit | `results/gas_fit_lines/` |
| continuum subtraction per band | `results/continuum_summary.csv`; point spectra in `data/emission/` |
| what was skipped and why | `results/skipped_groups.csv`, `results/not_fitted.csv` |
| phase (observing-state) groups | `results/phase_map.csv`, `results/phase_cuts.csv`, `data/phase_assignment.csv` |
| aperture per comet | `results/apertures.csv` |
| the run that produced them | `results/run.meta.json` (config hash), `results/logs/` |
| figures | `fig/emission_model/` (one per fit), `fig/cont_subtract/` (validation grid + raw spectrum per group), `fig/phase_group/`, `fig/summary_Q_vs_rhel.png`, `fig/summary_mixing_ratios.png` |
| stacked images per (comet, phase) | `fig/phase_images/<stem>.png` and `results/phase_stacks/<stem>_<band>.fits`: the ZTF r-band frames nearest the SPHEREx window and the SPHEREx exposures of the phase in the continuum (1.2–2.5 µm) and the H₂O / CO₂ / CO windows, comet-centred, north up (`scripts/phase_images.py`) |
| the review deck | `doc/figures.pptx` (`scripts/make_figure_slides.py`): three slides per fitted group — spectra and fits, the ZTF Afρ trend with the phase highlighted, the stacked images — with Q and the ZTF Afρ on every header |
| the studies that fixed the pipeline choices | `results/studies/`, `fig/studies/` (flag policy, `badphot`, distance correction, previous vs revised photometry) |
| placeholders still open | `results/placeholders.csv` |

## How to run

```bash
cd spherex-comspec
~/miniconda3/envs/spherex/bin/python main.py all          # group -> run -> analyze -> figures
~/miniconda3/envs/spherex/bin/python main.py run --variants dc_main --targets 24P 2P
~/miniconda3/envs/spherex/bin/python -m pytest -q tests
```

`main.py all` runs the main variant and its study partners; only the main
variant gets per-group figures (`--study-figs` for the rest).

The `afrho_*` columns come from the sibling `ztf-comet` project: its
`notebooks/afrho_trends.py` writes `results/afrho/spherex_afrho.csv` (and marks
the values in red on `fig/afrho/trend/<comet>.png` there), and

```bash
~/miniconda3/envs/spherex/bin/python scripts/attach_afrho_ztf.py   # [--ztf-root PATH]
```

pivots it to `results/afrho_ztf.csv` and attaches the columns to `gas_fit.csv`
and `phase_map.csv`.  The pipeline re-attaches them on every rerun while that
file exists and blanks rows whose phase geometry moved (`stale`); after a
regroup, rerun the ZTF script first, then this one.

The review deck and the stacked images:

```bash
~/miniconda3/envs/spherex/bin/python scripts/phase_images.py --workers 4   # fig/phase_images/, results/phase_stacks/
~/miniconda3/envs/spherex/bin/python scripts/make_figure_slides.py         # doc/figures.pptx; 24P:3 for one group
```

`phase_images.py` reads the SPHEREx cutouts through the `spherex-comet-apphot` package
(sibling checkout or `SPHEREX_APPHOT_ROOT`; the cutouts on its T7 root) and the ZTF
frames through `ztf-comet` (`ZTF_COMET_ROOT`); the trend slides come from
`ztf-comet/fig/afrho/trend_slide/`, drawn by its `afrho_trends.py`.

## Documents

`doc/README.md` is the index.  Read `doc/pipeline_decisions.md` first: what the
pipeline does, which alternatives were tried, why the current choice is the most
robust one, and which placeholders remain.  `spherex-comspec/README.md` is the
tool guide.  `handoff.md` is the working state for the next session.

## Layout

```
data/apphot/            revised aperture photometry (input; never written to)
data/reference/         orbit classes, literature production rates, the previous phase map
data/emission/          continuum summaries + point spectra (main variant)
results/, fig/          main products (above); studies/ one level down
notebooks/              rcparams.py, the photometry-comparison scripts, literature-table builders
spherex-comspec/        the pipeline package (git: github.com/bumhooLIM/spherex-comspec)
doc/                    method, decisions, comparison, literature
_archive/               everything superseded (previous photometry, notebooks, earlier runs); safe to delete
```
