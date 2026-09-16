# `notebooks/` — validation notebooks, study scripts and the figure style

`rcparams.py` is the project's matplotlib style (base font 20 pt, bold; see `CLAUDE.md`).
`spherex_apphot.plotting` and `spherex_comspec.plotting` import it for its side effects;
a notebook does `import rcparams` after putting this directory on `sys.path`.  Do not
restate the settings, and do not override them downward — enlarge `figsize` instead.
(`ztfcomet.rcparams` is the older 15-pt style of the ZTF figures.)

Every notebook walks up from its working directory to the project root, so it runs from
any kernel start directory.  Outputs saved in a notebook are kept only where they document
a validated run; strip them before committing (`nbstripout`).

## `ztf/`
| notebook | validates |
|---|---|
| `query.ipynb` | the query stage on 24P: path resolution, Horizons ephemeris, IRSA search, footprint check, download |
| `afrho.ipynb` | the photometry stage on 24P, and the four corrections that separate `ztfcomet.phot` from the prototype |
| `figure.ipynb` | the core figures (cutouts, Afρ) — still writes the old per-target figure names (open item) |

## `apphot/`
| notebook | purpose |
|---|---|
| `main.ipynb` | end-to-end walkthrough of the SPHEREx photometry on `data/spherex_sample/` (or the SSD when mounted), calling the same `pipeline.run_target` as `scripts/apphot/main.py`; every figure type |
| `build_db_filtered.ipynb` | data preparation: `db.parq` → `db_filtered.parq` (V < 22, r_h < 20 au, WCS attached); run once per data release |

## `comspec/`
| notebook | purpose |
|---|---|
| `figures.ipynb` | the catalog figures from `results/comspec/gas_fit.csv`: Q vs ⟨r_h⟩ per species with the 3σ limits, the comets with a value at ≥ 2 phases (tracks labelled by designation), and Q / A(0°)fρ vs ⟨r_h⟩ at ρ = 10 000 and 20 000 km → `fig/comspec/Q_*.png` |

Standalone study scripts of the catalog stage (run from the project root):
`method_matrix.py` (the 54-run continuum/fit method matrix in isolated result trees),
`apphot_comparison.py` / `_figures.py` / `_summary.py` (previous vs revised photometry —
the previous set was deleted on 2026-09-15; the products stay in
`results/comspec/studies/apphot_comparison/`), `build_gicquel2023_table.py` and
`build_harrington_pinto2022_table.py` (the literature Q tables in `data/reference/`), and
`fluorescence_gfm/` (the General Fluorescence Model reconstruction of the GSFC g-factor
database: `build_fluorescence_db.py` → `data/fluorescence/`, `validate_gfm.py`,
`plot_fluorescence_profiles.py`; physics in `doc/comspec/fluorescence_database.md`).
