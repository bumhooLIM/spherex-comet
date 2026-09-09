# `spherex-comspec` — gas production rates from SPHEREx comet spectrophotometry

An end-to-end model from the revised aperture photometry of the SPHEREx comet
catalog (`data/apphot/`, produced by `spherex_apphot`) to production
rates of H₂O, CO₂ and CO per comet and observing phase.  It reconstructs, as a
single package, what the project notebooks `phase_group_update.ipynb`,
`continuum_subtraction.ipynb` and `gas_emission_fit.ipynb` (with `emission-fitter/`)
did as three separate steps, and adds the two studies the revised photometry
makes possible: the effect of the **distance correction** on the continuum, and
of the **Gaia source flags** on contamination.

The scientific method is unchanged: the AKARI/IRC band-fluorescence technique
of Ootsubo et al. (2012) on a Haser coma with the Yamamoto (1981) aperture
filling factor, solved as one **joint linear least-squares** fit across the
2.7 / 4.3 / 4.7 µm bands.  [`doc/model_concept.md`](doc/model_concept.md) and
[`doc/fitting_methodology.md`](doc/fitting_methodology.md) describe the physics
and the formalism; this file describes how to run the tool and how to read what
it writes.  [`doc/spherex_comspec_reconstruction.md`](doc/spherex_comspec_reconstruction.md)
records what changed in the reconstruction, what the distance-correction and
source-flag studies found, and which placeholders to update first.

The package expects to live inside the `spherex-comet-catalog` project tree
(`data/apphot/` two levels up); set `COMSPEC_ROOT` / `COMSPEC_APPHOT_DIR`
to run it anywhere else.

---

## Layout

```
spherex-comspec/
├── main.py                    end-to-end driver: group | run | analyze | figures | all
├── tests/test_comspec.py      unit + reproduction tests (no pytest needed)
└── spherex_comspec/
    ├── config.py              every constant, window, threshold; the four dataclasses;
    │                          the PLACEHOLDERS registry; the DEFAULT_VARIANTS
    ├── directory.py           paths (environment-overridable)
    ├── dataio.py              photometry reading, aperture choice, flag policy, all file I/O
    ├── grouping.py            28-day epochs -> single-state phase groups (dynamic programme)
    ├── continuum.py           local polynomial continuum, validation, subtraction
    ├── gasmodel.py            Haser + Yamamoto + fluorescence  (Layers 0-4, unchanged)
    ├── instrument.py          SPHEREx channel bandpass          (Layer 5, unchanged)
    ├── fitting.py             design matrix, weighted solve, limits, coverage (Layer 6)
    ├── pipeline.py            orchestration: run_grouping, run_variant
    ├── analysis.py            flag-contamination and distance-correction studies
    └── plotting.py            every figure -- imported only by `main.py figures`
```

Nothing is installed; `main.py` and the tests put the package on `sys.path`.
Requires the `spherex` conda environment (numpy, pandas, astropy, scipy, matplotlib).

## Running

```bash
cd spherex-comspec
python main.py all                       # everything: group, six variants, analysis, figures
python main.py group                     # (1) regroup the epochs, once
python main.py run --variants dc_main    # (2) continuum + fit for one variant
python main.py analyze                   # (3) cross-variant tables
python main.py figures --variants dc_main # (4) figures for one variant
python main.py run --variants dc_main --targets 24P 2P   # a quick look
python tests/test_comspec.py
```

`group` writes `data/phase_assignment.csv` (one row per exposure) and
`results/phase_map.csv` (one row per group) and is shared by every
variant.  **The photometry tables are never modified.**

## The pipeline

### 1. Phase regrouping (`grouping.py`)

The `epoch` column of the photometry is the 28-day temporal grouping.  It has no
notion of where the comet was, so one epoch can span a large range of
heliocentric distance.  Each epoch is cut into `phase` groups under four rules in
priority order: (1) r_h range < 10 % inside a group; (2) inbound and outbound
never share a group where perihelion is resolved; (3) an emission band sampled
by epochs within tolerance is kept whole; (4) fewest groups, cuts placed in real
gaps.  Four targets are grouped by hand (`GroupingConfig.manual_edges`,
`manual_no_arc`).  The solver is a dynamic programme over
`(band breaks, groups, cut penalty, spread)` — a greedy sweep cannot build it.

Result on the revised photometry: **141 epochs → 174 phases**, 21 epochs
subdivided, 7 manual groups — the same count as before, and identical
group-by-group for 66 of 68 targets.  The two that differ (2025 K1, 2025 L1) do
so because rule 3 reads which exposures sampled a band with `badphot`, whose
meaning changed (below).

### 2. Aperture (`dataio.aperture_for`)

One aperture per target: 20 000 km when the mean r_h is inside 3 au, 40 000 km
beyond.  The revised photometry refuses an aperture below the PSF or beyond the
sky annulus, so the rule aperture can be missing for a distant target — 2014 UN271
at 14.6 au has no 40 000 km measurement at all.  When the rule aperture covers
fewer than 95 % of a target's exposures, the smallest larger `km` aperture that
does is used and logged: 2014 UN271 → 80 000 km; 2019 U5, 2022 R3, 2023 RS61
→ 60 000 km.  `results/apertures.csv` records the choice.

### 3. Flag policy (`config.FlagPolicy`, `dataio.select_spectrum`)

Which rows of that aperture enter a spectrum.  `badphot` rows (a bad pixel
inside the aperture) are dropped, optionally only above a `frac_badpix_ap`
threshold; source flags in `drop_flags` are dropped (`a` bright blend within
r_ap + 2 FWHM, `b` Gaia flux > 20 % of the comet's, `c` any source, `d`
SNR < 1).  The photometry's own flag precedence is `a > b > c > d > 0`.

### 4. Continuum subtraction (`continuum.py`)

Per band, a weighted polynomial in λ − λ_c fitted to the continuum window with
**all three** emission windows punched out, sigma-clipped, order capped at 1
when the sample does not bracket the band, with a 2-point fallback for sparse
continua.  Validated on bracketing and positivity (hard) and residual shape and
cross-validation (soft) → `PASS` / `WARN` / `FAIL`.  A group with fewer than 5
usable points or fewer than 2 inside the emission windows is skipped, not
`FAIL`ed.  Every band-row is written, including `FAIL`; filter on `verdict`.

### 5. Fitting (`fitting.py`)

`F_ν(λ_i) = Σ_X Q_X A_iX` solved by weighted linear least squares on the
`role == "emission"` channels of bands with an accepted verdict.  Each channel
keeps its own r_h and Δ and its own bandpass.  Coverage is decided on the
diagnostic `KEY_RANGES` (H₂O ≥ 3 channels in 2.60–2.80 µm, CO₂ ≥ 2 in
4.20–4.30, CO ≥ 2 in 4.60–4.70) — *not covered* is categorically different from
*not detected*.  `n_eff` (the participation ratio of the per-channel Fisher
information) says how many channels really carry a species; require ≥ 2.


**H₂O coverage rule.**  Q(H₂O) is anchored by the 2.7 µm main band (≥ 3 channels
in 2.60–2.80 µm).  Only when that band is *not covered* — the bright, close comets
such as 10P and 24P lose those channels to saturation and flags — do the 4.63 and
4.85 µm hot bands carry Q(H₂O) (`config.H2O_HOT_RANGE`: ≥ 3 channels in
4.55–4.90 µm, at least one beyond 4.75 µm so the value is separable from CO).  The
row then carries `h2o_source = "hot"`, the caveat says so, and the figures draw it
with an open marker; `FitConfig.h2o_hot_fallback=False` switches the rule off.
Hot-band g-factors are placeholders, so such values are provisional.

**Errors.**  The fit uses `source_sum_err_empirical_mjy` (`Variant.error_column`,
placeholder 10 applied): the formal error under-reports the annulus scatter.
## Distance correction — how Q stays physical

The revised photometry carries `flux_distcorr_mjy = F × r_h² × Δ²`, the flux the
comet would show at 1 au from both Sun and observer, with the factor stored per
row as `distcorr_factor`.  The `dc_*` variants fit the continuum on that column.

Q is a property of the comet, so it must not depend on which space the
continuum was fitted in.  The continuum-subtracted emission is therefore divided
back by **each channel's own** `distcorr_factor` before the design matrix is
built (`emis_raw_mjy`), and the solve runs in physical flux space.  There is no
separate "re-correction" step: that per-channel division *is* the
re-correction, and because the model is linear it is exact — the alternative
route (fit in corrected space with the design matrix scaled row-wise by the same
factor) gives identical Q to 3 × 10⁻¹⁶ on real data (`fitting.py`,
`tests/test_comspec.py`).  What the space *does* change is the continuum
polynomial, and that is what the distance-correction study measures.

Why fit the continuum in corrected space at all: SPHEREx scans a moving target
non-simultaneously, so the channels of one group were taken at different
geometries.  Inside one 28-day epoch of 24P the factor r_h²Δ² varies by 3×; a
polynomial through the raw spectrum would absorb part of that geometric gradient
into the band.  Even after regrouping to 10 % in r_h, the spread of r_h²Δ²
across the channels of a fitted group is 13 % (median), 55 % (90th percentile)
and up to 166 %, because Δ is not constrained by the grouping rules — and the
change in Q between the two continuum spaces correlates with it (r = 0.40 over
the groups detected in both spaces; 0.48 for the 2026-09-08 baseline).

## Variants (`config.DEFAULT_VARIANTS`)

| name | role | flag policy | continuum space | purpose |
|---|---|---|---|---|
| `dc_main` | main | drop rows with `frac_badpix_ap > 0.05`, drop flag `a` (`BASELINE_FLAGS`) | distance-corrected | **main result** (`MAIN_VARIANT`) |
| `dc_main_keep_a` | flags | baseline but flag `a` kept | distance-corrected | contamination study: flag `a` |
| `dc_main_no_b` | flags | baseline + drop flag `b` | distance-corrected | contamination study: flag `b` |
| `raw_main` | distcorr | baseline | physical | distance-correction study |
| `dc_main_strict` | badphot | baseline but any bad pixel drops the row | distance-corrected | badphot-policy study |
| `dc_all` | previous | strict `badphot`, every flag kept (the 2026-09-08 baseline) | distance-corrected | before/after the placeholder switch |

Each variant is a complete, independent run under `data/emission/`, `results/` and `fig/` for the main variant and under
`studies/<variant>/` for every other one.  A variant's `hash`
(written to `run.meta.json`) identifies its full parameter set.  `config.VARIANTS`
maps each name to its `Variant`, so a script drives one directly:
`run_variant(VARIANTS["dc_main"])`, `save_variant_figures(VARIANTS["dc_main_strict"], assignment)`.
The driver picks each study's partners by `Variant.role` (`main`, `flags`, `distcorr`,
`badphot`, `previous`), so a new variant only needs a row in `DEFAULT_VARIANTS`.

## Outputs

| path | content |
|---|---|
| `data/phase_assignment.csv` | `target, filename, obsid, jd_utc, epoch, arc, phase, manual` per exposure |
| `results/phase_map.csv` | one row per phase group: r_h / r_obs / timing / band-sampling statistics |
| `data/emission/<target>_<ap>km.csv` | per-band continuum summary (verdict, provenance, coefficients + covariance) |
| `data/emission/<target>_<ap>km_points.csv` | point-level spectrum: `flux`/`err` in the fit space, `emis_*` and `emis_raw_*`, `role`, flags |
| `results/gas_fit.csv` | **one row per (target, phase): Q, errors, limits, coverage, n_eff, mixing ratios, caveats** |
| `results/gas_fit_lines/` | dense model curves and per-channel residuals per fit |
| `results/continuum_summary.csv`, `skipped_groups.csv`, `not_fitted.csv`, `apertures.csv`, `run.meta.json` | provenance |
| `data/studies/<v>/`, `results/studies/<v>/`, `fig/studies/<v>/` | the same products for every study variant |
| `results/studies/flag_policy_*.csv`, `distcorr_effect_*.csv`; `fig/studies/*.png` | the cross-variant studies |
| `results/placeholders.csv` | the placeholder registry as a table |

`phase` in every file is the regrouped phase; `epoch` is its 28-day parent.
Pin `dtype={"target": str}` when reading: `2022E2` and `2024E1` are valid
scientific notation.

### Key columns of `gas_fit.csv`

| column | meaning |
|---|---|
| `Q_X`, `Q_X_err` | production rate [s⁻¹] and 1σ error rescaled by √χ²_ν when χ²_ν > 1 — **use this error** |
| `Q_X_fit` | the raw least-squares value, sign preserved (a negative one is a non-detection) |
| `Q_X_status` | `detected` / `upper_limit` / `negative_fit` / `not_covered` |
| `Q_X_upper_limit` | the k-σ limit where the status is not `detected` |
| `Q_X_n_eff` | effective channel count — **read before believing a small error** |
| `h2o_source`, `h2o_anchored` | `main`: Q(H₂O) rests on the 2.7 µm band; `hot`: on the 4.6–4.9 µm hot bands only (2.7 µm not covered; provisional, open markers in the figures); `none`: not covered |
| `CO2_H2O`, `CO_H2O` (+ `_err`) | mixing ratios with the covariance term |
| `fit_space`, `distcorr_factor_mean` | the space of the solve (always physical) and the group's mean factor |
| `n_flag_a`, `n_flag_b` | flagged channels that entered the fit |
| `caveats` | the interpretive warnings, joined by `;` |

## Placeholders

`config.PLACEHOLDERS` lists, in priority order, every value that is a stand-in
or an unattributed convention rather than a measurement — the SPHEREx LSF,
the band profiles, the expansion-velocity law, the 2.7 µm blue edge, the
`badphot` policy (applied on 2026-09-09 as `BASELINE_FLAGS`), the aperture rule,
the polynomial orders, the 1σ detection
tier, the negative-channel cut, the error column, the grouping thresholds, the
sufficiency gates, the opacity calibration, T_rot and the upstream flag
thresholds.  `results/placeholders.csv` is the same table; the
reconstruction document says which ones must be updated first.

## Tests

```bash
python tests/test_comspec.py
```

Filling-factor limits and linearity in Q; the grouping DP on synthetic
sequences; Q identical in physical and corrected space on a synthetic spectrum;
and, when the data are present, reproduction of the old `phase_update_map.csv`
for four targets, the aperture promotion of 2014 UN271, the nesting of the flag
policies, and the round trip of every saved continuum from its own columns.
