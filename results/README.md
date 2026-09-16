# `results/`

Main, clean products of each stage; gitignored (regenerable).  Studies and variants sit one
level down so they can never be mistaken for the catalog result.

## `comspec/` — the catalog (`scripts/comspec/main.py all`)

`phase` is the single-observing-state group of `data/comspec/phase_assignment.csv`.

| file | content |
|---|---|
| **`gas_fit.csv`** | **one row per (comet, phase): the production rates.**  Per species `Q_X`, `Q_X_err` (χ²-rescaled), `Q_X_err_formal`, `Q_X_nsig`, `Q_X_status` (`detected` ≥ 3σ, `marginal` 1–3σ, `upper_limit`, `negative_fit`, `not_covered`, `rejected` = a detection the 2026-09-15 review rejected as spurious: value withheld, `Q_X_fit` kept), `Q_X_upper_limit` (3σ), `Q_X_n_eff` (effective channels — read before believing a small error), `Q_X_n_key`; `h2o_source` (`main` = 2.7 µm band, `hot` = 4.6–4.9 µm hot bands), `h2o_anchored`; `CO2_H2O`, `CO_H2O` with errors; `chi2_red`, `n_points`, `bands_used`, geometry means incl. `v_hel_mean_kms` and `swings_CO`, `n_flag_a`, `n_flag_b`, `caveats`, `revision` (the case revision applied to the group, `spherex_comspec.revisions`); the ZTF dust context `afrho_rh_au`, `afrho_10k_cm`, `afrho_10k_err_cm`, `afrho_10k_method`, `afrho_20k_*`, `afrho_note` |
| `gas_fit_lines/` | model curves and residuals per fit (`<stem>_curve.csv`, `<stem>_points.csv`) |
| `continuum_summary.csv` | one row per (comet, phase, band): continuum fit, `verdict` PASS/WARN/FAIL with reasons, band emission flux, peak S/N, `revision` (the band's case directive of the 2026-09-15 memo; a waived check says so in `notes`) |
| `skipped_groups.csv`, `not_fitted.csv` | groups that could not be analysed or fitted, with the reason |
| `apertures.csv` | the aperture per (comet, phase): `r_ap_km`, `reason` (`near` 20 000 km inside 3 au, `far` 40 000 km, `far->enlarged` 60 000 km under 1.5 px), coverage, median geometry |
| `phase_map.csv`, `phase_cuts.csv` | the phase groups (r_h range, arc *index*, epochs; the `afrho_*` columns) and every cut with its cost |
| `afrho_ztf.csv` | the ZTF A(0°)fρ per (comet, phase) at ρ = 10 000 and 20 000 km with method, range, grade and note (`scripts/comspec/attach_afrho_ztf.py`) |
| `placeholders.csv` | the stand-in values still to be replaced, in priority order |
| `run.meta.json`, `logs/` | the configuration hash of the run (`a6463b4c610b`) and its logs |
| `phase_stacks/` | `<stem>_<band>.fits` comet-centred north-up median stacks (ZTF r, continuum, H₂O, CO₂, CO) per fitted phase (`scripts/comspec/phase_images.py`) |
| `studies/` | the same tables per study variant (`<variant>/`), the cross-variant tables (`flag_policy_*.csv`, `distcorr_effect_*.csv`), `apphot_comparison/` (previous vs revised photometry), `annulus_previous/` (the pipeline on the old-annulus photometry), `method_matrix/matrix.csv` (the 54-run method matrix of 2026-09-14) |

## `apphot/` — SPHEREx aperture photometry (`scripts/apphot/main.py`)

| file | content |
|---|---|
| **`photometry/<slug>.csv`** | one row per (exposure × aperture), 68 comets, 455 245 rows; **the input of the catalog stage**.  Identity (`filename`, `objdesig`, `obsid`, `detector`, `epoch`, `date_obs`, `jd_utc`, `jd_tdb`), channel (`wl`, `wlwidth`, `sun_jy`), aperture (`ap_label`, `r_ap_pix/km/arcsec`, `r_in_pix`, `r_out_pix`, `psf_fwhm_pix`, areas), bad pixels (`n_badpix_ap`, `frac_badpix_ap`, `badphot`), sky (`sky_median_mjy_per_pix`, `sky_std…`, `sky_excess_ratio`), flux (`source_sum_mjy` — negatives kept, `err_pix_mjy`, `err_skylevel_mjy`, `err_skyscatter_mjy`, `source_sum_err_mjy`, **`source_sum_err_empirical_mjy`** — the column the fits use, `snr`, `abmag`), distance correction (`flux_distcorr_mjy`, `distcorr_factor`), Gaia (`n_gaia`, `gmag_brightest`, `gmag_eff`, `gmag_nearest`, `dist_gmag_nearest`, `sourceflag` a/b/c/d/0), geometry (`r_hel`, `r_obs`, `alpha`, ecliptic coordinates and rates, `vmag`, `xcen`, `ycen`, …).  Read with `dtype={"sourceflag": str}` |
| `photometry/<slug>.meta.json` | the full `Config`, its hash (`9fcc7ea3871a`), counters and the epoch table that produced the CSV |
| `status.csv` | one row per target: status, counts, timing, config hash, output path — the resume table |
| `stacks/` | `<slug>_epoch<N>_<band>.fits` per-epoch band stacks (`--stack`) |
| `logs/` | the run logs (`apphot_<target>_<utc>.log`) |
| `contamination_by_target.csv`, `flag_effectiveness_survey.csv` | the growth-curve contamination metric and the source-flag recall/lift per target (`spherex_apphot.diagnostics`) |

## `ztf/` — ZTF dust photometry (`scripts/ztf/`)

| folder | content |
|---|---|
| `photometry/<slug>.csv` | one row per (frame × aperture) with calibrated magnitudes, Afρ, A(0°)fρ, the `flags` column and one boolean per check, `quality_ok`, `rho_km` — 56 comets |
| `profile/` | `<slug>_profile.csv` (comet radial profile per frame), `<slug>_stars.csv` (the field-star stacks), `<slug>_summary.csv` (slopes per frame) |
| `afrho/` | `trends.csv`, `peaks.csv`, `breaks.csv`, `outbursts.csv`, `colour.csv`, `elements.csv` (q, e, T_p), `spherex_windows.csv` (JD / r_h range of every SPHEREx phase), **`spherex_afrho.csv`** (the A(0°)fρ at every SPHEREx phase and aperture with method and note — what `attach_afrho_ztf.py` reads), `systematics_2024E1.csv` |
