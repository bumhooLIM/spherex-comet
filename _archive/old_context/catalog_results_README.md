# `results/`

Products of `spherex-comspec` (`python main.py all`).  `phase` is the
single-observing-state group of `data/phase_assignment.csv`.

| file | content |
|---|---|
| `gas_fit.csv` | **one row per (comet, phase): the production rates.**  Per species `Q_X`, `Q_X_err` (χ²-rescaled), `Q_X_err_formal`, `Q_X_nsig`, `Q_X_status` (`detected` ≥ 3σ, `marginal` 1–3σ with the value reported, `upper_limit` < 1σ, `negative_fit`, `not_covered`), `Q_X_upper_limit` (3σ), `Q_X_n_eff` (effective channels — read before believing a small error), `Q_X_n_key`; `h2o_source` (`main` = 2.7 µm band, `hot` = 4.6–4.9 µm hot bands only, used when 2.7 µm is not covered), `h2o_anchored`; `CO2_H2O`, `CO_H2O` mixing ratios with errors; `chi2_red`, `n_points`, `bands_used`, geometry means including `v_hel_mean_kms` (heliocentric radial velocity, positive receding) and `swings_CO` (the g(CO) velocity factor applied), `n_flag_a`, `n_flag_b`, `caveats` |
| `gas_fit_lines/` | model curves and residuals per fit |
| `continuum_summary.csv` | one row per (comet, phase, band): continuum fit, validation (`verdict` PASS/WARN/FAIL with reasons), band emission flux and peak S/N |
| `skipped_groups.csv`, `not_fitted.csv` | groups that could not be analysed or fitted, with the reason |
| `apertures.csv` | the aperture analysed per (comet, phase): `r_ap_km`, `reason` (`near` = 20 000 km inside 3 au, `far` = 40 000 km beyond, `far->enlarged` = 60 000 km because the rule radius was under 1.5 px), `coverage` of the phase's exposures, the phase's median `r_hel_med`, `r_obs_med`, `pixel_scale_km` and `r_ap_pix` (rule of 2026-09-14; `studies/dc_ap_snr/apertures.csv` keeps the S/N-rule evidence columns of 2026-09-12) |
| `phase_map.csv`, `phase_cuts.csv` | the phase groups (r_h range, arc, epochs) and every cut with its cost |
| `placeholders.csv` | the stand-in values still to be replaced, in priority order |
| `run.meta.json`, `logs/` | configuration hash and run logs |
| `studies/` | the same tables for each study variant (`studies/<variant>/`), the cross-variant tables (`flag_policy_*.csv`, `distcorr_effect_*.csv`), the previous-vs-revised photometry comparison (`apphot_comparison/`), the same pipeline on the archived photometry (`annulus_previous/`) and the 2026-09-14 continuum/fit method matrix (`method_matrix/matrix.csv`, one row per run with the robust, clean and negative-tail counts; `notebooks/method_matrix.py`) |
