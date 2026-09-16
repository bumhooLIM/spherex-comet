# Case revisions of the review memo (2026-09-15)

_Applied 2026-09-16 · memo `doc/notes_ver260915.xlsx` (38 rows) · code `spherex_comspec/revisions.py`,
`spherex_comspec.config.GroupingConfig`, `ztfcomet.config.AFRHO_EPOCH_OVERRIDES` · main run
`5741d6611f60` (`spherex_comspec` 1.3.0) against the 2026-09-14 run `a6463b4c610b` · report generator
`notebooks/comspec/case_revision_report.py` · deck `doc/figures.pptx` carries the previous slides of
every revised group directly before the revised ones._

The memo went through the per-phase figures and asked, case by case, for a different continuum
window, a fixed polynomial order, the exclusion of star-contaminated points, a waiver of the
negative-continuum or one-sided-continuum rejection, a fit on fewer channels than the coverage
rule demands, the rejection of a spurious detection, one regrouping (240P) and three Afρ
adjustments.  None of these is expressible as a global rule, so they are stored as data — one
`CaseRevision` per (target, phase) with a `BandRevision` per band — and applied wherever a group
is processed (batch run, every study variant, the figures), so nothing can drift apart.  The
variant `dc_main_norev` is the same configuration without them: it reproduces the 2026-09-14
census exactly, so every difference below is the memo's doing.

## 1. How each kind of directive was implemented

| memo wording | implementation | where it acts |
|---|---|---|
| "Move the left/right window to X µm" | `cont_lo` / `cont_hi` (continuum) or `em_lo` / `em_hi` (emission) of that band; the punch-out every other band applies follows the revised emission window | `continuum.fit_continuum` through `revisions.effective_windows` |
| "Increase / Reduce the order a → b" | `order = b`, fixed: no cross-validation for that band; a fixed order is fitted on as few as order + 2 points instead of falling back to two points below six | `fit_continuum` |
| "Exclude top N data points in the left/right window" | the N brightest points of that side of the continuum window are ranked over every point the validation figure shows there (those inside the neighbouring CO window included, since the reviewer counted them) and the continuum candidates among them are dropped before the sigma clipping | `fit_continuum` (`exclude_top_blue` / `exclude_top_red`) |
| "Ignore continuum to be negatived but fit it, and make PASS" | the positivity check is waived (`ignore_negative`); the fit and the subtraction are unchanged; the shape and CV soft checks still decide between PASS and WARN, both of which enter the fit | `continuum.validate_fit` |
| "Make 1D continuum fit only with left window from 2.0 to 2.5 µm" / "only with right window and make PASS" | `cont` window as given, `order = 1`, `one_sided_ok`: no extension to the empty side and the bracketing check is waived; the saved point neighbourhood now spans the emission window as well (it used to span the continuum window only, which lost the band's channels) | `fit_continuum`, `validate_fit`, `process_group` |
| "Exclude flag b points" | rows with that source flag leave the band (continuum and emission) before anything is fitted | `process_group` (`drop_flags`) |
| "Exclude the top N data points before emission fitting" | the N brightest emission channels of the band get `role = "excluded"`: out of the production-rate fit and of the band flux | `process_group` (`exclude_top_emission`) |
| "Make emission fit with N channels" / "single channel + nearby 2 continuum points as a baseline" | coverage of the species on ≥ N emission channels of its band instead of the `KEY_RANGES` rule (H₂O ≥ 3 in 2.55–2.80 µm, CO₂ ≥ 2 in 4.20–4.30, CO ≥ 2 in 4.60–4.70); the two-point continuum through the nearest points is what the pipeline already does below six continuum points | `fitting.fit_production_rates` (`min_channels`) |
| "This is fault detection. Reject the detection." | `Q_X_status = "rejected"`: the fitted value stays in `Q_X_fit` / `Q_X_nsig`, `Q_X` and the limit are blank, the species is outside every census and mixing ratio | `fit_production_rates` (`reject`) |
| "Regroup phase 2-3" (240P) | a manual r_h bin edge at 2.38 au (`manual_edges`) and the target exempted from rule 1b (`delta_tol_exempt`): the former phases 2 (94 exposures) and 3 (6) are one phase whose Δ spread is 26 % | `grouping.subdivide` |
| "Extrapolate the Afrho trend to S1" (47P) | the grade-D rising law may be extrapolated for that epoch (`extrap_grades` includes D; S1 lies 0.005 dex beyond its 2.82–2.90 au range) | `activity.afrho_at_epoch` via `afrho_trends.py` |
| "Refit the Afrho trend to S1–S4" (210P) | ZTF observed 210P after perihelion only; the four inbound SPHEREx phases are read off the outbound law (`leg="outbound"`), S1 as a 0.05 dex extrapolation of its 0.61–2.19 au range — a symmetric-activity assumption that the note states | `afrho_at_epoch` |
| "Overestimated. Refit the Afrho to the nearby point." (217P S1) | the two of the five in-window frames that sit at the start of the 2.30 au outburst are excluded from the direct mean (`exclude_outburst`) | `afrho_at_epoch` |

Every directive is repeated in the products: `continuum_summary.csv` and `gas_fit.csv` carry a
`revision` column and the text in `notes` / `caveats`, the validation figures print it under
the band, the ZTF `note` names the memo.

## 2. Case by case

Q entries: value ± error (σ, n_eff), `marginal` for 1–3σ, `< limit (σ)` for limits and negative
fits, `rejected (σ)` for a withheld detection; Afρ at 10 000 km.  Only the species that changed
are listed; `unchanged` means the memo row was applied and nothing moved.  Band verdicts give the
verdict, the polynomial order, the number of emission channels and the continuum window.

| comet | phase | directive (memo) | band verdicts before → after | Q before → after |
|---|---|---|---|---|
| 2P | 1 | 2.7um: continuum window 2.40-.. um; 4.7um: continuum window ..-4.95 um | 2.7um: FAIL (o2, 4 em, 2.30–3.00) → WARN (o1, 4 em, 2.40–3.00); 4.7um: PASS (o2, 11 em, 4.40–5.00) → PASS (o1, 11 em, 4.40–4.95) | H₂O: not covered → < 3.7e+26 (+0.6σ); CO₂: 9.31e+25 ± 6.5e+24 (14.3σ, n_eff 1.3) → 9.31e+25 ± 6.7e+24 (13.9σ, n_eff 1.3); CO: < 6.9e+26 (-1.6σ) → < 4.2e+26 (-2.7σ) |
| 2P | 2 | 4.3um: negative continuum accepted | 4.3um: FAIL (o1, 5 em, 3.90–4.65) → PASS (o1, 5 em, 3.90–4.65) | H₂O: 1.80e+26 ± 6.0e+25 (3.0σ, n_eff 1.8) marginal → 1.79e+26 ± 9.8e+25 (1.8σ, n_eff 1.8) marginal; CO₂: not covered → 7.83e+25 ± 8.1e+24 (9.6σ, n_eff 1.5); CO: < 3.1e+26 (-0.1σ) → < 5.0e+26 (-0.1σ) |
| 10P | 1 | 2.7um: flag b rows dropped | 2.7um: PASS (o1, 6 em, 2.30–3.00) → WARN (o1, 6 em, 2.30–3.00) | H₂O: < 2.9e+26 (-0.1σ) → < 2.4e+26 (-0.8σ); CO₂: 2.56e+25 ± 5.5e+24 (4.7σ, n_eff 2.0) → 2.56e+25 ± 6.1e+24 (4.2σ, n_eff 2.0) |
| 24P | 7 | 4.7um: order 2 fixed | 4.7um: PASS (o1, 29 em, 4.40–5.00) → PASS (o2, 29 em, 4.40–5.00) | H₂O: 3.14e+28 ± 1.2e+28 (2.6σ, n_eff 16.1) marginal → < 6.1e+28 (+0.7σ); CO: < 4.9e+26 (-2.0σ) → < 6.2e+26 (-2.6σ) |
| 24P | 8 | 4.7um: order 2 fixed | 4.7um: PASS (o1, 46 em, 4.40–5.00) → PASS (o2, 46 em, 4.40–5.00) | H₂O: 2.12e+28 ± 1.3e+28 (1.6σ, n_eff 32.4) marginal → 3.30e+28 ± 2.3e+28 (1.5σ, n_eff 34.4) marginal; CO: < 5.5e+26 (-3.8σ) → < 7.4e+26 (-2.3σ) |
| 43P | 1 | 4.7um: order 1 fixed | 4.7um: PASS (o3, 19 em, 4.40–5.00) → WARN (o1, 19 em, 4.40–5.00) | H₂O: 2.37e+26 ± 2.4e+25 (9.7σ, n_eff 3.6) → 2.37e+26 ± 2.4e+25 (10.1σ, n_eff 3.6); CO₂: 5.62e+25 ± 3.2e+24 (17.8σ, n_eff 1.3) → 5.62e+25 ± 3.1e+24 (18.3σ, n_eff 1.3); CO: < 2.5e+26 (+0.4σ) → < 1.8e+26 (+0.2σ) |
| 47P | 1 | ztf_afrho: Extrapolate the Afrho trend to S1 | — | Afρ(10k): none (rising law is unconstrained (grade D) and <r_h> = 2.93 au is) → 119 ± 8 cm (trend_extrap, rising) |
| 47P | 2 | 2.7um: continuum window ..-2.95 um; emission window ..-2.85 um | 2.7um: PASS (o1, 4 em, 2.30–3.00) → PASS (o1, 5 em, 2.30–2.95) | H₂O: 8.32e+25 ± 3.8e+25 (2.2σ, n_eff 2.4) marginal → 9.01e+25 ± 4.8e+25 (1.9σ, n_eff 2.5) marginal; CO₂: 2.31e+25 ± 3.2e+24 (7.3σ, n_eff 1.8) → 2.31e+25 ± 3.1e+24 (7.5σ, n_eff 1.8); CO: < 2.3e+26 (-0.1σ) → < 2.2e+26 (-0.1σ) |
| 63P | 2 | 2.7um: order 1 fixed | 2.7um: PASS (o2, 5 em, 2.30–3.00) → PASS (o1, 5 em, 2.30–3.00) | H₂O: 2.02e+27 ± 9.4e+26 (2.1σ, n_eff 4.5) marginal → 1.99e+27 ± 3.7e+26 (5.4σ, n_eff 4.3); CO₂: 1.63e+26 ± 1.1e+25 (14.3σ, n_eff 5.3) → 1.63e+26 ± 1.2e+25 (14.1σ, n_eff 5.3) |
| 124P | 1 | 4.3um: 3 brightest red continuum point(s) excluded; 4.7um: 3 brightest blue continuum point(s) excluded | 4.3um: PASS (o1, 4 em, 3.90–4.65) → PASS (o1, 4 em, 3.90–4.65); 4.7um: PASS (o1, 16 em, 4.40–5.00) → PASS (o1, 16 em, 4.40–5.00) | H₂O: < 7.6e+27 (-0.8σ) → < 6.9e+27 (-0.6σ); CO₂: 5.55e+24 ± 2.0e+24 (2.8σ, n_eff 1.4) marginal → 6.60e+24 ± 1.9e+24 (3.4σ, n_eff 1.4); CO: < 1.4e+26 (+0.2σ) → < 1.5e+26 (+0.5σ) |
| 124P | 3 | 4.7um: continuum window ..-4.95 um; order 1 fixed | 4.7um: PASS (o3, 23 em, 4.40–5.00) → PASS (o1, 23 em, 4.40–4.95) | H₂O: 1.47e+26 ± 3.4e+25 (4.4σ, n_eff 1.8) → 1.47e+26 ± 3.5e+25 (4.2σ, n_eff 1.8); CO₂: 3.85e+24 ± 2.0e+24 (2.0σ, n_eff 1.2) marginal → 3.85e+24 ± 2.1e+24 (1.9σ, n_eff 1.2) marginal; CO: < 1.1e+26 (-1.4σ) → < 9.7e+25 (-1.3σ) |
| 210P | 2 | ztf_afrho: Refit the Afrho trend to S1, S2, S3, and S4. | — | Afρ(10k): none (ZTF fits only the outbound phase; the SPHEREx epoch is inbou) → 47 ± 6 cm (trend, outbound) |
| 217P | 1 | ztf_afrho: Overestimated. Refit the Afrho to the nearby point. | — | Afρ(10k): 124 ± 25 cm (direct, outbound) → 54 ± 4 cm (direct, outbound) |
| 217P | 3 | 4.3um: negative continuum accepted | 4.3um: FAIL (o3, 3 em, 3.90–4.65) → PASS (o3, 3 em, 3.90–4.65) | CO₂: not covered → 5.66e+25 ± 8.8e+24 (6.4σ, n_eff 1.6) |
| 229P | 1 | 4.3um: negative continuum accepted | 4.3um: FAIL (o1, 6 em, 3.90–4.65) → PASS (o1, 6 em, 3.90–4.65) | H₂O: < 1.1e+27 (-1.5σ) → < 9.6e+26 (-1.7σ); CO₂: not covered → < 3.6e+25 (+0.8σ) |
| 240P | 2 | phase_note: Regroup phase 2-3 (GroupingConfig.manual_edges / delta_tol_exempt) | — | H₂O: < 2.1e+29 (+0.1σ) / no fit → < 2.0e+29 (+0.1σ); CO₂: not covered / no fit → not covered; CO: < 2.0e+27 (-0.0σ) / no fit → < 1.9e+27 (-0.0σ) |
| 306P | 4 | 2.7um: order 1 fixed; 4.3um: 3 brightest emission channel(s) excluded from the fit | 2.7um: PASS (o2, 10 em, 2.30–3.00) → PASS (o1, 10 em, 2.30–3.00); 4.3um: PASS (o1, 11 em, 3.90–4.65) → PASS (o1, 8 em, 3.90–4.65) | H₂O: 6.84e+25 ± 7.0e+24 (9.8σ, n_eff 7.0) → 6.27e+25 ± 4.7e+24 (13.3σ, n_eff 7.0); CO₂: < 9.7e+23 (-0.1σ) → 3.16e+23 ± 2.5e+23 (1.3σ, n_eff 4.4) marginal; CO: < 1.9e+25 (+0.4σ) → < 1.5e+25 (+0.6σ) |
| 2019U5 | 1 | 4.3um: negative continuum accepted | 4.3um: FAIL (o3, 8 em, 3.90–4.65) → PASS (o3, 8 em, 3.90–4.65) | H₂O: 7.15e+26 ± 3.7e+26 (1.9σ, n_eff 2.3) marginal → 7.14e+26 ± 3.2e+26 (2.2σ, n_eff 2.3) marginal; CO₂: not covered → 4.28e+25 ± 3.0e+25 (1.4σ, n_eff 1.4) marginal |
| 2019U5 | 2 | 2.7um: 2 brightest blue continuum point(s) excluded; 2 brightest red continuum point(s) excluded; negative continuum accepted | 2.7um: FAIL (o2, 6 em, 2.30–3.00) → WARN (o1, 6 em, 2.30–3.00) | H₂O: not covered → 7.99e+26 ± 3.2e+26 (2.5σ, n_eff 3.1) marginal |
| 2022N2 | 2 | 4.3um: emission window 4.13-.. um; fitted with >= 2 emission channel(s) | 4.3um: PASS (o1, 4 em, 3.90–4.65) → PASS (o1, 6 em, 3.90–4.65) | H₂O: < 4.6e+27 (+0.4σ) → 4.78e+26 ± 1.4e+26 (3.4σ, n_eff 1.4); CO₂: not covered → 3.30e+26 ± 8.7e+24 (38.2σ, n_eff 1.0); CO: < 4.4e+27 (-0.0σ) → < 4.6e+26 (-0.2σ) |
| 2022QE78 | 1 | 2.7um: fitted with >= 2 emission channel(s) | 2.7um: PASS (o1, 2 em, 2.30–3.00) → PASS (o1, 2 em, 2.30–3.00) | H₂O: not covered → 1.97e+27 ± 2.6e+26 (7.7σ, n_eff 1.3); CO₂: 1.53e+26 ± 7.2e+25 (2.1σ, n_eff 1.4) marginal → 1.53e+26 ± 2.8e+25 (5.5σ, n_eff 1.4); CO: < 4.7e+27 (+0.1σ) → < 1.8e+27 (+0.2σ) |
| 2022R6 | 2 | 4.3um: 4 brightest blue continuum point(s) excluded | 4.3um: PASS (o1, 2 em, 3.90–4.65) → PASS (o1, 2 em, 3.90–4.65) | CO₂: 7.42e+25 ± 6.3e+25 (1.2σ, n_eff 1.0) marginal → 8.56e+25 ± 3.6e+25 (2.4σ, n_eff 1.0) marginal |
| 2023A3 | 2 | 4.3um: 1 brightest red continuum point(s) excluded | 4.3um: PASS (o1, 4 em, 3.90–4.65) → PASS (o1, 4 em, 3.90–4.65) | CO₂: 3.77e+25 ± 1.2e+25 (3.1σ, n_eff 1.9) → 4.29e+25 ± 1.3e+25 (3.3σ, n_eff 1.9); CO: < 6.3e+26 (-0.7σ) → < 6.8e+26 (-0.7σ) |
| 2023R1 | 1 | 2.7um: continuum window 2.40-.. um | 2.7um: PASS (o2, 5 em, 2.30–3.00) → WARN (o1, 5 em, 2.40–3.00) | H₂O: 6.34e+26 ± 2.3e+26 (2.8σ, n_eff 2.3) marginal → 5.07e+26 ± 1.9e+26 (2.7σ, n_eff 2.2) marginal; CO₂: 1.60e+26 ± 1.1e+25 (14.0σ, n_eff 2.5) → 1.60e+26 ± 1.1e+25 (14.4σ, n_eff 2.5); CO: < 9.2e+26 (-1.2σ) → < 8.9e+26 (-1.3σ) |
| 2023R1 | 4 | 2.7um: continuum window 2.00-2.50 um; order 1 fixed; one-sided continuum accepted, no extension; fitted with >= 3 emission channel(s) | 2.7um: FAIL (o1, 3 em, 2.30–3.80) → PASS (o1, 3 em, 2.00–2.50) | H₂O: not covered → 3.97e+27 ± 1.7e+27 (2.4σ, n_eff 1.7) marginal; CO: 1.04e+27 ± 7.4e+26 (1.4σ, n_eff 3.4) marginal → 1.01e+27 ± 7.1e+26 (1.4σ, n_eff 3.4) marginal |
| 2023RS61 | 1 | rejected: H2O | — | H₂O: 1.14e+29 ± 2.4e+28 (4.7σ, n_eff 1.9) → rejected (+4.7σ) |
| 2023U1 | 2 | 4.7um: emission window 4.50-.. um | 4.7um: WARN (o1, 12 em, 4.40–5.00) → WARN (o1, 13 em, 4.40–5.00) | H₂O: < 2.7e+27 (+0.6σ) → < 2.8e+27 (+0.6σ); CO: < 4.3e+27 (-0.8σ) → < 4.4e+27 (-0.4σ) |
| 2023V1 | 2 | 4.3um: 3 brightest blue continuum point(s) excluded | 4.3um: PASS (o3, 4 em, 3.90–4.65) → PASS (o1, 4 em, 3.90–4.65) | H₂O: < 3.6e+26 (+0.4σ) → < 4.0e+26 (+0.4σ); CO₂: 5.54e+25 ± 1.4e+25 (3.9σ, n_eff 2.0) → 5.22e+25 ± 9.8e+24 (5.3σ, n_eff 1.9); CO: < 8.4e+26 (+0.8σ) → < 9.0e+26 (+0.8σ) |
| 2024A1 | 1 | 4.3um: 5 brightest blue continuum point(s) excluded | 4.3um: PASS (o1, 5 em, 3.90–4.65) → PASS (o1, 5 em, 3.90–4.65) | H₂O: < 3.4e+26 (-0.9σ) → < 5.7e+26 (-0.5σ); CO₂: 4.55e+25 ± 1.0e+25 (4.4σ, n_eff 2.0) → 6.85e+25 ± 1.5e+25 (4.5σ, n_eff 1.9); CO: < 7.9e+26 (-0.2σ) → < 1.3e+27 (-0.1σ) |
| 2024E1 | 2 | 4.3um: order 1 fixed; negative continuum accepted | 4.3um: FAIL (o2, 10 em, 3.90–4.65) → PASS (o1, 10 em, 3.90–4.65) | H₂O: 3.09e+27 ± 2.4e+26 (12.8σ, n_eff 2.0) → 3.09e+27 ± 4.2e+26 (7.4σ, n_eff 2.0); CO₂: not covered → 1.61e+26 ± 2.4e+25 (6.8σ, n_eff 2.7); CO: < 1.6e+27 (-0.4σ) → < 2.8e+27 (-0.2σ) |
| 2024J2 | 1 | 4.3um: fitted with >= 3 emission channel(s) | 4.3um: PASS (o1, 3 em, 3.90–4.65) → PASS (o1, 3 em, 3.90–4.65) | H₂O: < 1.6e+27 (+0.2σ) → < 1.5e+27 (+0.2σ); CO₂: not covered → 3.18e+26 ± 2.8e+26 (1.1σ, n_eff 1.1) marginal; CO: 6.01e+27 ± 8.3e+26 (7.2σ, n_eff 5.4) → 6.01e+27 ± 8.2e+26 (7.3σ, n_eff 5.4) |
| 2024L5 | 1 | 4.3um: fitted with >= 5 emission channel(s) | 4.3um: WARN (o1, 5 em, 3.90–4.65) → WARN (o1, 5 em, 3.90–4.65) | H₂O: < 3.4e+26 (+0.3σ) → < 1.6e+26 (+0.6σ); CO₂: not covered → 3.94e+25 ± 4.2e+24 (9.4σ, n_eff 1.1); CO: < 6.4e+26 (-0.1σ) → < 2.7e+26 (-0.4σ) |
| 2024L5 | 4 | 4.3um: negative continuum accepted | 4.3um: FAIL (o1, 4 em, 3.90–4.65) → PASS (o1, 4 em, 3.90–4.65) | H₂O: 1.94e+26 ± 1.3e+26 (1.5σ, n_eff 1.9) marginal → 1.94e+26 ± 1.2e+26 (1.6σ, n_eff 1.9) marginal; CO₂: not covered → 1.16e+25 ± 9.4e+24 (1.2σ, n_eff 1.5) marginal; CO: < 5.6e+26 (-1.4σ) → < 5.0e+26 (-1.5σ) |
| 2025L1 | 1 | 4.3um: fitted with >= 1 emission channel(s) | 4.3um: WARN (o1, 1 em, 3.90–4.65) → WARN (o1, 1 em, 3.90–4.65) | H₂O: 2.15e+26 ± 9.9e+25 (2.2σ, n_eff 4.8) marginal → 2.10e+26 ± 8.6e+25 (2.4σ, n_eff 4.8) marginal; CO₂: not covered → 4.30e+26 ± 2.2e+26 (2.0σ, n_eff 1.0) marginal |
| 2025R1 | 1 | 4.3um: fitted with >= 2 emission channel(s) | 4.3um: PASS (o1, 2 em, 3.90–4.65) → PASS (o1, 2 em, 3.90–4.65) | H₂O: 2.56e+26 ± 5.3e+25 (4.9σ, n_eff 1.9) → 2.56e+26 ± 3.4e+25 (7.6σ, n_eff 1.9); CO₂: not covered → 1.53e+25 ± 3.4e+24 (4.6σ, n_eff 1.0); CO: < 1.1e+26 (-1.2σ) → < 7.4e+25 (-1.8σ) |
| 2025R1 | 4 | 4.3um: 4 brightest red continuum point(s) excluded; rejected: H2O | 4.3um: PASS (o1, 5 em, 3.90–4.65) → PASS (o1, 5 em, 3.90–4.65) | H₂O: 2.48e+27 ± 1.8e+27 (1.4σ, n_eff 2.5) marginal → rejected (+1.4σ); CO₂: < 1.2e+26 (+0.4σ) → < 1.2e+26 (+0.5σ) |
| 2025R2 | 2 | 4.3um: order 1 fixed; one-sided continuum accepted, no extension | 4.3um: FAIL (o1, 5 em, 3.18–4.65) → PASS (o1, 5 em, 3.90–4.65) | H₂O: < 2.2e+28 (-0.7σ) → < 7.4e+28 (-0.2σ); CO₂: not covered → 7.51e+25 ± 6.3e+24 (11.9σ, n_eff 2.7); CO: < 3.3e+26 (+0.1σ) → < 1.1e+27 (+0.0σ) |
| 2025W2 | 2 | 4.3um: fitted with >= 2 emission channel(s) | 4.3um: PASS (o1, 2 em, 3.90–4.65) → PASS (o1, 2 em, 3.90–4.65) | H₂O: 1.48e+26 ± 1.5e+25 (10.1σ, n_eff 3.2) → 1.48e+26 ± 1.5e+25 (10.0σ, n_eff 3.2); CO₂: not covered → < 2.3e+25 (+0.7σ) |

## 3. What the memo changed overall

| | 2026-09-14 run (`a6463b4c610b`) | `dc_main_norev` (new grouping, no revisions) | **2026-09-16 run (`5741d6611f60`)** |
|---|---|---|---|
| phases / fits | 193 / 147 | 192 / 147 | 192 / 147 |
| robust ≥ 3σ (n_eff ≥ 2) H₂O / CO₂ / CO | 17 / 26 / 4 | 17 / 26 / 4 | **18 / 28 / 4** |
| clean (carrier band PASS) | 17 / 25 / 2 | 17 / 25 / 2 | 18 / 27 / 2 |
| marginal 1–3σ | 25 / 23 / 9 | 25 / 23 / 9 | 24 / 26 / 9 |
| rejected | — | — | H₂O 2 (2023 RS61 S1, 2025 R1 S4) |
| covered H₂O / CO₂ / CO | 113 / 98 / 92 | 113 / 98 / 92 | 117 / 111 / 92 |
| comets with a robust / clean detection | 35 / 34 | 35 / 34 | **37 / 36** |
| median χ²_ν | 4.56 | 4.56 | 4.24 |
| continuum verdicts 2.7 / 4.3 / 4.7 µm (PASS-WARN-FAIL) | 95-14-52 / 113-8-40 / 82-30-49 | same | 94-18-49 / 120-8-33 / 81-31-49 |
| Afρ values at the SPHEREx phases (10k / 20k) | 116 / 120 of 165 | — | 120 / 124 of 164 |
| fitted phases with an Afρ value (10k / 20k) | 89 / 95 | — | 91 / 97 |

The 34 pipeline directives touch 36 band rows; 9 groups change their H₂O status and 16 their CO₂
status, none their CO status.  CO₂ gains the most: the seven waived negative continua and the
seven reduced coverage requirements turn "not covered" into six detections (2P S2, 217P S3,
2022 N2 S2, 2024 E1 S2, 2024 L5 S1, 2025 R1 S1, 2025 R2 S2) and five marginal values; two H₂O
detections are withdrawn as spurious, two more appear (2022 N2 S2 from the CO₂ column being
freed, 2022 QE78 S1 on two channels) and 63P S2 becomes a detection at the linear continuum.
The regrouping of 240P changes nothing in its rates (H₂O hot-band limits either way).  The
three Afρ overrides add four values (47P S1, 210P S1–S4 minus the one 240P phase that vanished)
and lower 217P S1 from 124 to 54 cm.

## 4. Caveats

1. **These are judgements, not rules.**  Each entry rests on the reviewer's reading of one
   figure; `dc_main_norev` keeps the rule-only result so the effect is always visible, and the
   registry must be revisited whenever the windows, the LSF or the photometry change (placeholder
   18).  Detections that exist only because of a waiver carry the directive in `caveats`.
2. **Waived negative continua** (2P S2, 217P S3, 229P S1, 2019 U5 S1/S2, 2024 E1 S2, 2024 L5 S4)
   subtract a continuum below zero under the band; 217P S3's is 7σ negative.  The CO₂ values on
   them are upper bounds of what a physical continuum would give.
3. **Reduced coverage** puts single channels in charge: 2022 N2 S2 (38σ, n_eff 1.0), 2024 L5 S1,
   2025 R1 S1, 2025 L1 S1 and 2022 QE78 S1 have n_eff ≤ 1.4 and stay outside the robust census
   (n_eff ≥ 2) by construction; they are detections of one channel each.
4. **Accepting a band can worsen the joint fit**: 2024 E1 S2 gains CO₂ (6.8σ) at the price of
   χ²_ν 33 → 97 and H₂O 12.8σ → 7.4σ, because the 10 CO₂ channels disagree with the band model
   more than the errors allow; 2025 R2 S2 (one-sided linear continuum) has χ²_ν 4 → 46.
5. **The 240P merged phase spans Δ by 26 %** (rule 1b exempted): the distance-corrected
   continuum absorbs an r_h²Δ² spread that every other group is protected from.
6. **210P's inbound Afρ assumes symmetric activity**; the outbound law has index 3.3 and S1 is
   0.05 dex beyond its range.  47P S1 rests on a grade-D law (index 11.7 over 0.09 au): 119 ± 8 cm
   is a formal number.
7. Two memo rows were applied through the emission-window count rather than the key range
   (2022 N2 S2, "+0.05 µm"; 2023 U1 S2, "to 4.5 µm"): the moved edge adds the channels the
   reviewer wanted in the fit, and the key range (which judges coverage) is untouched.
