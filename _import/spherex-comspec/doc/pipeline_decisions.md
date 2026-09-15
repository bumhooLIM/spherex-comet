# Pipeline decisions — SPHEREx comet gas production rates

**Date:** 2026-09-09 · **Package:** `spherex-comspec/` (`spherex_comspec` 1.0.0, run hash
`bb0944596af3`) · **Input:** `data/apphot/` (revised `spherex_apphot` photometry, config
`5502194856bc`; 68 comets, 27 797 exposures) · **Result:** `results/gas_fit.csv`.

This is the record of what the pipeline does, which alternatives were tried, why the
current configuration was kept as the most robust one, what remains a concern, and
which placeholders still have to be replaced.  The physics is in `model_concept.md`
and `fitting_methodology.md`; the photometry comparison is in `apphot_comparison.md`.

---

## 0. The pipeline, and the conclusion

Revised photometry → one aperture per comet → phases of a single observing state →
per-band polynomial continuum in distance-corrected flux → joint weighted linear
least squares for Q(H₂O), Q(CO₂), Q(CO) in physical flux space, with the 2.7 µm band
anchoring H₂O and the 4.6–4.9 µm hot bands taking over only where 2.7 µm is not
covered.  Six configurations were run end to end and compared group by group
(§7).  **The configuration kept as the main result (`dc_main`) is:**

| choice | value | why it won |
|---|---|---|
| photometry | revised set (`data/apphot/`) | the previous set zeroed bad and star-adjacent pixels inside the aperture, which biased Q(CO₂) low by 13 % and cost 8 robust CO detections (`apphot_comparison.md`) |
| `badphot` rule | drop a row only when `frac_badpix_ap > 0.05` | the strict rule (any bad pixel) removed 52 % of 24P's rows for no change in Q; the lenient rule analyses 154 instead of 150 groups and fits 138 instead of 135, with Q identical (§7.2) |
| source flags | drop `a` (G < 13 star within r_ap + 2 FWHM), keep `b`, `c`, `d` | dropping `a` costs 3 band-rows and nothing in Q; dropping `b` removes 28 fits and a third of the channels, again without moving Q (§7.1) |
| continuum flux space | distance-corrected (F × r_h² Δ²) | the channels of one group span up to 166 % in r_h² Δ²; the corrected continuum removes that gradient and changes Q by ≤ 1.35σ against the physical-space fit (§7.3) |
| error column | `source_sum_err_empirical_mjy` | the formal error under-reports the annulus scatter; the empirical one leaves Q unchanged and lowers χ²_ν from 2.44 to 2.27 (§7.4) |
| H₂O coverage | 2.7 µm main band, hot-band fallback | the bright, close comets (24P at perihelion, 306P) lose the 2.7 µm channels; the fallback recovers their Q(H₂O) and removes the CO bias those groups carried (§7.5) |

It analyses every one of the 68 comets, fits 138 (comet, phase) groups, and yields a
robust detection (≥ 1σ, n_eff ≥ 2) of at least one species for 41 comets: H₂O 40
groups (23 from the main band, 17 from the hot bands), CO₂ 23, CO 18.  The
alternatives remain runnable as study variants (`config.DEFAULT_VARIANTS`,
products under `results/studies/`), so every number here can be regenerated.

## 1. Inputs

`data/apphot/` carries one row per (exposure, aperture) with the revised
treatment: predicted ephemeris centres, a fixed 15–20 px sky annulus, bad pixels
masked with an effective-area correction, Gaia sources flagged (`sourceflag`
`a`/`b`/`c`/`d`) rather than masked, a distance-corrected flux column and both a
formal and an empirical error.  Two properties matter downstream: `badphot` means
*a bad pixel in the aperture*, not a negative flux; and `jd_utc` is rounded to
0.1 day, so exposures are keyed on `filename` and timed from `date_obs`.

## 2. Phase grouping (`grouping.py`)

The 28-day `epoch` blocks are cut into phases that hold one physical state: r_h
spread < 10 %, a perihelion passage split (resolved when both sides move by
> 0.001 au and have ≥ 2 epochs), emission bands never cut, the fewest groups with
cuts placed in gaps (link 0.05 au), manual edges for 24P and 2024 E1.  Result: 141
epochs → **174 phases**, identical to the previous notebook's map for 66 of 68
comets (2025 K1 and 2025 L1 differ because rule 3 reads band sampling through
the new `badphot`).  `results/phase_map.csv`, `phase_cuts.csv`.

## 3. Aperture (`dataio.aperture_for`)

20 000 km inside 3 au, 40 000 km beyond.  The revised photometry refuses apertures
below the PSF FWHM or beyond the annulus, so the rule aperture can be absent; the
smallest larger `km` aperture present for ≥ 95 % of the exposures is used instead
(2014 UN271 → 80 000 km; 2019 U5, 2022 R3, 2023 RS61 → 60 000 km).  The model's
`rho_ap_km` follows.  33 comets at 20 000 km, 31 at 40 000 km, 4 promoted.

## 4. Row selection (`config.BASELINE_FLAGS`)

A row enters the spectrum unless `frac_badpix_ap > 0.05` or `sourceflag = a`.
Nothing is written into the photometry; the selection and its rejection counts
travel with every product.

## 5. Continuum subtraction (`continuum.py`)

Per band a weighted polynomial (orders 3 / 2 / 2 at 2.7 / 4.3 / 4.7 µm) in λ − λ_c
over the continuum window with all three emission windows punched out, 3σ
clipping, order capped at 1 when the continuum is one-sided, a two-point fallback
below six points.  Validation: bracketing and positivity are hard, shape and
cross-validation soft → PASS / WARN / FAIL.  Main run: 462 band-rows, 227 PASS /
97 WARN / 138 FAIL (70 % saved); 20 groups skipped for insufficiency.  The emission
is divided back by each channel's own `distcorr_factor` before the fit
(`emis_raw_mjy`), which is what keeps Q at the comet.

## 6. Fitting (`fitting.py`)

Optically thin Haser coma with the Yamamoto filling factor, g-factors ∝ r_h⁻², so
the model is linear in Q and the three species are solved jointly by weighted
least squares over every accepted emission channel, each at its own geometry and
bandpass.  Coverage is decided on the diagnostic ranges (H₂O ≥ 3 channels in
2.60–2.80 µm, CO₂ ≥ 2 in 4.20–4.30, CO ≥ 2 in 4.60–4.70); *not covered* is
categorically different from *not detected*.  The covariance is rescaled by χ²_ν
when > 1; `n_eff`, the participation ratio of the per-channel Fisher information,
says how many channels really carry a species.  Statuses: `detected` (≥ 1σ),
`upper_limit`, `negative_fit`, `not_covered`; "robust" throughout this document
means detected with n_eff ≥ 2.

**H₂O hot-band rule (2026-09-09).**  When the 2.7 µm band is not covered and at least
three channels lie in 4.55–4.90 µm with one beyond 4.75 µm (the 4.85 µm ν₁−ν₂ band,
outside CO v(1−0)), the ν₃−ν₂ and ν₁−ν₂ hot bands carry Q(H₂O).  The row is
labelled `h2o_source = "hot"`, the caveat says so, the figures draw it open.  Such
values are provisional: hot-band g-factors are placeholder 2.

## 7. What the studies showed

All comparisons are group by group against `dc_main`; "Q ratio" is the median
over groups detected in both runs with its 16–84 % range; "max" is the largest
shift in units of the larger error.  Tables: `results/studies/flag_policy_*.csv`,
`distcorr_effect_*.csv`, `apphot_comparison/`.

### 7.1 Source flags

| run | groups / fits | robust H₂O / CO₂ / CO | Q ratio | max |
|---|---|---|---|---|
| `dc_main_keep_a` (flag `a` kept) | 155 / 138 | 40 / 21 / 18 | 1.000, 1.000–1.000 (all) | 1.7σ |
| `dc_main_no_b` (flag `b` dropped too) | 138 / 110 | 29 / 13 / 11 | 1.000; H₂O 0.64–1.03, CO₂ 0.90–1.06, CO 1.00–1.21 | 3.0σ |

Flag `a` (342 of 12 438 channels) changes nothing; it is dropped because it marks
the one contamination mechanism the photometry project demonstrated.  Flag `b` is
too common to cut on: it removes a third of the channels and 28 fits and leaves
the surviving Q unchanged.

### 7.2 `badphot` rule

| run | groups / fits | robust H₂O / CO₂ / CO | Q ratio | max |
|---|---|---|---|---|
| `dc_main_strict` (any bad pixel drops the row) | 150 / 135 | 36 / 25 / 10 | 1.000; 0.98–1.14 / 0.99–1.10 / 0.96–1.15 | 4.1σ |
| `dc_all` (strict + every flag; the 2026-09-08 baseline) | 150 / 134 | 37 / 23 / 11 | 1.000; 0.98–1.16 / 0.98–1.10 / 0.94–1.15 | 5.1σ |

The strict rule discards 52 % of 24P's rows (a 50-pixel aperture almost always
contains one bad pixel) and costs four groups, three fits and eight robust CO
detections at unchanged Q.

### 7.3 Continuum flux space

`dc_main` versus `raw_main` (same policy, physical flux): Q ratio 0.998 (0.956–1.044)
for H₂O, 1.000 (0.990–1.003) for CO₂, 1.000 (0.976–1.078) for CO; nothing beyond
1.35σ; 444 of 462 band-rows keep their verdict.  The in-group spread of r_h² Δ²
has a median of 13 %, a 90th percentile of 58 % and a maximum of 166 %, and the size
of the Q change correlates with it (r = 0.23).  The two spaces are equivalent for
Q; the corrected space is kept because it removes a known gradient instead of
leaving it to the polynomial.

### 7.4 Error column

Empirical against formal errors (`_archive/results_comspec_2026-09-09/dc_main`):
Q ratio 0.999 / 1.003 / — for H₂O / CO₂, Q errors unchanged (the χ² rescaling had
been absorbing the difference), χ²_ν median 2.44 → 2.27.  The fits still run at
χ²_ν ≈ 2.3: the errors remain optimistic by ~1.5 in scatter, which the rescaling
carries into the quoted Q errors.

### 7.5 H₂O hot-band fallback

Before the rule, 16 fitted groups had CO detected with H₂O not in the model, so the
4.63 µm hot band was absorbed into Q(CO); now 1.  Robust H₂O rose from 28 to 40
(17 from hot bands: 24P phases 3–5 at 1.19–1.21 au, 306P phases 2–3, 2024 E1,
2025 L1, …), CO fell from 21 to 18 and the CO values shared by both runs dropped
by 9 % (0.73–1.03), the removed contamination.  The hot-band values are weaker than
the main-band ones (median 1.8σ against 2.9σ): ten of the 17 lie inside 3 au and two
of those reach 3σ (24P phases 4 and 5); the seven beyond 3 au (up to 29P at 6.3 au
and 2019 U5 at 9 au) are all 1–2σ with errors comparable to the value, and are
not to be quoted as water production rates without a distance cap or the 3σ tier.

### 7.6 Previous photometry

`apphot_comparison.md`: on identical channels the revised photometry changes
Q(H₂O) by 0.98, Q(CO₂) by 1.13 (1.00–1.22) and Q(CO) by 1.08; the previous set's
zeroed bad and star-adjacent pixels (29 % of emission-window channels) depressed
the 4.3 µm band of the brightest comets by 10–40 %.

## 8. Concerns

1. **The revised annulus** (fixed 15–20 px) sits at 40 000 km for 24P and 84 000 km
   for 3–6 px apertures, where the coma is in the sky: −3.5 to −9 % on clean
   channels.  A km-based inner radius in `spherex_apphot` is the fix.
2. **`jd_utc` is rounded to 0.1 day** upstream (`float_format="%.8g"`).
3. **Δ is unconstrained by the grouping** (r_h² Δ² spread up to 166 %).
4. **Errors are optimistic** (χ²_ν ≈ 2.3 even with the empirical column) and the
   continuum uncertainty is correlated across a band but propagated as diagonal
   (the coefficient covariance is saved; generalised least squares is a contained
   change).
5. **`detected` means ≥ 1σ**; the robust counts add n_eff ≥ 2, and 23 comets have a
   ≥ 3σ robust detection.
6. **Hot-band Q(H₂O) is provisional** (placeholder g-factors, shared feature with
   CO) and, beyond 3 au, weak (seven values, none ≥ 3σ).  It is labelled everywhere
   and easy to exclude (`h2o_source == "main"`); a `max_rh_au` in `H2O_HOT_RANGE`
   or the 3σ tier is the next decision.
7. **FAIL is 30 % of band-rows**, mostly one-sided continua; the 2.7 µm blue edge
   at 2.60 µm loses 6.7 % of the H₂O complex.
8. **A group can lose every channel to the 1σ negative cut** (499P phase 3, 2024 L5
   phase 3): reported as not fitted with the reason.
9. **Rule 3 of the grouping reads the raw `badphot`**, so the assignment does not
   move with the row policy (deliberate; 2025 K1 is the case to check).

## 9. Placeholders

Applied: 5 (`badphot` → `frac_badpix_ap > 0.05`), 10 (empirical error column).
Open, in priority order (`results/placeholders.csv`, `config.PLACEHOLDERS`):

| priority | placeholder | current value | evidence | update |
|---|---|---|---|---|
| 1 | SPHEREx LSF | Gaussian, FWHM = `wlwidth` | CO / H₂O-hot separability at 4.7 µm depends on the wings | as-built LVF R(λ) and LSF |
| 2 | band profiles | Gaussian, 0.02–0.10 µm | hot-band Q(H₂O) rests on them | PSG / GSFC templates |
| 3 | expansion velocity | 0.8 r_h⁻⁰·⁵ km/s | ±18 % in Q per 20 % in v_g | cited species-specific law |
| 4 | 2.7 µm blue edge | 2.60 µm | loses 6.7 % of the band | 2.50 µm |
| 6 | aperture rule | 20 k / 40 k km, promoted | 4 comets promoted | S/N-driven choice |
| 7 | polynomial orders | 3 / 2 / 2 | `cv_best_order` disagrees in most WARNs | re-tune |
| 8 | detection tier | 1σ | 23 of 41 comets survive 3σ | 3σ, or rename the tier |
| 9 | negative-channel cut | 1σ one-sided | biases Q up ~0.3σ per channel | noise-injection test |
| 11–15 | grouping thresholds, sufficiency gates, opacity calibration, T_rot, upstream flag thresholds | notebook / `spherex_apphot` values | no new evidence | see registry |
| new | sky annulus (upstream) | 15–20 px | concern 1 | ≥ 150 000 km or a multiple of the aperture |

## 10. Products

`results/gas_fit.csv` (one row per comet and phase; columns in `results/README.md`),
`results/gas_fit_lines/`, `results/continuum_summary.csv`, `skipped_groups.csv`,
`not_fitted.csv`, `apertures.csv`, `phase_map.csv`, `phase_cuts.csv`,
`placeholders.csv`, `run.meta.json`; `data/phase_assignment.csv`, `data/emission/`;
`fig/emission_model/`, `fig/cont_subtract/`, `fig/phase_group/`,
`fig/summary_Q_vs_rhel.png`, `fig/summary_mixing_ratios.png`.  Every study variant
has the same set under `results/studies/<variant>/`, `data/studies/`, `fig/studies/`;
the cross-variant tables and figures sit directly in `results/studies/` and
`fig/studies/`.  Everything older is in `_archive/`.
