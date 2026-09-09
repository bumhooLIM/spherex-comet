# Reconstruction of the gas-production-rate pipeline on the revised photometry

**Package:** `spherex-comspec/` (`spherex_comspec` 1.0.0) · **Date:** 2026-09-08 ·
**Input:** `data/apphot_revised/` (68 comets, 27 797 exposures, `spherex_apphot` config `5502194856bc`)

This document records what changed when the three-notebook pipeline
(`phase_group_update` → `continuum_subtraction` → `gas_emission_fit` + `emission-fitter/`)
was rebuilt as one package on the revised aperture photometry, what the two
requested studies found — the distance correction and the Gaia source flags —
and which placeholders must be updated before the numbers are quoted.  How to
run the tool is in `spherex-comspec/README.md`; the physics is in
`model_concept.md` and `fitting_methodology.md`.

---

## 1. What the reconstruction reproduces

| stage | reproduction on the revised data |
|---|---|
| phase regrouping | **141 epochs → 174 phases**, 21 subdivided, 7 manual — the same count as before; group-by-group identical (r_h mean, n_meas, n_epoch) for **66 of 68** targets |
| continuum subtraction | 150 groups analysed, 24 skipped, **450 band-rows**: PASS 202 / WARN 108 / FAIL 140 (68.9 % saved; previously 152 / 456 / 73.0 %) |
| fitting | **136 target-phases fitted**, 14 with no band surviving validation (previously 139 / 13) |
| detections (`detected`, then `n_eff ≥ 2`) | H₂O 33 → 21 robust, CO₂ 71 → 25, CO 19 → 12 (previously 45 → 25, 77 → 34, 17 → 12) |

The method is unchanged.  The differences trace to two properties of the revised
photometry — the meaning of `badphot` and the aperture validity gate — both
discussed below.  The two targets whose grouping differs (2025 K1, 2025 L1) do so
because rule 3 counts which exposures sampled a band using `badphot`, whose
meaning changed; for 2025 L1 the new grouping removes a 1-epoch fragment the old
one carried.

---

## 2. The distance correction

### 2.1 What was done

`flux_distcorr_mjy = F × r_h² × Δ²` is the flux the comet would show at 1 au from
both Sun and observer.  The `dc_*` variants fit the continuum on it; `raw_all`
fits the same groups on the physical flux.

Q is a property of the comet and must not depend on the space the continuum was
fitted in.  The emission is therefore divided back by **each channel's own**
`distcorr_factor` before the design matrix is built, and the solve always runs
in physical flux space.  That division is the "re-correction" — there is no
separate step — and because the model is linear it is exact: fitting the
corrected emission with the design matrix scaled row-wise by the same factor
gives identical Q to 3 × 10⁻¹⁶ on real data (24P phase 2) and is asserted in the
tests.  Every retrieved Q is the production rate at the comet, at the geometry of
the observation.

### 2.2 Why the continuum space matters at all

SPHEREx scans a moving target non-simultaneously, so the channels of one group
were taken at different geometries.  Inside one 28-day epoch of 24P, r_h²Δ²
varies by 3×; the raw spectrum carries that geometric gradient across
wavelength, and a polynomial continuum absorbs part of it into the band.  The
regrouping bounds r_h to 10 % but says nothing about Δ, so **even inside the
fitted groups the per-channel spread of r_h²Δ² is 13 % (median), 53 % (90th
percentile) and up to 157 %**.

### 2.3 What it changes

Over the 150 groups both runs analysed (`results/comspec/distcorr_effect_*.csv`):

| | H₂O | CO₂ | CO |
|---|---|---|---|
| groups detected in both | 30 | 71 | 19 |
| lost / gained by the correction | 2 / 3 | 0 / 0 | 1 / 0 |
| median Q(dc) / Q(raw) | 1.000 | 1.000 | 1.000 |
| 16–84 % range of the ratio | 0.95 – 1.03 | 0.99 – 1.01 | 0.98 – 1.03 |
| fraction changed by > 1σ | 3 % | 0 % | 0 % |
| largest change | 1.06 σ | 0.68 σ | 0.49 σ |

Continuum verdicts: 439 of 450 band-rows keep their verdict; 7 PASS → WARN,
2 → FAIL, 2 FAIL → PASS/WARN.

The largest individual changes are 499P ph 3 CO₂ (×1.31), 2024 E1 ph 1 CO
(×1.28), 47P ph 2 H₂O (×0.73), 24P ph 2 CO (×0.74) — and |ΔQ/Q| correlates with
the group's geometric spread at r = 0.48, so the effect is real and is the one
the correction is meant to remove.  But no group moves by more than 1.1σ.

**Reading.**  At the catalog level the choice of continuum space is below the
noise; the correction changes individual Q by up to 30 % exactly where the
geometry inside a group is most heterogeneous, and does so in the direction the
physics predicts.  The corrected space is therefore the right default — it is
the one in which the channels of a group are on a common footing — and the
physical-space run is retained as the check that this is not an artefact.

---

## 3. Source-flag contamination

### 3.1 Setup

Gaia sources are never masked in the revised photometry; each row carries a
flag (`a` bright blend within r_ap + 2 FWHM, G_eff < 13; `b` Gaia flux > 20 % of
the comet's within r_ap + FWHM; `c` any source; `d` SNR < 1; precedence
a > b > c > d).  Three policies were run against the baseline, all in
distance-corrected space and sharing the same 174 groups
(`results/comspec/flag_policy_*.csv`):

| variant | rows dropped | groups analysed | band-rows saved | fitted | robust H₂O / CO₂ / CO |
|---|---|---|---|---|---|
| `dc_all` (baseline) | `badphot` | 150 | 68.9 % | 136 | 21 / 25 / 12 |
| `dc_no_a` | + flag `a` | 150 | 68.7 % | 136 | 21 / 24 / 13 |
| `dc_no_b` | + flag `b` | 140 | 56.2 % | 111 | 18 / 14 / 9 |
| `dc_no_ab` | + `a` and `b` | 134 | 55.7 % | 106 | 18 / 14 / 9 |

In the emission windows the baseline carries 3 815 channels, of which **89 are
flagged `a` and 1 274 are flagged `b`**.

### 3.2 Excluding flag `a`

Changes nothing measurable.  The same 136 groups are fitted; over the groups
detected under both policies the median Q ratio is 1.000 for all three species
with a 16–84 % range of 1.000–1.000 (H₂O, CO₂) and 0.97–1.00 (CO); one CO₂
detection is lost and one CO gained.  Only 2.3 % of emission-window channels
carry the flag at the 20 000–40 000 km apertures, and those that do are mostly
in groups the fit already down-weights.  **Excluding `a` is free and removes the
one contamination class that is unambiguous; there is no reason not to.**

### 3.3 Excluding flag `b`

Costs a quarter of the sample — 25 fewer fitted groups, a third of the
emission-window channels, robust CO₂ detections 25 → 14 — and for the groups
kept by both policies the Q ratios are consistent with unity: median 1.000, 16 %
quantile 0.69 (H₂O), 0.87 (CO₂), 0.83 (CO), no group changed by more than 3σ,
and only 14 % of H₂O / 4 % of CO₂ / 0 % of CO groups changed by more than 1σ.
Detections lost and gained are balanced (3 / 1 for CO₂, 3 / 3 for CO), which is
what a cut that removes points at random with respect to contamination looks
like.  `b` fires on 33 % of channels because the criterion (Gaia flux > 20 % of
the comet's predicted V-band flux) is met by nearly any catalogued star near a
faint comet.  **Excluding `b` shrinks the sample without changing the answer;
keep it as a column, do not cut on it.**  This is the same conclusion the
survey-wide growth-curve test reached in the photometry project (lift 5.6 for
`a` alone, 1.7 for `a + b`).

### 3.4 The `badphot` policy is the cut that actually matters

`badphot` in the revised photometry marks *any* bad pixel inside the aperture;
in the old tables it was `source_sum ≤ 0`, a flux-sign cut.  The strict rule
removes a median 48 rows per target and **52 % of 24P's rows at 20 000 km**,
where any ~50-pixel aperture contains a bad pixel — yet 85 % of those rows lose
less than 5 % of their area (median `frac_badpix_ap` = 0.021).  A `dc_all_lenient`
variant that drops a row only above `frac_badpix_ap > 0.05` was run alongside:

| | strict | lenient |
|---|---|---|
| groups analysed / band-rows saved | 150 / 68.9 % | 155 / 69.5 % |
| fitted | 136 | 139 |
| detected H₂O / CO₂ / CO | 33 / 71 / 19 | 43 / 80 / 24 |
| robust (n_eff ≥ 2) | 21 / 25 / 12 | **28 / 25 / 19** |
| Q ratio lenient / strict, 16–84 % | 0.96–1.02 / 0.90–1.02 / 0.83–1.03 | |
| median Q_err ratio | 0.97 – 0.99 | |
| median χ²_ν | 2.9 | 2.8 |

Ten more H₂O and five more CO detections, seven more robust H₂O and CO each,
Q unchanged within errors for the groups both policies fit, and a slightly
*better* χ²_ν.  **This is the cut to revisit first** (placeholder 5): the strict
rule is implemented as specified and remains the default, but the evidence is
that it discards a third of the usable signal on the brightest comets for no
gain in accuracy.

---

## 4. Concerns found during the reconstruction

1. **`jd_utc` in the revised photometry is rounded to 0.1 day.**  The tables were
   written with `float_format="%.8g"`, which truncates a Julian date of ~2.46 × 10⁶
   to one decimal; 765 exposures of 24P collapse onto 135 distinct times, each
   spanning several distinct r_h.  The package rebuilds the exact time from the
   `date_obs` string and keys exposures on `obsid`, and carries the exact time
   into every product.  *The fix belongs upstream:* `spherex_apphot` should
   exempt the JD columns from `float_format` (or use `%.12g`).
2. **`badphot` changed meaning** (§3.4).  Any downstream code that treated it as a
   flux-sign cut now removes half of a bright comet's rows.
3. **The aperture rule can select an aperture that does not exist.**  The revised
   photometry refuses apertures below the PSF FWHM or beyond the 15 px sky
   annulus, so 2014 UN271 (14.6 au) has no 40 000 km measurement at all and
   2019 U5, 2022 R3, 2023 RS61 have it for only 79–88 % of exposures.  The
   package promotes to the smallest larger `km` aperture with ≥ 95 % coverage
   (80 000 / 60 000 km) and logs it.  The model's `rho_ap_km` follows.
4. **Δ is unconstrained by the grouping rules**, so the geometry inside a
   "single-state" group is heterogeneous by up to 157 % in r_h²Δ² (§2.2).  The
   distance-corrected continuum is the mitigation; a Δ-aware rule would be the
   cure.
5. **Formal photometric errors are a lower bound.**  The revised photometry reports
   `sky_excess_ratio` ≈ 1.2 (the VARIANCE plane under-reports the annulus
   scatter), and the continuum `err_scale` and the fit's χ²_ν ≈ 2.9 both absorb
   it after the fact.  `source_sum_err_empirical_mjy` is available as an
   alternative (`Variant.error_column`).
6. **The continuum uncertainty is correlated across a band but propagated as
   diagonal** (fitting_methodology §10.2a).  The coefficient covariance is saved
   in every summary row, so a generalised least-squares solve is a contained
   change to `fitting.py`.
7. **`detected` means ≥ 1σ** (`FitConfig.upper_limit_sigma = 1.0`) and gates the
   mixing ratios; the robust counts quoted here additionally require n_eff ≥ 2.
8. **The 2.7 µm blue edge at 2.60 µm** loses 6.7 % of the H₂O complex and leaves a
   `cont_used` channel at 2.59 µm carrying 23 % of the peak (project handoff);
   the window is a config value and was left at the notebook value for a faithful
   reproduction.
9. **FAIL is 31 % of band-rows** (27 % before).  As before, ~88 % of failures
   involve a continuum sampled on one side of the band; the stricter `badphot`
   removes continuum points as readily as emission points.
10. **Flag `a` is too rare at these apertures to test** (89 channels); its
    effectiveness was established in the photometry project's growth-curve
    study, not here.
11. **A group can lose every channel to the negative-channel cut.**  499P phase 3
    in `dc_all_lenient` keeps only its 4.7 µm band (2-point continuum, WARN), and
    that continuum sits above all five of its emission channels, so the 1σ cut
    leaves nothing to fit.  The fitter now reports such a group as not fitted,
    with the reason in `not_fitted.csv`, instead of failing inside the filling
    factor, and the figure stage skips it.  It is a continuum failure, not a
    non-detection, and it is the one group that stopped the figure batch twice.

---

## 5. Placeholders — what to update, in order

`config.PLACEHOLDERS` (also `results/comspec/placeholders.csv`) lists 15.  The
ones this run gives direct evidence on:

| priority | placeholder | current value | evidence from this run | update |
|---|---|---|---|---|
| 1 | SPHEREx LSF | Gaussian of FWHM = `wlwidth` | CO / H₂O-hot separability at 4.7 µm depends on the wings | as-built LVF R(λ) and LSF |
| 2 | band profiles Φ_b | Gaussian, 0.02–0.10 µm | neutral for Q (unit area); wrong for any band-shape claim | PSG / GSFC templates |
| 3 | expansion velocity | 0.8 r_h⁻⁰·⁵ km/s, unattributed | ±18 % in Q(CO₂) per 20 % in v_g; does not cancel in ratios | cited species-specific law |
| 4 | 2.7 µm blue edge | 2.60 µm | loses 6.7 % of the band; 1 comb offset in 5 leaves < 3 channels | 2.50 µm |
| **5** | **`badphot` policy** | drop any bad pixel | **strict rule discards 52 % of 24P; lenient (frac > 0.05) recovers +7 robust H₂O, +7 CO at unchanged Q** | `max_frac_badpix = 0.05` |
| 6 | aperture rule | 20 k / 40 k km at 3 au, promoted for coverage | 4 targets promoted; `rho_ap_km` is physics, not a label | S/N-driven per target |
| 7 | polynomial orders | 3 / 2 / 2 | `cv_best_order` disagrees in most WARNs | re-tune from `cv_best_order` |
| 8 | detection tier | 1σ | see concern 7 | 3σ or rename |
| 9 | negative-channel cut | 1σ one-sided | biases Q up ~0.29σ/channel (noise) | noise-injection test |
| 10 | error column | `source_sum_err_mjy` | `sky_excess_ratio` ≈ 1.2; χ²_ν ≈ 2.9 | empirical error |
| 11–15 | grouping thresholds, sufficiency gates, opacity calibration, T_rot, upstream flag thresholds | as in the notebooks / `spherex_apphot` | no new evidence | see registry |

Left as defaults deliberately: everything the user's specification fixed
(strict `badphot`, the 2.60 µm edge, the 20 k / 40 k rule, the 1σ tier).  Each
is a one-line change in `config.py`, and the variant hash written to every
`run.meta.json` records which value produced which file.

---

## 6. Products

```
data/comspec/phase_assignment.csv                    per-exposure phase labels (never written into apphot)
results/comspec/phase_map.csv, phase_cuts.csv        174 groups; every cut and what it cost
data/comspec/<variant>/emission/                     continuum summaries + point spectra, per target
results/comspec/<variant>/gas_fit.csv                Q per (target, phase)            <- main: dc_all
results/comspec/<variant>/gas_fit_lines/             model curves and residuals per fit
results/comspec/<variant>/{continuum_summary,skipped_groups,not_fitted,apertures}.csv, run.meta.json
results/comspec/flag_policy_{census,pairs,paired_Q}.csv
results/comspec/distcorr_effect_{summary,paired_Q,band_rows,verdicts}.csv
results/comspec/placeholders.csv
fig/comspec/phase_group/                             68 before/after grouping figures + summary
fig/comspec/<variant>/cont_subtract/                 raw spectrum + 3-band validation grid per group
fig/comspec/<variant>/emission_model/                model-over-data per fit
fig/comspec/<variant>/summary_{Q_vs_rhel,mixing_ratios}.png
fig/comspec/flag_policy_comparison.png, distcorr_effect.png
```
