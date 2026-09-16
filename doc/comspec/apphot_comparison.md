# Previous versus revised aperture photometry near the gas-emission windows

**Date:** 2026-09-09 · **Inputs:** `results/apphot/photometry/` (previous, 68 targets) and
`results/apphot/photometry/` (revised, `spherex_apphot` config `5502194856bc`) ·
**Scripts:** `notebooks/comspec/apphot_comparison.py` (matching and the four continuum
runs), `apphot_comparison_summary.py` (tables), `apphot_comparison_figures.py` ·
**Products:** `results/comspec/apphot_comparison/`, `fig/comspec/apphot_comparison/`.

Every one of the 27 797 exposures exists in both sets with identical wavelengths
and identical aperture centres (neither pipeline re-centres), so rows are matched
one to one on file name and aperture.  The comparison is made at each target's
rule aperture (20 000 km inside 3 au, 40 000 km beyond, promoted where the revised
set has no such aperture) inside the three emission windows (2.60–2.80,
4.18–4.35, 4.55–4.90 µm) and their continuum windows, and then propagated through
the pipeline's own continuum subtraction and production-rate fit in physical flux
space, in four runs:

| run | photometry | rows |
|---|---|---|
| a | previous | the rows the revised baseline policy keeps (`dc_main`: `frac_badpix_ap ≤ 0.05`, no flag `a`) |
| b | revised | the same rows |
| c | previous | the previous policy (`badphot` = flux ≤ 0 dropped, nothing else) |
| d | revised | the revised baseline policy |

a versus b isolates the photometry; c versus d is what each pipeline would have
delivered.

---

## 1. What actually changed between the two sets

Read from the data, not from the code:

1. **Sky annulus.**  Previous: fixed 150 000–300 000 km (7.5–15 × the 20 000 km
   aperture), scaled per exposure.  Revised: fixed 15–20 px, i.e. 105 000–300 000 km
   for 68 % of the sample and only 40 000–55 000 km for 24P.
2. **Bad pixels.**  Previous: counted (`nbadpix`) but *summed as zero*, so a bad
   pixel removed its share of the comet flux.  Revised: masked, and the sum is
   scaled by the effective aperture area (`aperture_area_eff_pix2`).
3. **Gaia stars.**  Previous: pixels near catalogued stars were *masked like bad
   pixels* — 82 % of the rows that only the previous set calls bad have a star
   within r_ap + 2 px (19 % for clean rows), and every flag-`a` row of the revised
   set has `nbadpix` > 0 in the previous set (median 8.4 pixels) while only 26 %
   have a bad pixel in the revised set.  Revised: pixels kept, row flagged
   `a`/`b`/`c`/`d`.
4. **`badphot`.**  Previous: flux ≤ 0 (3–37 % of rows per target).  Revised: any
   bad pixel in the aperture (up to 55 %), now relaxed to `frac_badpix_ap > 0.05`.
5. **Aperture availability.**  Revised refuses apertures below the PSF FWHM
   (0.86 px) or beyond the 15 px annulus: the 20 000 km aperture is absent for nine
   targets beyond ~4.5 au (promoted to 40 000–80 000 km) and 40 000 km is partly
   absent for the closest comets.  Raw emission-window channels at 20 000 km fall
   from 4 842 to 2 933 for that reason alone; it is a selection change, not a
   photometric one.
6. **Errors.**  Revised formal error = 0.71 × previous (16–84 %: 0.67–0.80), the
   double counting removed; revised empirical error = 0.76 × previous; the sky
   scatter itself agrees (1.03 ×); `sky_excess_ratio` ≈ 1.11.

Sample sizes per target barely move between c and d (153 groups, 459 band-rows
in both; per-target channel counts on the identity line in `per_target.png`).

## 2. Channel level inside the emission windows

`summary_row_level.csv`, channels with previous flux > 3σ, revised / previous:

| class of channel (share of rows) | emission windows | continuum windows |
|---|---|---|
| clean in both (56 %) | 0.989 (0.939–1.030) | 0.995 (0.941–1.066) |
| star-masked in the previous set only (29 %) | 1.106 (0.967–2.210), shifts up to +12σ | 1.289 (0.994–3.823) |
| bad pixels in both (15 %) | 0.976 (0.910–1.245) | 0.983 (0.914–1.394) |

The clean class shows the annulus: the ratio falls with aperture size in pixels
because the revised 15 px annulus moves into the coma of the close comets:

| aperture radius | 0.8–1.5 px | 1.5–3 px | 3–6 px | 6–13 px |
|---|---|---|---|---|
| revised inner annulus | 417 000 km | 170 000 km | 84 000 km | 40 000 km |
| flux ratio (clean) | 1.022 | 0.988 | 0.965 | 0.908 |

The sky level itself rises by +0.08σ per aperture (median; 16–84 % −0.20 to
+0.46σ) — small per channel, but systematic and largest for 24P, 306P, 2025 W2,
499P and 508P, the comets that carry most of the detections.

The star-masked class is the large effect: a bright star near the aperture cost
the previous set a median 10 % of the flux inside the emission windows and 29 %
in the continuum windows, with a tail to factors of several.  Emission and
continuum windows behave alike within every class (`windows_per_band.png`): the
change is a pixel-treatment effect, not a band effect.

## 3. Continuum subtraction on identical channels (a vs b)

* Verdicts agree for 401 of 459 band-rows (87 %); 21 accepted rows become FAIL and
  6 FAIL rows become accepted; 37 swap PASS and WARN.
* Band emission flux (both accepted, previous flux > 2σ), revised / previous:
  2.7 µm 0.985 (0.77–1.21, n = 37), 4.3 µm 1.025 (0.89–1.18, n = 63), 4.7 µm
  1.018 (0.72–1.18, n = 13).  In units of the previous error the change exceeds
  1σ for 29–34 % of the rows and 3σ for 7–13 %.
* The peak S/N rises by 1.12 / 1.23 / 1.27 × — almost entirely the smaller errors.

## 4. Production rates

`summary_Q.csv`, groups robust in both (detected, n_eff ≥ 2), revised / previous:

| species | same rows (a→b) | own selection (c→d) | robust previous → revised |
|---|---|---|---|
| H₂O | 0.983 (0.92–1.34), n = 23 | 0.976 (0.78–1.19) | 30 → 27 (same rows), 29 → 28 |
| CO₂ | **1.129 (1.005–1.216)**, n = 24 | **1.165 (1.02–1.41)** | 28 → 26, 34 → 25 |
| CO | 1.075 (0.87–1.58), n = 12 | 1.092 (0.81–2.08) | **13 → 21**, 12 → 21 |

Nothing moves beyond 4σ and only 4–8 % of the same-row pairs move beyond 3σ, but
the CO₂ shift is systematic and follows the affected channels: groups with fewer
than 20 % star-masked or bad-pixel channels in their emission windows have a
median ratio of 1.05, groups with 20 % or more have 1.16 (Spearman 0.29); the
largest are 24P phase 2 (1.40, 3.0σ; 61 % star-masked, 18 % bad pixels), 2025 L1
phase 2 (1.16, 3.4σ) and 2023 R1 phase 1 (1.16, 2.7σ).  The previous set's zeroed
pixels sat inside the 4.3 µm emission window often enough to depress the strongest
band of the brightest comets by 10–40 %.  CO gains nine robust detections on the
same rows: with 29 % smaller errors, marginal 4.7 µm detections cross the n_eff and
1σ tiers — an error-model effect, not new flux.  H₂O is unchanged.

Under each pipeline's own selection the picture is the same, plus a mild loss of
CO₂ robustness (34 → 25) because the revised set's FAIL count is higher (111 →
141): its smaller errors make the shape and χ² tests of the continuum validation
stricter at fixed thresholds.

## 5. Targets

`per_target_band_ratio.csv` (same rows, accepted, previous flux > 2σ; one or two
band-rows per target, so indicative only): the largest deviations are 40P (0.37)
and 2025 A6 (0.44), both with ≥ 50 % star-masked channels, and 2024 W1 (4.1),
161P (2.6), 2023 V1 (1.6) on single band-rows.  The bad-pixel-in-both class is
concentrated in 210P, 24P, 306P, 499P and 508P (35–50 % of channels), the
star-masked class in 491P, 2023 C2, 2024 J3, 2025 A6 and 229P (60–75 %).

## 6. Conclusions and caveats

1. **The revised set is the better photometry** and the previous set should not be
   used for gas production rates: zeroing bad and star-adjacent pixels removed
   comet flux inside the emission windows, which biased Q(CO₂) low by 13 % on
   average and by up to 40 % for 24P, and cost eight robust CO detections.
2. **The revised annulus is a new, smaller systematic in the other direction.**  A
   15 px inner radius is 40 000 km for 24P and 84 000 km at 3–6 px apertures, where
   the clean-channel flux drops by 3.5–9 % because the coma is in the sky.  For
   comets with extended CO₂ or H₂O comae this subtracts gas emission from the band.
   A km-based inner radius (the previous 150 000 km), or the annulus placed at a
   fixed multiple of the aperture, is the placeholder to update in
   `spherex_apphot`; the effect is measurable here because the previous set gives
   the same channels with a far annulus.
3. **Errors.**  The revised formal errors are 29 % smaller and the fits' χ²_ν of
   2.4–2.9 says they are still optimistic; the empirical column
   (`source_sum_err_empirical_mjy`, `Variant.error_column`) is the safer choice
   for detection tiers, which this comparison shows to be sensitive at the 1σ level.
4. **What does not change:** H₂O rates, the verdict of 87 % of band-rows, the
   sample size per target under each policy, and the aperture centres.

## 7. The annulus effect, measured (2026-09-12)

The photometry was regenerated with the sky annulus at 150 000 km (floored at 15 px,
capped at 40 px; `spherex_apphot` config `9fcc7ea3871a`).  Running the 2026-09-12
pipeline on both sets (`results/comspec/studies/annulus_previous/` holds the run on the archived
`5502194856bc` photometry) and pairing the 133 groups whose aperture did not change gives
Q ratios of 1.000 (CO₂, CO) and 1.001 (H₂O): the ring moves only for comets inside
~1.6 au, and there the four groups inside 1 au gain 4 % in Q(H₂O).  At the channel
level the 20 000 km fluxes of 306P and 508P rise by 2 and 6 %; the continuum subtraction
removes most of that from the band.  The concern of §6.2 is therefore closed at the
few-per-cent level; the larger effect of the new annulus is indirect — more valid large
apertures for close comets (the S/N-driven rule of 2026-09-12 used them; since 2026-09-14 the
aperture is fixed per phase at 20 000 / 40 000 / 60 000 km, `pipeline_decisions.md` §3).
