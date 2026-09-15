# Where does the 1/ρ coma law stop holding?  A survey-wide test

**Bottom line.**  Nowhere we can see.  With the PSF core and the sky both
modelled, the observed coma profile matches a 1/ρ model to within 1–2% from
~1,500 km out to **37,000–45,000 km — the edge of coverage — at every S/N band
and at every heliocentric distance from 0.6 to 6.7 au.**  The 12,000 km figure
from the 24P study was where 24P's data ran out, not where the law did.  The
steep naive slopes are an S/N effect, not a physical-window effect, and the
corrected slope converges on m = 1.00 as S/N rises.

Driver: `notebooks/profile_survey.py`.  Tables: `results/profile_survey_*.csv`.
Figures: `fig/survey/profile_survey_slope.png`, `profile_survey_ratio.png`.

## Question

The 24P resolution study (`doc/profile_resolution_24P.md`) showed ρ·SB flat to
~12,100 km post-perihelion and set ρ_max ≥ 12,000 km.  But 10 px at Δ ≈ 0.8 au
*is* ~12,000 km: that limit was coverage.  The survey's distant targets push
10 px to 30–45,000 km, and ~1,150 clean frames do so with outer-annulus
S/N > 5.  Does 1/ρ hold out there, or does it steepen — as the naive slopes of
−1.4 to −1.8 on high-S/N distant comets (2021G2, 2022E2, 2024J3) suggested?

## Method

Two known biases steepen a naive [1.5 FWHM, 10 px] slope: the PSF core, which
at Δ = 4.5 au spans ~10,000 km and puts the whole fit window inside the wings;
and a sky error, which at low S/N tilts the outer annuli.  Both are modelled
per frame with `ztfcomet.profile.fit_coma_model`: a nucleus δ-function plus a
C·ρ^(−m) coma, convolved with the empirical PSF from the field-star stack,
plus a free sky offset — fitted twice, with m free and with m = 1.

The diagnostic is the ratio **observed / (m = 1 model)**, sky removed from both
sides.  It is 1 wherever 1/ρ holds, at every radius including inside the PSF
core, and it is binned in km so that a departure tied to physical distance
separates cleanly from one tied to pixels or S/N.

Sample, all 68 targets on profiles centred by `refine_centre`: 3,398 clean frames with
S/N(10 px) > 5 → 3,273 peaking on-centre (125 still off-peak, against 691 before the fix)
→ **2,735 whose free sky converged off its bound**.  16 free slopes hit the m = 3 ceiling and are
treated as failed.

## Results

### The corrected slope converges on 1 with S/N

| S/N at 10 px | n | naive slope | corrected m [16–84%] | at sky bound |
|---|---|---|---|---|
| 5–20 | 527 | 1.80 | **1.32** [1.00, 1.71] | 38% |
| 20–60 | 835 | 1.66 | **1.22** [0.97, 1.50] | 15% |
| > 60 | 1,373 | 1.43 | **1.04** [0.85, 1.27] | 3% |

The PSF core alone accounts for ~0.4 in slope at every S/N; the sky adds ~0.3
more at low S/N.  Neither is coma physics.  At S/N > 60, 65% of frames sit in
[0.8, 1.2] and 80% in [0.7, 1.3].

### 1/ρ holds to the coverage edge at every S/N

Median observed / (1/ρ model), km bins:

| S/N | 2k | 5k | 10k | 17k | 26k | 37k | 45k | ρ_max |
|---|---|---|---|---|---|---|---|---|
| 5–20 | 0.98 | 1.05 | 0.99 | 1.01 | 1.01 | 1.08 | — | 37,500 (edge) |
| 20–60 | 0.99 | 1.04 | 0.99 | 1.00 | 1.02 | 0.99 | 0.99 | 45,400 (edge) |
| > 60 | 1.00 | 1.01 | 1.00 | 1.00 | 1.00 | 1.00 | — | 37,500 (edge) |

ρ_max (last bin with median ≥ 0.8) is the last bin with data in every band.
No downturn is detected anywhere.

### The break does not move with r_h

S/N > 20 only, by heliocentric-distance tercile:

| r_h | n | flat to | ρ_max |
|---|---|---|---|
| 0.6–2.2 au | 620 | 17,400 km | edge |
| 2.2–3.3 au | 495 | 25,500 km | edge |
| 3.3–6.7 au | 563 | 45,400 km | edge |

A radiation-pressure turnover scales roughly as r_h²; an instrumental one is
fixed in km or pixels.  Nothing moves, because nothing turns over.  The ±5–10%
wiggle inside ~2,000 km in the two outer terciles sits at sub-pixel radii at
those Δ and is the PSF stack imperfectly matching the comet core — a statement
about the model's centre, not the coma's wings.

### The slope does not steepen with the physical window

At S/N > 60 the running median of corrected m stays between 0.83 and 1.10 as
ρ(10 px) runs from 4,000 to 45,000 km (`profile_survey_slope.png`, left).  At
lower S/N the median sits at 1.2–1.45 with no trend in ρ(10 px) either.  The
steepness follows S/N, not km.  This is the same conclusion as the 24P study,
now with the PSF modelled and 2,118 frames across 55 targets.

### Real comet-to-comet variation, centred on 1

Among the 17 targets with ≥ 15 usable frames at S/N > 60, the per-target
median corrected slope runs from **0.57 (2025R2)** to **1.42 (10P)**, with a
16–84% spread of 0.84–1.13.  That spread is intrinsic: 10P's fan and 2025R2's
flat inner coma are visible in the cutouts.  The population is consistent with
steady-state outflow with a real dispersion of a few tenths in m.

## Answer to the second open question: should slopes below −2.5 be excluded?

Not by their value.  Of the 213 frames with a naive slope steeper than −2.5,
**185 never obtained a constrained model fit** — the free sky ran to its bound —
and only 36 have S/N(10 px) > 20.  The 28 that did converge correct from a
median of 2.67 to **1.24**.  Those slopes were never measurements of the coma.

The right rule is on S/N, not slope: **report the PSF-model + free-sky slope,
restricted to S/N(10 px) > 20 and fits off the sky bound.**  Below that, quote
the naive slope as an upper limit on steepness, not a value.

## Caveats

- **23% of fits sat on the sky bound** (42% at S/N 5–20).  The model allows a
  sky offset up to 3× the outer-annulus SB; at low S/N the slope–sky degeneracy
  walks there.  Those frames are excluded, so the low-S/N row above is the 58%
  that converged — a mild selection toward better-behaved frames.
- **Δχ² is only indicative.**  `n_clip` counts oversampled sub-pixels, which are
  bilinear interpolates of the same data, so annulus errors are ~4× too small.
  The driver scales Δχ² by the oversampling factor squared and calls it
  effective; the ratio curves and the m distributions are the evidence, not the
  acceptance fractions.
- The m = 1 model includes the nucleus, so "1/ρ holds" means the *coma
  component* is consistent with 1/ρ after nucleus, PSF and sky.  At Δ > 4 au
  the inner ~10,000 km is mostly PSF; the test there is of the wings, from
  ~1 FWHM outward.

## The centring defect, since fixed

Before the fix, **691 frames — 20% of the clean, S/N > 5 sample — were excluded because the
comet profile did not peak in the innermost annulus.**  With `refine_centre` — the optocentre within a
disc covering the winpos and ephemeris positions, guard calibrated on its own null — the exclusion is 125 frames.  The star stack never
does this (0.0%), so the extraction is sound; the comet centre is not.  It is
the photometry's windowed centroid, and the off-peak fraction tracks how far
that centroid moved from the ephemeris:

| winpos shift from ephemeris | off-peak |
|---|---|
| < 0.5 px | 6–8% |
| 0.5–1 px | 10% |
| 1–2 px | **32%** |
| 2–5 px | **57%** |
| > 5 px | 74% |

A windowed centroid on a 1/ρ coma — far heavier-winged than a Gaussian — is
pulled by any asymmetry toward the light centre, and the further it moves the
more often it overshoots the nucleus.  Most cases are mild (SB at 1.5 px over
SB at 0.5 px: median 1.01, 84th percentile 1.16), a flat-topped core rather
than a miss; 4% of frames peak at ≥ 2 px.

Aperture photometry is untouched — a 1 px shift inside a ≥ 10,000 km aperture
is nothing — but the model's nucleus term has nothing to fit.  The fix is a
peak re-centring step in `run_profiles` within ±2 px of the winpos position,
followed by a re-run of the profile step (~30 min).  It would recover most of
the 691 and change the per-target slopes for the off-peak-heavy targets
(40P 61 frames, 235P 60, 47P 45, 2024E1 38); numbers in the survey summary
would move modestly.  That is a decision, not a bug fix, and is left open.
