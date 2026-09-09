# ZTF 68-comet survey — process and flag summary

Run: 2026-09-07 17:25 → 2026-09-08 22:32 (one resume, one restart).
Window 2025-03-01 → present, Vmag < 20, cutouts 5′ (10′ when Vmag_min < 14
or Δ_min < 1 au).  Apertures ρ = 10/15/20/30/40 ×10³ km, each kept only where
`PSF_FWHM < r_ap < 60″`.

Numbers below are computed from `results/*/photometry_*.csv` and
`results/*/profile_summary_*.csv` on the local disk.  The authoritative
`survey_status.csv` and `survey_summary.md` live on the T7, which was
unmounted when this was written.

## Coverage

| | |
|---|---|
| designations processed | 68 / 68 |
| with photometry | 55 |
| no frames / no epochs | 13 |
| unique frames measured | 8,336 |
| measurement rows | 39,591 |
| filters | ZTF_r 4,388 · ZTF_g 2,978 · ZTF_i 970 |

Targets that produced no photometry — either nothing reached Vmag 20, or ZTF
never covered the positions that did:

`161P · 229P · 508P · 2P · 2014UN271 · 2019U5 · 2022R3 · 2023U1 · 2024G7 ·
2024T5 · 2024W1 · 2025O2 · 2025UX109`

`2P` is a special case: its 166 frames were queried earlier at Vmag < 21 and
are on disk, but the survey ran it at Vmag < 20 and found nothing.

## Aperture scale test

| ρ (km) | measured | skipped |
|---|---|---|
| 10,000 | 7,340 | 996 (11.9%) |
| 15,000 | 8,104 | 232 (2.8%) |
| 20,000 | 8,218 | 118 (1.4%) |
| 30,000 | 8,045 | 291 (3.5%) |
| 40,000 | 7,884 | 452 (5.4%) |

The 10,000 km aperture is the one that most often falls below the seeing disc;
40,000 km is the one that most often exceeds the 1′ ceiling.  20,000 km is the
best-sampled aperture and is the sensible default for cross-target comparison.

## Flag statistics

53.0% of all measurement rows are clean (`quality_ok`, no critical flag).

| CRITICAL flag | rows | share |
|---|---|---|
| contaminated | 7,970 | 20.1% |
| lowsnr | 6,167 | 15.6% |
| sky_edge | 4,881 | 12.3% |
| centroid | 3,754 | 9.5% |
| undersampled | 2,924 | 7.4% |
| negative_flux | 1,322 | 3.3% |
| aperture_edge | 1,207 | 3.0% |
| outside | 893 | 2.3% |

| ADVISORY flag | rows | share |
|---|---|---|
| color_default | 13,776 | 34.8% |
| nan_pixels | 678 | 1.7% |

Advisory flags qualify a measurement without invalidating it.  `color_default`
dominates simply because most frames are single-band, so no g−r colour is
available for the colour term.

Clean rate falls off at both aperture extremes — 53.8% at 10,000 km, peaking at
59.6% at 15,000 km, down to 42.0% at 40,000 km, where the aperture more often
reaches a frame edge or swallows a field star.

**Contamination is the largest single cause of rejection.**  One row in five
fails the Gaia DR3 test (G_eff within r_ap + FWHM brighter than 30% of the
comet's expected V).  That is the intended behaviour for a survey along the
ecliptic, but it means a fifth of the raw measurements are unusable rather than
merely noisy.

## Coma radial profiles

Available for 9 targets only (see the gap noted below): 942 frames, 524 clean.

| target | n | comet slope | star slope | excess at 3·FWHM | ρ_max (km) |
|---|---|---|---|---|---|
| 24P | 240 | −1.07 | −4.43 | 70.7 | 6,179 |
| 2025K1 | 83 | −1.60 | −4.55 | 53.0 | 7,696 |
| 2025L1 | 18 | −0.22 | −4.27 | 141.3 | 8,805 |
| 2025L2 | 15 | −2.49 | −4.12 | 16.2 | 20,320 |
| 2025M2 | 28 | −3.52 | −4.27 | 26.0 | 34,569 |
| 2025Q3 | 54 | −2.30 | −4.25 | 10.6 | 13,232 |
| 2025R1 | 36 | −2.49 | −4.16 | 18.8 | 13,593 |
| 2025R2 | 43 | −0.98 | −4.15 | 99.5 | 8,653 |
| 2025W2 | 7 | −1.81 | −4.57 | 28.9 | 7,039 |

Median comet slope **−1.23** against a median field-star slope of **−4.41**, and
**501 of 524** clean frames show the comet shallower than the stars in the same
image.  The comae are unambiguously extended, and the ensemble slope sits close
to the −1 expected for steady-state dust outflow.

The steep outliers (2025M2 at −3.52, 2025L2 and 2025R1 at −2.49) are the faint,
distant frames where the resolution study already identified sky subtraction as
the dominant bias, not a genuinely compact coma.

## Known gaps

1. **Radial profiles are missing for 46 of the 55 targets with data.**  The
   profile step was added to the package after the batch process had already
   started, and a running Python process keeps the modules it imported.  Only
   24P (run by hand) and the eight targets processed after the 20:01 restart
   have profiles.  Re-running `--steps profile figures` fills this in.
2. **2024E1 and 240P were queried over too short a window.**  `--end` was
   applied only when passed, so a stale per-target `end_date` in `config.py`
   won: 2024E1 stopped at 2025-10-30, missing its 2026-01-20 perihelion and the
   entire post-perihelion leg.  Fixed in `survey.py`; both need re-querying.
3. **Incomplete downloads.**  Roughly 300 frames failed, ~131 of them on 235P
   during a three-hour IRSA outage on the night of 09-07/08.  404s are permanent
   archive gaps; the 5xx failures are recoverable and IRSA has since recovered.
4. **Targets 1–12 were processed before the perihelion-axis change** and still
   carry signed-r_h Afρ plots.
