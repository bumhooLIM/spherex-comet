# ZTF 68-comet survey — process and flag summary

Survey 2026-09-07 17:25 → 2026-09-08 22:32.  Repair and re-reduction
2026-09-09.  Window 2025-03-01 → present, Vmag < 20, cutouts 5′ (10′ when
Vmag_min < 14 or Δ_min < 1 au).  Apertures ρ = 10/15/20/30/40 ×10³ km, each
kept only where `PSF_FWHM < r_ap < 60″`.

## Coverage

| | |
|---|---|
| designations processed | 68 / 68 |
| with photometry | 56 |
| no data | 12 |
| frames on disk | **8,730** |
| download completeness | **8,734 / 8,845 = 98.75%** |
| corrupt files | **0** |
| cutout PNGs (`fig/ztf/photometry/<T>/cutout/`) | 7,783 |
| measurement rows | 41,533 |
| clean rows | 51.2% (after the single-frame anomaly flag: 311 rows) |

### The 12 targets with no data

Three had **no epoch reaching Vmag 20** — the directory is empty by design:
`2014UN271` (Bernardinelli-Bernstein, ~17 au), `2022R3`, `229P`, `2025UX109`.

Nine had **epochs that passed the cut but no ZTF coverage** at those positions:
`161P` (17/112 epochs), `508P` (26), `2024W1` (26), `2025O2` (40), `2019U5`
(99), `2024G7` (109), `2023U1` (112), `2024T5` (112).

`2P` was recovered: its 166 frames were queried earlier at Vmag < 21, and the
repair pass adopted them, so it now carries full photometry despite the survey
scoring it `no_frames` at Vmag < 20.

### The 111 frames that could not be downloaded

Both kinds are permanent; a second repair pass recovered **zero** of them.

| cause | frames |
|---|---|
| HTTP 404 — not in the archive | 80 |
| HTTP 500 — `Cutout does not overlap image` | 31 |

The 500s are **not** server faults.  IRSA raises `ibe::HttpException` at
`stream_subimage.cxx:331` when the requested centre falls outside the image
footprint: the IBE metadata search matched the frame, but the pixel data does
not reach the comet.  Those frames would have failed `flag_outside` anyway.

## Aperture scale test

| ρ (km) | measured | skipped | clean |
|---|---|---|---|
| 10,000 | 7,731 | 999 (11.4%) | 53.6% |
| 15,000 | 8,498 | 232 (2.7%) | 58.8% |
| 20,000 | 8,612 | 118 (1.4%) | 57.9% |
| 30,000 | 8,439 | 291 (3.3%) | 49.2% |
| 40,000 | 8,253 | 477 (5.5%) | 40.8% |

10,000 km is the aperture most often below the seeing disc; 40,000 km most
often exceeds the 1′ ceiling or reaches a frame edge.  **20,000 km is the
best-sampled aperture** (98.6% of frames) and the right default for
cross-target comparison; 15,000 km has the highest clean fraction.

## Flag statistics

| CRITICAL flag | rows | share |
|---|---|---|
| contaminated | 8,279 | 19.9% |
| lowsnr | 7,012 | 16.9% |
| sky_edge | 5,207 | 12.5% |
| centroid | 3,867 | 9.3% |
| undersampled | 2,957 | 7.1% |
| negative_flux | 1,470 | 3.5% |
| aperture_edge | 1,268 | 3.1% |
| outside | 923 | 2.2% |

| ADVISORY flag | rows | share |
|---|---|---|
| color_default | 14,373 | 34.6% |
| nan_pixels | 696 | 1.7% |

Advisory flags qualify a measurement without invalidating it.  `color_default`
dominates because most frames are single-band, leaving no g−r for the colour
term.

**Contamination is the single largest cause of rejection**: one row in five
fails the Gaia DR3 test (G_eff within r_ap + FWHM brighter than 30% of the
comet's expected V).  For a survey along the ecliptic that is expected, but it
means a fifth of the raw measurements are unusable rather than merely noisy.

## Coma radial profiles

All targets, profiles centred on the optocentre (`refine_centre`, 52% of clean frames moved): **7,735 frames, 4,092 clean.**

| | |
|---|---|
| median comet slope | **−1.68** |
| median field−star slope | **−4.39** |
| comet shallower than stars | **96.5%** |

The comae are unambiguously extended, and the ensemble slope brackets the −1
expected for steady-state dust outflow.

| target | n | comet | stars | ρ_max (km) |
|---|---|---|---|---|
| 2025M2 | 28 | −3.52 | −4.27 | 34,569 |
| 2023V1 | 18 | −2.81 | −4.30 | 37,323 |
| 2024G4 | 157 | −2.70 | −4.51 | 32,765 |
| 63P | 23 | −2.59 | −4.67 | 13,712 |
| 491P | 144 | −2.54 | −4.58 | 23,057 |
| … | | | | |
| 24P | 240 | −1.07 | −4.43 | 6,179 |
| 2P | 54 | −1.07 | −4.40 | 23,453 |
| 210P | 97 | −1.06 | −4.55 | 9,203 |
| 2025R2 | 43 | −0.98 | −4.15 | 8,653 |
| 2025L1 | 18 | −0.22 | −4.27 | 8,805 |

The steep end is dominated by faint, distant targets, where the resolution
study (`doc/ztf/profile_resolution_24P.md`) showed sky subtraction — not a compact
coma — is the controlling bias.  Slopes near or below −2.5 should be read as
upper limits on steepness, not as measurements of coma structure.

## Caveats

- **2024E1 has mixed cutout sizes.**  Its re-query found Vmag_min = 13.62, so
  the adaptive rule now selects 10′; the 116 frames already on disk are 5′ and
  only the 2 new ones are 10′.  Harmless — the largest aperture spans ~32″ at
  Δ = 1.7 au, well inside a 300″ box — but worth knowing before comparing sky
  annuli across its frames.
- **2024E1's post-perihelion leg is two frames.**  Its window had been
  truncated at 2025-10-30 by a stale `end_date`; re-querying the full range
  took it from 49 to 112 passing epochs but added only 2 frames.  At
  q = 0.566 AU the comet was in solar conjunction from Palomar through
  perihelion (2026-01-20) and was caught only twice on the way out.
- **240P gained nothing** from the same fix: 112/112 epochs now pass instead of
  98/98, but ZTF observed none of the extra window.
