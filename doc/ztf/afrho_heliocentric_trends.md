# Heliocentric dependence of Afρ across the survey

Driver `scripts/ztf/afrho_trends.py`; science in `ztfcomet/activity.py`.  Tables
in `results/ztf/afrho/` (`trends.csv`, `peaks.csv`, `breaks.csv`, `colour.csv`,
`outbursts.csv`, `spherex_windows.csv`, `spherex_afrho.csv`); figures
`fig/ztf/afrho/trend/<comet>.png` (fits, peak, breaks, outbursts, colour, SPHEREx
windows, and in red the Afρ at each SPHEREx phase), `fig/ztf/afrho/trend_slide/`
(one figure per SPHEREx phase for the catalog's review deck), `fig/ztf/afrho/rh/`,
`fig/ztf/afrho/apertures/` and `fig/ztf/afrho/survey_overview.png`.  56 comets;
r-band is the main series and g is analysed separately; ρ = 10,000 and
20,000 km.  The aperture question raised by C/2024 E1 is answered in
`doc/ztf/afrho_aperture_systematics.md`.

## Method

- **Quantity**: phase-corrected A(0°)fρ.
- **Points**: `quality_ok` rows.  Two rejections were added for this
  analysis.  A **single frame far above its neighbours in time** (more than
  0.15 dex and 5× the series' robust scatter above the median of its
  nearest clean frames, with neither adjacent frame sharing half the
  excess) is flagged `anomalous_bright` in the photometry table itself and
  is no longer clean — 309 rows across 47
  comets, e.g. the 2022 E2 frame 50% above its own night.  Frames whose
  photometry centre sits more than half an aperture radius from the
  ephemeris are dropped here as before.  Survey clean rate is now
  51.4%.
- **Phases are defined by the peak of activity, not perihelion.**  A comet
  sampled on both sides of T_p gets a smoothed peak (Gaussian kernel,
  max(10 d, 2 × median gap)); when the maximum is *interior* to the coverage
  the series splits into *rising* and *fading* there — a plateau before the
  peak still splits (240P), and `bracketed` records whether the curve falls
  on both sides.  A maximum at the edge leaves one phase.  A one-sided comet
  is one leg named by its orbital direction.
- **Outbursts** are found first and set aside: a jump of ≥ 0.3 dex above
  the recent trend *extrapolated* to the frame (not the recent median — a
  steep smooth rise sits above its own median at every step), confirmed by
  the next frame, closing when the running median returns to the
  pre-outburst line.  A window must show a decay: a jump that keeps rising
  is an onset (261P) and stays in the trend.  Outburst frames are excluded
  from the primary fits and shaded in the figures.
- **Isolated points** at either end of the r_h range — at most max(4, 10%
  of N) frames beyond a gap of ≥ 0.1 dex — are set aside as a *tail* and
  fitted separately: a few frames far from the rest carry most of the
  leverage (10P's four frames at 3.6 au).
- **Fit**: `log Afρ = a − x log r_h`, weighted; the *scaled* error
  (× √χ²_red) is quoted.  Grades A–D as before; **grade-D fits are not
  drawn**, since a line through 27 points spanning 0.02 dex (2023 RS61) is
  arithmetic, not a trend.
- **Broken law**: legs with ≥ 10 points over ≥ 0.1 dex are fitted with a
  continuous two-slope law (free break) and with a break at 3 au, compared
  with the single law by ΔBIC on errors inflated so the single law has
  χ²_red = 1.  Where a break is preferred the two sides are refitted as
  plain power laws — the **divided trend** — and drawn as segments.
- **Smoothed trend**: a kernel smooth of log Afρ over log r_h through every
  quiet frame of every phase, drawn in black as a guide to the overall shape.
- **SPHEREx windows**: the JD and r_h ranges of each SPHEREx phase group, from
  the catalog's `phase_assignment.csv` and `apphot/`, are shaded on every
  panel (S1, S2, …): 165 windows over 56 comets (the catalog grouping of
  2026-09-14; `spherex_windows.py` is rebuilt by every full run).
- **Afρ at the SPHEREx epochs** (`activity.afrho_at_epoch`,
  `results/ztf/afrho/spherex_afrho.csv`, the red markers and the table strip of
  the trend figures): for each SPHEREx phase, at the mean r_h of the channels
  that carried its production rate (the exposure mean where no rate was
  fitted), r-band, both apertures.  Three estimates in order of what they
  assume — *direct*: at least two clean frames within 5 d of the window, each
  moved to ⟨r_h⟩ along the most local law of grade A–C and averaged with
  weights, the error that of the mean (scaled by χ²_ν) plus the slope term;
  *trend*: no frames then, but ⟨r_h⟩ lies inside a fitted law of the phase
  the epoch falls in (the rising/fading split at the peak, else
  inbound/outbound at T_p) — the most local of the primary law, its segments
  and its tail, any grade, since the *level* of a flat well-sampled series is
  measured even where its slope is not; the error carries the law's level and
  slope uncertainties and its RMS scatter, because the question is the
  comet's brightness at that epoch, not the mean law; *trend_extrap*: ⟨r_h⟩
  at most 0.1 dex beyond the primary law, grade A–C only.  Otherwise no value,
  and the note says why: no frames at that aperture, ZTF fitted only the other
  phase of the orbit, or ⟨r_h⟩ too far outside the data.  A single frame in
  the window is used only when no law can be evaluated (error floor 0.1 dex).
  Frames inside an outburst are used as they are — the value at the time is
  the question — and counted.
- **Colour**: same-night g/r pairs, excess over solar
  −2.5 log₁₀(Afρ_g / Afρ_r), with a two-level step as the change-point model.

## Reading the results

The generated section below is rebuilt by `scripts/ztf/afrho_trends_report.py`
after every run of `afrho_trends.py`.  What it shows, in brief:

- **Only four comets have both a rising and a fading phase with enough
  points to fit** (24P, 2023 R1, 2025 K1, 2024 G4); only 24P has both at
  grade B (r_h^−8.8 up, r_h^−4.0 down).  One-phase comets keep the earlier
  picture: outbound-only fading is uniform (median x ≈ 3, 16–84% within
  0.7–4.2), inbound-only indices scatter from −0.2 to 7.
- **The survivors of the broken-law test are regime changes** — 10P flat
  from 3.6 to 1.85 au then r_h^−10 (an onset); 2024 E1 an inbound maximum
  at 3.4 au; 2023 R1 flat beyond 3.9 au then steep — plus curvature within
  weeks of a peak (24P) and end-of-range features.  **3 au is not
  special**: only four comets straddle it and they disagree in sense.
- **Outbursts**: 217P's event at 3.10 au (×10, decaying over forty days)
  and a handful of smaller ones; the 217P outbound index without them is
  3.35 ± 0.09 (B).
- **Colour**: the dust is +0.10 mag redder than the Sun with no r_h trend;
  24P's change at 1.5 au doubles with aperture and is C₂ in g, not dust.
- **Sparse targets**: re-querying at V < 21 added 0–3 clean points at
  10,000 km, because the magnitude cut was never the limit: the distant
  ones (29P, 2023 RS61, 2022 R6, 2023 F3) fail the 10,000 km scale test on
  every frame and have 23–85 clean points at 20,000 km instead; the rest
  have little ZTF coverage.
- **Caveats**: a phase split at the peak can straddle perihelion; a 0.04 dex
  r_h range (240P) gives grade D on both sides of a real split; the anomaly
  rule cannot tell an artefact from a one-night brightening under sparse
  sampling (29P's flagged rows are its outbursts one frame at a time).

## Results (generated)

## Results (generated)

<!-- generated -->
### Headline numbers (r-band, primary phases, grades A/B)

| phase | 10,000 km: median x [16–84%] | 20,000 km |
|---|---|---|
| rising (to the peak) | 6.64 [3.91, 8.08] (n = 3) | 3.75 [3.02, 6.76] (n = 3) |
| fading (from the peak) | 3.86 [3.73, 3.99] (n = 2) | 3.62 [3.61, 3.63] (n = 2) |
| inbound only | 5.14 [3.19, 7.43] (n = 6) | 4.70 [0.80, 5.76] (n = 8) |
| outbound only | 3.18 [0.44, 4.22] (n = 11) | 2.64 [0.98, 3.77] (n = 12) |

Survey clean rate 51.2% after the anomaly flag (311 rows).  Of 14 two-sided comets, 8 have an interior maximum and are split there (6 bracketed on both sides).  11 primary legs at 10,000 km prefer a broken law and are divided: 10P, 145P, 2022E2, 2023A3, 2023R1, 2024E1, 2024J3, 2024L5, 2025A6, 2025K1, 2025M2, 2025Q3, 24P, 40P.  7 outburst windows on 5 comets are excluded from the fits; 8 comets have an isolated tail set aside.  SPHEREx: 165 windows over 56 comets.

### Activity peaks

| comet | ρ (km) | n rising / fading | peak T−T_p (d) | 16–84% | bracketed |
|---|---|---|---|---|---|
| 171P | 20k | 11 / 19 | **-22** | [-50, -15] | plateau on one side |
| 2022N2 | 10k | 13 / 48 | **+9** | [-5, +17] | plateau on one side |
| 2022N2 | 20k | 72 / 12 | **+168** | [+7, +169] | plateau on one side |
| 2022R6 | 20k | 4 / 49 | **-67** | [-97, +54] | plateau on one side |
| 2023R1 | 10k | 44 / 35 | **+17** | [+6, +23] | yes |
| 2023R1 | 20k | 80 / 31 | **+19** | [+13, +26] | yes |
| 2023V1 | 20k | 10 / 24 | **-56** | [-115, -55] | plateau on one side |
| 2024G4 | 10k | 30 / 27 | **-59** | [-238, -58] | yes |
| 2024G4 | 20k | 41 / 47 | **-84** | [-87, -8] | yes |
| 2024N1 | 20k | 25 / 8 | **+36** | [-156, +38] | plateau on one side |
| 2025K1 | 10k | 40 / 7 | **+38** | [-1, +40] | yes |
| 2025K1 | 20k | 28 / 7 | **+36** | [-3, +40] | yes |
| 235P | 10k | 8 / 71 | **-26** | [-90, -20] | yes |
| 235P | 20k | 11 / 63 | **-7** | [-28, -4] | yes |
| 240P | 10k | 25 / 23 | **-22** | [-129, -21] | plateau on one side |
| 240P | 20k | 22 / 21 | **-20** | [-22, -18] | yes |
| 24P | 10k | 27 / 66 | **+11** | [+6, +13] | yes |
| 24P | 20k | 20 / 53 | **+13** | [+7, +16] | yes |
| 47P | 10k | 38 / 17 | **-34** | [-35, -31] | yes |
| 47P | 20k | 41 / 20 | **-17** | [-19, -14] | yes |

### Where a single power law fails

| comet | ρ (km) | phase | r_h range | break (au) [16–84%] | x inside → outside | ΔBIC |
|---|---|---|---|---|---|---|
| 24P | 10k | rising | 1.18–2.05 | **1.24** [1.24, 1.31] | 21.9 ± 1.3 → 6.8 ± 0.2 | 14 |
| 2025K1 | 10k | rising | 0.69–2.82 | **1.32** [1.17, 1.50] | 3.5 ± 0.1 → 1.2 ± 0.2 | 15 |
| 24P | 20k | rising | 1.18–2.05 | **1.39** [1.39, 1.47] | 12.0 ± 0.5 → 4.4 ± 0.5 | 9 |
| 24P | 20k | fading | 1.24–2.39 | **1.47** [1.45, 1.65] | 5.4 ± 0.4 → 3.3 ± 0.1 | 7 |
| 2025A6 | 20k | inbound | 1.00–1.81 | **1.53** [1.50, 1.53] | 4.9 ± 0.1 → 5.7 ± 0.2 | 8 |
| 2025A6 | 10k | inbound | 1.00–1.81 | **1.53** [1.53, 1.56] | 4.6 ± 0.1 → 5.4 ± 0.2 | 10 |
| 10P | 20k | inbound | 1.45–2.63 | **1.59** [1.59, 2.07] | 11.2 ± 0.6 → 4.3 ± 0.2 | 17 |
| 10P | 10k | inbound | 1.44–2.63 | **1.85** [1.70, 1.89] | 9.4 ± 0.3 → -0.1 ± 0.4 | 31 |
| 145P | 10k | inbound | 1.89–2.48 | **2.05** [2.03, 2.16] | 11.8 ± 0.6 → 6.4 ± 0.4 | 11 |
| 145P | 20k | inbound | 1.89–2.48 | **2.08** [1.96, 2.16] | 12.3 ± 0.5 → 6.4 ± 0.6 | 13 |
| 40P | 10k | outbound | 1.83–2.87 | **2.41** [2.40, 2.45] | 5.5 ± 0.2 → 1.0 ± 0.5 | 22 |
| 40P | 20k | outbound | 1.83–2.87 | **2.45** [2.40, 2.54] | 5.2 ± 0.2 → 1.0 ± 0.6 | 10 |
| 2025Q3 | 10k | inbound | 2.13–2.99 | **2.66** [2.26, 2.68] | 8.6 ± 0.3 → -1.0 ± 1.1 | 12 |
| 2025Q3 | 20k | inbound | 2.13–2.99 | **2.66** [2.44, 2.66] | 10.1 ± 0.6 → -10.4 ± 1.6 | 12 |
| 2024E1 | 20k | inbound | 2.04–4.66 | **3.34** [3.34, 3.37] | -0.7 ± 0.1 → 1.0 ± 0.1 | 36 |
| 2024E1 | 10k | inbound | 2.01–4.66 | **3.37** [3.30, 3.49] | -0.7 ± 0.0 → 0.6 ± 0.1 | 29 |
| 2023R1 | 20k | rising | 3.57–4.85 | **3.87** [3.84, 3.88] | 13.6 ± 0.4 → -0.5 ± 0.2 | 59 |
| 2023R1 | 10k | rising | 3.57–4.60 | **3.91** [3.88, 3.91] | 15.4 ± 0.4 → -1.1 ± 0.4 | 31 |
| 2024L5 | 20k | outbound | 3.43–5.04 | **4.46** [3.46, 4.56] | 3.2 ± 0.1 → 7.1 ± 0.7 | 18 |
| 2024L5 | 10k | outbound | 3.43–5.04 | **4.47** [3.45, 4.55] | 3.7 ± 0.1 → 7.5 ± 0.7 | 13 |
| 2022E2 | 20k | outbound | 3.95–6.16 | **4.71** [4.71, 5.07] | 1.0 ± 0.2 → 3.0 ± 0.1 | 17 |
| 2023A3 | 10k | outbound | 3.03–5.18 | **4.82** [4.55, 4.82] | 2.3 ± 0.2 → 6.9 ± 1.1 | 10 |
| 2023A3 | 20k | outbound | 3.08–6.38 | **5.68** [5.18, 5.68] | 2.7 ± 0.1 → -5.2 ± 1.2 | 14 |
| 2024J3 | 20k | inbound | 4.41–6.03 | **5.74** [5.69, 5.76] | 0.3 ± 0.1 → 4.0 ± 0.7 | 8 |
| 2025M2 | 20k | inbound | 5.51–8.17 | **6.08** [6.08, 6.15] | -0.2 ± 0.5 → 3.8 ± 0.2 | 25 |

#### The 3 au hypothesis

| comet | ρ (km) | phase | x (< 3 au) | x (> 3 au) | split preferred |
|---|---|---|---|---|---|
| 2023C2 | 10k | outbound | -0.41 ± 0.36 | -0.53 ± 0.09 | no (ΔBIC -4) |
| 2023C2 | 20k | outbound | -0.64 ± 0.31 | -0.47 ± 0.07 | no (ΔBIC -4) |
| 2024E1 | 10k | inbound | -1.04 ± 0.07 | +0.33 ± 0.07 | yes (ΔBIC 27) |
| 2024E1 | 20k | inbound | -1.16 ± 0.09 | +0.82 ± 0.07 | yes (ΔBIC 35) |
| 217P | 10k | outbound | +3.37 ± 0.10 | -1.02 ± 3.88 | no (ΔBIC -3) |
| 43P | 10k | outbound | +4.27 ± 0.12 | +8.29 ± 1.58 | no (ΔBIC 1) |
| 43P | 20k | outbound | +3.90 ± 0.14 | +5.07 ± 1.73 | no (ΔBIC -4) |
| 48P | 10k | outbound | +2.77 ± 0.35 | +9.94 ± 3.22 | no (ΔBIC 1) |
| 48P | 20k | outbound | +2.63 ± 0.33 | +1.08 ± 2.69 | no (ΔBIC -3) |
| 78P | 10k | inbound | +2.74 ± 0.62 | +5.30 ± 0.57 | no (ΔBIC 1) |
| 78P | 20k | inbound | +3.46 ± 0.91 | +5.35 ± 0.94 | no (ΔBIC -2) |

### Outbursts

| comet | ρ (km) | r_h at onset (au) | T−T_p (d) | rise (dex) | frames | decayed within coverage |
|---|---|---|---|---|---|---|
| 10P | 10k | 2.21 | -177 … -149 | +0.40 | 6 | yes |
| 2023T3 | 20k | 5.69 | +514 … +531 | +0.38 | 8 | no |
| 2024G4 | 10k | 4.93 | +64 … +71 | +0.37 | 5 | yes |
| 217P | 10k | 1.80 | +109 … +110 | +0.37 | 4 | yes |
| 217P | 10k | 2.30 | +169 … +185 | +0.66 | 3 | yes |
| 217P | 10k | 3.10 | +273 … +317 | +0.98 | 26 | yes |
| 217P | 20k | 2.30 | +169 … +185 | +0.51 | 3 | yes |
| 217P | 20k | 3.10 | +273 … +318 | +0.71 | 28 | no |
| 235P | 10k | 1.98 | -3 … -0 | +1.25 | 3 | yes |
| 235P | 20k | 1.98 | -7 … -0 | +1.17 | 4 | yes |
| 29P | 20k | 6.27 | +2153 … +2172 | +0.43 | 5 | yes |
| 29P | 20k | 6.31 | +2496 … +2527 | +1.61 | 12 | yes |
| 40P | 10k | 1.95 | +70 … +93 | +1.01 | 16 | yes |
| 40P | 20k | 1.95 | +70 … +93 | +0.85 | 15 | yes |

### Dust colour

| comet | ρ (km) | pairs | r_h (au) | median excess (mag) | slope (mag/dex) | change at r_h | Δ colour |
|---|---|---|---|---|---|---|---|
| 10P | 10k | 23 | 1.44–3.58 | -0.051 | +0.49 ± 0.79 | none | — |
| 2023A3 | 10k | 12 | 3.35–4.97 | +0.116 | -0.49 ± 1.47 | none | — |
| 2023A3 | 20k | 24 | 3.35–5.49 | +0.121 | +0.36 ± 1.01 | none | — |
| 2023C2 | 10k | 22 | 3.12–4.96 | +0.092 | -0.15 ± 0.53 | none | — |
| 2023C2 | 20k | 21 | 3.10–4.96 | +0.115 | +0.15 ± 0.62 | none | — |
| 2023T3 | 20k | 9 | 3.69–5.78 | +0.054 | -0.68 ± 1.65 | none | — |
| 2024E1 | 10k | 24 | 2.01–4.32 | +0.090 | -0.05 ± 0.13 | none | — |
| 2024E1 | 20k | 41 | 2.04–4.66 | +0.098 | +0.09 ± 0.13 | none | — |
| 2024L5 | 10k | 29 | 3.43–5.04 | +0.131 | -0.62 ± 0.69 | none | — |
| 2024L5 | 20k | 52 | 3.43–5.04 | +0.110 | -0.48 ± 0.57 | none | — |
| 2025K1 | 10k | 22 | 0.87–2.70 | +0.067 | -0.18 ± 0.26 | none | — |
| 2025K1 | 20k | 15 | 0.87–2.66 | +0.060 | +0.01 ± 0.28 | none | — |
| 2025M2 | 20k | 41 | 5.51–8.16 | +0.112 | -0.72 ± 0.65 | none | — |
| 2025R2 | 10k | 15 | 1.09–2.60 | -0.051 | +1.45 ± 0.55 | none | — |
| 2025R2 | 20k | 13 | 1.28–2.62 | -0.177 | +1.59 ± 0.57 | none | — |
| 210P | 10k | 17 | 1.14–2.19 | -0.044 | -0.30 ± 0.42 | none | — |
| 217P | 10k | 39 | 2.03–3.42 | +0.181 | -0.29 ± 0.44 | none | — |
| 217P | 20k | 35 | 2.03–3.42 | +0.157 | -0.07 ± 0.46 | none | — |
| 24P | 10k | 63 | 1.19–2.39 | -0.082 | +1.52 ± 0.20 | **1.53** [1.49, 1.58] | **+0.24 ± 0.03** |
| 24P | 20k | 43 | 1.19–2.37 | -0.214 | +2.70 ± 0.21 | **1.50** [1.49, 1.50] | **+0.44 ± 0.04** |
| 40P | 10k | 20 | 1.95–2.83 | +0.202 | -0.27 ± 0.62 | none | — |
| 40P | 20k | 10 | 1.95–2.76 | +0.165 | -2.17 ± 1.14 | none | — |

### Per-comet indices, r-band, primary phases and segments

Scaled error, (N), grade; A/B in bold.  Segments (`r_h<` / `r_h>`) are the divided trend where a break is preferred.

| comet | r_h (au) | phase | 10,000 km | 20,000 km |
|---|---|---|---|---|
| 10P | 1.44–2.63 | inbound | **+5.57 ± 0.34 (50) B** | **+5.73 ± 0.38 (34) B** |
| 10P | 1.44–2.63 | inbound r_h<1.85 | +10.31 ± 0.31 (29) C | — |
| 10P | 1.44–2.63 | inbound r_h<1.59 | — | +11.36 ± 0.49 (11) D |
| 10P | 1.44–2.63 | inbound r_h>1.85 | **+0.66 ± 0.60 (21) B** | — |
| 10P | 1.44–2.63 | inbound r_h>1.59 | — | **+4.34 ± 0.33 (23) A** |
| 124P | 2.11–2.48 | inbound | +1.78 ± 1.43 (16) C | +3.02 ± 4.27 (6) D |
| 131P | 2.48–2.67 | inbound | +11.87 ± 6.63 (12) D | +21.24 ± 5.86 (4) D |
| 145P | 1.89–2.48 | inbound | **+8.62 ± 0.32 (38) B** | **+9.41 ± 0.36 (37) B** |
| 145P | 1.89–2.48 | inbound r_h<2.05 | +11.81 ± 1.02 (15) D | — |
| 145P | 1.89–2.48 | inbound r_h<2.08 | — | +12.35 ± 0.87 (18) D |
| 145P | 1.89–2.48 | inbound r_h>2.05 | +6.40 ± 0.44 (23) C | — |
| 145P | 1.89–2.48 | inbound r_h>2.08 | — | +6.45 ± 0.61 (19) C |
| 164P | 1.88–2.78 | outbound | **-0.62 ± 0.18 (75) B** | **-0.72 ± 0.19 (67) B** |
| 171P | 1.77–2.09 | rising | — | +33.42 ± 26.52 (11) D |
| 171P | 1.77–2.09 | fading | +9.02 ± 0.87 (40) C | +7.64 ± 1.00 (19) C |
| 172P | 3.40–3.51 | inbound | -1.59 ± 5.20 (18) D | +13.21 ± 14.97 (10) D |
| 2021G2 | 5.18–6.89 | outbound | **+1.07 ± 0.29 (17) B** | **+1.45 ± 0.06 (83) B** |
| 2022E2 | 3.95–6.16 | outbound | +2.25 ± 0.17 (36) C | **+2.29 ± 0.07 (87) B** |
| 2022E2 | 3.95–6.16 | outbound r_h<4.71 | — | +0.19 ± 0.64 (40) D |
| 2022E2 | 3.95–6.16 | outbound r_h>4.71 | — | **+2.91 ± 0.14 (47) B** |
| 2022N2 | 3.83–4.25 | rising | +3.14 ± 9.77 (13) D | +0.19 ± 0.27 (72) D |
| 2022N2 | 3.83–4.25 | fading | +0.17 ± 0.13 (48) D | +3.15 ± 1.47 (12) D |
| 2022QE78 | 5.48–5.64 | fading | — | -4.84 ± 0.77 (81) D |
| 2022QE78 | 5.48–5.64 | outbound | -6.06 ± 6.65 (5) D | — |
| 2022R6 | 6.58–6.74 | rising | — | +325.20 ± 414.97 (4) D |
| 2022R6 | 6.58–6.74 | fading | — | +1.74 ± 1.26 (49) D |
| 2023A3 | 3.03–6.38 | outbound | **+2.75 ± 0.15 (26) A** | **+2.46 ± 0.12 (45) B** |
| 2023A3 | 3.03–6.38 | outbound r_h<4.82 | **+2.30 ± 0.19 (21) A** | — |
| 2023A3 | 3.03–6.38 | outbound r_h<5.68 | — | **+2.66 ± 0.09 (38) A** |
| 2023A3 | 3.03–6.38 | outbound r_h>4.82 | +7.33 ± 5.63 (5) D | — |
| 2023A3 | 3.03–6.38 | outbound r_h>5.68 | — | -9.19 ± 3.65 (7) D |
| 2023C2 | 2.75–6.02 | outbound | **-0.51 ± 0.05 (56) A** | **-0.49 ± 0.05 (65) A** |
| 2023F3 | 5.71–6.10 | outbound | — | +1.96 ± 1.01 (23) D |
| 2023H5 | 4.32–4.43 | inbound | -6.81 ± 2.39 (8) D | -4.10 ± 1.54 (23) D |
| 2023R1 | 3.57–4.85 | rising | **+6.64 ± 0.35 (44) B** | **+3.75 ± 0.19 (80) B** |
| 2023R1 | 3.57–4.85 | rising r_h<3.91 | +15.50 ± 1.11 (17) D | — |
| 2023R1 | 3.57–4.85 | rising r_h<3.87 | — | +14.19 ± 1.29 (31) D |
| 2023R1 | 3.57–4.85 | rising r_h>3.91 | -1.04 ± 0.33 (27) C | — |
| 2023R1 | 3.57–4.85 | rising r_h>3.87 | — | -0.39 ± 0.13 (49) C |
| 2023R1 | 3.57–4.85 | fading | +15.89 ± 1.30 (35) D | +10.95 ± 1.27 (31) D |
| 2023RS61 | 8.54–8.97 | inbound | — | -43.20 ± 5.71 (27) D |
| 2023T3 | 3.68–5.67 | outbound | +3.23 ± 4.37 (6) D | +1.90 ± 0.31 (14) C |
| 2023V1 | 5.09–5.33 | rising | — | -4.70 ± 11.33 (10) D |
| 2023V1 | 5.09–5.33 | fading | +3.92 ± 4.90 (15) D | +7.74 ± 2.43 (24) D |
| 2024A1 | 3.92–4.37 | fading | +2.31 ± 0.54 (56) D | +2.19 ± 0.39 (84) D |
| 2024E1 | 2.01–4.66 | inbound | **-0.34 ± 0.05 (43) B** | **-0.06 ± 0.07 (50) B** |
| 2024E1 | 2.01–4.66 | inbound r_h<3.37 | **-0.74 ± 0.08 (17) A** | — |
| 2024E1 | 2.01–4.66 | inbound r_h<3.34 | — | **-0.67 ± 0.11 (16) A** |
| 2024E1 | 2.01–4.66 | inbound r_h>3.37 | **+0.59 ± 0.09 (26) B** | — |
| 2024E1 | 2.01–4.66 | inbound r_h>3.34 | — | **+1.05 ± 0.09 (34) B** |
| 2024G4 | 4.90–5.73 | rising | +1.17 ± 0.83 (30) D | +2.83 ± 0.68 (41) C |
| 2024G4 | 4.90–5.73 | fading | +5.17 ± 7.57 (27) D | -0.63 ± 15.32 (47) D |
| 2024J3 | 4.41–6.03 | inbound | +3.61 ± 1.18 (7) D | **+0.51 ± 0.07 (56) B** |
| 2024J3 | 4.41–6.03 | inbound r_h<5.74 | — | **+0.31 ± 0.10 (44) B** |
| 2024J3 | 4.41–6.03 | inbound r_h>5.74 | — | +4.19 ± 0.78 (12) D |
| 2024L5 | 3.43–5.04 | outbound | **+4.13 ± 0.11 (60) B** | **+3.71 ± 0.10 (90) B** |
| 2024L5 | 3.43–5.04 | outbound r_h<4.47 | **+3.54 ± 0.28 (27) B** | — |
| 2024L5 | 3.43–5.04 | outbound r_h<4.46 | — | **+3.22 ± 0.17 (48) B** |
| 2024L5 | 3.43–5.04 | outbound r_h>4.47 | +7.28 ± 0.71 (33) C | — |
| 2024L5 | 3.43–5.04 | outbound r_h>4.46 | — | +7.04 ± 0.82 (42) C |
| 2024N1 | 4.40–4.58 | rising | — | -10.36 ± 4.61 (25) D |
| 2024N1 | 4.40–4.58 | fading | -3.79 ± 4.97 (30) D | +90.70 ± 41.95 (8) D |
| 2025A6 | 1.00–1.81 | inbound | **+4.72 ± 0.05 (53) A** | **+5.02 ± 0.05 (52) A** |
| 2025A6 | 1.00–1.81 | inbound r_h<1.53 | **+4.56 ± 0.07 (35) B** | **+4.86 ± 0.07 (33) B** |
| 2025A6 | 1.00–1.81 | inbound r_h>1.53 | +5.53 ± 0.25 (18) C | +5.65 ± 0.23 (19) C |
| 2025K1 | 0.69–2.82 | rising | **+2.62 ± 0.09 (40) B** | **+2.68 ± 0.07 (28) B** |
| 2025K1 | 0.69–2.82 | rising r_h<1.32 | **+3.43 ± 0.49 (17) B** | — |
| 2025K1 | 0.69–2.82 | rising r_h>1.32 | **+1.21 ± 0.16 (23) B** | — |
| 2025K1 | 0.69–2.82 | fading | **+3.67 ± 0.48 (7) B** | **+3.63 ± 0.60 (7) B** |
| 2025L1 | 1.71–1.89 | outbound | +30.94 ± 2.43 (7) D | +15.27 ± 4.61 (8) D |
| 2025L2 | 3.22–3.54 | inbound | +2.63 ± 1.51 (9) D | +3.31 ± 2.07 (8) D |
| 2025M2 | 5.51–8.17 | inbound | +0.92 ± 2.11 (11) D | **+2.94 ± 0.14 (74) B** |
| 2025M2 | 5.51–8.17 | inbound r_h<6.08 | — | +1.27 ± 0.93 (28) D |
| 2025M2 | 5.51–8.17 | inbound r_h>6.08 | — | **+3.93 ± 0.19 (46) B** |
| 2025Q3 | 2.13–2.99 | inbound | **+7.14 ± 0.47 (27) B** | **+5.76 ± 0.97 (23) B** |
| 2025Q3 | 2.13–2.99 | inbound r_h<2.66 | +8.66 ± 0.35 (20) C | +10.15 ± 0.44 (19) C |
| 2025Q3 | 2.13–2.99 | inbound r_h>2.66 | -0.81 ± 1.68 (7) D | -8.62 ± 8.16 (4) D |
| 2025R1 | 1.98–2.41 | fading | -0.63 ± 0.71 (18) C | — |
| 2025R1 | 1.98–2.41 | inbound | — | -0.11 ± 1.31 (9) C |
| 2025R2 | 1.68–2.62 | outbound | **+2.75 ± 0.25 (16) B** | **+2.77 ± 0.18 (15) B** |
| 2025W2 | 1.46–1.59 | inbound | -4.79 ± 10.72 (3) D | -21.69 ± 6.63 (4) D |
| 210P | 0.61–2.19 | outbound | **+3.33 ± 0.04 (56) B** | **+3.07 ± 0.04 (44) B** |
| 217P | 1.54–3.42 | outbound | **+3.35 ± 0.09 (87) B** | **+3.30 ± 0.18 (74) B** |
| 235P | 1.98–2.48 | rising | +8.87 ± 5.54 (8) D | +10.07 ± 7.26 (11) D |
| 235P | 1.98–2.48 | fading | +7.46 ± 0.15 (71) C | +7.22 ± 0.19 (63) C |
| 240P | 2.12–2.34 | rising | +0.25 ± 0.64 (25) D | -0.28 ± 0.61 (22) D |
| 240P | 2.12–2.34 | fading | +20.13 ± 2.82 (23) D | +18.93 ± 3.66 (21) D |
| 24P | 1.18–2.39 | rising | **+8.76 ± 0.33 (27) B** | **+8.18 ± 0.72 (20) B** |
| 24P | 1.18–2.39 | rising r_h<1.24 | +26.79 ± 8.74 (6) D | — |
| 24P | 1.18–2.39 | rising r_h<1.39 | — | +11.91 ± 0.52 (9) C |
| 24P | 1.18–2.39 | rising r_h>1.24 | **+6.83 ± 0.21 (21) A** | — |
| 24P | 1.18–2.39 | rising r_h>1.39 | — | **+4.32 ± 0.66 (11) B** |
| 24P | 1.18–2.39 | fading | **+4.05 ± 0.10 (66) B** | **+3.60 ± 0.09 (53) B** |
| 24P | 1.18–2.39 | fading r_h<1.47 | — | +5.44 ± 0.38 (11) C |
| 24P | 1.18–2.39 | fading r_h>1.47 | — | **+3.28 ± 0.13 (42) B** |
| 261P | 2.01–2.42 | inbound | +15.72 ± 0.55 (56) C | +16.35 ± 0.99 (30) D |
| 29P | 6.27–6.31 | outbound | — | +13.90 ± 22.86 (54) D |
| 2P | 4.02–4.09 | inbound | -5.92 ± 12.77 (25) D | -17.39 ± 27.15 (8) D |
| 302P | 3.32–3.87 | outbound | +2.41 ± 0.61 (40) C | +1.58 ± 0.54 (44) C |
| 306P | 1.29–1.48 | outbound | -0.24 ± 3.35 (5) D | — |
| 40P | 1.83–2.87 | outbound | **+4.34 ± 0.14 (63) B** | **+4.13 ± 0.14 (55) B** |
| 40P | 1.83–2.87 | outbound r_h<2.41 | **+5.50 ± 0.18 (42) B** | — |
| 40P | 1.83–2.87 | outbound r_h<2.45 | — | **+5.12 ± 0.26 (42) B** |
| 40P | 1.83–2.87 | outbound r_h>2.41 | +0.96 ± 0.86 (21) C | — |
| 40P | 1.83–2.87 | outbound r_h>2.45 | — | +0.02 ± 1.29 (13) C |
| 43P | 2.45–3.10 | outbound | **+4.40 ± 0.10 (95) B** | **+3.95 ± 0.12 (88) B** |
| 47P | 2.81–3.14 | rising | +11.69 ± 0.94 (38) D | +16.13 ± 1.40 (41) D |
| 47P | 2.81–3.14 | fading | +24.61 ± 7.31 (17) D | +2.55 ± 1.65 (20) D |
| 486P | 2.33–2.86 | outbound | +0.49 ± 0.41 (57) C | -0.03 ± 0.48 (55) C |
| 48P | 2.37–3.15 | outbound | **+3.18 ± 0.29 (34) B** | **+2.52 ± 0.25 (34) B** |
| 491P | 3.91–4.64 | outbound | -3.41 ± 0.72 (62) C | -4.07 ± 0.49 (69) C |
| 493P | 3.82–4.26 | rising | +6.22 ± 1.35 (39) D | +6.14 ± 1.84 (60) D |
| 499P | 1.33–1.76 | outbound | -1.88 ± 2.37 (6) D | — |
| 63P | 2.45–2.76 | inbound | +6.72 ± 2.53 (10) D | +9.05 ± 3.13 (8) D |
| 78P | 2.58–3.54 | inbound | **+4.08 ± 0.27 (50) B** | **+4.38 ± 0.40 (42) B** |

### Afρ at the SPHEREx epochs

Per SPHEREx phase group, the ZTF r-band A(0°)fρ at the group's mean r_h (`spherex_afrho.csv`; `activity.afrho_at_epoch`).  `direct`: clean frames within 5 d of the window, moved to ⟨r_h⟩ along the local law and averaged; `trend`: the fitted law of the orbital phase the epoch falls in, at ⟨r_h⟩; `trend_extrap`: that law extended by at most 0.1 dex; `none`: no estimate, with the reason.  Errors are 1σ in log space (the scatter about the law included for `trend`); the linear error quoted is the symmetric approximation, the 1σ range is in the table.

| ρ (km) | phases | direct | trend | trend_extrap | none | with a value |
|---|---|---|---|---|---|---|
| 10k | 165 | 79 | 25 | 12 | 49 | 116 |
| 20k | 165 | 79 | 33 | 8 | 45 | 120 |

**Per-epoch overrides (2026-09-15 memo, `ztfcomet.config.AFRHO_EPOCH_OVERRIDES`; applied by
`afrho_trends.py` through keyword arguments of `afrho_at_epoch`, and named in the `note`).**
47P S1: the rising law is grade D (a 0.09 au baseline) and S1 lies 0.005 dex beyond it — the
extrapolation is allowed for this epoch (`extrap_grades` includes D).  210P S1–S4: ZTF observed
the comet after perihelion only; the four inbound SPHEREx phases are read off the outbound law
(`leg="outbound"`; S1 at 2.47 au is a 0.05 dex extrapolation of its 0.61–2.19 au range), a
symmetric-activity assumption the note states.  217P S1: two of the five in-window frames sat
at the start of the 2.30 au outburst and lifted the direct mean to 124 cm; they are excluded
(`exclude_outburst=True`) and the quiescent frames give the value.  The tables below and the
per-phase figures were regenerated for these three comets and for 240P (regrouped) after the
override; `doc/comspec/case_revisions.md` lists the before → after values.

Why there is no estimate:

| reason | 10,000 km | 20,000 km |
|---|---|---|
| the law is grade D and ⟨r_h⟩ lies outside its data | 19 | 20 |
| ZTF fitted only the other phase of the orbit | 17 | 15 |
| no clean ZTF r-band frames at this aperture | 6 | 0 |
| no fitted trend at this aperture | 4 | 7 |
| ⟨r_h⟩ more than 0.1 dex beyond the fitted range | 3 | 3 |

Per phase (value ± error in cm, then how: method, points, grade of the law used):

| comet | phase | ⟨r_h⟩ (au) | T−T_p (d) | 10,000 km | how | 20,000 km | how |
|---|---|---|---|---|---|---|---|
| 10P | S1 | 3.63 | -442 | 32.5 ± 4.4 | direct (4) B | 20.9 ± 8.8 | direct (1) |
| 10P | S2 | 2.03 | -147 | 31.0 ± 11 | trend (21) B | 31.2 ± 12 | trend (23) A |
| 10P | S3 | 1.91 | -128 | 32.5 ± 4.8 | direct (4) B | 40.5 ± 16 | trend (23) A |
| 10P | S4 | 1.66 | -84 | 80.1 ± 4.1 | direct (4) C | 73.7 ± 4.6 | direct (2) A |
| 124P | S1 | 2.67 | -247 | 7.02 ± 2.7 | trend_extrap (16) C | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 124P | S2 | 2.54 | -224 | 7.63 ± 2.6 | trend_extrap (16) C | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 124P | S3 | 1.89 | -83 | 13.0 ± 4.8 | trend_extrap (16) C | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 131P | S1 | 2.88 | -255 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 131P | S2 | 2.46 | -78 | 15.4 ± 5.1 | direct (1) | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 145P | S1 | 2.90 | -262 | 2.33 ± 0.45 | trend_extrap (38) B | 1.32 ± 0.28 | trend_extrap (37) B |
| 145P | S2 | 2.84 | -250 | 2.84 ± 0.53 | trend_extrap (38) B | 1.64 ± 0.34 | trend_extrap (37) B |
| 145P | S3 | 2.03 | -81 | 45.0 ± 1.5 | direct (9) B | 34.2 ± 1.3 | direct (10) B |
| 164P | S1 | 2.40 | +189 | 77.0 ± 7.7 | direct (4) B | 54.9 ± 5.3 | direct (5) B |
| 164P | S2 | 3.26 | +341 | 82.3 ± 16 | trend_extrap (75) B | 62.1 ± 13 | trend_extrap (67) B |
| 171P | S1 | 1.84 | +55 | 17.2 ± 2.8 | direct (4) C | 16.8 ± 4.4 | direct (3) C |
| 171P | S2 | 2.55 | +217 | 0.748 ± 0.29 | trend_extrap (40) C | 1.04 ± 0.41 | trend_extrap (19) C |
| 171P | S3 | 2.66 | +238 | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range |
| 172P | S1 | 3.41 | -138 | 31.1 ± 4.7 | direct (3) | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 172P | S2 | 3.39 | -106 | 30.3 ± 6.1 | trend (18) D | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 172P | S3 | 3.42 | +158 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2021G2 | S1 | 5.49 | +301 | 4,292 ± 337 | trend (17) B | 4,197 ± 129 | direct (2) B |
| 2021G2 | S2 | 6.46 | +546 | 3,609 ± 338 | trend (17) B | 3,351 ± 105 | direct (4) B |
| 2021G2 | S3 | 6.69 | +595 | 3,459 ± 752 | direct (2) B | 3,157 ± 120 | direct (5) B |
| 2022E2 | S1 | 4.64 | +324 | 2,055 ± 141 | direct (2) C | 2,098 ± 118 | direct (2) B |
| 2022E2 | S2 | 5.26 | +436 | 1,576 ± 113 | trend (36) C | 1,532 ± 120 | trend (47) B |
| 2022N2 | S1 | 3.83 | +17 | 1,672 ± 16 | direct (6) | 1,589 ± 31 | direct (5) |
| 2022N2 | S2 | 4.17 | +189 | 1,625 ± 17 | direct (9) | 1,533 ± 20 | direct (9) |
| 2022QE78 | S1 | 5.49 | +56 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 2,925 ± 75 | direct (4) |
| 2022QE78 | S2 | 5.70 | +218 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2022R6 | S1 | 6.59 | +86 | — | no clean ZTF r-band frames at this aperture | 1,999 ± 46 | direct (12) |
| 2022R6 | S2 | 6.75 | +236 | — | no clean ZTF r-band frames at this aperture | 1,960 ± 142 | trend (49) D |
| 2023A3 | S1 | 3.65 | +219 | 1,144 ± 93 | trend (21) A | 957 ± 136 | direct (2) A |
| 2023A3 | S2 | 4.83 | +324 | 596 ± 60 | direct (3) A | 511 ± 34 | direct (4) A |
| 2023A3 | S3 | 5.29 | +368 | 426 ± 49 | trend_extrap (26) A | 387 ± 74 | trend (38) A |
| 2023A3 | S4 | 7.14 | +564 | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range | 201 ± 72 | trend_extrap (45) B |
| 2023C2 | S1 | 3.45 | +233 | 916 ± 32 | direct (5) A | 783 ± 34 | direct (5) A |
| 2023C2 | S2 | 3.75 | +274 | 1,026 ± 59 | direct (3) A | 798 ± 52 | direct (3) A |
| 2023C2 | S3 | 4.56 | +378 | 1,201 ± 64 | direct (3) A | 945 ± 38 | direct (6) A |
| 2023F3 | S1 | 5.26 | +114 | — | no fitted trend at this aperture | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2023F3 | S2 | 5.72 | +322 | — | no fitted trend at this aperture | 762 ± 32 | direct (3) |
| 2023F3 | S3 | 6.21 | +458 | — | no fitted trend at this aperture | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2023H5 | S1 | 4.55 | +177 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2023R1 | S1 | 4.51 | -309 | 1,136 ± 66 | direct (2) C | 1,052 ± 22 | direct (7) C |
| 2023R1 | S2 | 4.00 | -200 | 1,132 ± 96 | direct (3) C | 979 ± 59 | trend (49) C |
| 2023R1 | S3 | 3.57 | -5 | 4,015 ± 308 | trend (17) D | 2,886 ± 308 | trend (31) D |
| 2023R1 | S4 | 3.57 | +6 | 4,163 ± 63 | direct (15) B | 3,224 ± 49 | direct (14) B |
| 2023RS61 | S1 | 8.97 | -1207 | — | no clean ZTF r-band frames at this aperture | 4,000 ± 145 | direct (10) |
| 2023RS61 | S2 | 8.74 | -1036 | — | no clean ZTF r-band frames at this aperture | 1,181 ± 111 | direct (3) |
| 2023T3 | S1 | 3.79 | +146 | 475 ± 20 | direct (2) | 398 ± 17 | direct (2) C |
| 2023T3 | S2 | 5.13 | +422 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 226 ± 76 | trend (14) C |
| 2023V1 | S1 | 5.10 | -40 | 357 ± 15 | direct (9) | 310 ± 11 | direct (17) |
| 2023V1 | S2 | 5.33 | +208 | 248 ± 77 | direct (1) | 222 ± 26 | direct (3) |
| 2024A1 | S1 | 4.05 | +137 | 484 ± 24 | direct (10) | 417 ± 18 | direct (10) |
| 2024A1 | S2 | 4.52 | +271 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2024E1 | S1 | 3.70 | -234 | 804 ± 8.3 | direct (15) B | 678 ± 7.8 | direct (19) B |
| 2024E1 | S2 | 3.18 | -191 | 813 ± 14 | direct (6) A | 730 ± 15 | direct (8) A |
| 2024G4 | S1 | 5.22 | -232 | 320 ± 27 | direct (4) | 241 ± 19 | direct (5) C |
| 2024G4 | S2 | 4.91 | +34 | 329 ± 53 | direct (3) | 245 ± 17 | direct (9) |
| 2024J2 | S1 | 2.46 | +135 | 82.8 ± 20 | direct (1) | 81.2 ± 20 | direct (1) |
| 2024J2 | S2 | 3.52 | +258 | — | no fitted trend at this aperture | — | no fitted trend at this aperture |
| 2024J3 | S1 | 5.38 | -442 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 2,370 ± 177 | trend (44) B |
| 2024J3 | S2 | 4.29 | -214 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 2,570 ± 233 | trend_extrap (56) B |
| 2024L5 | S1 | 3.50 | +72 | 486 ± 41 | direct (4) B | 367 ± 33 | direct (4) B |
| 2024L5 | S2 | 3.54 | +91 | 548 ± 72 | trend (27) B | 414 ± 45 | trend (48) B |
| 2024L5 | S3 | 4.30 | +279 | 274 ± 32 | direct (2) B | 230 ± 13 | direct (7) B |
| 2024L5 | S4 | 5.03 | +405 | 96.5 ± 12 | trend (33) C | 82.9 ± 9.7 | trend (42) C |
| 2024N1 | S1 | 4.44 | -73 | 84.4 ± 4.3 | direct (24) | 63.9 ± 3.7 | direct (22) |
| 2024N1 | S2 | 4.41 | -33 | 84.4 ± 14 | trend (30) D | 59.7 ± 15 | trend (25) D |
| 2024N1 | S3 | 4.50 | +114 | 104 ± 14 | direct (4) | 104 ± 17 | direct (8) |
| 2025A6 | S1 | 2.14 | +112 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025A6 | S2 | 2.47 | +136 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025A6 | S3 | 3.23 | +194 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025K1 | S1 | 2.49 | -127 | 93.1 ± 8.9 | direct (5) B | 56.3 ± 5.9 | direct (4) B |
| 2025K1 | S2 | 2.26 | -111 | 133 ± 9.1 | direct (3) B | 88.1 ± 5.9 | direct (3) B |
| 2025K1 | S3 | 1.36 | -56 | 246 ± 72 | trend (23) B | 314 ± 90 | trend (28) B |
| 2025K1 | S4 | 1.08 | +41 | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range |
| 2025K1 | S5 | 1.74 | +78 | 49.9 ± 3.3 | direct (2) B | 60.0 ± 4.0 | direct (2) B |
| 2025K1 | S6 | 1.83 | +84 | 40.7 ± 6.8 | trend (7) B | 49.4 ± 10 | trend (7) B |
| 2025K1 | S7 | 1.99 | +94 | 29.5 ± 2.3 | direct (5) B | 36.2 ± 3.4 | direct (5) B |
| 2025L1 | S1 | 3.06 | -206 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025L1 | S2 | 2.79 | -177 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025L1 | S3 | 2.63 | -160 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 2025L1 | S4 | 1.77 | +42 | 13.4 ± 1.2 | direct (3) | 20.5 ± 4.0 | direct (2) |
| 2025L1 | S5 | 2.18 | +108 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2025L1 | S6 | 2.33 | +126 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2025L2 | S1 | 3.53 | -211 | 46.2 ± 7.0 | trend (9) D | 28.5 ± 6.2 | trend (8) D |
| 2025L2 | S2 | 2.96 | -85 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 2025M2 | S1 | 8.54 | -923 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 315 ± 46 | trend_extrap (74) B |
| 2025M2 | S2 | 7.49 | -773 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 458 ± 70 | trend (46) B |
| 2025M2 | S3 | 6.02 | -569 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 882 ± 123 | trend (74) B |
| 2025Q3 | S1 | 2.89 | -173 | 37.1 ± 8.2 | direct (4) B | 48.0 ± 16 | direct (3) B |
| 2025Q3 | S2 | 2.69 | -145 | 31.6 ± 3.3 | direct (3) B | 21.6 ± 16 | trend (4) D |
| 2025Q3 | S3 | 2.18 | -46 | 176 ± 8.2 | direct (3) C | 125 ± 5.8 | direct (4) C |
| 2025R1 | S1 | 2.17 | -75 | 26.1 ± 3.2 | direct (3) C | 17.9 ± 2.9 | trend (9) C |
| 2025R1 | S2 | 2.02 | -33 | 29.1 ± 1.3 | direct (6) C | 18.2 ± 1.4 | direct (5) C |
| 2025R1 | S3 | 1.99 | -21 | 26.3 ± 4.1 | direct (2) C | 17.6 ± 2.8 | direct (2) C |
| 2025R1 | S4 | 2.56 | +137 | 33.8 ± 7.4 | trend_extrap (18) C | — | ZTF fitted only the other phase of the orbit |
| 2025R1 | S5 | 2.76 | +164 | 35.5 ± 9.3 | trend_extrap (18) C | — | ZTF fitted only the other phase of the orbit |
| 2025R2 | S1 | 1.10 | +45 | 833 ± 77 | trend (3) D | — | ⟨r_h⟩ more than 0.1 dex beyond the fitted range |
| 2025R2 | S2 | 2.00 | +101 | 94.5 ± 2.9 | direct (3) B | 112 ± 7.9 | trend (15) B |
| 2025R2 | S3 | 2.23 | +117 | 70.8 ± 4.6 | trend (16) B | 82.9 ± 6.1 | trend (15) B |
| 2025R2 | S4 | 2.33 | +124 | 67.6 ± 7.8 | direct (2) B | 77.5 ± 4.0 | direct (2) B |
| 2025W2 | S1 | 1.57 | -44 | — | the law is grade D and ⟨r_h⟩ lies outside its data | 23.3 ± 12 | trend (4) D |
| 2025W2 | S2 | 1.53 | -34 | 7.27 ± 2.8 | direct (1) | 15.0 ± 9.2 | direct (3) |
| 2025W2 | S3 | 1.50 | -25 | 6.75 ± 1.7 | trend (3) D | 5.60 ± 1.1 | direct (2) |
| 2025W2 | S4 | 1.46 | +9 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 210P | S1 | 2.47 | -163 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 210P | S2 | 1.38 | -70 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 210P | S3 | 1.12 | -52 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 210P | S4 | 1.08 | -49 | — | ZTF fitted only the other phase of the orbit | — | ZTF fitted only the other phase of the orbit |
| 210P | S5 | 1.54 | +82 | 32.4 ± 4.0 | trend (56) B | 32.9 ± 4.3 | trend (44) B |
| 210P | S6 | 1.67 | +92 | 24.2 ± 1.2 | direct (6) B | 26.3 ± 2.2 | direct (6) B |
| 210P | S7 | 1.88 | +109 | 16.2 ± 2.0 | direct (2) B | 17.8 ± 2.4 | trend (44) B |
| 217P | S1 | 2.26 | +165 | 124 ± 25 | direct (5) B | 111 ± 21 | direct (4) B |
| 217P | S2 | 3.32 | +304 | 51.7 ± 3.8 | direct (5) B | 47.6 ± 3.8 | direct (4) B |
| 217P | S3 | 3.42 | +319 | 37.2 ± 3.8 | direct (4) B | 37.4 ± 3.6 | direct (5) B |
| 235P | S1 | 2.06 | +73 | 455 ± 7.2 | direct (5) C | 369 ± 6.0 | direct (4) C |
| 235P | S2 | 2.15 | +109 | 356 ± 37 | trend (71) C | 305 ± 35 | trend (63) C |
| 240P | S1 | 2.39 | -143 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 240P | S2 | 2.35 | -129 | 442 ± 18 | direct (7) | 364 ± 14 | direct (7) |
| 240P | S3 | 2.26 | -98 | 393 ± 13 | direct (6) | 322 ± 11 | direct (4) |
| 240P | S4 | 2.16 | +52 | 216 ± 7.0 | direct (9) | 192 ± 7.0 | direct (6) |
| 24P | S1 | 1.91 | -121 | 8.15 ± 1.5 | direct (2) A | 12.0 ± 3.9 | trend (11) B |
| 24P | S2 | 1.82 | -111 | 10.6 ± 1.2 | direct (2) A | 14.8 ± 4.6 | trend (11) B |
| 24P | S3 | 1.74 | -102 | 12.9 ± 0.69 | direct (5) A | 13.1 ± 0.88 | direct (3) B |
| 24P | S4 | 1.58 | -82 | 24.4 ± 3.8 | trend (21) A | 26.7 ± 7.5 | trend (11) B |
| 24P | S5 | 1.34 | -49 | 74.9 ± 12 | trend (21) A | 72.5 ± 10 | trend (9) C |
| 24P | S6 | 1.21 | -20 | 188 ± 8.2 | direct (5) B | 208 ± 11 | direct (4) C |
| 24P | S7 | 1.19 | -7 | 293 ± 19 | direct (4) B | 296 ± 23 | direct (3) C |
| 24P | S8 | 1.20 | +12 | 356 ± 8.9 | direct (9) B | 366 ± 41 | direct (7) C |
| 261P | S1 | 2.53 | -204 | 2.82 ± 0.76 | trend_extrap (56) C | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 261P | S2 | 2.14 | -93 | 29.6 ± 1.0 | direct (10) C | 21.5 ± 1.1 | direct (6) |
| 261P | S3 | 2.05 | -50 | 56.4 ± 4.4 | direct (4) C | 40.9 ± 3.9 | direct (3) |
| 29P | S1 | 6.30 | +2429 | — | no clean ZTF r-band frames at this aperture | 2,621 ± 1,205 | trend (54) D |
| 29P | S2 | 6.31 | +2579 | — | no clean ZTF r-band frames at this aperture | 2,738 ± 410 | direct (2) |
| 2P | S1 | 4.10 | -593 | 31.0 ± 5.7 | trend (25) D | 21.4 ± 7.9 | trend (8) D |
| 2P | S2 | 3.93 | -456 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 302P | S1 | 3.33 | +121 | 238 ± 5.2 | direct (5) C | 185 ± 3.8 | direct (6) C |
| 302P | S2 | 3.52 | +284 | 204 ± 6.2 | direct (9) C | 171 ± 5.0 | direct (9) C |
| 306P | S1 | 1.63 | -90 | — | ZTF fitted only the other phase of the orbit | — | no fitted trend at this aperture |
| 306P | S2 | 1.52 | -73 | — | ZTF fitted only the other phase of the orbit | — | no fitted trend at this aperture |
| 306P | S3 | 1.42 | -55 | — | ZTF fitted only the other phase of the orbit | — | no fitted trend at this aperture |
| 306P | S4 | 1.30 | -17 | — | ZTF fitted only the other phase of the orbit | 2.87 ± 0.88 | direct (1) |
| 306P | S5 | 1.28 | +7 | 2.63 ± 0.66 | direct (2) | 2.87 ± 0.88 | direct (1) |
| 40P | S1 | 2.23 | +130 | 52.9 ± 2.5 | direct (5) B | 39.1 ± 2.1 | direct (4) B |
| 43P | S1 | 2.61 | +128 | 126 ± 5.2 | direct (6) B | 98.6 ± 4.4 | direct (6) B |
| 43P | S2 | 3.13 | +287 | 54.8 ± 4.7 | trend_extrap (95) B | 47.8 ± 4.1 | trend_extrap (88) B |
| 47P | S1 | 2.93 | -149 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 47P | S2 | 2.81 | +22 | 167 ± 4.3 | direct (9) | 136 ± 4.7 | direct (8) |
| 486P | S1 | 2.39 | +92 | 36.9 ± 1.3 | direct (10) C | 26.8 ± 1.2 | direct (9) C |
| 486P | S2 | 2.84 | +258 | 27.4 ± 3.0 | direct (3) C | 26.3 ± 3.5 | direct (3) C |
| 48P | S1 | 2.38 | +168 | 107 ± 3.4 | direct (10) B | 106 ± 3.2 | direct (10) B |
| 48P | S2 | 3.10 | +337 | 39.4 ± 4.0 | direct (5) B | 56.6 ± 6.9 | direct (6) B |
| 491P | S1 | 4.08 | +383 | 165 ± 7.9 | direct (6) C | 115 ± 7.5 | direct (5) C |
| 493P | S1 | 3.90 | -128 | 252 ± 15 | direct (7) | 186 ± 14 | direct (6) |
| 493P | S2 | 3.83 | +34 | 290 ± 34 | trend (39) D | 260 ± 28 | direct (2) |
| 499P | S1 | 1.33 | +68 | 2.40 ± 1.5 | trend (6) D | — | no fitted trend at this aperture |
| 499P | S2 | 1.43 | +78 | 2.54 ± 0.61 | direct (4) | 3.35 ± 1.7 | direct (2) |
| 499P | S3 | 1.58 | +94 | 2.07 ± 0.35 | direct (3) | — | no fitted trend at this aperture |
| 499P | S4 | 1.71 | +109 | 4.69 ± 2.0 | direct (2) | — | no fitted trend at this aperture |
| 63P | S1 | 3.25 | -286 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 63P | S2 | 2.25 | -112 | — | the law is grade D and ⟨r_h⟩ lies outside its data | — | the law is grade D and ⟨r_h⟩ lies outside its data |
| 78P | S1 | 3.54 | -418 | 58.7 ± 9.6 | trend (50) B | 40.4 ± 10 | trend (42) B |
| 78P | S2 | 3.38 | -383 | 70.3 ± 11 | trend (50) B | 49.1 ± 12 | trend (42) B |
| 78P | S3 | 2.87 | -271 | 130 ± 3.2 | direct (15) B | 93.1 ± 3.0 | direct (12) B |
<!-- /generated -->
