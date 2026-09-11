# Heliocentric dependence of Afρ across the survey

Driver `notebooks/afrho_trends.py`; science in `ztfcomet/activity.py`.  Tables
in `results/afrho/` (`trends.csv`, `peaks.csv`, `breaks.csv`, `colour.csv`),
figures `fig/afrho/<comet>_trend.png` and `fig/afrho/survey_overview.png`.
56 comets; r-band is the main series and g is analysed separately; apertures
ρ = 10,000 and 20,000 km.

## Method

- **Quantity**: phase-corrected A(0°)fρ, so the phase angle cannot masquerade
  as a heliocentric trend.
- **Points**: `quality_ok` rows, minus frames whose photometry centre sits more
  than half an aperture radius from the ephemeris (the `flag_centroid` gate is
  a 3-FWHM star-capture rule and lets those through).
- **Phases are defined by the peak of activity, not perihelion.**  A comet
  sampled on both sides of T_p gets a smoothed peak (Gaussian kernel,
  max(10 d, 2 × median gap)); when the peak is *bracketed* — interior to the
  coverage and ≥ 2 × RMS below the maximum at both ends — the series splits
  into *rising* and *fading* there.  When the maximum sits at an edge, the
  whole series is one phase (still rising, or already fading).  A comet
  sampled on one side only is one leg named by its orbital direction; its
  peak is not measurable and no split is attempted.  Perihelion-split rows
  remain in `trends.csv` under `split == "perihelion"` for reference only.
- **Fit**: `log Afρ = a − x log r_h`, weighted.  Every slope carries a formal
  error, a *scaled* error (× √χ²_red where the scatter exceeds the errors —
  the one quoted), and a bootstrap interval.
- **Broken law**: every leg with ≥ 10 points over ≥ 0.2 dex is also fitted
  with a continuous two-slope law, the break searched over every position
  leaving ≥ 4 points per side, and compared with the single law by ΔBIC
  **on errors inflated so the single law has χ²_red = 1**.  On raw χ² the
  intrinsic scatter makes any kink in the last four points look decisive
  (17 of 22 legs "preferred" a break; 13 do once scaled).  A break is claimed
  at ΔBIC > 6.  The same legs are also split at a fixed 3 au.
- **Colour**: g and r frames of the same aperture within one day are paired;
  the colour is the excess over solar, −2.5 log₁₀(Afρ_g / Afρ_r), in which no
  solar constants enter and the phase correction cancels.  A change point is
  a *step* between two levels (a kink cannot hold two plateaus), tested the
  same way against one level.
- **Grade**: A needs ≥ 10 points, ≥ 0.2 dex in r_h, scaled slope error ≤ 0.5,
  χ²_red ≤ 3; B ≥ 5 points, ≥ 0.1 dex, error ≤ 1; C ≥ 0.05 dex, error ≤ 2;
  D is a number, not a measurement.  The baseline leads: several comets span
  < 0.02 dex, where any slope is arithmetic.

## Headline numbers (r-band, primary phases, grades A/B)

| phase | 10,000 km: median x [16–84%] | 20,000 km |
|---|---|---|
| rising (to the peak) | 6.64 [3.91, 8.08] (n = 3) | 5.96 [3.19, 8.56] (n = 4) |
| fading (from the peak) | 4.04 [4.04, 4.04] (n = 1) | 3.60 [3.60, 3.60] (n = 1) |
| inbound only | 4.35 [-0.22, 6.85] (n = 8) | 3.92 [0.40, 6.47] (n = 6) |
| outbound only | 2.87 [0.69, 4.16] (n = 12) | 2.50 [1.11, 3.71] (n = 11) |

Only **four comets** have a bracketed peak with enough points on both sides
to fit both phases (24P, 2023 R1, 2025 K1, 2024 G4), and only 24P has both
phases at grade B: it rises as r_h^−8.8 and fades as r_h^−4.0.  2025 K1
rises as r_h^−2.6 (B) but its fading side has eight points.  The rest of the
sample is one-phase, and there the earlier picture stands: **outbound-only
comets fade uniformly** (12 reliable indices, median 2.9, 16–84% 0.7–4.2)
while **inbound-only ones scatter from −0.2 to 6.9**.

## Activity peaks (r-band, bracketed)

| comet | ρ (km) | n rising / fading | peak T−T_p (d) | 16–84% |
|---|---|---|---|---|
| 2023R1 | 10k | 44 / 35 | **+17** | [+6, +23] |
| 2023R1 | 20k | 80 / 31 | **+19** | [+13, +26] |
| 2024G4 | 10k | 30 / 32 | **-59** | [-238, -58] |
| 2024G4 | 20k | 41 / 48 | **-84** | [-87, +79] |
| 2025K1 | 10k | 40 / 8 | **+38** | [-1, +40] |
| 2025K1 | 20k | 28 / 8 | **+37** | [-3, +40] |
| 235P | 10k | 0 / 81 | **-29** | [-101, -26] |
| 235P | 20k | 9 / 70 | **-26** | [-45, -10] |
| 24P | 10k | 27 / 67 | **+11** | [+6, +13] |
| 24P | 20k | 20 / 53 | **+13** | [+7, +16] |
| 47P | 10k | 38 / 17 | **-34** | [-35, -31] |
| 47P | 20k | 41 / 20 | **-17** | [-19, -14] |

Three peak after perihelion (24P +11 d, 2023 R1 +17 d, 2025 K1 +38 d), two
before (47P −34 d, 235P −29 d); 2024 G4's is nominal (±100 d at 5 au).  The
10,000 and 20,000 km peaks agree within days.  Nine two-sided comets have
their maximum at the edge of coverage and are treated as one phase.

## Where a single power law fails

| comet | ρ (km) | phase | grade | r_h range | break (au) [16–84%] | x inside → outside | ΔBIC |
|---|---|---|---|---|---|---|---|
| 24P | 10k | rising | B | 1.18–2.05 | **1.24** [1.24, 1.31] | 21.9 ± 1.3 → 6.8 ± 0.2 | 14 |
| 2025K1 | 10k | rising | B | 0.69–2.82 | **1.32** [1.17, 1.50] | 3.5 ± 0.1 → 1.2 ± 0.2 | 15 |
| 24P | 20k | rising | B | 1.18–2.05 | **1.39** [1.39, 1.47] | 12.0 ± 0.5 → 4.4 ± 0.5 | 9 |
| 24P | 20k | fading | B | 1.24–2.39 | **1.47** [1.45, 1.65] | 5.4 ± 0.4 → 3.3 ± 0.1 | 7 |
| 2025A6 | 20k | inbound | A | 0.95–1.81 | **1.53** [1.50, 1.53] | 4.9 ± 0.1 → 5.7 ± 0.2 | 8 |
| 2025A6 | 10k | inbound | A | 0.95–1.81 | **1.53** [1.53, 1.56] | 4.6 ± 0.1 → 5.4 ± 0.2 | 10 |
| 10P | 10k | inbound | B | 1.44–3.63 | **1.85** [1.85, 1.89] | 9.6 ± 0.3 → -0.2 ± 0.3 | 42 |
| 217P | 10k | outbound | B | 1.54–3.42 | **2.90** [2.81, 2.90] | 2.8 ± 0.2 → -4.2 ± 1.3 | 10 |
| 217P | 20k | outbound | B | 1.54–3.42 | **2.90** [2.83, 2.90] | 2.9 ± 0.2 → -3.7 ± 1.2 | 10 |
| 2024E1 | 20k | inbound | B | 2.04–4.66 | **3.34** [3.34, 3.37] | -0.7 ± 0.1 → 1.0 ± 0.1 | 36 |
| 2024E1 | 10k | inbound | B | 2.01–4.66 | **3.37** [3.30, 3.49] | -0.7 ± 0.0 → 0.6 ± 0.1 | 29 |
| 2023A3 | 10k | outbound | A | 3.03–5.18 | **4.82** [4.55, 4.82] | 2.3 ± 0.2 → 6.9 ± 1.1 | 10 |
| 2023A3 | 20k | outbound | B | 3.08–6.38 | **5.68** [5.18, 5.68] | 2.7 ± 0.1 → -5.2 ± 1.2 | 14 |

Thirteen breaks survive the scaled test, and they sort into three kinds:

1. **Genuine changes of regime.**  10P is flat from 3.6 to 1.85 au and then
   rises as r_h^−9.6 — an activity onset at **1.85 au** (ΔBIC 42).  2024 E1's
   inbound activity *peaks at 3.4 au* and declines toward perihelion (ΔBIC
   29–36 at both apertures).  2025 K1's rise steepens inside 1.3 au.
2. **Curvature near the peak.**  24P's breaks at 1.24–1.47 au sit four points
   from the end of each phase: the rise flattening into the peak and the
   fading steep just after it.  A power law is the wrong local model within a
   few weeks of a peak; the break is telling you that.
3. **End-of-range features.**  217P beyond 2.9 au and 2023 A3 at 20,000 km
   beyond 5.7 au are the last few points turning *up* — an outburst or bad
   frames, not a heliocentric law; the two apertures disagree for 2023 A3.

## The 3 au hypothesis

Only four comets straddle 3 au with a usable leg:

| comet | ρ (km) | phase | x (< 3 au) | x (> 3 au) | split preferred |
|---|---|---|---|---|---|
| 10P | 10k | inbound | +5.60 ± 0.29 | -11.76 ± 2.70 | yes (ΔBIC 19) |
| 2023C2 | 10k | outbound | -0.41 ± 0.36 | -0.53 ± 0.09 | no (ΔBIC -4) |
| 2023C2 | 20k | outbound | -0.64 ± 0.31 | -0.47 ± 0.07 | no (ΔBIC -4) |
| 2024E1 | 10k | inbound | -1.04 ± 0.07 | +0.33 ± 0.07 | yes (ΔBIC 27) |
| 2024E1 | 20k | inbound | -1.16 ± 0.09 | +0.82 ± 0.07 | yes (ΔBIC 35) |
| 217P | 10k | outbound | +2.70 ± 0.19 | -6.16 ± 1.82 | yes (ΔBIC 13) |
| 217P | 20k | outbound | +2.75 ± 0.18 | -5.37 ± 1.73 | yes (ΔBIC 12) |

The split is "preferred" for three of them, but in **opposite senses**: 10P and
217P are steeper inside 3 au, 2024 E1 is the reverse (its maximum is at
3.4 au), and 2023 C2 shows no change.  Where the free search finds a break it
lands at 1.3, 1.5, 1.85, 2.9, 3.4 or 4.8 au — the transition is
comet-specific, not a property of 3 au.  The survey's ecliptic-comet windows
mostly lie entirely inside or entirely outside 3 au, so this is a four-comet
verdict, not a null result for the population.

## Dust colour with r_h

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
| 217P | 10k | 40 | 2.03–3.42 | +0.176 | -0.24 ± 0.43 | none | — |
| 217P | 20k | 36 | 2.03–3.42 | +0.149 | +0.02 ± 0.45 | none | — |
| 24P | 10k | 65 | 1.19–2.39 | -0.065 | +1.45 ± 0.19 | **1.53** [1.47, 1.58] | **+0.23 ± 0.03** |
| 24P | 20k | 43 | 1.19–2.37 | -0.214 | +2.70 ± 0.21 | **1.50** [1.49, 1.50] | **+0.44 ± 0.04** |
| 40P | 10k | 21 | 1.95–2.83 | +0.181 | -0.07 ± 0.58 | none | — |
| 40P | 20k | 10 | 1.95–2.76 | +0.165 | -2.17 ± 1.14 | none | — |

Eleven comets have enough same-night pairs.  Across all 22 series the dust is
**+0.10 mag redder than the Sun** (16–84%: −0.05 to +0.14) with no overall
trend in r_h (median slope −0.06 mag/dex).  One change is preferred, **24P at
1.50–1.53 au**: bluer inside, by +0.23 ± 0.03 mag at 10,000 km and
**+0.44 ± 0.04 at 20,000 km**.  That the step doubles with aperture is the
signature of gas, not dust: ZTF g carries the C₂ Swan bands, C₂ extends
further from the nucleus than the dust that dominates r, so the larger
aperture is more gas-affected and looks bluer near perihelion.  24P's
fading-phase colour slope (+1.1 to +2.3 mag/dex) is the same effect; the
survey does not show a dust colour that changes with r_h.

## Per-comet indices, r-band, primary phases

Scaled error, (N), grade; A/B in bold.

| comet | r_h (au) | phase | 10,000 km | 20,000 km |
|---|---|---|---|---|
| 10P | 1.44–3.63 | inbound | **+4.75 ± 0.33 (60) B** | +3.12 ± 1.71 (37) C |
| 124P | 2.11–2.48 | inbound | +1.75 ± 1.88 (17) C | +3.02 ± 4.27 (6) D |
| 131P | 2.48–2.67 | inbound | +11.87 ± 6.63 (12) D | +21.24 ± 5.86 (4) D |
| 145P | 1.89–2.48 | inbound | **+8.62 ± 0.32 (38) B** | **+9.29 ± 0.41 (38) B** |
| 164P | 1.88–2.78 | outbound | **-0.62 ± 0.18 (75) B** | **-0.72 ± 0.19 (67) B** |
| 171P | 1.77–2.09 | fading | +9.02 ± 0.87 (40) C | +12.66 ± 2.38 (30) D |
| 172P | 3.40–3.51 | inbound | -1.59 ± 5.20 (18) D | +13.21 ± 14.97 (10) D |
| 2021G2 | 5.18–6.89 | outbound | **+1.07 ± 0.29 (17) B** | +1.34 ± 4.00 (84) D |
| 2022E2 | 3.95–6.16 | outbound | +2.16 ± 0.18 (37) C | **+2.27 ± 0.08 (88) B** |
| 2022N2 | 3.83–4.25 | rising | — | +0.04 ± 0.14 (84) D |
| 2022N2 | 3.83–4.25 | fading | +0.20 ± 0.12 (61) D | — |
| 2022QE78 | 5.48–5.64 | fading | — | -4.01 ± 0.87 (80) D |
| 2022QE78 | 5.48–5.64 | outbound | -6.06 ± 6.65 (5) D | — |
| 2022R6 | 6.58–6.74 | outbound | — | +1.74 ± 1.26 (49) D |
| 2023A3 | 3.03–6.38 | outbound | **+2.75 ± 0.15 (26) A** | **+2.46 ± 0.12 (45) B** |
| 2023C2 | 2.75–6.02 | outbound | **-0.51 ± 0.05 (56) A** | **-0.49 ± 0.05 (65) A** |
| 2023F3 | 5.71–6.10 | outbound | — | +1.37 ± 0.99 (26) D |
| 2023H5 | 4.32–4.43 | inbound | -6.81 ± 2.39 (8) D | -4.10 ± 1.54 (23) D |
| 2023R1 | 3.57–4.85 | rising | **+6.64 ± 0.35 (44) B** | **+3.75 ± 0.19 (80) B** |
| 2023R1 | 3.57–4.85 | fading | +15.89 ± 1.30 (35) D | +10.95 ± 1.27 (31) D |
| 2023RS61 | 8.54–8.97 | inbound | — | -43.20 ± 5.71 (27) D |
| 2023T3 | 3.68–5.80 | outbound | **+2.61 ± 0.28 (9) B** | +1.52 ± 0.20 (22) C |
| 2023V1 | 5.09–5.33 | fading | +3.92 ± 4.90 (15) D | +5.54 ± 2.37 (35) D |
| 2024A1 | 3.92–4.37 | fading | +2.14 ± 0.54 (57) D | +2.06 ± 0.38 (87) D |
| 2024E1 | 2.01–4.66 | inbound | **-0.34 ± 0.05 (43) B** | **-0.06 ± 0.07 (50) B** |
| 2024G4 | 4.90–5.73 | rising | +1.17 ± 0.83 (30) D | +2.83 ± 0.68 (41) C |
| 2024G4 | 4.90–5.73 | fading | +2.04 ± 7.14 (32) D | -8.39 ± 14.73 (48) D |
| 2024J3 | 4.41–6.03 | inbound | **+0.66 ± 0.11 (10) B** | **+0.51 ± 0.07 (56) B** |
| 2024L5 | 3.43–5.04 | outbound | **+4.08 ± 0.12 (62) B** | **+3.66 ± 0.10 (92) B** |
| 2024N1 | 4.40–4.58 | rising | — | -5.59 ± 6.55 (33) D |
| 2024N1 | 4.40–4.58 | fading | -3.86 ± 4.97 (31) D | — |
| 2025A6 | 0.95–1.81 | inbound | **+4.72 ± 0.05 (54) A** | **+5.02 ± 0.05 (53) A** |
| 2025K1 | 0.69–2.82 | rising | **+2.62 ± 0.09 (40) B** | **+2.68 ± 0.07 (28) B** |
| 2025K1 | 0.69–2.82 | fading | +0.66 ± 1.96 (8) C | +1.26 ± 1.74 (8) C |
| 2025L1 | 1.71–1.89 | outbound | +30.94 ± 2.43 (7) D | +15.27 ± 4.61 (8) D |
| 2025L2 | 3.22–3.54 | inbound | +2.63 ± 1.51 (9) D | +3.31 ± 2.07 (8) D |
| 2025M2 | 5.51–8.17 | inbound | +0.92 ± 2.11 (11) D | **+2.82 ± 0.15 (76) B** |
| 2025Q3 | 2.13–2.99 | inbound | **+7.14 ± 0.47 (27) B** | **+5.76 ± 0.97 (23) B** |
| 2025R1 | 1.98–2.41 | fading | -0.63 ± 0.71 (18) C | — |
| 2025R1 | 1.98–2.41 | inbound | — | +0.57 ± 1.27 (10) C |
| 2025R2 | 1.09–2.62 | outbound | **+3.36 ± 0.12 (19) A** | **+2.94 ± 0.12 (18) A** |
| 2025W2 | 1.46–1.51 | inbound | -4.93 ± 11.18 (3) D | +6.95 ± 10.78 (3) D |
| 210P | 0.61–2.19 | outbound | **+3.33 ± 0.04 (56) B** | **+3.07 ± 0.04 (44) B** |
| 217P | 1.54–3.42 | outbound | **+2.08 ± 0.14 (120) B** | **+2.18 ± 0.13 (105) B** |
| 235P | 1.98–2.83 | rising | — | **+8.91 ± 0.96 (9) B** |
| 235P | 1.98–2.83 | fading | +7.86 ± 0.16 (81) C | +7.39 ± 0.18 (70) C |
| 240P | 2.12–2.34 | rising | -0.94 ± 2.66 (49) D | -0.33 ± 2.74 (44) D |
| 24P | 1.18–2.39 | rising | **+8.76 ± 0.33 (27) B** | **+8.18 ± 0.72 (20) B** |
| 24P | 1.18–2.39 | fading | **+4.04 ± 0.10 (67) B** | **+3.60 ± 0.09 (53) B** |
| 261P | 2.01–2.42 | inbound | +15.35 ± 0.73 (57) C | +16.35 ± 0.99 (30) D |
| 29P | 6.27–6.31 | outbound | — | -134.84 ± 54.26 (72) D |
| 2P | 3.01–4.09 | inbound | **-0.42 ± 0.84 (26) B** | -17.39 ± 27.15 (8) D |
| 302P | 3.32–3.87 | outbound | +2.41 ± 0.61 (40) C | +1.24 ± 0.53 (45) C |
| 306P | 1.29–1.48 | outbound | -0.15 ± 3.39 (5) D | — |
| 40P | 1.83–2.87 | outbound | **+4.91 ± 0.64 (81) B** | **+3.79 ± 0.78 (72) B** |
| 43P | 2.45–3.10 | outbound | **+4.40 ± 0.10 (95) B** | **+3.95 ± 0.12 (88) B** |
| 47P | 2.81–3.14 | rising | +11.69 ± 0.94 (38) D | +16.13 ± 1.40 (41) D |
| 47P | 2.81–3.14 | fading | +24.61 ± 7.31 (17) D | +2.55 ± 1.65 (20) D |
| 486P | 2.33–2.86 | outbound | +0.73 ± 0.43 (59) C | -0.05 ± 0.60 (57) C |
| 48P | 2.37–3.15 | outbound | **+3.00 ± 0.29 (35) B** | **+2.50 ± 0.25 (35) B** |
| 491P | 3.91–4.64 | outbound | -3.51 ± 0.76 (64) C | -4.46 ± 0.87 (71) C |
| 493P | 3.82–4.26 | rising | +6.22 ± 1.35 (39) D | +6.14 ± 1.84 (60) D |
| 499P | 0.93–1.47 | outbound | +2.88 ± 0.58 (4) D | — |
| 63P | 2.45–2.76 | inbound | +6.72 ± 2.53 (10) D | +9.05 ± 3.13 (8) D |
| 78P | 2.55–3.54 | inbound | **+3.98 ± 0.29 (51) B** | +21.65 ± 15.59 (43) D |

## Caveats

- Splitting at the peak means a phase can straddle perihelion (2025 K1's
  rising phase runs to +38 d): r_h is not monotonic in time within it, and
  the power law is then a net slope over two orbital directions.
- A power law assumes monotonic activity; near a peak or across an outburst
  it is the wrong local model, and the break test says so rather than the
  index.
- The 3 au test rests on four comets; the colour result on eleven.
- g-band Afρ carries C₂ emission and is reported separately.
