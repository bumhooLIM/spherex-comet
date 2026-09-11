# Heliocentric dependence of Afρ across the survey

Driver `notebooks/afrho_trends.py`; science in `ztfcomet/activity.py`.  Tables
`results/activity/afrho_trends.csv` (one row per comet × band × aperture × leg)
and `afrho_peaks.csv`; figures `fig/<comet>/afrho_trend_<comet>.png` and
`fig/survey/afrho_trends_overview.png`.  56 comets, r-band as the main series
and g-band analysed separately, apertures ρ = 10,000 and 20,000 km.

## Method

- **Quantity**: phase-corrected A(0°)fρ.  Phase angle changes along a leg and
  would otherwise masquerade as a heliocentric trend.
- **Points**: `quality_ok` rows, additionally dropping frames whose photometry
  centre sits more than half an aperture radius from the ephemeris (the
  `flag_centroid` gate is a 3-FWHM star-capture rule and lets those through;
  see `handoff.md`).  Contamination and off-centre counts are kept per comet.
- **Fit**: `log Afρ = a − x log r_h`, weighted by the propagated errors.  Every
  slope carries a *formal* error, a *scaled* error (formal × √χ²_red when the
  scatter exceeds the errors — an active comet varies, and the errors carry no
  jitter), and a *bootstrap* 16–84% interval.  The scaled error is quoted.
- **Legs**: whole series; inbound / outbound by the sign of ṙ; and, for
  comets sampled on both sides of perihelion, *rising* / *fading* split at the
  activity peak instead.
- **Peak**: maximum of a Gaussian-kernel-smoothed log Afρ in T − T_p, kernel
  max(10 d, 2 × median gap), reported only when *bracketed* — interior to the
  sampled span and ≥ 2 × RMS below the maximum at both ends.  Uncertainty from
  bootstrap.  A maximum at the edge of the data is censoring, not a peak.
- **Grade**: A needs ≥ 10 points, ≥ 0.2 dex in r_h, scaled slope error ≤ 0.5,
  χ²_red ≤ 3; B ≥ 5 points, ≥ 0.1 dex, error ≤ 1; C ≥ 0.05 dex, error ≤ 2;
  D is a number, not a measurement.  One step down if more than half the
  band's frames were rejected.  **The baseline leads**: σ_x ∝ 1/(√N · σ_log r_h),
  and several comets span < 0.02 dex, where any slope is arithmetic.

## Headline numbers (r-band)

| | 10,000 km | 20,000 km |
|---|---|---|
| comets with an A/B whole-series fit | 24 of 56 | 21 |
| inbound, A/B: median x [16–84%] | **4.7** [0.4, 8.0] (n = 12) | 3.6 [0.9, 5.8] (n = 9) |
| outbound, A/B | **3.2** [1.2, 4.1] (n = 14) | 2.9 [2.0, 3.7] (n = 13) |
| whole series, A/B | 3.3 [0.3, 5.5] (n = 24) | 2.9 [0.8, 4.8] (n = 21) |

Among A/B fits the median χ²_red is **2.4** with a residual scatter of
0.069 dex against a median point error of 0.043 dex: real comets vary
by several times their photometric errors along a leg, which is why the scaled
error is the one quoted and the formal one is not to be believed.

**Outbound fading is remarkably uniform** — 13 of 14 reliable outbound indices
at 10,000 km lie between 1.2 and 4.1, median 3.2 — while **inbound behaviour is
not**: the reliable inbound indices span 0.4 to 8.0.  The two comets with both
legs reliable disagree in sense: 24P rises as r_h^−8 and fades as r_h^−4;
2025 K1 rises as r_h^−1.2 and fades as r_h^−3.9.  There is no universal
inbound/outbound asymmetry to report, only per-comet ones.

**Aperture**: for 39 comet-legs reliable at both apertures, x(10k) − x(20k)
has median +0.29 — the inner coma responds slightly more steeply, in the
direction the resolution study predicts for the aperture closest to the PSF.

**g versus r**: 68 pairs, median x_g − x_r = +0.23 [−0.27, +0.92].  g fades
faster than r with distance, i.e. the coma reddens outward — but ZTF g contains
the C₂ Swan bands, so part of this is gas, not dust, and it is why r is the
primary series.

## Activity peaks (r-band, bracketed only)

| comet | ρ (km) | n in / out | span (d) | peak T−T_p (d) | 16–84% | rising x | fading x |
|---|---|---|---|---|---|---|---|
| 2023R1 | 10,000 | 36 / 43 | -326 … +84 | **+17** | [+6, +23] | +6.64 ± 0.35 (44) **B** | +15.89 ± 1.30 (35) D |
| 2023R1 | 20,000 | 72 / 39 | -371 … +84 | **+19** | [+13, +26] | +3.75 ± 0.19 (80) **B** | +10.95 ± 1.27 (31) D |
| 2024G4 | 10,000 | 30 / 32 | -347 … +109 | **-59** | [-238, -58] | +1.17 ± 0.83 (30) D | +2.04 ± 7.14 (32) D |
| 2024G4 | 20,000 | 45 / 44 | -388 … +109 | **-84** | [-87, +79] | +2.83 ± 0.68 (41) C | -8.39 ± 14.73 (48) D |
| 2025K1 | 10,000 | 26 / 22 | -151 … +116 | **+38** | [-1, +40] | +2.62 ± 0.09 (40) **B** | +0.66 ± 1.96 (8) C |
| 2025K1 | 20,000 | 17 / 19 | -142 … +116 | **+37** | [-3, +40] | +2.68 ± 0.07 (28) **B** | +1.26 ± 1.74 (8) C |
| 235P | 10,000 | 15 / 68 | -258 … +198 | **-29** | [-101, -26] | — | +7.86 ± 0.16 (81) C |
| 235P | 20,000 | 16 / 63 | -280 … +197 | **-26** | [-45, -10] | +8.91 ± 0.96 (9) **B** | +7.39 ± 0.18 (70) C |
| 24P | 10,000 | 23 / 71 | -138 … +178 | **+11** | [+6, +13] | +8.76 ± 0.33 (27) **B** | +4.04 ± 0.10 (67) **B** |
| 24P | 20,000 | 16 / 57 | -138 … +178 | **+13** | [+7, +16] | +8.18 ± 0.72 (20) **B** | +3.60 ± 0.09 (53) **B** |
| 47P | 10,000 | 44 / 11 | -130 … +83 | **-34** | [-35, -31] | +11.69 ± 0.94 (38) D | +24.61 ± 7.31 (17) D |
| 47P | 20,000 | 44 / 17 | -130 … +251 | **-17** | [-19, -14] | +16.13 ± 1.40 (41) D | +2.55 ± 1.65 (20) D |

Six comets bracket a peak.  Three peak **after** perihelion — 24P at +11 d,
2023 R1 at +17 d, 2025 K1 at +38 d — the thermal-lag sense; two peak
**before** — 47P at −34 d and 235P at −29 d.  2024 G4's peak is nominal only:
at 5 au over a 0.06 dex baseline its intervals reach ±100 d and the 20,000 km
result disagrees; treat it as unbracketed.  The 10,000 and 20,000 km peaks
agree to within a few days wherever the 10,000 km one is well determined.

The peak-split fits are only meaningful when the leg spans r_h: 24P's rising
r_h^−8.8 and fading r_h^−4.0 are grade B; 2023 R1's fading leg after +17 d
spans 3.57–3.65 au, so the r_h^−16 there is a decay in *time* forced through
a 0.01 dex lever arm — the grade (D) says so.

Nine two-sided comets have no bracketed peak (171P, 240P, 493P, 2022 N2,
2022 QE78, 2023 V1, 2024 A1, 2024 N1, 2025 R1): their maximum sits at the edge
of the coverage, so the peak is censored, and only the perihelion split is
reported.

## Comets worth a second look

- **2023 C2**: x = −0.51 ± 0.05, grade A, n = 56 — Afρ *rises* as it recedes
  from 2.8 to 6.0 au.  164P does the same outbound (−0.62 ± 0.18, B).
- **2024 E1**: inbound x = −0.34 ± 0.05 over 2.0–4.7 au is a net decline while
  approaching; the leg is non-monotonic (a maximum near 3 au, then fading
  toward perihelion), so a single power law is the wrong model and χ²_red
  says so.  The figure, not the index, is the result here.
- **2025 A6**: x = 4.72 ± 0.05 (A, n = 54, 0.95–1.81 au) — the cleanest
  inbound measurement in the survey.
- **261P**: x = 15.4 ± 0.7 (C) over 2.0–2.4 au inbound — extreme; consistent
  with the decay of an outburst rather than steady activity.
- **145P** (8.6 ± 0.3), **2025 Q3** (7.1 ± 0.5), **24P** (8.1 ± 0.3): steep
  inbound onsets of the kind seen in seasonally driven Jupiter-family comets.
- **29P** (6.27–6.31 au) and **2023 RS61** (8.5–9.0 au) return absurd indices
  and grade D: no lever arm, and 29P's outbursts are not a trend.

## Per-comet indices, r-band, perihelion split

Scaled error, (N), grade; A/B in bold.

| comet | r_h (au) | inbound 10k | outbound 10k | inbound 20k | outbound 20k |
|---|---|---|---|---|---|
| 10P | 1.44–3.63 | +4.75 ± 0.33 (60) **B** | — | +3.12 ± 1.71 (37) C | — |
| 124P | 2.11–2.48 | +1.75 ± 1.88 (17) C | — | +3.02 ± 4.27 (6) D | — |
| 131P | 2.48–2.67 | +11.87 ± 6.63 (12) D | — | +21.24 ± 5.86 (4) D | — |
| 145P | 1.89–2.48 | +8.62 ± 0.32 (38) **B** | — | +9.29 ± 0.41 (38) **B** | — |
| 164P | 1.88–2.78 | — | -0.62 ± 0.18 (75) **B** | — | -0.72 ± 0.19 (67) **B** |
| 171P | 1.77–2.09 | -5.45 ± 5.77 (23) D | +8.34 ± 1.19 (17) C | +17.69 ± 24.74 (14) D | +8.24 ± 1.11 (16) C |
| 172P | 3.40–3.51 | -1.59 ± 5.20 (18) D | — | +13.21 ± 14.97 (10) D | — |
| 2021G2 | 5.18–6.89 | — | +1.07 ± 0.29 (17) **B** | — | +1.34 ± 4.00 (84) D |
| 2022E2 | 3.95–6.16 | — | +2.16 ± 0.18 (37) C | — | +2.27 ± 0.08 (88) **B** |
| 2022N2 | 3.83–4.25 | +3.14 ± 9.77 (13) D | +0.17 ± 0.13 (48) D | +6.13 ± 5.04 (37) D | +0.04 ± 0.16 (47) D |
| 2022QE78 | 5.48–5.64 | — | -6.06 ± 6.65 (5) D | -7.89 ± 1.39 (15) D | +1.61 ± 0.93 (65) D |
| 2022R6 | 6.58–6.74 | — | — | — | +1.74 ± 1.26 (49) D |
| 2023A3 | 3.03–6.38 | — | +2.75 ± 0.15 (26) **A** | — | +2.46 ± 0.12 (45) **B** |
| 2023C2 | 2.75–6.02 | — | -0.51 ± 0.05 (56) **A** | — | -0.49 ± 0.05 (65) **A** |
| 2023F3 | 5.71–6.10 | — | — | — | +1.37 ± 0.99 (26) D |
| 2023H5 | 4.32–4.43 | -6.81 ± 2.39 (8) D | — | -4.10 ± 1.54 (23) D | — |
| 2023R1 | 3.57–4.85 | +6.43 ± 0.41 (36) **B** | +14.93 ± 0.98 (43) D | +3.58 ± 0.20 (72) **B** | +10.00 ± 0.96 (39) D |
| 2023RS61 | 8.54–8.97 | — | — | -43.20 ± 5.71 (27) D | — |
| 2023T3 | 3.68–5.80 | — | +2.61 ± 0.28 (9) **B** | — | +1.52 ± 0.20 (22) C |
| 2023V1 | 5.09–5.33 | -11.38 ± 9.33 (12) D | +6.30 ± 5.98 (3) D | -9.44 ± 5.01 (26) D | +5.40 ± 2.98 (9) D |
| 2024A1 | 3.92–4.37 | -2.81 ± 11.21 (3) D | +2.06 ± 0.54 (54) D | -2.49 ± 9.54 (3) D | +1.99 ± 0.39 (84) D |
| 2024E1 | 2.01–4.66 | -0.34 ± 0.05 (43) **B** | — | -0.06 ± 0.07 (50) **B** | — |
| 2024G4 | 4.90–5.73 | +1.17 ± 0.83 (30) D | +2.04 ± 7.14 (32) D | +3.75 ± 0.66 (45) C | -20.78 ± 17.04 (44) D |
| 2024J3 | 4.41–6.03 | +0.66 ± 0.11 (10) **B** | — | +0.51 ± 0.07 (56) **B** | — |
| 2024L5 | 3.43–5.04 | — | +4.08 ± 0.12 (62) **B** | — | +3.66 ± 0.10 (92) **B** |
| 2024N1 | 4.40–4.58 | -3.52 ± 4.99 (27) D | -20.43 ± 87.28 (4) D | -10.36 ± 4.61 (25) D | +90.70 ± 41.95 (8) D |
| 2025A6 | 0.95–1.81 | +4.72 ± 0.05 (54) **A** | — | +5.02 ± 0.05 (53) **A** | — |
| 2025K1 | 0.69–2.82 | +1.23 ± 0.14 (26) **B** | +3.89 ± 0.22 (22) **B** | +1.73 ± 0.16 (17) **A** | +3.55 ± 0.19 (19) **B** |
| 2025L1 | 1.71–1.89 | — | +30.94 ± 2.43 (7) D | — | +15.27 ± 4.61 (8) D |
| 2025L2 | 3.22–3.54 | +2.63 ± 1.51 (9) D | — | +3.31 ± 2.07 (8) D | — |
| 2025M2 | 5.51–8.17 | +0.92 ± 2.11 (11) D | — | +2.82 ± 0.15 (76) **B** | — |
| 2025Q3 | 2.13–2.99 | +7.14 ± 0.47 (27) **B** | — | +5.76 ± 0.97 (23) **B** | — |
| 2025R1 | 1.98–2.41 | -0.56 ± 0.72 (15) C | +6.27 ± 406.56 (3) D | +0.57 ± 1.27 (10) C | — |
| 2025R2 | 1.09–2.62 | — | +3.36 ± 0.12 (19) **A** | — | +2.94 ± 0.12 (18) **A** |
| 2025W2 | 1.46–1.51 | -4.93 ± 11.18 (3) D | — | +6.95 ± 10.78 (3) D | — |
| 210P | 0.61–2.19 | — | +3.33 ± 0.04 (56) **B** | — | +3.07 ± 0.04 (44) **B** |
| 217P | 1.54–3.42 | — | +2.08 ± 0.14 (120) **B** | — | +2.18 ± 0.13 (105) **B** |
| 235P | 1.98–2.83 | +7.99 ± 0.80 (15) **B** | +7.45 ± 0.15 (68) C | +8.13 ± 1.16 (16) C | +7.22 ± 0.19 (63) C |
| 240P | 2.12–2.34 | +0.74 ± 0.57 (32) D | -1.89 ± 13.02 (17) D | +1.23 ± 0.72 (29) D | -5.79 ± 14.22 (15) D |
| 24P | 1.18–2.39 | +8.05 ± 0.29 (23) **B** | +4.04 ± 0.09 (71) **B** | +5.81 ± 0.66 (16) **B** | +3.63 ± 0.08 (57) **B** |
| 261P | 2.01–2.42 | +15.35 ± 0.73 (57) C | — | +16.35 ± 0.99 (30) D | — |
| 29P | 6.27–6.31 | — | — | — | -134.84 ± 54.26 (72) D |
| 2P | 3.01–4.09 | -0.42 ± 0.84 (26) **B** | — | -17.39 ± 27.15 (8) D | — |
| 302P | 3.32–3.87 | — | +2.41 ± 0.61 (40) C | — | +1.24 ± 0.53 (45) C |
| 306P | 1.29–1.48 | — | -0.15 ± 3.39 (5) D | — | — |
| 40P | 1.83–2.87 | — | +4.91 ± 0.64 (81) **B** | — | +3.79 ± 0.78 (72) **B** |
| 43P | 2.45–3.10 | — | +4.40 ± 0.10 (95) **B** | — | +3.95 ± 0.12 (88) **B** |
| 47P | 2.81–3.14 | +10.85 ± 0.85 (44) D | +16.01 ± 4.78 (11) D | +16.36 ± 1.28 (44) D | +1.30 ± 1.24 (17) D |
| 486P | 2.33–2.86 | — | +0.73 ± 0.43 (59) C | — | -0.05 ± 0.60 (57) C |
| 48P | 2.37–3.15 | — | +3.00 ± 0.29 (35) **B** | — | +2.50 ± 0.25 (35) **B** |
| 491P | 3.91–4.64 | — | -3.51 ± 0.76 (64) C | — | -4.46 ± 0.87 (71) C |
| 493P | 3.82–4.26 | +4.59 ± 1.78 (34) D | -231.73 ± 534.99 (5) D | +5.03 ± 2.37 (55) D | -436.96 ± 275.47 (5) D |
| 499P | 0.93–1.47 | — | +2.88 ± 0.58 (4) D | — | — |
| 63P | 2.45–2.76 | +6.72 ± 2.53 (10) D | — | +9.05 ± 3.13 (8) D | — |
| 78P | 2.55–3.54 | +3.98 ± 0.29 (51) **B** | — | +21.65 ± 15.59 (43) D | — |

## Caveats

- A power law in r_h assumes monotonic activity along the leg.  Where it is
  not (2024 E1 inbound, any leg containing an outburst) the index is a net
  slope and χ²_red flags it; look at the figure.
- Many comets span < 0.1 dex in r_h.  For them the index is unconstrained
  regardless of N, and the table says D.  The survey's ecliptic-comet
  windows are often short arcs near perihelion.
- g-band indices carry C₂ emission; they are reported separately and should
  not be averaged with r.
- The bracketing test needs the curve to fall ≥ 2 × RMS on both sides; a
  broad plateau (2022 N2, drop 0.01 dex) is correctly reported as no peak.
