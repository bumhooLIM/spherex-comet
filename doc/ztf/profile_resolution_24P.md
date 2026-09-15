# Why does 24P's coma profile look steeper inbound? — a resolution study

**Date:** 2026-09-08  **Data:** 24P/Schaumasse, 295 ZTF frames (240 clean), 2025-03 → 2026-09
**Code:** `ztfcomet.profile`, `notebooks/profile_resolution.py`
**Figures:** `fig/24P/profile_resolution_24P.png`, `profile_validity_24P.png`, `profile_correction_24P.png`
**Tables:** `results/24P/profile_resolution_24P.csv`, `results/profile_feasibility_survey.csv`

## The question

The power-law slope of the coma's radial surface brightness, fitted from 1.5 FWHM
to 10 px, is −1.05 post-perihelion (a steady-state 1/ρ coma) but −1.25 with a
tail to −2.9 on the inbound leg. Hypothesis: at large observer distance Δ a fixed
10 px window spans a much larger cometocentric distance ρ, reaching beyond where
1/ρ holds, so the fit steepens. Is that right, is it an inherent limitation of
ZTF's 1″ pixels, and can it be corrected?

## Scales

| | Δ (au) | km/px | ρ at 1.5 FWHM | ρ at 10 px |
|---|---|---|---|---|
| near perihelion | 0.6 | 437 | 1 070 km | 4 500 km |
| far inbound | 2.1 | 1 509 | 5 970 km | 15 500 km |

The fit window in km changes by 3.5× across the apparition, so the hypothesis is
at least plausible on its face.

## Test 1 — does the slope track ρ(10 px)?

Spearman correlations of the naive slope, clean frames:

| | ρ(10 px) | Δ | Afρ | SB(10 px)/σ_sky |
|---|---|---|---|---|
| pre-perihelion (n=51) | **−0.63** (p 7×10⁻⁷) | −0.63 | **+0.71** (p 5×10⁻⁹) | +0.56 |
| post-perihelion (n=189) | −0.02 | −0.02 | — | +0.02 |

Inbound, the slope does correlate with ρ(10 px) — but ρ(10 px) ∝ Δ, so that is
the same correlation, and brightness correlates *more* strongly. The decisive
control is the post-perihelion leg: those frames span the same ρ(10 px) range
(4 500–13 500 km) and show **no steepening at all**. If a large ρ window alone
caused it, they would be steep too.

## Test 2 — where does 1/ρ actually hold?

ρ·SB(ρ), normalised per frame, is flat where the profile is 1/ρ.

- **Post-perihelion: flat (0.99–1.05) from 1 300 km out to 12 100 km**, the largest
  bin with enough frames. 1/ρ holds to **ρ_max ≥ 12 000 km** — a lower bound set
  by coverage, not a detected break.
- Pre-perihelion: flat to ~5 400 km, drops to 0.71 at 8 100 km, then **recovers to
  0.96 at 12 100 km**. A physical break does not recover; noise does.
- In pixel units, by Δ tercile: flat to 10 px in all three. No departure tied to a
  pixel radius either.

Only 7 % of 24P frames have ρ(10 px) beyond 12 000 km. Accordingly, **refitting in
a fixed km window changes nothing** (inbound median −1.25 → −1.25). The
hypothesis in its literal form is not supported for 24P.

## Test 3 — what does cause it?

The six steepest inbound frames all have Afρ 15–33 cm (the rest: 60–600 cm) and
SB(10 px) only 7–14× the formal sky uncertainty. Post-perihelion frames with low
outer-annulus S/N are steep too (`profile_correction` right panel): **S/N is the
common cause on both legs.**

Nucleus/PSF leakage is ruled out. The forward model puts the nucleus at 1–2 % of
the flux, and synthetic tests show that even a nucleus carrying 70 % of the flux
under a 4.5 px PSF only biases the naive slope to −1.23. That mechanism cannot
reach −2.

The mechanism is sky subtraction. Letting the PSF-convolved model fit a free sky
offset:

| inbound, clean (n=51) | median [16–84 %] | within [−1.3, −0.7] |
|---|---|---|
| naive power law | −1.25 [−1.62, −0.95] | 26 / 51 |
| fixed km window | −1.25 [−1.62, −0.95] | — |
| PSF-convolved model | −1.08 [−1.36, −0.87] | 35 / 51 |
| **PSF model + free sky** | **−1.01 [−1.19, −0.76]** | **40 / 51** |
| post-perihelion, same | −0.94 [−1.10, −0.74] | — |

The fitted offsets are **negative in 66 % of frames** (sky over-estimated), median
−4 DN. Two contributions:

1. **Coma light in the sky annulus.** The photometry's sky annulus sits at
   60–75 px; if 1/ρ continues, the coma there is 2–5 DN/px. The median ratio of
   fitted offset to that expectation is 0.93 — the right size and sign on
   average — but with a wide spread (Spearman +0.27).
2. **Background structure.** On the faintest frames the offset (−3 to −6 DN) is
   10× any coma at the annulus. The *formal* sky uncertainty (σ_sky/√N ≈ 0.07 DN)
   understates the real systematic by ~50×; "SB(10 px)/σ_sky ≈ 10" is really an
   outer S/N near 1 once background structure at the few-DN level is counted.

Either way: a few DN of sky error on a profile whose outer annuli are at 1–2 DN
tilts the whole tail.

## Test 4 — can it be corrected?

| technique | effect | why |
|---|---|---|
| **oversampling** | **none** (native − oversampled slope = +0.009, 16–84 % [−0.02, +0.05]) | adds no spatial information; only improves annulus sampling geometry |
| fixed km window | none for 24P | ρ_max ≥ 12 000 km exceeds the window on 93 % of frames |
| PSF-convolved model | partial (−1.25 → −1.08) | may legitimately use the high-S/N inner annuli, which the naive fit excludes |
| PSF model + free sky | **most of it** (→ −1.01) | absorbs the sky systematic the inner annuli are immune to |

The remaining scatter in the sky-fitted slope (some frames at −0.3 to −0.5) is
the price: a free offset is partly degenerate with a shallow coma. It is the
correct tool for faint frames, not a replacement for the naive fit on bright ones.

Also implied for the Afρ photometry: coma in the sky annulus biases the aperture
sum low by ~0.7 % on bright frames and up to ~5 % on the faintest. Not fixed here.

## The inherent limit, and other comets

With ρ_max ≥ 12 000 km (a lower bound) the usable window
[1.5 FWHM, min(ρ_max, 10 px)] shrinks with Δ:

| target | Δ (au) | km/px | ρ_max in px | usable window (px) | frames feasible |
|---|---|---|---|---|---|
| 24P | 0.84 | 615 | 19.7 | 6.8 | 100 % |
| 10P | 1.40 | 1 031 | 11.7 | 6.0 | 81 % |
| 43P, 78P | 2.3 | 1 680 | 7.2 | 3.9 | 86–89 % |
| 172P | 2.64 | 1 940 | 6.2 | 2.7 | 81 % |
| 2024 E1 | 3.33 | 2 447 | 4.9 | 1.9 | 45 % |
| **29P** | **5.48** | **4 019** | **3.0** | **0** | **0 %** |

For 29P the whole region where 1/ρ is *known* to hold lies inside the 1.5 FWHM
core. **That is the inherent ZTF limitation:** no annulus outside the PSF measures
the coma, and no amount of oversampling changes it. The forward model is the only
route, and there the coma and nucleus terms become strongly degenerate.

Two caveats keep this from being the final word:

- ρ_max = 12 000 km is a lower bound from 24P's coverage. Comets at Δ ≥ 1.8 au
  have their whole 10 px window beyond it (`frac_rho10_beyond_rhomax = 1.0`).
  Whether 1/ρ holds at 15 000–40 000 km is untested — and testable with the
  survey's own bright distant targets once their profiles exist
  (`profile_resolution.py --target 29P`).
- The sky term should be validated against an independent sky estimate (a wider,
  coma-free annulus) before it is trusted on faint frames in bulk.

## Bottom line

The inbound steepness of 24P is **not** a pixel-scale/ρ-range effect and **not**
nucleus leakage. It is sky-subtraction error on faint frames, and it is largely
correctable by a PSF-convolved model with a free sky term — not by oversampling.
The genuine ZTF limitation appears only when the 1/ρ-valid region falls inside
the PSF core, which for ρ_max ≈ 12 000 km begins around Δ ≈ 3 au and is total by
Δ ≈ 5 au.
