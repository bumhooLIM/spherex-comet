# What the background annulus does to Afρ

Two ways to define the background of an aperture measurement of a coma:

1. **a sky annulus far from the aperture** — what `ztfcomet.phot` does:
   inner radius 3ρ, outer radius 4ρ + 20 px, sigma-clipped median;
2. **an annulus just outside the aperture edge.**

They are not two estimates of the same thing.  They measure different
integrals of the coma, and only one of them is Afρ.

## Setup

Let the coma surface brightness be a power law on a constant sky,

    I(r) = C r^(−m) + S,        m = 1 for steady-state outflow.

The true coma flux inside ρ, against the sky alone, is

    F_true(ρ) = ∫₀^ρ 2πr · C r^(−m) dr = 2πC ρ^(2−m) / (2 − m)       (m < 2),

which for m = 1 is 2πCρ.  Afρ is defined from this flux
(A'Hearn et al. 1984): Afρ = (2 Δ r_h / ρ)² (F_true / F_⊙) ρ, and for m = 1
it is independent of ρ.  It is an *integral* quantity — the cross-section
of dust inside ρ — and for a steady-state coma it scales with the dust
production rate, Afρ ∝ Q_dust / v.  That is the "activity level".

Any annulus estimate of the background is B = S + ⟨I_coma⟩_annulus, and the
measured flux is F = F_true − πρ² ⟨I_coma⟩_annulus.  The two methods differ
only in what ⟨I_coma⟩ is.

## Method 1: far annulus

For a median over an annulus [r₁, r₂] the estimator returns the coma
brightness at the radius that halves the annulus area,
r_med² = (r₁² + r₂²)/2, so B₁ = S + C r_med^(−m) and

    F₁ = F_true · [ 1 − (2 − m)/2 · (ρ / r_med)^m ].

For m = 1 the bias is ½ ρ/r_med.  In our geometry r_med = 6.4ρ at
10,000 km (r₂ = 8.6ρ) and 4.2ρ at 40,000 km (the fixed 20 px pad matters
less for a large aperture), so the far annulus removes **7.8%** of the coma
flux at 10k and **11.8%** at 40k.  Small, constant for a given geometry and
slope, and exactly correctable: multiply by 1/[1 − (2−m)/2 (ρ/r_med)^m]
with m from the radial profile (≈ 1 for this survey to 40,000 km,
`doc/profile_survey_1rho.md`).  Method 1 approximates the definition of
Afρ and its error is a known factor of 1.08–1.13.

Its weakness is not the coma in the annulus; it is the accuracy of the sky
*level* itself when the aperture is sky-dominated
(`doc/afrho_aperture_systematics.md`).  That is a separate problem, and
the local annulus does not solve it.

## Method 2: local annulus at the aperture edge

A thin annulus at the edge returns B₂ = S + I(ρ), so

    F₂ = F_true − πρ² C ρ^(−m)
       = ∫₀^ρ 2πr [ I(r) − I(ρ) ] dr
       = − ∫₀^ρ πr² I′(r) dr                (integrating by parts).

The last form is the interpretation: **the local-annulus flux is the
gradient-weighted integral of the profile** — how much the coma inside ρ
stands *above its own edge*.  For a power law it is a fixed fraction of the
true flux,

    F₂ / F_true = m / 2 .

| coma | F₂ / F_true |
|---|---|
| uniform disc (m = 0) | **0** — however bright |
| m = 0.5 | 0.25 |
| steady-state 1/ρ | **0.50** |
| m = 1.5 | 0.75 |
| point source (nucleus) | 1 |

So a local annulus returns half of Afρ for a steady-state coma, nothing at
all for a uniform bright coma, and the full flux for a bare nucleus.  It
is a measure of **central concentration**, not of dust content: it
multiplies Afρ by m/2 and so entangles the profile slope with the activity
level, changes from comet to comet and from epoch to epoch as m changes,
and — because the nucleus flux passes through it untouched while the coma
is halved — *doubles the nucleus fraction* relative to the truth.  For the
"pure dust activity of the inner coma" it moves in exactly the wrong
direction.

## The survey confirms the factor

From the radial-profile tables (background: the far annulus), the flux
inside 8 px against the far background versus against the level at 8 px,
for 4,086 clean frames of 55 comets:

| profile slope m (naive) | n | F₂/F₁ measured | m/2 |
|---|---|---|---|
| 0.75–1.0 | 193 | **0.46** | 0.45 |
| 1.0–1.25 | 395 | **0.51** | 0.56 |
| 1.25–1.5 | 787 | 0.60 | 0.70 |
| 1.5–2.0 | 1367 | 0.68 | 0.87 |

Where the profile is cleanly 1/ρ the ratio is ½ (24P 0.49, 210P 0.46,
48P 0.51, 2025 R2 0.43).  At steeper *naive* slopes the measured ratio
falls below m/2 because the naive slope is itself biased steep by sky
error and the PSF core (the resolution study); the inner 8 px are closer to
1/ρ than the naive fit says.  Survey-wide the local annulus would have
returned 0.67 [0.54, 0.83] of the far-annulus flux — a scatter driven by
morphology, not by activity.

## What the local annulus *is* good for

It is a differential measurement: a background gradient, scattered light
or a crowded field is subtracted locally, and the ratio F₂/F₁ is a
one-number **concentration index** of the coma.  For a time series of one
comet with a stable profile it tracks F_true to within a constant.  It
cannot be compared across comets, against literature Afρ, or turned into
Q_dust.

## Which is correct for the dust activity level of the inner coma

**Method 1**, the far sky annulus, with two refinements:

- correct the coma removed by the annulus with the factor above
  (1.08 at 10,000 km, 1.13 at 40,000 km for m = 1; per row it is
  1/[1 − ½ ρ_pix/r_med] with `rho_pix`, `sky_in_pix` and `sky_out_pix` already in
  the photometry table — a calibration change of 8–12% to every Afρ, so a
  decision rather than a bug fix);
- treat sky-dominated apertures (sky/flux > 10) as unreliable, since the
  sky *level* error, not the coma in the annulus, is what limits them.

If the goal is the dust *alone*, the remaining contaminant is the nucleus,
and the way to remove it is the PSF-modelled profile fit
(`ztfcomet.profile.fit_coma_model`, which returns `F_nuc`):
Afρ_dust = Afρ − Afρ_nucleus.  The local annulus does not remove the
nucleus; it removes half the dust.
