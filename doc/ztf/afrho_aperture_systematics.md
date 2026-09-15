# When Afρ rises with distance: real, or an aperture effect?  (C/2024 E1)

Driver `scripts/ztf/aperture_systematics.py`; figure `fig/ztf/afrho/systematics_2024E1.png`;
table `results/ztf/afrho/systematics_2024E1.csv`.

## The question

In the 10,000 km aperture C/2024 E1 brightens as it approaches from 4.7 to
~3.4 au and then *fades* toward perihelion; the 40,000 km aperture declines
monotonically over the same arc.  The broken-law test puts the change at
3.37 au [3.30, 3.49] with x = +0.59 ± 0.09 outside and −0.74 ± 0.08 inside
(both A/B).  Is the reversal a property of the comet, or of the photometry?
Two systematics were suspected: pixel resolution, and an effective aperture
too small to hold the coma.

## What the tables say

The diagnostic is the ratio of the two apertures.  For a 1/ρ coma
Afρ(10k)/Afρ(40k) = 1 at every distance; anything else is either a coma
steeper than 1/ρ or a bias that depends on Δ.

| r_h (au) | Δ (au) | Afρ(10k)/Afρ(40k) | m_eff from the curve of growth | m from the PSF-modelled profile (≤ 10 px) | sky/flux at 10k | sky/flux at 40k | apcor at 10k (mag) | ρ(10k)/FWHM |
|---|---|---|---|---|---|---|---|---|
| 2.1 | 2.44 | 1.33 | 1.20 | **0.92** (n = 13) | 4.9 | 24.8 | −0.011 | 2.6 |
| 2.8 | 2.69 | 1.48 | 1.28 | **0.83** (n = 12) | 1.2 | 7.2 | −0.014 | 2.7 |
| 3.4 | 3.13 | 1.71 | 1.36 | — | 1.0 | 6.6 | −0.035 | 2.3 |
| 3.6 | 3.31 | 1.81 | 1.42 | 1.10 (n = 1) | 2.3 | 16.5 | −0.038 | 2.3 |
| 3.9 | 3.63 | 1.91 | 1.45 | — | 2.4 | 11.9 | −0.057 | 2.2 |
| 4.3 | 4.15 | 2.07 | 1.43 | — | 2.2 | 15.7 | −0.065 | 1.8 |

The ratio climbs from 1.3 to 2.1 with distance.  Read as a coma slope through
the curve of growth (m_eff = 1 − dlog Afρ / dlog ρ across the five apertures)
it is 1.20 at 2 au and 1.45 at 4 au.  But the PSF-modelled radial profile
(`doc/ztf/profile_survey_1rho.md`), which fits the sky as a free parameter, gives
m = 0.83–0.92 on the same frames inside 3.3 au: a 1/ρ coma, or slightly
shallower, within 30,000 km.  The apertures disagree with the resolved
profile.  The difference between the two measurements is how they treat
the sky: the profile fits it; the aperture photometry subtracts a fixed
annulus.

## The four candidates

1. **Coma light in the sky annulus.**  The annulus runs from 3ρ to
   4ρ + 20 px and takes the median, which for a 1/ρ coma removes
   ½ ρ/r_med of the coma flux: 7.8% at 10,000 km, 11.8% at 40,000 km
   (`doc/ztf/afrho_background_annulus.md`).  A 4% differential, in the
   direction of *lowering* the ratio.  Not this.
2. **Aperture correction.**  A point-source correction applied to an
   extended source over-brightens the small aperture; the applied value
   grows from −0.011 to −0.065 mag with distance.  Six per cent at most,
   and in the direction of *raising* the far 10k points — removing it makes
   the rise to 3.4 au slightly steeper, not flatter.  Not this.
3. **Resolution.**  ρ(10k) is 1.8–2.7 FWHM.  PSF scattering removes flux
   from a small aperture on an extended source, which would *lower* the
   10k values at large Δ and push the ratio toward 1 — the opposite of what
   is seen.  Not this.
4. **The sky level.**  The aperture sum is sky-dominated: at 40,000 km the
   sky inside the aperture is 7–25 times the source flux, at 10,000 km 1–5
   times.  A systematic error of 1% in the sky level therefore removes
   7–25% of the 40k flux and 1–5% of the 10k flux.  To turn an intrinsic
   ratio of ~1.3 into the observed 2.1 at 4.3 au the 40k flux must be
   depressed by 40%, i.e. the sky over-estimated by 1.4% — the level at
   which a background gradient across a 90-pixel annulus, or coma light
   beyond 1/ρ, is ordinary.  **This is the one.**

The survey confirms it is not peculiar to 2024 E1: across all clean r-band
rows the median sky/flux at 40,000 km is 23–79 depending on Δ, **80% of
clean 40k rows exceed 10**, and even at 10,000 km 55% do.  The 30,000 and
40,000 km apertures are sky-systematics-limited for most of the survey; the
10,000–20,000 km apertures are the robust ones.

## Verdict

The concern was that the small aperture is the unreliable one.  It is the
reverse.  The 10,000 km series — Afρ rising to a maximum near 3.4 au and
then declining by ~25% toward perihelion — is the measurement least exposed
to the sky and unaffected by the two resolution-type systematics; the
decline inside 3.4 au appears in every aperture (10k −25%, 20k −23%, 30k
−14%, 40k −5%), weakening outward exactly as a sky bias on the larger
apertures would make it.  **The reversal is real in the inner coma.**  The
monotonic 40k decline is the artefact.

The aperture dependence of the decline may itself be physical — old dust
at 40,000 km lags a drop in production by weeks for slow grains — but that
cannot be separated from the sky bias with a fixed annulus.

## Recommendations

- Quote 10,000–20,000 km for faint or distant comets; treat 30,000 and
  40,000 km as indicative only where sky/flux > 10.
- Add an advisory flag `sky_dominated` at sky/flux > 10 (a table-only
  change, like `flag_anomalous_bright`), so the aperture figures can show it.
- The profile study's free-sky fit is the right estimator of the outer-coma
  slope; the aperture curve of growth is not, at these S/N.
