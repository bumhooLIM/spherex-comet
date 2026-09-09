# Handoff

## Current State
The 68-comet survey is complete and fully reduced; PR #2 merged to `master`.
Statistics: `doc/survey_summary_68comets.md`.

**The 1/ρ question is answered** (`doc/profile_survey_1rho.md`, driver
`notebooks/profile_survey.py`): with the PSF core and sky modelled, the coma
matches 1/ρ to 1–2% from ~1,500 km to the edge of coverage — 37,000–45,000 km —
at every S/N band and every r_h from 0.6 to 6.7 au.  No turnover anywhere.
The 24P "ρ_max ≥ 12,000 km" was coverage.  Corrected slope converges on
m = 1.00 at S/N(10 px) > 60 (naive 1.42); per-target medians span 0.57–1.42,
which is real comet-to-comet dispersion.

Slopes steeper than −2.5 are not measurements: 185 of 213 such frames never
got a constrained fit, and the 28 that did correct from 2.67 to 1.24.  Rule:
use the PSF-model + free-sky slope at S/N(10 px) > 20, off the sky bound.

This is on branch `analysis/profile-survey-1rho`, not yet pushed.

## Next Steps
1. **Decide on the profile centring fix.**  20% of clean S/N > 5 frames were
   excluded because the comet profile peaks off-centre; the off-peak rate
   tracks the winpos shift from ephemeris (6% at < 0.5 px, 57% at 2–5 px).
   Fix = peak re-centring within ±2 px in `run_profiles`, then re-run
   `--steps profile figures` (~30 min).  Changes numbers in the merged
   summary for 40P, 235P, 47P, 2024E1 most.  User's call.
2. Open a PR for `analysis/profile-survey-1rho` if the note is wanted on master.
3. Optional: tighten the sky bound in `fit_coma_model` (now 3× outer SB) — 23%
   of fits sit on it, 42% at S/N 5–20.  A bound tied to the sky uncertainty
   would keep more low-S/N frames constrained.

## Blind Spots / Dead Ends
- **Oversampled sub-pixels are correlated**: `n_clip` overstates independent
  samples by ~oversample², so any χ² from the profile tables is ~16× too large.
  Scale before thresholding, or don't threshold.
- **The free sky can run to its bound** and the fit then reports a slope that
  means nothing; always check `at_bound`.  16 slopes also hit the m = 3
  ceiling.
- **`fit_coma_model` needs an on-peak profile.**  Off-peak frames drive the
  nucleus term to zero and the sky to its bound (seen on 2022E2).
- A running process keeps its imported modules — restart batches after
  package changes.  `--no-resume` rewrites every status row.  AppleDouble
  `._` sidecars on the exFAT SSD match FITS globs; filtered in `phot.py`.
