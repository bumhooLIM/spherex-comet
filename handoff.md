# Handoff

## Current State
Three threads.  The survey is complete and reduced (PR #2, #4 merged).

1. **Heliocentric Afρ trends** — done, on branch
   `analysis/afrho-heliocentric-trends`: `ztfcomet/activity.py`,
   `notebooks/afrho_trends.py`, `doc/afrho_heliocentric_trends.md`.
   56 comets, r main / g separate, ρ = 10k and 20k.  24 comets have an A/B
   whole-series index at 10k.  Outbound fading is uniform (median x = 3.2,
   13/14 within 1.2–4.1); inbound spans 0.4–8.0.  Six bracketed activity
   peaks: 24P +11 d, 2023R1 +17 d, 2025K1 +38 d, 47P −34 d, 235P −29 d,
   2024G4 nominal only.  2023C2 brightens outbound (x = −0.51 ± 0.05, A).
   Also fixed here: `orbit.fetch_elements` lacked `id_type="smallbody"`
   (bare designations returned None silently); regression test added.
2. **Profile centring fix** — PR #5 open.  Code complete (141 tests); the
   profile re-run is 22/68 because the T7 dropped out.  46 targets still
   have pre-fix profiles; the 46-target list and resume command are in the
   PR.  Survey-wide profile numbers in the docs are pre-fix until then.
3. **Aperture off-centre follow-up** (not started): `flag_centroid` is a
   3-FWHM star-capture rule; at ρ = 10k, 5.3% of clean rows have the centre
   > 0.5 r_ap off the nucleus.  The trend analysis drops those rows itself.

## Next Steps
1. PR for `analysis/afrho-heliocentric-trends` (one commit, 149 tests).
2. Mount the T7 → finish the profile re-run (PR #5), refit 1/ρ, update docs.
3. Decide the aperture-centring fix: `refine_centre` in phot, or a flag.

## Blind Spots / Dead Ends
- `fetch_elements` returns None on failure rather than raising; log it.
  Resolve comets from the survey-list *designation* ("2022 E2"), never the
  results-directory slug ("2022E2") — Horizons cannot match the latter.
- A power law in r_h needs a lever arm: σ_x ∝ 1/(√N σ_log r_h).  Several
  comets span < 0.02 dex; grade on the baseline before N.
- Scaled errors, not formal: A/B fits have median χ²_red ≫ 1.
- A sigma cut cannot guard a max-minus-sample statistic (17.5% false
  positives at 3σ); use the null quantile through the same code path.
- Gate on log lines, not file mtimes; `tail -1` after pytest hides its exit.
