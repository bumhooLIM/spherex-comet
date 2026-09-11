# Handoff

## Current State
Three threads.  The survey is complete and reduced (PR #2, #4 merged).

1. **Heliocentric Afρ trends** — on `feat/afrho-peak-fits-and-reorg`
   (which merges the unmerged `analysis/afrho-heliocentric-trends` and
   `fix/profile-peak-centring`; rebases to nothing once both PRs land).
   Phases split at the *activity peak*, not perihelion; broken power law
   with a free break and a fixed 3 au split (ΔBIC on scaled errors — raw
   χ² called 17/22 legs broken, scaled 13); g−r colour as the excess over
   solar with a step change point.  `doc/afrho_heliocentric_trends.md`.
   Results: only 4 comets have both phases fittable; outbound-only comets
   fade uniformly (x ≈ 2.9), inbound-only scatter −0.2…6.9.  Breaks: 10P
   onset at 1.85 au, 2024E1 maximum at 3.4 au inbound; the 3 au hypothesis
   is comet-specific, not universal (4 comets, opposite senses).  Colour:
   +0.10 mag redder than solar, no r_h trend; 24P bluer inside 1.5 au by an
   amount that doubles with aperture — C₂ in g, not dust.
   **Outputs are now organised by subject**: results/{photometry,profile,
   afrho}/, fig/{photometry,profile,afrho}/; `result_dir`/`fig_dir` take a
   subject and reject target names; `notebooks/reorganize_outputs.py`
   migrated 8,277 files.  `figure.ipynb` still writes the old per-target
   figure names and needs updating to `fig_path()`.
2. **Profile centring fix** — PR #5 open.  Code complete (141 tests); the
   profile re-run is 22/68 because the T7 dropped out.  46 targets still
   have pre-fix profiles; the 46-target list and resume command are in the
   PR.  Survey-wide profile numbers in the docs are pre-fix until then.
3. **Aperture off-centre follow-up** (not started): `flag_centroid` is a
   3-FWHM star-capture rule; at ρ = 10k, 5.3% of clean rows have the centre
   > 0.5 r_ap off the nucleus.  The trend analysis drops those rows itself.

## Next Steps
1. Merge PR #5 and the activity branch, then PR `feat/afrho-peak-fits-and-reorg`.
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
