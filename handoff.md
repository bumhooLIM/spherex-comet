# Handoff

## Current State
**Blocked on hardware, mid-task.**  The T7 dropped off the bus with
`OSError: [Errno 5] Input/output error` during the profile re-run (second
time in a day).  `results/` is a MIX: 22 targets (2P … 240P in survey order)
have profiles from the new optocentre code; the other 46 still have pre-fix
profiles.  **Do not compute survey-wide profile statistics until the 46 are
re-run.**  Photometry tables are untouched and valid.

Branch `fix/profile-peak-centring` (4 commits, not pushed):
- `refine_centre()` in `ztfcomet/profile.py`: comet profiles now start at
  the optocentre, not the winpos centroid.  Search disc covers winpos and
  the ephemeris (midpoint, half-separation + 2 px, cap 4); the move guard is
  the 99.5th percentile of the statistic's own null distribution, measured
  per box shape through the same code path (false-positive rate 0.9%).
  Live: 10P off-peak 21% → 0.0%, 40P 40% → 4.6%.  141 tests.
- `survey.py` refuses to start when `ZTFCOMET_DATA` points at a missing dir.

Merged today: PR #2 (survey), PR #4 (1/ρ holds to the coverage edge).

## Next Steps
1. Mount the T7, confirm `/Volumes/T7/data/ztf-comet` exists, then:
   `ZTFCOMET_DATA=/Volumes/T7/data/ztf-comet python notebooks/survey.py
   --steps profile figures --no-resume --vmag-max 20 --start 2025-03-01
   --targets $(cat <scratch>/remaining_targets.txt)` — the 46 are 261P
   through 2025W2; plain resume would skip them (every row is already `ok`).
2. `notebooks/profile_survey.py`, then `profile_resolution.py --target 24P`.
3. Before/after vs the pre-fix numbers in `doc/survey_summary_68comets.md`
   (7,698 frames, 4,145 clean, −1.73 vs −4.39, 95.6%) and
   `doc/profile_survey_1rho.md` (2,744 fitted, 2,118 good, 691 off-peak).
4. Push and open the PR for `fix/profile-peak-centring`.
5. **Follow-up to raise:** `flag_centroid` fires at 3 FWHM from the
   ephemeris — a star-capture rule, not an aperture rule.  Among clean rows
   the centroid shift is 5.4 px at the 99th percentile; at ρ = 10,000 km
   (r_ap ≈ 6 px) 5.3% of clean rows are > 0.5 r_ap off and 1.4% > 1 r_ap
   (nucleus at the aperture edge; ~30–40% low bias on a 1/ρ coma).  The
   per-frame `recentre_shift_pix` now measures this directly.  Proposed
   fix: centre the aperture photometry with `refine_centre` too, or add a
   critical flag at shift > 0.5 r_ap.  Changes Afρ for a few % of rows.

## Blind Spots / Dead Ends
- A sigma cut cannot guard a "max minus one sample" statistic (3σ moved
  17.5% of pure-noise frames, calibrated 4σ still 4%); only the null
  quantile through the same code path holds.  A fixed 2 px disc lands
  confidently on its edge when the nucleus is further away — hence the
  ephemeris-covering disc.
- Gate on the survey log line, not a file mtime; `tail -1` hides pytest's exit.
- χ² from the profile tables is ~16× too large (correlated sub-pixels).
  The free sky can run to its bound; check `at_bound`.
