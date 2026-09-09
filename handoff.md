# Handoff

## Current State
The 68-comet survey is **complete and fully reduced**.  Full statistics:
`doc/survey_summary_68comets.md`.

- 56 targets with photometry, 12 with no data (3 never reach Vmag 20, 9 have
  no ZTF coverage at the positions that do).
- Downloads **98.75% complete** (8,734 / 8,845), 0 corrupt.  The 111 missing
  are permanent: 80 archive 404s and 31 that IRSA rejects with "Cutout does not
  overlap image" — the comet falls off the quadrant, so they would have tripped
  `flag_outside` anyway.  Two repair passes recovered zero of them; do not
  retry a third time.
- 41,533 measurement rows over five apertures, 52.1% clean.  Contamination is
  the largest rejection cause at 19.9%.
- Radial profiles for **all 56** targets: median comet slope −1.73 against
  −4.39 for field stars, 95.6% of clean frames shallower than the stars.
- Afρ figures for every target use the r_h − q abscissa with the Kepler date
  axis; the targets 1–12 that were stuck on signed-r_h have been regenerated.

## Next Steps
Nothing is outstanding from the survey itself.  Open scientific questions:

1. Extend `profile_resolution.py` beyond 24P — does 1/ρ hold past ~12,000 km
   for the targets with large ρ_max (2023V1 37,000 km; 2025M2 34,600 km)?
   The steep-slope end of the summary table is the sample to test.
2. Decide whether slopes below about −2.5 should be excluded from any
   population analysis, given they track frame depth rather than coma
   structure.
3. `fix/horizons-query-correctness` is well ahead of `master` and has never
   been merged.  Consider opening the PR.

## Blind Spots / Dead Ends
- A **running process keeps its imported modules**.  The profile step and the
  perihelion axis were both added mid-run and silently did not apply to
  targets already in flight; both needed a full re-reduction afterwards.
  Restart the batch after changing the package.
- **`--no-resume` rewrites every status row**, including targets it cannot
  act on.  A download-only pass relabelled the 12 no-data targets `no_urls`,
  destroying the `no_frames` / `no_epochs` distinction; it was recovered from
  the run logs.  Read the logs, not just the status CSV, when reconstructing.
- The run **hung 54 min on one Horizons call** at 0% CPU while the service
  answered fresh requests in 0.8 s.  `socket.setdefaulttimeout(300)` in
  `survey.py` now bounds this.
- macOS writes `._name` AppleDouble sidecars on the exFAT T7; they matched the
  FITS glob and were reported as half the frames being unreadable.  Filtered in
  `phot.py` — do not "fix" that by re-downloading.
