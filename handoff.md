# Handoff

## Current State
Branch `feat/afrho-peak-fits-and-reorg` (merges the unmerged PR #5 and the
activity branch; rebases to nothing once they land).  163 tests.

**Afρ trend analysis** (`doc/afrho_heliocentric_trends.md`): phases split at
the *activity peak* (interior maximum, plateau allowed); outbursts detected
(jump above the extrapolated trend, confirmed, must decay) and excluded;
isolated r_h tails set aside; segments where a broken law is preferred;
smoothed black guide curve; SPHEREx phase windows shaded; grade-D fits not
drawn.  `flag_anomalous_bright` is a new CRITICAL photometry flag (single
frame > 0.15 dex and 5σ above its time-neighbours, alone); 309 rows.
**C/2024 E1's Afρ reversal is real** and the large-aperture decline is the
artefact: the 30–40k apertures are sky-dominated (80% of clean 40k rows have
sky/flux > 10) — `doc/afrho_aperture_systematics.md`.

**Outputs by kind**: fig/afrho/{rh,apertures,trend,lightcurve,colour}/<T>.png
via `fig_kind_path`; `reorganize_outputs.py` migrated them.

**A chained job is running** (`<scratch>/chain.sh`, status in
`<scratch>/chain_status`): (1) profile+figures for all 68 (the recentring
fix, PR #5) → (2) re-query of the 12 sparse targets at **V < 21** (the survey
was already V < 20; the request said 20) → (3) cutout PNGs for every frame
(`survey.py --steps cutouts`, fig/photometry/<T>/cutout/).  ~5 h.

## Next Steps
1. When the chain ends: rerun `afrho_trends.py` (sparse targets changed),
   refresh the numbers in the trends note, commit; then push and PR.
2. Merge PR #5 and the activity branch first so this PR shows its own diff.
3. Add an advisory `sky_dominated` flag (sky/flux > 10, table-only) and show
   it on the aperture figures — recommended in the systematics note.
4. `figure.ipynb` still writes the old per-target figure names.

## Blind Spots / Dead Ends
- An outburst test against the recent *median* opens on any steep smooth
  rise and never closes; compare with the extrapolated trend, and require
  the window to decay (a jump that keeps rising is an onset: 261P).
- BIC on raw χ² calls every kink decisive when scatter ≫ errors; inflate
  the errors to χ²_red = 1 on the simpler model first.
- A kinked line cannot represent a colour *step*; use two levels.
- `df.flags` is a pandas attribute — use `df["flags"]`.
- `--targets` in survey.py takes designations with spaces ("2022 QE78");
  the results-file slug fails at Horizons.
- The T7 drops off the bus under sustained reads; `survey.py` refuses a
  pinned-but-missing root, and the chain records each stage's exit code.
