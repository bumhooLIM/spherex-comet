# Handoff

## Current State
Branch `feat/afrho-peak-fits-and-reorg` (holds the unmerged PR #5 and the activity
branch).  170 tests.  The chain finished cleanly (2026-09-11): all 68 profiles
recentred, the 12 sparse targets re-queried at V < 21, 7,783 cutout PNGs.

Afρ analysis (`doc/afrho_heliocentric_trends.md`; tables by
`notebooks/afrho_trends_report.py`): phases at the activity peak, outbursts
excluded, tails set aside, segments where a break is preferred, smoothed guide
curve, SPHEREx windows, grade-D fits not drawn.  Physics notes:
`doc/afrho_aperture_systematics.md`, `doc/afrho_background_annulus.md`.

**Afρ at the SPHEREx epochs (2026-09-14)**: `activity.afrho_at_epoch` gives, per
SPHEREx phase and aperture, the r-band A(0°)fρ at the phase's ⟨r_h⟩ — `direct`
(≥ 2 clean frames within 5 d of the window, moved to ⟨r_h⟩ along the local law),
`trend` (the fitted law of the orbital phase the epoch falls in, any grade inside
its data), `trend_extrap` (≤ 0.1 dex, grade A–C), else `none` with the reason.
`results/afrho/spherex_afrho.csv`: 165 phases of 56 comets, values for 116 (10k)
/ 120 (20k); red markers and a table strip on every `fig/afrho/trend/<T>.png`;
`fig/afrho/trend_slide/<T>_S<k>.png` per phase (highlighted) for the catalog deck.
`afrho_trends.py --figures-only` redraws from the saved tables; `--targets`
merges into them; a full run rebuilds `spherex_windows.csv`.  The SPHEREx
catalog attaches the table with its `scripts/attach_afrho_ztf.py`.

## Next Steps
1. Merge PR #5 and the activity branch, then PR this branch.
2. Apply the coma-in-annulus factor (1/[1 − ½ ρ_pix/r_med], 8–12%) and an
   advisory `sky_dominated` flag (sky/flux > 10) — calibration decisions
   spelled out in the two notes.
3. `figure.ipynb` still writes the old per-target figure names.
4. Phases without an Afρ (49 at 10k / 45 at 20k): 19–20 outside a grade-D law,
   17–15 on the orbital phase ZTF never covered; `extrap_grades` and
   `max_extrap_dex` in `afrho_at_epoch` are the knobs if the catalog needs more.

## Blind Spots / Dead Ends
- **The catalog's `arc` is an arc index, not a direction**: `in` until a
  perihelion resolved inside the SPHEREx coverage, so a comet observed only after
  perihelion (217P, 2025 A6) reads `in` throughout.  Decide the orbital phase from
  T_p (`elements.csv`), never from that label.
- A line drawn in axis-fraction y on an otherwise empty axes leaves the y view at
  ±0.05 and a log axis then fails in `tight_layout`; use `ax.text`.
- An outburst test against the recent *median* opens on any steep smooth rise;
  compare with the extrapolated trend and require decay (261P is an onset).
- BIC on raw χ² calls every kink decisive; inflate errors to χ²_red = 1 on the
  simpler model first.  A kinked line cannot fit a colour step.
- The coma-in-annulus fraction for a median estimator is ½ ρ/r_med with
  r_med² = (r₁² + r₂²)/2, not ρ/(2r₁).
- `df.flags` is a pandas attribute; `--targets` needs designations with spaces;
  a survey run that empties a target's epochs now says which cut.
- The T7 drops under sustained reads; the chain script records each stage's exit code.
