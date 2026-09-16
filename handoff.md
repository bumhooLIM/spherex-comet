# Handoff — `spherex-comet`
_2026-09-16: the review memo of 2026-09-15 applied case by case (`doc/comspec/case_revisions.md`); gas fits, Afρ epochs, figures and the deck reproduced.  Git history restarted the same day: the old 80-commit history is archived at `~/Desktop/claude/spherex-comet-git-history-backup-2026-09-16` (the former `.git`; `git clone` it to browse) and, as of 2026-09-16, still on GitHub `bumhooLIM/spherex-comet`._
## Current State
**Main result:** `results/comspec/gas_fit.csv` — 147 fits over 68 comets (run `5741d6611f60`,
`spherex_comspec` 1.3.0): robust ≥ 3σ (n_eff ≥ 2) H₂O 18 / CO₂ 28 / CO 4 phases → 37 comets
(36 clean), 19 more with marginal values only, 2 H₂O detections `rejected`; χ²_ν median 4.2;
192 phases (240P's 2–3 merged).  Input photometry unchanged (`spherex_apphot` 2.2.0, `9fcc7ea3871a`).
ZTF: Afρ at 120 / 124 of 164 phases (10k / 20k), attached to 91 / 97 fitted phases.  Deck
`doc/figures.pptx` (2026-09-16): 3 slides per fitted phase, with the previous state's 3 slides
placed before them for every revised group; the old deck is `doc/figures_before_rev260915.pptx`.
`pytest tests/` — 234 pass (7 new: case-revision mechanics, Afρ overrides).
**Case revisions** (`spherex_comspec/revisions.py`, `CASE_REVISIONS`, 34 groups + memo-only
entries): per-band window edges, fixed orders, top-N point exclusions, negative / one-sided
waivers, flag drops, excluded emission channels, reduced coverage, rejected species; applied in
`continuum.process_group` and `fitting.fit_production_rates` wherever a group is processed, and
recorded in the `revision` column, `notes`, `caveats` and the figures.  240P: `GroupingConfig`
(`manual_edges`, `delta_tol_exempt`).  Afρ: `ztfcomet.config.AFRHO_EPOCH_OVERRIDES` (47P S1
grade-D extrapolation, 210P S1–S4 from the outbound law, 217P S1 without outburst frames).
`dc_main_norev` (role `rules`) is the rule-only run and reproduces the 2026-09-14 census exactly.
Snapshots for the comparison: `fig/comspec/previous_rev260915/` (figures + `gas_fit.csv` of the
revised groups); the before/after report is `notebooks/comspec/case_revision_report.py BEFORE_DIR`.
`notebooks/comspec/figures.ipynb` (2026-09-16): Q vs ⟨r_h⟩ with the 3σ limits, the ≥ 2-phase comets as
labelled tracks, Q / A(0°)fρ at 10k / 20k → `fig/comspec/Q_*.png`; paths via `spherex_comspec.directory`.

## Next Steps
1. Review the deck pairs and `case_revisions.md` §4 (waived negative continua are upper bounds;
   single-channel CO₂ detections at n_eff ≈ 1; 240P's 26 % Δ spread; 210P's symmetric-activity
   assumption); drop or keep each entry.  Re-read the registry whenever windows / LSF change.
2. Placeholder 1 (LSF / band shape) for the bright comets (χ²_ν 10²–10³ at 20 000 km).
3. `scripts/apphot/make_figures.py … --no-cutouts` to refresh `fig/apphot/`; `pip install -e .`.
4. Flag `b` re-test at the fixed apertures; annulus 100 000 / 200 000 km; aperture thresholds.
5. Decide: quote marginal values or limits only; unify the ZTF (15 pt) and SPHEREx (20 pt) styles.

## Blind Spots / Dead Ends
- **A revised continuum window need not bracket the band**: the saved point neighbourhood must
  span the emission window too, or the band's channels vanish from the fit (2023 R1 S4).
- **"Exclude top N on a side" counts what the figure shows**, including points inside the
  neighbouring CO window; ranking only continuum candidates emptied 2025 R1 S4's red side.
- **Rule 1b would re-split 240P's merged phase** — the target is exempted, and the regroup
  assertion skips `delta_tol_exempt` targets.  After any regroup: `scripts/ztf/spherex_windows.py`
  → `afrho_trends.py --targets …` (a partial run does not rebuild the windows) → `attach_afrho_ztf.py`;
  stale figure files of vanished phases (240P ph4) must be deleted by hand.
- `phase_images.py` embeds Q and Afρ in its header text: regenerate it after the attach step.
- **The deck builder's figure cache must be keyed on the source directory**: keyed on the file
  name it served the previous snapshot's cropped figures on the revised slides (fixed 2026-09-16).
- **Emission files are keyed on (target, aperture)**; `arc` in `phase_map.csv` is an arc index.
- Isolated comspec runs set `COMSPEC_{APPHOT,DATA,RESULT,FIG}_DIR` before import; band windows
  are module constants to monkeypatch (`notebooks/comspec/method_matrix.py`).
- Never mask Gaia sources; never project the whole Gaia subset through the WCS; FLAG bit 21
  fires on the comet; `pgrep -fc` is invalid on macOS; macOS writes `._` sidecars on exFAT.
