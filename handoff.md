# Handoff — `spherex-comet`
_2026-09-15: the three projects merged into this tree, every driver, notebook and test re-pointed and verified; science state unchanged from 2026-09-14._
## Current State
**Layout:** packages `ztfcomet/`, `spherex_apphot/`, `spherex_comspec/` at the root; drivers in
`scripts/<stage>/`, notebooks in `notebooks/<stage>/` (+ shared `rcparams.py`), tests in `tests/`,
products namespaced `results|fig|doc/<stage>/`.  One git repo (ztfcomet history + spherex-comspec
history via tag `spherex-comspec-import`; apphot added).  Verified: `pytest tests/` 226 pass
(169 + 34 + 23); an isolated comspec run reproduces `gas_fit.csv` for 24P exactly.
**Main result:** `results/comspec/gas_fit.csv` — 147 fits over 68 comets (run `a6463b4c610b`,
`spherex_comspec` 1.2.0): robust ≥ 3σ (n_eff ≥ 2) H₂O 17 / CO₂ 26 / CO 4 phases → 35 comets
(34 clean), 19 more with marginal values only; χ²_ν median 4.6.  Input photometry
`results/apphot/photometry/` (`spherex_apphot` 2.2.0, config `9fcc7ea3871a`, 2026-09-12: 68/68,
455 245 rows).  ZTF: 56 comets with photometry, Afρ trends, Afρ at 165 SPHEREx phases
(`results/ztf/afrho/spherex_afrho.csv`, values for 116 / 120 at 10k / 20k) attached as `afrho_*`
to `gas_fit.csv` (89 / 95 of 147 phases).  Deck `doc/figures.pptx` (2026-09-14).
Deleted on 2026-09-15: the stale 7.6 GB cutout PNGs, the duplicate `data/apphot`, the archived
old photometry, obsolete notebooks; `fig/apphot/<T>/` summaries still date from 2026-09-07.
Pre-merge backup of code/docs: `../_backup_pre-merge_2026-09-15.tar.gz`.

## Next Steps
1. `python scripts/apphot/make_figures.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 6 --no-cutouts`
   to refresh `fig/apphot/` for the current photometry; then `scripts/comspec/make_figure_slides.py`.
2. Run the suites once in the `spherex` env (`pytest tests/ -q`; `pip install -e .`), commit,
   add a remote for the monorepo; retire the two old repos (nothing was pushed by the merge).
3. Placeholder 1: the SPHEREx LSF / band shape for the bright comets (χ²_ν 10²–10³ at 20 000 km).
4. Re-test flag `b` at the fixed apertures (`dc_main_no_b`); the 2.42 µm band-3/4 step;
   annulus 100 000 / 200 000 km on 24P, 306P, 508P; aperture thresholds (3 au, 1.5 px, 60 000 km).
5. Decide: quote marginal values or limits only; 10 000 vs 20 000 km as the dust-context aperture.
6. ZTF: apply the coma-in-annulus factor (1/[1 − ½ ρ/r_med], 8–12 %) and the `sky_dominated`
   flag; `notebooks/ztf/figure.ipynb` still uses the old per-target figure names; merge any
   `analysis/*`, `fix/*` branch content not on the current branch.
7. Unify the figure style (`ztfcomet.rcparams` 15 pt vs `notebooks/rcparams.py` 20 pt).

## Blind Spots / Dead Ends
- **`arc` in `phase_map.csv` is an arc index** (`in` until a resolved perihelion): 217P and
  2025 A6 read `in` outbound.  Decide the orbital phase from T_p (`results/ztf/afrho/elements.csv`).
- **After a regroup: rerun `scripts/ztf/afrho_trends.py`, then `scripts/comspec/attach_afrho_ztf.py`**
  (moved phases are blanked as `stale`).  Emission files are keyed on (target, aperture).
- **Isolated comspec runs** set `COMSPEC_{APPHOT,DATA,RESULT,FIG}_DIR` before import (see
  `notebooks/comspec/method_matrix.py`); band windows are module constants to monkeypatch.
- **Never mask Gaia sources; never project the whole Gaia subset through the WCS** (use
  `NeighborIndex`); **FLAG bit 21 fires on the comet**; the FLAG plane misses ~16 % NaN pixels.
- Diagonal errors add 5–11 "detections" but 3–6× the Gaussian negative tail — GLS stays; an S/N
  score over flagged channels picks the most contaminated aperture (2024 N1); census
  differences of ±1–2 between methods are noise.  BIC on raw χ² calls every kink decisive
  (inflate to χ²_red = 1 first); an outburst test against the recent *median* opens on any
  steep rise (compare with the extrapolated trend).
- `test_measure_photometry_actually_refines_the_centroid` needs > 3 GB in a bare Linux VM
  (sep/photutils build there) and passes in the `spherex` env.  The T7 drops under sustained
  reads; macOS writes `._` sidecars on exFAT; `pgrep -fc` is invalid on macOS.
