# Handoff — `spherex-comet-catalog`
_2026-09-14. Fixed per-phase apertures; continuum/fit methods chosen by a 54-run matrix; full run reproduced._

## Current State
**Main result: `results/gas_fit.csv`** (147 fits over all 68 comets, run hash `a6463b4c610b`,
`spherex_comspec` 1.2.0; tiers `detected` ≥ 3σ, `marginal` 1–3σ, limits at 3σ). **Robust ≥ 3σ
(n_eff ≥ 2): H₂O 17, CO₂ 26, CO 4 → 35 comets, 34 of them "clean" (carrier band PASS: H₂O 17,
CO₂ 25, CO 2); marginal-only a further 19 comets** (H₂O 25 / CO₂ 23 / CO 9 groups).
`python main.py all` = group → 10 variants → analyze → figures, 10.4 min; 22 tests pass.
Input `data/apphot/` unchanged (`spherex_apphot` `9fcc7ea3871a`, 2026-09-12; apphot **not** rerun).

**2026-09-14 rules** (`doc/pipeline_decisions.md` §3, §5, §7.9): one fixed aperture per *phase*
(20 000 km inside 3 au, 40 000 beyond, 60 000 when under 1.5 px; 99 / 73 / 21 phases;
`results/apertures.csv`, `dataio.aperture_table`); windows 2.55–2.80 µm over 2.30–3.00 (H₂O),
CO₂ continuum 3.90–4.65; one-sided extension 1 µm; order by LOO-CV (10 %); GLS; no negative cut;
hot-band cap 3 au. Chosen by `notebooks/method_matrix.py` (`results/studies/method_matrix/matrix.csv`):
the score is clean robust detections, the negative tail (≤ −2σ) the false-positive gauge — diagonal errors add
5–11 "detections" but 3–6× the Gaussian tail, so GLS stays.  Variant `dc_ap_snr` (S/N rule): Q within 1 %.

**Consequences to know:** Q unchanged against 2026-09-12 (H₂O 1.01, CO₂ 1.01, CO 1.00 matched) but the census
tightened (H₂O 20 → 17, CO 8 → 4); χ²_ν median 4.6 (10.4 at 20 000 km: bright comets' LSF); CO at the noise level.

**ZTF dust context (2026-09-14):** `results/afrho_ztf.csv` + the `afrho_*` columns of `gas_fit.csv` /
`phase_map.csv`: r-band A(0°)fρ at each phase's ⟨r_h⟩ (10 000 / 20 000 km, method, note), from
`ztf-comet/notebooks/afrho_trends.py` → `scripts/attach_afrho_ztf.py`; values for 89 / 95 of the 147 fitted
phases (12 comets have no ZTF data).  The pipeline re-attaches them on every rerun, blanking moved phases
(`stale`): **after a regroup, rerun the ZTF script first, then the attach script.**  Review deck
`doc/figures.pptx` (`scripts/make_figure_slides.py`, python-pptx in `spherex`): 3 slides per fitted group —
spectra/fits, the ZTF trend with the phase highlighted (`ztf-comet/fig/afrho/trend_slide/`), and
`fig/phase_images/<stem>.png` (`scripts/phase_images.py`: ZTF r frames nearest the window + SPHEREx cutouts
of the phase in the continuum / H₂O / CO₂ / CO windows, comet-centred north-up medians; `results/phase_stacks/`).

## Next Steps
1. **LSF / band shape for the bright comets** (placeholder 1): χ²_ν 10²–10³ at 20 000 km.
2. Re-test flag `b` with `dc_main_no_b` (120 of 147 fits have flag-b channels); the 2.42 µm band-3/4 step.
3. Quote `marginal` values or limits only?  Which ZTF aperture (10 000 km literature / 20 000 km sampled)?
4. Conventions to vary: aperture thresholds (3 au, 1.5 px, 60 000 km), annulus 100 000 / 200 000 km
   (upstream), Δ tolerance.  Commit `spherex-comspec` (and `spherex-apphot`) to their repos.

## Blind Spots / Dead Ends
- **Emission files are keyed on (target, aperture)**: `save_emission` splits by `r_ap_km`; one write loses the rest.
- **`arc` in `phase_map.csv` is an arc index** (`in` until a resolved perihelion): 217P, 2025 A6 read `in` outbound.
- **Stacks need both T7 roots mounted**; judge a ZTF frame's quality over all its apertures (10 000 km fails far out).
- **Isolated runs**: set `COMSPEC_{APPHOT,DATA,RESULT,FIG}_DIR` before import.  **Windows are module
  constants** (`config.BAND_WINDOWS`): a windows study must monkeypatch them; `dc_rules_previous` shares them.
- **Census differences of ±1–2 detections between methods are noise**; decide on physics among ties.  **An S/N
  score over flagged channels picks the most star-contaminated aperture** (2024 N1).
- **A wide one-sided quadratic fails near perihelion** for v_h (local weighted quartic validated).  **Band
  width drives Q at SPHEREx sampling**; compare integrals, not peaks, when changing profiles/LSF.
- HITRAN labels are not unique (`gfm_core.py`); no PDF library (PDFKit via `osascript`); never read `_archive/`.
