# handoff — 2026-09-12

## Current state
Package `spherex-apphot/spherex_apphot/`, driven by `main.py` (photometry only),
`make_figures.py`, `extract_targets.py`, `build_gaia_cache.py`.
`notebooks/main.ipynb` validates the same code end to end. Prototype in
`legacy/`, not runnable. Disposition: `doc/pipeline_upgrade_notes.md`. Tests 34/34.

**Production run 2026-09-12, config hash `9fcc7ea3871a`** — the 68 comets of
`doc/sx_comet_list_ver2607.xlsx`: 27 797 exposures → **455 245 rows, 68/68 ok**,
25.3 min at 4 workers (`--force`). Two changes from the 2026-09-06 run
(`5502194856bc`, 447 942 rows, archived in the catalog project as
`spherex-comet-catalog/_archive/apphot_5502194856bc/`):
1. **Sky annulus is physical**: inner radius 150 000 km at the comet
   (`Config.annulus_r_in_km`), floored at 15 px and capped at 40 px inside the
   91-px cutout, 5 px wide (`apertures.annulus_radii_pix`; `r_in_pix`/`r_out_pix`
   per row). The fixed 15–20 px ring sat at 40 000 km for 24P, inside the coma,
   and removed 3.5–9 % of the clean-channel flux (catalog `doc/apphot_comparison.md`
   §6). Close comets gain 2–6 % flux at 20 000 km and more valid large apertures;
   comets beyond ~1.6 au keep the near ring. `annulus_r_in_km=None` restores the
   fixed ring.
2. **`jd_utc`/`jd_tdb` at full precision** (`Config.jd_float_format="%.9f"`,
   `pipeline.format_jd_columns`); `"%.8g"` had rounded them to 0.1 day.

Staging on T7 (exFAT, `comets_v5_fits_filtered/<target>/`, 27 797 files) unchanged.
Outputs: `results/apphot/<slug>.{csv,meta.json}`, `results/status.csv`, `results/logs/`;
`fig/` still holds the 2026-09-07 figures (cutouts not regenerated: 8 GB).

## Next steps
1. `make_figures.py --target-list doc/sx_comet_list_ver2607.xlsx --workers 6 --no-cutouts`
   to refresh the per-target summaries for the new run (cutouts only on request).
2. **Aperture correction (review S3) still open** — the catalog now requires
   `r_ap >= 2 PSF FWHM`, which sidesteps it for the analysed apertures only.
3. Confirm `config.FLAG_BIT_MEANINGS` against the SPHEREx Explanatory Supplement.
4. The physical annulus radius (150 000 km) is a convention: test 100 000 and
   200 000 km on the bright comets (24P, 306P, 508P) and compare the band fluxes.

## Blind spots / dead ends
- **Do not mask Gaia sources** — the prototype's discs deleted real coma flux.
- **Do not project the whole Gaia subset through the WCS**: use `NeighborIndex`.
- **macOS writes `._name` sidecars onto exFAT**: `COPYFILE_DISABLE=1`; `--clean-appledouble`.
- **`pgrep -fc` is invalid on macOS**; use `pgrep -f … | wc -l`. No `setsid`.
- **Sigma clipping does not save a mean stack** — `stack_combine="median"`.
- **photutils suffixes columns for any *list* of apertures, including length 1.**
- **LaTeX in non-raw f-strings breaks** (`\r`): use `rf"..."`.
- **The FLAG plane misses bad pixels** (~16 % NaN science, no flag bit). **Bit 21 fires on the comet.**
- **`sky_excess_ratio` ~1.2**: VARIANCE under-reports scatter; `source_sum_err_mjy` is a lower bound.
- A ring beyond 45 px leaves the cutout: the caps in `Config` are not tunable past it.
- `notebooks/rcparams.py` disagrees with `CLAUDE.md` (font.size 15 vs 20); left untouched.
