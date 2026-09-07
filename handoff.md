# Handoff

## Current State
`ztfcomet` v0.3.0 — merged package, plus Gaia contamination flagging and
fragment-safe Horizons resolution.

- `directory`, `config`, `query`, `horizons`, `cutout`, `phot`, `gaia`, `plotting`.
- `notebooks/main.py`: query → download → photometry → figures.
- `query/afrho/figure.ipynb` all execute clean end-to-end (nbconvert-verified).
- 89 offline tests pass (`pytest tests/ -q`). Review: `doc/primitive_code_analysis.md`.
- C/2024 E1 (111 frames): 85 unflagged, A(0)fρ 439–831 cm at ρ=15000 km, ~7% scatter.

**Contamination:** Gaia sources within `rho_pix + FWHM` summed to
`G_eff = -2.5·log10(Σ10^-0.4G)`, flagged at ≥30% of the JPL `Tmag` flux
(`G_eff ≤ Tmag + 1.31`). Catches 5/111 frames — 5 of the 6 highest Afρ points
(contaminated median 1365 cm vs 654 cm clean).

**Fragments:** `240P` → parent per epoch (2018→90001211, 2025→90001212), never
`240P-B`. `get_target("240P-B")` requests the fragment explicitly.

## Next Steps
1. **Publish.** Committed locally, no remote: `git remote add origin <url> && git push -u origin main`.
2. Re-reduce 240P and 24P. Beyond the C1 zeropoint bug, earlier 240P results may
   be of the *fragment* or of 233P/234P via stale record numbers.
3. Repair the 2 corrupt 2024E1 downloads (`--overwrite`); both URLs are fine today.
4. `pip install -e .` not yet run — notebooks reach the package via `sys.path`.
5. Consider per-target `phase_beta`; `phase_beta_err` is 0, so it adds no error.

## Blind Spots / Dead Ends
- **`sep.winpos` returns `(x, y, flag)`.** Unpacking two raised, was swallowed by
  the surrounding `except`, and left every aperture on the raw ephemeris position
  — ~1.5 px (25% of aperture radius) off on 79/111 frames. Symptom: an all-zero
  centroid-shift histogram.
- **Horizons answers questions you did not ask, and the answer looks fine.**
  Three variants, all silent, all caught only by `verify_targetname`:
  (a) a bare designation may not be read as a comet at all — `id="2P"` returns
  **Styx (905)**, a moon of Pluto; always pass `id_type="smallbody"`.
  (b) record numbers are renumbered — `90001203`/`90001204` meant 240P once and
  now mean 233P/234P; never hardcode one.
  (c) a designation can resolve to a *fragment* — `240P` lists `240P-B`, 2.9 mag
  fainter. Orbit epoch matters too: 240P's parent solutions differ by ~50 arcsec
  in 2025, more than the aperture.
- **A silent "success" with zero output is a failure.** After targets moved to
  designations, `search_frames` still read the now-`None` `horizons_id`; every
  query aborted, was handled as "no coverage", and a two-comet run exited 0 with
  no data. Check the counts, not the exit code.
- **Gaia is complete only to G < 18.5** — "uncontaminated" means "no *catalogued*
  source". `gaiadr3_all.npy` (11 GB) is unsorted; use `gaiadr3_deccache/`
  (dec-sorted, 3–6 ms per cone search).
- **`nbformat.validate()` does not catch broken notebooks.** Building `source`
  with `split("\n")` drops newlines, collapsing every cell to one unusable line
  while validation passes. Check that cells *compile*. Notebook `sys.path`
  bootstraps must search upward for `ztfcomet/`. Do not label pixel axes
  "RA"/"Dec" — subplots need `projection=wcs`.
- **Horizons epoch lists ride in the URL**: 75 JDs work, 100 give 502; chunk at 50.
  `elong` needs quantity 23, not 16. Join on JD, never by row position (tol 1e-5 d).
- **IRSA returns HTTP 200 with HTML error bodies.** Check the FITS signature; never
  use bare `if path.exists(): continue` as a resume check. `np.sqrt(data/gain+...)`
  NaNs on negative pixels — clip and flag. `flag_color_default`/`flag_nan_pixels`
  stay advisory: folding them into `quality_ok` rejected 9 good single-band frames.
