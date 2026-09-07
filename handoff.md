# Handoff

## Current State
`ztfcomet` v0.2.0 — the merge of `ztf-sso-query` + `ztf-comet` is complete and working.

- Package: `directory`, `config`, `query`, `cutout`, `phot`, `plotting`, `rcparams`.
- `notebooks/main.py` runs query → download → photometry → figures for one or many targets.
- Validation notebooks: `query.ipynb`, `afrho.ipynb`, `figure.ipynb` (test target 24P).
- 38 offline tests pass; `test_regression_c*` pin the reviewed defects.
- Verified end-to-end on the 111-frame C/2024 E1 sample: 91/111 unflagged,
  A(0)fρ = 448–1552 cm at ρ = 15000 km, g and r tracking each other.
- Review: `doc/primitive_code_analysis.md` (22 concerns, C1–C22).

Four photometric corrections applied vs. the old notebook: per-frame zeropoint (C1),
colour term (C4), aperture correction (C5), uncertainties through to Afρ (C8).

## Next Steps
1. **`pytest` is not installed in `spherex`.** Tests currently run only via a shim.
   `conda install -n spherex pytest`, then `pytest tests/ -q`.
2. Re-reduce 240P and 24P with the fixed pipeline — every pre-merge Afρ value is
   invalid (C1). Old figures in `notebooks/legacy/afrho_240P.ipynb` must not be reused.
3. Repair the 2 corrupt downloads in the 2024E1 set: rerun with `--overwrite` or
   let `download_urls(..., repair=True)` handle them; both URLs serve valid FITS today.
4. Consider a per-target `phase_beta` in `config.py` — 0.03 mag/deg is a generic
   dust value, and `phase_beta_err` is currently 0 so it contributes no uncertainty.

## Blind Spots / Dead Ends
- **Horizons epoch lists ride in the request URL.** 75 JDs succeed, 100 return
  HTTP 502. `max_epochs_per_call` is 50. Do not raise it.
- **`elong` needs Horizons quantity 23**, not 16. The old notebook only got it by
  calling `ephemerides()` with no quantities and taking service defaults.
- **Do not join Horizons output by row position.** It returns epochs *sorted*, so
  `pd.concat(axis=1)` was right only by accident. Merge on JD (tolerance 1e-5 d;
  Horizons rounds its JDs).
- **IRSA returns HTTP 200 with HTML error bodies.** Status alone is not enough —
  check the FITS signature. And never use bare `if path.exists(): continue` as a
  resume check; that made corrupt files permanent in the old code.
- **`ccdproc.ImageFileCollection` silently skips unreadable files**, hiding bad
  downloads. Replaced by `phot.build_frame_table`, which reports them.
- **`np.sqrt(data/gain + ...)` NaNs on negative pixels** — 97.5% of one real frame.
  Clip the Poisson term at zero and flag instead.
- **Do not hardcode N/E as (0,−1)/(−1,0).** ZTF frames carry up to ~5° rotation;
  `plotting._sky_axes` derives them from the WCS.
- **`flag_color_default` must stay advisory.** Folding it into `quality_ok` rejected
  9 good single-band frames for a ~0.015 mag systematic.
