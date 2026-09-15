# Project directives for the SPHEREx aperture photometry

_Bumhoo Lim's instructions to the pipeline rewrite (2026-09-06; converted from
`note_to_claude.rtf` on 2026-09-15).  These are the decisions that overruled or
qualified the code review (`code_review_primitive.md`); their disposition is in
`pipeline_upgrade_notes.md`._

## Structure

- Elaborate the overall code set-up into a clean, structured format: individual
  modules such as `__init__`, `config`, `directory`, `phot`, `main`, so the code
  can be run in a clean and understandable way.
- Besides `main.py` as the batch execution, write a Jupyter notebook `main.ipynb`
  that follows and validates the pipeline end-to-end (`notebooks/`).
- All figure plotting belongs in the notebooks.  `main.py` executes and saves
  only the main photometric results and does not plot.
- Generate, update and save a status log in the pipeline.

## Answers to the review's "concerns" (section 4 of `code_review_primitive.md`)

- **4.1** Make the script run as committed.
- **4.2 Scientific correctness**
  - **S1** The FITS data unit is actually mJy/pixel, so no conversion is needed.
    The information about the MJy/sr unit is outdated.
  - **S2** Fix the sky annulus size for all aperture sizes (`r_in = 15 px`,
    `r_out = 20 px`).  Since the SPHEREx central wavelength varies for every
    single pixel, it is important to set the annulus as close to the centre
    pixel as possible.  Aperture radius set-up:
    - fixed km: 2 000 km steps from 10 000 km to 40 000 km, then 60 000,
      80 000, 100 000 km; skip an aperture if its radius is too small
      (< PSF FWHM) or too large (> `r_in`);
    - fixed pixels: 2 px, 5 px.

    _(Superseded on 2026-09-12 by the physical annulus: inner radius 150 000 km
    at the comet, floored at 15 px and capped at 40 px, 5 px wide — see
    `pipeline_upgrade_notes.md` and `doc/comspec/apphot_comparison.md` §6.)_
  - **S4** Do not adjust the photocentre with centroiding; `xcen` and `ycen`
    are already corrected.
  - **S5** Do not drop negative source detections.
  - **S6** Correct the double-counting problem in the sky error estimation.
    Make sure the final error includes not only the pixel-based Poisson noise
    but also the sky fluctuations via `sky_stddev`.
  - **S8** Conduct median stacking in a more robust way.
  - **S11** Conduct the reflectance estimation in a more robust way.
- **4.3, 4.4, 4.5** Adopt all of them to improve code robustness.

## Data masking and flagging

There are two types of masking, `flag_to_mask` and `gaia_to_mask`, and they are
treated differently.

1. **Flag masking.**  A flagged pixel is literally a "bad" pixel.  Treat these
   pixels as `badpix` and mask them out: during photometry ignore the pixel as
   if it were not part of the photometric calculation.  If one bad pixel falls
   inside a photometric aperture (`r_ap`), ignore it, estimate the total flux by
   summing all other pixels, and treat the effective aperture area as reduced
   by one pixel.  Record the number of bad pixels and mark `badphot` only if a
   bad pixel is inside the photometric aperture (not in the sky annulus).

2. **Source masking.**  Do not build a data mask for Gaia sources and do not
   flag out the pixels with Gaia contamination.  Naively conduct the aperture
   photometry without considering the source contamination, and instead record
   the `sourceflag` of the measurement:
   - `0`: no corresponding flag detected;
   - `a`: the photometric aperture + 3 PSF FWHM from the centre contains any
     *bright* Gaia source, bright meaning an effective G < 15 mag.  The effective
     G magnitude sums all Gaia fluxes within the region (two G = 15 sources give
     G_eff ≈ 14.25);
   - `b`: the photometric aperture + PSF FWHM contains Gaia sources with
     V_mag > 0.2 × G_eff, where V is the expected coma brightness from
     JPL/Horizons, already given as `vmag` in `db_filtered.parq`;
   - `c`: the photometric aperture + PSF FWHM contains any Gaia source;
   - `d`: S/N < 1.

   Record the columns `gmag_brightest` (brightest Gaia magnitude within
   `r_ap_eff = r_ap + PSF_FWHM`), `gmag_eff` (effective Gaia magnitude within
   `r_ap_eff`), `n_gaia` (number of Gaia sources within `r_ap_eff`),
   `gmag_nearest`, `dist_gmag_nearest` (pixel distance to the nearest Gaia
   source) and `sourceflag`, with priority `a` > `b` > `c` > `d` > `0` — if
   both `a` and `c` are true, record only `a`.

   _(As implemented: `a` = combined G inside `r_ap + 2 FWHM` brighter than 13;
   `b` = Gaia flux inside `r_ap + FWHM` above 20 % of the comet's flux, i.e.
   `G_eff < vmag + 1.75`; see `spherex_apphot/README.md`.)_
