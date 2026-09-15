# `data/`

## `apphot/` — revised aperture photometry (input, never written to)

One CSV per comet (`<designation>.csv`, 68 comets, 27 797 exposures) plus a
`.meta.json` with the `spherex_apphot` configuration (`config_hash`
`5502194856bc`).  One row per (exposure, aperture).  Fluxes are in mJy; the
photometry was made on the SPHEREx Level-2 images with the predicted ephemeris
position (no re-centring), a sky annulus whose inner radius is 150 000 km at the comet
(floored at 15 px, capped at 40 px inside the 91-px cutout, 5 px wide; since 2026-09-12 —
the fixed 15–20 px ring used before sat inside the coma of close comets), bad pixels masked
with an effective-area correction, and Gaia sources flagged rather than masked.  `jd_utc`
and `jd_tdb` are written at full precision since 2026-09-12; the 2026-09-08 set (rounded to
0.1 d, near annulus, config hash `5502194856bc`) is kept in `_archive/apphot_5502194856bc/`.

| group | columns |
|---|---|
| identity | `filename`, `objdesig`, `obsid`, `detector`, `epoch` (28-day observing block), `date_obs`, `jd_utc`, `jd_tdb` (full precision since 2026-09-12) |
| channel | `wl`, `wlwidth` (µm), `sun_jy` |
| aperture | `ap_label` (`km20000`, `pix05.0`, …), `ap_kind`, `r_ap_pix`, `r_ap_km`, `r_ap_arcsec`, `r_in_pix`, `r_out_pix`, `pixel_scale_km`, `psf_fwhm_pix`, `aperture_area_pix2`, `aperture_area_eff_pix2` |
| bad pixels | `n_badpix_ap`, `badpix_area_ap_pix2`, `frac_badpix_ap`, `badphot` (any bad pixel in the aperture — *not* a flux-sign cut) |
| sky | `sky_median_mjy_per_pix`, `sky_std_mjy_per_pix`, `sky_area_pix2`, `n_badpix_sky`, `sky_excess_ratio` (annulus scatter / variance-plane prediction) |
| flux | `aperture_sum_mjy`, `source_sum_mjy`, `err_pix_mjy`, `err_skylevel_mjy`, `err_skyscatter_mjy`, `source_sum_err_mjy` (formal), **`source_sum_err_empirical_mjy`** (used by the pipeline), `snr`, `abmag`, `abmag_err` |
| distance correction | `flux_distcorr_mjy` = `source_sum_mjy` × r_h² × Δ², `flux_distcorr_err_mjy`, `distcorr_factor` |
| Gaia | `n_gaia`, `gmag_brightest`, `gmag_eff`, `gmag_eff_bright`, `gmag_nearest`, `dist_gmag_nearest`, `sourceflag` (`a` G < 13 star within r_ap + 2 FWHM; `b` Gaia flux > 0.2 × comet flux within r_ap + FWHM; `c` any source; `d` S/N < 1; `0` clean) |
| geometry | `r_hel`, `r_obs`, `alpha`, ecliptic longitudes/latitudes, rates, `sky_motion`, `sky_motion_pa`, `vmag`, `xcen`, `ycen`, `ltv1`, `ltv2`, `cutout_size`, `pix_scale`, `ra`, `dec` |

The pipeline analyses one aperture per comet — since 2026-09-12 the smallest km aperture
(present for ≥ 95 % of the exposures, ≥ 2 PSF FWHM, ≤ ⅓ of the annulus inner radius) within
10 % of the best median emission-window S/N (`results/apertures.csv`) — and drops a row
only when `frac_badpix_ap > 0.05` or `sourceflag = a`.

## `phase_assignment.csv`, `emission/`

Written by the pipeline: the phase (single observing state) of every exposure,
and per comet the continuum summary (`<comet>_<aperture>.csv`) and the
continuum-subtracted point spectrum (`<comet>_<aperture>_points.csv`) of the main
variant.  `studies/<variant>/emission/` holds the same for the study variants.

## `reference/`

`comet_orbit_classes.csv` (dynamical classes), `gicquel2023_neowise_production_rates.csv`
and `harrington_pinto2022_gas_production_rates.csv` (literature Q values, built by
`notebooks/build_*_table.py`), `phase_update_map_previous.csv` (the phase grouping of
the previous pipeline, used by the regression test of the regrouping).

## `fluorescence/` — reconstructed GSFC-style fluorescence g-factors (0.7–5.0 µm)

Line lists, band g-factors and spectral profiles of H2O, CO2, CO and eleven organic / N- / S-bearing
species at five rotational temperatures, rebuilt from HITRAN 2020 and public solar spectra with the
General Fluorescence Model of Villanueva et al.; validated against the published values to ~10 %.
Formats in `fluorescence/README.md`, method in `../doc/fluorescence_database.md`, builder
`../notebooks/fluorescence_gfm/build_fluorescence_db.py`. Used by the pipeline since 2026-09-11 (`spherex-comspec/spherex_comspec/fluorescence.py`): the species templates of the fit and the velocity-dependent g(CO) come from `profiles/` and `co_swings.csv`.
