# `fig/`

All figures, namespaced by stage; gitignored (regenerable).  Style: `notebooks/rcparams.py`
(SPHEREx stages) and `ztfcomet.rcparams` (ZTF); `savefig.dpi` 200 for one-off figures, 50 in
batch.

## `comspec/` (`scripts/comspec/main.py figures`, main variant only)
`emission_model/<stem>.png` one per fit (data, continuum, the three species, residuals);
`cont_subtract/<stem>_validation.png` and `_raw.png` per group and band;
`phase_group/<comet>_phase_group.png` the epoch → phase cuts; `phase_images/<stem>.png` the
stacked ZTF and SPHEREx images per fitted phase (`phase_images.py`); `summary_Q_vs_rhel.png`,
`summary_mixing_ratios.png`; `fluorescence_profiles.png` (the g-factor database);
`studies/` the cross-variant figures (`flag_policy_comparison.png`, `distcorr_effect.png`),
`studies/<variant>/` per-group figures of a study variant (`--study-figs`),
`studies/apphot_comparison/`.  `<stem> = <comet>_<aperture>km_ph<phase>`.

## `apphot/` (`scripts/apphot/make_figures.py`)
`<slug>/`: `spectrum_<slug>_epoch<N>.png` (flux by source flag), `reflectance_*`,
`radialprofile_*`, `stacks_*`, `contamination_<slug>.png`, `sourceflags_<slug>.png`,
`apertures_<slug>.png` (aperture acceptance); `flag_effectiveness_survey.png`.  The
per-exposure cutout PNGs (`<slug>/cutouts/`, two per exposure, ~8 GB for the survey) are
not kept — `make_figures.py` without `--no-cutouts` regenerates them on request.
**The current summaries date from the 2026-09-07 photometry**; regenerate for `9fcc7ea3871a`.

## `ztf/` (`scripts/ztf/`)
`afrho/trend/<comet>.png` (fits, peak, breaks, outbursts, colour, SPHEREx windows, the Afρ
at each SPHEREx phase in red), `afrho/trend_slide/<comet>_S<k>.png` (one per SPHEREx phase,
for the review deck), `afrho/rh/`, `afrho/apertures/`, `afrho/lightcurve/`, `afrho/colour/`,
`afrho/survey_overview.png`, `afrho/systematics_2024E1.png`; `photometry/<comet>/cutout/`
annotated cutouts per frame (7 783 PNGs); `profile/<comet>/` the radial profile per frame
and `profile/<comet>_summary.png`.
