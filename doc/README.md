# `doc/` — index

**Start with [`summary.md`](summary.md)**: the scientific background, the concepts, the
three-stage pipeline, the main results and the open items in one document.  The review deck
is `figures.pptx` (three slides per fitted phase; `scripts/comspec/make_figure_slides.py`).

## `comspec/` — the catalog stage (gas production rates)

| document | read it for |
|---|---|
| `pipeline_decisions.md` | **the pipeline as it stands**, every alternative tried (flag policies, `badphot`, flux space, error column, windows, apertures, the method matrix), why the current configuration is the most robust, the concerns, the placeholders in priority order |
| `model_concept.md` | the physics of the synthetic coma emission model: fluorescence, Haser coma, aperture integration, opacity, instrument convolution |
| `fitting_methodology.md` | the formalism of the fitter layer by layer (geometry → column density → excitation → spectrum → instrument → linear inverse problem) and its verification; the code is `spherex_comspec/{gasmodel,instrument,fitting}.py` |
| `fluorescence_database.md` | solar-pumped fluorescence (pumps, cascade, Swings effect) and the reconstruction of the GSFC g-factor database for 14 species over 0.7–5.0 µm, with its validation; products `data/fluorescence/`, code `notebooks/comspec/fluorescence_gfm/` |
| `apphot_comparison.md` | what changed between the previous and the revised aperture photometry inside the emission windows, and what it did to the band fluxes and Q |
| `comspec/case_revisions.md` | the review memo of 2026-09-15 (`notes_ver260915.xlsx`) applied case by case: what each directive became in `spherex_comspec.revisions` / `GroupingConfig` / `ztfcomet.config`, and what it changed, before → after |

## `apphot/` — the SPHEREx photometry stage

| document | read it for |
|---|---|
| `pipeline_upgrade_notes.md` | the disposition of every code-review item and the decisions that overruled the review; the measured effects (error budget, stacking, distance correction, the physical annulus) |
| `project_directives_apphot.md` | the project owner's instructions to the rewrite (units, annulus, no centroiding, negative fluxes, masking vs flagging rules) |
| `code_review_primitive.md` | the read-only review of the single-file prototype (now `_archive/legacy/apphot/`); historical paths |

## `ztf/` — the dust stage

| document | read it for |
|---|---|
| `afrho_heliocentric_trends.md` | method and generated results of the Afρ trend analysis: peaks, broken laws, the 3 au test, outbursts, colour, per-comet indices, and the Afρ at every SPHEREx phase |
| `afrho_aperture_systematics.md` | why Afρ can rise with distance in a small aperture (C/2024 E1): a sky systematic, not coma structure |
| `afrho_background_annulus.md` | what the far sky annulus removes from Afρ (½ ρ/r_med for a 1/ρ coma) and the correction |
| `profile_survey_1rho.md` | the survey-wide test of the 1/ρ coma law: it holds to the edge of coverage; slopes are an S/N effect |
| `profile_resolution_24P.md` | the 24P resolution study that motivated the survey test |
| `survey_summary_68comets.md` | coverage, download completeness, flag statistics and the profile census of the 68-comet ZTF survey |
| `primitive_code_analysis.md` | the review of the pre-merge ZTF projects (now `_archive/legacy/ztf/`); the defects pinned by `tests/test_ztfcomet.py` |

## `literature/`

The papers behind the model with `.summary.md` files beside the PDFs where one was written:
`fluorescence-emission/` (Crovisier & Encrenaz 1983, Ootsubo et al. 2012, Debout et al. 2016,
Gicquel et al. 2023, Harrington Pinto et al. 2022, Villanueva et al. 2011, 2018, the PSG
handbook), `fluorescence-db/` (the five Villanueva et al. 2011–2013 papers that define the
GSFC fluorescence models), `dust-continuum/` (Harker et al. 2002, 2007, 2023; Kelley et al.
2016), `ice-absorption/` (Protopapa et al. 2014).  When reading one, write a lean `.summary.md`
next to it rather than loading the full text into a session.

The tool guides are the package READMEs (`ztfcomet/`, `spherex_apphot/`,
`spherex_comspec/`); the working state is `../handoff.md`; conventions are in `../CLAUDE.md`.
