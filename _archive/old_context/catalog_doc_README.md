# `doc/` — index

| document | read it for |
|---|---|
| `pipeline_decisions.md` | **start here.**  The pipeline as it stands, the alternatives that were tried (flag policies, `badphot` rule, continuum flux space, error column, previous photometry), why the current configuration is the most robust one, the concerns found, and the placeholders that remain, in priority order. |
| `model_concept.md` | the physics of the synthetic coma emission model: fluorescence, Haser coma, aperture integration, opacity, instrument convolution. |
| `fitting_methodology.md` | the formalism of the fitter layer by layer (geometry → column density → excitation → spectrum → instrument → linear inverse problem), with the verification that was performed. The code it describes lives in `spherex-comspec/spherex_comspec/{gasmodel,instrument,fitting}.py`. |
| `fluorescence_database.md` | the physics of solar-pumped fluorescence (pumps, cascade, Swings effect, spectroscopic inputs) and the procedure used to rebuild the GSFC fluorescence line database for 14 species over 0.7–5.0 µm from HITRAN and public solar spectra, with its validation against the published g-factors. Products in `../data/fluorescence/`, code in `../notebooks/fluorescence_gfm/`. |
| `apphot_comparison.md` | what changed between the previous and the revised aperture photometry inside the emission windows, and what it did to the band fluxes and Q. |
| `literature/` | the papers behind the model (Ootsubo et al. 2012 and the fluorescence literature) with their `.summary.md` files; `literature/fluorescence-db/` holds the five Villanueva et al. papers that define the GSFC fluorescence models. |

The tool guide (installation, stages, variants, output columns) is
`../spherex-comspec/README.md`; the working state is `../handoff.md`.
