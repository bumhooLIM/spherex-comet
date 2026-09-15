# Claude Instructions for project folder: `spherex-comet-catalog`

## Project Context
- **Purpose:** Data analysis, paper review/writing, and project managing of SPHEREx-observed comet spectrum catalog.
- **Audience/Users:** Professional academic astronomers, data scientists

## Tech Stack & Dependencies
- **Languages:** Mainly Python 3.10
- **File Extenstions:** .py, .ipynb
- **Virtual Conda Environment:** `spherex`
- **Core Libraries:** numpy, pandas, astropy, photutils

## Coding Conventions
- **Documentation:** Use NumPy-style docstrings. Document the *why* of complex algorithms, not just the *what*.
- **rcparams:** For all figures, follow the rcparams set up saved in `notebooks/rcparams.py`. Import it (`import rcparams`) rather than restating the settings.
- **Figure font sizes (default, set in `notebooks/rcparams.py`):** base `font.size` **20**, `axes.titlesize` 22, `axes.labelsize` 20, `xtick/ytick.labelsize` 18, `legend.fontsize` 16, `figure.titlesize` 26. Weights are bold (`font.weight`, `axes.titleweight`, `axes.labelweight`).
    - These were doubled on 2026-09-01 from the ~10 pt that notebooks had been rendering at.
    - **Do not override them downward** to make a dense multi-panel figure fit. Enlarge `figsize` instead, and scale marker sizes and any hard-coded `fontsize=` arguments to match — a 3-panel row wants roughly `figsize=(24, 7)`.
    - `savefig.dpi` is 200 for one-off figures; drop to **50** when writing hundreds of files in a batch (`notebooks/figures.ipynb` writes ~840 PNGs and stays under 80 MB at 50 dpi, vs 240 MB at 110).

## Architecture & Structure
- `spherex-comspec/`: **the pipeline** (`spherex_comspec` package, its own git repo): phase grouping → continuum subtraction → production-rate fitting → cross-variant studies → figures. `main.py all` runs everything. Tool guide in its `README.md`.
- `data/apphot/`: revised aperture photometry (input; never written to). `data/reference/`: orbit classes, literature Q tables, the previous phase map. `data/emission/`, `data/phase_assignment.csv`: pipeline intermediates.
- `results/`: **`gas_fit.csv` is the main result**; continuum, skipped/not-fitted, phase and aperture tables beside it; `results/studies/` holds every study variant and cross-variant table one level down.
- `fig/`: figures of the main variant (`emission_model/`, `cont_subtract/`, `phase_group/`, `summary_*.png`); `fig/studies/` the rest.
- `doc/`: `README.md` index → `pipeline_decisions.md` (start here), `model_concept.md`, `fitting_methodology.md`, `apphot_comparison.md`, `literature/`.
- `notebooks/`: `rcparams.py`, standalone study scripts (`apphot_comparison*.py`), literature-table builders.
- `_archive/`: everything superseded (previous photometry, notebooks, earlier runs). Never read it for current numbers; safe to delete.
- Keep the layout flat: main products at the top of `data/`, `results/`, `fig/`; anything experimental under `studies/`.

## Strict Constraints (Do NOT Do These)
- Never introduce new third-party dependencies without asking for confirmation first.
- **No Bulk Loading:** Do NOT attempt to read entire FITS files, massive CSV catalogs from `data/`, or full PDFs into the active context window. 
- **Data Inspection:** To examine large datasets, write and execute transient Python snippets to print headers, `df.describe()`, or array shapes instead of loading the data directly into chat.
- **Document Summaries:** When reading long papers, output a lean `.summary.md` file next to the original document and drop the full text from active memory.

## State Management & Handoffs
- **No Chronological Logs:** Do NOT append every execution, file modification, or debugging step to a running log.
- **Use `handoff.md`:** Update the `handoff.md` file in the root directory ONLY at natural breakpoints (e.g., completing a feature, closing a session).
- **Handoff Structure (Max 50 lines):**
  1. **Current State:** What is currently working (e.g., "Photometry pipeline extracts fluxes successfully").
  2. **Next Steps:** The exact immediate task for the next session.
  3. **Blind Spots/Dead Ends:** Specific bugs fixed or failed approaches (e.g., "Do not use `scipy.optimize.curve_fit` for the background model; it fails to converge. Use Astropy's `Background2D` instead").
- Always read `handoff.md` at the start of a new task to resume context instantly.

## Literature Review & Manuscript Writing
- **Data Extraction:** When parsing PDFs or text in the `doc/` directory, focus strictly on extracting methodology, spectral features, and numerical results relevant to cometary science and SPHEREx. Do not generate generic, high-level summaries.
- **Academic Tone:** When drafting manuscript sections, use the objective, precise tone standard for journals like *ApJ* or *MNRAS*. 
- **Citations:** Never hallucinate citations, authors, publication years, or DOIs. If a citation is required but not present in the loaded context, insert a distinct placeholder (e.g., `[CITATION NEEDED: Author/Topic]`).
- **LaTeX/Markdown:** Maintain clean formatting. Do not wrap standard text in LaTeX equations; reserve LaTeX strictly for variables (e.g., $R_h$, $\Delta$) and mathematical formulas.

