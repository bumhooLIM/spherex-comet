# Claude Instructions for project folder: `spherex-comet-apphot`

## Project Context
- **Purpose:** Conduct astronomical aperture photometry for the SPHEREx comet cutout dataset.
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
    - `savefig.dpi` is 200 for one-off figures; drop to **50** when writing hundreds of files in a batch.

## Architecture & Structure
- `doc/`: Technical guidebooks
- `data-sample/`: Sample FITS data to be processed in testbed. 
- `results/`: Main, clean results of processing
- `fig/` : All figures
- `notebooks/`: Project notebooks, Unit and integration tests, Standalone data processing, Figure plotting. Mostly `.ipynb` and `.py` extensions.
- `emission-fitter/`: Main project scripts for SPHEREx emission line fitter.
- `README.md` (if exists) in the subfolder describe the detailed instruction when dealing with the corresponding data/files.

## Strict Constraints (Do NOT Do These)
- Never introduce new third-party dependencies without asking for confirmation first.
- **Data Managing:** The file size of raw dataset to be processed is very large and saved separately in out SSD (defined in `notebooks\directory.py`). Do NOT attempt to read entire FITS files, massive CSV catalogs from `data/`, or full PDFs into the active context window. 
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

