# Claude Instructions for project folder: `ztf-comet`

## Project Context
- **Purpose:** Retrieve and analyse FITS data from the Zwicky Transient Facility (ZTF), Palomar Observatory (MPC code: I41), for comet science — cutout querying, Afρ photometry, figures.
- **Audience/Users:** Professional academic astronomers, data scientists.
- **Repository:** `ztfcomet` — merges the former `ztf-sso-query` (querying) and `ztf-comet` (photometry) projects.

## Tech Stack & Dependencies
- **Languages:** Python >3.10
- **File Extensions:** `.py`, `.ipynb`
- **Virtual Conda Environment:** `spherex`
- **Core Libraries:** numpy, pandas, astropy, astroquery, photutils, sep, matplotlib, requests, tqdm
- Declared in `pyproject.toml`. **Never add a third-party dependency without asking first.**

## Coding Conventions
- **Documentation:** NumPy-style docstrings. Document the *why* of complex algorithms, not just the *what* — especially where a choice guards against a past bug.
- **rcparams:** For all figures follow `ztfcomet/rcparams.py`. Import it (`from ztfcomet import rcparams`) rather than restating the settings.
  - `savefig.dpi` is 200 for one-off figures; drop to **50** when writing hundreds of files in a batch.
- **Directory:** Every path comes from `ztfcomet/directory.py`. Do **not** hardcode paths or derive them from `Path.cwd()` — that is what broke the pre-merge notebooks whenever they ran from anywhere but `notebooks/`.
- **Configuration:** Per-target constants (orbit records, dates, ρ, β) belong in `ztfcomet/config.py`, never pasted into a notebook.

## Architecture & Structure
- `ztfcomet/`: The package. `directory` (paths), `config` (targets and tunables), `query` (Horizons + IRSA), `cutout` (URLs + validated download), `phot` (photometry + Afρ), `plotting`, `rcparams`.
- `notebooks/`: `main.py` (end-to-end driver), `query.ipynb`, `afrho.ipynb`, `figure.ipynb`. Validation and figure work only — **science logic belongs in the package**, so that what runs in batch is what the notebooks validate.
- `notebooks/legacy/`: Pre-merge notebooks, outputs stripped, provenance only. They do not run against this package and contain known bugs. **Do not copy code out of them.**
- `doc/`: Technical guidebooks. `primitive_code_analysis.md` is the review that motivated this structure.
- `data/`, `results/`, `fig/`: Gitignored outputs. Raw FITS live on the SSD (see `directory.py`).
- `tests/`: Offline test suite. `test_regression_c*` pin defects from the review.
- `README.md` (if present) in a subfolder gives detailed instructions for that folder.

## Scientific Invariants (do not regress)
These are fixed bugs with regression tests. Changing them silently corrupts results.
- **Zeropoint is per frame.** `MAGZP` varies by >3 mag across a run. Never reuse a scalar zeropoint across rows.
- **Ephemerides join on JD**, never by row position. Positional joins depend on undocumented Horizons ordering.
- **Downloads are validated.** IRSA returns HTTP 200 with HTML error bodies; a file is only accepted if it carries the FITS signature.
- **Frames are flagged, never dropped.** `quality_ok` covers `CRITICAL_FLAGS` only; advisory flags qualify a measurement without invalidating it.
- **Afρ is aperture-dependent.** Always report `rho_km` alongside a value.
- Horizons epoch lists travel in the request URL: chunk at ≤50 (75 works, 100 returns HTTP 502).

## Strict Constraints (Do NOT Do These)
- Never introduce new third-party dependencies without asking for confirmation first.
- **Data Managing:** Raw datasets are very large and live on the SSD (defined in `ztfcomet/directory.py`). Do NOT read entire FITS files, massive CSV catalogs, or full PDFs into the active context window.
- **Data Inspection:** To examine large datasets, write and execute transient Python snippets that print headers, `df.describe()`, or array shapes instead of loading data into chat.
- **Document Summaries:** When reading long papers, write a lean `.summary.md` next to the original and drop the full text from active memory.
- Never commit FITS files or notebook outputs. The predecessor repository tracked 36 MB of FITS; do not repeat it.

## State Management & Handoffs
- **No Chronological Logs:** Do NOT append every execution, file modification, or debugging step to a running log.
- **Use `handoff.md`:** Update it in the root directory ONLY at natural breakpoints (completing a feature, closing a session).
- **Handoff Structure (Max 50 lines):**
  1. **Current State:** What is currently working.
  2. **Next Steps:** The exact immediate task for the next session.
  3. **Blind Spots/Dead Ends:** Specific bugs fixed or failed approaches.
- Always read `handoff.md` at the start of a new task to resume context instantly.
