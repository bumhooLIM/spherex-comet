# Claude Instructions for project folder: `spherex-comet`

## Project Context
- **Purpose:** Production of the SPHEREx comet spectrum catalog — gas production rates
  (H₂O, CO₂, CO) of 68 comets from SPHEREx spectrophotometry with the ZTF dust activity
  context — plus the paper review/writing and project management around it.
- **Audience/Users:** Professional academic astronomers, data scientists.
- **Three stages, one tree** (merged on 2026-09-15 from `ztf-comet`, `spherex-comet-apphot`
  and `spherex-comet-catalog`): `ztfcomet` (ZTF query, Afρ photometry, coma profiles,
  heliocentric trends) → `spherex_apphot` (SPHEREx aperture photometry) → `spherex_comspec`
  (phase grouping, continuum subtraction, production-rate fitting).  `doc/summary.md` is the
  scientific overview; read it once.

## Tech Stack & Dependencies
- **Languages:** Python ≥ 3.10 (`.py`, `.ipynb`).  **Conda environment:** `spherex`
  (`~/miniconda3/envs/spherex/bin/python`).
- **Core libraries:** numpy, pandas, scipy, astropy, astroquery, photutils, sep, matplotlib,
  requests, tqdm, pyarrow, openpyxl; python-pptx + Pillow for the review deck.  All declared in
  `pyproject.toml` (one project, three packages).  `pip install -e .` once, or run from the
  project root — every driver and notebook puts the root on `sys.path`.
- **Never add a third-party dependency without asking for confirmation first.**

## Coding Conventions
- **Documentation:** NumPy-style docstrings.  Document the *why* of complex algorithms, not
  just the *what* — especially where a choice guards against a past bug.
- **Paths come from the package `directory.py` modules** (`ztfcomet.directory`,
  `spherex_apphot.directory`, `spherex_comspec.directory`), each with environment overrides.
  Never hardcode a path or derive one from `Path.cwd()`; scripts locate the root with
  `Path(__file__).resolve().parents[2]`, notebooks walk up to the package.
- **Configuration:** per-target constants in `ztfcomet/config.py`; every SPHEREx tunable in
  `spherex_apphot.config.Config` (frozen, hashed into every `.meta.json`) and
  `spherex_comspec.config` (variants, windows, placeholders).  Never paste constants into a
  notebook.
- **rcparams:** figures of the SPHEREx stages import `notebooks/rcparams.py` for its side
  effects (base `font.size` **20**, `axes.titlesize` 22, `axes.labelsize` 20, tick labels 18,
  `legend.fontsize` 16, `figure.titlesize` 26, bold weights).  Do not override them downward
  to make a dense figure fit — enlarge `figsize` (a 3-panel row wants ~`(24, 7)`) and scale
  markers and hard-coded `fontsize=` to match.  `ztfcomet.rcparams` still carries the older
  15-pt style its figures were tuned at; unifying the two is an open housekeeping item.
- `savefig.dpi` 200 for one-off figures; **50** when writing hundreds of files in a batch.
- **Notebooks validate, packages compute.**  Science logic belongs in the packages so that
  what runs in batch is what the notebooks validate; study scripts live in `notebooks/<stage>/`.
- **Products are namespaced by stage**: `results/{ztf,apphot,comspec}/`, `fig/{ztf,apphot,comspec}/`,
  `doc/{ztf,apphot,comspec}/`.  The SPHEREx photometry lives once, in
  `results/apphot/photometry/`, and stage 3 reads it there — never copy it.

## Architecture & Structure
- `ztfcomet/`, `spherex_apphot/`, `spherex_comspec/`: the packages (each has a `README.md`).
- `scripts/{ztf,apphot,comspec}/`: batch drivers and standalone processing scripts.
- `notebooks/`: `rcparams.py` + `{ztf,apphot,comspec}/` validation notebooks and study scripts.
- `tests/`: the three offline suites (`pytest tests/`), regression tests pinning fixed defects.
- `data/`: `reference/` (working list, orbit classes, literature Q tables), `fluorescence/`
  (the g-factor database), `comspec/` (phase assignment, point spectra, study intermediates),
  `spherex_sample/` (2P + 24P cutouts for the end-to-end test), `ztf/` (local fallback).
  Raw data live on the SSD: `/Volumes/T7/data/spherex-comet`, `/Volumes/T7/data/ztf-comet`,
  Gaia DR3 in `~/Desktop/data/gaia_dr3`.
- `results/`: main clean products per stage; `results/comspec/gas_fit.csv` is **the** result.
- `fig/`: all figures per stage.  `doc/`: `summary.md`, `README.md` index, per-stage notes,
  `literature/`, `figures.pptx`.  `_archive/`: superseded code kept for the code reviews only —
  never read it for current numbers.
- A `README.md` in every one-step subfolder gives the detailed instructions for that folder.

## Scientific Invariants (do not regress)
Fixed bugs with regression tests; changing them silently corrupts results.
- **ZTF** — the zeropoint is per frame (`MAGZP` varies by > 3 mag across a run); ephemerides
  join on JD, never by row; downloads are validated for the FITS signature; frames are
  flagged, never dropped (`quality_ok` = no critical flag); Afρ is aperture-dependent, always
  quote `rho_km`; Horizons record numbers are never hardcoded and fragments are excluded
  unless `allow_fragment=True`; `sep.winpos` returns three values; Afρ figures plot r_h − q,
  never signed r_h; profiles are normalised to the mean SB within 1.5 px, the comet
  oversampled ×4, and centred on the optocentre; a long run must not fall back to the local
  `data/ztf` (`survey.py` refuses).
- **SPHEREx photometry** — the data unit is mJy per pixel (the MJy/sr note is outdated);
  `xcen`/`ycen` are used as given, no centroiding; negative fluxes are kept (`abmag` is NaN
  there); bad pixels are masked with an effective-area correction and `badphot` means *a bad
  pixel in the aperture*, never a flux sign; Gaia sources are **never masked**, only flagged
  (`sourceflag` a > b > c > d > 0); the sky term is not double-counted (`sky_noise_mode="level"`);
  FLAG bit 21 fires on the comet and must never be masked; the annulus is physical
  (150 000 km, floored 15 px, capped 40 px); `jd_utc`/`jd_tdb` at full precision.
- **Gas fits** — one fixed aperture per phase (20 000 km inside 3 au, 40 000 beyond, 60 000
  under 1.5 px); rows enter unless `frac_badpix_ap > 0.05` or `sourceflag = a`; the continuum
  is fitted in distance-corrected flux and divided back per channel; GLS with the continuum
  covariance; `detected` ≥ 3σ, `marginal` 1–3σ, robust = n_eff ≥ 2; no negative-channel cut;
  hot-band H₂O only inside 3 au; coverage is judged before significance (*not covered* ≠
  *not detected*); band shapes and g-factors from `data/fluorescence/`, g(CO) at the comet's
  v_h.  Emission files are keyed on (target, aperture) — `save_emission` splits by `r_ap_km`.
- **Cross-stage** — `arc` in `phase_map.csv` is an arc *index*, not a direction: decide the
  orbital phase from T_p.  After a regroup, rerun `scripts/ztf/afrho_trends.py` first, then
  `scripts/comspec/attach_afrho_ztf.py`; the pipeline blanks moved phases as `stale`.

## Strict Constraints (Do NOT Do These)
- Never introduce new third-party dependencies without asking for confirmation first.
- **No bulk loading:** do not read entire FITS files, the photometry CSVs, `db_filtered.parq`
  or full PDFs into the context window.  Inspect with transient Python snippets that print
  headers, `df.describe()`, shapes or `usecols` slices.
- **Document summaries:** when reading a long paper, write a lean `.summary.md` beside the
  PDF and drop the full text from active memory.
- Never commit FITS, Parquet, notebook outputs, `results/`, `fig/` or `data/` beyond
  `data/reference/` and the READMEs (`.gitignore` enforces it); never mask Gaia sources; never
  read `_archive/` for current numbers; never let a survey run fall back to the local data dir.
- macOS on exFAT writes `._name` sidecars: `COPYFILE_DISABLE=1`, `--clean-appledouble`;
  `pgrep -fc` is invalid on macOS (use `pgrep -f … | wc -l`); no `setsid`.

## State Management & Handoffs
- **No chronological logs.**  Do not append every execution or debugging step anywhere.
- **`handoff.md`** in the root is updated ONLY at natural breakpoints (a feature completed, a
  session closed), max 50 lines: **Current State**, **Next Steps**, **Blind Spots/Dead Ends**.
  Read it at the start of every task to resume context instantly.
- Decisions that change the science go into the stage note (`doc/comspec/pipeline_decisions.md`,
  `doc/ztf/afrho_heliocentric_trends.md`, `doc/apphot/pipeline_upgrade_notes.md`) and into
  `results/comspec/placeholders.csv` / `config.PLACEHOLDERS`; `doc/summary.md` is refreshed
  when the headline numbers move.
- **Git:** one repository at the root, remote `bumhooLIM/spherex-comet`; the history was
  restarted on 2026-09-16 (the earlier merged `ztfcomet` / `spherex-comspec` history, tag
  `spherex-comspec-import` included, is archived in the local backup named in `handoff.md`).
  Commit at breakpoints with a message that says *why*; never push without being asked.

## Literature Review & Manuscript Writing
- **Data extraction:** when parsing PDFs or text in `doc/`, extract methodology, spectral
  features and numerical results relevant to cometary science and SPHEREx; no generic summaries.
- **Academic tone:** the objective, precise tone of *ApJ* / *MNRAS* / *PSJ*.
- **Citations:** never hallucinate citations, authors, years or DOIs.  If a citation is
  needed but not in the loaded context, insert `[CITATION NEEDED: Author/Topic]`.
- **LaTeX/Markdown:** reserve LaTeX for variables and formulas (e.g. $r_h$, $\Delta$,
  $Af\rho$); never wrap ordinary text in equations.
