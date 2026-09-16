# `spherex_apphot` — SPHEREx comet aperture photometry

Aperture photometry for the SPHEREx comet cutout dataset — the second stage of the
`spherex-comet` project (`ztfcomet` → `spherex_apphot` → `spherex_comspec`).  A rewrite
of the single-file prototype (now in `_archive/legacy/apphot/`) following the review in
[`doc/apphot/code_review_primitive.md`](../doc/apphot/code_review_primitive.md); what changed
and why is in [`doc/apphot/pipeline_upgrade_notes.md`](../doc/apphot/pipeline_upgrade_notes.md).

---

## Layout

```
spherex-comet/
├── scripts/apphot/main.py              batch driver — photometry only, never plots
├── scripts/apphot/make_figures.py      all figures for already-processed targets
├── scripts/apphot/extract_targets.py   build the per-target FITS tree from the flat archive
├── scripts/apphot/build_gaia_cache.py  one-off Gaia declination cache builder
├── scripts/apphot/h5cut2fits.py        the HDF5 cutout archive -> individual FITS (data preparation)
├── notebooks/apphot/main.ipynb         end-to-end validation on the sample data
├── notebooks/apphot/build_db_filtered.ipynb   db.parq -> db_filtered.parq (data preparation)
├── tests/test_apphot_pipeline.py       unit tests (pytest, or run the file directly)
├── _archive/legacy/apphot/             the superseded prototype, kept for the review only
└── spherex_apphot/
    ├── config.py            every tunable number, frozen and hashable
    ├── directory.py         paths, with environment overrides
    ├── logging_utils.py     run logging
    ├── status.py            results/apphot/status.csv and the target-slug rule
    ├── fitsio.py            cutout FITS + Parquet index reading
    ├── wcsutil.py           WCS from a header, or rebuilt from an index row
    ├── masking.py           bad pixels: flags, NaN science, bad variance
    ├── catalog.py           Gaia DR3 access + declination cache + KD-trees
    ├── apertures.py         per-exposure aperture set and validity
    ├── sourceflag.py        Gaia contamination flags (recorded, never masked)
    ├── phot.py              the aperture photometry itself
    ├── epochs.py            grouping exposures into observing epochs
    ├── stacking.py          robust band stacking (sigma-clipped median)
    ├── reflectance.py       reflectance spectra
    ├── diagnostics.py       growth-curve metrics + source-flag effectiveness
    ├── targetlist.py        reading a target list from .xlsx / .csv / .txt
    ├── pipeline.py          orchestration: run_target / stack_target
    └── plotting.py          figures — imported by notebooks only
```

`notebooks/apphot/main.ipynb` walks the whole pipeline end to end and validates it.
It calls the *same* `pipeline.run_target` that `scripts/apphot/main.py` runs in batch,
so the notebook and the production path cannot drift apart.

---

## Quick start

```bash
# one-off, ~5 min, ~8 GB: makes every later Gaia query ~100x faster
python scripts/apphot/build_gaia_cache.py

# one-off per working list: copy just those targets out of the 139k-file
# archive into comets_v5_fits_filtered/<target>/  (exFAT scans that flat
# directory linearly, which is why this is worth doing once)
python scripts/apphot/extract_targets.py --target-list data/reference/sx_comet_list_ver2607.xlsx

# one target
python scripts/apphot/main.py --objdesig 24P

# the whole working list, four at a time, resuming from results/apphot/status.csv
python scripts/apphot/main.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 4

# end-to-end check against the sample data in the repository
python scripts/apphot/main.py --sample --objdesig 2P 24P --stack

# figures for targets already processed (main.py itself never plots)
python scripts/apphot/make_figures.py --target-list data/reference/sx_comet_list_ver2607.xlsx --workers 6

# tests
python -m pytest tests/test_apphot_pipeline.py -q      # or: python tests/test_apphot_pipeline.py
```

Run everything from the project root (`spherex-comet/`).  There is nothing to
install: the drivers and the notebook put the project root on `sys.path`
(`pip install -e .` once does the same for every shell).

## Outputs

| path | contents |
|---|---|
| `results/apphot/photometry/<slug>.csv` | one row per (exposure × aperture) |
| `results/apphot/photometry/<slug>.meta.json` | the full config, counters and epoch table that produced it |
| `results/apphot/status.csv` | one row per target: status, counts, timing, config hash |
| `results/apphot/logs/apphot_<utc>.log` | full run log |
| `results/apphot/stacks/<slug>_epoch<N>_<band>.fits` | band stacks (`--stack`) |
| `fig/apphot/<slug>/` | per-target figures; `fig/apphot/<slug>/cutouts/` two per exposure |

`<slug>` is the designation with spaces removed (`2021 G2` → `2021G2`), defined
once in `status.slugify` and used by both the writer and the resume check.

Read a result with `sourceflag` kept as text:

```python
df = pd.read_csv("results/apphot/photometry/24P.csv", dtype={"sourceflag": str})
```

## Configuration

Everything lives in `spherex_apphot.config.Config`, which is frozen, hashes to a
12-character digest, and is written into every `.meta.json`.

```python
from spherex_apphot import Config, run_target
cfg = Config(sky_noise_mode="empirical", annulus_r_in_km=150000.0)   # None -> fixed r_in_pix/r_out_pix ring
res = run_target("24P", cfg)
```

Paths can be overridden without editing code:
`SPHEREX_T7_DIR`, `SPHEREX_DB`, `SPHEREX_FITS_DIR`, `SPHEREX_GAIA_DIR`,
`SPHEREX_RESULT_DIR`, `SPHEREX_FIG_DIR`.

## What the pipeline does

1. Load the target's rows from `db_filtered.parq` (predicate pushdown).
2. Sort chronologically and split into **epochs** at a 28-day gap.
3. One Gaia cone search covering every pointing; one KD-tree for the target.
4. Per exposure:
   - read `IMAGE` / `VARIANCE` / `FLAG` **by name**;
   - build the bad-pixel mask — flag bits, non-finite science, non-positive variance;
   - resolve the aperture set for this geometry, keeping `PSF_FWHM ≤ r_ap < r_in`;
   - measure the background in the sky annulus (sigma-clipped): inner radius
     **150 000 km** at the comet, floored at 15 px and capped at 40 px so the ring
     stays inside the 91-px cutout, 5 px wide (`Config.annulus_r_in_km`; since
     2026-09-12 — the fixed 15–20 px ring sat inside the coma of close comets and
     removed 3.5–9 % of their flux; `r_in_pix`/`r_out_pix` record the ring per row);
   - measure every valid aperture in one photutils call;
   - record Gaia contamination rather than masking it.
5. Concatenate, write the CSV, sidecar and status row.

Optional: `--stack` writes per-epoch band stacks — background-subtracted,
sub-pixel registered, sigma-clipped mean, with a per-pixel count map.

## Key columns

| column | meaning |
|---|---|
| `source_sum_mjy` | background-subtracted aperture flux [mJy] — **negatives are kept** |
| `source_sum_err_mjy` | uncertainty in the configured `sky_noise_mode` (default: no sky double-count) |
| `source_sum_err_empirical_mjy` | the same from the measured annulus scatter alone |
| `err_pix_mjy`, `err_skylevel_mjy`, `err_skyscatter_mjy` | components, so the budget can be recombined |
| `sky_excess_ratio` | annulus scatter ÷ scatter predicted by the variance plane; **> 1 means the errors are a lower bound** |
| `aperture_area_pix2` / `aperture_area_eff_pix2` | geometric vs bad-pixel-corrected area |
| `n_badpix_ap`, `badphot` | bad pixels *inside the aperture* (annulus ones do not flag) |
| `flux_distcorr_mjy` | flux normalised to `r_hel = r_obs = 1 au`: `source_sum_mjy × r_hel² × r_obs²` |
| `sourceflag` | `a` > `b` > `c` > `d` > `0`, see below |
| `gmag_eff`, `gmag_brightest`, `n_gaia`, `gmag_nearest`, `dist_gmag_nearest` | contamination detail |
| `abmag` | `nan` where the flux is not positive — the magnitude is undefined, the flux is not |

## Source flags

| flag | condition |
|---|---|
| `a` | combined G inside `r_ap + 2·FWHM` brighter than 13 — a bright blend |
| `b` | Gaia flux inside `r_ap + FWHM` exceeds 20 % of the comet's own flux, i.e. `G_eff < vmag + 1.75` |
| `c` | any source inside `r_ap + FWHM` |
| `d` | `SNR < 1` |
| `0` | none of the above |

`a` and `b` are not nested.  `a` is absolute (is there a bright star nearby at
all); `b` is relative (do the blended stars matter compared with *this* comet).
Only the highest-priority flag is stored.

`diagnostics.py` scores the flags against an independent growth-curve symptom of
contamination, per target and across the survey.  Over 67 targets (median 7.9 %
of exposures contaminated):

| | median recall | median lift | median data kept |
|---|---|---|---|
| `cut a` | 0.35 | **5.60** | **0.95** |
| `cut a+b` | 1.00 | 1.72 | 0.44 |

`a` alone beats `a+b` on lift in **56 of 67 targets**.  **Cut on `a`; keep `b`
and `c` as columns**, applying them per target where completeness matters more
than sample size.  `badphot` is *anti*-correlated with contamination (lift 0.8)
and must not be used as a contamination cut — it is about detector pixels, not
neighbours.

    python -c "from spherex_apphot.diagnostics import survey_flag_effectiveness as s; \
               print(s('results/apphot'))"

## Figures

All plotting lives in `plotting.py` and is called only from the notebook.

```bash
python scripts/apphot/make_figures.py --objdesig 2P 24P              # everything
python scripts/apphot/make_figures.py --objdesig 24P --no-cutouts    # summaries only
python scripts/apphot/make_figures.py --target-list LIST --workers 6 # a whole list
```

Two figures per exposure is a lot at survey scale — the 68-comet working list is
27 797 exposures, so ~55 600 per-cutout PNGs and ~8 GB.  Use `--no-cutouts` when
only the per-target summaries are wanted, or `--max-cutouts N` to sample.

or from a notebook:

```python
from spherex_apphot import plotting
plotting.save_cutout_figures(phot, resolver, cfg, gaia_subset=subset,
                             objdesig="24P", dpi=50)     # 2 figures per exposure
plotting.save_target_figures(phot, cfg, objdesig="24P",
                             stacks=stacks, metrics=metrics)
```

Per-cutout: the three-panel view (raw science, FLAG plane, masked science with
every Gaia source marked at a size proportional to its flux) plus the growth
curve.  Per-target: spectrum by flag, reflectance, band stacks, radial profiles,
contamination analysis, flag census, aperture acceptance.

## Caveats that remain open

* **No aperture (encircled-energy) correction.**  Apertures below the PSF FWHM
  are rejected, which removes the worst of it, but the smallest retained
  apertures are still ~1 pixel across against a ~1-pixel PSF, and the SPHEREx PSF
  broadens with wavelength.  A residual wavelength-dependent aperture loss
  therefore remains in every spectrum.  This is review item **S3** and is the
  largest known systematic in the current products.
* **`sky_excess_ratio` is typically 1.1–1.5** in the sample, so the variance
  plane under-reports the real scatter and `source_sum_err_mjy` should be read as
  a lower bound.
* **Gaia positions carry no proper motion** (epoch J2016.0).  Negligible at
  6.2 arcsec pixels, but it is an approximation.
* **Flag-bit meanings** in `config.FLAG_BIT_MEANINGS` were recovered empirically
  by correlating cutout bit counts against the frame-level `N_*` counters.
  Confirm against the SPHEREx Explanatory Supplement before publishing.
