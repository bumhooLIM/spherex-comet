# Pipeline upgrade notes

**Date:** 2026-09-07 · **Package:** `spherex_apphot` 2.2.0 · **Supersedes:** `spherex-apphot/legacy/apphot.py`

This is the disposition record for every item in
[`code_review_primitive.md`](code_review_primitive.md), plus the decisions taken
where the project overruled the review.  Read it alongside
[`../spherex-apphot/README.md`](../spherex-apphot/README.md).

---

## 1. Structure

The 975-line single file became a package with one responsibility per module
(layout in the README).  Three consequences worth naming:

* **`main.py` produces photometry and nothing else.**  No figure is drawn on the
  batch path; `plotting.py` is imported only by `notebooks/main.ipynb`.  Keeping
  the drawing code in a module rather than pasted into the notebook means the
  figures are version-controlled and reviewable, while the notebook stays a
  narrative that calls them.
* **The notebook runs the production code.**  `notebooks/main.ipynb` calls the
  same `pipeline.run_target`.  The prototype's notebook and script were separate
  copies of the same logic and had already diverged.
* **Provenance is written, not remembered.**  Every result carries a
  `.meta.json` with the full `Config`, its hash, per-run counters and the epoch
  table.

## 2. Status logging

`results/status.csv` holds one row per target — status, timing, exposure and row
counts, missing/unreadable files, apertures skipped, `badphot` count, Gaia count,
config hash, output path, error message.  Written atomically (temp file +
`os.replace`), keyed on the slug so a re-run replaces rather than appends, and
marked `running` *before* the work starts so an interrupted run is visible as
such.  `main.py --all` consults it to resume.

Full logs go to `results/logs/apphot_<utc>.log`.

---

## 3. Review dispositions

### 4.1 Blockers — all fixed

| | |
|---|---|
| **B1** `APPHOT_DIR` / `COMBFITS_DIR` undefined | Defined in `directory.py`, along with `LOG_DIR` and `STATUS_CSV`. |
| **B2** `import rcparams` fails outside `notebooks/` | `plotting.apply_rcparams()` puts `notebooks/` on `sys.path` and imports it for its side effects, keeping one source of truth for the style. |
| **B3** no output directory created | `directory.ensure_dirs()` runs at the top of `main.py` and in the notebook. |
| **B4** `apphot_all.py` assumed its own cwd | `main.py` is cwd-independent; it resolves the package from `__file__` and runs targets in-process. |

### 4.2 Scientific correctness

| item | disposition |
|---|---|
| **S1** units | **Overruled by the project: the data is mJy/pixel; the `MJy/sr` note in the exporter is outdated.** No conversion is applied. The assumption is no longer implicit — `Config.flux_unit` is written into every CSV sidecar and every stack header, so the next reader does not have to re-derive it. |
| **S2** annulus geometry | Fixed **`r_in = 15 px`, `r_out = 20 px`** for every aperture and every exposure. It can no longer leave the 91-pixel cutout, and it stays close to the target — which matters because SPHEREx's central wavelength varies across the detector, so a distant annulus samples a different bandpass than the aperture it corrects. Apertures are now 10 000–40 000 km in 2 000 km steps plus 60 000/80 000/100 000 km, and two fixed-pixel apertures (2 px, 5 px). An aperture is dropped when `r_ap < PSF_FWHM` or `r_ap ≥ r_in`. |
| **S3** aperture correction | **Still open.** Rejecting sub-PSF apertures removes the worst of it, but the smallest retained apertures are ~1 pixel against a ~1-pixel PSF and the SPHEREx PSF broadens with wavelength, so a residual wavelength-dependent aperture loss remains in every spectrum. This is the largest known systematic in the current products. |
| **S4** centroiding | **Overruled by the project: `xcen`/`ycen` are already corrected.** Used as given; no centroiding step exists. |
| **S5** negative fluxes | **Kept.** `badphot` no longer looks at flux sign. `abmag` is `nan` where the flux is not positive — the magnitude is undefined there — but the flux and its uncertainty are always reported. |
| **S6** double-counted sky noise | Fixed. Default `sky_noise_mode="level"`: `σ² = σ_pix² + A_eff²·σ_sky²/N_sky`. The variance plane already contains the background's photon noise, so DAOPHOT's `A·σ_sky²` term counted it twice; `σ_sky` still enters through the uncertainty on the *subtracted sky level*. All three components are reported separately, plus `source_sum_err_empirical_mjy` (annulus scatter only) and `sky_excess_ratio`. Measured effect on the sample: the default error is **24 % smaller** than the double-counted form. |
| **S7** flux lost to masked pixels | Addressed by the project's masking rule: the aperture sums good pixels over a correspondingly reduced effective area, no flux is invented, and `badphot` marks it. Both `aperture_area_pix2` and `aperture_area_eff_pix2` are reported so a reader can rescale if they choose. |
| **S8** stacking | Rewritten. Per-frame background subtraction, NaN-safe sub-pixel registration, optional north alignment, sigma-clipped **median** (`Config.stack_combine`, `"mean"` available), and a three-extension output (`IMAGE`/`COUNT`/`ERROR`) with a WCS when — and only when — the frames were rotated to a common orientation. The median is not a cosmetic choice: sigma clipping alone does not remove a star that lands at the same stamp position in a *minority* of frames. With 2 contaminated frames out of 10, a clipped mean returns ~100 where the truth is 0; the median returns ~0 (unit-tested). In a field as crowded as 24P's that is the difference between stacking the comet and stacking its neighbours. The reported `ERROR` plane carries the √(π/2) factor that converts the standard error of the mean into that of the median. |
| **S9** `neargaia_*` computed but unused | Superseded: the Gaia information now drives `sourceflag` and five reported columns. |
| **S10** inconsistent `radius_scale` | Moot: Gaia sources are never masked, so there is no mask radius to be inconsistent about. |
| **new** distance-corrected flux | Every row gains `flux_distcorr_mjy`, `flux_distcorr_err_mjy` and `distcorr_factor`. `Config.distcorr_mode="multiply"` (default) computes `F × r_hel² × r_obs²` — the flux the comet would have shown at `r_hel = r_obs = 1 au`. Reflected flux falls as `1/(r_hel² r_obs²)` (the illuminating sunlight with heliocentric distance, the reflected light again with observer distance), so this exactly cancels the geometric dilution and leaves a quantity constant for an unchanging coma. Measured on 24P: the raw median flux differs by a factor of 64 between its two epochs (r_hel 1.79 → 1.20 au); after correction the ratio is 5.2, and that residual is genuine brightening toward perihelion rather than geometry. `"divide"` is retained for comparison only. |
| **S11** reflectance | Rewritten. `R = F_obs/F_⊙ × r_hel² × r_obs²` — the distance factors were missing, without which epochs at different distances are not on a common scale. Normalisation uses an inverse-variance weighted, sigma-clipped mean over a continuum window (default 1.0–2.5 µm, avoiding the 2.7 µm water band) and propagates the normalisation's own uncertainty. `bin_spectrum()` collapses the single-exposure scatter into a weighted spectrum. An optional linear phase correction exists and is off by default. |

### 4.3 Robustness

| | |
|---|---|
| **R1** resume check vs. filename mismatch | One `status.slugify()` used by both the writer and the check. Covered by a test. |
| **R2** positional HDU access | Extensions looked up by `EXTNAME` (`IMAGE`/`VARIANCE`/`FLAG` with aliases); positional access only as an explicitly logged fallback. |
| **R3** bare `except` returning NaN | Removed. Failures are logged with the filename; the target-level handler records the exception type and message in `status.csv` and the traceback in the log. |
| **R4** unguarded `fits.open` in the main loop | A missing or unreadable cutout costs one exposure; counts land in `n_missing_fits` / `n_read_errors`. |
| **R5** NaN/zero variance unmasked | `mask_bad_variance` rejects non-finite and non-positive variance. In the sample these coincide exactly with the NaN science pixels — which carry **no flag bit at all**, so the flag plane alone was never sufficient. |
| **R6** `int(round(nan))` in the stacker | `extract_stamp` returns an all-NaN stamp for a non-finite centre. |
| **R7** meaningless error fallback | Deleted; the error model is explicit and configurable. |
| **R8** stale index columns | Confirmed a non-issue: the raw `db_filtered.parq` (240 columns) carries neither `phase` nor `r_ap_*`; the prototype notebook added them itself. |
| **R9** magic flag-bit list | `Config.bad_flag_bits`, documented against `FLAG_BIT_MEANINGS`. Bit meanings were **recovered empirically** — per-bit cutout pixel counts correlated against the frame-level `L1/L2 N_*` counters, scaled by the cutout/frame area ratio. Confident: bit 0 = L1 TRANSIENT, 1 = L1 OVERFLOW, 2 = L1 SUR_ERROR, 11 = L2 COLD, 17 = L2 PERSIST, 19 = L2 OUTLIER, 21 = L2 SOURCE. Note bit 21 fires on the comet itself and must never be masked. Verify against the Explanatory Supplement before publishing. |
| **R10** half-pixel `imshow` offset | `extent=(-0.5, nx-0.5, -0.5, ny-0.5)` throughout, so overlaid apertures register with the data. |
| **R11** dead code / contradictions | Gone. |
| **R12** duplicated Gaia rows | One cone search de-duplicates across pointings via the KD-tree. |

### 4.4 Performance

Measured on the sample (2P: 226 exposures; 24P: 765 exposures, 59 427 Gaia sources):

| | |
|---|---|
| **P1** 11.7 GB Gaia read per target | `build_gaia_cache.py` writes a declination-sorted cache; a query binary-searches a memory map and reads only the relevant band. **21 s → 0.2 s per target**, identical results. The cache builder is an external bucket sort, so peak memory is set by one 0.5° bucket rather than by the 419 M-row catalogue. Note a magnitude cut alone would *not* have helped: G < 18 keeps 77 % of this already-limited catalogue. |
| **P2** `O(N_pointings × N_gaia)` intersection | One `cKDTree` on 3-vectors; a second per-target tree (`NeighborIndex`) serves the per-exposure queries. |
| **P3** repeated WCS/mask work | Each exposure is read and masked once for photometry; stacking re-reads only the frames it needs. |
| **P4** unbounded figure sizes | `plot_cutout_grid` caps panels; `plot_coverage` subsamples footprints. |
| **P5** subprocess per target | `main.py --all` runs in-process with one shared catalogue handle; `--workers N` uses a process pool over targets, each memory-mapping the shared cache. |

Whole-target wall time on the sample: 2P ≈ 10 s, 24P ≈ 21 s (765 exposures × up to 21 apertures = 14 526 rows).

### 4.5 Conventions

* **Dependencies.** `skimage` is gone — no Gaia discs are drawn, so
  `skimage.draw.disk` is no longer needed. The package uses numpy, pandas,
  astropy, photutils, scipy and matplotlib. `scipy` is new relative to
  `CLAUDE.md`'s core list; it is used for `cKDTree` and `ndimage` and is already
  present in the `spherex` environment. Flagging it rather than assuming.
* **Naming.** `phase` → `epoch`, so it no longer collides with the solar phase
  angle `alpha` in the same table.
* **Tests.** `tests/test_pipeline.py`, 29 tests, no pytest required. They already
  caught one real bug: photutils suffixes its output columns (`aperture_sum_0`)
  whenever a *list* of apertures is passed — including a list of one — so a
  single-aperture configuration raised `KeyError`.
* **`rcparams` discrepancy (unresolved).** `CLAUDE.md` documents
  `font.size = 20`, `axes.titlesize = 22`, tick labels 18, `figure.titlesize` 26
  and a `legend.fontsize`; `notebooks/rcparams.py` actually sets 15/15/12/20 and
  no `legend.fontsize`. The file was left untouched — changing it would alter
  every existing figure — but one of the two should be corrected.

---

## 4. Masking and flagging, as specified

**Bad pixels are excluded, not masked over.**  A flagged pixel (or a non-finite
science value, or a non-positive variance) is removed from the sum and the
effective aperture area shrinks by exactly its overlap.  Verified to machine
precision: `A_geom − A_eff − badpix_area < 1.2e-13 pix²` across the sample.
`badphot` fires only for bad pixels **inside the photometric aperture**; in the
24P sample 6 228 rows have a bad pixel in the annulus alone and are correctly
left unflagged.

**Gaia sources are recorded, not masked.**  Photometry is performed naively and
each measurement carries `n_gaia`, `gmag_brightest`, `gmag_eff`,
`gmag_eff_bright`, `gmag_nearest`, `dist_gmag_nearest` and one `sourceflag`:

| flag | condition |
|---|---|
| `a` | combined G inside `r_ap + 2·FWHM` brighter than **13** |
| `b` | Gaia flux inside `r_ap + FWHM` exceeds **20 %** of the comet's own predicted flux |
| `c` | any source inside `r_ap + FWHM` |
| `d` | `SNR < 1` |
| `0` | none of the above |

Priority `a` > `b` > `c` > `d` > `0`; only the winner is stored.  Blends are
summed in flux, so two G = 15 sources give `gmag_eff = 14.25` (tested).

`a` and `b` are deliberately **not nested** — they ask different questions.
`a` is absolute: is there a bright star nearby at all, within a generous
`r_ap + 2·FWHM`, whose wings and ghosts contaminate the aperture regardless of
how bright the comet is.  `b` is relative: inside the tighter `r_ap + FWHM`, do
the blended stars matter *compared with this comet*.  `b` compares fluxes,

```
10**(-0.4·G_eff)  >  0.2 · 10**(-0.4·vmag)      i.e.   G_eff < vmag + 1.75
```

which replaces the earlier `vmag > 0.2·gmag_eff` form.  That form compared two
magnitudes on a scale where multiplying one by 0.2 has no photometric meaning,
and for any faint comet it was satisfied whenever *any* source was present — so
`'b'` absorbed everything `'c'` was meant to catch and `'c'` never fired.  With
the flux-ratio form all five flags are populated (24P: a=1890, b=6069, c=1246,
d=392, 0=4929).

---

## 5. Do the flags work?

A flag that never fires is useless; one that fires on everything is equally
useless.  Neither failure is visible unless the flags are tested against an
**independent** symptom of contamination, so `diagnostics.py` builds one from the
growth curve.

Differentiating the growth curve gives each annulus's mean surface brightness.
For a steady-state coma this falls monotonically outward (Σ ∝ 1/ρ); a field star
entering the aperture makes it **rise** instead — a signature no coma model
produces.  `growth_metrics()` measures the most significant such rise per
exposure in units of its own uncertainty, and `flag_effectiveness()` scores each
cut against it.  On synthetic data the separation is complete: a pure 1/ρ coma
gives −2.3 σ, the same coma plus a star gives +107 σ.

`lift` is precision ÷ base rate: above 1 the cut preferentially removes
contaminated exposures, at 1 it is uncorrelated with contamination and is
discarding good data for nothing.

**24P** — bright (V ≈ 18), nearby, Galactic-plane field, 59 427 Gaia sources.
409 usable growth curves, 62 (15.2 %) contaminated:

| cut | removes | recall | precision | lift | purity kept |
|---|---|---|---|---|---|
| `a` | 93 | 0.82 | 0.55 | **3.62** | 0.965 |
| `a+b` | 313 | 1.00 | 0.20 | 1.31 | 1.000 |
| `a+b+c` | 330 | 1.00 | 0.19 | 1.24 | 1.000 |
| `badphot` | 386 | 0.76 | 0.12 | 0.80 | 0.348 |

**2P** — faint (V ≈ 21), distant, sparse field, 1 060 Gaia sources.
215 usable growth curves, 24 (11.2 %) contaminated:

| cut | removes | recall | precision | lift | purity kept |
|---|---|---|---|---|---|
| `a` | 9 | 0.21 | 0.56 | **4.98** | 0.908 |
| `a+b` | 106 | 1.00 | 0.23 | **2.03** | 1.000 |
| `badphot` | 97 | 0.08 | 0.02 | 0.18 | 0.814 |

Reading them:

* **`a` is always precise but its reach depends on the field.**  Lift is ~4–5 in
  both targets, so when it fires it is nearly always right.  In 24P's crowded
  field it also catches most of the contamination (recall 0.82) at a cost of
  23 % of the sample — clearly the efficient cut.  In 2P's sparse field there
  are few G < 13 stars, so `a` reaches only 21 % of the contamination.
* **The right additional cut is target-dependent.**  For 24P, adding `b` removes
  three quarters of the sample at lift 1.3 — most of those extra removals are
  not contaminated.  For 2P the same addition has lift 2.0, genuinely
  informative, because a faint comet is contaminated by stars far below the
  G < 13 threshold that only a flux-ratio test can see.  **Cut on `a` for bright
  targets; consider `a+b` for faint ones**, and check the lift rather than
  assuming.
* **`badphot` is anti-correlated with contamination** in both (lift 0.80 and
  0.18), which is the right answer: it is a data-quality flag about detector
  pixels, not about neighbours, and must not be used as a contamination cut.
* The ranking is stable between 2 σ and 10 σ thresholds (notebook §9).

---

## 6. Validation performed---

## 5. Validation performed

`notebooks/main.ipynb` runs against `data-sample/` and checks, with output:

| check | result |
|---|---|
| growth curve monotonic within errors | most significant decrease −0.26 σ |
| area accounting `A_geom − A_eff = badpix_area` | max residual 1.1e-13 pix² |
| `badphot` ⇔ bad pixel inside the aperture | exact |
| annulus-only bad pixels not flagged | 6 228 rows, correctly unflagged |
| negative fluxes retained, `abmag` undefined there | 524 rows, errors still finite |
| fully-masked apertures all flagged | 73 rows |
| header WCS vs. index-reconstructed WCS | max 5.2e-07 pixels over 25 cutouts |
| default vs. double-counted error | 24 % smaller |
| growth curve separates a star from a coma | −2.3 σ vs +107 σ (synthetic) |
| median vs mean stack, 2/10 frames contaminated | −0.3 vs +99.8 (truth 0) |
| distance correction | max relative error 1.9e-07 (CSV `%.8g` round-trip) |
| unit tests | 29/29 pass |

Both sample targets complete with `status = ok`; nine band stacks written for 24P.


---

## 7. Production run — the 68-comet working list

`doc/sx_comet_list_ver2607.xlsx` lists 69 designation cells; one of them is the
sheet's `Total` row, whose `desig` cell holds `68` (the *count*).
`targetlist.read_target_list` drops summary rows for exactly that reason, so the
run covers **68 comets**, all present in the index.

### Getting the data off the archive

`comets_v5_fits` is a single flat directory of ~139 000 files on an **exFAT**
volume, and exFAT scans a directory linearly — which is why it is slow to work
with.  `extract_targets.py` builds `comets_v5_fits_filtered/<target>/` holding
only the working-list exposures.  The file list comes from the Parquet index
rather than from listing the source directory, so the slow directory is never
enumerated.

* 27 797 files, 3.60 GB, **0 missing, 0 errors**, 23.8 min at 8 threads.
* macOS wrote an AppleDouble `._name` sidecar beside every copied file (exFAT has
  no native extended attributes).  At a 128 kB exFAT cluster each, 27 797 4 kB
  sidecars cost **3.4 GB of slack** and doubled every directory's entry count.
  `--clean-appledouble` removes them, and `COPYFILE_DISABLE=1` is now set in the
  script so they are not created again.  `du` fell from 6.8 GB to 3.4 GB.
* `directory.FITS_FILTERED_DIR` is first in `FITS_ROOTS`, so the pipeline prefers
  the split tree automatically with no argument.

### Photometry

`main.py --target-list ... --workers 4 --stack`

| | |
|---|---|
| targets | 68 / 68 `ok`, **0 failed** |
| exposures | 27 797 (0 missing, 0 read errors) |
| photometry rows | 447 942 |
| epochs | 141 |
| band stacks | 595 FITS |
| apertures skipped by the validity window | 135 795 |
| wall clock | 29.2 min (93.4 min of worker time across 4 processes) |
| slowest target | 2025 L1, 847 s |
| config hash | `5502194856bc` on every row — one parameter set throughout |

Sanity of the combined table: λ 0.74–5.01 µm, r_hel 1.06–15.19 au (2014 UN271 at
~15 au), r_obs 0.30–15.28 au, AB 10.7–21.8 (1st–99th pct), 49.5 % of rows at
SNR > 5, median `sky_excess_ratio` 1.19, 11.4 % negative fluxes retained, median
18 apertures per exposure.

Flag census over all 447 942 rows: `0` 188 042, `b` 157 689, `d` 74 210,
`a` 19 989, `c` 8 012 — all five populated, which the pre-2.1 magnitude-ratio
form of `b` could not achieve.

### Configuration changes carried into this run

Made outside the pipeline and picked up automatically:

* `gaia_gmag_limit` 18.0 → **18.5**
* `stack_bands` → `Continuum` 1.3–2.0, `H2O` 2.55–2.8, **`H2O_ICE` 2.8–3.3**,
  `CO2` 4.15–4.35, `CO` 4.6–4.8 µm

The spectrum plots previously carried a hard-coded band list that would have
silently disagreed with the new definitions; `plotting.bands_from_config()` now
derives the shading from `Config.stack_bands`, so the two cannot drift.

### Fixed in passing

* `status.py` built its table with `pd.concat` of a single all-None row, which
  raised pandas' empty/all-NA dtype `FutureWarning` once per target.  Rebuilt
  from records; a test now fails on *any* warning from `StatusLog`.
* `Path.glob("*.fits")` matches dot-files, unlike shell globbing — which is how
  the AppleDouble sidecars were first noticed as a doubled file count.


### Figures

`make_figures.py --target-list ... --workers 6`

**56 351 figures, 7.8 GB, 211 min**, no failures: two per exposure (55 594 — the
three-panel cutout with Gaia sources marked, and the growth curve) plus 757
per-target summaries.  Every exposure of every target has both figures.

`--skip-existing` now checks per exposure, *before* the FITS is read, so an
interrupted render resumes in seconds rather than redoing hours of work.

### Flag effectiveness across the whole survey

A single comet answers "did the flags work here"; 67 targets answer whether the
definitions generalise.  `diagnostics.survey_flag_effectiveness` scores every
target against the growth-curve contamination test
(`results/flag_effectiveness_survey.csv`, `fig/flag_effectiveness_survey.png`):

| | median recall | median lift | median data kept |
|---|---|---|---|
| `cut a` | 0.35 | **5.60** | **0.95** |
| `cut a+b` | 1.00 | 1.72 | 0.44 |
| `cut a+b+c` | 1.00 | 1.72 | 0.43 |

Contaminated fraction per target: median 7.9 %, range 0–29.8 %.

**`a` alone beats `a+b` on lift in 56 of 67 targets.**  Flag `a` is five and a
half times more likely than chance to be removing genuinely contaminated data
while discarding only 5 % of the sample.  `a+b` does catch everything, but at a
lift of 1.7 and the cost of *more than half the survey* — for most targets those
extra removals are not contaminated.

The recommendation is therefore **cut on `a`; keep `b` and `c` as columns** and
apply them per target where the science demands completeness over sample size.
The earlier two-target result (24P favouring `a`, 2P favouring `a+b`) was
directionally right but understated how strongly `a` wins in general: 2P is one
of the minority of sparse-field targets where the extra recall of `b` is worth
its cost.
