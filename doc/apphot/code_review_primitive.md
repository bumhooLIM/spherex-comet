# Code Review — Primitive SPHEREx Comet Aperture Photometry (`spherex-apphot/`)

> **Provenance note (2026-09-15).** This review describes the single-file prototype as it
> stood on 2026-09-05; every path below refers to that layout.  The prototype is kept, for
> these line references only, in `_archive/legacy/apphot/` (`apphot.py`, `apphot_all.py`, the
> three prototype notebooks with outputs stripped).  Its replacement is the `spherex_apphot`
> package; the disposition of every item below is `pipeline_upgrade_notes.md` beside this file.

**Reviewed:** 2026-09-05 · **Status:** read-only review, nothing executed or modified.

**Files reviewed**

| File | Lines | Role |
|---|---|---|
| `spherex-apphot/apphot.py` | 975 | All science logic + the per-target driver (`__main__`) |
| `spherex-apphot/apphot_all.py` | 69 | Batch orchestrator: one subprocess per target |
| `spherex-apphot/directory.py` | 11 | Path constants |

**Supporting context consulted (not reviewed in depth):** `notebooks/h5cut2fits.py` (the upstream
HDF5 → FITS exporter that defines the FITS layout and header keywords),
`notebooks/apphot_sample_2P.ipynb` and `apphot_sample_UN271.ipynb` (the prototypes `apphot.py`
was lifted from, including their stored outputs), `notebooks/reflectance.ipynb` (the downstream
consumer of the photometry CSVs), `notebooks/rcparams.py`, `notebooks/clean_double.py`.

---

## 1. Data model

### 1.1 Inputs

| Input | Path (via `directory.py`) | Notes |
|---|---|---|
| Cutout index DB | `DB_DIR/db_filtered.parq` | ~250-column denormalized Parquet: ephemeris, frame WCS, SIP + PV coefficients, L1/L2 DQA counters, and `filename` pointing at the exported FITS |
| Cutout FITS | `FITS_DIR/<filename>` | One file per (target, obsid, detector) |
| Gaia DR3 | `REFCAT_DIR/gaiadr3_all.npy` | ~11 GB structured `.npy`, opened `mmap_mode='r'` |

### 1.2 FITS layout (verified from a stored `hdul.info()` in `apphot_sample_2P.ipynb`)

```
0  PRIMARY    PrimaryHDU    7 cards
1  IMAGE      ImageHDU    202 cards   (91, 91) float32   <- full header lives here only
2  VARIANCE   ImageHDU      8 cards   (91, 91) float32
3  FLAG       ImageHDU      8 cards   (91, 91) int32
```

Cutouts are 91×91 px; SPHEREx pixels are ~6.2″, so the field is ~9.4′ (half-width 45.5 px).
Only the `IMAGE` extension carries the WCS/ephemeris header (`h5cut2fits._export_one` applies the
full header to the first extension only) — which is why `gaia_to_mask` is always handed `hdul[1]`.

### 1.3 Key header/DB quantities used

`xcen`, `ycen` (ephemeris position of the target **in cutout pixel coordinates**, quantized to
1/32 px, values ≈ 44.5–45.5 on a 91-px cutout → consistent with 0-indexed, i.e. correct for
photutils), `ltv1/ltv2` (cutout origin offset), `pix_scale` / `PIX-SCL` (″/px), `PSF_FWHM`,
`r_hel`, `r_obs`, `wl`, `wlwidth`, `sun_jy`, `vmag`, `DATE-OBS`, `jd_utc`.

---

## 2. Processing concept

The scientific goal is a **1–5 μm reflectance/emission spectrum of a comet's coma measured inside
a fixed *physical* aperture** (a fixed cometocentric radius ρ in km, not a fixed angular radius),
so that spectra taken at different geocentric distances sample the same volume of coma. SPHEREx's
spectral dimension is obtained by *scanning*: each exposure of a given sky position records one
narrow wavelength per detector pixel, so a comet's spectrum is assembled from **many separate
exposures at many epochs**, each contributing one (λ, flux) point. The pipeline therefore:

1. groups exposures into observing **epochs** (called `phase` in the code),
2. measures one flux per exposure per aperture size,
3. and assembles the per-epoch scatter plot of flux vs λ as "the spectrum".

---

## 3. Procedure, step by step (`apphot.py::__main__`)

**[1/7] Configuration.** `--objdesig` from argv. Nine physical aperture radii are hard-coded:
`[10, 15, 20, 25, 30, 40, 60, 80, 100] × 10³ km`. `gmag_limit = 18.0`, Gaia search box
`±0.2°`.

**[2/7] Load + epoch grouping.** `pd.read_parquet(..., filters=[("objdesig","==",objdesig)])` uses
predicate pushdown so only the target's rows enter RAM. `spherex_phase_grouping()` sorts by
`DATE-OBS`, takes consecutive differences, and starts a new `phase` wherever the gap exceeds
**28 days** (`is_new_phase.cumsum() + 1`).

**[3/7] Gaia subsetting.** The 11 GB memmap is boolean-masked on `phot_g_mean_mag < 18` to produce
an in-RAM "light" array. Then, per epoch, `create_gaia_subset()` narrows further:
a global Dec bounding box, then a per-pointing box loop where the RA half-width is inflated by
`1/cos δ` (clipped to 180° above |δ| > 89°) with modular-arithmetic RA wrap handling. Pointings are
deduplicated by snapping to a `0.5 × min(Δα, Δδ)` grid, with the box padded by half a grid cell to
compensate. A coverage figure is written.

**[4/7] Cross-match.** All epoch subsets are concatenated and `SkyCoord.match_to_catalog_sky` finds
the nearest Gaia star to each target position; `neargaia_gmag` and `neargaia_dist_pixel`
(= `sep2d.arcsec / pix_scale`) are attached to every row.

**[5/7] Aperture parameterization.** The plate scale is converted to a physical scale,

```
pixel_scale_km = pix_scale["] × (1" in rad) × r_obs[au] × (1 au in km)
r_ap_NN_pix    = r_ap_NN_km / pixel_scale_km
```

and the shared sky annulus is pinned to the **largest** aperture:
`r_in = 1.5 × r_ap_max_pix`, `r_out = 3.0 × r_ap_max_pix`.

**[6/7] Photometry loop.** For every row:

- `sci  = hdul[1].data` (IMAGE), `err = sqrt(hdul[2].data)` (√VARIANCE), `flag = hdul[3].data`.
- **Masking** — three sources OR-ed into a `master_mask`:
  - `flag_to_mask(flag, [2,6,7,9,10,11,12,14,15,17,22,24,26,27,28,29])` — bitwise DQ rejection;
  - `gaia_to_mask(...)` — circular stamps on Gaia stars, radius
    `radius_scale × (PSF_FWHM/PIX-SCL) × (gmag_limit − G)`, i.e. linearly larger for brighter
    stars, drawn with `skimage.draw.disk`; `radius_scale = 0.6` in the photometry path;
  - `np.isnan(sci)`.
- **`perform_aperture_photometry`** — for all 9 radii at once:
  1. one `CircularAnnulus` + `ApertureStats` with `SigmaClip(3σ, 5 iters)` → `msky` (median),
     `ssky` (std), `nsky` (area);
  2. `aperture_photometry(sci, [apertures], error=err, mask=master_mask)` → raw sums;
  3. per aperture, exact fractional unmasked area `ap_area` and `nbadpix`
     (= `ApertureStats(master_mask.astype(float), ap).sum`);
  4. `source_sum = ap_sum − ap_area × msky`;
  5. DAOPHOT-style error
     `σ = sqrt(σ_map² + ap_area·ssky² + ap_area²·ssky²/nsky)`, `snr = source_sum/σ`;
  6. `abmag = −2.5 log10(source_sum × 1e−3) + 8.90`, `abmag_err = 1.0857 σ/source_sum`
     (only when `source_sum > 0`);
  7. `badphot = (ap_area == 0) or (source_sum <= 0) or (nsky < 10)`.
  The whole body is wrapped in `try/except Exception` returning an all-NaN frame.
- 9 rows per exposure are emitted, metadata columns broadcast in, appended to a list.

All results are concatenated once and written to `OUT_DIR/<objdesig>.csv`.

**[7/7] Per-epoch products.** For each epoch: pick the *smallest* aperture whose minimum pixel
radius exceeds 2.0 px; filter to `badphot == False & snr > 1`; then
`plot_spec` (flux vs λ, colour-coded by `r_hel`, with H₂O/CO₂/CO bands shaded),
`plot_cutout` (a 5-column grid of masked postage stamps with aperture/annulus overlays), and
`combine_fits` (NaN-median stacks in five hand-defined bands — `dust_cont_1` 1.3–2.0,
`h2o` 2.6–2.8, `dust_cont_2` 3.0–4.0, `co2` 4.2–4.4, `co` 4.6–4.8 μm — written as single-HDU FITS
plus a summary figure).

**Batch layer.** `apphot_all.py` reads only the `objdesig` column, then for each unique target
skips if `OUT_DIR/<target>.csv` exists, else `subprocess.run([sys.executable, "apphot.py",
"--objdesig", target], check=True)`, logging failures to `failed_targets.log`.

---

## 4. Concerns

### 4.1 Blockers — the script cannot run as committed

| # | Issue | Location |
|---|---|---|
| B1 | **`directory.APPHOT_DIR` and `directory.COMBFITS_DIR` do not exist.** `directory.py` defines only `WORK_DIR, T7_DIR, DB_DIR, FITS_DIR, REFCAT_DIR, FIG_DIR, RESULT_DIR, DOC_DIR`. Both copies (`spherex-apphot/` and `notebooks/`) are missing them → `AttributeError` on the first line of `__main__`. `apphot_all.py` hits the same. | `apphot.py:737,738`, `apphot_all.py:14` |
| B2 | **`import rcparams` fails.** `rcparams.py` lives in `notebooks/`, not `spherex-apphot/`. The script only imports cleanly if launched with `notebooks/` on `sys.path`. | `apphot.py:28` |
| B3 | **No output directories are created.** There is no `FIG_DIR.mkdir()` / `OUT_DIR.mkdir()`; `fig/` does not exist in the repo and `results/` is empty. Only `COMBFITS_DIR` is created (inside `combine_fits`). First `savefig` → `FileNotFoundError`. | `apphot.py` §3, §6 |
| B4 | **`apphot_all.py` assumes `cwd == spherex-apphot/`** (relative `"apphot.py"`), yet the existing `failed_targets.log` sits in `notebooks/`, implying it was actually run from there. Use `Path(__file__).parent / "apphot.py"`. | `apphot_all.py:44` |

### 4.2 Scientific correctness — these change the numbers

**S1 — Units are wrong: the images are MJy/sr, not mJy/pixel.**
This is the most consequential finding. `perform_aperture_photometry`'s docstring asserts
"Values should be in milliJanskys (mJy) per pixel", every output column is suffixed `_mjy`, and
`abmag` is computed as if `source_sum` were already in mJy. But `h5cut2fits.EPH_KEYWORDS` documents
the upstream sums as `[MJy/sr*pixscale^2]`, and the earlier prototype cells in
`apphot_sample_2P.ipynb` do the conversion explicitly:

```python
apsum_mjy = apsum * 1e9 * (fits_summary["PIX-SCL"]/206265)**2
```

That conversion was **dropped** when the code was refactored into `apphot.py`. The missing factor is

```
1e9 × (pix_scale / 206264.8)²  =  0.9035  at 6.2″/px      (0.8462 at 6.0″/px)
```

so every `source_sum_mjy`, `aperture_sum_mjy`, `source_sum_err_mjy`, `annulus_median_mjy_per_pix`
and `bkg_std_mjy_per_pix` is **too high by ×1.107**, and every `abmag` is **0.110 mag too bright**
(0.18 mag at 6.0″). It is close enough to unity to have passed unnoticed. Because `pix_scale` varies
per row, this is not even a constant offset. `combine_fits` propagates the error by stamping
`BUNIT = 'mJy/pixel'` on the stacked FITS. `reflectance.ipynb` also divides `source_sum_mjy` (mJy)
by `sun_jy` (Jy) — harmless after normalization, but dimensionally inconsistent and worth fixing
at the same time.

**S2 — The sky annulus leaves the cutout for any nearby target, and samples the coma for all of them.**
`r_in/r_out` are pinned to the *largest* aperture (100 000 km), giving annulus radii of 150 000 and
300 000 km. In pixels, with a 45.5 px cutout half-width:

| `r_obs` | km/px | `r_ap`(10⁴ km) | `r_ap`(10⁵ km) | `r_in` | `r_out` |
|---|---|---|---|---|---|
| 0.5 au | 2 248 | 4.45 px | 44.5 px | **66.7 px** | **133.4 px** |
| 1.0 au | 4 497 | 2.22 px | 22.2 px | **33.4 px** | **66.7 px** |
| 1.5 au | 6 745 | 1.48 px | 14.8 px | 22.2 px | 44.5 px |
| 4.07 au (2P) | 18 301 | 0.55 px | 5.46 px | 8.2 px | 16.4 px |

`r_out ≤ 45.5 px` requires **`r_obs ≳ 1.47 au`**; `r_in ≤ 45.5 px` requires `r_obs ≳ 0.73 au`. Below
~1.5 au the annulus is progressively truncated by the cutout edge — the surviving sky pixels are an
asymmetric corner sample, and `nsky` silently shrinks; below ~0.73 au the annulus is entirely off-image
and `msky` becomes NaN or is drawn from nothing. There is **no check that the annulus fits inside the
array.** Separately, and independent of geometry: a coma commonly extends to 10⁵–10⁶ km, so an annulus
at 150 000–300 000 km is measuring *comet*, not sky. Subtracting it removes exactly the extended
signal the project is trying to measure — a bias that grows with activity level.

**S3 — No aperture (encircled-energy) correction, with sub-PSF apertures.**
For 2P at 4.07 au the 10 000 km aperture is **r = 0.55 px** (area 0.94 px²), i.e. comparable to or
smaller than one pixel and smaller than the PSF half-width (`PSF_FWHM` ≈ 6–7″ ≈ 1.0–1.1 px). Such an
aperture captures only a fraction of a point-like inner coma, no correction is applied, and — because
the SPHEREx PSF broadens toward long wavelength — **the fraction captured is itself a function of λ**.
The resulting flux-vs-λ curve therefore contains a wavelength-dependent aperture-loss slope that is
indistinguishable from a real reflectance slope. Since `r_ap_pix = r_ap_km / (pix_scale × r_obs)`,
the loss also varies with `r_obs`, i.e. across epochs. This is a first-order threat to the science
product and should be the highest-priority fix after S1.

**S4 — No centroiding.** `xcen, ycen` come straight from the ephemeris. Comet astrometry
(non-gravitational forces, stale orbit solutions) is routinely off by ≳1″; at 6.2″/px that is a
sizeable fraction of a pixel, and with r ≈ 0.5–2 px apertures a fractional-pixel centering error
translates directly into flux error. There is no re-centroiding, no fallback, and no diagnostic on
the offset. Related: `sky_motion` / `sky_motion_pa` are carried in the DB but never used, so
trailing during the exposure is unmodelled and circular apertures are assumed.

**S5 — Positive-flux truncation biases faint spectra high.**
`badphot` sets `source_sum <= 0` as *bad*, `abmag` is only computed for `source_sum > 0`, and the
plotting/stacking stage further imposes `snr > 1`. For low-S/N points the noise distribution is
symmetric about the true flux; discarding the negative half and then averaging/median-stacking the
survivors produces a systematically **positive** bias. For faint volatile bands — the science target —
this can manufacture a detection. Negative fluxes are physically meaningful measurements and must be
retained through any averaging; only mark them, do not drop them.

**S6 — Sky noise is double-counted in the error budget.**
`σ² = σ_map² + ap_area·ssky² + ap_area²·ssky²/nsky`. The DAOPHOT formula's `ap_area·ssky²` term is
the per-pixel sky *noise*, appropriate when there is no error map. Here `σ_map²` comes from the
VARIANCE extension, which already includes the background's photon noise. The second term therefore
adds it a second time, inflating `source_sum_err_mjy` and deflating `snr` (and the `snr > 1` cut then
throws away real data). Either drop the term or use `ssky` only for the extra systematic
(confusion/flat) component.

**S7 — Flux lost to masked pixels inside the aperture is never corrected.**
`aperture_photometry(mask=...)` zeroes masked pixels, and `ap_area` is the *unmasked* fractional
area, so the sky subtraction is self-consistent — but the source flux under a masked pixel is simply
gone. `nbadpix` is recorded and then never used: there is no threshold on `nbadpix/ap_area`, no
rescaling, and `badphot` does not consider it. An exposure with a cosmic ray or a masked star on the
nucleus silently yields a low flux.

**S8 — Median stacking is done on un-background-subtracted, un-registered, un-rotated frames.**
`combine_fits` stacks `sci_masked` directly. Three problems: (i) no background subtraction, and the
zodiacal foreground varies strongly between visits, so the median is dominated by the sky level
rather than the comet; (ii) alignment is `int(round(xcen))` only — up to 0.5 px of registration
jitter, which is large next to a ~1 px PSF; (iii) SPHEREx roll angle differs between visits and no
rotation to a common frame is applied, so any real coma morphology (tail direction, sunward jet) is
smeared. Additionally the output FITS carries **no WCS**, no `NCOMBINE` weighting, and a band with a
single contributing frame produces a "median" of one with no minimum-frame guard.

**S9 — `neargaia_dist_pixel` is computed and never used.** The whole cross-match step (§4) produces
a star-contamination metric that no filter, flag, or output decision consumes. Either use it (e.g.
flag rows where a bright star falls within `r_out`) or drop the step.

**S10 — `radius_scale` is inconsistent between paths.** The photometry loop masks Gaia stars with
`radius_scale=0.6`; `plot_cutout` and `combine_fits` both use the function default `1.0`. The stamps
shown to the user and the frames that get stacked are therefore masked *differently* from the frames
that were actually measured. Also, the radius law `∝ (18 − G)` is linear and unbounded below G ≈ 18
and has no floor/ceiling; and `psf_fwhm_pixel = PSF_FWHM / PIX-SCL` silently assumes `PSF_FWHM` is
stored in arcsec — worth confirming against the header, since a pixel-valued `PSF_FWHM` would inflate
every mask radius by ~6×.

**S11 — Reflectance conversion is incomplete (downstream).** `reflectance.ipynb` computes
`refl = source_sum_mjy / sun_jy` with no `r_hel²` scaling and no `r_obs²`/phase-angle term. Within one
epoch these are near-constant and normalize away, but comparing epochs or targets requires them.

### 4.3 Robustness & correctness of the plumbing

| # | Issue |
|---|---|
| R1 | **`apphot_all.py`'s resume check is broken for multi-word designations.** `apphot.py` strips spaces before naming the output (`"2021 G2"` → `2021G2.csv`), but the orchestrator tests `OUT_DIR / f"{target}.csv"` with the *unstripped* name. Every target whose designation contains a space is reprocessed on every run — which is most comets. (The one entry in `failed_targets.log` is `2020 R7`.) |
| R2 | **HDUs are addressed positionally** (`hdul[1]/[2]/[3]`) rather than by name. It happens to be correct for this dataset, but `h5cut2fits` writes extensions in the order `IMAGE, MASK, FLAGS, VARIANCE, FLAG` filtered to whatever exists in the HDF5 group — so a group carrying `MASK` or naming its DQ plane `FLAGS` instead of `FLAG` silently swaps VARIANCE and FLAG. That failure mode is *silent and catastrophic* (√FLAG as the error map). Use `hdul['IMAGE']`, `hdul['VARIANCE']`, `hdul['FLAG']`. |
| R3 | **`except Exception: return all-NaN`** in `perform_aperture_photometry` swallows every bug with no logging. Across ~10⁵ exposures this converts programming errors into quiet data loss that is indistinguishable from genuine bad pixels. Log the exception with the filename at minimum. |
| R4 | **The main photometry loop does not guard `fits.open`.** `plot_cutout` and `combine_fits` both check `fpath.exists()`; the loop that actually matters does not. One missing file raises, `apphot_all` catches the non-zero exit, and **the entire target is discarded** after all its work — the failure mode that likely produced `failed_targets.log`. |
| R5 | **NaNs in VARIANCE are not masked.** `master_mask` includes `np.isnan(sci_data)` but not `np.isnan(err_data)`; a negative or NaN variance yields `sqrt` → NaN → a NaN `aperture_sum_err` and a NaN `snr` that is not flagged by `badphot`. |
| R6 | **`get_padded_cutout` will raise on a NaN centroid** (`int(round(nan))` → `ValueError`), inside `combine_fits` which has no try/except. |
| R7 | **The no-error-map fallback is dimensionally meaningless**: `ap_sum_err_sq = np.abs(ap_sum - ap_area*msky)`, i.e. the absolute source flux used as a *variance*. That is Poisson reasoning applied to a surface-brightness unit. It is unreachable in the current main path but is a trap for anyone reusing the function. |
| R8 | **Stale columns from the DB.** `db_filtered.parq` already carries `phase` and `r_ap_01..06_{km,pix}` (visible in the prototype notebook's dataframe). The script overwrites `phase` and `r_ap_01..06`, then *adds* `07..09`, so the working frame mixes fresh and stale values under identical names. Worth verifying against the actual Parquet schema and either dropping or renaming on load. |
| R9 | **The magic flag-bit list `[2,6,7,9,10,11,12,14,15,17,22,24,26,27,28,29]` is duplicated verbatim in three places** with no comment on what any bit means. Any future change will be applied inconsistently. Promote to a documented module constant. |
| R10 | **`plot_cutout`'s `imshow(extent=[x_min,x_max,y_min,y_max])` is off by half a pixel** relative to the `Circle((xcen,ycen), ...)` overlays, because `extent` refers to pixel *edges* while the aperture centre is in pixel-*centre* coordinates. Cosmetic, but this is the figure a human uses to judge centering — a systematic half-pixel misregistration in the diagnostic is a bad property. |
| R11 | **Dead code / contradictions.** `df_cutout` is built with the identical mask as `df_spec` while its comment says "Cutouts don't strictly need the SNR filter"; `chosen_ap_km` is printed as "Selected …" but `combine_fits` is handed `median_ap_pixel` instead; `hasattr(phot_summary_filtered,'iloc')` is always true so the `'all'` fallback is unreachable; `math`, `RESULT_DIR`, and several `typing` imports are unused; `numpy`, `pandas`, and `astropy.io.fits` are each imported twice. |
| R12 | **`np.concatenate` over the per-epoch Gaia subsets duplicates stars** that fall in overlapping epochs. Harmless for a nearest-neighbour query, but it inflates the array and would corrupt any count-based use. |

### 4.4 Performance & scale

| # | Issue |
|---|---|
| P1 | **The 11 GB Gaia catalog is fully read from disk once per target.** `gaia_all_mmap[valid_mag_mask]` is a boolean fancy-index over a memmap — it touches every byte of the file. Because `apphot_all.py` spawns a *fresh subprocess per target*, this is `11 GB × N_targets` of I/O, and the in-code comment calls it a "speed hack". The fix is trivial and large: pre-filter `G < 18` **once** and cache the result (a few hundred MB) to disk; every run then loads the small file. |
| P2 | **`create_gaia_subset` is O(N_pointings × N_gaia)** — a full boolean pass over the catalog per unique pointing, in a Python loop. `SkyCoord.search_around_sky` or a `cKDTree` on unit vectors would collapse this. |
| P3 | **`gaia_to_mask` rebuilds `WCS(hdr)` and re-projects the whole star list for every exposure**, and is called again (redundantly) inside `plot_cutout` and `combine_fits`, so the same file is re-opened and re-masked up to three times per epoch. |
| P4 | **Figure sizes are unbounded.** `plot_coverage` uses `figsize=(8 × num_phases, 8)`; `plot_cutout` uses `nrows = ceil(N/5)` with 4-inch rows, so an epoch with 200 exposures produces a 20 × 160 inch canvas. At the `savefig.dpi = 200` set in `rcparams.py` that is a ~19 000 px image. `CLAUDE.md` specifies dropping to **dpi 50** for batch output — the pipeline never does. |
| P5 | **`apphot_all.py` is strictly sequential** with full interpreter + import + Gaia-load startup cost per target. Combined with P1 this dominates wall time. Per-target work is embarrassingly parallel. |

### 4.5 Conventions & project hygiene

- **Undeclared dependencies.** `CLAUDE.md` lists numpy/pandas/astropy/photutils as core and forbids new
  third-party deps without confirmation. `apphot.py` imports **`skimage`** (scikit-image, for
  `draw.disk`); `reflectance.ipynb` additionally uses **`sep`**, and `h5cut2fits.py` uses **`click`**,
  **`h5py`**, and a local **`pqfilt`**. `disk` is ~5 lines of numpy and could remove the scikit-image
  dependency outright.
- **Naming collision on "phase".** In this code `phase` means *observing epoch*, while the DB's
  `alpha` is the *solar phase angle* — the standard meaning of "phase" in comet photometry. For an
  audience of professional astronomers this is actively misleading. Rename to `epoch` or `visit`.
- **`rcparams.py` does not match `CLAUDE.md`.** The file sets `font.size = 15`, `axes.titlesize = 15`,
  `xtick/ytick.labelsize = 12`, `figure.titlesize = 20`; `CLAUDE.md` documents 20 / 22 / 18 / 26. One
  of the two should be corrected. `legend.fontsize` is documented but not set in the file at all.
- **Docstring drift.** `gaia_to_mask`'s docstring documents a `dec_col` parameter that does not exist
  and lists `radius_scale` twice with two different defaults (1.0 and 0.3).
- **No provenance in the outputs.** The CSV records no configuration — aperture list, flag bits,
  `gmag_limit`, `radius_scale`, annulus scaling, code version. For a science product intended for
  publication this should be written into a sidecar or as CSV header comments.
- **No tests.** There is no unit test for the pieces that are cheap to test and easy to get wrong:
  `flag_to_mask` bit semantics, `reconstruct_cutout_wcs` round-tripping against the on-disk header,
  the km↔pixel conversion, and the sky/error arithmetic on a synthetic frame with known answers.
- **WCS note (minor, likely benign).** `reconstruct_cutout_wcs` scoops SIP terms with
  `^[AB]P?(_ORDER|_\d_\d)$`, which correctly captures `A_*`, `B_*`, `AP_*`, `BP_*`. The DB also
  carries `PV1_0..10` / `PV2_0..10`, which are **not** copied. If `CTYPE` is `RA---TAN-SIP` those are
  vestigial and ignoring them is right; the assumption is just worth confirming once. Note also that
  this function is used only for the coverage plot — the photometry path reads the WCS straight from
  the `IMAGE` header, so the two are independent and could drift apart.

---

## 5. Suggested order of work

1. **Make it run**: define `APPHOT_DIR`/`COMBFITS_DIR`, `mkdir(parents=True, exist_ok=True)` all
   outputs, fix the `rcparams` import path and the `cwd` assumption (B1–B4), and fix the
   space-stripping mismatch in the resume check (R1).
2. **Fix the units** (S1) — a one-line factor, but it invalidates every number currently produced,
   so it must land before anything is regenerated. Add a `BUNIT`-driven assertion so it cannot
   silently regress.
3. **Fix the aperture/annulus geometry** (S2): decouple the annulus from `r_ap_max`, and validate
   that `r_out` fits inside the array — reject or fall back explicitly when it does not. Decide
   deliberately how to estimate a background that is not itself coma.
4. **Add centroiding + an encircled-energy correction** (S3, S4). Without these the extracted
   spectra carry a λ-dependent slope of instrumental origin.
5. **Stop truncating at zero flux** (S5) and **fix the error budget** (S6). Keep negative fluxes;
   flag rather than drop.
6. **Harden the loop** (R2–R5): address HDUs by name, log exceptions with filenames, guard missing
   files per-row so one bad file costs one row and not one target.
7. **Cache the G < 18 Gaia subset once** (P1) and parallelize `apphot_all` (P5). These are the two
   changes that make a full-catalog run practical.
8. **Then** revisit stacking (S8) — background-subtract, register sub-pixel, rotate to a common
   frame, and write a real WCS.
