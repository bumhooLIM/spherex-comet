# Primitive Code Analysis — `ztf-comet` + `ztf-sso-query`

> **Provenance note (2026-09-15).** This review describes the two pre-merge projects as they
> stood on 2026-09-07; every path below refers to that layout.  The code it reviews is
> superseded by the `ztfcomet` package (now the ZTF stage of `spherex-comet`); the pre-merge
> notebooks it cites are kept, outputs stripped, in `_archive/legacy/ztf/`.  The defects are
> pinned by `tests/test_ztfcomet.py` (`test_regression_c*`).

**Date:** 2026-09-07
**Scope:** `/Users/bumhoo7/Desktop/claude/ztf-comet` and `/Users/bumhoo7/Desktop/repo/ztf-sso-query`
**Purpose:** Baseline review before merging both trees into a single `ztfcomet` repository.

All quantitative claims below were measured on the sample dataset
`ztf-sso-query/_data/2024E1` (113 URLs, 111 readable FITS cutouts) in the `spherex`
conda environment. Nothing here is inferred from reading alone.

---

## 1. What exists today

### 1.1 `ztf-sso-query` — the query half

A real, installed Python package (`ztfssoquery` 0.1.0, editable in `spherex`), tracked
in git with remote `https://github.com/bumhooLIM/ztf-query-comet.git` (5 commits, branch `main`).

```
ztf-sso-query/
├── pyproject.toml            setuptools; deps declared = astropy, numpy only
├── README.md                 one line: "# ztfssoquery"
├── .gitignore                ignores _doc/, _data/, __pycache__, .DS_Store
├── ztfssoquery/
│   ├── __init__.py           re-exports 7 names from construct_fitsurl
│   ├── construct_fitsurl.py  285 lines — THE live module
│   ├── _construct_fitsurl.py 12 KB legacy script (mode 0700)
│   └── ztf-fitsurl.py        12 KB — byte-identical to the above
├── notebooks/ztfquery_2019Y3.ipynb    4 code cells, thin driver
├── example-data/2019Y3/      36 MB, 35 files — TRACKED IN GIT
├── _data/{2019Y3,2024E1}/    172 MB — gitignored working copy
└── _doc/ztf_pipelines_deliverables.pdf   gitignored
```

**Public API** (`ztfssoquery/__init__.py`): `generate_fits_urls`, `import_fitsurl`,
`query_sso_ephemeris`, `query_ztf_metadata`, `construct_fitsurl`, `save_fitsurl`,
`extract_lastrecnum`.

**Pipeline implemented by `generate_fits_urls`:**

1. Coarse JPL Horizons ephemeris for the target at observatory `I41`, stepped every
   `interval_days`, quantities `1,3,9,19`.
2. Cut on `r < rh_max` and `Tmag < vmag_max`; write `eph.csv`.
3. For each surviving ephemeris row, size a search box from the apparent rate
   (`RA_rate`/`DEC_rate`, arcsec/hr → deg/day × `interval_days`) and query the IRSA
   IBE ZTF science-image metadata service over `obsjd ± 0.5·interval_days`.
4. Re-query Horizons at the exact `obsjd` of every returned frame (quantities
   `1,3,4,8,9,16,18,19,20,24,27`), `pd.concat(axis=1)` onto the ZTF metadata.
5. Keep frames whose footprint corners (`ra1..ra4`, `dec1..dec4`) bracket the target;
   write `ztf.csv`.
6. Build cutout URLs (`?center=RA,DEC&size=10arcmin&gzip=false`); write `fits_urls.txt`.
7. `import_fitsurl` streams each URL to disk, skipping paths that already exist.

This works. On 2024 E1 it produced 113 URLs from 2025-01-01 → 2025-10-30.

### 1.2 `ztf-comet` — the analysis half

**Not a git repository.** No package, no `pyproject.toml`, no `README.md`.

```
ztf-comet/
├── CLAUDE.md                 3.2 KB — describes a layout that does not exist
├── handoff.md                0 bytes
├── data/                     EMPTY
├── fig/                      EMPTY
├── ztf-comet/                EMPTY nested dir
├── 무제 폴더/                 EMPTY ("Untitled Folder")
└── notebooks/
    ├── rcparams.py           matplotlib rcParams (savefig.dpi=200)
    ├── ztfquery_24P.ipynb    1.14 MB — query + cutout figures
    ├── ztfquery_2P.ipynb     1.18 MB — same, for 2P
    ├── ztfquery_240P.ipynb   2.50 MB — same, for 240P
    ├── afrho_240P.ipynb      0.14 MB — THE photometry / Afρ pipeline
    └── 240P.ipynb            0.02 MB — scratch: fragment separation of 240P
```

**`afrho_240P.ipynb` is the only Afρ implementation in either tree.** Its 21 cells are
the science pipeline, and none of it is packaged:

| Step | Implementation |
|---|---|
| Frame inventory | `ccdproc.ImageFileCollection(...).summary` → DataFrame |
| Per-frame ephemeris | Horizons at `obsjd`, `location='I41'`, by `orb_id` |
| Aperture scale | fixed `rho_km = 15000`; `rho_pix = rho_km / (pixscale·Δ)` small-angle |
| Sky annulus | `r_in = 3·rho_pix`, `r_out = 4·rho_pix + 20` |
| Centroid | `sep.winpos` seeded from `WCS.world_to_pixel` of the Horizons RA/Dec |
| Photometry | `photutils.aperture_photometry`, 3σ-clipped median sky, 5 iterations |
| Calibration | `filter_mag = inst_mag + MAGZP` |
| Solar flux ratio | Willmer (2018) PS1 AB solar mags: g −26.54, r −26.93, i −27.05 |
| Afρ | `4·Δ²·r²/ρ · (F_c/F_☉)`, `.to_value(u.cm)` |
| Phase correction | linear, `Φ = 10^(0.4·β·α)`, β = 0.03 mag/deg |

The Afρ formula itself is correct — `Δ` carries AU units, `r` enters as a bare number
in AU, `ρ` in km, and astropy resolves AU²/km → cm properly. The problems are upstream
of it.

### 1.3 How the halves connect

Only by import: `ztf-comet/notebooks/*.ipynb` do `import ztfssoquery`, resolved from
the editable install of the *other* repository. There is no declared dependency, no
version pin, and nothing in `ztf-comet` records that the coupling exists.

---

## 2. Concerns

Ordered by severity. Each is labelled with how it was established.

### 2.1 Blocking — wrong science

---

#### C1. The zeropoint of the last image is applied to every image

`afrho_240P.ipynb`, cell 8:

```python
phot_target["filter_mag"] = phot_target["inst_mag"] + row.zpmag
```

`row` is the loop variable left over from cell 7. It is a scalar — the final frame's
`MAGZP` — broadcast across the whole table. The per-row column `phot_target["zpmag"]`
was populated in cell 3 and then never used.

Measured `MAGZP` spread in the 2024 E1 sample:

| filter | mean | std | min | max |
|---|---|---|---|---|
| ZTF_g | 26.117 | 0.520 | 23.033 | 26.353 |
| ZTF_r | 26.206 | 0.103 | 25.767 | 26.368 |

The bug also silently mixes filters: whichever frame happens to be last sets the
zeropoint for g and r alike.

**Measured impact**, from running both versions over the same 111 frames
(old = last frame's `MAGZP` for all, new = per-frame):

| | median error | max error | Afρ factor |
|---|---|---|---|
| all 111 frames | 0.138 mag | 3.075 mag | **×17.0** |
| 82 frames passing quality cuts | 0.142 mag | 0.341 mag | ×1.37 |

The ×17 case is a single ZTF_g frame (`ztf_20250826162176_...`) whose own
`MAGZP` is 23.03 — and which the new quality flags reject anyway as low-SNR. So
the honest statement is not "everything is wrong by ×21": it is that **55 of the
82 usable frames are wrong by more than 10%**, with a worst usable case of 37%.
That is not a normalisation offset that divides out; it is epoch-dependent
scatter injected straight into the lightcurve, so it distorts shape — exactly
the signal an Afρ study is trying to measure.

Every `afrho_cm` and `afrho_cm_corr` value produced before the fix, and both
lightcurve figures in that notebook, are affected.

**Fix:** `phot_target["inst_mag"] + phot_target["zpmag"]`.

---

#### C2. HTTP error pages are written to disk as `.fits` and never repaired

`construct_fitsurl.py:274-281` writes the response body whenever `status_code == 200`.
IRSA returns 200 with an HTML error body for some cutout requests.

Measured in `_data/2024E1`: **2 of 113 files are 245-byte HTML documents** named
`*_sciimg.fits`:

```
ztf_20250114536794_000724_zg_c10_o_q4_sciimg.fits   245 bytes
ztf_20250114558461_000724_zr_c10_o_q4_sciimg.fits   245 bytes
```

Content: `<title>404 Not Found</title>`.

Re-fetching the first of those URLs today returns HTTP 200, `application/fits`,
1 437 120 bytes — the data is fine, the download was not. But line 271:

```python
if save_path.exists():
    continue
```

means **no re-run will ever fix it**. The corrupt files are permanent until deleted by
hand. Downstream, `ccdproc.ImageFileCollection` silently drops them, so the frames
vanish from the analysis without any message.

**Fix:** validate `Content-Type` and the FITS `SIMPLE` magic before committing the
file; write to a `.part` file and rename on success; add a `--force`/`overwrite`
path and a manifest recording per-URL status.

---

#### C3. 37% of cutouts are edge-clipped, and nothing checks it

A 10 arcmin cutout at 1.012″/px should be ~593 px square. Measured shapes in
`_data/2024E1`:

- **41 of 111 readable cutouts are truncated** by the CCD-quadrant boundary
  (e.g. 252×594, 530×303, 593×457).

Projecting the Horizons RA/Dec into each frame's WCS and comparing against
`rho_pix` and `sky_out = 4·rho_pix + 20`:

| condition | frames |
|---|---|
| target falls **outside** the cutout entirely | 6 / 111 |
| aperture runs off the array edge | 6 / 111 |
| sky annulus runs off the array edge | **12 / 111** |
| min. distance from target to nearest edge | −45.9 px |

Truncated apertures lose flux; truncated annuli bias the sky estimate. Neither is
flagged, so the affected points enter the Afρ lightcurve indistinguishable from good
ones.

**Fix:** compute the edge distance at photometry time; reject or flag frames where the
annulus is not fully contained; ideally request the cutout so the target is centred
and reject frames whose footprint cannot supply the full box.

---

### 2.2 Serious — silently degrades results

---

#### C4. The ZTF colour term is never applied

Every header carries `CLRCOEFF` (with `CLRCOUNC`, `ZPCLRCOV`, `CLRMED`). The correct
ZTF calibration is `mag = −2.5·log₁₀(DN) + MAGZP + CLRCOEFF·colour`. The notebook stops
at `MAGZP`.

Measured range in the sample: `CLRCOEFF` ∈ [−0.2028, +0.1250], median |CLRCOEFF| = 0.082.
For a comet whose g−r differs from the median calibrator colour (`CLRMED` ≈ 0.48) by a
few tenths, this is a systematic of a few percent in Afρ — and it is
filter-and-epoch-dependent, so it distorts the *shape* of the lightcurve, not just its
normalisation.

Comets are redder than the stellar calibrators; this is not a random error.

---

#### C5. `MAGZP` is the PSF zeropoint, used on aperture sums with no aperture correction

The header comment is explicit: `MAGZP = Magnitude zero point for PSF-photometry [mag]`.
The headers also supply `APCOR1..APCOR6` (`PSF_mag − AP_mag` for 2–14 px diameters).

For the large apertures used here the correction is small — `APCOR6` (14 px diam.) is
−0.0095 mag — but it is unexamined and undocumented, and it grows fast for the small
apertures noted in C6.

---

#### C6. The aperture can be smaller than the seeing disc

`rho_km` is hardwired to 15000 km for all epochs, so `rho_pix` shrinks as Δ grows.
Measured over the sample:

- `rho_pix` ∈ [3.8, 8.0] px
- `rho_pix / FWHM` ∈ [0.89, 5.31]
- **5 of 111 frames have `rho/FWHM < 1.5`** — the aperture no longer contains the PSF,
  so the "15000 km aperture" is not measuring 15000 km of coma.

Telling detail: `phot_target["rho_fwhm"]` is computed in cell 6 and then never
referenced again. The diagnostic exists; the cut was never wired up.

---

#### C7. `sqrt` of negative pixel values produces a NaN error array

Cell 7: `err = np.sqrt(data/row.egain + (hdr['READNOI']/row.egain)**2)`.

ZTF science images are not background-subtracted, so this is usually fine — but not
always. Measured: **1 of 111 frames has 97.5% NaN in its error array** (and emits
`RuntimeWarning: invalid value encountered in sqrt`). The NaN propagates into
`source_sum_err`, `snr`, and `mag_err` without ever raising.

Related: `inst_mag = −2.5·log10(source_sum)` is NaN whenever a faint or mis-centred
source gives a negative background-subtracted sum. No detection flag distinguishes
"non-detection" from "processing failure".

---

#### C8. No uncertainty is propagated to Afρ

`source_sum_err`, `mag_err` and `zpmagerr` are all computed and then dropped. There is
no `afrho_err` column. `MAGZPRMS` (0.047 mag in the frame inspected) is ignored, as is
the β = 0.03 mag/deg phase-coefficient uncertainty. The lightcurve plots have no error
bars — not because the errors are small, but because they were never carried through.

---

#### C9. Wrong target ID pasted into the 2P notebook

`ztfquery_2P.ipynb` cell 9, in a notebook whose `targetname = "2P"`:

```python
target_id = 90000355  # 24P/Schaumasse
```

Copy-pasted from `ztfquery_24P.ipynb`. Every ephemeris annotation on the 2P cutout
figures describes 24P. This is the failure mode that hardcoded per-target constants
scattered through notebooks invites.

---

### 2.3 Fragile — works now, will break later

---

#### C10. `pd.concat(axis=1)` joins ephemeris to metadata by position

`construct_fitsurl.py:206-211` concatenates the Horizons block onto the ZTF block by
row position. That is correct **only** because Horizons happens to return epochs in
ascending order and the code sorts `ztf` by `obsjd` first.

I verified the behaviour directly: sending `[dup, dup, jd₀, jd₁, jd₂]` returns five
rows in ascending JD order (Horizons sorts, and does not de-duplicate). The current
code is therefore right today — by coincidence of an undocumented service behaviour,
not by construction. One upstream change to Horizons row ordering and every ephemeris
in the table silently attaches to the wrong image.

**Fix:** merge explicitly on JD with a tolerance.

---

#### C11. Duplicate `airmass` column after the concat

Both IRSA metadata and Horizons supply `airmass`. After `pd.concat(axis=1)` the
DataFrame carries the label twice — confirmed in the raw `ztf.csv` header (pandas
mangles it to `airmass.1` on re-read, hiding the problem). In memory, `df["airmass"]`
returns a two-column DataFrame, so any downstream scalar use raises or misbehaves.

The same collision pattern applies to `ra`/`dec` (ZTF) vs `RA`/`DEC` (Horizons) —
distinguished only by case, which is a bug waiting for a `str.lower()`.

---

#### C12. `query_sso_ephemeris` can raise `UnboundLocalError`

```python
try:
    eph = obj.ephemerides(...)
except ValueError as e:
    id_record = extract_lastrecnum(str(e))
    if id_record:
        try:
            eph = obj.ephemerides(...)
        except Exception as e:
            print(f"Failed to get ephemerides ...")   # eph never assigned
    else:
        return pd.DataFrame()
return eph.to_pandas()                                # UnboundLocalError
```

If the record-number retry fails, the function prints a message and then crashes on
the return line with a confusing traceback that hides the real cause.

---

#### C13. Bare `except Exception` hides coverage gaps

`generate_fits_urls` wraps the whole per-epoch body in `except Exception as e:
print(f"Skipping ephemeris row; {e}")`. A transient network failure and a genuine
"no ZTF coverage" are indistinguishable, and neither is recorded. The final `ztf.csv`
gives no way to tell whether an epoch was absent or merely dropped.

Likewise `query_ztf_metadata` swallows every parse error and returns `None`, which the
caller cannot tell apart from an empty result.

---

#### C14. No retry, backoff, or rate limiting

Every Horizons and IRSA call is a bare `requests.get` / `astroquery` call. A long
multi-target run will hit a transient failure and lose that epoch (via C13). Note that
`import time` at `construct_fitsurl.py:6` is unused — a leftover from an earlier
throttling attempt.

---

### 2.4 Structural — the reason for the merge

---

#### C15. `CLAUDE.md` documents an architecture that does not exist

| `CLAUDE.md` claims | Reality |
|---|---|
| `doc/` — technical guidebooks | did not exist until this document |
| `data-sample/` | does not exist |
| `results/` | does not exist |
| `ztf-api-call/` | does not exist |
| `afrho/` | does not exist |
| `scripts/directory.py` | `scripts/` does not exist |
| `notebooks/directory.py` | **does not exist** |
| `import rcparams` | the notebooks all `import _rcparams` |
| `fig/`, `data/` | exist but are **empty** |

`directory.py` is named as the single source of truth for every path in the project.
It has never been written. Instead each notebook re-derives `WORKDIR = Path.cwd() / ".."`,
which silently breaks the moment a notebook is executed from anywhere but `notebooks/`.

The rcparams mismatch is worse than cosmetic: `import _rcparams` fails on a clean
checkout, so **every plotting notebook is broken as committed**.

---

#### C16. `handoff.md` is empty

0 bytes. The state-management protocol in `CLAUDE.md` — read it at the start of every
task, update at natural breakpoints — has never been exercised.

---

#### C17. 36 MB of FITS binaries are tracked in git

`git ls-files example-data` → 35 files, 36 MB. `.git` is already 28 MB on 5 commits.
`_data/` (172 MB) is correctly gitignored, but `example-data/2019Y3` duplicates
`_data/2019Y3` almost exactly and is committed. Binaries never delta-compress; this
grows monotonically and will follow the history into the new `ztfcomet` repo unless it
is dropped at the merge.

---

#### C18. Dead and duplicated source shipped inside the package

`ztfssoquery/_construct_fitsurl.py` and `ztfssoquery/ztf-fitsurl.py` are
**byte-identical** (verified by `diff`) 12 KB legacy scripts with module-level
side effects (`OUTDIR.mkdir()` at import), mode `0700`, and an old CamelCase API
(`Query_Ephemerides`, `Extract_LastRecNum`). Neither is imported by `__init__.py`, both
are listed in `SOURCES.txt`, and `ztf-fitsurl.py` cannot even be imported as a module
because of the hyphen.

---

#### C19. Declared dependencies do not match the imports

`pyproject.toml` declares `astropy` and `numpy`. `construct_fitsurl.py` imports
`pandas`, `requests`, `tqdm`, `astroquery` — none declared. `requires-python = ">=3.9"`
contradicts `CLAUDE.md`'s ">3.10".

The analysis notebooks additionally import `photutils`, `sep`, `ccdproc`, `matplotlib`,
plus `sxobsplan` and `astrometry` — the latter two **imported but never used** in both
`afrho_240P.ipynb` and `240P.ipynb`. `sxobsplan` is a local package from an unrelated
project; carrying it forward would make the new repo depend on the user's private
environment.

---

#### C20. The cutout-figure block is copy-pasted four times

The ~70-line annotated-cutout plot (N/E arrows, −V and −☉ vectors, ephemeris text box)
appears in `ztfquery_24P.ipynb` twice, `ztfquery_2P.ipynb` once, and `afrho_240P.ipynb`
once. The variants have already drifted — the `afrho_240P` copy adds an `Elong` line,
the 24P copies do not; the single-frame copies keep a blank line the loop copies drop.

---

#### C21. Notebook outputs are committed

`ztfquery_240P.ipynb` is 2.50 MB, `ztfquery_2P.ipynb` 1.18 MB, `ztfquery_24P.ipynb`
1.14 MB — almost entirely base64 PNGs. Combined with C17, the new repo starts heavy.

---

#### C22. Miscellaneous

- `무제 폴더/` ("Untitled Folder"), `ztf-comet/ztf-comet/`, `data/`, `fig/` — four empty
  directories with no `.gitkeep` and no stated purpose. A non-ASCII directory name is a
  portability hazard.
- `.DS_Store` files present in both trees.
- `.vscode/settings.json` is untracked but not gitignored in `ztf-sso-query`.
- `eph.to_csv()` / `df_ztf.to_csv()` are called without `index=False`, so every CSV has
  a leading unnamed index column and must be re-read with `index_col=0` — a trap the
  analysis notebooks currently sidestep only by not reading those files back.
- `README.md` in `ztf-sso-query` is a single heading. `ztf-comet` has none.
- No tests of any kind in either tree.
- Outputs (FITS, `eph.csv`, `ztf.csv`, `photometry_*.csv`, figures) all land in one
  `data/<target>/` directory, contradicting the `results/` + `fig/` split in `CLAUDE.md`.

---

## 3. What should carry forward unchanged

The review is not all negative. These are sound and should be preserved through the merge:

- **The rate-based search-box sizing** in `generate_fits_urls` — deriving the query box
  from the target's own apparent motion over the step interval is the right approach,
  and correctly converts arcsec/hr → deg/day.
- **The footprint containment test** using `ra1..ra4`/`dec1..dec4` — a cheap, effective
  filter on IRSA's "overlaps" results.
- **The record-number retry** in `query_sso_ephemeris` — parsing Horizons' ambiguity
  error to recover the right apparition record is genuinely useful for periodic comets,
  and is the kind of hard-won detail that must not be lost. (Fix C12, keep the logic.)
- **`sep.winpos` centroid refinement** seeded from the WCS-projected ephemeris — correct
  handling of ephemeris uncertainty.
- **The Afρ formula and its unit handling** — dimensionally correct, and the
  `.to_value(u.cm)` conversion via astropy is the right way to do it.
- **The linear phase correction** — the standard reduction to A(0°)fρ.
- **The annotated-cutout figure design** — N/E, −V and −☉ vectors plus an ephemeris
  panel is exactly the right diagnostic plot. It needs to become one function (C20),
  not be redesigned.
- **`rcparams.py`** — a single shared style module is the right pattern; only the name
  and import path need fixing.

---

## 4. Proposed target structure

```
ztf-comet/                        →  github.com/<user>/ztfcomet
├── pyproject.toml                complete dependency set, requires-python >=3.10
├── README.md
├── CLAUDE.md                     rewritten to match reality
├── handoff.md
├── .gitignore                    excludes data/, results/, fig/, *.fits
├── ztfcomet/                     the package
│   ├── __init__.py               public API
│   ├── directory.py              THE path authority (C15)
│   ├── config.py                 targets, orbit IDs, solar mags, β, ρ_km, cutout size
│   ├── query.py                  Horizons + IRSA metadata  (C10, C12, C13, C14)
│   ├── cutout.py                 URL construction + validated download (C2)
│   ├── phot.py                   apertures, sky, calibration, Afρ (C1, C4–C8)
│   ├── plotting.py               annotated cutout + Afρ lightcurve (C20)
│   └── rcparams.py               (C15)
├── notebooks/
│   ├── main.py                   single/multi-target end-to-end driver
│   ├── query.ipynb               validate querying — test target 24P
│   ├── afrho.ipynb               validate photometry — test target 24P
│   └── figure.ipynb              cutouts + Afρ, single and multiple targets
├── doc/
│   └── primitive_code_analysis.md   (this file)
├── data/                         raw FITS + eph.csv/ztf.csv   (gitignored)
├── results/                      photometry tables            (gitignored)
└── fig/                          figures                      (gitignored)
```

---

## 5. Decisions needed before the upgrade

Each of these either invalidates results already produced, changes the scientific
output, or discards history. They are listed for explicit sign-off.

| # | Question | Recommendation |
|---|---|---|
| D1 | Fix C1 (zeropoint)? It invalidates every existing Afρ number and both lightcurve figures. | **Yes** — the current values are wrong by up to ×21. |
| D2 | Apply the `CLRCOEFF` colour term (C4)? Requires a comet colour, or per-epoch g−r where both filters exist. | Yes, with a configurable default and a flag when colour is unmeasured. |
| D3 | Reject or merely flag frames failing the edge/`rho_fwhm` cuts (C3, C6)? | Flag in the table, exclude from plots by default — keeps the data, protects the figures. |
| D4 | Start `ztfcomet` as a fresh repo, or graft `ztf-sso-query`'s 5 commits? | Fresh — it drops the 36 MB of tracked FITS (C17) cleanly. |
| D5 | Keep `example-data/` in the new repo? | No. Ship `fits_urls.txt` + `eph.csv` and let `main.py` re-fetch. |
| D6 | Drop `sxobsplan` and `astrometry` (C19, unused)? | Yes. |
| D7 | Delete `_construct_fitsurl.py` / `ztf-fitsurl.py` (C18)? | Yes — identical, dead, superseded. |
| D8 | Keep the legacy notebooks (`240P.ipynb`, `ztfquery_*.ipynb`)? | Archive under `notebooks/legacy/` outputs-stripped, or delete. |
| D9 | `rho_km` fixed at 15000, or per-target in `config.py`? | Per-target in config, defaulting to 15000. |
| D10 | Propagate uncertainties to Afρ (C8)? | Yes — `MAGZPRMS` + photon noise at minimum. |


---

## 6. Addendum — post-merge findings (2026-09-07)

Three defects found *after* the merge, while adding background-source
contamination flagging and fragment-safe target resolution. All are fixed and
carry regression tests; they are recorded here because each was silent.

### C23. `sep.winpos` returns three values, so centroiding never ran

`sep.winpos` returns `(x, y, flag)`. The merged `phot.measure_photometry`
unpacked two:

```python
xw, yw = sep.winpos(...)          # ValueError, every frame
```

The `ValueError` was caught by the surrounding `except Exception`, which fell
back to the unrefined WCS position. Centroid refinement was therefore **disabled
on every frame**, and the only visible symptom was a `centroid_shift_pix`
histogram that was exactly zero everywhere — easy to read as "the ephemeris is
excellent" rather than "the code never ran".

Measured on the 111-frame C/2024 E1 sample once fixed:

| | value |
|---|---|
| median centroid shift | 1.51 px (1.53″) |
| 90th percentile | 2.72 px |
| frames shifted > 1 px | 79 / 111 |
| shift as a fraction of aperture radius | ~25% |

A quarter of the aperture radius is not a rounding error: it moves flux out of
the aperture and sky into it.

**Lesson.** A bare `except` around a library call turns an API mismatch into a
silent no-op. The fallback path must be counted and reported, not merely taken —
`measure_photometry` now escalates to `log.error` when more than half the frames
fail to centroid.

### C24. Horizons record numbers are not stable

The pre-merge notebooks hardcoded `90001203` and `90001204`, both labelled
"240P/NEAT", and those numbers were carried into `config.py` during the merge.
Queried today:

| record | resolves to |
|---|---|
| 90001203 | **233P/La Sagra** |
| 90001204 | **234P/LINEAR** |
| 90001211 | 240P/NEAT (2014 solution) |
| 90001212 | 240P/NEAT (2024 solution) |
| 90001213 | 240P-B/NEAT (fragment) |

JPL renumbers its small-body records, so a pinned number is a cache key with no
expiry check. Any 240P reduction run with the old notebooks may be of a
different comet entirely.

**Fix.** `Target` now carries a *designation*; `ztfcomet.horizons` resolves it
per epoch and `verify_targetname` checks what came back. `orbit_records` remains
as an override but warns.

### C25. An ambiguous designation resolves to a fragment

Horizons answers `240P` with three candidates — two parent orbit solutions and
the fragment `240P-B`. The predecessor's `extract_lastrecnum` took the **last**
row, which is the fragment: 2.9 mag fainter than the parent and a different
object.

`select_record` now drops fragments unless `allow_fragment=True`, restricts
candidates to the requested parent designation, and prefers the orbit solution
nearest the observation. That last part is not cosmetic — 240P's two parent
solutions differ by ~50 arcsec in 2025, several times the aperture radius, so
the wrong one would put the aperture on empty sky.

### Also fixed

- **Broken notebooks passed validation.** The three notebooks were written with
  `source` lists built by `split("\n")`, which drops the newlines; Jupyter joins
  the list with `""`, so every cell collapsed to one unusable line.
  `nbformat.validate()` checks schema only and reported them valid. All three
  now execute clean under `nbconvert --execute`, which is the check that
  actually means something.
- **Mislabelled figure axes.** `plot_cutout` labelled axes "RA"/"Dec"
  unconditionally, including on subplots created without a WCS projection, which
  show pixels. It now labels by projection.

### C26. Contamination flagging (new capability, not a defect)

Background stars inside the aperture inflate Afρ indistinguishably from
activity. `phot.flag_contamination` sums Gaia DR3 sources within
`rho_pix + FWHM` into `G_eff = -2.5 log₁₀ Σ 10^(-0.4 Gᵢ)` and flags the frame at
≥30% of the comet's predicted `Tmag` flux.

On C/2024 E1 it flags 5 of 111 frames — and **5 of the 6 highest Afρ points**.
Contaminated frames have median A(0)fρ = 1365 cm against 654 cm for clean ones.
The sixth high point is separately flagged as an edge-truncated cutout, so every
outlier in the lightcurve is now accounted for.

Two limits belong with any result: the catalogue stops at `G = 18.5`, so
"uncontaminated" means "no *catalogued* source"; and Gaia `G` is compared
directly with a visual `Tmag`, a ~0.1–0.2 mag mismatch against a 0.28 mag
threshold.

### C27. A bare designation is not necessarily read as a comet

Horizons guesses the object class of an id it is given. Asking for `2P` with no
class hint returns **Styx (905)** — a moon of Pluto — with a complete,
plausible-looking ephemeris and no error, no warning, and no ambiguity listing.
The wrong-object failure is entirely silent; the only clue is that `Tmag` is
absent, because satellites do not carry cometary magnitude parameters.

```
Horizons(id="2P", ...)                        -> Styx (905)
Horizons(id="2P", id_type="smallbody", ...)   -> ambiguity listing, 61 records
                                                 for 2P/Encke (1786 .. present)
```

Every Horizons call in the package now passes `id_type="smallbody"`. Record
numbers are unaffected — they resolve identically with or without the hint.

This is the third variant of the same failure mode, after the fragment
substitution (C25) and stale record numbers (C24): **Horizons will answer a
question you did not ask, and the answer looks fine.** `verify_targetname` is
the backstop for all three, and it is the reason this one was caught.

### C28. `search_frames` read the retired `horizons_id`

When targets moved from record numbers to designations (C24), `Target.horizons_id`
became `None` for every entry — but `query.search_frames` still read it directly
and passed `None` to Horizons. Every query aborted with `'id' parameter not set`,
and because the failure was handled as "no ephemeris returned", a full run over
two comets completed successfully with **zero frames** and only a warning.

`search_frames` now resolves the designation for the mid-point of the query
window, logs the record it chose, and verifies `targetname` before searching.

A related gap closed at the same time: `query_sso_ephemeris` had its own
ambiguity-retry path that still used `extract_lastrecnum` — "take the last
record" — so a designation resolved *inside* that function would have picked the
fragment for 240P even though the resolver did the right thing everywhere else.
It now goes through `horizons.select_record` like every other path.
