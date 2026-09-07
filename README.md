# ztfcomet

Query and analyse comet observations from the **Zwicky Transient Facility**
(Palomar Observatory, MPC code `I41`): cutout retrieval from IRSA, aperture
photometry, and Afρ.

This repository merges two earlier projects — `ztf-sso-query` (cutout FITS
retrieval via the IRSA API) and `ztf-comet` (Afρ photometry and figures) — into
one package, and fixes the defects catalogued in
[`doc/primitive_code_analysis.md`](doc/primitive_code_analysis.md).

---

## Install

Requires Python ≥ 3.10. Development happens in the `spherex` conda environment.

```bash
conda activate spherex
pip install -e .
```

Raw FITS are **not** stored in the repository. Paths resolve in this order:

1. `$ZTFCOMET_DATA` if set,
2. the external SSD at `/Volumes/T7/data/ztf-comet` when mounted,
3. `<project>/data` otherwise.

```python
from ztfcomet import directory
print(directory.describe())
```

---

## Quick start

Query, download, reduce and plot a target in one command:

```bash
python notebooks/main.py 24P
```

Several targets, reusing data already on disk:

```bash
python notebooks/main.py 24P 240P 2P --steps phot figures
```

`--help` lists every override (date window, `rho_km`, cutout size, and so on).

From Python:

```python
import ztfcomet as zc

target = zc.get_target("24P")
eph, frames, report = zc.search_frames(target)          # Horizons + IRSA
urls = zc.build_urls(frames, cutout_size="10arcmin")
zc.download_urls(urls, zc.data_dir(target.name))        # validated download
phot = zc.run_photometry(target)                        # apertures -> Afrho
zc.plot_afrho({target.name: phot}, filters=["ZTF_r"], x="rh")
```

---

## Layout

```
ztfcomet/            the package
├── directory.py     THE path authority — nothing else builds a path
├── config.py        targets, orbit records, solar magnitudes, tunables
├── query.py         JPL Horizons ephemerides + IRSA image search
├── cutout.py        URL construction + validated download
├── phot.py          aperture photometry, calibration, Afrho
├── plotting.py      annotated cutouts, Afrho lightcurves
└── rcparams.py      shared matplotlib style

notebooks/
├── main.py          end-to-end driver, single or multiple targets
├── query.ipynb      validate the query stage (test target 24P)
├── afrho.ipynb      validate the photometry stage (test target 24P)
├── figure.ipynb     the core figures
└── legacy/          pre-merge notebooks, outputs stripped — do not reuse

doc/                 technical documents and the code review
data/                raw FITS + eph.csv/ztf.csv      (gitignored)
results/             photometry tables               (gitignored)
fig/                 figures                         (gitignored)
tests/               offline test suite
```

---

## Pipeline

**Query** — a coarse Horizons ephemeris decides where and when to look; the
IRSA search box is sized from the target's own apparent motion over the step
interval; each returned frame gets an exact-time ephemeris, joined **on JD**;
frames are kept only if the target falls inside the footprint corners.

**Download** — cutout URLs are built per frame and fetched with validation. IRSA
answers some requests with HTTP 200 and an HTML error body, so the payload is
checked for the FITS signature, written through a `.part` file, and logged in
`download_manifest.csv`. Re-running repairs anything corrupt.

**Photometry** — a circular aperture of physical radius `rho_km` at the comet,
centroided with `sep.winpos` from the WCS-projected ephemeris, on a
sigma-clipped annulus sky. Calibration:

```
filter_mag = inst_mag + MAGZP + CLRCOEFF·(g−r) + APCOR(ρ)
```

with uncertainties propagated from photon noise, `MAGZPRMS`, `CLRCOUNC` and
`ZPCLRCOV`.

**Afρ** — A'Hearn et al. (1984):

$$Af\rho = \frac{4\Delta^2 r_h^2}{\rho}\,10^{-0.4(m_c-m_\odot)},\qquad
A(0°)f\rho = Af\rho\cdot 10^{0.4\beta\alpha}$$

Solar magnitudes are PS1 AB (Willmer 2018), matching ZTF's calibration.

> Afρ is aperture-dependent by construction. Always quote `rho_km` with a value.

---

## Quality flags

Frames are **flagged, never dropped**. Every measurement stays in the table with
the reason attached, in a `flags` column plus one boolean per check.

| flag | meaning |
|---|---|
| `outside` | target falls outside the cutout |
| `aperture_edge` | aperture clipped by the array edge |
| `sky_edge` | sky annulus clipped by the array edge |
| `undersampled` | aperture narrower than `min_rho_fwhm` × seeing |
| `centroid` | `winpos` moved too far from the ephemeris position |
| `negative_flux` | background-subtracted sum ≤ 0 |
| `lowsnr` | SNR < 3 |
| `color_default` | *(advisory)* colour assumed, not measured |
| `nan_pixels` | *(advisory)* negative pixels clipped in the variance model |

`quality_ok` is the negation of the **critical** flags only. The two advisory
flags qualify a good measurement rather than invalidating it — excluding a whole
single-band night because its colour had to be assumed (≈0.015 mag) would throw
away usable data.

Plots draw flagged points as open symbols instead of hiding them; pass
`only_good=True` once you are satisfied the flags are right.

---

## Adding a target

Edit `ztfcomet/config.py` rather than pasting constants into a notebook — that
is how the old `ztfquery_2P.ipynb` came to annotate its figures with 24P's
ephemeris.

```python
"81P": Target(
    name="81P",
    horizons_id=90000664,
    start_date="2025-01-01",
    end_date="2025-12-31",
    perihelion_jd={"2028": 2461800.5},
    note="81P/Wild 2.",
),
```

An unknown name still works — `get_target` passes it straight to Horizons — but
the record number is what disambiguates apparitions for periodic comets.

---

## Tests

```bash
pytest tests/ -q
```

Offline and deterministic. The `test_regression_c*` cases pin the specific
defects from the review so they cannot return.

---

## Notes

- Do not read whole FITS files or large catalogues into an interactive session;
  print headers and shapes instead.
- Import `ztfcomet.rcparams` for figure style rather than restating rcParams.
  `savefig.dpi` is 200 for one-off figures; use 50 for batch output.
- `notebooks/legacy/` is kept for provenance only. Those notebooks do not run
  against this package and contain the bugs the review documents.
