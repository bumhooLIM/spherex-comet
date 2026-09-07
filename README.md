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

The Gaia DR3 catalogue used for contamination flagging resolves from
`$ZTFCOMET_GAIA`, defaulting to `~/Desktop/data/gaia_dr3`. It needs
`gaiadr3_all.npy`; the `gaiadr3_deccache/` sibling (sorted by declination) turns
a cone search from a full 11 GB scan into a 3–6 ms lookup and is used
automatically when present. `GaiaCatalog.build_dec_cache()` creates it.

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
├── gaia.py          Gaia DR3 cone search + contamination test
├── horizons.py      designation -> orbit record, fragment-aware
├── orbit.py         perihelion elements, Kepler t(r_h) for the date axis
├── profile.py       radial SB profile of the comet vs field stars (1/rho test)
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

**Contamination** — a comet drifts across the star field, so on some frames a
background star lands inside the aperture and inflates Afρ in a way the image
alone cannot distinguish from activity. Every catalogued Gaia DR3 source within
`rho_pix + FWHM` is summed into one effective magnitude

$$G_\mathrm{eff} = -2.5\log_{10}\sum_i 10^{-0.4G_i}$$

and compared with the comet's predicted `Tmag` from JPL. The frame is flagged
when the background carries ≥30% of the comet's flux ($G_\mathrm{eff} \le
T_\mathrm{mag} + 1.31$). On the C/2024 E1 sample this flags 5 of 111 frames —
and 5 of the 6 highest Afρ points.

> The catalogue is complete only to `G = 18.5`, so "uncontaminated" means "no
> *catalogued* source". Gaia `G` is also compared directly with a visual `Tmag`;
> the passbands differ by ~0.1–0.2 mag, small against a 0.28 mag threshold.

**Coma radial profile** — aperture photometry says how much light is in the
coma, not how it is distributed. `ztfcomet.profile` measures the comet's
surface brightness in annuli at 0.5–10 px (0.5 px step, sigma-clipped mean,
with the plain mean kept as a check on the clipping), on a ×4 bilinear
oversampled cutout so 0.5 px annuli on 1″ pixels are well sampled. Up to 20
unsaturated, isolated field stars with S/N > 10 are profiled on the same frame
and stacked (sigma-clipped median) as the PSF reference. A power law is fitted
outside the core: **−1 is a steady-state coma; stars give ≈ −4.** Output:
`results/<target>/profile_*.csv`, `fig/<target>/profile/` and a per-target
summary. On 24P the clean-frame median slope is −1.07.

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
| `contaminated` | catalogued background source(s) inside the aperture |
| `color_default` | *(advisory)* colour assumed, not measured |
| `nan_pixels` | *(advisory)* negative pixels clipped in the variance model |

`quality_ok` is the negation of the **critical** flags only. The two advisory
flags qualify a good measurement rather than invalidating it — excluding a whole
single-band night because its colour had to be assumed (≈0.015 mag) would throw
away usable data.

Plots draw flagged points as open symbols instead of hiding them; pass
`only_good=True` once you are satisfied the flags are right.

---

## Target resolution: fragments and stale records

Two traps sit between a designation and the right ephemeris, and both silently
return the **wrong object**.

**Fragments share the parent's designation.** Asking Horizons for `240P` returns
an ambiguity listing containing two parent solutions *and* `240P-B`, a fragment
about 2.9 mag fainter. Taking the last record — the obvious rule, and the one
the predecessor used — selects the fragment.

**Record numbers are not stable.** Horizons renumbers its small-body records.
The numbers the pre-merge notebooks hardcoded for 240P (`90001203`, `90001204`)
today resolve to **233P/La Sagra** and **234P/LINEAR**.

`ztfcomet.horizons` therefore resolves by designation, drops fragments, picks
the orbit solution nearest the observation, and verifies `targetname` on the way
back. For 240P the two parent solutions differ by ~50 arcsec in 2025 — more than
the aperture — so the epoch choice is not cosmetic:

```python
zc.get_target("240P").resolve_orbit_record(2458300.5)   # 2018 obs -> 90001211
zc.get_target("240P").resolve_orbit_record(2460900.5)   # 2025 obs -> 90001212
zc.get_target("240P-B").resolve_orbit_record(2460900.5) # the fragment, explicitly
```

Pass `--allow-fragment` (or use the `240P-B` target) when the fragment really is
the object of interest.

## Adding a target

Edit `ztfcomet/config.py` rather than pasting constants into a notebook — that
is how the old `ztfquery_2P.ipynb` came to annotate its figures with 24P's
ephemeris.

```python
"81P": Target(
    name="81P",
    designation="81P",          # never a record number: they get renumbered
    start_date="2025-01-01",
    end_date="2025-12-31",
    perihelion_jd={"2028": 2461800.5},
    note="81P/Wild 2.",
),
```

An unknown name still works — `get_target` resolves it against Horizons — but an
entry here is what makes the choice reviewable.

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
