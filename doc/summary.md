# SPHEREx comet catalog — scientific summary

_Project `spherex-comet` · written 2026-09-15 at the merge of the three former projects
(`ztf-comet`, `spherex-comet-apphot`, `spherex-comet-catalog`) · state of the results:
ZTF chain 2026-09-11/14, SPHEREx photometry 2026-09-12 (config `9fcc7ea3871a`), gas fits
2026-09-14 (run `a6463b4c610b`, `spherex_comspec` 1.2.0)._

This document is the one-stop description of the project: what question it answers, the
physics and the definitions it relies on, how the three-stage pipeline works, what the
numbers currently are, and what still has to be added, upgraded or questioned.  The
detailed technical notes it condenses are indexed in `doc/README.md`; the working state
for the next session is `handoff.md`; the conventions for code and documents are in
`CLAUDE.md`.

---

## 1. Purpose

SPHEREx surveys the whole sky in spectrophotometry between 0.75 and 5 µm.  Comets cross its
fields by chance, so every comet brighter than the survey limit is observed repeatedly, at
one wavelength per exposure, over weeks to months.  The project turns those cutouts into a
**catalog of gas production rates** — Q(H₂O), Q(CO₂), Q(CO) — for the 68 comets of the
working list (`data/reference/sx_comet_list_ver2607.xlsx`), each per *observing phase*
(a stretch of a few weeks in one physical state), together with the **dust activity
context** of every phase from ground-based ZTF imaging (the phase-corrected A(0°)fρ and
its heliocentric trend).  The products are the tables under `results/`, the figures under
`fig/`, and the review deck `doc/figures.pptx` (three slides per fitted phase).

The scientific goal is a homogeneous CO₂/H₂O and CO/H₂O census across dynamical classes
and heliocentric distances, with the dust production of the same epochs measured
independently, so that gas-to-dust behaviour and the CO₂-driven activity of distant
comets can be read from one consistent dataset.

## 2. Scientific background

**SPHEREx spectra of comets.**  SPHEREx disperses with linear variable filters, so a pixel
sees one narrow bandpass whose central wavelength depends on where the source falls on the
detector; a comet's spectrum is not one exposure but the collection of all exposures of a
phase, each a photometric point at its own wavelength (`wl`, width `wlwidth`), taken at its own geometry (r_h, Δ, phase angle); the resolving power is R ≈ 40–100.  The three stages
therefore have to (i) measure a flux per exposure inside a physically defined aperture,
(ii) group exposures into phases whose geometry is nearly constant, and (iii) fit the
assembled point spectrum with a model that knows each channel's own bandpass and geometry.
Pixels are 6.2″, so comae are barely resolved; the point-spread function is about one pixel
and broadens with wavelength.

**Gas: fluorescence emission bands.**  The coma is optically thin and collisionless over
almost all of its volume; vibrational levels are populated by absorption of solar photons
and re-emission (resonance fluorescence), so the emission per molecule is a fluorescence
efficiency, the g-factor [photons s⁻¹ molecule⁻¹], that scales with r_h⁻² and does not
depend on the local density.  The bands SPHEREx reaches are the H₂O ν₃ band at 2.7 µm
(2.55–2.9 µm), the CO₂ ν₃ band at 4.26 µm, the CO v(1–0) band at 4.67 µm, and the H₂O hot
bands ν₃−ν₂ (4.63 µm) and ν₁−ν₂ (4.85 µm) that sit on top of the CO band and carry ~3 % of
the water emission.  A Haser coma (constant outflow velocity v_g, photodissociation
lifetime τ ∝ r_h²) integrated over the aperture (the Yamamoto filling factor) links the
number of molecules in the aperture to the production rate, so the band flux is linear in
Q and the three species can be solved jointly.  This is the AKARI/IRC band technique of
Ootsubo et al. (2012), extended here with band *shapes* from a reconstruction of the GSFC
fluorescence line database (Villanueva et al. 2011–2013) and with the velocity dependence
of g(CO) (the Swings effect: pump lines shift against solar Fraunhofer lines with the
comet's heliocentric velocity).

**Dust: Afρ.**  The dust cross-section in the coma is measured on ZTF r- and g-band
images through the A'Hearn et al. (1984) quantity Afρ = (2Δr_h/ρ)² · 10^(−0.4(m_c − m_⊙)) · ρ,
which for a steady-state coma whose surface brightness falls as 1/ρ is independent of the
aperture radius ρ and proportional to the dust production rate divided by the outflow
speed.  Values are corrected to zero phase angle, A(0°)fρ, and are compared between comets
at fixed physical apertures (10 000 and 20 000 km).  How Afρ changes with r_h (the power-law
index, activity peaks relative to perihelion, breaks, outbursts, colour) is the dust
activity trend; evaluated at the epoch of a SPHEREx phase it is the dust context of that
phase's gas production rates.

**Star contamination.**  A comet drifts across the star field, so on some exposures a
background star lands inside the aperture.  Both the ZTF and the SPHEREx photometry test
every aperture against Gaia DR3 and *flag* rather than mask the contamination: masking
stars out of the aperture deletes real coma flux (the prototype did that, and it biased
Q(CO₂) low by 13 %).

## 3. Concepts and definitions

| term | meaning in this project |
|---|---|
| **epoch** | a block of SPHEREx exposures separated by ≥ 28 days from the next (the survey returns to a field roughly every six months) |
| **phase** | a subdivision of the epochs into a single observing state: r_h spread < 10 %, Δ spread < 20 %, no perihelion passage inside, emission bands never cut, cuts placed in gaps (`data/comspec/phase_assignment.csv`, `results/comspec/phase_map.csv`) — 193 phases over the 68 comets |
| **aperture** | circular, defined in km at the comet (10 000–100 000 km set, plus 2 and 5 px); the gas fit uses one fixed aperture per phase: 20 000 km inside 3 au, 40 000 km beyond, 60 000 km when that is under 1.5 px |
| **sky annulus** | SPHEREx: inner radius 150 000 km at the comet, floored at 15 px and capped at 40 px, 5 px wide; ZTF: 3ρ–(4ρ + 20 px), sigma-clipped median |
| **`badphot`** | a bad (flagged or non-finite) pixel inside the aperture; the flux is summed over the good pixels and the effective area recorded — not a flux-sign cut |
| **`sourceflag`** | Gaia contamination class of a SPHEREx measurement: `a` a G < 13 star within r_ap + 2 FWHM, `b` Gaia flux above 20 % of the comet's inside r_ap + FWHM, `c` any source inside r_ap + FWHM, `d` S/N < 1, `0` clean; priority a > b > c > d |
| **`quality_ok`** (ZTF) | none of the critical flags (`outside`, `aperture_edge`, `sky_edge`, `undersampled`, `centroid`, `negative_flux`, `lowsnr`, `contaminated`, `anomalous_bright`); advisory flags (`color_default`, `nan_pixels`) qualify without invalidating |
| **distance-corrected flux** | F × r_h² Δ², the flux the comet would show at 1 au from both Sun and observer; the continuum is fitted in this space and the emission divided back by each channel's own factor before the fit |
| **continuum verdict** | PASS / WARN / FAIL of the per-band polynomial continuum (bracketing and positivity hard, shape and cross-validation soft); "clean" detections sit on a PASS continuum |
| **tiers** | `detected` ≥ 3σ, `marginal` 1–3σ (value reported, 3σ limit quoted), `upper_limit` < 1σ, `negative_fit`, `not_covered`; **robust** = detected with n_eff ≥ 2 |
| **n_eff** | participation ratio of the per-channel Fisher information — how many channels really carry a species |
| **GLS** | generalised least squares: the continuum polynomial's coefficient covariance is propagated into a full data covariance, correlated across the channels of a band |
| **hot-band fallback** | when the 2.7 µm band is not covered and ≥ 3 channels lie in 4.55–4.90 µm with one beyond 4.75 µm, the H₂O hot bands carry Q(H₂O) (`h2o_source = hot`); offered only inside 3 au |
| **Afρ methods at a SPHEREx phase** | `direct` (≥ 2 clean ZTF frames within 5 d of the window, moved to ⟨r_h⟩ along the local law), `trend` (the fitted law of that orbital phase), `trend_extrap` (≤ 0.1 dex beyond it, grade A–C), `none` with the reason |
| **trend grades** | A–D reliability of a fitted power law `log Afρ = a − x log r_h`; grade-D laws are not drawn and are not extrapolated |
| **arc** | in `phase_map.csv` an arc *index* (`in` until a perihelion resolved inside the SPHEREx coverage), not the orbital direction — decide inbound/outbound from T_p |

## 4. Data

- **SPHEREx cutouts** (external SSD `/Volumes/T7/data/spherex-comet`): the flat archive
  `comets_v5_fits/` (~139 000 cutout FITS of 91 px with IMAGE / VARIANCE / FLAG planes, in
  mJy per pixel), its index `comets_v5/db_filtered.parq` (predicted V < 22, r_h < 20 au,
  with ephemeris geometry, `vmag`, corrected centres `xcen`/`ycen`), and the per-target tree
  `comets_v5_fits_filtered/<target>/` that `scripts/apphot/extract_targets.py` copies out of
  it for the working list — 27 797 exposures for 68 comets.  `data/spherex_sample/` holds
  2P and 24P (991 files) for the end-to-end test without the SSD.
- **ZTF frames** (`/Volumes/T7/data/ztf-comet/<target>/`): 5′ (10′ for bright or close
  comets) cutouts of every ZTF science image containing the comet since 2025-03-01 at
  V < 20 (V < 21 for the sparse targets), queried through JPL Horizons and IRSA — 8 730
  frames, download completeness 98.75 %; the 111 missing frames are permanent (archive 404s
  and cutouts that do not overlap the image).
- **Gaia DR3** (`~/Desktop/data/gaia_dr3/gaiadr3_all.npy`, 11 GB, complete to G ≈ 18.5)
  with the declination-sorted cache that turns a cone search into milliseconds.
- **Reference tables** (`data/reference/`): the working list, dynamical classes, literature
  production rates (Harrington Pinto et al. 2022; Gicquel et al. 2023 NEOWISE), the previous
  phase map for the regrouping regression test.
- **Fluorescence database** (`data/fluorescence/`): the project's reconstruction of the GSFC
  g-factors and band profiles for 14 species over 0.7–5.0 µm at five rotational
  temperatures, from HITRAN 2020 and public solar spectra, validated to ~10 % against the
  published values (`doc/comspec/fluorescence_database.md`).

## 5. The pipeline

```
 ZTF frames ──► ztfcomet ──────────────► results/ztf/photometry, profile, afrho
 (IRSA, Horizons)   query · download · aperture photometry · Gaia flags · radial
                    profiles · heliocentric Afρ trends · Afρ at the SPHEREx epochs
                                                          │  spherex_afrho.csv
 SPHEREx cutouts ─► spherex_apphot ────► results/apphot/photometry/<T>.csv  (input of stage 3)
 (T7, db_filtered)  masks · physical annulus · 21 apertures · Gaia sourceflags · stacks
                                                          │
                    spherex_comspec ───► results/comspec/gas_fit.csv  ◄── attach_afrho_ztf.py
                    phase grouping · continuum subtraction · joint H₂O/CO₂/CO fit
                    · study variants · figures · phase_images.py · make_figure_slides.py
```

### 5.1 Stage 1 — ZTF dust photometry (`ztfcomet`, `scripts/ztf/`)

*Query and download.*  A coarse Horizons ephemeris decides where and when to look; the
IRSA search box is sized from the comet's own motion; each frame gets an exact-time
ephemeris joined **on JD**; cutout URLs are fetched with validation (IRSA answers some
requests with HTTP 200 and an HTML body).  Designations are resolved by name — Horizons
record numbers are renumbered and fragments share the parent's designation.

*Photometry.*  Circular apertures of 10/15/20/30/40 × 10³ km at the comet, centroided with
`sep.winpos` from the WCS-projected ephemeris, sigma-clipped annulus sky, calibration with
the **per-frame** zeropoint, colour term and aperture correction, uncertainties from photon
noise and the calibration terms.  Every frame stays in the table with its flags; Gaia
DR3 contamination is flagged when the catalogued flux inside r_ap + FWHM exceeds 30 % of the
comet's expected V.  Afρ follows A'Hearn et al. (1984) with PS1 AB solar magnitudes and a
linear phase correction.

*Coma radial profiles.*  Surface brightness in 0.5 px annuli on a ×4 oversampled cutout,
centred on the optocentre, against a stack of field stars on the same frame; a nucleus + C ρ^(−m)
coma convolved with the empirical PSF plus a free sky is fitted per frame.

*Heliocentric trends and the SPHEREx epochs* (`afrho_trends.py`, `ztfcomet.activity`).
Outbursts and single-frame anomalies are set aside, phases are split at the activity peak
rather than at perihelion, `Afρ = A r_h^(−x)` is fitted per phase, aperture and band, a
broken law is tested by ΔBIC on errors inflated to χ²_ν = 1, dust colour is measured from
same-night g/r pairs, and for every SPHEREx phase the A(0°)fρ at the phase's ⟨r_h⟩ is
estimated (`direct` / `trend` / `trend_extrap`) and written to
`results/ztf/afrho/spherex_afrho.csv`, with one review figure per phase in
`fig/ztf/afrho/trend_slide/`.

### 5.2 Stage 2 — SPHEREx aperture photometry (`spherex_apphot`, `scripts/apphot/`)

Per target: the index rows are loaded with predicate pushdown, sorted and cut into 28-day
epochs; one Gaia cone search covers every pointing and feeds a KD-tree.  Per exposure the
IMAGE / VARIANCE / FLAG planes are read by name; the bad-pixel mask combines the flag bits,
non-finite science and non-positive variance (~16 % of pixels are NaN with no flag bit);
the aperture set is resolved for the exposure's pixel scale, keeping PSF FWHM ≤ r_ap < r_in;
the sky is the sigma-clipped annulus at 150 000 km; every valid aperture is measured in one
photutils call with the reduced effective area; the Gaia neighbourhood is summarised into
`sourceflag` and five columns.  Errors are reported as components (pixel, sky level, sky
scatter), as the formal total without the double-counted sky term, and as the empirical
annulus scatter (`source_sum_err_empirical_mjy`, the column the fits use); `sky_excess_ratio`
≈ 1.1–1.5 says the variance plane under-reports the scatter.  A distance-corrected flux, a
reflectance and optional per-epoch band stacks (sigma-clipped median, `results/apphot/stacks/`)
complete the table.  The production run of 2026-09-12: 68/68 comets, 27 797 exposures,
455 245 rows, 25 min at four workers, `results/apphot/status.csv` as the resume table.

### 5.3 Stage 3 — gas production rates (`spherex_comspec`, `scripts/comspec/`)

`main.py all` runs grouping → the main variant and nine study variants → cross-variant
analysis → figures (10 min).  Rows enter a spectrum unless `frac_badpix_ap > 0.05` or
`sourceflag = a`.  Per band a weighted polynomial continuum in distance-corrected flux is
fitted over the continuum window with the three emission windows punched out (H₂O
2.55–2.80 µm over 2.30–3.00; CO₂ over 3.90–4.65; CO 4.55–4.90), the order chosen by
leave-one-out cross-validation, an empty one-sided window extended by 1 µm, and validated
(PASS / WARN / FAIL).  The emission channels of the three bands are then fitted jointly by
GLS for Q(H₂O), Q(CO₂), Q(CO) with the Haser–Yamamoto model, the reconstructed g-factors and
band shapes, and g(CO) at the comet's heliocentric velocity; coverage is judged before
significance, the covariance is rescaled by χ²_ν when above 1, and every value carries its
tier, n_eff, caveats and the geometry it belongs to.  `attach_afrho_ztf.py` pivots the ZTF
table into `results/comspec/afrho_ztf.csv` and adds the `afrho_*` columns to `gas_fit.csv`
and `phase_map.csv` (the pipeline re-attaches them on every rerun and blanks moved phases
as `stale`); `phase_images.py` builds comet-centred, north-up median stacks of the ZTF frames
nearest each window and of the SPHEREx exposures in the continuum and the three emission
windows (`fig/comspec/phase_images/`, `results/comspec/phase_stacks/`); `make_figure_slides.py`
assembles the review deck.

## 6. Main results

### 6.1 ZTF (56 of the 68 comets have ZTF photometry)

Twelve comets have no usable ZTF data: four never reached V = 20 (2014 UN271, 2022 R3, 229P,
2025 UX109) and eight had no ZTF coverage at their positions.  The 41 739 measurement rows
(8 730 frames × apertures) are 51.2 % clean; contamination is the largest single cause of
rejection (one row in five), then low S/N and sky annuli reaching a frame edge.  20 000 km is
the best-sampled aperture (98.6 % of frames), 10 000 km the one most often below the seeing
disc.

*The 1/ρ law holds to the edge of coverage.*  With the PSF core and the sky modelled, the
observed coma profile matches a 1/ρ coma within 1–2 % from ~1 500 km out to 37 000–45 000 km
at every S/N band and from 0.6 to 6.7 au; the naive slopes of −1.4 to −1.8 on distant comets
were an S/N effect, and the corrected slope converges on m = 1.00 (median 1.04 at S/N > 60)
with a real comet-to-comet dispersion of a few tenths (0.57 for 2025 R2 to 1.42 for 10P).
Field stars give ≈ −4.4.  Aperture-ratio "steepening" of C/2024 E1 with distance is a sky
systematic of fixed-annulus photometry, not coma structure; the far annulus removes a
constant 8–12 % of the coma flux (½ ρ/r_med for m = 1) and is exactly correctable.

*Heliocentric Afρ trends.*  Only four comets have both a rising and a fading phase with
enough points (24P: r_h^−8.8 rising, r_h^−4.0 fading, both grade B); outbound-only fading
is uniform (median index ≈ 3, 16–84 % within 0.7–4.2) while inbound-only indices scatter from
−0.2 to 7.  Eleven legs prefer a broken law, and the survivors are regime changes (10P flat
from 3.6 to 1.85 au then r_h^−10; C/2024 E1 an inbound maximum at 3.4 au; 2023 R1 flat beyond
3.9 au then steep); 3 au is not special.  217P's outburst at 3.10 au (×10, decaying over forty
days) is the largest of seven windows on five comets.  The dust is +0.10 mag redder than the
Sun with no r_h trend.  Of the 14 comets sampled on both sides of perihelion, eight have an
interior activity maximum and are split there; the bracketed peaks fall within about a month
of perihelion on either side (24P +11 to +13 d, 2023 R1 +17 to +19 d, 2025 K1 +36 to +38 d,
47P −17 to −34 d, 235P −7 to −26 d, 240P −20 d), and 2024 G4 peaked 60–85 d before.

*Afρ at the SPHEREx epochs.*  Over the 165 SPHEREx phases of the 56 comets, an A(0°)fρ
value exists for 116 at 10 000 km and 120 at 20 000 km (79 direct in each, the rest from a
law); the missing ones fall outside a grade-D law (19–20), on the orbital leg ZTF never
covered (15–17), or beyond 0.1 dex of the data.  Attached to the catalog, 89 / 95 of the 147
fitted phases carry a dust value.

### 6.2 SPHEREx photometry

All 68 targets processed (`9fcc7ea3871a`).  Median 7.9 % of exposures per target are
contaminated by the growth-curve criterion; against that independent symptom flag `a`
alone has the best lift (5.6, keeping 95 % of the data) while `a + b` has full recall at
44 % of the data kept — hence *cut on `a`, keep `b` and `c` as columns*.  `badphot` is
anti-correlated with contamination and must not be used as a contamination cut.  The
physical annulus of 2026-09-12 changed the flux of close comets by 2–6 % at 20 000 km (the
fixed 15–20 px ring had sat inside the coma of 24P at 40 000 km and removed 3.5–9 % of the
clean-channel flux) and left comets beyond ~1.6 au unchanged.  The known systematic that
remains is the absence of an aperture (encircled-energy) correction: the smallest retained
apertures are ~1 px against a ~1 px PSF that broadens with wavelength.

### 6.3 Gas production rates (`results/comspec/gas_fit.csv`, run `a6463b4c610b`)

161 of the 193 phases were analysable (32 skipped for insufficient channels), 147 were fitted
(14 not fittable), over 67 comets.  **Robust ≥ 3σ detections (n_eff ≥ 2): H₂O in 17 phases
(14 comets, all from the 2.7 µm band), CO₂ in 26 (25 comets), CO in 4 (3 comets)** — 35
comets with at least one robust detection, 34 of them clean (on a PASS continuum); a further
19 comets have only marginal (1–3σ) values (H₂O 25, CO₂ 23, CO 9 phases).  The apertures in
use are 20 000 km for 70 fitted phases, 40 000 for 61 and 60 000 for 16.  The median χ²_ν is
4.6 (10.4 at 20 000 km, 4.0 at 40 000, 2.3 at 60 000): the bright, close comets at small
apertures in pixels are limited by the line-spread function and the band shape, not by
noise, and since the errors are rescaled by √χ²_ν every tier count is conservative.
Q(H₂O) from the hot bands is provisional (23 rows, 1 detected, 6 marginal, all inside 3 au).

The values are stable against the alternatives: Q ratios of 1.01 / 1.01 / 1.00 (H₂O / CO₂ /
CO) against the 2026-09-12 run and 0.99 / 1.01 / 1.00 against the S/N-driven apertures;
the flag-`b` cut, the strict `badphot` rule and the physical-vs-distance-corrected
continuum space each leave Q unchanged (nothing beyond 1.4σ for the flux space).  Two
choices did move the numbers, and both are physics: the reconstructed band shapes raised
Q(CO₂) by 1.27 and Q(CO) by 1.51 over the Gaussian placeholders (which were 2–3× too narrow
and biased the peak-sampled bands low) and improved χ²_ν in 95 of 138 groups; the revised
photometry raised Q(CO₂) by 1.13 over the prototype's masked photometry.  The census, by
contrast, depends on the error model: diagonal errors add 5–11 "detections" per configuration
but produce three to six times the Gaussian expectation of fits at ≤ −2σ, which no coma can
produce, so GLS stays and its census is the honest one.  The 54-run method matrix
(`results/comspec/studies/method_matrix/matrix.csv`) chose the 2.55 µm emission edge, the
2.30–3.00 µm continuum, cross-validated orders and the one-sided extension by the same
criterion (clean robust detections against the negative tail).

### 6.4 Fluorescence database

Line lists, band g-factors and spectral profiles of H₂O, CO₂, CO and eleven organic,
nitrogen- and sulphur-bearing species at 30–130 K, reconstructed with the General
Fluorescence Model from HITRAN 2020 and the Kurucz/Toon solar spectra: H₂O ν₃ 3.14 × 10⁻⁴,
CO₂ ν₃ 2.71 × 10⁻³, CO v(1–0) 1.92 × 10⁻⁴ photons s⁻¹ molecule⁻¹ at 1 au and v_h = 0;
only CO responds to the heliocentric velocity (+25–31 % beyond 10 km s⁻¹).  Validated to
~10 % line by line against Villanueva et al. (2011, 2012, 2013); NH₃ is the least accurate
(−20 %), CH₃OH and HC₃N are band models.

## 7. Known systematics and concerns

1. **Line-spread function and band shape** (placeholder 1).  A Gaussian of FWHM = `wlwidth`
   stands in for the as-built LVF response; the bright comets with 2–4 channels on the CO₂
   band reach χ²_ν of 10²–10³, and the CO / H₂O-hot separation at 4.7 µm depends on it.
2. **No aperture correction** in the SPHEREx photometry (review item S3); the catalog
   sidesteps it by requiring r_ap ≥ 2 PSF FWHM for the analysed apertures only, so a
   residual wavelength-dependent aperture loss remains in every spectrum.
3. **Star contamination inside the fits**: 120 of the 147 fits contain flag-`b` channels
   (median 4); the finding that dropping `b` leaves Q unchanged was made at the earlier
   apertures and is re-tested by `dc_main_no_b`.  A star-free selection is also what
   decides the aperture scoring: an S/N score over flagged channels picks the most
   contaminated aperture (2024 N1).
4. **The H₂O continuum spans the band-3/band-4 detector boundary at 2.42 µm**; a
   calibration step between detectors would bend the continuum under the 2.7 µm band, but a
   band-4-only blue side is too short to validate.
5. **CO is at the noise level**: 4 robust detections against 5 of 92 covered fits at ≤ −2σ;
   the CO hump is degenerate with a common-mode offset of the continuum interpolated across
   4.55–4.90 µm.
6. **Conventions not yet varied**: the aperture thresholds (3 au, 1.5 px, 60 000 km), the
   annulus radius (150 000 km; 100 000 and 200 000 km untested), the Δ tolerance (20 %), the
   expansion-velocity law v_g = 0.8 r_h^−0.5 km s⁻¹ (±18 % in Q per 20 % in v_g), T_rot = 70 K,
   the opacity calibration.
7. **Errors remain optimistic** by ~1.5 in scatter even with the empirical column
   (`sky_excess_ratio` ≈ 1.2; `source_sum_err_mjy` is a lower bound); the χ² rescaling carries
   this into the quoted Q errors.
8. **Flag-bit meanings** of the SPHEREx FLAG plane were recovered empirically (bit 21 fires
   on the comet itself and must never be masked); confirm against the Explanatory
   Supplement before publication.  Gaia positions carry no proper motion (negligible at
   6.2″ pixels).
9. **ZTF side**: the coma-in-annulus factor (8–12 %) and an advisory `sky_dominated` flag
   (sky/flux > 10) are decided but not yet applied; a single-frame anomaly cannot be told
   from a one-night brightening under sparse sampling (29P); 23 % of profile fits sit on
   the sky bound at low S/N; `figure.ipynb` still writes the old per-target figure names.
10. **Provenance of the figures**: the per-comet SPHEREx summary figures in `fig/apphot/`
    date from the 2026-09-07 photometry (the per-exposure cutout PNGs, 7.6 GB, were deleted
    on 2026-09-15); regenerate with `make_figures.py --no-cutouts` for the current run.

## 8. To be added, upgraded or modified

**Highest priority — the physics placeholders.**  Replace the Gaussian LSF with the
as-built SPHEREx spectral response (the largest systematic of the bright comets); cite and
adopt a species-specific expansion-velocity law; decide, on the bright comets (24P, 306P,
508P), the annulus radius and the aperture thresholds by varying them; confirm the flag-bit
table; implement the aperture (encircled-energy) correction as a function of wavelength
and pixel scale, or demonstrate its size on the retained apertures.

**Catalog decisions still open.**  Whether marginal (1–3σ) values are quoted or only
limits; which ZTF aperture defines the dust context (10 000 km for literature comparison,
20 000 km for sampling); the treatment of the 2.42 µm detector step; the re-test of flag `b`
at the fixed apertures; the case 2025 K1 (grouping rule 3 reads the raw `badphot`).

**Extensions.**  Literature comparison of Q against the NEOWISE (Gicquel et al. 2023) and
Harrington Pinto et al. (2022) tables already in `data/reference/`; mixing ratios versus
dynamical class and r_h with the dust context (gas-to-dust from Afρ and Q); the coma-in-annulus
correction and `sky_dominated` flag in `ztfcomet.phot`; the additional species of the
fluorescence database (CH₄, C₂H₆, CH₃OH, …) as upper limits where SPHEREx covers them;
reflectance spectra of the dust continuum per phase (the stage-2 product exists,
unexploited downstream); a per-phase figure of the fitted spectrum with the ZTF frames for
the manuscript (the deck already holds it).

**Housekeeping.**  Regenerate `fig/apphot/` summaries and the review deck after the next
run; commit the monorepo to a remote and retire the two old repositories (their histories are
merged: `git log spherex-comspec-import`); unify the two figure styles (ztfcomet 15 pt
versus the project's 20 pt `notebooks/rcparams.py`); merge the outstanding ZTF branches
(`analysis/*`, `fix/*`) if anything in them is not on the current branch; keep
`results/comspec/placeholders.csv` current with every decision.

## 9. Provenance

| product | run | date | hash / version |
|---|---|---|---|
| ZTF photometry, profiles, Afρ trends | `scripts/ztf/survey.py`, `profile_survey.py`, `afrho_trends.py` | 2026-09-11 (chain), 2026-09-14 (SPHEREx epochs) | `ztfcomet` 0.3.0, branch `feat/afrho-peak-fits-and-reorg` |
| SPHEREx photometry | `scripts/apphot/main.py --target-list … --workers 4 --force` | 2026-09-12 | `spherex_apphot` 2.2.0, config `9fcc7ea3871a` |
| gas production rates | `scripts/comspec/main.py all` | 2026-09-14 | `spherex_comspec` 1.2.0, run `a6463b4c610b` |
| fluorescence database | `notebooks/comspec/fluorescence_gfm/build_fluorescence_db.py` | 2026-09-11 | `gfm-2026-09-11` |
| review deck | `scripts/comspec/make_figure_slides.py` | 2026-09-14 | `doc/figures.pptx` |
