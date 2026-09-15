# Pipeline decisions — SPHEREx comet gas production rates

**Date:** 2026-09-14 (first written 2026-09-09) · **Package:** `spherex_comspec/`
(`spherex_comspec` 1.2.0, run hash `a6463b4c610b`) · **Input:** `results/apphot/photometry/` (revised
`spherex_apphot` photometry, config `9fcc7ea3871a`; 68 comets, 27 797 exposures) ·
**Result:** `results/comspec/gas_fit.csv`.

This is the record of what the pipeline does, which alternatives were tried, why the
current configuration was kept as the most robust one, what remains a concern, and
which placeholders still have to be replaced.  The physics is in `model_concept.md`
and `fitting_methodology.md`; the photometry comparison is in `apphot_comparison.md`.

---

## 0. The pipeline, and the conclusion

Revised photometry → one aperture per comet → phases of a single observing state →
per-band polynomial continuum in distance-corrected flux → joint weighted linear
least squares for Q(H₂O), Q(CO₂), Q(CO) in physical flux space, with the 2.7 µm band
anchoring H₂O and the 4.6–4.9 µm hot bands taking over only where 2.7 µm is not
covered.  Six configurations were run end to end and compared group by group
(§7).  **The configuration kept as the main result (`dc_main`) is:**

| choice | value | why it won |
|---|---|---|
| photometry | revised set (`results/apphot/photometry/`) | the previous set zeroed bad and star-adjacent pixels inside the aperture, which biased Q(CO₂) low by 13 % and cost 8 robust CO detections (`apphot_comparison.md`) |
| `badphot` rule | drop a row only when `frac_badpix_ap > 0.05` | the strict rule (any bad pixel) removed 52 % of 24P's rows for no change in Q; the lenient rule analyses 154 instead of 150 groups and fits 138 instead of 135, with Q identical (§7.2) |
| source flags | drop `a` (G < 13 star within r_ap + 2 FWHM), keep `b`, `c`, `d` | dropping `a` costs 3 band-rows and nothing in Q; dropping `b` removes 28 fits and a third of the channels, again without moving Q (§7.1) |
| continuum flux space | distance-corrected (F × r_h² Δ²) | the channels of one group span up to 166 % in r_h² Δ²; the corrected continuum removes that gradient and changes Q by ≤ 1.35σ against the physical-space fit (§7.3) |
| error column | `source_sum_err_empirical_mjy` | the formal error under-reports the annulus scatter; the empirical one leaves Q unchanged and lowers χ²_ν from 2.44 to 2.27 (§7.4) |
| H₂O coverage | 2.7 µm main band, hot-band fallback | the bright, close comets (24P at perihelion, 306P) lose the 2.7 µm channels; the fallback recovers their Q(H₂O) and removes the CO bias those groups carried (§7.5) |
| emission model | reconstructed fluorescence database + velocity-dependent g(CO) (2026-09-11) | physically validated band strengths and shapes; the Gaussian placeholders were 2–3× too narrow and biased Q(CO₂) and Q(CO) low by 20–50 % at SPHEREx sampling; χ²_ν improves in 95 of 138 groups (§7.7) |
| photometry (2026-09-12) | regenerated with a physical sky annulus (150 000 km, 15–40 px) and full-precision `jd_utc` (`9fcc7ea3871a`) | the fixed 15–20 px ring sat inside the coma of close comets; Q changes by +4 % for the groups inside 1 au and not at all beyond ~1.6 au (§7.8) |
| grouping (2026-09-12) | Δ spread < 20 % (rule 1b) | 193 phases instead of 174; r_h²Δ² no longer spreads by up to 183 % inside a group |
| aperture (2026-09-14) | one fixed aperture per phase: 20 000 km inside 3 au, 40 000 km beyond, 60 000 km when the rule radius is under 1.5 px | 99 / 73 / 21 phases; Q within 1 % of the S/N-driven apertures of 2026-09-12, census within a few groups (§3, §7.9) |
| continuum (2026-09-14) | 2.55–2.80 µm emission over a 2.30–3.00 µm continuum, CO₂ continuum 3.90–4.65 µm, one-sided windows extended to 1 µm, order by CV (10 %) | the method matrix's best clean census; 71 % of the band-rows saved, 78 % of the continua linear (§5, §7.9) |
| fit (2026-09-12/14) | GLS with the continuum covariance, ≥ 3σ tier + marginal tier, no negative cut, hot-band cap at 3 au | the negative tail of the fits (≤ −2σ: 5 of 113 H₂O, 5 of 92 CO) matches the Gaussian expectation only with GLS; diagonal errors are 3–6× too optimistic (§6, §7.9) |

It analyses every one of the 68 comets, fits 147 (comet, phase) groups, and yields a
**robust ≥ 3σ detection** (n_eff ≥ 2) of at least one species for 35 comets — 34 with a
*clean* one, i.e. on a PASS continuum: H₂O 17 groups (all clean, all from the main band;
the hot-band fallback stops at 3 au), CO₂ 26 (25 clean), CO 4 (2 clean); a further 19
comets have only marginal (1–3σ) values (H₂O 25, CO₂ 23, CO 9 groups).
The alternatives remain runnable as study variants (`config.DEFAULT_VARIANTS`, products
under `results/comspec/studies/`), so every number here can be regenerated.

## 1. Inputs

`results/apphot/photometry/` carries one row per (exposure, aperture) with the revised
treatment: predicted ephemeris centres, a fixed 15–20 px sky annulus, bad pixels
masked with an effective-area correction, Gaia sources flagged (`sourceflag`
`a`/`b`/`c`/`d`) rather than masked, a distance-corrected flux column and both a
formal and an empirical error.  Two properties matter downstream: `badphot` means
*a bad pixel in the aperture*, not a negative flux; and `jd_utc` is rounded to
0.1 day, so exposures are keyed on `filename` and timed from `date_obs`.

## 2. Phase grouping (`grouping.py`)

The 28-day `epoch` blocks are cut into phases that hold one physical state: r_h
spread < 10 %, a perihelion passage split (resolved when both sides move by
> 0.001 au and have ≥ 2 epochs), emission bands never cut, the fewest groups with
cuts placed in gaps (link 0.05 au), manual edges for 24P and 2024 E1.  Result: 141
epochs → **174 phases**, identical to the previous notebook's map for 66 of 68
comets (2025 K1 and 2025 L1 differ because rule 3 reads band sampling through
the new `badphot`).  `results/comspec/phase_map.csv`, `phase_cuts.csv`.

**Rule 1b (2026-09-12).**  The observer distance is bounded too: `(max − min)/mean`
of Δ below 20 % inside every group (`GroupingConfig.delta_tol`), including the
hand-set r_h bins of 24P and 2024 E1, which are subdivided by Δ alone.  Before,
r_h² Δ² spread by up to 183 % inside a group (concern 3); 20 of the 174 groups
of the 2026-09-11 run violated the new rule.  `delta_tol=None` reproduces the old map.

## 3. Aperture (`dataio.aperture_table`)

**One aperture per phase (2026-09-14, `ApertureConfig.rule = "fixed"`).**  20 000 km when
the phase's median r_h is inside 3 au, 40 000 km beyond; when that radius is smaller than
1.5 pixels at the phase's median pixel scale it is enlarged to 60 000 km.  Nothing else
enters the choice.  The upstream photometry refuses only apertures below the PSF FWHM
(0.86–1.08 px), so the rule aperture is measured for every exposure of a phase except for
2014 UN271 (Δ ≈ 14–15 au, where 60 000 km is 0.9 px; coverage 50 %).  Over the 193 phases:
99 at 20 000 km, 73 at 40 000 km, 21 enlarged to 60 000 km (all beyond 3 au, Δ > 5.9 au);
22 comets use two apertures across their phases.  The model's `rho_ap_km` follows the
phase.  `results/comspec/apertures.csv` records the radius, the reason (`near`, `far`,
`far->enlarged`), the coverage, the phase's median geometry and the radius in pixels.

The rule replaces the S/N-driven per-target aperture of 2026-09-12 (kept as
`rule="snr"`, study variant `dc_ap_snr`), which chose 10 000–100 000 km and coincides
with the fixed rule for 12 of the 193 phases; §7.9 gives what the change did (Q within
1 %, two CO and two CO₂ robust detections lost, two H₂O and one CO₂ gained).  The r_h
rule before that (per target, promoted for coverage) survives as `rule="rh"`
(`dc_rules_previous`).

## 4. Row selection (`config.BASELINE_FLAGS`)

A row enters the spectrum unless `frac_badpix_ap > 0.05` or `sourceflag = a`.
Nothing is written into the photometry; the selection and its rejection counts
travel with every product.

## 5. Continuum subtraction (`continuum.py`)

Per band a weighted polynomial (order by cross-validation up to 3; 3 / 2 / 2 before 2026-09-12) in λ − λ_c
over the continuum window with all three emission windows punched out, 3σ
clipping, order capped at 1 when the continuum is one-sided, a two-point fallback
below six points.  Validation: bracketing and positivity are hard, shape and
cross-validation soft → PASS / WARN / FAIL.  Main run (2026-09-14): 483 band-rows, 290 PASS /
52 WARN / 141 FAIL (71 % saved); 32 groups skipped for insufficiency.  The emission
is divided back by each channel's own `distcorr_factor` before the fit
(`emis_raw_mjy`), which is what keeps Q at the comet.

**2026-09-12 / 2026-09-14 revisions.**  (i) Windows: the 2.7 µm emission window is
2.55–2.80 µm, where the H₂O template starts (the 2.60 µm edge cut 6.7 % of the band; the
2.50 µm edge of 2026-09-12 added signal-free channels to the fit), over the 2.30–3.00 µm
continuum (the 2.20–3.10 µm window of 2026-09-12 cost one to two clean H₂O detections);
the CO₂ continuum is 3.90–4.65 µm (4.00–4.55 µm gives the same census).  Both were chosen
by the method matrix of §7.9.
(ii) When the window holds continuum on one side of the band only, the empty side
is extended to 1 µm from the band edge before fitting (every emission window is
still punched out); the fit is then bracketed and the order is no longer capped
at 1, and the summary records `cont_extended` and the window used.  (iii) The
polynomial order of every bracketed fit is chosen by leave-one-out
cross-validation among orders 1–3, the lowest order within 10 % of the best CV
error (`ContinuumConfig.order_mode = "cv"`): in the 2026-09-11 run the CV
preferred order 1 in half of the fits that used the fixed 3/2/2.

## 6. Fitting (`fitting.py`)

Optically thin Haser coma with the Yamamoto filling factor, g-factors ∝ r_h⁻², so
the model is linear in Q and the three species are solved jointly by weighted
least squares over every accepted emission channel, each at its own geometry and
bandpass.  **Since 2026-09-11 the g-factors and band shapes are the project's
reconstruction of the GSFC fluorescence database** (`data/fluorescence/`,
`doc/comspec/fluorescence_database.md`; validated to ~10 % line by line against the published
values) and **g(CO) follows the comet's heliocentric velocity** (Swings effect: 1.92 × 10⁻⁴ s⁻¹
at v_h = 0, 25–31 % more for |v_h| ≳ 10 km/s), with v_h of every pointing derived from the
ephemeris r_h(t) of the photometry (agrees with JPL Horizons to 0.05 km/s at worst).  The
previous model — eight Gaussian bands with the Ootsubo et al. (2012) g-factors and a constant
g(CO) = 2.6 × 10⁻⁴ — is the `dc_main_gauss` study variant; §7.7 gives what changed.  Coverage is decided on the diagnostic ranges (H₂O ≥ 3 channels in
2.60–2.80 µm, CO₂ ≥ 2 in 4.20–4.30, CO ≥ 2 in 4.60–4.70); *not covered* is
categorically different from *not detected*.  The covariance is rescaled by χ²_ν
when > 1; `n_eff`, the participation ratio of the per-channel Fisher information,
says how many channels really carry a species.  Statuses: `detected` (≥ 3σ since 2026-09-12; `marginal` 1–3σ),
`upper_limit`, `negative_fit`, `not_covered`; "robust" throughout this document
means detected with n_eff ≥ 2.

**H₂O hot-band rule (2026-09-09).**  When the 2.7 µm band is not covered and at least
three channels lie in 4.55–4.90 µm with one beyond 4.75 µm (the 4.85 µm ν₁−ν₂ band,
outside CO v(1−0)), the ν₃−ν₂ and ν₁−ν₂ hot bands carry Q(H₂O).  The row is
labelled `h2o_source = "hot"`, the caveat says so, the figures draw it open.  Such
values are provisional: the hot bands carry ~3 % of the water emission, and the
reconstruction makes the 4.85 µm band 40 % weaker than the harmonic estimate used before.

**2026-09-12 revisions.**  (i) *Generalised least squares*: the continuum-model
uncertainty is correlated across the channels of a band through the polynomial
coefficients; their saved covariance now builds the full data covariance
(`fitting.data_covariance`, `FitConfig.gls`), whose diagonal equals the errors
used before.  (ii) *Tiers*: `detected` means ≥ 3σ, `marginal` 1–3σ (value
reported with its error, 3σ limit quoted, `Q_X_nsig` in the table), `upper_limit`
below 1σ; limits are 3σ.  (iii) *No negative-channel cut*: every channel enters the
solve (the one-sided 1σ cut biased Q upward and could empty a group).  (iv) *Hot-band
cap*: the hot-band fallback for Q(H₂O) is offered only inside 3 au
(`FitConfig.h2o_hot_max_rh_au`).  `dc_main_diag` reruns the main configuration
with diagonal errors; `dc_rules_previous` with every 2026-09-11 rule a variant
can carry (r_h aperture rule, fixed 3/2/2 orders without the extension, the 1σ cut
and tier, diagonal errors, no cap).

## 7. What the studies showed

All comparisons are group by group against `dc_main`; "Q ratio" is the median
over groups detected in both runs with its 16–84 % range; "max" is the largest
shift in units of the larger error.  Tables: `results/comspec/studies/flag_policy_*.csv`,
`distcorr_effect_*.csv`, `apphot_comparison/`.

### 7.1 Source flags

| run | groups / fits | robust H₂O / CO₂ / CO | Q ratio | max |
|---|---|---|---|---|
| `dc_main_keep_a` (flag `a` kept) | 155 / 138 | 40 / 21 / 18 | 1.000, 1.000–1.000 (all) | 1.7σ |
| `dc_main_no_b` (flag `b` dropped too) | 138 / 110 | 29 / 13 / 11 | 1.000; H₂O 0.64–1.03, CO₂ 0.90–1.06, CO 1.00–1.21 | 3.0σ |

Flag `a` (342 of 12 438 channels) changes nothing; it is dropped because it marks
the one contamination mechanism the photometry project demonstrated.  Flag `b` is
too common to cut on: it removes a third of the channels and 28 fits and leaves
the surviving Q unchanged.

### 7.2 `badphot` rule

| run | groups / fits | robust H₂O / CO₂ / CO | Q ratio | max |
|---|---|---|---|---|
| `dc_main_strict` (any bad pixel drops the row) | 150 / 135 | 36 / 25 / 10 | 1.000; 0.98–1.14 / 0.99–1.10 / 0.96–1.15 | 4.1σ |
| `dc_all` (strict + every flag; the 2026-09-08 baseline) | 150 / 134 | 37 / 23 / 11 | 1.000; 0.98–1.16 / 0.98–1.10 / 0.94–1.15 | 5.1σ |

The strict rule discards 52 % of 24P's rows (a 50-pixel aperture almost always
contains one bad pixel) and costs four groups, three fits and eight robust CO
detections at unchanged Q.

### 7.3 Continuum flux space

`dc_main` versus `raw_main` (same policy, physical flux): Q ratio 0.998 (0.956–1.044)
for H₂O, 1.000 (0.990–1.003) for CO₂, 1.000 (0.976–1.078) for CO; nothing beyond
1.35σ; 444 of 462 band-rows keep their verdict.  The in-group spread of r_h² Δ²
has a median of 13 %, a 90th percentile of 58 % and a maximum of 166 %, and the size
of the Q change correlates with it (r = 0.23).  The two spaces are equivalent for
Q; the corrected space is kept because it removes a known gradient instead of
leaving it to the polynomial.

### 7.4 Error column

Empirical against formal errors (`_archive/results_comspec_2026-09-09/dc_main`):
Q ratio 0.999 / 1.003 / — for H₂O / CO₂, Q errors unchanged (the χ² rescaling had
been absorbing the difference), χ²_ν median 2.44 → 2.27.  The fits still run at
χ²_ν ≈ 2.3: the errors remain optimistic by ~1.5 in scatter, which the rescaling
carries into the quoted Q errors.

### 7.5 H₂O hot-band fallback

Before the rule, 16 fitted groups had CO detected with H₂O not in the model, so the
4.63 µm hot band was absorbed into Q(CO); now 1.  Robust H₂O rose from 28 to 40
(17 from hot bands: 24P phases 3–5 at 1.19–1.21 au, 306P phases 2–3, 2024 E1,
2025 L1, …), CO fell from 21 to 18 and the CO values shared by both runs dropped
by 9 % (0.73–1.03), the removed contamination.  The hot-band values are weaker than
the main-band ones (median 1.8σ against 2.9σ): ten of the 17 lie inside 3 au and two
of those reach 3σ (24P phases 4 and 5); the seven beyond 3 au (up to 29P at 6.3 au
and 2019 U5 at 9 au) are all 1–2σ with errors comparable to the value, and are
not to be quoted as water production rates without a distance cap or the 3σ tier.

### 7.6 Previous photometry

`apphot_comparison.md`: on identical channels the revised photometry changes
Q(H₂O) by 0.98, Q(CO₂) by 1.13 (1.00–1.22) and Q(CO) by 1.08; the previous set's
zeroed bad and star-adjacent pixels (29 % of emission-window channels) depressed
the 4.3 µm band of the brightest comets by 10–40 %.

### 7.7 Fluorescence database and the CO Swings factor (2026-09-11)

`dc_main_gauss` re-runs the main configuration with the previous emission model
(eight Gaussian bands, Ootsubo et al. 2012 g-factors, constant g(CO) = 2.6 × 10⁻⁴)
and reproduces the 2026-09-09 census exactly.  Against it, the main run with the
reconstructed database (`doc/comspec/fluorescence_database.md`; H₂O ν₃ 3.14 × 10⁻⁴, CO₂ ν₃
2.71 × 10⁻³, CO 1.92 × 10⁻⁴ s⁻¹ at v_h = 0 and 25–31 % more beyond 10 km/s) changes
the production rates of the groups detected in both by (median, 16–84 %):

| species | Q_new / Q_old | n | what drives it |
|---|---|---|---|
| H₂O, main band | 1.09 (0.89–1.16) | 38 | ν₃ 11 % stronger but the real band (2.55–2.9 µm) is 2× wider than the 0.05 µm Gaussians |
| H₂O, hot bands | 1.22 (0.90–2.23) | 24 | 4.63/4.85 µm hot bands 20 % weaker in total, 4.85 µm 40 % weaker |
| CO₂ | 1.27 (1.09–1.38) | 75 | ν₃ 7 % weaker, but the band is 0.05 µm wide against a 0.02 µm Gaussian |
| CO | 1.51 (1.22–2.27) | 20 | Swings factor (×1.35 at v_h ≈ 0, ×1.03 beyond 10 km/s), a resolved 0.13 µm-wide P/R band against a 0.045 µm Gaussian, and less H₂O hot-band flux subtracted |

Why the shapes matter this much: with dense channel sampling the integrated model flux
per unit Q follows the g-factors exactly (ratios 0.97, 1.00, 1.11 for CO₂, CO, H₂O),
but the *peak* flux per unit Q is 33 %, 47 % and 14 % lower for the real bands, and a
SPHEREx phase samples each band with a few channels, usually near the peak.  The
narrow placeholders therefore biased Q(CO₂) and Q(CO) low.  The database templates fit
the data better — χ²_ν median 1.91 against 2.27, lower in 95 of 138 groups; 45 of the
56 CO₂ groups above 3σ — and the fitted spectra now show the resolved P/R structure of
CO at 4.6–4.75 µm.  The census grows because the wider templates spread the
information over more channels (n_eff rises): robust H₂O 40 → 41 (main-band 23 → 29,
hot-band 17 → 12), CO₂ 23 → 34, CO 18 → 22; comets with a robust detection 41 → 52,
with a robust ≥ 3σ one 23 → 35.  The velocity of every pointing comes from the
ephemeris r_h(t) in the photometry and agrees with JPL Horizons to 0.05 km/s at
worst (118 checks on 10 comets); four of the 22 CO detections sit at |v_h| < 3 km/s
where the factor is < 1.15.

### 7.8 The 2026-09-12 revision (`936ed47a7a81`)

**What changed.** Photometry regenerated with the physical annulus and full-precision
times; Δ rule 1b; S/N-driven apertures; the 2.50 µm edge, wider continua, one-sided
extension and CV orders; GLS, 3σ/marginal tiers, no negative cut, the 3 au hot-band cap.
Groups 174 → 193, fits 138 → 151, apertures changed for 64 of 68 targets (60 000 km for
14 targets, 100 000 km for one; the PSF bound relaxed for 7 distant targets, three
crowded fields at their smallest bounded aperture).

**Annulus alone** (`results/comspec/studies/annulus_previous/`: the same pipeline on the archived
photometry, paired on the 133 groups whose aperture did not change): Q ratios of
1.000 for CO₂ and CO and 1.001 for H₂O overall — the ring only moves for comets inside
~1.6 au, and there the four groups inside 1 au gain 4 % in Q(H₂O).  The continuum
subtraction removes most of the coma flux the near ring had taken from both the band
and its continuum.

**Rule by rule** (each switched on the new photometry with the Δ grouping and the 3σ
tiers; the first-version aperture rule scored every non-`a` channel, the released rule
star-free channels only): median χ²_ν and ≥ 3σ robust counts (H₂O / CO₂ / CO):

| configuration | χ²_ν | ≥ 3σ, n_eff ≥ 2 |
|---|---|---|
| all previous rules, diagonal errors, 1σ cut | 2.19 | 21 / 26 / 13 |
| … without the 1σ negative cut | 3.16 | 19 / 26 / 12 |
| … with generalised least squares | 4.39 | 14 / 26 / 5 |
| + S/N-driven apertures (first version) | 5.88 | 20 / 29 / 6 |
| + new windows | 7.37 | 17 / 27 / 6 |
| + cross-validated orders | 9.88 | 20 / 30 / 6 |
| **released run** (apertures rescored on star-free channels) | **7.69** | **20 / 26 / 8** |

The χ²_ν increase is the sum of honest changes, not a bug: the channels the old cut
discarded, continuum errors the diagonal treatment had averaged down, and apertures and
windows that expose more of the systematics (star contamination in wide apertures, a
possible band-3/band-4 continuum step, the LSF/band-shape mismatch of the bright comets,
whose χ²_ν reach 10³–10⁴ with 2–4 channels on the CO₂ band).  GLS changes the values
little (`dc_main_diag` vs main: 0.96 / 1.00 / 0.99) but halves the ≥ 3σ CO count: the CO
hump is degenerate with a common-mode offset of a continuum interpolated across
4.55–4.90 µm, which the diagonal errors had priced as N independent channels.

**Against the 2026-09-11 run** (groups matched by target and time, ≥ 3σ in both):
Q(H₂O) 1.20 (1.10–1.40, n = 17), Q(CO₂) 0.98 (0.88–1.06, n = 43), Q(CO) 1.05
(0.72–1.14, n = 5).  Against the previous rules on the same photometry
(`dc_rules_previous`): H₂O 0.80, CO₂ 1.005, CO 0.84, i.e. the new apertures and
windows raise Q(H₂O) by 25 % and Q(CO) by 19 %.  Continuum verdicts: 2.7 µm 99 PASS /
14 WARN / 49 FAIL, 4.3 µm 110 / 8 / 44, 4.7 µm 94 / 24 / 44; 77 windows extended
(22 blue, 55 red), of which 9 reach PASS or WARN; the CV keeps 319 of 367 polynomial
continua linear.  Hot-band Q(H₂O): 22 rows, 2 detected and 3 marginal, all inside
2.84 au; 13 groups beyond 3 au carry the "fallback withheld" caveat.

### 7.9 The 2026-09-14 revision (`a6463b4c610b`): fixed apertures and the method matrix

**Apertures.**  One fixed aperture per phase (§3) in place of the S/N-driven per-target
choice: 20 000 / 40 000 / 60 000 km for 99 / 73 / 21 phases; 12 phases keep the aperture
the S/N rule had chosen, 103 are smaller, 78 larger.  Against `dc_ap_snr` (the S/N rule
on the otherwise identical pipeline) the values are unchanged — Q(H₂O) 0.99, Q(CO₂) 1.01,
Q(CO) 1.00 for the groups ≥ 3σ in both — and the census moves by a few groups: two CO
(2025 A6, 210P) and two CO₂ robust detections are lost, two H₂O and one CO₂ gained;
comets with a clean detection 33 → 34.  χ²_ν scales with the aperture in pixels: medians
10.4 at 20 000 km, 4.0 at 40 000, 2.3 at 60 000 — the small apertures of the bright, close
comets carry the LSF/band-shape systematic (concern 1).  The upstream photometry was not
rerun: the three radii are in every table since 2026-09-12.

**Method matrix** (`notebooks/comspec/method_matrix.py`; `results/comspec/studies/method_matrix/matrix.csv`,
54 runs of the main variant on the fixed apertures, each in an isolated tree).  Factors:
windows (the 2026-09-12 set; the 2026-09-11 set; eight intermediates over the 2.7 µm
emission edge 2.50 / 2.55 / 2.60 µm, the H₂O continuum 2.20–3.10 / 2.30–3.00 / 2.45–3.00 µm
and the CO₂ continuum 4.00–4.55 / 3.90–4.65 µm), polynomial order (CV at 5 / 10 / 20 %;
fixed 3/2/2, 2/2/2, 1/1/1), the one-sided extension (1 µm / off), the error model (GLS /
diagonal) and the accepted verdicts (PASS+WARN / PASS).  Score: *clean* robust detections
— ≥ 3σ, n_eff ≥ 2 and the carrier band's continuum verdict PASS — with the negative tail
of the fits (Q_fit/Q_err ≤ −2, which no coma can produce) as the false-positive gauge.

| factor | effect on the clean census (GLS) |
|---|---|
| one-sided extension | +0 to +2 clean and +7 fits in every pair; kept |
| order rule | CV 10 % = CV 5 % = CV 20 % = fixed linear within one count (44–45); fixed 3/2/2 and 2/2/2 lose 6–9; CV kept (78 % of the continua end up linear) |
| 2.7 µm emission edge | 2.55 µm (17 clean H₂O) ≥ 2.50 (16) ≥ 2.60 (15), each over 2.30–3.00 µm |
| H₂O continuum | 2.30–3.00 µm beats 2.20–3.10 by one to two clean H₂O in both pairs; 2.45–3.00 (band 4 only) halves the PASS count (56 vs 95): the blue side is too short |
| CO₂ continuum | 4.00–4.55, 3.90–4.55, 4.00–4.65, 3.90–4.65 µm: identical census |
| accepted verdicts | PASS only: the same clean count, 9 fewer fits and 2 fewer robust CO |
| error model | diagonal errors add 5–11 "robust" detections in every pair (H₂O +4 to +8, CO₂ +2 to +4, CO +1 to +2) |

**Why the diagonal model's extra detections are not taken.**  The negative tail
calibrates the errors: with the chosen windows 5 of 113 H₂O fits and 5 of 92 CO fits sit
at ≤ −2σ under GLS (Gaussian expectation 2.6 and 2.1), against 9 and 12 under diagonal
errors — three to six times the expectation.  The GLS error is 1.21× the diagonal one for
H₂O (IQR 1.10–1.39), 1.10× for CO₂ and 1.13× for CO, and every diagonal-only "detection"
sits at 0.8–2.8σ under GLS, i.e. in the marginal tier.  GLS stays (§6, concern 4), and its
census is the honest one.

**Result** (147 fits of 161 analysed groups; 32 skipped, 14 not fittable): robust ≥ 3σ
H₂O 17 (all clean, all from the main band), CO₂ 26 (25 clean), CO 4 (2 clean) → 35 comets
with a robust and 34 with a clean detection; marginal 25 / 23 / 9 groups → 19 further
comets.  χ²_ν median 4.6 (7.7 on 2026-09-12; 2.1 under the old rules).  Continuum
verdicts 2.7 µm 95 / 14 / 52, 4.3 µm 113 / 8 / 40, 4.7 µm 82 / 30 / 49 (PASS / WARN /
FAIL); 77 windows extended (20 blue, 57 red), 14 of them reaching PASS or WARN.  Hot-band
Q(H₂O): 23 rows, 1 detected and 6 marginal; 19 groups beyond 3 au carry the "fallback
withheld" caveat.  **Against the 2026-09-12 run** (same phases, ≥ 3σ in both): Q(H₂O) 1.01
(0.93–1.23, n = 16), Q(CO₂) 1.01 (0.97–1.09, n = 22), Q(CO) 1.00 (n = 4) — the values are
stable, the census tightened (H₂O 20 → 17, CO 8 → 4).  Against `dc_rules_previous` (old
aperture rule, orders, no extension, 1σ cut and tier, diagonal errors, no hot-band cap, on
the same photometry and windows): H₂O 1.19, CO₂ 0.97, CO 1.27; its counts at ≥ 3σ by
nsig are 22 / 32 / 10, with diagonal errors and the 1σ cut.  `dc_main_diag` vs main:
1.07 / 0.99 / 1.01.

## 8. Concerns

Resolved on 2026-09-12 (§7.8): **1** the sky annulus is physical (150 000 km, upstream);
**2** `jd_utc` is written at full precision upstream; **3** Δ is bounded by the grouping
(rule 1b); **4** the continuum covariance is propagated in full (GLS) — and on 2026-09-14
the negative tail of the fits showed the diagonal errors to be three to six times too
optimistic (§7.9); **5** `detected` means ≥ 3σ, with a `marginal` tier for 1–3σ; **6** the
hot-band fallback stops at 3 au; **8** no channel is dropped for being negative.  On
2026-09-14 the aperture returned to a fixed rule, one per phase (§3), and the windows,
order rule, extension and error model were chosen by the method matrix (§7.9).

Open:

1. **χ²_ν of the fits is ~4.6** (median; 10.4 at 20 000 km, 4.0 at 40 000, 2.3 at 60 000;
   2.1 under the old rules on the same photometry).  The bright comets at small apertures
   in pixels (χ²_ν of 10²–10³ with few channels on a narrow band) are limited by the
   instrument line-spread function and the band shape, not by noise (placeholder 1).  The
   errors are rescaled by √χ²_ν, so every tier count is conservative.
2. **Star contamination.**  Flag-`b` channels enter the fit under `BASELINE_FLAGS`; 120 of
   the 147 fits contain some (median 4 channels).  The §7.1 finding ("dropping `b` does not
   move Q") was made at the apertures of the time; `dc_main_no_b` is the re-test at
   20 000–60 000 km.
3. **The H₂O continuum window spans the band-3/band-4 boundary at 2.42 µm.**  A
   band-4-only blue side (2.45–2.60 µm) is too short to validate (PASS 56 instead of 95,
   §7.9), so the window keeps 2.30–3.00 µm; a calibration step between the detectors
   would bend the continuum under the 2.7 µm band.
4. **The one-sided extension saves few fits** (14 of 77 extended windows reach PASS/WARN);
   most one-sided cases lack emission channels altogether.
5. **CO is the weakest species**: 4 robust detections (2 clean) against 5 of 92 covered fits
   at ≤ −2σ; the CO hump is degenerate with a common-mode offset of the continuum
   interpolated across 4.55–4.90 µm, and the S/N-driven apertures had found 8 (§7.9).
6. **Conventions**: the aperture thresholds (3 au, 1.5 px, 60 000 km), the physical annulus
   radius (150 000 km), the Δ tolerance (20 %) and the CV margin (10 %; 5 and 20 % give the
   same census).
7. **Rule 3 of the grouping reads the raw `badphot`**, so the assignment does not move with
   the row policy (deliberate; 2025 K1 is the case to check).

## 9. Placeholders

Applied: 2 (fluorescence database and CO Swings factor, 2026-09-11), 7, 8, 9, 10, 16, 17
(the 2026-09-12 revision: CV orders, 3σ tiers, no negative cut, GLS, the physical annulus,
the Δ rule), 4 and 6 (windows and the aperture rule, re-decided on 2026-09-14 by the method
matrix: 2.55–2.80 µm over 2.30–3.00 µm; one fixed aperture per phase), 5 (`badphot` →
`frac_badpix_ap > 0.05`).
Open, in priority order (`results/comspec/placeholders.csv`, `config.PLACEHOLDERS`):

| priority | placeholder | current value | evidence | update |
|---|---|---|---|---|
| 1 | SPHEREx LSF | Gaussian, FWHM = `wlwidth` | χ²_ν of 10³–10⁴ for bright comets with 2–4 channels on the CO₂ band; the CO / H₂O-hot separability at 4.7 µm | as-built LVF R(λ) and LSF |
| 3 | expansion velocity | 0.8 r_h⁻⁰·⁵ km/s | ±18 % in Q per 20 % in v_g | cited species-specific law |
| 11–15 | grouping thresholds, sufficiency gates, opacity calibration, T_rot, upstream flag thresholds | notebook / `spherex_apphot` values | no new evidence | see registry |
| conv. | aperture thresholds (3 au, 1.5 px, 60 000 km), annulus radius 150 000 km, Δ tolerance 20 %, CV margin 10 % | 2026-09-12/14 conventions | §7.8, §7.9 (the CV margin is insensitive) | vary the aperture thresholds and the annulus on the bright comets |

## 10. Products

`results/comspec/gas_fit.csv` (one row per comet and phase; columns in `results/README.md`),
`results/comspec/gas_fit_lines/`, `results/comspec/continuum_summary.csv`, `skipped_groups.csv`,
`not_fitted.csv`, `apertures.csv`, `phase_map.csv`, `phase_cuts.csv`,
`placeholders.csv`, `run.meta.json`, `afrho_ztf.csv` (the ZTF dust context -- A(0°)fρ at
each phase's ⟨r_h⟩ from the `ztfcomet` (the ZTF stage) trends, attached to `gas_fit.csv` and
`phase_map.csv` as the `afrho_*` columns by `scripts/comspec/attach_afrho_ztf.py`);
`data/comspec/phase_assignment.csv`, `data/comspec/emission/`;
`fig/comspec/emission_model/`, `fig/comspec/cont_subtract/`, `fig/comspec/phase_group/`,
`fig/comspec/summary_Q_vs_rhel.png`, `fig/comspec/summary_mixing_ratios.png`.  Every study variant
has the same set under `results/comspec/studies/<variant>/`, `data/comspec/studies/`, `fig/comspec/studies/`;
the cross-variant tables and figures sit directly in `results/comspec/studies/` and
`fig/comspec/studies/`; `results/comspec/studies/method_matrix/matrix.csv` (+ `.jsonl`) is the 2026-09-14
method matrix of §7.9, regenerated by `notebooks/comspec/method_matrix.py`.  Everything older is in `_archive/`.
