# Cometary infrared fluorescence: concepts and the reconstruction of the GSFC line database

Scope: solar-pumped fluorescence g-factors of cometary parent molecules over 0.7–5.0 µm, for the
SPHEREx coma emission model. Sources: the five papers in `doc/literature/fluorescence-db/`
(Villanueva et al. 2011 Icarus, 2011 JGR, 2012 JQSRT, 2012 ApJ, 2013 JQSRT; summaries next to the
PDFs) and the PSG description of the GSFC Fluorescence Database (Villanueva et al. 2018, JQSRT 217,
86: 17 species, 530 281 non-LTE lines below 10 µm, derived from 2 billion ab initio/empirical lines).
Products: `data/fluorescence/` (formats in its README); code: `notebooks/comspec/fluorescence_gfm/`;
figure: `fig/comspec/fluorescence_profiles.png`.

## 1. Scientific concepts

**Why fluorescence.** Coma densities are too low for collisions to excite vibrations, and quenching
is far slower than radiative decay, so vibrational populations are set by the solar radiation field
alone (non-LTE). A molecule in a rotational level of the ground vibrational state absorbs a solar
photon into an excited vibrational level and re-emits within ~10 ms. Decay straight back to the
ground state is *resonant* fluorescence; decay through intermediate vibrational levels is
*non-resonant* fluorescence and produces the *hot bands* (e.g. H2O ν3−ν2 at 4.63 µm, ν1−ν2 at
4.85 µm, the L-band 2ν1−ν3 / ν1+ν3−ν1 / … bands at 2.8–3.0 µm). The emission rate per molecule of a
line or band is its **fluorescence efficiency** g [photons s⁻¹ molecule⁻¹], quoted at 1 au and
scaled by R_h⁻². The observed band flux is F = N_ap · g · hc/λ / (4πΔ²) with N_ap the number of
molecules in the aperture (`doc/comspec/model_concept.md`).

**The pump.** For a line l→u the absorption rate per molecule in level l is
(g_u/g_l) A_ul c² F_ν(ν) / (8π h ν³) — the Einstein-B form of Villanueva et al. (2011) Eqs. C3–C7 —
times the fractional population of l, w_l exp(−c2 E_l/T_rot)/Q(T_rot). Only ground-vibrational-state
levels are populated, so the line list must contain every band that starts from the ground state up
to the energy where the solar photon flux × band strength becomes negligible (for H2O 99 % of the
pumping goes through E_up < 7500 cm⁻¹, and levels above 12 000 cm⁻¹ contribute band rates < 10⁻⁷ s⁻¹).
T_rot enters only through the rotational populations; the *band* g-factors are nearly independent
of it, the *line* g-factors are not.

**The solar spectrum.** A 5778 K blackbody is adequate only near 3 µm. The real solar continuum
exceeds a blackbody by 25 % at 1.6 µm (H⁻ opacity minimum) and falls below it by 5 % at 4.3 µm and
7 % at 4.7 µm; Fraunhofer lines (CO fundamental at 4.4–5 µm, OH and CH at 3 µm, atomic lines below
1 µm) remove flux exactly where cometary lines of the same molecule absorb, unless the comet's
heliocentric velocity shifts them apart — the **Swings effect**. The GSFC recipe is an empirical
solar line list (Hase et al. 2006) multiplied by a theoretical continuum (Kurucz), broadened by the
solar rotation (2.5 km s⁻¹ FWHM) and evaluated at the Doppler-shifted frequency ν(1 + v_h/c).

**The cascade.** Each excited level u receives the direct pump plus the cascade from every higher
level that decays into it, and distributes the total among its decay lines with the branching ratios
A_i/A_tot(u) (Eq. C8). Processing levels in order of decreasing energy makes the network a single
pass. Branching ratios are temperature independent; pumps are not. Full cascade matters for H2O
(the 4.6–4.9 µm and L-band hot bands are fed by pumps at 1.38, 1.14 and 0.94 µm), for NH3 (many
cascades; the net ν1 g-factor is 7 % *below* the resonant value because 1000 also decays to 0001),
and hardly at all for HCN (cascade lowers ν1 slightly), CO, CO2, C2H6, CH3OH.

**Spectroscopic inputs.** Terrestrial databases (HITRAN, GEISA) are precise but complete only for
lines measurable at 296 K, hence they lack most hot-band decay routes of highly excited levels.
GSFC therefore takes Einstein A coefficients from self-consistent ab initio lists (BT2/VTT for
H2O/HDO, BYTe for NH3, HPT for HCN/HNC) and corrects their energies with semi-empirical levels
(SELP, HITRAN, GEISA, IUPAC) when the six quantum numbers match within 0.3 cm⁻¹ (1 cm⁻¹ for NH3,
HCN). Where no line list exists, an effective-Hamiltonian band model is built from laboratory
constants (C2H6 ν7: 17 266 lines; CH3OH ν3: 2 948 lines; HC3N ν1: 397 lines) with intensities from
laboratory band strengths and A coefficients from the Šimečková et al. (2006) conversion using the
*total* partition sum at 296 K (rotational × vibrational/torsional). Statistical weights carry the
nuclear-spin degeneracies, so ortho/para populations are in equilibrium unless a spin temperature is
imposed.

## 2. Procedure followed here

1. **Line lists.** HITRAN 2020 main isotopologue, 0–25 000 cm⁻¹, through `astroquery.hitran`
   (H2O 269 127 lines, CO2 173 129, CO 1 344, CH4 323 311, C2H6 63 516 including the ν7 and ν5
   bands, C2H2 81 454, C2H4 59 536, H2CO 346 546, NH3 92 117, HCN 128 215, OCS 530 134, H2S 90 507).
   CH3OH (HITRAN stops at 1407 cm⁻¹) and HC3N (absent) are regenerated as band models from the
   constants of the ApJ 2012 and JQSRT 2013 papers (`gfm_bandmodels.py`): symmetric-top energies
   with the torsional K ladders and Table 1 sub-band shifts for CH3OH ν3, Σ–Σ P/R structure for
   HC3N ν1; the CH3OH ν2 and ν9 bands use the Rogers (1980) intensity split quoted in the paper
   (28 % and 50 % of the C–H region, hot-band fraction 33 %) with approximate centres at 3000 and
   2970 cm⁻¹ and perpendicular-band Hönl–London factors — an envelope-level approximation.
2. **Level bookkeeping** (`gfm_core.load_hitran_csv`, `prepare`). HITRAN quantum-number strings are
   normalised per molecule class so that one level has one key (asymmetric tops → J Ka Kc; linear
   → J with e/f parity, Q branches flip parity; CH4 → J, symmetry, ranking index; NH3 → J K s/a).
   This matters: taken literally, the mixed-source spellings double-count the ground-state levels
   of H2O, HCN, CH4 and NH3 and halve every pump. Lines whose upper level is unlabelled are treated
   as isolated resonant lines. Ground-state labels are those whose lowest level lies below
   250 cm⁻¹ with at least 20 such lines (CH4 has two spellings, NH3 the s and a inversion
   components). Lines with unknown (negative) lower energy are dropped.
3. **Partition sums.** Explicit sums over the ground-state levels (Eq. C1), which reproduce the
   HITRAN TIPS values at 296 K to < 1 % for H2O, CO, HCN, C2H2, H2S, H2CO and to 1–2 % for CO2, OCS
   once the harmonic vibrational factor of the low-lying modes is included (Eq. C2 generalised).
   For CH4 and C2H4 the HITRAN labels still double count, so the sums are replaced by analytic
   rigid-rotor sums with the correct spin statistics (Td reduction for CH4; 586.7 vs TIPS 590.5).
4. **Solar pump** (`solar_pump.py`). Kurucz ATLAS9 solar model (`fsunallp`, R = 10 000 resampled),
   whose *line-free continuum column* is scaled to 1 au (× 4π(R☉/au)² = 2.720 × 10⁻⁴), multiplied
   by the disk-integrated Fraunhofer pseudo-transmittance of Toon (JPL, merged ATMOS/ACE/Kitt Peak
   atlas, 2024 release, 600–33 300 cm⁻¹ at 0.01 cm⁻¹). This is the Villanueva Appendix B recipe
   with a newer empirical line atlas. Check against the TSIS-1 Hybrid Solar Reference Spectrum:
   product/TSIS = 1.015 ± 0.004 over 0.7–2.7 µm in 100-nm bins; the Kurucz spectrum integrates to
   1368.6 W m⁻² (TSI 1361). No rescaling applied.
5. **Pumps and cascade** (`gfm_core.gfactors`): Eqs. C3–C8 with level-by-level cascade in
   descending energy order; solar flux read at ν(1 + v_h/c).
6. **Products** (`build_fluorescence_db.py`): line g-factors at T_rot = 30, 50, 70, 100, 130 K,
   band sums, spectral densities g_λ on an R = 5000 grid, Swings sensitivity at 70 K for
   v_h = −30 … +30 km s⁻¹, the pumping spectrum, and the validation table.

## 3. Validation (`data/fluorescence/validation.csv`)

Reconstructed / published g-factors (line level unless noted), at the published T_rot and, for the
Boattini values, v_h = +9.8 km s⁻¹:

| species | published values | median ratio | scatter |
|---|---|---|---|
| CO 1–0, 5 lines (Icarus 2011) | 12.6–18.4 × 10⁻⁶ | 1.00 | ± 0.01 |
| H2O, 16 hot-band lines at 2.9–3.0 and 4.65 µm (Icarus 2011) | 0.2–5.3 × 10⁻⁷ | 1.01 | ± 0.10 |
| CH4 ν3, 6 lines | 4.5–20.6 × 10⁻⁶ | 1.01 | ± 0.01 |
| HCN ν1, 12 lines | 8–26 × 10⁻⁶ | 0.99 | ± 0.01 |
| C2H2 ν3 / ν2+ν4+ν5, 11 lines | 5–19 × 10⁻⁶ | 0.95 / 1.02 | ± 0.01 |
| C2H6 ν7, 7 Q-branch windows | 16–40 × 10⁻⁶ | 1.07 | ± 0.07 |
| CH3OH ν3 band (ApJ 2012) at 75 / 150 K; Q branch at 80 K | 143.4, 128.4, 15.5 × 10⁻⁶ | 1.03, 1.08, 0.94 | |
| HC3N ν1 band at 100 K (JQSRT 2013) | 9.88 × 10⁻⁴ | 1.01 | |
| NH3 ν1, 3 lines | 0.6–1.5 × 10⁻⁶ | 0.81 | ± 0.02 |

The CO lines also confirm the Swings treatment: at v_h = 0 they would be 16–22 % weaker than at
+9.8 km s⁻¹, exactly the effect the papers describe. The 10 % scatter of the H2O ortho/para lines
around unity (mean 0.995) shows that the HITRAN cascade routes of the 200, 101 and 001 levels are adequate.
NH3 is 20 % low because HITRAN lacks most of the BYTe cascade network (and carries wrong weights
for 1 586 lines, as the 2013 paper notes); use NH3 with that caveat.

## 4. Results relevant to SPHEREx (1 au, 70 K, v_h = 0)

| band | g [s⁻¹] | reference value | note |
|---|---|---|---|
| H2O ν3 2.68 µm | 3.14 × 10⁻⁴ | 2.82 × 10⁻⁴ (Ootsubo et al. 2012 / Crovisier) | continuum +5 % over blackbody, cascade into 001 |
| H2O ν1 2.78 µm | 2.9 × 10⁻⁵ | 2.47 × 10⁻⁵ | |
| H2O ν2+ν3−ν2 2.70 µm; ν1+ν3−ν1 2.81 µm | 2.85, 1.94 × 10⁻⁵ | 2.81, 2.19 × 10⁻⁵ | under the 2.7 µm band |
| H2O 2ν1−ν3, 2ν1−ν1 and other L-band hot bands 2.8–3.0 µm | ≈ 1 × 10⁻⁵ total | — | the red wing of the 2.7 µm feature |
| H2O ν3−ν2 4.67 µm; ν1−ν2 4.85 µm; ν2+ν3−2ν2 4.63 µm | 6.6, 4.4, 1.1 × 10⁻⁶ | 7.66, 7.35 × 10⁻⁶ | the CO-region blend |
| H2O ν1+ν3 1.38 µm; ν2+ν3 1.88 µm | 6.5, 5.6 × 10⁻⁶ | — | 2 % of ν3 each; SPHEREx bands 2–3 |
| CO2 ν3 4.26 µm | 2.71 × 10⁻³ | 2.86–2.9 × 10⁻³ | continuum 5 % below blackbody |
| CO2 ν3 hot bands 4.30 µm | 1.0 × 10⁻⁴ | — | 3.8 % of ν3, unresolved |
| CO 1–0 4.67 µm | 1.92 × 10⁻⁴ at v_h = 0; 2.4–2.5 × 10⁻⁴ for |v_h| ≥ 10 km s⁻¹ | 2.6 × 10⁻⁴ (blackbody) | **the only band with a large Swings effect (−25 % at v_h = 0)** |
| CH4 ν3 3.30 µm | 3.5 × 10⁻⁴ | | plus ν1+ν4 at 2.34 µm, 1.7 × 10⁻⁵ |
| C2H6 ν7 3.35 µm; ν5 3.45 µm | 4.6, 1.2 × 10⁻⁴ | | ν7 5.0 × 10⁻⁴ at 30 K |
| CH3OH ν9 3.37; ν2 3.34 (approx.); ν3 3.52 µm | 3.5, 2.0, 1.5 × 10⁻⁴ | ν3 1.43 × 10⁻⁴ | |
| C2H2 ν2+ν4+ν5 3.05; ν3 3.04 µm | 2.0, 1.8 × 10⁻⁴ | | |
| H2CO ν5 3.53; ν1 3.60 µm | 4.0, 3.2 × 10⁻⁴ | | |
| HCN ν1 3.02 µm; HC3N ν1 3.01 µm | 3.5 × 10⁻⁴; 1.07 × 10⁻³ | | |
| NH3 ν1 3.0 µm and other bands 2.0–3.1 µm | ≈ 1.3 × 10⁻⁵ each, 1.3 × 10⁻⁴ total | | −20 % |
| OCS ν3 4.85 µm | 2.8 × 10⁻³ | | sits on the H2O ν1−ν2 / CO blend |
| H2S 2.65 µm | 1 × 10⁻⁵ | | negligible |

Use in the pipeline (wired in on 2026-09-11, `spherex_comspec/fluorescence.py`):
the species templates of the fit are the `profiles/` spectral densities at T_rot, and g(CO) is scaled
with the heliocentric velocity of every channel through `co_swings.csv` (0.25 km/s grid, ±60 km/s,
30/70/130 K). Relative to the previous eight-band model the H2O ν3 and 2.7 µm hot-band strengths
differ by 10–15 %, the 4.85 µm hot band is 40 % weaker, and g(CO) at v_h ≈ 0 is 25 % lower than
the constant 2.6 × 10⁻⁴ used before (comets near perihelion). The `dc_main_gauss` study variant keeps
the old model for comparison. OCS remains a possible contaminant of the 4.85 µm region for
sulphur-rich comets.

## 5. Limitations

- HITRAN hot-band completeness: decay routes of levels above ~8000 cm⁻¹ are truncated, so weak
  hot bands fed by short-wavelength pumps are lower limits (the validated 2.9–3.0 and 4.65 µm
  H2O lines are not affected at the 10 % level).
- Band models: CH3OH ν2/ν9 centres and rotational structure are approximate (envelope accuracy
  ~30 cm⁻¹, fine at R ≤ 130); no torsional hot bands; HC3N follows the paper's intensity
  convention (249.4 cm⁻² atm⁻¹ taken as the fundamental at 296 K).
- C2H6: the explicit ground-state sum exceeds the TIPS rotational sum by 33 % (mixed statistical
  weights in HITRAN), yet the line-level validation is fine (1.07 ± 0.07); C2H4 uses the analytic
  partition sum and is unvalidated.
- The solar pump is sampled at the line centre (no Doppler broadening by the coma expansion), and
  the Toon atlas is a pseudo-transmittance whose weak-line depths carry a few per cent uncertainty;
  the continuum is a model (no rescaling to TSIS was applied; the offset is +1.5 %).
- Optically thin everywhere; opacity is handled by the pipeline (`gasmodel.opacity_factor`).
- Not included: HDO, CH3D, HNC (in the GSFC list), NH2, OH prompt emission, and any electronic
  transitions below 0.7 µm.
