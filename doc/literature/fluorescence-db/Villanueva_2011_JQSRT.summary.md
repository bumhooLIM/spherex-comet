# Villanueva et al. (2012) — Water in planetary and cometary atmospheres: H2O/HDO transmittance and fluorescence models

*JQSRT 113, 202–220 (received 2011 July; the file is labelled 2011). Authors: G. L. Villanueva, M. J. Mumma, B. P. Bonev, R. E. Novak, R. J. Barber, M. A. DiSanti. doi:10.1016/j.jqsrt.2011.11.001*

**Role in this project:** the reference description of the GSFC *full non-resonance cascade* fluorescence model for H2O (the species that dominates SPHEREx comet spectra at 2.7 µm and whose hot bands contaminate the CO region at 4.6–4.9 µm). The generic machinery (pumps, branching ratios, solar spectrum) is Appendix C of the ethane paper (`Villanueva_2011_JGR.summary.md`); this paper adds the *cascade* and the database work.

## 1. Physics of cometary water emission (Sect. 1.1, 3)
- Coma collisions cannot excite vibrations and quenching is slower than radiative decay, so the vibrational manifold is far from LTE. Solar photons pump ground-state molecules into excited vibrational states; decay is radiative, either back to the ground state (**resonant fluorescence**) or through intermediate vibrational levels (**non-resonant fluorescence**, which produces the "hot bands"). All three modes of H2O are IR active.
- With a 5778 K pumping source, 99 % of the pumping of a 100 K coma goes through transitions with E_up < 7500 cm⁻¹; pumps and cascades from levels above 12 000 cm⁻¹ have band rates < 10⁻⁷ s⁻¹ and are negligible. Cutting the line list at that energy reduces BT2 from 505.8 million to 1.3 million lines (factor 380) and VTT (HDO) to ~8 million (factor 90) with no significant change in the g-factors.
- Hot bands matter observationally because their lower levels are unpopulated in Earth's atmosphere, so they are transmitted from the ground. Bands used historically: ν3 fundamental (2.7 µm, KAO/Halley), 011–010 (ν2+ν3−ν2), 100–010 and 001–010 (4.85 and 4.63 µm), 111–100 (~2 µm), 021–010, and a network of ≥ 10 hot bands in the L band (2.9–3.0 µm: 101–001, 101–100, 200–001, 200–100, 201–200, 300–101, …) that NIRSPEC/CRIRES retrievals rely on. Photons in the fundamental bands are absorbed by telluric water.
- Earlier models (Bockelée-Morvan & Crovisier) assumed uncoupled modes and the harmonic-oscillator approximation to scale hot-band strengths from cold bands → imprecise branching ratios, rates and positions. Dello Russo et al. (2004/2005) added BT2 rotational branching ratios for four bands (200–100, 200–001, 101–100, 101–001) but kept Born–Oppenheimer vibrational branching and a blackbody pump.
- **Solar pump:** a blackbody is adequate only for the continuum at 2900–3300 cm⁻¹. Elsewhere the continuum deviates and Fraunhofer lines shift in and out of pump transitions with the comet's heliocentric velocity (**Swings effect**); ignoring it biases the *relative* line intensities, hence T_rot and spin temperatures. They use an empirical solar line list (Hase et al. 2006/2010) scaled by a Kurucz continuum model (Appendix B of the ethane paper).

## 2. Spectral database construction (Sect. 2)
- Line counts: GEISA 41 447, HITRAN 37 432, BT2 (ab initio) 505 806 255, HITEMP 114 209 395 lines for H2O; HDO: GEISA 11 980, HITRAN 13 238, VTT 697 450 825. Terrestrial databases are complete only for lines strong at 296 K and become severely incomplete above ~4000 cm⁻¹ of level energy — inadequate for cascade calculations.
- **Rule:** adopt the self-consistent **Einstein A coefficients of BT2/VTT** (branching ratios must be internally consistent) but **replace the ab initio energies** by semi-empirical values (SELP → HITEMP, HITRAN, GEISA; IUPAC for HDO) when the six quantum numbers (v1 v2 v3 J Ka Kc) match and |ΔE| < 0.3 cm⁻¹. Corrected: 11 989 BT2 levels with SELP, 7 486 with HITRAN, 7 044 with GEISA; 5 287 VTT levels with IUPAC. Even a single known (usually lower) level is substituted, which narrows the error spread of the whole set.
- Consistency checks: ortho/para symmetry from the parity of Ka+Kc+v3 (84 699 BT2 lines flagged); "negative rotational energy" flags for mislabelled quanta; energy consistency of terrestrial lists (99 % of HITRAN levels consistent to 0.3 cm⁻¹, 86.5 % to 0.001 cm⁻¹). BT2 leaves 97 % of lines unlabelled (conservative), VTT 21 % unlabelled but with an 11 % error rate.
- Corrections to BT2 frequencies reach 0.3 cm⁻¹ = 7 NIRSPEC pixels at 0.04 cm⁻¹/pixel — essential for line identification.

## 3. Fluorescence computation (Sect. 3) — the four steps
1. Update the BT2/VTT energy tables with SELP/IUPAC values.
2. Synthesize the high-resolution solar spectrum at the comet's heliocentric velocity.
3. Compute the solar pumps for all 500 million H2O lines (700 million HDO).
4. Compute g-factors for 1200 million lines with **line-by-line and level-by-level branching ratios** — each pumped level decays with branching A_i/A_tot; the receiving level (if excited) decays in turn, until the ground state is reached (Figs. 3–4: cascade diagrams at 1 au, 100 K; only vibrational rates > 10⁻⁷ s⁻¹ shown).
Parallelized; a full set of cometary g-factors takes hours. Line lists and energy tables were posted at astrobiology.gsfc.nasa.gov/Villanueva/spec.html (now the GSFC Fluorescence Database served by PSG).

## 4. Validation on comet C/2007 W1 (Boattini) (Sect. 5.1)
- NIRSPEC/Keck II, 2008 July 9–10, KL1/KL2/MWA settings; R_h = 0.893–0.899 au, v_h = +9.77 to +10.36 km s⁻¹, Δ = 0.35 au. Nucleus-centred 9-row extract (±0.9″).
- 26 water lines (10 para, 16 ortho): T_rot = 80 ± 2 K, OPR = 2.87 ± 0.15 (T_spin > 34 K); D/H < 8.3 × VSMOW (3σ).
- New vs old model: Q(H2O) = 116.25 ± 9.47 vs 133.45 ± 20.15 (10²⁶ s⁻¹): −15 % in value, factor 2 better precision. Lines from the 200–221 level (3394.1 and 3372.8 cm⁻¹), previously predicted 1.5–2× too weak and excluded, are now reproduced.
- Line-by-line g-factors at 80 K for the L-band hot bands and for the 4.65 µm 001–010 band are tabulated in the Boattini paper (`Villanueva_2011_Icarus.summary.md`): typical values 10⁻⁸–10⁻⁷ s⁻¹ per line (L band) and up to 5 × 10⁻⁷ s⁻¹ per line at 2137–2151 cm⁻¹.

## 5. LTE part (Sect. 4, 5.2) — not used here
Integration of the corrected H2O/HDO/CO2/C2H6 lists into LBLRTM; Mars CO2-broadening via the complex Robert–Bonamy formalism; retrievals of Mars and Mauna Kea water columns and D/H (Table 4). Irrelevant to SPHEREx except as a reminder that the same list serves telluric modelling.

## 6. Implications for the SPHEREx reconstruction
- H2O g-factors require cascade: the 4.63/4.85 µm hot bands and the 2.9–3.0 µm hot bands are fed by pumps at 1.38, 1.14, 0.94 µm etc., not by their own (unpopulated) lower levels.
- A HITRAN-only reconstruction reproduces the strongest pumps (all ground-state lines up to 12 000 cm⁻¹ are present) but truncates the decay routes of high levels, because HITRAN lists hot-band lines only if they are measurable at 296 K. Expect the reconstructed hot-band rates to be lower limits; document the deficit against the tabulated Boattini values.
- Use a solar spectrum with Fraunhofer lines and a real continuum; a blackbody is acceptable only around 3 µm.
