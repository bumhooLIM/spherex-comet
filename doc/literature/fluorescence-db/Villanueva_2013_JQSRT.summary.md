# Villanueva, Magee-Sauer & Mumma (2013) — Modeling of nitrogen compounds in cometary atmospheres: fluorescence models of NH3, HCN, HNC and HC3N

*JQSRT 129, 158–168. doi:10.1016/j.jqsrt.2013.06.010*

**Role in this project:** full-cascade models for NH3, HCN and HNC (ab initio databases), a resonant band model for HC3N ν1, and quantitative statements about how much cascade changes band g-factors — the best guide to what a HITRAN-only reconstruction loses.

## 1. Databases (Sect. 2; Tables 1–2)
- NH3: HITRAN/GEISA 27 994 lines (main isotopologue), reduced to 23 596 lines between 3 305 levels after removing incomplete IDs, K > J, ℓ > 2 and mixed s/a labels; > 95 % of levels consistent to 0.3 cm⁻¹ (78 % to 0.001 cm⁻¹). Ab initio BYTe: 1 138 323 351 lines, 4 167 360 levels; 3 010 levels corrected with HITRAN/GEISA energies when |ΔE| < 1 cm⁻¹. Many HITRAN NH3 lines had **wrong Einstein A and statistical weights** (1 586 flagged) — use BYTe's self-consistent A values.
- HCN: HITRAN 2 955 lines (841 levels); GEISA 2011 79 957 lines (77 782 theoretical from Harris); HPT ab initio 34 433 190 lines, 18 723 HCN levels + 3 979 HNC levels (HNC ground state 5 186 cm⁻¹ above HCN; delocalized levels above 16 800 cm⁻¹). HPT energies corrected with HITRAN (834) and GEISA (636 HCN, 1 458 HNC). HNC: GEISA 5 619 lines.
- HC3N: absent from HITRAN, GEISA only below 760 cm⁻¹. Fundamentals ν1 3327 (CH stretch), ν2 2274, ν3 2079, ν4 862, ν5 663, ν6 499, ν7 222 cm⁻¹. ν1 band model: 397 lines in 3258–3378 cm⁻¹ from B_gs = 0.1574028, D_gs = 1.8165 × 10⁻⁸, H = 2.87 × 10⁻¹⁵; ν1: 3327.37085, B = 0.15149762, D = 1.8065 × 10⁻⁸ cm⁻¹; band intensity 249.4 cm⁻² atm⁻¹ (1 × 10⁻¹⁷ cm molecule⁻¹) at 296 K. Hot bands carry 67 % of the flux at 296 K but < 8 % below 100 K and nothing at 50 K. **g(ν1) = 9.88 × 10⁻⁴ s⁻¹ at 100 K, 1 au** (Crovisier 1987: 8.8 × 10⁻⁴).

## 2. Cascade results (Sect. 4.1; Figs. 6–7)
- Solar pump as in the ethane paper (Kurucz continuum × Hase line atlas, differential-rotation broadening); the published g-factors assume v_h = 0.
- HCN ν1 (3311 cm⁻¹) and HNC ν3 (2024 cm⁻¹): no significant difference from pure resonant fluorescence — in fact the ν1 emission of HCN *decreases* slightly because part of the ν1 population cascades to ν2. Hence cascade cannot explain the factor ~2 IR/radio disagreement in HCN production rates.
- NH3: many significant cascades. Total pump + cascade into the 1000 level is 4.08 × 10⁻⁵ s⁻¹, 18 % above the direct ν1 pump (3.45 × 10⁻⁵), but the level decays with branching 0.789 to ground (ν1) and 0.191 to 0001 (ν1−ν4), so the net ν1 g-factor is 7 % *lower* than the resonant value, not 17 % higher as Kawakita & Mumma (2011) estimated. ν3 is more complex still (numerous cascades to and from it). NH2 emission at 3 µm is removed with the Kawakita & Mumma model.
- Cascade maps (Fig. 6) shown for vibrational rates > 10⁻⁶ (NH3) and > 10⁻⁷ s⁻¹ (HCN, HNC) at 1 au, 100 K.

## 3. HCN ⇌ HNC radiative isomerization (Sect. 3)
Three mechanisms; the only trusted one (excitation to a labelled delocalized level) gives g(HCN→HNC) = 1.7 × 10⁻¹² and g(HNC→HCN) = 2.9 × 10⁻⁹ s⁻¹; over the 77 000 s lifetime only 0.1 ppm of HCN converts to HNC and 0.02 % of HNC to HCN — far too small to explain the HNC/HCN variation with R_h in Hale–Bopp.

## 4. Cometary application (Sect. 4.2; Fig. 8)
C/2007 W1 (Boattini), NIRSPEC KL2 order 25, 2008-07-09: T_rot(HCN) = 84 ± 5 K (H2O 80 ± 2 K in the same setting); NH3 detected at 3295.2 cm⁻¹ (aqP(2,1) + aqP(2,0) blend; sqP(1,0) at 3317.1 blended with HCN R1 at R = 25 000): NH3/H2O ≈ 0.7 %; HC3N < 0.017 % (3σ). Blended regions are fitted globally (all species simultaneously).

## 5. Implications for the SPHEREx reconstruction
- HCN from HITRAN (128 215 lines, ν1 3.0 µm, ν3 4.8 µm, ν2 overtone/combination pumps) is adequate: cascade corrections are at the percent level.
- NH3 from HITRAN (92 117 lines to 12 810 cm⁻¹) gives the pumps but only part of the cascade network; expect the ν1/ν3 band totals to differ by ~10 % from the BYTe-based values — acceptable at SPHEREx resolution, flag it.
- HC3N ν1 can be regenerated exactly from the constants above (Σ–Σ band, P/R branches only) and is a clean check of the pump normalization: 9.88 × 10⁻⁴ s⁻¹ at 100 K.
