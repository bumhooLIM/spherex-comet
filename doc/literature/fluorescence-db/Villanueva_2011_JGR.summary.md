# Villanueva, Mumma & Magee-Sauer (2011) — Ethane in planetary and cometary atmospheres: transmittance and fluorescence models of the ν7 band at 3.3 µm

*J. Geophys. Res. 116, E08012. doi:10.1029/2010JE003794*

**Role in this project:** the paper that defines the **General Fluorescence Model (GFM)** (Appendix C) and the **synthetic solar pumping spectrum** (Appendix B) used by every later GSFC fluorescence list, plus a complete line-by-line model of the strongest cometary ethane band. It is the recipe we reproduce for CO, CO2, CH4, C2H2, … from HITRAN.

## 1. The ν7 band model of C2H6 (Sect. 2)
- HITRAN (2008) had only some Q branches of ν7 (Brown et al. 1987; Pine & Rinsland 1999 for PQ3); no P/R lines. This work builds the full band: **17 266 lines** (8 680 in ν7, 8 586 in the torsional hot band ν7+ν4−ν4), K_max = 20, J_max = 40, 2900–3100 cm⁻¹, HITRAN-2008 format, G36⁺ point-group labels (vibrational tags V7, GROUND, V7+V4, V4; local quanta J, K, ℓ). Lines are listed per spin symmetry even when unresolved (e.g. A34s–A12s reported as A3s–A1s and A4s–A2s).
- Perpendicular band from C–H stretching; severely perturbed by overtones/combinations in Fermi and Coriolis resonance, and by the low-frequency torsion ν4 (289 cm⁻¹). The rotational structure is fitted **per K ladder** (individual constants for K″ ≤ 8 from 654 assigned lines, σ = 0.005 cm⁻¹; mean progressions otherwise). Reference Q-branch origins: RQ0 2986.72, RQ1 2990.08, RQ2 2993.46, RQ4 3000.21, PQ1 2983.38, PQ2 2980.08, PQ3 2976.79 cm⁻¹.
- Spin statistics: three spin modifications (A, E, G types with weights in Table 1c); E/A ratios in comets consistent with equilibrium (T_spin > 10 K). Torsional hot band: population of the first torsional level is ¼ of the ground state at 296 K; hot-band intensity 25 % of the fundamental at 296 K (Pine & Rinsland found ~20 % for the PQ3 replica); modelled only approximately (rotational constants of ν4 from Blass et al. 1990) and provided on request. In comets (T < 100 K) the hot band is not detected.
- Intensities: band intensity re-derived from PNNL absorption spectra; line intensities → Einstein A via Šimečková et al. (2006) with TIPS Q_tot(296 K) = 70 881 (Q_vib = 1.3732) and isotopic abundance I_a = 0.97699.

## 2. Solar spectrum for pumping (Appendix B) — reproduce this
- Flux-calibrated high-resolution solar spectrum = **empirical solar line list** (Hase et al. 2006: centre, strength, width, Gaussian/Lorentz shape per line, with centre-to-limb variation; validated against the ACE-FTS atlas of Hase et al. 2010 and ATMOS) × **theoretical continuum** (Kurucz irradiance model), 700–6400 cm⁻¹.
- Disk integration: convolve the disk-centre synthesis with a solar rotation kernel (differential rotation A = 14.713, B = −2.396, C = −1.787 deg/day; limb darkening u = 0.6; ~2.5 km s⁻¹ FWHM at the ecliptic). Called "far from optimum" but adequate until a disk-integrated line list exists.
- Why: a blackbody continuum is only right at 2900–3300 cm⁻¹; using it would produce g-factor errors up to 30 % for ν7 lines and up to 40 % for CO 1–0 lines, through the Swings effect (Fraunhofer lines Doppler-shifted by the comet's heliocentric velocity). Therefore g-factors must be recomputed for each (T_rot, v_h).

## 3. General Fluorescence Model (Appendix C) — the reconstruction recipe
Inputs per line (11 HITRAN parameters): molecule, isotopologue, ν (cm⁻¹), A21 (s⁻¹), E″ (cm⁻¹), global vibrational quanta V′, V″, local rotational quanta L′, L″, statistical weights w′, w″.
1. **Energy levels** rebuilt from the line list: E″(L″) = E_i″, E′(L′) = E_i″ + ν_i; L = J for linear molecules, (J, K, ℓ) for symmetric tops, (J, Ka, Kc) for asymmetric tops.
2. **Partition function**: HITRAN TIPS above 70 K; below 70 K an explicit sum over the rotational levels, Q_r = Σ w_i″ exp(−c2 E_i″/T) (C1), with Q_vib ≈ 1 except C2H6 (torsion at 289 cm⁻¹: Q_vib = [1 − exp(−c2 E_ν4/T)]⁻¹, C2).
3. **Pumping rate** into upper level L′ (C3–C7): g_pump(L′) = Σ J_s(ν_s) B12 w″ exp(−c2E″/T)/Q_tot(T), with B12 = B21 w′/w″, B21 = A21/(8πhν³), ν_s = (1 − u/c) ν the Doppler-shifted frequency at which the solar flux J_s is read (u = heliocentric radial velocity), and A_tot(L′) = Σ A21 over all lines leaving L′ (C4).
4. **Emission**: g_i = g_pump(L′) A21,i/A_tot(L′) (C8). Branching ratios are temperature independent; pumps are not, so g-factors are recomputed at each T_rot. Selection rules, Hönl–London and Herman–Wallis factors are implicit in the line list.
This resonant version was applied here to C2H6 ν7 (and to ν5, HCN ν1, CH4 ν3, CO/¹³CO 1–0, HDO ν1 and 2ν2, C2H2 ν3 and ν2+ν4+ν5 in follow-up work); the full-cascade extension is in the water paper.

## 4. Cometary validation (Sect. 4, 5.3; Tables 3–4)
- Three comets: 8P/Tuttle (CRIRES, 2008-01-26, T_rot = 60 K, v_h = −0.37 km s⁻¹), C/2004 Q2 Machholz (NIRSPEC, 2005-01-19, 93 K, v_h = −2.01), C/2007 W1 Boattini (NIRSPEC, 2008-07-10, 79 K, v_h = +9.70).
- Integrated Q-branch g-factors (10⁻⁵ s⁻¹ at 1 au, including underlying P/R lines) at 79 K: RQ4 1.585, RQ2 2.560, RQ1 3.071, RQ0 3.974, PQ1 3.067, PQ2 2.581, PQ3 2.196; at 93 K RQ0 = 3.654; at 60 K RQ0 = 4.484.
- Production rates: W1 (2.35 ± 0.02) × 10²⁶ s⁻¹ (1.957 ± 0.053 % of H2O); Q2 (13.3 ± 0.25) × 10²⁶ at 93 K (0.488 ± 0.016 %, −12 % vs earlier; consistent with the ν5-based 0.48 ± 0.06 %); 8P (1.74 ± 0.06) × 10²⁶ (0.291 ± 0.017 %). Confidence limits improve by a factor of ~3 over the Dello Russo et al. (2001) model.
- The old model failed because of (1) wrong ℓ-split symmetries → wrong branching ratios, (2) outdated constants and a T-independent "g-band" treatment, (3) missing P/R structure under the Q branches, (4) a blackbody pump.

## 5. Terrestrial/Mars part (Sect. 3, 5.1) — context only
LBLRTM transmittances against Mars with CSHELL/CRIRES/NIRSPEC; Mars ethane upper limits < 0.8 and < 1.5 ppb (3σ); telluric ethane 0.97 ppb (Mauna Kea) vs 0.30 ppb (Paranal), reproducing the hemispheric asymmetry. Appendix A lists HITRAN corrections (new CO2 isotopic bands at 3.3–3.7 µm).

## 6. Implications for the SPHEREx reconstruction
- Implement (C1)–(C8) exactly; they are what PSG's fluorescence lists are built from.
- HITRAN 2020 now contains the ν7 band (4 679 lines, V7–GROUND) and ν5 (4 596 lines) of C2H6, so the 3.3–3.5 µm ethane emission can be reconstructed directly; the torsional hot band is absent (negligible below 100 K).
- Expect resonance-only treatment for C2H6 (no cascade data), consistent with the original.
