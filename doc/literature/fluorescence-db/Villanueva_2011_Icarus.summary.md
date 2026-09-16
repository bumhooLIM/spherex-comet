# Villanueva et al. (2011) — The molecular composition of Comet C/2007 W1 (Boattini): evidence of a peculiar outgassing and a rich chemistry

*Icarus 216, 227–240. Authors: G. L. Villanueva, M. J. Mumma, M. A. DiSanti, B. P. Bonev, E. L. Gibb, K. Magee-Sauer, G. A. Blake, C. Salyk. doi:10.1016/j.icarus.2011.08.024*

**Role in this project:** the application paper. Table 2 tabulates **line-by-line g-factors at 80 K** for H2O (L band and 4.65 µm), C2H6, HCN, CH3OH, CO, CH4, C2H2 and NH3 computed with the new GSFC models — the numbers we use to validate the reconstruction. It also documents the retrieval chain (line flux → Q) that our low-resolution band fitting mirrors.

## 1. Observations
- NIRSPEC/Keck II, 2008 July 9 and 10 (R ≈ 25 000), settings KL1 (HDO, CH3OH, C2H6, H2O), KL2 (HCN, C2H2, NH3, CH4, C2H6, H2CO, H2O; order 25 is a dense blend) and MWA (CO 1–0, 5 lines; H2O 3 lines near 4.7 µm). R_h ≈ 0.90 au, v_h ≈ +9.8 to +10.4 km s⁻¹, Δ ≈ 0.35 au.
- Nucleus-centred extracts: 9 rows = 0.432″ × 1.782″ ≈ 110 × 460 km²; also spatial profiles along the Sun–comet line.

## 2. Retrieval method (Sect. 2.2–2.5)
- T_rot from excitation analysis (F/g diagrams: an optimum temperature flattens Q retrieved from lines of different E_rot) for H2O (26 lines: 80 ± 2 K, OPR 2.87 ± 0.15), HCN (84 ± 5 K), C2H6 (79 ± 3 K); other species assume 80 K. Spectrally blended regions are fitted with a global Levenberg–Marquardt model of all species.
- Nucleus-centred production rates Q_nc from line fluxes with the new g-factors, outflow velocity v_0 = 0.8 R_h^−0.5 km s⁻¹; growth factors from the spatial profiles correct slit losses.
- g-factors: C2H6 (Villanueva et al. 2011 JGR), H2O/HDO (2012 JQSRT, full non-resonant cascade of 500 million lines, BT2/GEISA/HITRAN/HITEMP), CH3OH (2012 ApJ), CH4/C2H2/CO (GFM with revised HITRAN parameters), HCN/NH3 (2013 JQSRT; cascade; three NH3 line atlases TROVE/GEISA/HITRAN), H2CO (DiSanti et al. 2006), NH2 (Kawakita & Mumma 2011). All use the realistic solar pumping spectrum (Kurucz continuum × solar line list) instead of a blackbody.
- Remaining systematics: telluric transmittance for Doppler-broadened cometary lines, Fermi/Coriolis perturbations in ab initio intensities, instrumental artefacts; hot-band intensities from ab initio models need laboratory checks (H2O, NH3 in particular).

## 3. Results (Tables 3–5)
| species | T_rot (K) | Q (10²⁶ s⁻¹) | X/H2O (%) |
|---|---|---|---|
| H2O (Jul 9) | 80 ± 2 | 120.32 ± 8.90 | 100 |
| CH4 | (80) | 1.89 ± 0.23 | 1.57 ± 0.16 |
| CH3OH | (80) | 4.42 ± 0.34 | 3.67 ± 0.11 |
| C2H6 | 79 ± 3 | 2.36 ± 0.15 | 1.96 ± 0.04 |
| H2CO | (80) | < 0.14 (3σ) | < 0.12 |
| HCN | 84 ± 5 | 0.61 ± 0.04 | 0.50 ± 0.01 |
| C2H2 | (80) | 0.34 ± 0.03 | 0.29 ± 0.02 |
| NH3 | (80) | 2.09 ± 0.24 | 1.74 ± 0.17 |
| H2O (Jul 10) | 79 ± 3 | 122.34 ± 9.11 | 100 |
| CO | (80) | 5.52 ± 0.71 | 4.50 ± 0.51 |
Mixing ratios drop by ~2.5× (e.g. CO 1.8 %, C2H6 0.78 %, HCN 0.20 %) if the extended water source (source II) is included. Cosmogonic indicators: D/H < 12.9 × 10⁻⁴ (8.3 VSMOW, 3σ; T_form > 14 K), OPR(H2O) 2.87 ± 0.15 (T_spin > 34 K), CH4 F/A = 1.57 ± 0.48 (T_spin > 19 K).
- Spatial profiles: polar species (H2O, CH3OH) extended anti-sunward, apolar (C2H6, HCN) symmetric → two ice moieties: mixed ice/dust clumps released sunward, and tiny pure polar-ice grains (source II). Abundance ratios among trace species are more robust than ratios to water.

## 4. Table 2 — g-factors used for validation (1 au, photons s⁻¹ molecule⁻¹, comet velocity v_h ≈ +9.8 km s⁻¹)
- H2O ortho, 80 K (10⁻⁹ s⁻¹): (101)000–(001)111 3458.12: 73.63; (200)110–(001)111 3450.29: 110.23; (200)212–(001)313 3382.10: 164.44; (200)221–(001)322 3372.75: 90.11; (200)303–(001)404 3358.92: 80.49; (201)000–(200)101 3388.77: 29.79; (200)321–(100)432 3374.40: 16.94; (001)111–(010)110 2151.20: 451.12; (100)221–(010)110 2148.19: 117.99; (001)220–(010)221 2144.81: 200.34; (001)000–(010)101 2137.38: 530.28.
- H2O para, 80 K (10⁻⁹): (200)220–(001)221 3445.89: 34.65; (101)212–(001)303 3434.40: 21.87; (200)202–(001)303 3378.48: 46.57; (200)211–(001)312 3366.35: 34.41; (200)313–(001)414 3360.99: 27.10; (001)212–(010)211 2139.93: 47.21.
- C2H6 ν7 Q-branch windows, 79 K (10⁻⁶): RQ4 15.85, RQ2 25.60, RQ1 30.71, RQ0 39.74, PQ1 30.67, PQ2 25.81, PQ3 21.96 (include underlying P/R lines).
- HCN ν1, 84 K (10⁻⁶): R4 20.35, R2 19.05, R1 14.71, R0 8.13, P2 16.29, P3 22.12, P4 25.49, P5 26.17, P6 24.54, P7 21.22, P8 17.21, P9 13.09.
- CH3OH ν3 Q-branch (2844.3 cm⁻¹), 80 K: 15.5 × 10⁻⁶.
- CO 1–0, 80 K (10⁻⁶): R3 15.45, R2 16.72, R1 12.58, P2 14.91, P3 18.36.
- CH4 ν3, 80 K (10⁻⁶): R2(E) 11.26, R2(F) 16.95, R1 17.21, R0 20.57, P2(E) 4.46, P2(F) 6.69.
- C2H2, 80 K (10⁻⁶): ν3 R9e 5.80, R4e 5.07, R3e 14.98, R1e 9.92, P3e 14.97, P5e 18.89; ν2+ν4+ν5 R9e 6.55, R7e 11.29, R5e 15.51, R3e 16.05, R1e 10.70.
- NH3 ν1, 80 K (10⁻⁶): sqP(1,0) 3317.21: 1.02; aqP(2,1) 3295.43: 0.56; aqP(2,0) 3295.38: 1.52 (Kawakita & Mumma 2011 values 1.25, 0.63, 1.91).

## 5. Implications for the SPHEREx reconstruction
- Line-level agreement with these values (within the Swings-effect uncertainty at v_h ≈ +10 km s⁻¹) validates pump normalization, partition functions and branching for each species.
- The 4.63/4.85 µm H2O hot-band lines (10⁻⁷–5 × 10⁻⁷ s⁻¹ each) are the direct check of the cascade implementation.
