# Villanueva, DiSanti, Mumma & Xu (2012) — A quantum band model of the ν3 fundamental of methanol (CH3OH) and its application to fluorescence spectra of comets

*ApJ 747, 37 (11 pp). doi:10.1088/0004-637X/747/1/37*

**Role in this project:** the only source for the 3.52 µm methanol band, which HITRAN does not contain (HITRAN CH3OH stops at 1407 cm⁻¹). Gives the constants needed to regenerate a ν3 line list and the band g-factors to validate it. Methanol dominates the C–H stretching region of many comets and blends with C2H6/CH4 at SPHEREx resolution.

## 1. Band model (Sect. 2)
- Methanol is a slightly asymmetric top with internal rotation; levels are labelled in the C3v/G6 scheme (J, K, T_s) with torsional symmetry T_s = A, E1, E2, A split into A± for K > 0. Torsional index τ = 1, 2, 3 maps to T_s via mod(K+τ, 3) = 0 (E1), 1 (A), 2 (E2).
- Ground-state energies from Mekhtiev et al. (1999) and Xu et al. (2008); rotational constants A = 4.2537233, B = 0.8236523, C = 0.7925575 cm⁻¹. Torsional energy per K ladder (Eq. 1): E_tors(K,τ) = E(J=K,τ) − [(B+C)/2 · J(J+1) + (A − (B+C)/2) K²]. Harmonic torsional pattern (Eq. 2): E_tors = c0 + c1 K + c2 K² + c3 cos(c4 K + 2π(τ−1)/3), with (Table 2) ground: c0 = 133.84356, c1 = 0.23722, c2 = 0.027307, c3 = −6.1352, c4 = 0.39735; ν3: 134.83284, −0.16811, 0.013409, −6.2039, 0.38327 (cm⁻¹, relative to the potential well, which lies 128.10687 cm⁻¹ below J = 0, K = 0, A+).
- Upper state (Eq. 3): E′(J,K,τ) = E″(J,K,τ) + ν0 + δν(K,τ) + δB(K,τ) J(J+1), ν0 = 2844.2 cm⁻¹; 19 sub-band origins δν and 19 δB for K′ ≤ 6 from 709 assignments (Table 1; K′ = 3 strongly perturbed); other ladders use the torsional diagram and the median δB = −0.0048371871 cm⁻¹.
- Partition functions (Sect. 2.2, Table 3, Fig. 4): explicit sum over 2.2 million generated levels; Z_tot(296 K) = 36 761, Z_rot = 24 558, Z_tors = 1.46800, Z_vib = 1.01971; agrees with JPL (Pearson & Xu 2010, revised) to 0.2 %; HITRAN's parsum.dat value (82 469) is inconsistent with its molparam value (35 314).
- Intensities (Sect. 2.3): total C–H stretch region 2650–3175 cm⁻¹ = 211 × 10⁻¹⁹ cm molecule⁻¹ (Rogers 1980, confirmed with PNNL/Sharpe & Sams), split ν2 : ν3 : ν9 = 28 % : 22 % : 50 %. Hot bands contribute f_hot = 33 % at 296 K (first torsional level 294.45 cm⁻¹), so the ν3 fundamental intensity is S_ν3(296 K) = 0.22 × 211 × (1 − 0.33) = 31.01 × 10⁻¹⁹ cm molecule⁻¹. Line intensities (Eq. 4) use parallel-band Hönl–London factors (Eq. 5), w″ = 4(2J″+1) (spin multiplicity 4) and Z_rot; selection rules ΔJ = 0, ±1, ΔK = 0, A± ↔ A± (ΔJ ≠ 0), A± ↔ A∓ (ΔJ = 0), E ↔ E. Einstein A from Šimečková et al. (2006) with I_a = 0.98593.
- Atlas (Sect. 2.4): 2 948 lines, J_max = 23, 2802–2892 cm⁻¹, HITRAN-2008 format (V3/GROUND; J, K, T_s); γ_air = 0.100, γ_self = 0.400 cm⁻¹ atm⁻¹, n_air = 0.75.

## 2. Fluorescence efficiencies (Sect. 3, Table 3; Fig. 3)
- Resonant fluorescence only (no data for cascades from higher levels); realistic solar pump (Kurucz continuum × Hase line list, differential rotation) via the GFM.
- Band g-factors (1 au, photons s⁻¹ molecule⁻¹): 143.4 × 10⁻⁶ at 75 K, 128.4 × 10⁻⁶ at 150 K, 97.3 × 10⁻⁶ at 296 K. The old Crovisier & Encrenaz-style band value was 1.47 × 10⁻⁴ (Bockelée-Morvan et al. 1995, valid only for T < 100 K).
- Q-branch (2843.5–2845.0 cm⁻¹) g-factor: 12.9, 14.4, 15.5, 15.5, 14.6 × 10⁻⁶ at 10, 40, 75, 100, 130 K. The Q branch is a compact ~1.5 cm⁻¹ feature with K ladders (K = 2, 3(A), 4, 5) resolved at R ≥ 50 000; the P and R branches extend over ~2810–2880 cm⁻¹ (Fig. 3).
- The "simplified model" (identical rotational structure in ground and ν3 states) that earlier authors used mislocates the sub-branches by several cm⁻¹ but gives a similar band envelope — the level of fidelity that matters at SPHEREx resolution.

## 3. Cometary validation (Sect. 3.1–3.4, Table 4)
- C/2001 A2 (LINEAR), C/2004 Q2 (Machholz; NIRSPEC 2005-01-19, T_rot = 75 K), 8P/Tuttle (NIRSPEC 2007-12-22 at 45 K; CRIRES 2008-01-29 and 02-03 at 55 K, R up to 98 000). Line-by-line agreement at all three resolutions.
- 103P/Hartley 2: two teams' CH3OH abundances differed by ~2×; with the new Q-branch g-factor 1.55 × 10⁻⁵ s⁻¹ both give ~1.8 % relative to water. The new model reconciles IR and radio production rates in C/2001 A2.

## 4. Implications for the SPHEREx reconstruction
- Regenerate ν3 with the constants above (symmetric-top rotational structure, torsional ladders, median δB) and S_ν3 = 31.01 × 10⁻¹⁹ cm molecule⁻¹; validate the band total against 143.4 × 10⁻⁶ s⁻¹ at 75 K.
- ν2 and ν9 (the other 78 % of the C–H region, at ~3000 and ~2970 cm⁻¹) are not modelled in this paper; a band-model treatment with the Rogers intensity split is the best available approximation and must be labelled as such.
