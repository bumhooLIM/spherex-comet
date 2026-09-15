"""Band-model line lists for species absent from HITRAN, following the GSFC papers.

CH3OH nu3 (Villanueva et al. 2012, ApJ 747, 37): symmetric-top rotational structure with torsional
K ladders (A, E1, E2), parallel-band Honl-London factors, S(nu3, 296 K) = 31.01e-19 cm/molecule,
Einstein A via Simeckova et al. (2006) with the *total* partition sum Z_tot(296 K) = 36761.
CH3OH nu2 / nu9: same rotational scheme, perpendicular-band Honl-London factors, intensities from the
Rogers (1980) split of the C-H stretch region (28 % / 50 % of 211e-19) quoted in the same paper.
Band centres of nu2 and nu9 are approximate (3000 and 2970 cm-1) -- see doc/fluorescence_database.md.
HC3N nu1 (Villanueva et al. 2013, JQSRT 129, 158): Sigma-Sigma linear band, constants of that paper,
band intensity 249.4 cm-2 atm-1 at 296 K treated as the fundamental with Q_vib(296 K) = 3.03.
"""
import numpy as np
import pandas as pd

C2 = 1.438776877
C_CM = 2.99792458e10
T0 = 296.0

# ---------------------------------------------------------------- CH3OH constants (Villanueva 2012)
A_ROT, B_ROT, C_ROT = 4.2537233, 0.8236523, 0.7925575
BBAR = 0.5 * (B_ROT + C_ROT)
TORS_G = (133.84356, 0.23722, 0.027307, -6.1352, 0.39735)      # Eq. 2, ground state
TORS_V3 = (134.83284, -0.16811, 0.013409, -6.2039, 0.38327)    # Eq. 2, nu3 state
NU0_V3 = 2844.2
DB_MEDIAN = -0.0048371871
Z_TOT_296 = 36761.0
Z_ROT_296_PAPER = 24558.0
IA_CH3OH = 0.98593
S_V3_296 = 31.01e-19                     # cm molecule^-1, fundamental after removing hot bands
S_CH_TOTAL = 211e-19                     # cm molecule^-1, whole C-H stretch region at 296 K
F_HOT_296 = 0.33
S_V2_296 = 0.28 * S_CH_TOTAL * (1 - F_HOT_296)
S_V9_296 = 0.50 * S_CH_TOTAL * (1 - F_HOT_296)
NU0_V2, NU0_V9 = 3000.0, 2970.0          # approximate band centres (not from the paper)
# Table 1 of Villanueva 2012: (delta_nu, delta_B) per K' for the columns E1, A+, A-, E2 (K' <= 7)
TABLE1 = {
    0: {"E1": (0.472427, -0.000163), "A": (0.515070, -0.002999)},
    1: {"E1": (0.232969, -0.002845), "A+": (1.644089, 0.001430), "A-": (1.648615, -0.012424), "E2": (0.240550, -0.016448)},
    2: {"E1": (0.283519, -0.010409), "A+": (-0.109288, -0.004837), "A-": (-0.117130, -0.004596), "E2": (0.495529, -0.006310)},
    3: {"E1": (0.053291, -0.004447), "A+": (-2.106990, -0.005146), "A-": (-2.111880, -0.005093), "E2": (1.663912, -0.008170)},
    4: {"E1": (0.189937, DB_MEDIAN), "A+": (0.730155, -0.007722), "A-": (0.730167, -0.007725), "E2": (0.225078, 0.007489)},
    5: {"E1": (-0.069808, DB_MEDIAN), "A+": (0.565761, DB_MEDIAN), "A-": (0.565761, DB_MEDIAN), "E2": (0.104760, DB_MEDIAN)},
    6: {"E1": (-0.046387, DB_MEDIAN), "A+": (-0.264249, -0.001436), "A-": (-0.281261, -0.001261), "E2": (0.911340, DB_MEDIAN)},
    7: {"E1": (1.104191, DB_MEDIAN), "A+": (-0.092171, DB_MEDIAN), "A-": (-0.092171, DB_MEDIAN), "E2": (-0.160811, DB_MEDIAN)},
}


def tors_energy(coef, K, tau):
    c0, c1, c2, c3, c4 = coef
    return c0 + c1 * K + c2 * K ** 2 + c3 * np.cos(c4 * K + 2 * np.pi / 3 * (tau - 1))


def tors_symmetry(K, tau):
    return {0: "E1", 1: "A", 2: "E2"}[(K + tau) % 3]


def ch3oh_levels(Kmax=15, Jmax=30):
    """Ground-state torsion-rotation levels (J, K, tau, sym, E'').  A levels are split into A+/A-
    (degenerate here) for K > 0.  Weight 4(2J+1) per level (paper convention)."""
    e00 = tors_energy(TORS_G, 0, 1)            # J=0,K=0,A reference (tau=1 -> A for K=0)
    rows = []
    for K in range(Kmax + 1):
        for tau in (1, 2, 3):
            sym = tors_symmetry(K, tau)
            comps = ["A+", "A-"] if (sym == "A" and K > 0) else [sym]
            et = tors_energy(TORS_G, K, tau) - e00
            for J in range(K, Jmax + 1):
                E = et + BBAR * J * (J + 1) + (A_ROT - BBAR) * K ** 2
                for c in comps:
                    rows.append((J, K, tau, c, E))
    return pd.DataFrame(rows, columns=["J", "K", "tau", "sym", "E"])


def _delta_v3(K, tau, sym):
    """(delta_nu, delta_B) of the nu3 K sub-band for torsional species sym."""
    if K in TABLE1:
        t = TABLE1[K]
        key = sym if sym in t else ("A" if sym in ("A+", "A-") else sym)
        if key in t:
            return t[key]
    dnu = tors_energy(TORS_V3, K, tau) - tors_energy(TORS_G, K, tau)
    return dnu, DB_MEDIAN


def _hl_parallel(J, K, dJ):
    if dJ == 1:
        return ((J + 1) ** 2 - K ** 2) / ((J + 1) * (2 * J + 1))
    if dJ == 0:
        return 0.0 if J == 0 else K ** 2 / (J * (J + 1))
    return 0.0 if J == 0 else (J ** 2 - K ** 2) / (J * (2 * J + 1))


def _hl_perp(J, K, dJ, dK):
    """Perpendicular-band Honl-London factors (Herzberg), K signed by dK; normalised later."""
    s = dK
    if dJ == 1:
        return (J + 2 + s * K) * (J + 1 + s * K) / (2 * (J + 1) * (2 * J + 1))
    if dJ == 0:
        return 0.0 if J == 0 else (J + 1 + s * K) * (J - s * K) / (2 * J * (J + 1))
    return 0.0 if J == 0 else (J - 1 - s * K) * (J - s * K) / (2 * J * (2 * J + 1))


def _einstein_A(nu, S296, El, wu, Ztot, Ia):
    return 8 * np.pi * C_CM * nu ** 2 * Ztot * S296 / ((1 - np.exp(-C2 * nu / T0)) * np.exp(-C2 * El / T0) * Ia * wu)


def ch3oh_band(which="nu3", Kmax=15, Jmax=30):
    lev = ch3oh_levels(Kmax, Jmax)
    lev["w"] = 4 * (2 * lev["J"] + 1)
    Zrot = float(np.sum(lev["w"] * np.exp(-C2 * lev["E"] / T0)))
    key = {(r.J, r.K, r.tau, r.sym): r.E for r in lev.itertuples()}
    rows = []
    if which == "nu3":
        S0, nu0 = S_V3_296, NU0_V3
        for r in lev.itertuples():
            J, K, tau, sym, El, wl = r.J, r.K, r.tau, r.sym, r.E, r.w
            dnu, dB = _delta_v3(K, tau, sym)
            for dJ in (-1, 0, 1):
                Ju = J + dJ
                if Ju < K or Ju > Jmax:
                    continue
                hl = _hl_parallel(J, K, dJ)
                if hl <= 0:
                    continue
                # A+ <-> A- for dJ = 0, A+ <-> A+ otherwise (degenerate here); E <-> E
                symu = sym
                if sym in ("A+", "A-") and dJ == 0:
                    symu = "A-" if sym == "A+" else "A+"
                Eu = key[(Ju, K, tau, symu)] + nu0 + dnu + dB * Ju * (Ju + 1)
                nu = Eu - El
                wu = 4 * (2 * Ju + 1)
                S = S0 * (nu / nu0) * hl * wl * np.exp(-C2 * El / T0) * (1 - np.exp(-C2 * nu / T0)) / Zrot
                rows.append(dict(nu=nu, a=_einstein_A(nu, S, El, wu, Z_TOT_296, IA_CH3OH), elower=El, gp=wu, gpp=wl,
                                 vup="V3", vlow="GROUND", lup=f"{Ju} {K} {tau} {symu}", llow=f"{J} {K} {tau} {sym}"))
    else:
        S0, nu0 = (S_V2_296, NU0_V2) if which == "nu2" else (S_V9_296, NU0_V9)
        vlab = "V2" if which == "nu2" else "V9"
        for r in lev.itertuples():
            J, K, tau, sym, El, wl = r.J, r.K, r.tau, r.sym, r.E, r.w
            cand = []
            for dK in (-1, 1):
                Ku = K + dK
                if Ku < 0 or Ku > Kmax:
                    continue
                for dJ in (-1, 0, 1):
                    Ju = J + dJ
                    if Ju < Ku or Ju < 0 or Ju > Jmax:
                        continue
                    hl = _hl_perp(J, K, dJ, dK)
                    if hl <= 0:
                        continue
                    # torsional species: keep tau, symmetry follows (K+tau) mod 3 of the upper ladder
                    symu = tors_symmetry(Ku, tau)
                    if symu == "A" and Ku > 0:
                        symu = "A+"
                    cand.append((Ju, Ku, dJ, dK, hl, symu))
            tot = sum(c[4] for c in cand)
            if tot <= 0:
                continue
            for Ju, Ku, dJ, dK, hl, symu in cand:
                Eu = key[(Ju, Ku, tau, symu)] + nu0
                nu = Eu - El
                if nu <= 0:
                    continue
                wu = 4 * (2 * Ju + 1)
                S = S0 * (nu / nu0) * (hl / tot) * wl * np.exp(-C2 * El / T0) * (1 - np.exp(-C2 * nu / T0)) / Zrot
                rows.append(dict(nu=nu, a=_einstein_A(nu, S, El, wu, Z_TOT_296, IA_CH3OH), elower=El, gp=wu, gpp=wl,
                                 vup=vlab, vlow="GROUND", lup=f"{Ju} {Ku} {tau} {symu}", llow=f"{J} {K} {tau} {sym}"))
    df = pd.DataFrame(rows)
    df.attrs["Zrot296"] = Zrot
    return df


# ---------------------------------------------------------------- HC3N nu1 (Villanueva 2013)
HC3N = dict(B0=0.1574028, D0=1.8165e-8, H0=2.87e-15, nu0=3327.37085, B1=0.15149762, D1=1.8065e-8, H1=2.87e-15,
            S296=1.006e-17, fundamentals=[(3327, 1), (2274, 1), (2079, 1), (862, 1), (663, 2), (499, 2), (222, 1 * 2)])


def qvib(fundamentals, T=T0):
    q = 1.0
    for nu, d in fundamentals:
        q *= (1 - np.exp(-C2 * nu / T)) ** (-d)
    return q


def hc3n_band(Jmax=120):
    p = HC3N
    E = lambda J, B, D, H: B * J * (J + 1) - D * (J * (J + 1)) ** 2 + H * (J * (J + 1)) ** 3
    Js = np.arange(0, Jmax + 1)
    w = 2 * Js + 1
    Eg = E(Js, p["B0"], p["D0"], p["H0"])
    Qrot = float(np.sum(w * np.exp(-C2 * Eg / T0)))
    Ztot = Qrot * qvib(p["fundamentals"])
    rows = []
    for J in Js:
        for dJ, hl in ((1, (J + 1) / (2 * J + 1)), (-1, J / (2 * J + 1))):
            Ju = J + dJ
            if Ju < 0 or hl <= 0:
                continue
            El, Eu = Eg[J], p["nu0"] + E(Ju, p["B1"], p["D1"], p["H1"])
            nu = Eu - El
            wl, wu = 2 * J + 1, 2 * Ju + 1
            S = p["S296"] * (nu / p["nu0"]) * hl * wl * np.exp(-C2 * El / T0) * (1 - np.exp(-C2 * nu / T0)) / Qrot
            rows.append(dict(nu=nu, a=_einstein_A(nu, S, El, wu, Ztot, 1.0), elower=El, gp=wu, gpp=wl,
                             vup="V1", vlow="GROUND", lup=f"J={Ju}", llow=f"J={J}"))
    df = pd.DataFrame(rows)
    df.attrs["Qrot296"] = Qrot; df.attrs["Qvib296"] = Ztot / Qrot
    return df


def as_gfm_frame(df):
    """Add the columns gfm_core.prepare expects (labelled both-given format)."""
    df = df.copy()
    df["eup"] = df["elower"] + df["nu"]
    df["low_key"] = df["vlow"] + "|" + df["llow"]
    df["up_key"] = df["vup"] + "|" + df["lup"]
    df["up_labeled"] = True
    df["low_labeled"] = True
    return df
