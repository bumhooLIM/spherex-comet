"""Reconstruction of the GSFC General Fluorescence Model (Villanueva et al. 2011 App. C; 2012 cascade).

Computes solar-pumped fluorescence efficiencies (g-factors, photons s^-1 molecule^-1 at 1 au) for
every line of a HITRAN-format line list, with level-by-level cascade: each upper level receives the
direct solar pump from the ground vibrational state plus the cascade from higher levels, and decays
with branching ratios A_i / A_tot.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

C2 = 1.438776877          # cm K   (second radiation constant)
H_PLANCK = 6.62607015e-34
C_LIGHT = 2.99792458e8    # m/s
C_CM = 2.99792458e10      # cm/s
KB = 1.380649e-23

LINEAR = {"CO", "CO2", "HCN", "C2H2", "OCS", "N2O", "HC3N"}
ASYM = {"H2O", "HDO", "H2CO", "H2S", "SO2", "C2H4"}
GROUND_EMAX = 250.0        # a vibrational label whose lowest level lies below this is a ground-state label
MIN_LINES = 20             # labels with fewer lines are not trusted to define vibrational states

_ws = re.compile(r"\s+")


def _norm(s):
    return _ws.sub(" ", str(s).strip())


def _asym_key(s):
    m = re.match(r"^\s*(\d+)\s+(\d+)\s+(\d+)", s)
    return f"{m.group(1)} {m.group(2)} {m.group(3)}" if m else ""


def _ch4_key(s):
    m = re.match(r"^\s*(\d+)\s*(A1|A2|E|F1|F2)\s*(\d*)", s)
    if not m:
        return ""
    n = m.group(3) or "1"
    return f"{m.group(1)}{m.group(2)} {n}"


def _nh3_key(s):
    m = re.match(r"^\s*(\d+)\s+(\d+)\s+([sa])\b", s)
    return f"{m.group(1)} {m.group(2)} {m.group(3)}" if m else ""


def _numeric_core(label):
    """Vibrational label with symmetry tokens removed ('0 0 0 0 1A1' -> '0 0 0 0')."""
    return " ".join(t for t in label.split() if not re.search(r"[A-Za-z]", t))


# low-lying fundamentals (cm-1, degeneracy) for the vibrational partition function at coma T
LOW_MODES = {
    "CO2": [(667.4, 2)], "OCS": [(520.4, 2), (859.0, 1)], "C2H2": [(612.9, 2), (730.3, 2)],
    "HCN": [(712.0, 2)], "C2H6": [(289.0, 1), (821.7, 2), (822.0, 1)], "N2O": [(588.8, 2)],
    "C2H4": [(826.0, 1), (943.0, 1), (949.0, 1)], "CH3OH": [(294.45, 1)],
    "HC3N": [(222.0, 2), (499.0, 2), (663.0, 2), (862.0, 1)],
}


def q_vib(species, T):
    q = 1.0
    for nu, d in LOW_MODES.get(species, []):
        q *= (1.0 - np.exp(-C2 * nu / T)) ** (-d)
    return q


def q_rot_ch4(T, B=5.2410, D=1.10e-4, Jmax=60):
    """Rigid spherical-top partition sum of CH4 with Td spin statistics (A:5, E:2, F:3) from the
    reduction of D^J onto the octahedral group; HITRAN's CH4 level labels are not consistent enough
    to sum directly."""
    Q = 0.0
    classes = [(1, 0.0), (8, 2 * np.pi / 3), (3, np.pi), (6, np.pi / 2), (6, np.pi)]      # (N_c, angle)
    chars = {"A1": [1, 1, 1, 1, 1], "A2": [1, 1, 1, -1, -1], "E": [2, -1, 2, 0, 0],
             "F1": [3, 0, -1, 1, -1], "F2": [3, 0, -1, -1, 1]}
    gns = {"A1": 5, "A2": 5, "E": 2, "F1": 3, "F2": 3}
    for J in range(Jmax + 1):
        chi = [(2 * J + 1) if th == 0 else np.sin((2 * J + 1) * th / 2) / np.sin(th / 2) for _, th in classes]
        E = B * J * (J + 1) - D * (J * (J + 1)) ** 2
        for G_, ch in chars.items():
            n = sum(N * c * x for (N, _), c, x in zip(classes, ch, chi)) / 24.0
            n = int(round(n))
            Q += n * gns[G_] * (2 * J + 1) * np.exp(-C2 * E / T)
    return Q


def q_rot_asym(T, A, B, C, sigma, gns_total):
    """Classical rigid asymmetric-top partition sum (used for C2H4, whose HITRAN labels double count)."""
    kT = T / C2
    return gns_total / sigma * np.sqrt(np.pi * kT ** 3 / (A * B * C))


Q_OVERRIDE = {"CH4": lambda T: q_rot_ch4(T),
              "C2H4": lambda T: q_rot_asym(T, 4.8646, 1.0011, 0.8282, 4, 16)}


def load_hitran_csv(path, species):
    """Read the astroquery HITRAN CSV and build level keys.

    Returns a DataFrame with nu, a, elower, gp, gpp, vup, vlow, up_key, low_key, up_labeled,
    low_labeled, eup.  Keys are normalised per molecule class so that one physical level has one
    key whatever source list the line came from (HITRAN 2020 mixes formats: flag characters on
    H2O levels, missing sub-level indices for CH4, hyperfine suffixes for HCN, inversion-symmetry
    suffixes on the NH3 vibrational label).
    """
    df = pd.read_csv(path, dtype={"global_upper_quanta": str, "global_lower_quanta": str,
                                  "local_upper_quanta": str, "local_lower_quanta": str},
                     keep_default_na=False,
                     usecols=["nu", "a", "elower", "gp", "gpp", "global_upper_quanta",
                              "global_lower_quanta", "local_upper_quanta", "local_lower_quanta"])
    df = df.rename(columns={"global_upper_quanta": "vup", "global_lower_quanta": "vlow",
                            "local_upper_quanta": "lup", "local_lower_quanta": "llow"})
    for c in ("vup", "vlow", "lup", "llow"):
        df[c] = df[c].map(_norm)
    df = df[df["elower"] >= 0].reset_index(drop=True)      # negative = unknown lower energy
    df["eup"] = df["elower"] + df["nu"]
    if species in LINEAR:
        m = df["llow"].str.extract(r"^([PQR])\s*(\d+)\s*([ef]?)")
        br, jl, sym = m[0], pd.to_numeric(m[1], errors="coerce"), m[2].fillna("")
        ok = br.notna() & jl.notna()
        dj = br.map({"P": -1, "Q": 0, "R": 1})
        ju = jl + dj
        symu = pd.Series(np.where(dj == 0, sym.map({"e": "f", "f": "e", "": ""}), sym), index=df.index)
        jl_s = jl.astype("Int64").astype(str)
        ju_s = ju.astype("Int64").astype(str)
        df["low_key"] = np.where(ok, df["vlow"] + "|J=" + jl_s + sym, "")
        df["up_key"] = np.where(ok, df["vup"] + "|J=" + ju_s + symu, "")
        df["up_labeled"] = ok & (df["vup"] != "")
        df["low_labeled"] = ok & (df["vlow"] != "")
    else:
        if species in ASYM:
            lup, llow = df["lup"].map(_asym_key), df["llow"].map(_asym_key)
        elif species in ("CH4", "CH3D"):
            lup, llow = df["lup"].map(_ch4_key), df["llow"].map(_ch4_key)
        elif species == "NH3":
            lup, llow = df["lup"].map(_nh3_key), df["llow"].map(_nh3_key)
        else:
            lup, llow = df["lup"], df["llow"]
        df["low_key"] = df["vlow"] + "|" + llow
        df["up_key"] = df["vup"] + "|" + lup
        df["up_labeled"] = (df["vup"] != "") & (lup != "")
        df["low_labeled"] = (df["vlow"] != "") & (llow != "")
    df.loc[~df["up_labeled"], "up_key"] = ""
    df.loc[~df["low_labeled"], "low_key"] = ""
    return df


def ground_labels(df):
    """Vibrational labels of the ground state: every labelled vlow whose lowest level is below
    GROUND_EMAX (several spellings of the ground state can coexist in HITRAN)."""
    d = df.loc[df["low_labeled"]]
    lab = d.groupby("vlow")["elower"].agg(emin="min", nlow=lambda e: int(np.sum(e < GROUND_EMAX)), n="size")
    ground = set(lab.index[(lab["emin"] < GROUND_EMAX) & (lab["nlow"] >= MIN_LINES)])
    cores = {_numeric_core(g) for g in ground}
    ground |= {l for l in lab.index if _numeric_core(l) in cores and lab.loc[l, "emin"] < GROUND_EMAX}
    minor = set(lab.index[lab["n"] < MIN_LINES]) - ground
    return ground, minor


def prepare(df, species):
    """Attach ground-state flags, level tables and pseudo-levels for unlabeled upper states."""
    gset, minor = ground_labels(df)
    lab = df.loc[df["low_labeled"] & ~df["vlow"].isin(gset | minor)].groupby("vlow")["elower"].min()
    lab = lab[lab >= GROUND_EMAX]
    e_vib1 = float(lab.min()) if len(lab) else np.inf
    df = df.copy()
    untrusted = (df["vlow"] == "") | df["vlow"].isin(minor)
    df["low_ground"] = df["vlow"].isin(gset) | (untrusted & (df["elower"] < e_vib1))
    df["up_ground"] = df["vup"].isin(gset)
    # one key per ground level whatever the spelling of the ground vibrational label
    for g in gset:
        m = df["low_ground"] & (df["vlow"] == g) & df["low_labeled"]
        loc = df.loc[m, "low_key"].str.split("|").str[1]
        if species in LINEAR:
            loc = loc.str.replace(r"[ef]$", "", regex=True)
        df.loc[m, "low_key"] = "GROUND|" + loc
    unl = ~df["up_labeled"]
    df.loc[unl, "up_key"] = ["__L%d" % i for i in np.flatnonzero(unl.to_numpy())]
    df = df[~df["up_ground"]].reset_index(drop=True)          # drop pure-rotational ground lines
    lev = df.groupby("up_key").agg(E=("eup", "median"), Atot=("a", "sum"), n=("a", "size"))
    lev = lev.sort_values("E", ascending=False)
    lev["idx"] = np.arange(len(lev))
    df["up_idx"] = df["up_key"].map(lev["idx"]).astype(int)
    # cascade target: the lower level as an upper level of other lines; try e/f variants when the
    # parity label is missing; -1 for ground or unknown
    idx = df["low_key"].map(lev["idx"])
    for suffix in ("e", "f"):
        alt = (df["low_key"] + suffix).map(lev["idx"])
        idx = idx.fillna(alt)
    df["low_idx"] = idx.fillna(-1).astype(int)
    df.loc[df["low_ground"], "low_idx"] = -1
    gl = df.loc[df["low_ground"] & df["low_labeled"]].groupby("low_key").agg(E=("elower", "median"), w=("gpp", "max"))
    gl = gl[gl["w"] > 0]
    meta = dict(species=species, ground=sorted(gset), e_vib1=e_vib1, ground_levels=gl, levels=lev)
    return df, meta


def partition_function(meta, T, vib=False):
    """Rotational partition sum of the ground vibrational state at T (explicit level sum, or the
    analytic override for species whose HITRAN labels double count); times q_vib(T) if vib."""
    sp = meta["species"]
    if sp in Q_OVERRIDE:
        q = float(Q_OVERRIDE[sp](T))
    else:
        gl = meta["ground_levels"]
        q = float(np.sum(gl["w"].to_numpy() * np.exp(-C2 * gl["E"].to_numpy() / T)))
    return q * q_vib(sp, T) if vib else q


def gfactors(df, meta, T, solar_fnu, v_h_kms=0.0):
    """Line g-factors [photons s^-1 molecule^-1] at 1 au.

    solar_fnu(nu_cm) -> F_nu [W m^-2 Hz^-1] at 1 au, sampled at the rest wavenumber the solar
    photons had in the Sun's frame: nu_sun = nu_line * (1 + v_h/c), v_h > 0 receding.
    """
    Q = partition_function(meta, T, vib=True)
    nu = df["nu"].to_numpy(float)
    A = df["a"].to_numpy(float)
    gp = df["gp"].to_numpy(float)
    gpp = df["gpp"].to_numpy(float)
    El = df["elower"].to_numpy(float)
    pumpable = df["low_ground"].to_numpy() & (nu > 0)
    nu_sun = nu * (1.0 + v_h_kms * 1e3 / C_LIGHT)
    fnu = solar_fnu(nu_sun)
    nu_hz = nu * C_CM
    ok = pumpable & (gpp > 0) & (A > 0) & np.isfinite(fnu)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = (gp / gpp) * A * C_LIGHT ** 2 * fnu / (8.0 * np.pi * H_PLANCK * nu_hz ** 3)
        pop = gpp * np.exp(-C2 * El / T) / Q
    g_pump = np.where(ok, rate * pop, 0.0)
    g_pump = np.nan_to_num(g_pump, nan=0.0, posinf=0.0, neginf=0.0)
    # cascade: process upper levels in descending energy
    nlev = len(meta["levels"])
    up_idx = df["up_idx"].to_numpy()
    low_idx = df["low_idx"].to_numpy()
    Atot = meta["levels"]["Atot"].to_numpy(float)
    P = np.bincount(up_idx, weights=g_pump, minlength=nlev)      # direct pumps per level
    Cin = np.zeros(nlev)                                          # cascade inflow
    order = np.argsort(up_idx, kind="stable")                     # lines grouped by level index
    bounds = np.searchsorted(up_idx[order], np.arange(nlev + 1))
    g = np.zeros(len(df))
    A_ord = A[order]; low_ord = low_idx[order]
    for L in range(nlev):                                         # levels already sorted by E desc
        s, e = bounds[L], bounds[L + 1]
        if s == e:
            continue
        flow = P[L] + Cin[L]
        if flow <= 0.0 or Atot[L] <= 0.0:
            continue
        gl = flow * A_ord[s:e] / Atot[L]
        g[order[s:e]] = gl
        tgt = low_ord[s:e]
        m = tgt >= 0
        if m.any():
            np.add.at(Cin, tgt[m], gl[m])
    return g, g_pump, Q


def bin_profile(lam_um, g, grid_um):
    """Sum line g-factors into a wavelength grid; returns g per micron (photons s^-1 molec^-1 um^-1)."""
    edges = np.sqrt(grid_um[:-1] * grid_um[1:])
    edges = np.concatenate([[grid_um[0] ** 2 / edges[0]], edges, [grid_um[-1] ** 2 / edges[-1]]])
    hist, _ = np.histogram(lam_um, bins=edges, weights=g)
    return hist / np.diff(edges)


def planck_fnu_1au(nu_cm, T=5772.0):
    """Blackbody solar irradiance at 1 au, W m^-2 Hz^-1 (for tests only)."""
    nu = np.asarray(nu_cm, float) * C_CM
    R_SUN, AU = 6.957e8, 1.495978707e11
    B = 2 * H_PLANCK * nu ** 3 / C_LIGHT ** 2 / np.expm1(H_PLANCK * nu / (KB * T))
    return np.pi * B * (R_SUN / AU) ** 2


def linear_band_lines(B0, D0, H0, nu0, B1, D1, H1, S_band_296, Jmax=80, iso_abund=1.0, weight_even=1.0, weight_odd=1.0):
    """Σ–Σ band of a linear molecule (P and R branches) as a HITRAN-like DataFrame.

    S_band_296 in cm molecule^-1 at 296 K.  Hönl–London: R: (J''+1)/(2J''+1), P: J''/(2J''+1).
    """
    T0 = 296.0
    E = lambda J, B, D, H: B * J * (J + 1) - D * (J * (J + 1)) ** 2 + H * (J * (J + 1)) ** 3
    rows = []
    Js = np.arange(0, Jmax + 1)
    w = np.where(Js % 2 == 0, weight_even, weight_odd) * (2 * Js + 1)
    Eg = E(Js, B0, D0, H0)
    Q0 = np.sum(w * np.exp(-C2 * Eg / T0))
    for J in Js:
        for dJ, hl in ((1, (J + 1) / (2 * J + 1)), (-1, J / (2 * J + 1))):
            Ju = J + dJ
            if Ju < 0:
                continue
            wl = w[J]; wu = (weight_even if Ju % 2 == 0 else weight_odd) * (2 * Ju + 1)
            El = Eg[J]; Eu = nu0 + E(Ju, B1, D1, H1)
            nu = Eu - El
            if nu <= 0:
                continue
            S = S_band_296 * (nu / nu0) * hl * wl * np.exp(-C2 * El / T0) * (1 - np.exp(-C2 * nu / T0)) / Q0
            A = 8 * np.pi * C_CM * nu ** 2 * Q0 * S / ((1 - np.exp(-C2 * nu / T0)) * np.exp(-C2 * El / T0) * iso_abund * wu)
            rows.append(dict(nu=nu, a=A, elower=El, gp=wu, gpp=wl, vup="V1", vlow="GROUND",
                             lup="", llow=("R" if dJ > 0 else "P") + " %2d" % J))
    return pd.DataFrame(rows)
