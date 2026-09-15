"""Build the reconstructed GSFC-style cometary fluorescence database for SPHEREx (0.7-5.0 um).

Outputs (data/fluorescence/):
  band_gfactors.csv          band-integrated g-factors [photons s^-1 molecule^-1 at 1 au] per T_rot
  species_summary.csv        per-species totals in 0.7-5.0 um, partition sums, provenance
  swings_sensitivity.csv     band g-factors at 70 K versus heliocentric velocity (Swings effect)
  co_swings.csv              CO v(1-0) g-factor on a fine velocity grid (0.25 km/s, -60..+60) at 30/70/130 K
  solar_flux_1au.csv         the pumping spectrum, 1 cm-1 sampling
  profiles/<species>_gprofile.csv   g per micron on a common R = 5000 grid, one column per T_rot
  profiles/all_species_T070K.csv    the 70 K profiles of every species side by side
  lines/<species>_lines.csv         line list (0.7-5.0 um) with g at every T_rot

Run:  python notebooks/fluorescence_gfm/build_fluorescence_db.py [--species H2O CO2 ...]
Inputs are read from data/fluorescence/inputs/ (HITRAN subsets are fetched with astroquery when
missing).  See doc/fluorescence_database.md for the method and data/fluorescence/README.md for
the formats.
"""
from __future__ import annotations
import argparse, os, sys, time, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # notebooks/comspec/fluorescence_gfm/ -> project
sys.path.insert(0, HERE)
import gfm_core as G                      # noqa: E402
import gfm_bandmodels as BM               # noqa: E402
import solar_pump as SP                   # noqa: E402

OUT = os.path.join(ROOT, "data", "fluorescence")
INPUTS = os.path.join(OUT, "inputs")
SPECIES = ["H2O", "CO2", "CO", "CH4", "C2H6", "C2H2", "C2H4", "CH3OH", "H2CO", "NH3", "HCN", "HC3N", "OCS", "H2S"]
HITRAN_ID = {"H2O": 1, "CO2": 2, "CO": 5, "CH4": 6, "NH3": 11, "HCN": 23, "C2H2": 26, "C2H6": 27,
             "H2CO": 20, "OCS": 19, "H2S": 31, "C2H4": 38}
BAND_MODEL = {"CH3OH", "HC3N"}
TEMPS = [30, 50, 70, 100, 130]
T_REF = 70
VH_GRID = [-30, -20, -10, -5, 0, 5, 10, 20, 30]
LAM_MIN, LAM_MAX, R_GRID = 0.70, 5.00, 5000
G_LINE_MIN = 1e-11            # lines below this at every T are not written to lines/
G_BAND_MIN = 1e-9             # bands below this at 70 K are not tabulated
# readable names for the bands of the main species (upper-lower vibrational labels)
BAND_NAMES = {
    "H2O": {"0 0 1-0 0 0": "nu3", "1 0 0-0 0 0": "nu1", "0 1 1-0 1 0": "nu2+nu3-nu2", "1 1 0-0 1 0": "nu1+nu2-nu2",
            "1 0 1-1 0 0": "nu1+nu3-nu1", "1 0 1-0 0 1": "nu1+nu3-nu3", "0 0 1-0 1 0": "nu3-nu2",
            "1 0 0-0 1 0": "nu1-nu2", "0 2 0-0 0 0": "2nu2", "1 1 1-1 0 0": "nu1+nu2+nu3-nu1",
            "2 0 0-0 0 1": "2nu1-nu3", "2 0 0-1 0 0": "2nu1-nu1", "0 1 1-0 0 0": "nu2+nu3", "0 2 1-0 1 0": "2nu2+nu3-nu2",
            "0 2 1-0 2 0": "2nu2+nu3-2nu2", "1 1 0-1 0 0": "nu1+nu2-nu1", "0 1 1-0 0 1": "nu2+nu3-nu3",
            "1 0 1-0 0 0": "nu1+nu3 (1.38um)", "0 1 1-0 2 0": "nu2+nu3-2nu2", "2 0 1-2 0 0": "2nu1+nu3-2nu1",
            "1 1 1-0 1 0": "nu1+nu2+nu3-nu2 (1.9um)", "0 2 1-0 0 0": "2nu2+nu3 (1.14um)"},
    "CH4": {"0 0 1 0 1F2-0 0 0 0 1A1": "nu3", "1 0 0 1 1F2-0 0 0 0 1A1": "nu1+nu4 (2.3um)", "0 0 1 1 1F2-0 0 0 1 1F2": "nu3+nu4-nu4"},
    "C2H6": {"V7-GROUND": "nu7", "V5-GROUND": "nu5"}, "CH3OH": {"V3-GROUND": "nu3", "V2-GROUND": "nu2 (approx.)", "V9-GROUND": "nu9 (approx.)"},
    "C2H2": {"001 0 0 0 0+ u-000 0 0 0 0+ g": "nu3", "010 1 1 0 0+ u-000 0 0 0 0+ g": "nu2+nu4+nu5"},
    "H2CO": {"0 0 0 0 1 0-0 0 0 0 0 0": "nu5", "1 0 0 0 0 0-0 0 0 0 0 0": "nu1"},
    "HCN": {"1 0 0 0-0 0 0 0": "nu1", "0 0 0 1-0 0 0 0": "nu3"}, "HC3N": {"V1-GROUND": "nu1"},
    "OCS": {"0 0 0 1-0 0 0 0": "nu3"}, "NH3": {"1000 00 0 A2\"-0000 00 0 A1'": "nu1 (a-s)", "1000 00 0 A1'-0000 00 0 A2\"": "nu1 (s-a)"},
    "H2S": {"0 1 1-0 0 0": "nu2+nu3", "1 1 0-0 0 0": "nu1+nu2"},
    "CO2": {"0 0 0 11-0 0 0 01": "nu3", "1 0 0 11-1 0 0 01": "nu3 hot (10011-10001)", "1 0 0 12-1 0 0 02": "nu3 hot (10012-10002)",
            "0 1 1 11-0 1 1 01": "nu3 hot (01111-01101)", "1 0 0 11-0 0 0 01": "nu1+nu3 (2.7um)", "1 0 0 12-0 0 0 01": "2nu2+nu3 (2.7um)",
            "0 1 1 01-0 0 0 01": "nu2"},
    "CO": {"1-0": "v(1-0)", "2-1": "v(2-1)", "2-0": "v(2-0)", "3-2": "v(3-2)"},
}


def hitran_cache(species):
    path = os.path.join(INPUTS, "hitran", f"{species}.csv.gz")
    if not os.path.exists(path):
        from astroquery.hitran import Hitran
        import astropy.units as u
        tbl = Hitran.query_lines(molecule_number=HITRAN_ID[species], isotopologue_number=1,
                                 min_frequency=0 * u.cm ** -1, max_frequency=25000 * u.cm ** -1)
        cols = ["nu", "sw", "a", "elower", "gp", "gpp", "global_upper_quanta", "global_lower_quanta",
                "local_upper_quanta", "local_lower_quanta"]
        tbl[cols].to_pandas().to_csv(path, index=False, compression="gzip")
    return path


def load_species(species):
    if species == "CH3OH":
        df = BM.as_gfm_frame(pd.concat([BM.ch3oh_band("nu3"), BM.ch3oh_band("nu2"), BM.ch3oh_band("nu9")], ignore_index=True))
        source = "band model (Villanueva et al. 2012 ApJ nu3; nu2/nu9 approximate)"
    elif species == "HC3N":
        df = BM.as_gfm_frame(BM.hc3n_band())
        source = "band model (Villanueva et al. 2013 JQSRT nu1)"
    else:
        df = G.load_hitran_csv(hitran_cache(species), species)
        source = "HITRAN 2020 main isotopologue via astroquery"
    df, meta = G.prepare(df, species)
    meta["source"] = source
    return df, meta


def lam_grid():
    n = int(np.log(LAM_MAX / LAM_MIN) * R_GRID) + 1
    return LAM_MIN * np.exp(np.arange(n) / R_GRID)


def wq(x, w, q):
    """Weighted quantile."""
    o = np.argsort(x); cw = np.cumsum(w[o]); return float(np.interp(q * cw[-1], cw, x[o]))


def co_swings_table(sp, v_max=60.0, dv=0.25, temps=(30, 70, 130)):
    """CO v(1-0) g-factor versus heliocentric velocity on a fine grid.

    CO is the only species whose band g-factor depends on v_h at the > 1 % level (its lines coincide
    with the solar CO fundamental absorption lines at v_h = 0), so the pipeline scales g(CO) with
    ``ratio_T070K`` (g(v_h) / g(0) at 70 K).  Positive v_h = receding from the Sun.
    """
    df, meta = load_species("CO")
    inrange = ((df["nu"] >= 1e4 / LAM_MAX) & (df["nu"] <= 1e4 / LAM_MIN)).to_numpy()
    fund = ((df["vup"] == "1") & (df["vlow"] == "0")).to_numpy()
    v = np.round(np.arange(-v_max, v_max + dv / 2, dv), 3)
    out = {"v_h_kms": v}
    for T in temps:
        g10, gtot = [], []
        for vh in v:
            g, _, _ = G.gfactors(df, meta, float(T), sp, v_h_kms=float(vh))
            g10.append(g[fund].sum()); gtot.append(g[inrange].sum())
        out[f"g_1-0_T{T:03d}K"] = np.array(g10)
        if T == T_REF:
            out[f"g_total_0.7-5um_T{T:03d}K"] = np.array(gtot)
    d = pd.DataFrame(out)
    i0 = int(np.argmin(np.abs(v)))
    for T in temps:
        d[f"ratio_T{T:03d}K"] = d[f"g_1-0_T{T:03d}K"] / d[f"g_1-0_T{T:03d}K"].iloc[i0]
    d.to_csv(os.path.join(OUT, "co_swings.csv"), index=False, float_format="%.6g")
    print(f"CO Swings table: g(1-0, 70 K) = {d['g_1-0_T070K'].iloc[i0]:.4e} at v=0, "
          f"{d['g_1-0_T070K'].max():.4e} max; ratio at +-20 km/s = "
          f"{d['ratio_T070K'].iloc[int(np.argmin(np.abs(v-20)))]:.3f} / {d['ratio_T070K'].iloc[int(np.argmin(np.abs(v+20)))]:.3f}")
    return d


def main(species_list, quick=False):
    t0 = time.time()
    os.makedirs(os.path.join(OUT, "profiles"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "lines"), exist_ok=True)
    sp = SP.SolarPump(INPUTS)
    grid = lam_grid()
    nu_grid = 1e4 / grid
    profiles_70 = pd.DataFrame({"lam_um": grid})
    band_rows, summary_rows, swings_rows = [], [], []
    temps = [T_REF] if quick else TEMPS
    for s in species_list:
        ts = time.time()
        df, meta = load_species(s)
        inrange = (df["nu"] >= 1e4 / LAM_MAX) & (df["nu"] <= 1e4 / LAM_MIN)
        gcols = {}
        for T in temps:
            g, gp, Q = G.gfactors(df, meta, float(T), sp, v_h_kms=0.0)
            gcols[T] = g
        lam = 1e4 / df["nu"].to_numpy()
        # ---- profiles
        prof = pd.DataFrame({"lam_um": grid, "nu_cm": nu_grid})
        for T in temps:
            prof[f"g_lam_T{T:03d}K"] = G.bin_profile(lam[inrange], gcols[T][inrange], grid)
        prof.to_csv(os.path.join(OUT, "profiles", f"{s}_gprofile.csv"), index=False, float_format="%.5g")
        profiles_70[s] = prof[f"g_lam_T{T_REF:03d}K"]
        # ---- lines
        L = df.loc[inrange, ["nu", "elower", "a", "gp", "gpp", "vup", "vlow", "up_key", "low_key"]].copy()
        L.insert(1, "lam_um", 1e4 / L["nu"])
        for T in temps:
            L[f"g_T{T:03d}K"] = gcols[T][inrange.to_numpy()]
        keep = L[[f"g_T{T:03d}K" for T in temps]].max(axis=1) >= G_LINE_MIN
        L = L[keep].sort_values("nu")
        L["upper"] = L["up_key"].str.split("|").str[1]
        L["lower"] = L["low_key"].str.split("|").str[1]
        L = L.drop(columns=["up_key", "low_key"]).rename(columns={"nu": "nu_cm", "a": "A_s-1", "elower": "E_low_cm",
                                                                    "vup": "vib_upper", "vlow": "vib_lower"})
        L.to_csv(os.path.join(OUT, "lines", f"{s}_lines.csv"), index=False, float_format="%.6g")
        # ---- bands
        d = df.loc[inrange, ["nu", "vup", "vlow"]].copy()
        for T in temps:
            d[f"g_T{T:03d}K"] = gcols[T][inrange.to_numpy()]
        gref = f"g_T{T_REF:03d}K"
        tot_range = d[gref].sum()
        for (vu, vl), b in d.groupby(["vup", "vlow"]):
            if b[gref].sum() < G_BAND_MIN:
                continue
            w = b[gref].to_numpy(); nu = b["nu"].to_numpy()
            key = f"{vu}-{vl}"
            row = dict(species=s, band=key, band_name=BAND_NAMES.get(s, {}).get(key, ""),
                       nu_center_cm=float(np.average(nu, weights=w)), lam_center_um=1e4 / float(np.average(nu, weights=w)),
                       lam_05_um=1e4 / wq(nu, w, 0.95), lam_95_um=1e4 / wq(nu, w, 0.05), n_lines=len(b),
                       frac_of_species_70K=b[gref].sum() / tot_range)
            for T in temps:
                row[f"g_T{T:03d}K"] = b[f"g_T{T:03d}K"].sum()
            band_rows.append(row)
        # ---- Swings sensitivity at 70 K
        if not quick:
            main_bands = [k for k in d.groupby(["vup", "vlow"])[gref].sum().sort_values(ascending=False).index[:6]]
            for v in VH_GRID:
                g, _, _ = G.gfactors(df, meta, float(T_REF), sp, v_h_kms=float(v))
                dd = df.loc[inrange, ["vup", "vlow"]].assign(g=g[inrange.to_numpy()])
                sums = dd.groupby(["vup", "vlow"])["g"].sum()
                for key in main_bands:
                    swings_rows.append(dict(species=s, band=f"{key[0]}-{key[1]}", v_h_kms=v, g_T070K=sums.get(key, 0.0)))
        summary_rows.append(dict(species=s, source=meta["source"], n_lines_total=len(df), n_lines_in_range=int(inrange.sum()),
                                 n_ground_levels=len(meta["ground_levels"]), ground_labels=";".join(meta["ground"]),
                                 E_vib1_cm=meta["e_vib1"], Q_296K=G.partition_function(meta, 296.0, vib=True),
                                 **{f"Q_T{T:03d}K": G.partition_function(meta, float(T), vib=True) for T in temps},
                                 **{f"g_total_0.7-5um_T{T:03d}K": float(gcols[T][inrange.to_numpy()].sum()) for T in temps},
                                 **{f"g_total_all_T{T:03d}K": float(gcols[T].sum()) for T in temps}))
        print(f"{s:6s} {len(df):7d} lines, {int(inrange.sum()):6d} in range, g(0.7-5um, 70K) = "
              f"{gcols[T_REF][inrange.to_numpy()].sum():.3e}, {time.time() - ts:.1f} s", flush=True)
    bands = pd.DataFrame(band_rows).sort_values(["species", f"g_T{T_REF:03d}K"], ascending=[True, False])
    bands.to_csv(os.path.join(OUT, "band_gfactors.csv"), index=False, float_format="%.5g")
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT, "species_summary.csv"), index=False, float_format="%.5g")
    profiles_70.to_csv(os.path.join(OUT, "profiles", f"all_species_T{T_REF:03d}K.csv"), index=False, float_format="%.5g")
    if swings_rows:
        sw = pd.DataFrame(swings_rows)
        ref = sw[sw.v_h_kms == 0].set_index(["species", "band"])["g_T070K"]
        r0 = ref.reindex(pd.MultiIndex.from_frame(sw[["species", "band"]])).to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            sw["ratio_to_v0"] = np.where(r0 > 0, sw["g_T070K"].to_numpy() / r0, np.nan)
        sw.to_csv(os.path.join(OUT, "swings_sensitivity.csv"), index=False, float_format="%.5g")
    if not quick:
        co_swings_table(sp)
    # ---- solar spectrum table (1 cm-1 bins)
    edges = np.arange(600.0, 33301.0, 1.0)
    nu_c = 0.5 * (edges[:-1] + edges[1:])
    fine = np.arange(600.0, 33300.0, 0.01)
    tr = sp.transmittance(fine)
    idx = np.clip(((fine - 600.0) / 1.0).astype(int), 0, len(nu_c) - 1)
    tr_mean = np.bincount(idx, weights=tr, minlength=len(nu_c)) / np.bincount(idx, minlength=len(nu_c))
    cont = sp.continuum_flam(nu_c)
    pd.DataFrame({"nu_cm": nu_c, "lam_um": 1e4 / nu_c, "continuum_W_m-2_um-1": cont, "transmittance_mean": tr_mean,
                  "flux_W_m-2_um-1": cont * tr_mean}).to_csv(os.path.join(OUT, "solar_flux_1au.csv"), index=False, float_format="%.6g")
    meta_out = dict(built=time.strftime("%Y-%m-%d %H:%M"), temps_K=temps, v_h_grid_kms=VH_GRID, lam_range_um=[LAM_MIN, LAM_MAX],
                    grid_resolving_power=R_GRID, g_line_min=G_LINE_MIN, g_band_min=G_BAND_MIN, species=species_list,
                    solar=dict(continuum="Kurucz ATLAS9 fsunallp.10000resam25 x 2.720e-4", lines="Toon JPL solar_merged_20240731 (disk-integrated)"))
    json.dump(meta_out, open(os.path.join(OUT, "build_meta.json"), "w"), indent=1)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="*", default=SPECIES)
    ap.add_argument("--quick", action="store_true", help="70 K only, no Swings grid")
    a = ap.parse_args()
    main(a.species, quick=a.quick)
