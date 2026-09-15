"""Compare reconstructed line g-factors with the values published by Villanueva et al.
(validation_targets.csv) and write data/fluorescence/validation.csv."""
import os, re, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # -> project
sys.path.insert(0, HERE)
import gfm_core as G, gfm_bandmodels as BM, solar_pump as SP
from build_fluorescence_db import load_species, INPUTS, OUT

V_BOATTINI = 9.8     # km/s, heliocentric velocity of C/2007 W1 on 2008 July 9-10


def main():
    sp = SP.SolarPump(INPUTS)
    tg = pd.read_csv(os.path.join(HERE, "validation_targets.csv"))
    tg["g_pub"] = tg["g_published"] * tg["unit_scale"]
    out = []
    for species, grp in tg.groupby("species", sort=False):
        if species == "CH3OH":
            df, meta = G.prepare(BM.as_gfm_frame(BM.ch3oh_band("nu3")), species)
        else:
            df, meta = load_species(species)
        for T, g2 in grp.groupby("T_K"):
            vh = 0.0 if g2["source"].str.contains("ApJ|2013").any() else V_BOATTINI
            res = {v: G.gfactors(df, meta, float(T), sp, v_h_kms=v)[0] for v in sorted({vh, 0.0})}
            for r in g2.itertuples():
                vals = {}
                for v, g in res.items():
                    if r.line == "band":
                        m = np.ones(len(df), bool)
                    elif species == "CH3OH":
                        m = ((df.nu >= 2843.5) & (df.nu <= 2845.0)).to_numpy()
                    elif species == "C2H6":
                        m = (np.abs(df.nu - r.nu_cm) <= 1.0).to_numpy()
                    elif species == "H2O":
                        mm = re.match(r"\((\d)(\d)(\d)\)(\d+)-\((\d)(\d)(\d)\)(\d+)", r.line)
                        vup = " ".join(mm.group(1, 2, 3)); vlow = " ".join(mm.group(5, 6, 7))
                        m = ((np.abs(df.nu - r.nu_cm) < 0.15) & (df.vup == vup) & (df.vlow == vlow)).to_numpy()
                    else:
                        m = (np.abs(df.nu - r.nu_cm) < 0.03).to_numpy()
                        if species == "CH4" and "(" in r.line:           # e.g. R2(F1): upper J=3, symmetry F1
                            jup = int(r.line[1]) + (1 if r.line[0] == "R" else -1)
                            sym = r.line[r.line.index("(") + 1:r.line.index(")")]
                            loc = df.up_key.str.split("|").str[1].fillna("")
                            m &= loc.str.startswith(f"{jup}{sym} ").fillna(False).to_numpy(dtype=bool)
                    vals[v] = g[m].sum() if m.any() else np.nan
                out.append(dict(species=species, band=r.band, line=r.line, nu_cm=r.nu_cm, T_K=T, source=r.source,
                                g_published=r.g_pub, g_model=vals.get(vh), g_model_vh0=vals.get(0.0), v_h_kms=vh))
    res = pd.DataFrame(out)
    res["ratio"] = res.g_model / res.g_published
    res["ratio_vh0"] = res.g_model_vh0 / res.g_published
    res.to_csv(os.path.join(OUT, "validation.csv"), index=False, float_format="%.4g")
    print(res.groupby("species")["ratio"].agg(["median", "mean", "std", "count"]).to_string(float_format=lambda x: f"{x:.3f}"))
    return res


if __name__ == "__main__":
    main()
