#!/usr/bin/env python
"""Rebuild the tables in doc/afrho_heliocentric_trends.md from results/ztf/afrho/.

The narrative in the note is written by hand; everything between the
``<!-- generated -->`` markers is produced here, so a rerun of
``afrho_trends.py`` is followed by this script and nothing is transcribed.
"""

from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd

import ztfcomet as zc

S = {"target": str}
DOC = zc.DOC_ROOT / "afrho_heliocentric_trends.md"


def cell(s):
    if s.empty or not np.isfinite(s["x"].iloc[0]):
        return "—"
    r = s.iloc[0]
    core = f"{r['x']:+.2f} ± {r['x_err_scaled']:.2f} ({int(r['n'])})"
    return f"**{core} {r['grade']}**" if r["grade"] in "AB" else f"{core} {r['grade']}"


def _num(v, sig=3):
    """*sig* significant figures without an exponent."""
    if not np.isfinite(v):
        return "—"
    if v == 0:
        return "0"
    d = max(0, sig - 1 - int(np.floor(np.log10(abs(v)))))
    return f"{v:,.{d}f}"


def _reason(note):
    """The note of a phase without an estimate, reduced to its kind."""
    if "no clean ZTF" in note:
        return "no clean ZTF r-band frames at this aperture"
    if "fits only the" in note:
        return "ZTF fitted only the other phase of the orbit"
    if "no fitted ZTF trend" in note:
        return "no fitted trend at this aperture"
    if "unconstrained" in note:
        return "the law is grade D and ⟨r_h⟩ lies outside its data"
    if "beyond the ZTF range" in note:
        return "⟨r_h⟩ more than 0.1 dex beyond the fitted range"
    return note


def spherex_section(sx):
    """Census of the Afρ estimates at the SPHEREx epochs, then one line per phase."""
    out = ["\n### Afρ at the SPHEREx epochs\n",
           "Per SPHEREx phase group, the ZTF r-band A(0°)fρ at the group's mean r_h (`spherex_afrho.csv`; "
           "`activity.afrho_at_epoch`).  `direct`: clean frames within 5 d of the window, moved to ⟨r_h⟩ along the "
           "local law and averaged; `trend`: the fitted law of the orbital phase the epoch falls in, at ⟨r_h⟩; "
           "`trend_extrap`: that law extended by at most 0.1 dex; `none`: no estimate, with the reason.  Errors are "
           "1σ in log space (the scatter about the law included for `trend`); the linear error quoted is the "
           "symmetric approximation, the 1σ range is in the table.\n",
           "| ρ (km) | phases | direct | trend | trend_extrap | none | with a value |\n|---|---|---|---|---|---|---|"]
    for rho, g in sx.groupby("rho_km"):
        c = g["method"].value_counts()
        out.append(f"| {int(rho) // 1000}k | {len(g)} | {c.get('direct', 0)} | {c.get('trend', 0)} | "
                   f"{c.get('trend_extrap', 0)} | {c.get('none', 0)} | {int(np.isfinite(g['afrho_cm']).sum())} |")
    none = sx[sx["method"] == "none"].assign(reason=lambda d: d["note"].map(_reason))
    if len(none):
        out.append("\nWhy there is no estimate:\n")
        out.append("| reason | 10,000 km | 20,000 km |\n|---|---|---|")
        tab = none.groupby(["reason", "rho_km"]).size().unstack(fill_value=0)
        for reason, r in tab.sort_values(10000 if 10000 in tab else tab.columns[0], ascending=False).iterrows():
            out.append(f"| {reason} | {r.get(10000, 0)} | {r.get(20000, 0)} |")
    out.append("\nPer phase (value ± error in cm, then how: method, points, grade of the law used):\n")
    out.append("| comet | phase | ⟨r_h⟩ (au) | T−T_p (d) | 10,000 km | how | 20,000 km | how |\n|---|---|---|---|---|---|---|---|")
    by = {rho: sx[sx["rho_km"] == rho].set_index(["target", "phase"]) for rho in (10000, 20000)}

    def cell(rho, key):
        if key not in by[rho].index:
            return "—", ""
        r = by[rho].loc[key]
        if np.isfinite(r["afrho_cm"]):
            g = f" {r['grade']}" if isinstance(r["grade"], str) and r["grade"] else ""
            return f"{_num(r['afrho_cm'])} ± {_num(r['afrho_err_cm'], 2)}", f"{r['method']} ({int(r['n'])}){g}"
        return "—", _reason(str(r["note"]))

    for key in sorted(set(by[10000].index) | set(by[20000].index), key=lambda k: (k[0], k[1])):
        r = by[10000].loc[key] if key in by[10000].index else by[20000].loc[key]
        t = f"{r['t_tp']:+.0f}" if np.isfinite(r["t_tp"]) else "—"
        v1, h1 = cell(10000, key); v2, h2 = cell(20000, key)
        out.append(f"| {key[0]} | S{key[1]} | {r['rh']:.2f} | {t} | {v1} | {h1} | {v2} | {h2} |")
    return out


def leg_order(l):
    for i, k in enumerate(("rising", "fading", "inbound", "outbound")):
        if l.startswith(k):
            return (i, 0 if " r_h" not in l else (1 if "<" in l else 2))
    return (9, 0)


def build():
    R = zc.result_dir("afrho", create=False)
    t = pd.read_csv(R / "trends.csv", dtype=S); p = pd.read_csv(R / "peaks.csv", dtype=S)
    b = pd.read_csv(R / "breaks.csv", dtype=S); c = pd.read_csv(R / "colour.csv", dtype=S)
    o = pd.read_csv(R / "outbursts.csv", dtype=S); w = pd.read_csv(R / "spherex_windows.csv", dtype=S)
    sx = pd.read_csv(R / "spherex_afrho.csv", dtype=S) if (R / "spherex_afrho.csv").exists() else pd.DataFrame()
    pr = t[t["primary"] & (t["band"] == "r")]
    out = []
    ab = pr[pr["grade"].isin(["A", "B"]) & pr["split"].isin(["peak", "perihelion"])]

    def med(leg, rho):
        v = ab[(ab["leg"] == leg) & (ab["rho_km"] == rho)]["x"]
        return f"{v.median():.2f} [{v.quantile(.16):.2f}, {v.quantile(.84):.2f}] (n = {len(v)})" if len(v) else "—"

    a = pd.concat([pd.read_csv(f, dtype=S) for f in glob.glob(str(zc.result_dir("photometry", create=False) / "*.csv"))])
    q = p[(p["band"] == "r") & (p["rho_km"] == 10000)]
    bp = b[(b["band"] == "r") & b["primary"] & b["preferred"]].sort_values("r_break")
    oo = o[o["band"] == "r"].sort_values(["target", "rho_km", "rh_start"])
    seg_t = sorted(pr[pr["split"] == "segment"]["target"].unique())
    n_tail = int((pr[pr["n_tail"] > 0].groupby("target").size() > 0).sum())
    out.append("### Headline numbers (r-band, primary phases, grades A/B)\n")
    out.append("| phase | 10,000 km: median x [16–84%] | 20,000 km |\n|---|---|---|")
    for leg, lab in (("rising", "rising (to the peak)"), ("fading", "fading (from the peak)"),
                     ("inbound", "inbound only"), ("outbound", "outbound only")):
        out.append(f"| {lab} | {med(leg, 10000)} | {med(leg, 20000)} |")
    out.append(f"\nSurvey clean rate {100 * a['quality_ok'].mean():.1f}% after the anomaly flag "
               f"({int(a['flag_anomalous_bright'].sum())} rows).  Of {len(q)} two-sided comets, "
               f"{int(q['interior'].sum())} have an interior maximum and are split there "
               f"({int(q['bracketed'].sum())} bracketed on both sides).  "
               f"{int((bp['rho_km'] == 10000).sum())} primary legs at 10,000 km prefer a broken law and are divided: "
               f"{', '.join(seg_t)}.  {int((oo['rho_km'] == 10000).sum())} outburst windows on "
               f"{oo[oo['rho_km'] == 10000]['target'].nunique()} comets are excluded from the fits; "
               f"{n_tail} comets have an isolated tail set aside.  SPHEREx: {len(w)} windows over {w['target'].nunique()} comets.\n")
    out.append("### Activity peaks\n")
    out.append("| comet | ρ (km) | n rising / fading | peak T−T_p (d) | 16–84% | bracketed |\n|---|---|---|---|---|---|")
    for _, r in p[(p["band"] == "r") & p["interior"]].sort_values(["target", "rho_km"]).iterrows():
        k = pr[(pr["target"] == r["target"]) & (pr["rho_km"] == r["rho_km"]) & (pr["split"] == "peak")]
        out.append(f"| {r['target']} | {int(r['rho_km']) // 1000}k | {int(k[k['leg'] == 'rising']['n'].sum())} / "
                   f"{int(k[k['leg'] == 'fading']['n'].sum())} | **{r['t_peak']:+.0f}** | [{r['t_lo']:+.0f}, {r['t_hi']:+.0f}] | "
                   f"{'yes' if r['bracketed'] else 'plateau on one side'} |")
    out.append("\n### Where a single power law fails\n")
    out.append("| comet | ρ (km) | phase | r_h range | break (au) [16–84%] | x inside → outside | ΔBIC |\n|---|---|---|---|---|---|---|")
    for _, r in bp.iterrows():
        out.append(f"| {r['target']} | {int(r['rho_km']) // 1000}k | {r['leg']} | {r['rh_min']:.2f}–{r['rh_max']:.2f} | "
                   f"**{r['r_break']:.2f}** [{r['r_break_lo']:.2f}, {r['r_break_hi']:.2f}] | "
                   f"{r['x_inner']:.1f} ± {r['x_inner_err']:.1f} → {r['x_outer']:.1f} ± {r['x_outer_err']:.1f} | {r['dbic']:.0f} |")
    out.append("\n#### The 3 au hypothesis\n")
    out.append("| comet | ρ (km) | phase | x (< 3 au) | x (> 3 au) | split preferred |\n|---|---|---|---|---|---|")
    for _, r in b[(b["band"] == "r") & b["primary"] & b["at3_testable"]].iterrows():
        out.append(f"| {r['target']} | {int(r['rho_km']) // 1000}k | {r['leg']} | {r['x_lt3']:+.2f} ± {r['x_lt3_err']:.2f} | "
                   f"{r['x_gt3']:+.2f} ± {r['x_gt3_err']:.2f} | {'yes' if r['at3_dbic'] > 6 else 'no'} (ΔBIC {r['at3_dbic']:.0f}) |")
    out.append("\n### Outbursts\n")
    out.append("| comet | ρ (km) | r_h at onset (au) | T−T_p (d) | rise (dex) | frames | decayed within coverage |\n|---|---|---|---|---|---|---|")
    for _, r in oo.iterrows():
        out.append(f"| {r['target']} | {int(r['rho_km']) // 1000}k | {r['rh_start']:.2f} | {r['t_start']:+.0f} … {r['t_end']:+.0f} | "
                   f"+{r['rise_dex']:.2f} | {int(r['n_points'])} | {'yes' if r['ended'] else 'no'} |")
    out.append("\n### Dust colour\n")
    out.append("| comet | ρ (km) | pairs | r_h (au) | median excess (mag) | slope (mag/dex) | change at r_h | Δ colour |\n|---|---|---|---|---|---|---|---|")
    for _, r in c[c["leg"] == "all"].sort_values(["target", "rho_km"]).iterrows():
        ch = f"**{r['r_change']:.2f}** [{r['r_change_lo']:.2f}, {r['r_change_hi']:.2f}]" if r["change_preferred"] else "none"
        dc = f"**{r['delta_colour']:+.2f} ± {r['delta_colour_err']:.2f}**" if r["change_preferred"] else "—"
        out.append(f"| {r['target']} | {int(r['rho_km']) // 1000}k | {int(r['n_pairs'])} | {r['rh_min']:.2f}–{r['rh_max']:.2f} | "
                   f"{r['excess_median']:+.3f} | {r['slope']:+.2f} ± {r['slope_err']:.2f} | {ch} | {dc} |")
    out.append("\n### Per-comet indices, r-band, primary phases and segments\n\nScaled error, (N), grade; A/B in bold.  "
               "Segments (`r_h<` / `r_h>`) are the divided trend where a break is preferred.\n")
    out.append("| comet | r_h (au) | phase | 10,000 km | 20,000 km |\n|---|---|---|---|---|")
    for tg, s in pr.groupby("target"):
        for leg in sorted(s["leg"].unique(), key=leg_order):
            out.append(f"| {tg} | {s['rh_min'].min():.2f}–{s['rh_max'].max():.2f} | {leg} | "
                       f"{cell(s[(s['rho_km'] == 10000) & (s['leg'] == leg)])} | {cell(s[(s['rho_km'] == 20000) & (s['leg'] == leg)])} |")
    if len(sx):
        out += spherex_section(sx)
    return "\n".join(out) + "\n"


def main():
    doc = DOC.read_text()
    block = "<!-- generated -->\n" + build() + "<!-- /generated -->"
    if "<!-- generated -->" in doc:
        doc = re.sub(r"<!-- generated -->.*?<!-- /generated -->", lambda m: block, doc, flags=re.S)
    else:
        doc = doc.rstrip() + "\n\n## Results (generated)\n\n" + block + "\n"
    DOC.write_text(doc)
    print("tables regenerated in", DOC)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
