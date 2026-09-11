#!/usr/bin/env python
"""Rebuild the tables in doc/afrho_heliocentric_trends.md from results/afrho/.

The narrative in the note is written by hand; everything between the
``<!-- generated -->`` markers is produced here, so a rerun of
``afrho_trends.py`` is followed by this script and nothing is transcribed.
"""

from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

import ztfcomet as zc

S = {"target": str}
DOC = zc.PROJECT_ROOT / "doc" / "afrho_heliocentric_trends.md"


def cell(s):
    if s.empty or not np.isfinite(s["x"].iloc[0]):
        return "—"
    r = s.iloc[0]
    core = f"{r['x']:+.2f} ± {r['x_err_scaled']:.2f} ({int(r['n'])})"
    return f"**{core} {r['grade']}**" if r["grade"] in "AB" else f"{core} {r['grade']}"


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
