"""
Continuum / fitting method matrix on the fixed per-phase apertures (2026-09-14).

Every run is the main variant (``BASELINE_FLAGS``, distance-corrected continuum, fixed
apertures) with one combination of continuum windows, order rule, one-sided extension and
error model, executed in an isolated result tree (``results/studies/method_matrix/runs/<name>``)
so the catalog products are untouched.  The census of each run -- ``robust`` (>= 3 sigma,
n_eff >= 2), ``clean`` (robust *and* the carrier band's continuum verdict PASS), the
negative tail (fits at <= -2 and <= -3 sigma, the empirical false-positive proxy) -- is
appended to ``matrix.jsonl`` and tabulated in ``matrix.csv``.

    python notebooks/method_matrix.py run NAME '{"windows": {...}, "continuum": {...}, "fit": {...}}'
    python notebooks/method_matrix.py batch RUNS.txt      # NAME<TAB>JSON per line, four at a time
    python notebooks/method_matrix.py table               # matrix.jsonl -> matrix.csv (+ stdout)

``windows`` overrides ``config.BAND_WINDOWS`` per band (``em`` / ``cont`` pairs) and moves
``KEY_RANGES["H2O"]["lo"]`` to the 2.7 um emission edge; ``continuum`` / ``fit`` are keyword
arguments of ``ContinuumConfig`` / ``FitConfig``.  The runs that decided the 2026-09-14
configuration are listed in ``doc/pipeline_decisions.md`` section 7.9.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "comspec" / "studies" / "method_matrix"
PY = sys.executable
CARRIER = {"H2O": "2.7um", "CO2": "4.3um", "CO": "4.7um"}


def run(name: str, opts: dict) -> dict:
    """One experiment in an isolated tree; returns (and appends) its census."""
    exp = OUT / "runs" / name
    os.environ["COMSPEC_RESULT_DIR"] = str(exp)
    os.environ["COMSPEC_DATA_DIR"] = str(exp / "data")
    os.environ["COMSPEC_FIG_DIR"] = str(exp / "fig")
    sys.path.insert(0, str(ROOT))
    import shutil
    (exp / "data").mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "data" / "comspec" / "phase_assignment.csv", exp / "data" / "phase_assignment.csv")
    import logging
    logging.disable(logging.WARNING)
    import numpy as np
    import pandas as pd
    from spherex_comspec import config as C
    from spherex_comspec.config import (BASELINE_FLAGS, ApertureConfig, ContinuumConfig, FitConfig,
                                        Variant)
    from spherex_comspec.pipeline import run_variant

    for b, w in opts.get("windows", {}).items():
        C.BAND_WINDOWS[b].update({k: tuple(v) for k, v in w.items()})
    if opts.get("windows"):
        C.KEY_RANGES["H2O"]["lo"] = C.BAND_WINDOWS["2.7um"]["em"][0]
        C.EMISSION_WINDOWS.update({b: w["em"] for b, w in C.BAND_WINDOWS.items()})
        C.ALL_EM_WINDOWS[:] = [w["em"] for w in C.BAND_WINDOWS.values()]
    fit_kw = dict(opts.get("fit", {}))
    if "accept_verdicts" in fit_kw:
        fit_kw["accept_verdicts"] = tuple(fit_kw["accept_verdicts"])
    v = Variant("dc_main", BASELINE_FLAGS, use_distcorr=True,
                continuum=ContinuumConfig(**opts.get("continuum", {})), fit=FitConfig(**fit_kw),
                aperture=ApertureConfig(**opts.get("aperture", {})), role="main")
    res = run_variant(v, progress_every=0)
    f, c = res.fits, res.cont_summary
    out = dict(name=name, opts=opts, n_fits=int(len(f)), chi2_med=float(f.chi2_red.median()),
               n_skipped=int(len(res.skipped)))
    for b in CARRIER.values():
        vc = c[c.band == b].verdict.value_counts()
        out.update({f"pass_{b}": int(vc.get("PASS", 0)), f"warn_{b}": int(vc.get("WARN", 0)),
                    f"fail_{b}": int(vc.get("FAIL", 0))})
    any_robust = np.zeros(len(f), bool)
    any_clean = np.zeros(len(f), bool)
    for s, band0 in CARRIER.items():
        src = f["h2o_source"] if "h2o_source" in f else pd.Series("main", index=f.index)
        band = np.where((s == "H2O") & (src == "hot"), "4.7um", band0)
        key = f[["target", "phase", "r_ap_km"]].assign(band=band)
        verdict = key.merge(c[["target", "phase", "r_ap_km", "band", "verdict"]],
                            on=["target", "phase", "r_ap_km", "band"], how="left").verdict.to_numpy()
        robust = ((f[f"Q_{s}_status"] == "detected") & (f[f"Q_{s}_n_eff"] >= 2)).to_numpy()
        clean = robust & (verdict == "PASS")
        cov = f[f"Q_{s}_covered"].astype(bool).to_numpy()
        nsig = f[f"Q_{s}_nsig"].to_numpy()
        out.update({f"robust_{s}": int(robust.sum()), f"clean_{s}": int(clean.sum()),
                    f"marginal_{s}": int((f[f"Q_{s}_status"] == "marginal").sum()),
                    f"covered_{s}": int(cov.sum()),
                    f"neg2_{s}": int((cov & (nsig <= -2)).sum()), f"neg3_{s}": int((cov & (nsig <= -3)).sum()),
                    f"comets_robust_{s}": int(f.loc[robust, "target"].nunique())})
        any_robust |= robust
        any_clean |= clean
    out.update(robust_total=sum(out[f"robust_{s}"] for s in CARRIER),
               clean_total=sum(out[f"clean_{s}"] for s in CARRIER),
               neg3_total=sum(out[f"neg3_{s}"] for s in CARRIER),
               comets_any_robust=int(f.loc[any_robust, "target"].nunique()),
               comets_any_clean=int(f.loc[any_clean, "target"].nunique()))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "matrix.jsonl", "a") as fh:
        fh.write(json.dumps(out) + "\n")
    print(f"{name:26s} fits {out['n_fits']:3d} chi2 {out['chi2_med']:5.2f} | robust "
          f"{out['robust_H2O']:2d}/{out['robust_CO2']:2d}/{out['robust_CO']:2d} | clean "
          f"{out['clean_H2O']:2d}/{out['clean_CO2']:2d}/{out['clean_CO']:2d} = {out['clean_total']:2d} | "
          f"<=-3sig {out['neg3_H2O']}/{out['neg3_CO2']}/{out['neg3_CO']} | comets {out['comets_any_clean']}")
    return out


def batch(path: str, workers: int = 4) -> None:
    from concurrent.futures import ThreadPoolExecutor
    runs = [line.rstrip("\n").split("\t", 1) for line in open(path) if line.strip()]

    def one(r):
        p = subprocess.run([PY, __file__, "run", r[0], r[1]], capture_output=True, text=True)
        tail = (p.stdout.strip().splitlines() or [""])[-1]
        return tail if p.returncode == 0 else f"FAILED {r[0]}: {p.stderr.strip().splitlines()[-1]}"

    with ThreadPoolExecutor(workers) as ex:
        for line in ex.map(one, runs):
            print(line, flush=True)


def table() -> None:
    import pandas as pd
    rows = [json.loads(line) for line in open(OUT / "matrix.jsonl")]
    d = pd.DataFrame(rows).drop_duplicates("name", keep="last")
    d["net"] = d.clean_total - d.neg3_total
    cols = ["name", "n_fits", "chi2_med", "pass_2.7um", "warn_2.7um", "fail_2.7um", "pass_4.3um", "fail_4.3um",
            "pass_4.7um", "fail_4.7um", "robust_H2O", "robust_CO2", "robust_CO", "robust_total",
            "clean_H2O", "clean_CO2", "clean_CO", "clean_total", "neg2_H2O", "neg2_CO2", "neg2_CO",
            "neg3_H2O", "neg3_CO2", "neg3_CO", "neg3_total", "net", "marginal_H2O", "marginal_CO2", "marginal_CO",
            "comets_any_robust", "comets_any_clean", "opts"]
    d = d[[c for c in cols if c in d]].sort_values(["clean_total", "net"], ascending=False)
    d.to_csv(OUT / "matrix.csv", index=False)
    pd.set_option("display.width", 300)
    pd.set_option("display.max_rows", 200)
    print(d.drop(columns=["opts"]).round(2).to_string(index=False))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "run":
        run(sys.argv[2], json.loads(sys.argv[3]))
    elif cmd == "batch":
        batch(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
    elif cmd == "table":
        table()
    else:
        raise SystemExit(__doc__)
