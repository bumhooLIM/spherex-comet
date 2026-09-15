"""
Build a tabulated CSV of cometary gas (CO+CO2 proxy) and dust production rates
from Gicquel, Bauer, Kramer, Mainzer & Masiero (2023), PSJ 4, 3,
"CO and CO2 Production Rates of Comets Observed by NEOWISE within Year 1 of the
Reactivated Mission"  doi:10.3847/PSJ/aca8ac

Parses Table 1 (orbital properties) and Table 2 (fluxes, Q_CO2^proxy, Afrho)
straight out of the PDF text layer and joins them on (object, visit).

IMPORTANT PHYSICS.  NEOWISE has no spectral resolution across its W2 band, so
CO2 nu3 (4.23 um) and the CO fundamental (4.67 um) fall in one 4.6 um filter and
CANNOT be separated.  The published quantity is Q_CO2^proxy: the whole W2
infrared excess converted as if it were all CO2, and it serves as a proxy for
the TOTAL CO+CO2 production.  There is therefore no independent Q(CO) in this
paper, and q_co is empty in every row by construction.  Every row is `indirect`.

Usage:  python3 scripts/build_gicquel2023_table.py
"""
import csv, os, re, subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # notebooks/comspec/ -> project
PDF  = os.path.join(ROOT, "doc", "literature", "fluorescence-emission", "Gicquel_2023_PSJ.pdf")
OUT  = os.path.join(ROOT, "data", "reference", "gicquel2023_neowise_production_rates.csv")

txt = subprocess.run(["pdftotext", "-layout", PDF, "-"],
                     capture_output=True, text=True, check=True).stdout
lines = txt.split("\n")

SPLIT = re.compile(r"\s{2,}")
NAME  = re.compile(r"^(?:\d+P|[CPX]/\d{4}\s)")          # 16P... / C/2013 R1... / P/2014 L2...

def num(s):
    """'0.7 1' -> 0.71 ; '2 518.55' -> 2518.55 ; '120/17' -> 120.17 ; '3,43' -> 3.43"""
    if s is None: return None
    s = s.replace("/", ".").replace(" ", "")
    if re.match(r"^\d{1,3},\d{3}", s):      # 56,963.0468 thousands separator
        s = s.replace(",", "")
    else:
        s = s.replace(",", ".")            # OCR decimal comma, e.g. 3,43
    if re.match(r"^\d+\.\d+\.", s):        # 56,989.194,28 -> 56989.19428
        a, b, c = s.split(".", 2)
        s = a + "." + b + c
    try: return float(s)
    except ValueError: return None

def pm(s):
    """'0.11 ± 0.03' -> (0.11, 0.03);  '(4.60 ± 0.56)x10^25' -> (4.60e25, 5.6e24)"""
    if s is None: return (None, None)
    s = s.replace("±", "±")
    exp = 0
    m = re.search(r"[×x]\s*10\s*(\d+)", s)
    if m:
        exp = int(m.group(1)); s = s[:m.start()]
    s = s.strip().strip("()")
    parts = [p for p in s.split("±")]
    v = num(parts[0]); e = num(parts[1]) if len(parts) > 1 else None
    if exp:
        v = v * 10.0**exp if v is not None else None
        e = e * 10.0**exp if e is not None else None
    return (v, e)

def clean_name(n):
    """strip the Table 2 footnote marker 'a' (Rh<2 au nucleus caution); return (name, flag)"""
    n = re.sub(r"\s+", " ", n.strip())
    if n.endswith(")a"):
        return n[:-1].strip(), True
    return n, False

def norm(n):
    """join key: designation only, ignoring the parenthetical common name"""
    return re.sub(r"\s+", "", n.split("(")[0]).upper()

# ------------------------------------------------------------------ Table 1
t1 = {}
in_t1 = False
for ln in lines:
    if "Orbital Properties of Comets Observed by NEOWISE" in ln: in_t1 = True; continue
    if in_t1 and ln.startswith("Note."): in_t1 = False; continue
    if not in_t1: continue
    f = SPLIT.split(ln.strip())
    if len(f) != 9 or not NAME.match(f[0]): continue
    name, _ = clean_name(f[0])
    t1[(norm(name), f[1])] = dict(name=name, visit=f[1], incl=num(f[2]), ecc=num(f[3]),
                                  q_peri=num(f[4]), phase=num(f[5]), cls=f[6],
                                  nobs=int(num(f[7])), mjd=num(f[8]))

# ------------------------------------------------------------------ Table 2
t2 = []
in_t2 = False
for ln in lines:
    if "Analysis of Comets Observed by NEOWISE" in ln: in_t2 = True; continue
    if in_t2 and ln.startswith("Notes."): in_t2 = False; continue
    if not in_t2: continue
    f = SPLIT.split(ln.strip())
    if len(f) != 9 or not NAME.match(f[0]): continue
    name, caution = clean_name(f[0])
    w1, w1e = pm(f[4]); w2, w2e = pm(f[5]); q, qe = pm(f[6]); af, afe = pm(f[7])
    t2.append(dict(name=name, visit=f[1], rh=num(f[2]), delta=num(f[3]),
                   w1=w1, w1e=w1e, w2=w2, w2e=w2e, q=q, qe=qe,
                   af=af, afe=afe, teff=num(f[8]), caution=caution))

# Table 1 misprints "C/2013 E2 (Iwamoto)" as "C/2012 E2"; Iwamoto is C/2013 E2.
ALIAS = {("C/2013E2", "A"): ("C/2012E2", "A")}

PREV = {"P/2014L2": "also reported by Bauer et al. 2015 (same epoch, results agree)",
        "C/2006S3": "also reported by Bauer et al. 2015, but from the WISE cryogenic prime mission, not this epoch",
        "C/2014C3": "also reported by Bauer et al. 2015 (same epoch, results agree)",
        "C/2014N3": "also reported by Bauer et al. 2015 (same epoch, results agree)",
        "C/2014Q3": "also reported by Rosser et al. 2018 (same epoch, results agree)"}

CAUTION = ("Rh < 2 au: possible nucleus contribution, use with caution (Gicquel+2023 Table 2 footnote a)")
NONUC   = ("Rh < 2 au but the signal is coma-dominated (extended coma), per Gicquel+2023 Table 2 footnote a")
EXTENDED = {"C/2012K1", "C/2012X1", "C/2013A1", "C/2014Q2", "C/2014Q3"}
PROXY = ("Q_CO2^proxy: the full W2 (4.6 um) excess converted as if it were all CO2; "
         "it is a proxy for TOTAL CO+CO2, not a CO2-only rate. NEOWISE cannot separate "
         "CO2 nu3 (4.23 um) from CO v(1-0) (4.67 um), so no independent Q(CO) exists")

COLS = ["designation", "visit", "r_hel", "delta_au", "q_peri", "ecc", "incl_deg",
        "orbital_class", "dynamical_group", "q_co2", "q_co2_err", "q_co",
        "measured_quantity", "afrho", "afrho_err", "w1_flux_mjy", "w1_flux_err_mjy",
        "w2_flux_mjy", "w2_flux_err_mjy", "t_eff_k", "phase_angle_deg", "n_obs",
        "mjd_midpoint", "reference", "instrument", "regime", "methodology", "notes"]

rows, unmatched = [], []
for r in t2:
    key = (norm(r["name"]), r["visit"])
    o = t1.get(ALIAS.get(key, key))
    if o is None:
        unmatched.append(key); o = {}
    cls = o.get("cls", "")
    grp = "SPC" if cls in ("JFC", "Centaur", "HTC") else ("LPC" if cls else "")
    notes = [PROXY]
    if r["caution"]:
        notes.append(NONUC if norm(r["name"]) in EXTENDED else CAUTION)
    if norm(r["name"]) in PREV:
        notes.append(PREV[norm(r["name"])])
    if ALIAS.get(key):
        notes.append("Gicquel+2023 Table 1 misprints this object as C/2012 E2; orbital elements taken from that row")
    if norm(r["name"]).startswith("C/2011J2-"):
        notes.append("Fragment of the split comet C/2011 J2 (LINEAR); fragments B and C are tabulated separately")
    def f(v, p="%g"): return "" if v is None else p % v
    rows.append(dict(zip(COLS, [
        r["name"], r["visit"], f(r["rh"]), f(r["delta"]), f(o.get("q_peri")),
        f(o.get("ecc")), f(o.get("incl")), cls, grp,
        f(r["q"], "%.3e"), f(r["qe"], "%.3e"), "",
        "Q_CO2^proxy (total CO+CO2)",
        f(r["af"], "%.2f"), f(r["afe"], "%.2f"),
        f(r["w1"]), f(r["w1e"]), f(r["w2"]), f(r["w2e"]),
        f(r["teff"]), f(o.get("phase")), o.get("nobs", ""), f(o.get("mjd"), "%.5f"),
        "Gicquel et al. 2023", "NEOWISE (WISE) W1 3.4 um / W2 4.6 um",
        "near-IR", "indirect", "; ".join(notes)])))

rows.sort(key=lambda x: (x["dynamical_group"], x["designation"], x["visit"]))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)

print("Table 1 rows parsed:", len(t1))
print("Table 2 rows parsed:", len(t2))
print("unmatched joins    :", unmatched)
print("written            :", len(rows), "->", os.path.relpath(OUT, ROOT))
