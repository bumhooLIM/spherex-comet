"""
Build a skim deck of the per-phase figures (``doc/figures.pptx``).

Three slides per (target, phase) group of ``results/comspec/gas_fit.csv``:

  1/3  the spectra of that group at its adopted aperture --
         raw spectrum          ``fig/comspec/cont_subtract/<stem>_raw.png``
         continuum fit         ``fig/comspec/cont_subtract/<stem>_validation.png``
         emission model fit    ``fig/comspec/emission_model/<stem>.png``
  2/3  the ZTF Afρ trend of the comet with this phase highlighted
         ``fig/ztf/afrho/trend_slide/<target>_S<phase>.png``
         (``scripts/ztf/afrho_trends.py``; the r_h panels at 10 000 and 20 000 km
         and the per-phase table; absent for comets outside the ZTF survey)
  3/3  the stacked images at the phase: the ZTF r-band frames nearest the SPHEREx
       window and the SPHEREx exposures of the phase in the dust continuum
       (1.2-2.5 um) and the H2O, CO2 and CO emission windows
         ``fig/comspec/phase_images/<stem>.png``  (``scripts/comspec/phase_images.py``)

with ``<stem> = <target>_<aperture_label>km_ph<phase>``.  Every slide carries the same
two header lines -- aperture, r_h, channel count; the three production rates with
their tiers and chi2_nu; and the ZTF A(0°)fρ at ⟨r_h⟩ for both ZTF apertures with how
it was obtained (``afrho_*`` columns of ``gas_fit.csv``) -- so the deck can be paged
through without the tables open.  The raster-heavy phase images are embedded as JPEG
copies to keep the deck within reason; the files on disk stay PNG.

Before/after comparison (2026-09-16).  When a *previous* figure set exists (by default
``fig/comspec/previous_rev260915/``, the figures and ``gas_fit.csv`` of the groups the review
memo of 2026-09-15 revised, copied before the rerun), every group that has figures there gets
its three *previous* slides (header "BEFORE the 2026-09-15 revision") immediately before the
three current ones ("AFTER"), so the two states can be paged through side by side.  A regrouped
comet maps its new phases onto the old ones through ``PREVIOUS_PHASES``.  ``--no-previous``
builds the plain deck.

Usage
-----
    python scripts/comspec/make_figure_slides.py                 # every group
    python scripts/comspec/make_figure_slides.py 2P:1            # one group (example slide)
    python scripts/comspec/make_figure_slides.py 2P 10P:2        # a target, or single groups
    python scripts/comspec/make_figure_slides.py --no-previous   # without the before/after pairs
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "fig" / "comspec"
RESULTS = ROOT / "results" / "comspec"
OUT = ROOT / "doc" / "figures.pptx"
CACHE = Path()          # set in main(): cropped copies of the figures
TREND_FIG = Path(os.environ.get("ZTF_TREND_SLIDES", ROOT / "fig" / "ztf" / "afrho" / "trend_slide")).expanduser()
IMAGE_FIG = FIG / "phase_images"
#: the figure set of the state before the 2026-09-15 case revisions (see the module docstring)
PREVIOUS = Path(os.environ.get("COMSPEC_PREVIOUS_FIGS", FIG / "previous_rev260915")).expanduser()
#: new phase -> the old phases it replaced, for comets regrouped by the revision
PREVIOUS_PHASES = {("240P", 2): [2, 3], ("240P", 3): [4]}

SLIDE_W, SLIDE_H = 13.333, 7.5          # 16:9 inches

# --- palette -----------------------------------------------------------------
NAVY = RGBColor(0x1B, 0x2A, 0x41)
BAND = RGBColor(0xEE, 0xF1, 0xF5)
INK = RGBColor(0x1B, 0x2A, 0x41)
MUTED = RGBColor(0x5A, 0x66, 0x76)
GREEN = RGBColor(0x1B, 0x7F, 0x3B)
AMBER = RGBColor(0xB5, 0x6B, 0x00)
RED = RGBColor(0xB3, 0x2D, 0x2D)

TIER_COLOR = {"detected": GREEN, "marginal": AMBER}
AFRHO_COLOR = {"direct": GREEN, "trend": AMBER, "trend_extrap": AMBER}

SPECIES = [("H2O", "H₂O"), ("CO2", "CO₂"), ("CO", "CO")]


# --- helpers -----------------------------------------------------------------
def aperture_label(r_km: float) -> str:
    """``20000 -> "2e4"`` -- the file-name convention of ``spherex_comspec``."""
    e = int(np.floor(np.log10(float(r_km))))
    return f"{float(r_km) / 10 ** e:g}e{e}"


def _sci(x: float, nd: int = 2) -> str:
    """``9.31e+25 -> "9.3e25"`` (compact, for a one-line header)."""
    if not np.isfinite(x):
        return "--"
    e = int(np.floor(np.log10(abs(x))))
    return f"{x / 10 ** e:.{nd - 1}f}e{e}"


def q_phrase(row: pd.Series, key: str, pretty: str) -> tuple[str, RGBColor]:
    """One species' result as it should read on the header line."""
    status = str(row.get(f"Q_{key}_status", ""))
    if status == "detected" or status == "marginal":
        val, err = row[f"Q_{key}"], row[f"Q_{key}_err"]
        nsig = row.get(f"Q_{key}_nsig", np.nan)
        tail = f" ({nsig:.1f}σ)" if np.isfinite(nsig) else ""
        tag = "" if status == "detected" else " marginal"
        return f"Q({pretty}) = {_sci(val)} ± {_sci(err)}{tail}{tag}", TIER_COLOR[status]
    lim = row.get(f"Q_{key}_upper_limit", np.nan)
    if np.isfinite(lim):
        return f"Q({pretty}) < {_sci(lim)}", MUTED
    return f"Q({pretty}) {status.replace('_', ' ')}", MUTED


def afrho_reason(note: str, tag: str) -> str:
    """The ZTF-side reason a phase has no Afρ, compacted from ``afrho_note``."""
    if not isinstance(note, str) or not note.strip():
        return "not in the ZTF survey"
    part = ""
    for chunk in str(note).split(";"):
        if chunk.strip().startswith(f"{tag}:"):
            part = chunk.split(":", 1)[1].strip()
    if "no clean ZTF" in part:
        return "no clean ZTF frames at this aperture"
    if "fits only the" in part:
        return "ZTF covers only the other side of perihelion"
    if "no fitted ZTF trend" in part:
        return "no fitted ZTF trend"
    if "unconstrained" in part:
        return "grade-D law, ⟨r_h⟩ outside its data"
    if "beyond the ZTF range" in part:
        return "⟨r_h⟩ beyond the fitted range"
    if part.startswith("stale"):
        return "stale (phase regrouped)"
    return part[:40] if part else "no ZTF estimate"


def afrho_phrase(row: pd.Series, tag: str, label: str) -> tuple[str, RGBColor]:
    """One aperture's ZTF Afρ as it should read on the header line."""
    val = row.get(f"afrho_{tag}_cm", np.nan)
    if np.isfinite(val):
        err, method = row.get(f"afrho_{tag}_err_cm", np.nan), str(row.get(f"afrho_{tag}_method", ""))
        return (f"{label}: {_sig(val)} ± {_sig(err, 2)} cm ({method.replace('_', ' ')})",
                AFRHO_COLOR.get(method, MUTED))
    return f"{label}: — ({afrho_reason(row.get('afrho_note', ''), tag)})", MUTED


def _sig(x: float, sig: int = 3) -> str:
    """*sig* significant figures without an exponent: 4,030 / 150 / 8.15 / 0.69."""
    if not np.isfinite(x):
        return "--"
    if x == 0:
        return "0"
    d = max(0, sig - 1 - int(np.floor(np.log10(abs(x)))))
    return f"{x:,.{d}f}"


def textbox(slide, x, y, w, h, text, *, size=12, bold=False, color=INK,
            align=None, font="Arial"):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    if align is not None:
        p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = color
    run.font.name = font
    return tb


def rich_line(slide, x, y, w, h, parts, *, size=11, font="Arial"):
    """A single line built from ``(text, color, bold)`` runs."""
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    for text, color, bold in parts:
        run = p.add_run()
        run.text = text
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = color
        run.font.name = font
    return tb


def rect(slide, x, y, w, h, color):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Inches(h))
    sh.fill.solid()
    sh.fill.fore_color.rgb = color
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def _ink_rows(img: Image.Image, thresh: int = 200) -> np.ndarray:
    """Boolean per row: does this row hold any ink (dark pixel)?"""
    a = np.asarray(img.convert("L"))
    return (a < thresh).any(axis=1)


def _blocks(ink: np.ndarray, upto: int) -> list[tuple[int, int]]:
    """``[(start, end), ...]`` of the runs of inked rows in ``ink[:upto]``."""
    out, start = [], None
    for i in range(upto):
        if ink[i] and start is None:
            start = i
        elif not ink[i] and start is not None:
            out.append((start, i))
            start = None
    if start is not None:
        out.append((start, upto))
    return out


def strip_header(path: Path, n_lines: int, cache: Path) -> Path:
    """
    Drop the ``n_lines`` figure-level header lines, which repeat the slide header.

    The figures are saved with ``bbox_inches="tight"``, so each header line is its own
    run of inked rows at the top of the image and the cut is the start of the run that
    follows them -- the per-panel titles, which must survive.  The line count is a
    property of the generating function (``plot_fit`` always writes suptitle + Q line +
    chi2 line; ``plot_validation_grid`` writes a suptitle), not of the data.  A cut that
    would reach past 30 % of the height, or leave nothing behind, is abandoned and the
    figure used whole: a duplicated title costs a little space, a clipped panel costs
    information.
    """
    if n_lines <= 0:
        return path
    out = cache / path.name
    if out.exists():
        return out
    img = Image.open(path)
    ink = (np.asarray(img.convert("L")) < 200).any(axis=1)
    limit = int(0.30 * len(ink))
    blk = _blocks(ink, limit)
    if len(blk) <= n_lines:
        return path                                   # header and panels run together
    cut = max(0, blk[n_lines][0] - 2)
    if not 0 < cut <= limit:
        return path
    img.crop((0, cut, img.width, img.height)).save(out)
    return out


def as_jpeg(path: Path, cache: Path, quality: int = 88) -> Path:
    """A JPEG copy of a raster-heavy PNG (stacked images compress poorly as PNG)."""
    out = cache / (path.stem + ".jpg")
    if not out.exists():
        Image.open(path).convert("RGB").save(out, "JPEG", quality=quality, optimize=True)
    return out


def place(slide, path: Path, x, y, max_w, max_h, label, *, header_lines=0, jpeg=False,
          missing="not produced"):
    """Drop a figure into the box (x, y, max_w, max_h), preserving its aspect."""
    if label:
        textbox(slide, x, y, max_w, 0.24, label, size=12, bold=True, color=NAVY)
        y += 0.28
    if not path.exists():
        textbox(slide, x, y + 0.2, max_w, 0.3, f"[{path.name}: {missing}]",
                size=11, color=MUTED)
        return
    path = strip_header(path, header_lines, CACHE)
    if jpeg:
        path = as_jpeg(path, CACHE)
    iw, ih = Image.open(path).size
    scale = min(max_w / iw, max_h / ih)
    w, h = iw * scale, ih * scale
    slide.shapes.add_picture(str(path), Inches(x + (max_w - w) / 2), Inches(y),
                             width=Inches(w), height=Inches(h))


# --- the three slides of a group ----------------------------------------------
RESULT_TOP, RESULT_H = 0.80, 0.86            # the two-line result band
BODY_TOP, BODY_BOT = 1.75, 7.28


def _header(slide, row: pd.Series, kind: str, stem: str):
    """The navy title band, the two-line result band and the footer, shared by the
    three slides of a group so each one reads on its own."""
    target, phase, r_ap = row["target"], int(row["phase"]), float(row["r_ap_km"])
    rect(slide, 0, 0, SLIDE_W, 0.80, NAVY)
    textbox(slide, 0.35, 0.17, 8.0, 0.45,
            f"{target}   —   phase {phase}", size=26, bold=True,
            color=RGBColor(0xFF, 0xFF, 0xFF))
    textbox(slide, 6.2, 0.10, 6.75, 0.30, kind, size=12, bold=True,
            color=RGBColor(0xFF, 0xD7, 0x6E), align=3)                # 3 = right
    textbox(slide, 6.2, 0.42, 6.75, 0.30,
            f"r_ap = {r_ap:,.0f} km    ⟨r_h⟩ = {row['r_hel_mean']:.2f} au"
            f"    N = {int(row['n_points'])} ch",
            size=13, color=RGBColor(0xC8, 0xD3, 0xE0), align=3)

    rect(slide, 0, RESULT_TOP, SLIDE_W, RESULT_H, BAND)
    parts = []
    for key, pretty in SPECIES:
        txt, col = q_phrase(row, key, pretty)
        parts.append((txt, col, True))
        parts.append(("     ·     ", MUTED, False))
    parts.append((f"χ²_ν = {row['chi2_red']:.2f}", INK, False))
    parts.append((f"     bands: {row['bands_used']}", MUTED, False))
    rich_line(slide, 0.35, RESULT_TOP + 0.10, SLIDE_W - 0.7, 0.3, parts, size=11.5)
    parts = [("ZTF A(0°)fρ at ⟨r_h⟩:   ", INK, False)]
    for tag, label in (("10k", "10,000 km"), ("20k", "20,000 km")):
        txt, col = afrho_phrase(row, tag, label)
        parts.append((txt, col, True))
        parts.append(("     ·     ", MUTED, False))
    parts.pop()
    rich_line(slide, 0.35, RESULT_TOP + 0.46, SLIDE_W - 0.7, 0.3, parts, size=11.5)
    textbox(slide, 0.35, SLIDE_H - 0.30, 6.0, 0.22, stem, size=9, color=MUTED)


def _stem(row: pd.Series) -> str:
    return f"{row['target']}_{aperture_label(float(row['r_ap_km']))}km_ph{int(row['phase'])}"


class FigSet:
    """Where the three figures of a group live: the current set, or a previous snapshot."""

    def __init__(self, cont: Path, emission: Path, images: Path, trend: Path, tag: str = ""):
        self.cont, self.emission, self.images, self.trend, self.tag = cont, emission, images, trend, tag

    def has(self, row: pd.Series) -> bool:
        stem = _stem(row)
        return (self.cont / f"{stem}_validation.png").is_file() or (self.emission / f"{stem}.png").is_file()


CURRENT = FigSet(FIG / "cont_subtract", FIG / "emission_model", IMAGE_FIG, TREND_FIG)


def previous_set(root: Path) -> FigSet:
    return FigSet(root / "cont_subtract", root / "emission_model", root / "phase_images", root / "trend_slide",
                  tag="BEFORE the 2026-09-15 revision   ·   ")


def add_slide(prs, row: pd.Series, figs: FigSet = CURRENT):
    """1/3 -- the spectra: raw, continuum fit, emission-model fit."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])       # blank
    stem = _stem(row)
    _header(slide, row, figs.tag + "1 / 3   spectra  ·  continuum  ·  emission model", stem)
    top, bot = BODY_TOP, BODY_BOT
    body_h = bot - top
    left_w = 5.35
    place(slide, figs.cont / f"{stem}_validation.png",
          0.35, top, left_w, body_h - 0.28,
          "②  continuum fit  ·  subtraction  ·  residuals",
          header_lines=1)
    right_x, right_w = 6.45, 6.13
    place(slide, figs.cont / f"{stem}_raw.png",
          right_x, top, right_w, 2.00, "①  raw spectrum")
    place(slide, figs.emission / f"{stem}.png",
          right_x, top + 2.45, right_w, 3.05, "③  emission model fit",
          header_lines=3)
    return slide


def add_trend_slide(prs, row: pd.Series, figs: FigSet = CURRENT):
    """2/3 -- the ZTF Afρ trend of the comet, this phase highlighted."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    stem = _stem(row)
    _header(slide, row, figs.tag + "2 / 3   ZTF Afρ trend  ·  this phase highlighted", stem)
    path = figs.trend / f"{row['target']}_S{int(row['phase'])}.png"
    place(slide, path, 0.35, BODY_TOP, SLIDE_W - 0.7, BODY_BOT - BODY_TOP - 0.28,
          "ZTF A(0°)fρ against r_h at ρ = 10,000 and 20,000 km; red: the value at each SPHEREx phase; "
          "gold ring and tinted row: this phase",
          missing="comet not in the ZTF survey, no trend figure")
    return slide


def add_images_slide(prs, row: pd.Series, figs: FigSet = CURRENT):
    """3/3 -- the stacked images at the phase."""
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    stem = _stem(row)
    _header(slide, row, figs.tag + "3 / 3   stacked images at the phase", stem)
    place(slide, figs.images / f"{stem}.png", 0.35, BODY_TOP, SLIDE_W - 0.7, BODY_BOT - BODY_TOP - 0.28,
          "ZTF r nearest the SPHEREx window  ·  SPHEREx exposures of the phase in the continuum (1.2–2.5 µm), "
          "H₂O, CO₂ and CO windows  ·  north up, comet-centred medians",
          jpeg=True, missing="run scripts/phase_images.py")
    return slide


def previous_rows(row: pd.Series, prev_fits: pd.DataFrame) -> list:
    """The previous-state rows of this group: the same (target, phase), or the old phases a
    regrouped phase replaced (``PREVIOUS_PHASES``)."""
    t, ph = row["target"], int(row["phase"])
    phases = PREVIOUS_PHASES.get((t, ph), [ph])
    out = []
    for old in phases:
        r = prev_fits[(prev_fits.target == t) & (prev_fits.phase.astype(int) == int(old))]
        if len(r):
            out.append(r.iloc[0])
    return out


# --- driver ------------------------------------------------------------------
def select(fits: pd.DataFrame, args: list[str]) -> pd.DataFrame:
    if not args:
        return fits
    keep = pd.Series(False, index=fits.index)
    for a in args:
        if ":" in a:
            t, ph = a.split(":")
            keep |= (fits.target == t) & (fits.phase.astype(int) == int(ph))
        else:
            keep |= fits.target == a
    return fits[keep]


def sort_key(target: str):
    """``2P`` < ``10P`` < ``508P`` < ``2019U5`` -- numbered comets first, by number."""
    m = re.fullmatch(r"(\d+)P", target)
    return (0, int(m.group(1)), "") if m else (1, 0, target)


def main(args: list[str]) -> None:
    global CACHE
    use_previous = "--no-previous" not in args
    args = [a for a in args if not a.startswith("--")]
    fits = pd.read_csv(RESULTS / "gas_fit.csv", dtype={"target": str})
    sel = select(fits, args).copy()
    sel["_k"] = sel.target.map(sort_key)
    sel = sel.sort_values(["_k", "phase"])
    if sel.empty:
        raise SystemExit(f"no groups matched {args}")
    prev, prev_fits = None, None
    if use_previous and (PREVIOUS / "gas_fit.csv").is_file():
        prev = previous_set(PREVIOUS)
        prev_fits = pd.read_csv(PREVIOUS / "gas_fit.csv", dtype={"target": str})

    tmp = TemporaryDirectory(prefix="figslides_")
    CACHE = Path(tmp.name)

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(SLIDE_W), Inches(SLIDE_H)
    n_pairs = 0
    for _, row in sel.iterrows():
        if prev is not None:
            for old in previous_rows(row, prev_fits):
                if prev.has(old):
                    add_slide(prs, old, prev)
                    add_trend_slide(prs, old, prev)
                    add_images_slide(prs, old, prev)
                    n_pairs += 1
        add_slide(prs, row)
        add_trend_slide(prs, row)
        add_images_slide(prs, row)

    OUT.parent.mkdir(exist_ok=True)
    prs.save(OUT)
    tmp.cleanup()
    mb = OUT.stat().st_size / 1024 ** 2
    print(f"{OUT.relative_to(ROOT)}: {len(sel)} group(s), {len(prs.slides)} slides"
          + (f" ({n_pairs} previous-state group(s) placed before their revised slides)" if n_pairs else "")
          + f", {mb:.1f} MB")


if __name__ == "__main__":
    main(sys.argv[1:])
