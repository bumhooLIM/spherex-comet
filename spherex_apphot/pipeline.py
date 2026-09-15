"""
End-to-end orchestration for one target.

:func:`run_target` is the single entry point used by ``main.py`` and by
``notebooks/apphot/main.ipynb``, so the notebook validates exactly the code the batch
run executes rather than a parallel copy of it -- which is how the primitive
notebook and script drifted apart in the first place.

No figure is produced anywhere in this module.  Plotting lives in
:mod:`spherex_apphot.plotting` and is called only from the notebook.
"""

from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from . import directory as _dir
from .apertures import annulus_radii_pix, build_apertures, pixel_scale_km
from .catalog import GaiaCatalog, GaiaSubset, NeighborIndex
from .config import Config
from .epochs import epoch_table, group_epochs
from .fitsio import FitsResolver, load_database, read_cutout
from .logging_utils import get_logger, utcnow_iso
from .masking import build_badpix_mask
from .phot import (DISTCORR_COLUMNS, PHOT_COLUMNS, add_distance_corrected_flux,
                   estimate_sky, measure_exposure_records)
from .sourceflag import project_gaia
from .status import StatusLog, slugify
from .wcsutil import north_angle_deg, wcs_from_header

__all__ = ["TargetResult", "run_target", "stack_target", "META_COLUMNS", "format_jd_columns"]

log = get_logger("pipeline")

#: Gaia sources are kept this far outside the cutout, so a star just off
#: the edge still counts towards a large aperture.
_GAIA_PAD_PIX = 30.0

#: Index columns copied verbatim onto every photometry row.
#: ``DATE-OBS`` and ``PSF_FWHM`` are renamed because ``itertuples`` mangles
#: names that are not valid Python identifiers, and because a lower-case,
#: unit-suffixed name is clearer in the output file.
META_COLUMNS: Dict[str, str] = {
    "filename": "filename", "objdesig": "objdesig", "obsid": "obsid",
    "detector": "detector", "epoch": "epoch", "epoch_n": "epoch_n",
    "wl": "wl", "wlwidth": "wlwidth", "sun_jy": "sun_jy",
    "xcen": "xcen", "ycen": "ycen", "ltv1": "ltv1", "ltv2": "ltv2",
    "cutout_size": "cutout_size", "pix_scale": "pix_scale",
    "PSF_FWHM": "psf_fwhm_arcsec",
    "ra": "ra", "dec": "dec", "r_hel": "r_hel", "r_obs": "r_obs", "alpha": "alpha",
    "hel_ecl_lon": "hel_ecl_lon", "hel_ecl_lat": "hel_ecl_lat",
    "obs_ecl_lon": "obs_ecl_lon", "obs_ecl_lat": "obs_ecl_lat",
    "racosdec_rate": "racosdec_rate", "dec_rate": "dec_rate",
    "sky_motion": "sky_motion", "sky_motion_pa": "sky_motion_pa",
    "vmag": "vmag", "jd_utc": "jd_utc", "jd_tdb": "jd_tdb",
    "DATE-OBS": "date_obs",
}


@dataclass
class TargetResult:
    """Everything one target's run produced, plus the counters the log needs."""

    objdesig: str
    slug: str
    status: str = "ok"
    phot: Optional[pd.DataFrame] = None
    epochs: Optional[pd.DataFrame] = None
    gaia: Optional[GaiaSubset] = None
    output_csv: Optional[Path] = None
    message: str = ""
    elapsed_s: float = 0.0
    counters: Dict[str, int] = field(default_factory=dict)

    def status_fields(self) -> dict:
        """Flatten into the keyword arguments :meth:`StatusLog.update` expects."""
        c = self.counters
        return dict(
            elapsed_s=round(self.elapsed_s, 2),
            n_exposures=c.get("n_exposures", 0),
            n_exposures_ok=c.get("n_exposures_ok", 0),
            n_missing_fits=c.get("n_missing_fits", 0),
            n_read_errors=c.get("n_read_errors", 0),
            n_epochs=c.get("n_epochs", 0),
            n_rows_out=c.get("n_rows_out", 0),
            n_ap_skipped=c.get("n_ap_skipped", 0),
            n_badphot=c.get("n_badphot", 0),
            n_gaia=c.get("n_gaia", 0),
            output_csv=str(self.output_csv) if self.output_csv else None,
            message=self.message,
        )


def _psf_fwhm_pix(psf_fwhm_arcsec: float, pix_scale: float) -> float:
    """PSF FWHM in pixels.  ``PSF_FWHM`` is stored in arcsec (~5-7 for SPHEREx)."""
    if not (np.isfinite(psf_fwhm_arcsec) and np.isfinite(pix_scale)) or pix_scale <= 0:
        return float("nan")
    return float(psf_fwhm_arcsec) / float(pix_scale)


def _prepare_frame(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Sort into epochs and rename the index columns that ``itertuples`` mangles."""
    work = group_epochs(df, gap_days=cfg.epoch_gap_days)
    rename = {src: dst for src, dst in META_COLUMNS.items()
              if src in work.columns and src != dst}
    return work.rename(columns=rename)


def run_target(
    objdesig: str,
    cfg: Optional[Config] = None,
    *,
    db_path: Optional[Path] = None,
    fits_roots: Optional[Sequence[Path]] = None,
    gaia: Optional[GaiaCatalog] = None,
    out_dir: Optional[Path] = None,
    status: Optional[StatusLog] = None,
    write: bool = True,
    limit: Optional[int] = None,
    progress_every: int = 200,
) -> TargetResult:
    """
    Run the full photometry pipeline for one target.

    Parameters
    ----------
    objdesig : str
        Target designation exactly as it appears in the index (e.g. ``"24P"``,
        ``"2021 G2"``).
    cfg : Config, optional
        Parameter set; the defaults are used when omitted.
    db_path : Path, optional
        Cutout index.  Defaults to ``directory.DB_PATH``.
    fits_roots : sequence of Path, optional
        Roots searched for cutout files.  Defaults to ``directory.FITS_ROOTS``.
    gaia : GaiaCatalog, optional
        Shared catalogue handle.  Passing one in is what lets a batch run open
        the 11.7 GB reference catalogue once instead of once per target.
    out_dir : Path, optional
        Where the CSV is written.  Defaults to ``directory.APPHOT_DIR``.
    status : StatusLog, optional
        Status table to update.
    write : bool
        Write the CSV and its sidecar metadata.  The notebook sets this to
        ``False`` when it only wants the DataFrame.
    limit : int, optional
        Process at most this many exposures -- for quick interactive checks.
    progress_every : int
        Log a progress line every this many exposures.

    Returns
    -------
    TargetResult

    Notes
    -----
    A missing or unreadable cutout costs one exposure, not the whole target: the
    primitive loop let a single ``FileNotFoundError`` abort a run that had
    already processed hundreds of files (review item R4).  Counts of both are
    carried into the status table.
    """
    cfg = cfg or Config()
    cfg.validate()
    slug = slugify(objdesig)
    t0 = time.time()
    result = TargetResult(objdesig=objdesig, slug=slug)

    db_path = Path(db_path) if db_path else _dir.DB_PATH
    out_dir = Path(out_dir) if out_dir else _dir.APPHOT_DIR
    resolver = FitsResolver(fits_roots if fits_roots is not None else _dir.FITS_ROOTS)

    if status is not None:
        status.update(objdesig, "running", config_hash=cfg.hash)

    try:
        # --- 1. index ---------------------------------------------------
        df = load_database(db_path, objdesig=objdesig)
        if df.empty:
            result.status, result.message = "empty", f"no rows for objdesig={objdesig!r}"
            log.warning(result.message)
            return _finish(result, t0, status, cfg)

        work = _prepare_frame(df, cfg)
        if limit:
            work = work.head(int(limit))
        n_exp = len(work)
        result.epochs = epoch_table(work, key_dateobs="date_obs")
        result.counters["n_exposures"] = n_exp
        result.counters["n_epochs"] = int(work["epoch"].nunique())

        # --- 2. pre-flight: are the cutouts actually here? ---------------
        # Resolving every path first costs one stat per file and populates the
        # resolver's cache for the loop below, but it lets a target with no data
        # exit before paying for a Gaia cone search -- which matters when
        # running against a partial dataset.
        present = sum(1 for fn in work["filename"]
                      if resolver.find(fn, objdesig) is not None)
        if present == 0:
            result.status = "empty"
            result.message = f"none of the {n_exp} cutout file(s) could be found"
            result.counters["n_missing_fits"] = n_exp
            log.warning("%s: %s", objdesig, result.message)
            return _finish(result, t0, status, cfg)
        if present < n_exp:
            log.warning("%s: %d of %d cutout file(s) missing", objdesig, n_exp - present, n_exp)

        # --- 3. Gaia ----------------------------------------------------
        gaia = gaia or GaiaCatalog(cache_dir=_dir.GAIA_CACHE_DIR, src_npy=_dir.GAIA_ALL_NPY)
        # One cone search covers every pointing of every epoch; the KD-tree
        # inside handles the union, so there is no per-epoch catalogue pass.
        subset = gaia.cone_search(work["ra"].to_numpy(float), work["dec"].to_numpy(float),
                                  radius_deg=cfg.gaia_search_radius_deg,
                                  gmag_limit=cfg.gaia_gmag_limit)
        result.gaia = subset
        result.counters["n_gaia"] = len(subset)
        # One KD-tree for the target; each exposure then pulls only its own
        # neighbours, which is both far cheaper and the only way to keep the
        # iterative SIP inverse inside its convergent regime.
        neighbors = NeighborIndex(subset)
        log.info("%s: %d exposures, %d epoch(s), %d Gaia source(s) G<%.1f",
                 objdesig, n_exp, result.counters["n_epochs"], len(subset), cfg.gaia_gmag_limit)

        # --- 4. per-exposure photometry ---------------------------------
        records: List[dict] = []
        n_ok = n_missing = n_read_err = n_ap_skipped = 0

        for i, row in enumerate(work.itertuples(index=False)):
            path = resolver.find(getattr(row, "filename"), objdesig)
            if path is None:
                n_missing += 1
                if n_missing <= 10:
                    log.warning("missing cutout: %s", getattr(row, "filename"))
                continue
            try:
                cut = read_cutout(path)
            except Exception as exc:            # noqa: BLE001 - logged, then skipped
                n_read_err += 1
                log.warning("unreadable cutout %s: %s", path.name, exc)
                continue

            badpix, _ = build_badpix_mask(
                cut.sci, cut.var, cut.flag, cfg.bad_flag_bits,
                mask_nonfinite_sci=cfg.mask_nonfinite_sci,
                mask_bad_variance=cfg.mask_bad_variance)

            pix_scale = float(getattr(row, "pix_scale"))
            fwhm_pix = _psf_fwhm_pix(float(getattr(row, "psf_fwhm_arcsec", np.nan)), pix_scale)
            kmpp = pixel_scale_km(pix_scale, float(getattr(row, "r_obs")))
            specs = build_apertures(cfg, pix_scale_arcsec=pix_scale,
                                    r_obs_au=float(getattr(row, "r_obs")),
                                    psf_fwhm_pix=fwhm_pix)
            n_ap_skipped += sum(1 for s in specs if not s.valid)

            xcen, ycen = float(getattr(row, "xcen")), float(getattr(row, "ycen"))
            r_in, r_out = annulus_radii_pix(cfg, kmpp)
            sky = estimate_sky(cut.sci, cut.var, badpix, xcen, ycen, cfg, r_in, r_out)

            # Cone radius covering the cutout's half-diagonal plus the pad that
            # `project_gaia` keeps outside the array.
            ny, nx = cut.sci.shape
            search_deg = (0.5 * np.hypot(nx, ny) + _GAIA_PAD_PIX) * pix_scale / 3600.0
            near = neighbors.query(float(getattr(row, "ra")), float(getattr(row, "dec")),
                                   search_deg)
            field_px = project_gaia(near, wcs_from_header(cut.header),
                                    xcen, ycen, cut.sci.shape, pad_pix=_GAIA_PAD_PIX)

            recs = measure_exposure_records(
                cut.sci, cut.var, badpix, xcen, ycen, specs, cfg,
                pixel_scale_km=kmpp, psf_fwhm_pix=fwhm_pix,
                vmag=float(getattr(row, "vmag", np.nan)),
                gaia_field=field_px, sky=sky)
            for r in recs:
                r["_row"] = i
                r["psf_fwhm_pix"] = fwhm_pix
            records.extend(recs)
            n_ok += 1

            if progress_every and (i + 1) % progress_every == 0:
                log.info("  %s: %d/%d exposures", objdesig, i + 1, n_exp)

        result.counters.update(n_exposures_ok=n_ok, n_missing_fits=n_missing,
                               n_read_errors=n_read_err, n_ap_skipped=n_ap_skipped)

        if not records:
            result.status = "empty"
            result.message = (f"no usable exposures "
                              f"({n_missing} missing, {n_read_err} unreadable)")
            log.warning("%s: %s", objdesig, result.message)
            return _finish(result, t0, status, cfg)

        # --- 5. assemble ------------------------------------------------
        phot = pd.DataFrame.from_records(records)
        meta_cols = [c for c in META_COLUMNS.values() if c in work.columns]
        phot = phot.merge(work[meta_cols].reset_index(drop=True),
                          left_on="_row", right_index=True, how="left").drop(columns="_row")

        # Needs r_hel / r_obs, so it can only run once the index metadata is on.
        phot = add_distance_corrected_flux(phot, cfg)

        lead = ["filename", "objdesig", "obsid", "detector", "epoch", "wl"]
        keep = PHOT_COLUMNS + DISTCORR_COLUMNS
        ordered = ([c for c in lead if c in phot.columns]
                   + [c for c in keep if c in phot.columns]
                   + [c for c in phot.columns
                      if c not in lead and c not in keep])
        phot = phot[ordered]

        result.phot = phot
        result.counters["n_rows_out"] = len(phot)
        result.counters["n_badphot"] = int(phot["badphot"].sum())
        log.info("%s: %d rows (%d exposures x apertures), %d flagged badphot",
                 objdesig, len(phot), n_ok, result.counters["n_badphot"])

        # --- 6. write ---------------------------------------------------
        if write:
            result.output_csv = _write_outputs(phot, out_dir, slug, objdesig, cfg, result)

    except Exception as exc:                     # noqa: BLE001 - recorded, not hidden
        result.status = "failed"
        result.message = f"{type(exc).__name__}: {exc}"
        log.error("%s FAILED: %s", objdesig, result.message)
        log.debug("traceback:\n%s", traceback.format_exc())

    return _finish(result, t0, status, cfg)


def format_jd_columns(phot: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Copy of ``phot`` with the ``jd_*`` columns rendered at full precision.

    ``to_csv(float_format="%.8g")`` would round a Julian date of 2.46e6 to 0.1 day;
    the catalog pipeline groups exposures by time and needs the second.
    """
    out = phot.copy()
    for c in out.columns:
        if c.startswith("jd_") and pd.api.types.is_float_dtype(out[c]):
            out[c] = out[c].map(lambda v: (cfg.jd_float_format % v) if np.isfinite(v) else "")
    return out


def _write_outputs(phot: pd.DataFrame, out_dir: Path, slug: str,
                   objdesig: str, cfg: Config, result: TargetResult) -> Path:
    """Write the photometry CSV plus a JSON sidecar recording how it was made."""
    import json

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{slug}.csv"
    format_jd_columns(phot, cfg).to_csv(csv_path, index=False, float_format=cfg.float_format)

    meta = {
        "objdesig": objdesig, "slug": slug, "written_utc": utcnow_iso(),
        "n_rows": int(len(phot)), "flux_unit": cfg.flux_unit,
        "config_hash": cfg.hash, "config": cfg.to_dict(),
        "counters": result.counters,
        "epochs": (result.epochs.to_dict("records")
                   if result.epochs is not None and not result.epochs.empty else []),
    }
    (out_dir / f"{slug}.meta.json").write_text(
        json.dumps(meta, indent=2, default=str))
    log.info("wrote %s (%d rows) and %s.meta.json", csv_path, len(phot), slug)
    return csv_path


def _finish(result: TargetResult, t0: float,
            status: Optional[StatusLog], cfg: Config) -> TargetResult:
    result.elapsed_s = time.time() - t0
    if status is not None:
        status.update(result.objdesig, result.status,
                      config_hash=cfg.hash, **result.status_fields())
    return result


# ---------------------------------------------------------------------------
# Band stacking
# ---------------------------------------------------------------------------
def stack_target(
    phot: pd.DataFrame,
    cfg: Optional[Config] = None,
    *,
    objdesig: Optional[str] = None,
    ap_label: Optional[str] = None,
    fits_roots: Optional[Sequence[Path]] = None,
    out_dir: Optional[Path] = None,
    epochs: Optional[Sequence[int]] = None,
    write: bool = True,
) -> Dict[tuple, object]:
    """
    Build robust band stacks for every epoch of one target.

    Parameters
    ----------
    phot : pandas.DataFrame
        Photometry table from :func:`run_target`.
    cfg : Config, optional
    objdesig : str, optional
        Overrides the designation taken from ``phot``.
    ap_label : str, optional
        Aperture whose measurements supply the sky level and the label.
        Defaults to the smallest valid fixed-pixel aperture.
    fits_roots, out_dir : optional
        As in :func:`run_target`.
    epochs : sequence of int, optional
        Restrict to these epochs.
    write : bool
        Write each stack to ``COMBFITS_DIR``.

    Returns
    -------
    dict
        ``{(epoch, band): StackResult}``.  Bands with fewer than
        ``cfg.stack_min_frames`` usable frames are omitted and logged.

    Notes
    -----
    Frames are background-subtracted, registered to a common sub-pixel centre,
    optionally rotated to north-up and combined with a sigma-clipped mean --
    see :mod:`spherex_apphot.stacking` for why each of those matters.
    """
    from .stacking import StackInput, stack_frames, write_stack_fits

    cfg = cfg or Config()
    out_dir = Path(out_dir) if out_dir else _dir.COMBFITS_DIR
    resolver = FitsResolver(fits_roots if fits_roots is not None else _dir.FITS_ROOTS)
    objdesig = objdesig or str(phot["objdesig"].iloc[0])
    slug = slugify(objdesig)

    if ap_label is None:
        pix_aps = sorted(phot.loc[phot["ap_kind"] == "pix", "ap_label"].unique())
        ap_label = pix_aps[0] if pix_aps else str(phot["ap_label"].iloc[0])
    sub_all = phot[phot["ap_label"] == ap_label]
    if sub_all.empty:
        log.warning("no rows for aperture %r; nothing to stack", ap_label)
        return {}
    log.info("stacking %s using the %s aperture for sky levels", objdesig, ap_label)

    out: Dict[tuple, object] = {}
    epoch_list = sorted(sub_all["epoch"].unique()) if epochs is None else list(epochs)

    for epoch in epoch_list:
        ep = sub_all[sub_all["epoch"] == epoch]
        for band, (lo, hi) in cfg.stack_bands.items():
            rows = ep[(ep["wl"] >= lo) & (ep["wl"] <= hi)]
            if rows.empty:
                continue
            items: List[StackInput] = []
            for row in rows.itertuples(index=False):
                path = resolver.find(getattr(row, "filename"), objdesig)
                if path is None:
                    continue
                try:
                    cut = read_cutout(path)
                except Exception as exc:          # noqa: BLE001
                    log.warning("stack: skipping unreadable %s (%s)", path.name, exc)
                    continue
                badpix, _ = build_badpix_mask(
                    cut.sci, cut.var, cut.flag, cfg.bad_flag_bits,
                    mask_nonfinite_sci=cfg.mask_nonfinite_sci,
                    mask_bad_variance=cfg.mask_bad_variance)
                xcen, ycen = float(getattr(row, "xcen")), float(getattr(row, "ycen"))
                north = (north_angle_deg(wcs_from_header(cut.header), xcen, ycen)
                         if cfg.stack_align_north else 0.0)
                items.append(StackInput(
                    sci=cut.sci, badpix=badpix, xcen=xcen, ycen=ycen,
                    sky_median=float(getattr(row, "sky_median_mjy_per_pix", 0.0) or 0.0),
                    north_deg=north, label=getattr(row, "filename")))

            if len(items) < cfg.stack_min_frames:
                log.info("epoch %s band %s: only %d frame(s) (< %d); skipped",
                         epoch, band, len(items), cfg.stack_min_frames)
                continue

            res = stack_frames(
                items, radius=cfg.stack_radius_pix,
                subtract_sky=cfg.stack_subtract_sky,
                subpixel=cfg.stack_subpixel_shift,
                align_north=cfg.stack_align_north,
                sigma=cfg.stack_sigma, maxiters=cfg.stack_maxiters,
                combine=cfg.stack_combine,
                pix_scale_arcsec=float(rows["pix_scale"].median()))
            if res is None:
                log.warning("epoch %s band %s: no usable pixels", epoch, band)
                continue

            out[(int(epoch), band)] = res
            if write:
                write_stack_fits(
                    res, out_dir / f"{slug}_epoch{int(epoch)}_{band}.fits",
                    objdesig=objdesig, epoch=int(epoch), band=band,
                    wl_range=(lo, hi), flux_unit=cfg.flux_unit,
                    extra={"AP_LABEL": (ap_label, "Aperture supplying sky levels"),
                           "CFGHASH": (cfg.hash, "Config hash")})
    return out
