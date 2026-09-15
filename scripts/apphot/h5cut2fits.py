"""
Expand HDF5 cutout archive into individual FITS files.

This script reads the Parquet index and HDF5 data created by extract_ssos.py
and exports individual FITS files with full headers reconstructed from
the Parquet database.

"""

import os
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import click
import h5py
import numpy as np
import pqfilt
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning


# Dict mapping Parquet column -> (FITS keyword, comment)
# Single source of truth for ephemeris keyword translation
EPH_KEYWORDS = {
    "wl": ("WLEN-CEN", "[um] Central wavelength"),
    "objdesig": ("OBJDESIG", "Object designation"),
    "ra": ("RA-OBJ", "[deg] RA of object"),
    "dec": ("DEC-OBJ", "[deg] Dec of object"),
    "r_hel": ("R_HEL", "[au] Heliocentric distance"),
    "r_obs": ("R_OBS", "[au] Observer distance"),
    "alpha": ("ALPHA", "[deg] Phase angle"),
    "hel_ecl_lon": ("HECL-LON", "[deg] Heliocentric ecl. lon"),
    "hel_ecl_lat": ("HECL-LAT", "[deg] Heliocentric ecl. lat"),
    "obs_ecl_lon": ("OECL-LON", "[deg] Observer ecl. lon"),
    "obs_ecl_lat": ("OECL-LAT", "[deg] Observer ecl. lat"),
    "racosdec_rate": ("RA_RATE", "[arcsec/min] RA*cos(DEC) motion"),
    "dec_rate": ("DEC_RATE", "[arcsec/min] Dec motion"),
    "sky_motion": ("SKY-MOT", "[arcsec/min] Sky motion"),
    "sky_motion_pa": ("SKY-PA", "[deg] Motion PA"),
    "vmag": ("VMAG", "[mag] Visual magnitude"),
    "jd_utc": ("MIDJDUTC", "[day] JD midpoint UTC"),
    "jd_tdb": ("MIDJDTDB", "[day] JD midpoint TDB"),
    "pix_scale": ("PIX-SCL", "[arcsec/pix] Pixel scale"),
    "sun_jy": ("SUN-JY", "[Jy] Sun, convolved to the pixel bandpass"),
    # Photometry (aperture 1: 2 pix)
    "ap1_rad": ("AP-RAD1", "[pix] AP1 radius"),
    "ap1_sum": ("AP-SUM1", "[MJy/sr*pixscale^2] AP1 bkg-subtracted sum"),
    "ap1_err": ("AP-ERR1", "[MJy/sr*pixscale^2] AP1 sum error"),
    # Photometry (aperture 2: 12.3 arcsec)
    "ap2_rad": ("AP-RAD2", "[pix] AP2 radius"),
    "ap2_sum": ("AP-SUM2", "[MJy/sr*pixscale^2] AP2 bkg-subtracted sum"),
    "ap2_err": ("AP-ERR2", "[MJy/sr*pixscale^2] AP2 sum error"),
    # Sky (shared annulus)
    "ap_msky": ("AP-MSKY", "[MJy/sr*pixscale^2] Background median"),
    "ap_ssky": ("AP-SSKY", "[MJy/sr*pixscale^2] Background std"),
    "ap_nsky": ("AP-NSKY", "Num of original sky pixel count"),
    "ap_nrej": ("AP-NREJ", "Rejected sky pixels from sigclip"),
    # Annulus
    "an_rin": ("AN-RIN", "[pix] Sky annulus inner radius"),
    "an_rout": ("AN-ROUT", "[pix] Sky annulus outer radius"),
}

# Columns to skip when copying row fields to the IMAGE header
_SKIP_COLS = {
    "hdf5_file", "hdf5_path", "objdesig", "obsid", "detector",
    "xcen", "ycen", "ltv1", "ltv2", "cutout_size",
    "sky_obsids", "n_sky",
    *EPH_KEYWORDS.keys(),
}

# ---------------------------------------------------------------------------
# Per-worker globals (set by _worker_init, used by _export_one)
# ---------------------------------------------------------------------------
_worker_h5_handles: dict = {}
_worker_outdir: Path = Path(".")
_worker_keys: list = []
_worker_delimiter: str = "_"


def _worker_init(h5_paths: list[str], outdir: str, keys: list[str], delimiter: str):
    """Initialiser for each worker process: opens HDF5 files."""
    global _worker_h5_handles, _worker_outdir, _worker_keys, _worker_delimiter
    _worker_outdir = Path(outdir)
    _worker_keys = keys
    _worker_delimiter = delimiter
    _worker_h5_handles = {}
    for p in h5_paths:
        try:
            _worker_h5_handles[Path(p).name] = h5py.File(p, "r")
        except OSError:
            pass


def build_cutout_wcs(row):
    """Build WCS header keywords for a cutout.

    Parameters
    ----------
    row : dict
        Denormalized row containing both frame WCS and cutout geometry.

    Returns
    -------
    dict
        WCS keywords for the cutout
    """
    xcen = row["xcen"]
    ycen = row["ycen"]
    ltv1 = row["ltv1"]
    ltv2 = row["ltv2"]
    cutout_size = row["cutout_size"]

    wcs_keys = {
        "NAXIS1": cutout_size,
        "NAXIS2": cutout_size,
        "LTV1": ltv1,
        "LTV2": ltv2,
        "XCEN": xcen,
        "YCEN": ycen,
    }

    if row.get("CRPIX1") is not None:
        wcs_keys["CRPIX1"] = row["CRPIX1"] + ltv1
    if row.get("CRPIX2") is not None:
        wcs_keys["CRPIX2"] = row["CRPIX2"] + ltv2

    return wcs_keys


def build_eph_header(cutout_row: dict) -> dict:
    """Build ephemeris header keywords from cutout row.

    Parameters
    ----------
    cutout_row : dict
        Row from cutouts table

    Returns
    -------
    dict
        Ephemeris keywords for the header, mapping FITS keyword to (value, comment)
    """
    result = {}
    for db_col, (fits_key, comment) in EPH_KEYWORDS.items():
        val = cutout_row.get(db_col)
        if val is None or (isinstance(val, float) and np.isnan(val)):
            continue
        result[fits_key] = (val, comment)
    return result


def _export_one(row_dict):
    """Build and write a single FITS file from a cutout row.

    Uses module-level worker globals set by ``_worker_init``.

    Returns
    -------
    str or None
        Output filename on success, ``None`` on skip.
    """
    keys = _worker_keys
    delimiter = _worker_delimiter
    outdir = _worker_outdir
    h5_handles = _worker_h5_handles

    # Build filename
    parts = []
    for key in keys:
        key_lower = key.lower()
        if key_lower in row_dict:
            val = row_dict[key_lower]
            if key_lower == "wl" and val is not None:
                parts.append(f"{val:.4f}um")
            elif key_lower == "vmag" and val is not None:
                parts.append(f"{val:.2f}Vmag")
            elif val is not None:
                parts.append(str(val))

    if not parts:
        return None

    fname = delimiter.join(parts) + ".fits"

    # Find the HDF5 group
    hdf5_grp_path = row_dict["hdf5_path"]
    hdf5_file = row_dict.get("hdf5_file")

    h5f = h5_handles.get(hdf5_file)
    if h5f is None:
        for name, handle in h5_handles.items():
            if name == hdf5_file or name in hdf5_file:
                h5f = handle
                break
    if h5f is None:
        return None

    if hdf5_grp_path not in h5f:
        return None

    grp = h5f[hdf5_grp_path]

    # Read arrays
    ext_names = ["IMAGE", "MASK", "FLAGS", "VARIANCE", "FLAG"]
    arrays = {}
    for ext_name in ext_names:
        if ext_name in grp:
            arrays[ext_name] = grp[ext_name][:]
    if "SIM" in grp:
        arrays["SIM"] = grp["SIM"][:]

    # Build WCS and ephemeris headers
    cutout_wcs = build_cutout_wcs(row_dict)
    eph_header = build_eph_header(row_dict)

    # Create FITS HDU list
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header["OBJDESIG"] = row_dict.get("objdesig", "")
    primary_hdu.header["OBSID"] = row_dict.get("obsid", "")
    primary_hdu.header["DETECTOR"] = row_dict.get("detector", 0)

    hdu_list = [primary_hdu]

    for j, ext_name in enumerate(ext_names):
        if ext_name not in arrays:
            continue
        img_hdu = fits.ImageHDU(data=arrays[ext_name], name=ext_name)

        # Apply full header to IMAGE extension only
        if j == 0:
            for key, val in row_dict.items():
                if key.lower() in _SKIP_COLS or key in _SKIP_COLS:
                    continue
                if val is None:
                    continue
                if isinstance(val, float) and np.isnan(val):
                    continue
                fits_key = key.upper()
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", category=VerifyWarning)
                        img_hdu.header[fits_key] = val
                except (ValueError, TypeError):
                    pass

            for key, val in cutout_wcs.items():
                img_hdu.header[key] = val

            for key, (val, comment) in eph_header.items():
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=VerifyWarning)
                    img_hdu.header[key] = (val, comment)

        hdu_list.append(img_hdu)

    if "SIM" in arrays:
        hdu_list.append(fits.ImageHDU(data=arrays["SIM"], name="SIM"))

    hdul = fits.HDUList(hdu_list)
    hdul.writeto(outdir / fname, overwrite=True)
    return fname


def h5cut2fits(
    hdf5_paths,
    db_path=None,
    outdir="./FITS_OUT",
    *,
    filters=None,
    keys=("objdesig", "wl", "vmag", "obsid", "detector"),
    delimiter="_",
    limit=None,
    nworkers=4,
):
    """Expand HDF5 cutout archive into individual FITS files.

    Parameters
    ----------
    hdf5_paths : str, Path, or list thereof
        Path(s) to HDF5 cutout file(s). Shell globs are expanded.
    db_path : str, Path, or None
        Path to the cutout index Parquet file (denormalized).
        If ``None``, defaults to ``db.parq`` in the parent directory
        of the first HDF5 file.
    outdir : str or Path, optional
        Output directory for FITS files (default: ``"./FITS_OUT"``).
    filters : str, list, or pqfilt.ExprNode, optional
        Filter expression passed to ``pqfilt.read()``.  Accepts:

        - Expression string: ``"vmag < 20"``,
          ``"(wl > 1.0 & wl < 2.0) | detector == 3"``
        - List of 3-tuples (AND): ``[("wl", ">", 1.0), ("vmag", "<", 20)]``
        - Pre-parsed AST node (``pqfilt.FilterExpr`` etc.)
    keys : sequence of str, optional
        Column names used to build output filenames (default:
        ``("objdesig", "wl", "vmag", "obsid", "detector")``).
        Better not change this to include other columns...
    delimiter : str, optional
        Delimiter between filename parts (default: ``"_"``).
    limit : int or None, optional
        Maximum number of FITS files to export.
    nworkers : int, optional
        Number of worker processes for parallel FITS export (default: 4).
        Each worker opens its own HDF5 file handles for true parallelism.
        Set to 1 for sequential execution.

    Returns
    -------
    int
        Number of FITS files exported.

    Examples
    --------
    >>> from spherex_l1l2tools.fileio.h5cut2fits import h5cut2fits
    >>> h5cut2fits("cutouts.h5", filters="objdesig == '3200'")
    >>> h5cut2fits(["a.h5", "b.h5"],
    ...           filters="wl > 1.0 & vmag <= 20", nworkers=8)
    """
    import glob as glob_module

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Normalise hdf5_paths to a list of expanded strings
    if isinstance(hdf5_paths, (str, Path)):
        hdf5_paths = [str(hdf5_paths)]
    else:
        hdf5_paths = [str(p) for p in hdf5_paths]

    expanded_h5_paths = []
    for p in hdf5_paths:
        if "*" in p or "?" in p:
            expanded_h5_paths.extend(glob_module.glob(p))
        else:
            expanded_h5_paths.append(p)

    if not expanded_h5_paths:
        raise FileNotFoundError("No HDF5 files found from the provided paths.")

    # Validate HDF5 files (quick open/close to weed out bad files early)
    valid_h5_paths = []
    for p in expanded_h5_paths:
        try:
            with h5py.File(p, "r"):
                valid_h5_paths.append(p)
        except OSError as e:
            print(f"Warning: cannot open {p} as HDF5 ({e}), skipping.")

    if not valid_h5_paths:
        raise FileNotFoundError("No valid HDF5 files found.")

    # Default db_path: db.parq in the same directory as the first HDF5 file
    if db_path is None:
        db_path = Path(valid_h5_paths[0]).parent / "db.parq"
        print(f"Using default Parquet DB: {db_path}")
    db_path = Path(db_path)

    # Normalise keys
    if isinstance(keys, str):
        keys = [k.strip() for k in keys.split(",") if k.strip()]
    else:
        keys = list(keys)

    # Load Parquet database with predicate-pushdown filtering via pqfilt
    df = pqfilt.read(db_path, filters=filters)

    if limit:
        df = df.head(limit)

    if df.empty:
        print("No rows found matching the filters.")
        return 0

    rows = df.to_dict(orient="records")
    total = len(rows)
    print(f"Found {total} cutouts to export (nworkers={nworkers}).")

    exported_count = 0

    if nworkers <= 1:
        # Sequential: use in-process HDF5 handles (no fork overhead)
        _worker_init(valid_h5_paths, str(outdir), keys, delimiter)
        try:
            for i, row_dict in enumerate(rows):
                result = _export_one(row_dict)
                if result is not None:
                    exported_count += 1
                if exported_count % 500 == 0 and exported_count > 0:
                    print(f"Exported {exported_count}/{total} files...")
        finally:
            for h5f in _worker_h5_handles.values():
                h5f.close()
    else:
        # Parallel: each worker process opens its own HDF5 handles
        with ProcessPoolExecutor(
            max_workers=nworkers,
            initializer=_worker_init,
            initargs=(valid_h5_paths, str(outdir), keys, delimiter),
        ) as pool:
            for result in pool.map(_export_one, rows, chunksize=64):
                if result is not None:
                    exported_count += 1
                if exported_count % 500 == 0 and exported_count > 0:
                    print(f"Exported {exported_count}/{total} files...")

    print(f"Done. Exported {exported_count} FITS files to {outdir}")
    return exported_count


# ---------------------------------------------------------------------------
# CLI wrapper
# ---------------------------------------------------------------------------


@click.command()
@click.option(
    "--db",
    "db_path",
    default=None,
    help="Path to cutout index Parquet file (default: db.parq next to HDF5)",
)
@click.option(
    "--h5",
    "hdf5_paths",
    multiple=True,
    help="Path to HDF5 cutout file(s) (can repeat; use quotes for globs)",
)
@click.argument("hdf5_files", nargs=-1, type=click.Path(exists=True))
@click.option(
    "--outdir",
    "outdir",
    default="./FITS_OUT",
    show_default=True,
    help="Output directory for FITS files",
)
@click.option(
    "-f",
    "--filter",
    "filter_exprs",
    multiple=True,
    help=(
        'Filter expression (pqfilt syntax), e.g. "vmag < 20", '
        '"wl > 1.0 & wl < 2.0". Multiple -f flags are AND-ed together.'
    ),
)
@click.option(
    "--keys",
    "keys_str",
    default="objdesig,wl,vmag,obsid,detector",
    show_default=True,
    help="Column names for filename (comma separated)",
)
@click.option(
    "--delimiter",
    "delimiter",
    default="_",
    show_default=True,
    help="Delimiter between filename parts",
)
@click.option(
    "--limit", "limit", default=None, type=int, help="Limit number of files to export"
)
@click.option(
    "-j",
    "--nworkers",
    "nworkers",
    default=4,
    show_default=True,
    type=int,
    help="Number of worker processes for parallel FITS export",
)
def main(
    db_path,
    hdf5_paths,
    hdf5_files,
    outdir,
    filter_exprs,
    keys_str,
    delimiter,
    limit,
    nworkers,
):
    """Expand HDF5 cutout archive into individual FITS files with full headers.

    Examples::

        spherex-h5cut2fits cutouts.h5 -f "objdesig == '3200'"
        spherex-h5cut2fits cutouts.h5 -f "vmag < 20"
        spherex-h5cut2fits cutouts.h5 -f "wl > 1.0" -f "wl < 2.0" -j 8

    After generating files, maybe you can zip them up with something like::

        tar cf - folder/ | pv | pigz -p8 -v > folder.tar.gz

    later unzip with::

        pigz -p8 -d -c folder.tar.gz | tar xf -
    """
    all_h5_paths = list(hdf5_paths) + list(hdf5_files)

    # Combine multiple -f expressions with AND (same pattern as pqfilt CLI)
    combined_filter = None
    if filter_exprs:
        if len(filter_exprs) == 1:
            combined_filter = filter_exprs[0]
        else:
            combined_filter = " & ".join(f"({expr})" for expr in filter_exprs)

    h5cut2fits(
        hdf5_paths=all_h5_paths,
        db_path=db_path,
        outdir=outdir,
        filters=combined_filter,
        keys=keys_str,
        delimiter=delimiter,
        limit=limit,
        nworkers=nworkers,
    )


if __name__ == "__main__":
    main()
