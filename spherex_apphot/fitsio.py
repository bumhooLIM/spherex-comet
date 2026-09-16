"""
Reading cutout FITS files and the Parquet cutout index.

Two problems from the review are fixed here.

R2 -- the primitive code addressed extensions positionally (``hdul[1]``,
``hdul[2]``, ``hdul[3]``).  That happens to be correct for the current export,
but the writer emits whichever of ``IMAGE / MASK / FLAGS / VARIANCE / FLAG``
exist, in that order, so a product that also carried ``MASK`` would silently
hand ``sqrt(FLAG)`` to the code as its error array.  Extensions are now looked
up by ``EXTNAME``, with the positional order used only as an explicitly logged
fallback.

R4 -- the photometry loop opened files without guarding, so one missing cutout
aborted the whole target after all its work.  :meth:`FitsResolver.find` returns
``None`` for a missing file and the caller skips that one row.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional, Sequence

import numpy as np
import pandas as pd
from astropy.io import fits

from .logging_utils import get_logger
from .status import slugify

__all__ = ["Cutout", "FitsResolver", "read_cutout", "load_database", "list_targets"]

log = get_logger("fitsio")

_SCI_NAMES = ("IMAGE", "SCI", "SCIENCE")
_VAR_NAMES = ("VARIANCE", "VAR")
_FLAG_NAMES = ("FLAG", "FLAGS", "DQ")


class Cutout(NamedTuple):
    """One cutout: science, variance and flag planes plus the science header."""

    sci: np.ndarray          #: float64, mJy/pixel
    var: np.ndarray          #: float64, (mJy/pixel)^2
    flag: np.ndarray         #: uint32 bitmask
    header: fits.Header      #: header of the science extension (carries the WCS)
    path: Path


class FitsResolver:
    """
    Locate a cutout FITS file by name across several possible roots.

    Two layouts are supported and probed in order for each root:
    ``<root>/<filename>`` (the flat production export) and
    ``<root>/<objdesig>/<filename>`` (the per-target ``data/spherex_sample/`` tree).

    Successful lookups are memoised, so a second pass over the same target --
    the stacking step, say -- costs no further filesystem calls.
    """

    def __init__(self, roots: Sequence[Path]):
        self.roots: List[Path] = [Path(r) for r in roots]
        self._cache: Dict[str, Optional[Path]] = {}
        self.n_missing = 0

    def find(self, filename: str, objdesig: Optional[str] = None) -> Optional[Path]:
        """Return the path to ``filename``, or ``None`` if no root holds it."""
        key = f"{objdesig}/{filename}"
        if key in self._cache:
            return self._cache[key]

        candidates: List[Path] = []
        for root in self.roots:
            candidates.append(root / filename)
            if objdesig:
                candidates.append(root / str(objdesig) / filename)
                candidates.append(root / slugify(objdesig) / filename)

        found = next((c for c in candidates if c.is_file()), None)
        if found is None:
            self.n_missing += 1
        self._cache[key] = found
        return found

    def describe(self) -> str:
        return "FITS roots: " + ", ".join(str(r) for r in self.roots)


def _get_hdu(hdul: fits.HDUList, names: Iterable[str], fallback_index: int) -> fits.ImageHDU:
    """Fetch an extension by EXTNAME, falling back to position with a warning."""
    for name in names:
        try:
            return hdul[name]
        except KeyError:
            continue
    log.warning(
        "no extension named %s in %s; falling back to index %d (%r)",
        "/".join(names), hdul.filename(), fallback_index,
        hdul[fallback_index].name,
    )
    return hdul[fallback_index]


def read_cutout(path: Path) -> Cutout:
    """
    Read one cutout FITS file.

    Returns
    -------
    Cutout
        Planes are returned as float64 (science, variance) and uint32 (flag).
        Float64 costs nothing at 91x91 and avoids float32 round-off
        accumulating through the aperture sums.

    Raises
    ------
    OSError, ValueError
        Propagated from astropy for a truncated or malformed file; the caller
        records the failure against the individual exposure and continues.
    """
    with fits.open(path, memmap=False) as hdul:
        sci_hdu = _get_hdu(hdul, _SCI_NAMES, 1)
        var_hdu = _get_hdu(hdul, _VAR_NAMES, 2)
        flag_hdu = _get_hdu(hdul, _FLAG_NAMES, 3)

        sci = np.asarray(sci_hdu.data, dtype=np.float64)
        var = np.asarray(var_hdu.data, dtype=np.float64)
        # int32 -> uint32 keeps bit 31 addressable and matches the `1 << bit`
        # tests in masking.py.
        flag = np.asarray(flag_hdu.data).astype(np.uint32, copy=False)
        header = sci_hdu.header.copy()

    if not (sci.shape == var.shape == flag.shape):
        raise ValueError(
            f"{path.name}: plane shapes disagree "
            f"(sci {sci.shape}, var {var.shape}, flag {flag.shape})"
        )
    return Cutout(sci=sci, var=var, flag=flag, header=header, path=Path(path))


def load_database(
    db_path: Path,
    objdesig: Optional[str] = None,
    columns: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    Load the cutout index, optionally for a single target.

    Uses Parquet predicate pushdown so that only the requested target's row
    group data is materialised -- the full index is ~104k rows x 240 columns.

    Parameters
    ----------
    db_path : Path
        Path to ``db_filtered.parq``.
    objdesig : str, optional
        Restrict to one target designation.
    columns : sequence of str, optional
        Restrict to these columns.

    Returns
    -------
    pandas.DataFrame
        Row order is not meaningful; :func:`epochs.group_epochs` sorts by time.
    """
    db_path = Path(db_path)
    if not db_path.is_file():
        raise FileNotFoundError(f"cutout index not found: {db_path}")
    filters = [("objdesig", "==", objdesig)] if objdesig is not None else None
    df = pd.read_parquet(db_path, filters=filters, columns=list(columns) if columns else None)
    log.info("loaded %d rows from %s%s", len(df), db_path.name,
             f" for objdesig={objdesig!r}" if objdesig else "")
    return df


def list_targets(db_path: Path) -> List[str]:
    """Return every unique ``objdesig`` in the index, sorted."""
    df = pd.read_parquet(Path(db_path), columns=["objdesig"])
    return sorted(df["objdesig"].dropna().unique().tolist())
