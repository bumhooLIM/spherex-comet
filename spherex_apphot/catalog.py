"""
Gaia DR3 access.

The reference catalogue is a 419-million-row, 11.7 GB ``.npy``.  The primitive
pipeline boolean-indexed the whole memory map on magnitude once per target --
and, because the batch driver spawned one subprocess per target, paid that
11.7 GB read 445 times over (review item P1).  It then intersected the result
with each pointing in a Python loop, one full pass per pointing (review item P2).

Two changes remove both costs:

* **Declination cache.**  ``build_gaia_cache.py`` writes ``dec``, ``ra`` and
  ``gmag`` as three ``.npy`` files sorted by declination.  A query then
  binary-searches the declination array on the memory map and reads only the
  matching band -- a few MB instead of 11.7 GB.  A comet footprint spans
  ~0.4 deg, so this is a reduction of roughly three orders of magnitude.  A
  magnitude cut alone would not have helped: a G < 18 cut keeps 77 % of this
  already magnitude-limited catalogue.
* **KD-tree intersection.**  Pointings are matched against the declination band
  with a single ``cKDTree`` on 3-vectors, so cost scales as
  ``O((N_band + N_point) log N)`` instead of ``O(N_point x N_gaia)``.

Without a cache the class falls back to one chunked pass over the source file,
which is still a single pass per *run* rather than per target.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np

from .logging_utils import get_logger

__all__ = ["GaiaCatalog", "GaiaSubset", "NeighborIndex",
           "build_dec_cache", "radec_to_unit"]

log = get_logger("catalog")

_CACHE_FILES = ("dec.npy", "ra.npy", "gmag.npy")
_META = "meta.json"


def radec_to_unit(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    """
    Convert spherical coordinates to unit 3-vectors.

    Chord distance in this space is monotonic in angular separation, which is
    what lets a Euclidean KD-tree answer cone searches exactly and without the
    ``1/cos(dec)`` special-casing (and pole breakdown) of a box search.
    """
    ra = np.radians(np.asarray(ra_deg, dtype=np.float64))
    dec = np.radians(np.asarray(dec_deg, dtype=np.float64))
    cd = np.cos(dec)
    return np.column_stack((cd * np.cos(ra), cd * np.sin(ra), np.sin(dec)))


def _chord(radius_deg: float) -> float:
    """Chord length on the unit sphere subtended by ``radius_deg``."""
    return 2.0 * np.sin(np.radians(float(radius_deg)) / 2.0)


class GaiaSubset:
    """A small in-memory set of Gaia sources: ``ra``, ``dec``, ``gmag``."""

    __slots__ = ("ra", "dec", "gmag")

    def __init__(self, ra: np.ndarray, dec: np.ndarray, gmag: np.ndarray):
        self.ra = np.asarray(ra, dtype=np.float64)
        self.dec = np.asarray(dec, dtype=np.float64)
        self.gmag = np.asarray(gmag, dtype=np.float64)

    def __len__(self) -> int:
        return int(self.ra.size)

    def __repr__(self) -> str:
        if not len(self):
            return "GaiaSubset(empty)"
        return (f"GaiaSubset(n={len(self)}, "
                f"G={np.nanmin(self.gmag):.2f}-{np.nanmax(self.gmag):.2f})")

    def to_dataframe(self):
        """Return the subset as a DataFrame (for plotting in the notebook)."""
        import pandas as pd
        return pd.DataFrame({"ra": self.ra, "dec": self.dec, "phot_g_mean_mag": self.gmag})


def build_dec_cache(
    src_npy: Path,
    out_dir: Path,
    gmag_limit: Optional[float] = None,
    chunk: int = 20_000_000,
    bucket_deg: float = 0.5,
) -> Path:
    """
    Build the declination-sorted cache used by :class:`GaiaCatalog`.

    Parameters
    ----------
    src_npy : Path
        The full Gaia ``.npy`` with fields ``ra``, ``dec``, ``phot_g_mean_mag``.
    out_dir : Path
        Destination directory; ``dec.npy``, ``ra.npy``, ``gmag.npy`` and
        ``meta.json`` are written into it.
    gmag_limit : float, optional
        Drop sources fainter than this.  ``None`` keeps everything, which is
        usually right -- the source catalogue is already magnitude-limited and a
        G < 18 cut removes only ~23 % of it.
    chunk : int
        Rows read per pass over the source file.
    bucket_deg : float
        Width of the declination buckets used by the external sort.

    Returns
    -------
    Path
        ``out_dir``.

    Notes
    -----
    A single in-memory ``argsort`` of 419 million rows would need roughly 15 GB
    of working set, which is not safe on a laptop.  This is therefore an
    **external bucket sort**: pass one streams the catalogue and appends each
    row to a per-declination-bucket scratch file; pass two sorts one bucket at a
    time and writes it straight into an output memory map.  Peak memory is set
    by the largest single bucket (tens of MB), not by the catalogue size.

    Everything is written under a ``.tmp`` sibling and renamed on success, so an
    interrupted build cannot leave a partial cache that later looks valid.
    """
    src_npy, out_dir = Path(src_npy), Path(out_dir)
    if not src_npy.is_file():
        raise FileNotFoundError(f"Gaia catalogue not found: {src_npy}")

    src = np.load(src_npy, mmap_mode="r")
    n_total = len(src)
    tmp = out_dir.with_name(out_dir.name + ".tmp")
    if tmp.exists():
        shutil.rmtree(tmp)
    scratch = tmp / "buckets"
    scratch.mkdir(parents=True)

    n_buckets = int(np.ceil(180.0 / float(bucket_deg)))
    log.info("building declination cache from %s (%d rows, %d buckets)",
             src_npy.name, n_total, n_buckets)

    # -- pass 1: scatter into declination buckets --------------------------
    counts = np.zeros(n_buckets, dtype=np.int64)
    handles: dict = {}
    try:
        for i0 in range(0, n_total, chunk):
            i1 = min(i0 + chunk, n_total)
            dec = np.asarray(src["dec"][i0:i1], dtype=np.float64)
            g = np.asarray(src["phot_g_mean_mag"][i0:i1], dtype=np.float32)
            keep = np.isfinite(dec) & np.isfinite(g)
            if gmag_limit is not None:
                keep &= g < float(gmag_limit)
            if not keep.any():
                continue
            dec = dec[keep]
            g = g[keep]
            ra = np.asarray(src["ra"][i0:i1], dtype=np.float64)[keep]

            b = np.clip(((dec + 90.0) / bucket_deg).astype(np.int64), 0, n_buckets - 1)
            order = np.argsort(b, kind="stable")
            b, dec, ra, g = b[order], dec[order], ra[order], g[order]
            edges = np.searchsorted(b, np.arange(n_buckets + 1))
            for k in range(n_buckets):
                lo, hi = edges[k], edges[k + 1]
                if lo == hi:
                    continue
                fh = handles.get(k)
                if fh is None:
                    fh = handles[k] = open(scratch / f"b{k:04d}.bin", "wb")
                rec = np.empty(hi - lo, dtype=[("dec", "f8"), ("ra", "f8"), ("g", "f4")])
                rec["dec"], rec["ra"], rec["g"] = dec[lo:hi], ra[lo:hi], g[lo:hi]
                fh.write(rec.tobytes())
                counts[k] += hi - lo
            log.info("  pass 1: %d/%d rows scanned, %d kept", i1, n_total, int(counts.sum()))
    finally:
        for fh in handles.values():
            fh.close()

    n_kept = int(counts.sum())
    if n_kept == 0:
        shutil.rmtree(tmp)
        raise ValueError("no Gaia rows survived the magnitude cut")

    # -- pass 2: sort each bucket into the output memory maps --------------
    dec_out = np.lib.format.open_memmap(tmp / "dec.npy", mode="w+", dtype="f8", shape=(n_kept,))
    ra_out = np.lib.format.open_memmap(tmp / "ra.npy", mode="w+", dtype="f8", shape=(n_kept,))
    g_out = np.lib.format.open_memmap(tmp / "gmag.npy", mode="w+", dtype="f4", shape=(n_kept,))
    rec_dtype = np.dtype([("dec", "f8"), ("ra", "f8"), ("g", "f4")])
    pos = 0
    for k in range(n_buckets):
        if counts[k] == 0:
            continue
        rec = np.fromfile(scratch / f"b{k:04d}.bin", dtype=rec_dtype)
        rec.sort(order="dec", kind="stable")
        n = rec.size
        dec_out[pos:pos + n] = rec["dec"]
        ra_out[pos:pos + n] = rec["ra"]
        g_out[pos:pos + n] = rec["g"]
        pos += n
        del rec
    if pos != n_kept:
        raise RuntimeError(f"bucket sort wrote {pos} rows, expected {n_kept}")
    dec_min, dec_max = float(dec_out[0]), float(dec_out[-1])
    dec_out.flush(); ra_out.flush(); g_out.flush()
    del dec_out, ra_out, g_out

    (tmp / _META).write_text(json.dumps({
        "source": str(src_npy), "n_rows": n_kept, "n_source_rows": int(n_total),
        "gmag_limit": gmag_limit, "sorted_by": "dec",
        "dec_min": dec_min, "dec_max": dec_max,
    }, indent=2))
    shutil.rmtree(scratch)

    if out_dir.exists():
        shutil.rmtree(out_dir)
    tmp.rename(out_dir)
    log.info("cache written: %s (%d rows)", out_dir, n_kept)
    return out_dir


class GaiaCatalog:
    """
    Query interface over the Gaia reference catalogue.

    Parameters
    ----------
    cache_dir : Path, optional
        Declination cache from :func:`build_dec_cache`.  Used when present.
    src_npy : Path, optional
        Full catalogue, used when no cache is available.
    chunk : int
        Chunk size for the fallback scan.

    Raises
    ------
    FileNotFoundError
        If neither a usable cache nor the source catalogue exists.
    """

    def __init__(self, cache_dir: Optional[Path] = None,
                 src_npy: Optional[Path] = None, chunk: int = 20_000_000):
        self.chunk = int(chunk)
        self._dec = self._ra = self._g = None
        self.mode = "none"
        self.src_npy = Path(src_npy) if src_npy else None

        if cache_dir is not None and self._cache_ok(Path(cache_dir)):
            cache_dir = Path(cache_dir)
            self._dec = np.load(cache_dir / "dec.npy", mmap_mode="r")
            self._ra = np.load(cache_dir / "ra.npy", mmap_mode="r")
            self._g = np.load(cache_dir / "gmag.npy", mmap_mode="r")
            self.mode = "cache"
            log.info("Gaia: declination cache %s (%d rows)", cache_dir, len(self._dec))
        elif self.src_npy is not None and self.src_npy.is_file():
            self.mode = "scan"
            log.info("Gaia: no cache; falling back to a chunked scan of %s "
                     "(run build_gaia_cache.py once to make this fast)", self.src_npy.name)
        else:
            raise FileNotFoundError(
                f"no Gaia catalogue available (cache_dir={cache_dir}, src_npy={src_npy})")

    @staticmethod
    def _cache_ok(d: Path) -> bool:
        return d.is_dir() and all((d / f).is_file() for f in _CACHE_FILES)

    # ------------------------------------------------------------------
    def cone_search(
        self,
        ra_deg: Sequence[float],
        dec_deg: Sequence[float],
        radius_deg: float,
        gmag_limit: Optional[float] = None,
    ) -> GaiaSubset:
        """
        Return every catalogue source within ``radius_deg`` of any pointing.

        Parameters
        ----------
        ra_deg, dec_deg : sequence of float
            Pointing centres [deg].  Non-finite entries are ignored.
        radius_deg : float
            Cone radius about each pointing [deg].
        gmag_limit : float, optional
            Keep only sources brighter than this.

        Returns
        -------
        GaiaSubset
            De-duplicated union over all pointings, so a source seen by two
            overlapping pointings appears once (review item R12).
        """
        ra_deg = np.asarray(ra_deg, dtype=np.float64)
        dec_deg = np.asarray(dec_deg, dtype=np.float64)
        ok = np.isfinite(ra_deg) & np.isfinite(dec_deg)
        if not ok.any():
            return GaiaSubset(np.empty(0), np.empty(0), np.empty(0))
        ra_deg, dec_deg = ra_deg[ok], dec_deg[ok]

        dec_lo = float(np.min(dec_deg) - radius_deg)
        dec_hi = float(np.max(dec_deg) + radius_deg)
        ra_band, dec_band, g_band = self._read_dec_band(dec_lo, dec_hi, gmag_limit)
        if ra_band.size == 0:
            return GaiaSubset(ra_band, dec_band, g_band)

        keep = self._cone_filter(ra_band, dec_band, ra_deg, dec_deg, radius_deg)
        return GaiaSubset(ra_band[keep], dec_band[keep], g_band[keep])

    # ------------------------------------------------------------------
    @staticmethod
    def _cone_filter(ra_cat, dec_cat, ra_pt, dec_pt, radius_deg) -> np.ndarray:
        """Indices of catalogue rows within ``radius_deg`` of any pointing."""
        from scipy.spatial import cKDTree

        # De-duplicate pointings on a grid one tenth of the search radius:
        # SPHEREx revisits the same field many times and a comet barely moves
        # between consecutive exposures, so this typically collapses hundreds of
        # pointings to a handful without changing the result.
        grid = max(float(radius_deg) / 10.0, 1e-6)
        key = np.column_stack((np.round(ra_pt / grid), np.round(dec_pt / grid)))
        _, uniq = np.unique(key, axis=0, return_index=True)
        ra_pt, dec_pt = ra_pt[uniq], dec_pt[uniq]

        tree = cKDTree(radec_to_unit(ra_cat, dec_cat))
        # Grid snapping can displace a pointing by up to half a cell; pad by a
        # full cell so no genuine neighbour is lost.
        r = _chord(float(radius_deg) + grid)
        hits = tree.query_ball_point(radec_to_unit(ra_pt, dec_pt), r=r)
        parts = [np.asarray(h, dtype=np.intp) for h in hits if len(h)]
        if not parts:
            return np.empty(0, dtype=np.intp)
        return np.unique(np.concatenate(parts))

    def _read_dec_band(self, dec_lo: float, dec_hi: float,
                       gmag_limit: Optional[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Read all sources in a declination band, by whichever backend exists."""
        if self.mode == "cache":
            i0 = int(np.searchsorted(self._dec, dec_lo, side="left"))
            i1 = int(np.searchsorted(self._dec, dec_hi, side="right"))
            dec = np.asarray(self._dec[i0:i1], dtype=np.float64)
            ra = np.asarray(self._ra[i0:i1], dtype=np.float64)
            g = np.asarray(self._g[i0:i1], dtype=np.float64)
        else:
            ra, dec, g = self._scan_dec_band(dec_lo, dec_hi)

        if gmag_limit is not None and g.size:
            keep = np.isfinite(g) & (g < float(gmag_limit))
            ra, dec, g = ra[keep], dec[keep], g[keep]
        return ra, dec, g

    def _scan_dec_band(self, dec_lo: float, dec_hi: float):
        """Fallback: one chunked pass over the unsorted source catalogue."""
        src = np.load(self.src_npy, mmap_mode="r")
        n = len(src)
        ra_p, dec_p, g_p = [], [], []
        for i0 in range(0, n, self.chunk):
            i1 = min(i0 + self.chunk, n)
            d = np.asarray(src["dec"][i0:i1], dtype=np.float64)
            keep = (d >= dec_lo) & (d <= dec_hi)
            if not keep.any():
                continue
            dec_p.append(d[keep])
            ra_p.append(np.asarray(src["ra"][i0:i1], dtype=np.float64)[keep])
            g_p.append(np.asarray(src["phot_g_mean_mag"][i0:i1], dtype=np.float64)[keep])
        if not ra_p:
            return np.empty(0), np.empty(0), np.empty(0)
        return np.concatenate(ra_p), np.concatenate(dec_p), np.concatenate(g_p)


class NeighborIndex:
    """
    Cone-search index over one target's Gaia subset.

    Built once per target and queried once per exposure.  Two things make this
    necessary rather than merely faster:

    * A target crossing the Galactic plane can pull in tens of thousands of
      sources over its whole track, while any single 91-pixel cutout sees a
      handful.  Projecting the full subset through the WCS for every exposure
      would be hundreds of times more work than needed.
    * ``WCS.all_world2pix`` inverts the SIP distortion iteratively, and that
      iteration *diverges* for points far outside the field of view.  Feeding it
      the whole subset raises ``NoConvergence`` and aborts the target.  Handing
      it only genuine neighbours keeps every point inside the regime where the
      inverse is well behaved.
    """

    __slots__ = ("subset", "_tree")

    def __init__(self, subset: "GaiaSubset"):
        self.subset = subset
        self._tree = None
        if len(subset):
            from scipy.spatial import cKDTree
            self._tree = cKDTree(radec_to_unit(subset.ra, subset.dec))

    def __len__(self) -> int:
        return len(self.subset)

    def query(self, ra_deg: float, dec_deg: float, radius_deg: float) -> "GaiaSubset":
        """
        Return the sources within ``radius_deg`` of one pointing.

        Parameters
        ----------
        ra_deg, dec_deg : float
            Cone centre [deg].
        radius_deg : float
            Cone radius [deg].

        Returns
        -------
        GaiaSubset
            Empty if the index is empty or the centre is not finite.
        """
        empty = GaiaSubset(np.empty(0), np.empty(0), np.empty(0))
        if self._tree is None or not (np.isfinite(ra_deg) and np.isfinite(dec_deg)):
            return empty
        pt = radec_to_unit(np.array([ra_deg]), np.array([dec_deg]))[0]
        idx = self._tree.query_ball_point(pt, r=_chord(radius_deg))
        if not idx:
            return empty
        idx = np.asarray(idx, dtype=np.intp)
        return GaiaSubset(self.subset.ra[idx], self.subset.dec[idx], self.subset.gmag[idx])
