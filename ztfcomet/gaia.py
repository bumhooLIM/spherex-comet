"""Gaia DR3 cone searches, for flagging background-source contamination.

A comet moves against the star field, so on some frames a background star sits
inside the photometric aperture and inflates the measured flux — and therefore
Afrho.  Nothing in the image itself distinguishes that from a genuine brightening,
which makes it exactly the kind of artefact that survives into a lightcurve and
gets interpreted as activity.

This module answers "is anything bright enough inside the aperture?" from a
local Gaia DR3 catalogue of every source with ``G < 18.5``.

Catalogue layout
----------------
The source catalogue is ``gaiadr3_all.npy``: 419 445 290 rows of
``(ra, dec, phot_g_mean_mag)``, 11 GB, and **unsorted** — so a cone search
against it means scanning the whole file.

Alongside it, ``gaiadr3_deccache/`` holds the same data split into three
contiguous arrays (``ra.npy``, ``dec.npy``, ``gmag.npy``) **sorted by
declination**.  That turns a cone search into a ``searchsorted`` for the
declination band followed by a filter, which is why this module prefers it:
measured, a 30-arcsec cone takes 3-6 ms and touches ~40 000 rows instead of
419 million.  Everything is memory-mapped, so nothing large is ever read into
memory.

If only the monolithic file is present, :meth:`GaiaCatalog.build_dec_cache`
creates the cache; failing that the search falls back to a chunked scan, which
works but is far too slow for a whole run.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "GaiaCatalog", "effective_magnitude", "flux_ratio", "ContaminationResult",
    "GAIA_MAG_LIMIT",
]

log = logging.getLogger(__name__)

#: Completeness limit of the bundled catalogue.  Sources fainter than this are
#: absent, so a "clean" verdict only means clean down to here.
GAIA_MAG_LIMIT = 18.5

_CHUNK = 8_000_000          # rows per block in the slow fallback scan


def effective_magnitude(gmags):
    """Combine source magnitudes into one total-flux magnitude.

    .. math:: G_\\mathrm{eff} = -2.5 \\log_{10} \\sum_i 10^{-0.4 G_i}

    Parameters
    ----------
    gmags : array_like
        Magnitudes of the sources falling inside the aperture.

    Returns
    -------
    float
        Combined magnitude, or ``inf`` when there are no sources — the natural
        limit, since zero flux is infinitely faint.

    Examples
    --------
    Two equal sources are 0.75 mag brighter together:

    >>> round(float(effective_magnitude([15.0, 15.0])), 2)
    14.25
    """
    gmags = np.asarray(gmags, dtype=float)
    gmags = gmags[np.isfinite(gmags)]
    if gmags.size == 0:
        return np.inf
    return float(-2.5 * np.log10(np.sum(10.0 ** (-0.4 * gmags))))


def flux_ratio(contaminant_mag, source_mag):
    """Flux of *contaminant_mag* relative to *source_mag*.

    Returns ``0.0`` when there is no contaminant and ``inf`` when the comet
    magnitude is unknown but a contaminant is present.
    """
    if not np.isfinite(contaminant_mag):
        return 0.0
    if not np.isfinite(source_mag):
        return np.inf
    return float(10.0 ** (-0.4 * (float(contaminant_mag) - float(source_mag))))


@dataclass(frozen=True)
class ContaminationResult:
    """Outcome of one aperture contamination check."""

    n_sources: int
    g_eff: float                 # combined magnitude of sources in the aperture
    g_brightest: float           # brightest single source
    ratio: float                 # contaminant flux / comet flux
    comet_mag: float             # the comparison magnitude used
    radius_arcsec: float         # search radius actually used
    contaminated: bool


class GaiaCatalog:
    """Memory-mapped Gaia DR3 catalogue with a declination-sorted fast path.

    Parameters
    ----------
    path : path-like, optional
        Either the ``gaiadr3_all.npy`` file or the directory containing it.
        Defaults to :data:`ztfcomet.directory.GAIA_ROOT`.
    cache_dir : path-like, optional
        The declination-sorted cache.  Defaults to ``gaiadr3_deccache`` beside
        the source file.

    Notes
    -----
    Opening is lazy and costs nothing: the arrays are memory-mapped on first
    use, so constructing a catalogue that is never queried reads no data.
    """

    def __init__(self, path=None, cache_dir=None):
        if path is None:
            from . import directory as d
            path = d.GAIA_ROOT

        path = Path(path).expanduser()
        if path.is_dir():
            self.root = path
            self.source_path = path / "gaiadr3_all.npy"
        else:
            self.root = path.parent
            self.source_path = path

        self.cache_dir = (Path(cache_dir).expanduser() if cache_dir
                          else self.root / "gaiadr3_deccache")
        self._dec = self._ra = self._gmag = None
        self._source = None
        self._warned_slow = False

    # ------------------------------------------------------------------ status
    @property
    def has_cache(self) -> bool:
        """True when the declination-sorted cache is present and usable."""
        meta = self.cache_dir / "meta.json"
        if not meta.is_file():
            return False
        try:
            if json.loads(meta.read_text()).get("sorted_by") != "dec":
                return False
        except (OSError, ValueError):
            return False
        return all((self.cache_dir / f"{n}.npy").is_file() for n in ("ra", "dec", "gmag"))

    @property
    def available(self) -> bool:
        """True when a cone search can be served at all."""
        return self.has_cache or self.source_path.is_file()

    def describe(self) -> str:
        """One-line summary, for notebook sanity checks."""
        if self.has_cache:
            return f"Gaia DR3: dec-sorted cache at {self.cache_dir} (fast)"
        if self.source_path.is_file():
            return (f"Gaia DR3: {self.source_path} — no dec cache, cone searches "
                    f"will scan the whole file (SLOW; call build_dec_cache())")
        return f"Gaia DR3: NOT FOUND at {self.root} — contamination flagging disabled"

    # -------------------------------------------------------------------- open
    def _open_cache(self):
        if self._dec is None:
            self._dec = np.load(self.cache_dir / "dec.npy", mmap_mode="r")
            self._ra = np.load(self.cache_dir / "ra.npy", mmap_mode="r")
            self._gmag = np.load(self.cache_dir / "gmag.npy", mmap_mode="r")
        return self._ra, self._dec, self._gmag

    def _open_source(self):
        if self._source is None:
            self._source = np.load(self.source_path, mmap_mode="r")
        return self._source

    # ------------------------------------------------------------------ search
    def cone_search(self, ra, dec, radius_arcsec, mag_limit=None):
        """Sources within *radius_arcsec* of ``(ra, dec)``.

        Parameters
        ----------
        ra, dec : float
            Centre in degrees (ICRS).
        radius_arcsec : float
            Search radius in arcseconds.
        mag_limit : float, optional
            Discard sources fainter than this.  The catalogue itself already
            stops at ``G = 18.5``.

        Returns
        -------
        gmag, sep_arcsec : ndarray
            Magnitudes and angular separations, sorted brightest first.

        Notes
        -----
        Separations use the haversine formula rather than a flat-sky
        approximation, so the result stays correct near the poles where
        ``cos(dec)`` collapses.
        """
        if not self.available:
            return np.empty(0), np.empty(0)

        radius_deg = float(radius_arcsec) / 3600.0
        if self.has_cache:
            cat_ra, cat_dec, cat_g = self._cone_from_cache(ra, dec, radius_deg)
        else:
            cat_ra, cat_dec, cat_g = self._cone_from_source(ra, dec, radius_deg)

        if cat_ra.size == 0:
            return np.empty(0), np.empty(0)

        sep = _angular_separation(ra, dec, cat_ra, cat_dec)
        keep = sep <= radius_deg
        if mag_limit is not None:
            keep &= cat_g <= float(mag_limit)

        gmag, sep = cat_g[keep], sep[keep] * 3600.0
        order = np.argsort(gmag)
        return gmag[order], sep[order]

    def _cone_from_cache(self, ra, dec, radius_deg):
        """Declination band via binary search, then a flat pre-filter in RA."""
        cat_ra, cat_dec, cat_g = self._open_cache()

        lo = int(np.searchsorted(cat_dec, dec - radius_deg, side="left"))
        hi = int(np.searchsorted(cat_dec, dec + radius_deg, side="right"))
        if hi <= lo:
            return np.empty(0), np.empty(0), np.empty(0)

        band_dec = np.asarray(cat_dec[lo:hi])
        band_ra = np.asarray(cat_ra[lo:hi])
        band_g = np.asarray(cat_g[lo:hi])

        # Cheap RA pre-filter before the exact separation.  Guard the pole,
        # where the RA window widens without bound.
        cos_dec = np.cos(np.radians(np.clip(abs(dec) + radius_deg, 0.0, 89.999)))
        if cos_dec > 1e-6:
            half_width = min(radius_deg / cos_dec, 180.0)
            delta = np.abs((band_ra - ra + 180.0) % 360.0 - 180.0)
            keep = delta <= half_width
            band_ra, band_dec, band_g = band_ra[keep], band_dec[keep], band_g[keep]
        return band_ra, band_dec, band_g

    def _cone_from_source(self, ra, dec, radius_deg):
        """Chunked scan of the unsorted 11 GB file.  Correct, but very slow."""
        if not self._warned_slow:
            log.warning("No Gaia dec cache at %s — scanning the full catalogue "
                        "per frame. Run GaiaCatalog.build_dec_cache() once.",
                        self.cache_dir)
            self._warned_slow = True

        data = self._open_source()
        ras, decs, gs = [], [], []
        for start in range(0, len(data), _CHUNK):
            block = data[start:start + _CHUNK]
            block_dec = np.asarray(block["dec"])
            keep = np.abs(block_dec - dec) <= radius_deg
            if not keep.any():
                continue
            ras.append(np.asarray(block["ra"])[keep])
            decs.append(block_dec[keep])
            gs.append(np.asarray(block["phot_g_mean_mag"])[keep])

        if not ras:
            return np.empty(0), np.empty(0), np.empty(0)
        return np.concatenate(ras), np.concatenate(decs), np.concatenate(gs)

    # ------------------------------------------------------------------- check
    def check_aperture(self, ra, dec, radius_arcsec, comet_mag,
                       flux_ratio_threshold=0.30, mag_limit=None):
        """Test one aperture for background-source contamination.

        Sums the flux of every catalogued source inside the aperture into a
        single effective magnitude and compares it with the comet's expected
        brightness.

        Parameters
        ----------
        ra, dec : float
            Aperture centre, degrees.
        radius_arcsec : float
            Aperture radius **already widened by the seeing FWHM**, so that a
            star just outside the aperture whose PSF spills into it still counts.
        comet_mag : float
            Expected total magnitude of the comet, from the JPL ephemeris.
        flux_ratio_threshold : float
            Flag when contaminating flux reaches this fraction of the comet's.
            0.30 corresponds to ``G_eff <= comet_mag + 1.31``.

        Returns
        -------
        ContaminationResult
        """
        gmag, _ = self.cone_search(ra, dec, radius_arcsec, mag_limit=mag_limit)
        g_eff = effective_magnitude(gmag)
        ratio = flux_ratio(g_eff, comet_mag)
        return ContaminationResult(
            n_sources=int(gmag.size),
            g_eff=g_eff,
            g_brightest=float(gmag[0]) if gmag.size else np.inf,
            ratio=ratio,
            comet_mag=float(comet_mag) if comet_mag is not None else np.nan,
            radius_arcsec=float(radius_arcsec),
            contaminated=bool(ratio >= flux_ratio_threshold),
        )

    # ------------------------------------------------------------------- build
    def build_dec_cache(self, overwrite=False, chunk=_CHUNK):
        """Create the declination-sorted cache from the monolithic catalogue.

        Writes ``ra.npy``, ``dec.npy``, ``gmag.npy`` and ``meta.json`` into
        :attr:`cache_dir`.  Needs roughly the size of the source file again in
        free disk space and takes a few minutes for 419 million rows.
        """
        if self.has_cache and not overwrite:
            log.info("Dec cache already present at %s", self.cache_dir)
            return self.cache_dir
        if not self.source_path.is_file():
            raise FileNotFoundError(self.source_path)

        data = self._open_source()
        n = len(data)
        log.info("Sorting %d Gaia rows by declination -> %s", n, self.cache_dir)

        order = np.argsort(np.asarray(data["dec"]), kind="stable")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        for name, dtype in (("ra", "<f8"), ("dec", "<f8"), ("gmag", "<f4")):
            field = {"ra": "ra", "dec": "dec", "gmag": "phot_g_mean_mag"}[name]
            out = np.lib.format.open_memmap(
                self.cache_dir / f"{name}.npy", mode="w+", dtype=dtype, shape=(n,))
            for start in range(0, n, chunk):
                stop = min(start + chunk, n)
                out[start:stop] = np.asarray(data[field])[order[start:stop]]
            out.flush()
            del out

        (self.cache_dir / "meta.json").write_text(json.dumps({
            "source": str(self.source_path), "n_rows": int(n), "sorted_by": "dec",
        }, indent=2))
        self._dec = self._ra = self._gmag = None
        log.info("Dec cache built at %s", self.cache_dir)
        return self.cache_dir


def _angular_separation(ra1, dec1, ra2, dec2):
    """Great-circle separation in degrees (haversine, pole-safe)."""
    phi1, phi2 = np.radians(dec1), np.radians(dec2)
    dphi = phi2 - phi1
    dlam = np.radians(np.asarray(ra2) - ra1)
    hav = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0))))
