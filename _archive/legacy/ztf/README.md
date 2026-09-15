# Legacy notebooks

The pre-merge notebooks from `ztf-comet/notebooks/` and `ztf-sso-query/notebooks/`,
kept for provenance with their outputs stripped (4.9 MB of embedded PNGs removed).

**They do not run against this package** and are not maintained. They import
`ztfssoquery`, `_rcparams`, `sxobsplan` and `astrometry`, none of which are
dependencies of `ztfcomet`.

**Do not copy code out of them.** `afrho_240P.ipynb` in particular contains the
single-zeropoint bug (C1) that made every Af-rho value it produced wrong, and
`ztfquery_2P.ipynb` annotates its figures with 24P's ephemeris because the orbit
ID was pasted across from the 24P notebook (C9).

| notebook | superseded by |
|---|---|
| `ztfquery_24P.ipynb`, `ztfquery_2P.ipynb`, `ztfquery_240P.ipynb`, `ztfquery_2019Y3.ipynb` | `notebooks/query.ipynb` |
| `afrho_240P.ipynb` | `notebooks/afrho.ipynb`, `ztfcomet/phot.py` |
| the cutout-figure cells (duplicated 4x) | `ztfcomet.plot_cutout` |
| `240P.ipynb` | scratch work on 240P fragment separation; not ported |

See `doc/primitive_code_analysis.md` for the full review.
