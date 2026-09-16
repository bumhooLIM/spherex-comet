# Legacy SPHEREx aperture-photometry prototype

`apphot.py` and `apphot_all.py` are the original single-file prototype, and
`apphot_sample_2P.ipynb`, `apphot_sample_UN271.ipynb` and `reflectance.ipynb` the
notebooks it grew out of (outputs stripped on 2026-09-15: 94 MB of embedded
figures removed).  They are kept only because `doc/apphot/code_review_primitive.md`
refers to them line by line.

**Do not run them.**  They cannot run as written -- they reference
`directory.APPHOT_DIR` and `directory.COMBFITS_DIR`, which never existed in
their `directory.py`, and import `skyloc` and `skimage`, which the project does
not depend on -- and where they do work they carry the defects catalogued in the
review: an annulus that falls outside the cutout for nearby targets, Gaia stars
masked out of the science aperture, negative fluxes discarded, and a
double-counted sky term in the uncertainties.

The working pipeline is the `spherex_apphot` package at the project root;
start from `spherex_apphot/README.md`.
