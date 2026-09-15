# Superseded prototype

`apphot.py` and `apphot_all.py` are the original single-file prototype.  They are
kept only because `doc/code_review_primitive.md` refers to them line by line.

**Do not run them.**  They cannot run as written -- they reference
`directory.APPHOT_DIR` and `directory.COMBFITS_DIR`, which never existed -- and
where they do work they carry the defects catalogued in the review: an annulus
that falls outside the cutout for nearby targets, Gaia stars masked out of the
science aperture, negative fluxes discarded, and a double-counted sky term in
the uncertainties.

The working pipeline is the `spherex_apphot` package next to this directory;
start from `spherex-apphot/README.md`.
