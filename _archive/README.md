# `_archive/` — superseded material, kept for provenance only

Never read anything here for current numbers or copy code out of it.

| folder | content | why it is kept |
|---|---|---|
| `legacy/ztf/` | the pre-merge `ztf-comet` / `ztf-sso-query` notebooks, outputs stripped | `doc/ztf/primitive_code_analysis.md` cites them line by line; they carry the bugs the review documents |
| `legacy/apphot/` | the single-file SPHEREx photometry prototype (`apphot.py`, `apphot_all.py`) and the three notebooks it grew out of, outputs stripped (94 MB removed on 2026-09-15) | `doc/apphot/code_review_primitive.md` cites them line by line |

Deleted on 2026-09-15 and not recoverable from here: the previous photometry set
(`apphot_5502194856bc`, fixed 15–20 px annulus — its comparison products remain in
`results/comspec/studies/apphot_comparison/` and `annulus_previous/`), the duplicate copy of
the photometry the catalog used to read, and the per-exposure cutout PNGs.  A pre-merge
archive of all code and documents — including the `CLAUDE.md`, `README.md` and `handoff.md` of
the three former projects, from which the root documents were merged — is
`../_backup_pre-merge_2026-09-15.tar.gz`; the ZTF ones are also in the git history.
