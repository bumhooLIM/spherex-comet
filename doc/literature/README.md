# `doc/literature/` — the papers behind the model

Read a paper into a lean `.summary.md` beside it rather than loading the full text into a
session (`CLAUDE.md`).  The PDFs are gitignored (copyrighted); the summaries are tracked.

**Status 2026-09-15.**  During the merge, 18 PDFs were removed from this folder by a git
history rewrite (see `handoff.md`); the folder is synced by Google Drive, so they can be
restored in one step from Drive's Trash (drive.google.com → Trash → Restore).  Two were
recovered from other copies meanwhile.  ✔ = present, ✘ = restore from Drive Trash or
re-download from ADS.

| folder | file | status |
|---|---|---|
| `dust-continuum/` | `Harker_2002_ApJ.pdf` (Harker et al. 2002, ApJ 580, 579) | ✘ |
| | `Harker_2007_Icarus.pdf` (Harker et al. 2007, Icarus 190, 432) | ✘ |
| | `Harker_2023_PSJ.pdf` (Harker et al. 2023, PSJ) | ✘ |
| | `Kelley_2016_PASP.pdf` (Kelley et al. 2016, PASP 128, 018009) | ✔ (copy from `literature-review/`) |
| `fluorescence-db/` | `Villanueva_2011_Icarus.pdf` (Villanueva et al. 2011, Icarus 216, 227) + `.summary.md` | ✔ (copy from Drive) |
| | `Villanueva_2011_JGR.pdf` + `.summary.md` | ✘ (summary present) |
| | `Villanueva_2011_JQSRT.pdf` + `.summary.md` | ✘ (summary present) |
| | `Villanueva_2012_ApJ.pdf` + `.summary.md` | ✘ (summary present) |
| | `Villanueva_2013_JQSRT.pdf` + `.summary.md` | ✘ (summary present) |
| `fluorescence-emission/` | `Crovisier_1983_AA.pdf` (Crovisier & Encrenaz 1983, A&A 126, 170) | ✘ |
| | `Debout_2016_Icarus.pdf` (Debout et al. 2016, Icarus 265, 110) | ✘ |
| | `Gicquel_2023_PSJ.pdf` (Gicquel et al. 2023, PSJ; the NEOWISE Q table in `data/reference/`) | ✘ |
| | `Harrington_Pinto_2022_PSJ.pdf` (Harrington Pinto et al. 2022, PSJ 3, 247; the Q table in `data/reference/`) | ✘ |
| | `Ootsubo_2012_ApJ.pdf` (Ootsubo et al. 2012, ApJ 752, 15 — the AKARI band technique) | ✘ |
| | `PSG_handbook.pdf` (Planetary Spectrum Generator handbook; another copy is in the user's Drive) | ✘ |
| | `Villanueva_2011_Icarus.pdf` (a second copy of the paper above) | ✘ (not needed) |
| | `Villanueva_2018_JQRST.pdf` (Villanueva et al. 2018, JQSRT 217, 86 — PSG; a `Villanueva_2017_JQSRT.pdf` copy is in the user's Drive) | ✘ |
| `ice-absorption/` | `Protopapa_2014_Icarus.pdf` (Protopapa et al. 2014, Icarus 238, 191) | ✘ |

`Kelley_2016_PASP.pdf` used to exist in both `dust-continuum/` and `fluorescence-emission/`;
one copy is enough.  Delete this status table once the files are back.
