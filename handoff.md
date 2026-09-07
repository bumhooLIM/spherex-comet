# Handoff

## Current State
`ztfcomet` v0.3 — 68-comet survey **running unattended** (resumed 21:25 on
2026-09-07, `ZTFCOMET_DATA` pinned to `/Volumes/T7/data/ztf-comet`;
log `survey_log_20260907T212513.log`, ~18/68 at 23:20, ~13 min/target).

- Package: `directory`, `config`, `query`, `horizons`, `cutout`, `phot`, `gaia`,
  `orbit`, `profile`, `plotting`. 126 tests.
- `notebooks/survey.py`: query → adaptive cutout (5′; 10′ if Tmag<14 or Δ<1 au) →
  download → re-verify on disk → multi-aperture Afρ (10/15/20/30/40k km, only where
  FWHM < r_ap < 1′) → radial profile → figures. Resumable via `survey_status.csv`.
- Afρ vs r_h figures: x = r_h − q signed by leg, perihelion marked, T−Tp axis on top
  (Kepler, validated to 0.01 d on 24P/10P); clean frames only.
- `ztfcomet.profile`: comet SB profile (×4 oversampled; clipped vs plain mean) vs
  ≤20 field stars; power-law slope. 24P clean median −1.07 (coma −1; stars −4.3).
- Data on T7: 24P 295, 2P 166 (Vmag<21 exception), plus 10P…210P; ~35 GB.

## Next Steps (in order, once `Survey complete` appears in the log)
1. Repair pass for transient failures: `survey.py --steps download phot profile
   figures --no-resume` over targets with `complete == False` in `survey_status.csv`
   (31 failed files so far: 15×404 permanent, 16×500/504 transient).
2. `survey.py --steps profile figures --no-resume` — targets 1–12 ran on old code
   (signed-r_h plots, no profile step). Regenerates everything with current code.
3. 2P: `main.py 2P --steps phot profile figures` (its 166 files are Vmag<21 data;
   the survey marks it `no_frames` at Vmag<20). Note the exception in the summary.
4. Read `survey_summary.md`; report process + flag statistics to the user.
5. All post-PR commits are on `fix/horizons-query-correctness` (PR #1 open).

## Blind Spots / Dead Ends
- **A spun-down SSD answers `is_dir()` False once.** The first resume silently ran
  on the local fallback and restarted from target 1. Resolver now retries; survey
  refuses the fallback root. Always pin `ZTFCOMET_DATA` for unattended runs.
- **Horizons clips `targetname` at 31 chars** ("Bernardinelli-Bernstein (C/2014").
  `verify_targetname` prefix-matches clipped names; fragments still fail.
- **Horizons rejects an all-identical TLIST** ("Bad dates"): dedupe epochs; every
  caller joins on JD. astroquery raises ValueError for transient failures too —
  only an ambiguity *listing* means "resolve a record".
- **Signed r_h is wrong**: an unreachable band between −q and +q. Use r_h − q.
  Beyond aphelion must be NaN, not saturated at half a period.
- **Native-pixel annulus means are biased for steep profiles** (they sample only
  the radii that exist); that is why the comet is oversampled. Normalise both comet
  and stars by the mean SB within 1.5 px, else stars lose their 0.5 px bin (<3 px).
- **Contamination flagging scales with aperture**: 30k/40k km rows are mostly
  flagged on bright low-latitude targets; expect thin clean samples there.
- `sep.winpos` returns 3 values; IRSA returns 200+HTML; check FITS magic; never
  hardcode Horizons record numbers; `"2P"` bare → Styx (905): pass `id_type="smallbody"`.
