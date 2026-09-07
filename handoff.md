# Handoff

## Current State
`ztfcomet` v0.3 — 68-comet survey **running unattended** (resumed 21:25 on
2026-09-07, `ZTFCOMET_DATA` pinned to `/Volumes/T7/data/ztf-comet`;
log `survey_log_20260907T212513.log`, ~18/68 at 23:20, ~13 min/target).

- Package: `directory`, `config`, `query`, `horizons`, `cutout`, `phot`, `gaia`,
  `orbit`, `profile`, `plotting`. 135 tests. Commits on `fix/horizons-query-correctness`.
- `notebooks/survey.py`: query → adaptive cutout (5′; 10′ if Tmag<14 or Δ<1 au) →
  download → re-verify on disk → multi-aperture Afρ (10/15/20/30/40k km, only where
  FWHM < r_ap < 1′) → radial profile → figures. Resumable via `survey_status.csv`.
- Afρ vs r_h figures: x = r_h − q signed by leg, perihelion marked, T−Tp axis on
  top (Kepler, validated to 0.01 d); clean frames only.
- `ztfcomet.profile`: comet SB profile vs ≤20 field stars (both ×4 oversampled),
  naive slope + PSF-convolved nucleus/coma model with optional free sky.
  `doc/profile_resolution_24P.md`: inbound steepness is sky error on faint frames
  (not ρ-range, not nucleus); oversampling changes nothing; model+sky gives −1.01;
  1/ρ holds to ≥12 000 km (lower bound).
- Data on T7: 24P 295, 2P 166 (Vmag<21 exception), plus 10P…210P; ~35 GB.

## Next Steps (in order, once `Survey complete` appears in the log)
1. Repair pass over targets with `complete == False`: `survey.py --steps download
   phot profile figures --no-resume` (404s are permanent; 500/504s retry).
2. `survey.py --steps profile figures --no-resume` — targets 1–12 ran on old code
   (signed-r_h plots, no profile step). Regenerates everything with current code.
3. 2P: `main.py 2P --steps phot profile figures` (its 166 files are Vmag<21 data;
   the survey marks it `no_frames` at Vmag<20). Note the exception in the summary.
4. Read `survey_summary.md`; report process + flag statistics. Then
   `profile_resolution.py --target 29P` etc.: does 1/ρ hold beyond 12 000 km?

## Blind Spots / Dead Ends
- **A spun-down SSD answers `is_dir()` False once**; a resume silently restarted
  on the local fallback. Resolver retries, survey refuses fallback; pin `ZTFCOMET_DATA`.
- **Horizons clips `targetname` at 31 chars** ("Bernardinelli-Bernstein (C/2014").
  `verify_targetname` prefix-matches clipped names; fragments still fail.
- **Horizons rejects an all-identical TLIST** ("Bad dates"): dedupe, join on JD.
  astroquery's ValueError is transient unless it carries an ambiguity *listing*.
- **Signed r_h is wrong** (unreachable band between −q and +q): use r_h − q; beyond
  aphelion is NaN, not half a period.
- **Native annulus means are biased for steep profiles**; oversample comet AND
  stars, normalise within 1.5 px. Oversampling adds no information.
- **Formal sky error (0.07 DN) understates the real systematic ~50×**; few-DN
  background structure tilts faint outer profiles. Fit sky as a free term.
- **A cache keyed on a rounded parameter gives least_squares a zero gradient.**
- **pandas parses "2024E1" as a float** in a homogeneous column: every `read_csv`
  with a target column needs `dtype={"target": str}` or resume reprocesses it.
- **Contamination flagging scales with aperture**: 30k/40k km rows mostly flagged
  on bright low-latitude targets.
- `sep.winpos` returns 3 values; IRSA sends 200+HTML (check FITS magic); bare `"2P"` → Styx: `id_type="smallbody"`.
