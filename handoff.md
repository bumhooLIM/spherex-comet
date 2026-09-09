# Handoff

## Current State
The 68-comet survey is **complete** (2026-09-07 17:25 → 2026-09-08 22:32).
55 targets have photometry, 13 produced no frames or no epochs at Vmag < 20.
8,336 frames, 39,591 measurement rows across five apertures, 53.0% clean.
Process and flag statistics: `doc/survey_summary_68comets.md`.

Coma profiles confirm extended comae: median slope −1.23 vs −4.41 for field
stars, 501/524 clean frames shallower than the stars in the same image.

**The T7 is currently unmounted.** Every step below reads raw FITS from it, so
mount it first and confirm `ZTFCOMET_DATA=/Volumes/T7/data/ztf-comet` resolves
before starting anything. `survey.py` refuses the local fallback root only when
the env var is unset — with it set to a missing path it will not protect you.

## Next Steps
Run in this order, all from the project root with ZTFCOMET_DATA pinned:

1. `survey.py --steps download --no-resume` — repair pass. ~300 frames missing,
   131 of them 235P from the 09-07/08 IRSA outage. 404s are permanent; the 5xx
   are recoverable and IRSA has recovered. Re-check fields 000616 c02 q4 and
   000375 c10, which timed out repeatedly while neighbours served fine.
2. Re-query `2024E1` and `240P` — their windows were truncated by a stale
   `end_date` in `config.py`. 2024E1 is missing its 2026-01-20 perihelion and
   all post-perihelion data, so its Afρ figure currently shows one leg only.
3. `survey.py --steps profile figures --no-resume` — profiles exist for 9 of 55
   targets; this also regenerates the targets 1–12 figures still on signed-r_h.
4. `main.py 2P --steps phot profile figures` — 2P's 166 frames are Vmag < 21
   data; the survey scored it `no_frames` at Vmag < 20. Note the exception.

## Blind Spots / Dead Ends
- A **running process keeps its imported modules**: the profile step and the
  perihelion axis were both added mid-run and silently did not apply. Restart
  the batch after changing the package, or accept that only later targets get it.
- The run **hung 54 min on one Horizons call** at 0% CPU while the service
  answered fresh requests in 0.8 s. astroquery's 30 s timeout never fired.
  `socket.setdefaulttimeout(300)` in `survey.py` now bounds this.
- macOS writes `._name` AppleDouble sidecars on the exFAT T7. They matched the
  FITS glob and were reported as half the frames being unreadable. Filtered in
  `phot.py`; do not "fix" that by re-downloading.
- Contamination is the top rejection cause at 20.1% — expected along the
  ecliptic, not a bug.
