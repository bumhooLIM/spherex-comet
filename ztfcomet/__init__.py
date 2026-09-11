"""ztfcomet — query and analyse ZTF observations of comets.

Beyond query and Af-rho photometry, :mod:`ztfcomet.profile` compares the
comet's radial surface-brightness profile with field stars on the same frame
-- the test of whether the coma is extended and whether it falls as 1/rho.

Combines what used to be two separate projects: ``ztfssoquery`` (cutout FITS
retrieval from IRSA) and the ``ztf-comet`` notebooks (Af-rho photometry and
figures).

Typical use::

    import ztfcomet as zc

    target = zc.get_target("24P")
    eph, frames, report = zc.search_frames(target)
    urls = zc.build_urls(frames, cutout_size=target.query.cutout_size)
    zc.download_urls(urls, zc.directory.data_dir(target.name))
    phot = zc.run_photometry(target)
    zc.plot_afrho({target.name: phot}, filters=["ZTF_r"], x="rh")

``notebooks/main.py`` runs exactly that chain for one or many targets.

Paths come from :mod:`ztfcomet.directory`; nothing else should build one.
Per-target constants live in :mod:`ztfcomet.config`.

See ``doc/primitive_code_analysis.md`` for the review that motivated this
structure, in particular the photometric corrections in :mod:`ztfcomet.phot`.
"""

from __future__ import annotations

import logging

from . import (config, cutout, directory, gaia, horizons, orbit, phot, plotting,
               profile, query, rcparams)
from .config import (
    PhotConfig, QueryConfig, Target, TARGETS, SOLAR_APPMAG_AB, get_target,
)
from .cutout import (build_urls, choose_cutout_size, construct_fitsurl,
                     download_urls, save_urls, verify_downloads)
from .directory import (
    DATA_ROOT, FIG_ROOT, GAIA_ROOT, PROJECT_ROOT, RESULT_ROOT, SPHEREX_ROOT, SUBJECTS,
    fig_kind_path, fig_path, photometry_path, profile_paths, result_path,
    data_dir, fig_dir, result_dir, target_slug,
)
from .gaia import GaiaCatalog, effective_magnitude
from .horizons import (
    is_fragment_designation, parse_ambiguity_table, resolve_record,
    resolve_target_id, verify_targetname,
)
from .phot import (
    ADVISORY_FLAGS, CRITICAL_FLAGS, FLAG_COLUMNS, attach_ephemerides,
    aperture_scale_ok, build_frame_table, calibrate, compute_afrho,
    flag_contamination, measure_photometry, run_multi_aperture, run_photometry,
)
from .profile import (ProfileConfig, fit_powerlaw, radial_profile, run_profiles,
                      select_field_stars, stack_star_profiles)
from .orbit import fetch_elements, time_from_perihelion
from .plotting import (plot_afrho, plot_afrho_apertures, plot_afrho_vs_rh,
                       plot_cutout, plot_cutout_grid, save_all_cutouts, signed_rh)
from .query import (
    extract_lastrecnum, query_sso_ephemeris, query_ztf_metadata, search_frames,
)

__version__ = "0.3.0"

__all__ = [
    "__version__",
    # submodules
    "config", "cutout", "directory", "gaia", "horizons", "orbit", "phot", "plotting",
    "profile", "query", "rcparams",
    # config
    "Target", "QueryConfig", "PhotConfig", "TARGETS", "get_target", "SOLAR_APPMAG_AB",
    # directory
    "PROJECT_ROOT", "DATA_ROOT", "RESULT_ROOT", "FIG_ROOT", "GAIA_ROOT", "SPHEREX_ROOT", "SUBJECTS",
    "result_path", "fig_path", "fig_kind_path", "photometry_path", "profile_paths",
    "data_dir", "result_dir", "fig_dir", "target_slug",
    # query
    "search_frames", "query_sso_ephemeris", "query_ztf_metadata", "extract_lastrecnum",
    # cutout
    "construct_fitsurl", "build_urls", "save_urls", "download_urls",
    "choose_cutout_size", "verify_downloads",
    # phot
    "build_frame_table", "attach_ephemerides", "measure_photometry",
    "calibrate", "compute_afrho", "run_photometry", "run_multi_aperture",
    "flag_contamination", "aperture_scale_ok",
    "FLAG_COLUMNS", "CRITICAL_FLAGS", "ADVISORY_FLAGS",
    # profile
    "ProfileConfig", "radial_profile", "select_field_stars", "stack_star_profiles",
    "fit_powerlaw", "run_profiles",
    # orbit
    "fetch_elements", "time_from_perihelion",
    # gaia
    "GaiaCatalog", "effective_magnitude",
    # horizons
    "resolve_record", "resolve_target_id", "verify_targetname",
    "parse_ambiguity_table", "is_fragment_designation",
    # plotting
    "plot_cutout", "plot_cutout_grid", "plot_afrho", "plot_afrho_vs_rh",
    "plot_afrho_apertures", "signed_rh", "save_all_cutouts",
    # helpers
    "setup_logging",
]


def setup_logging(level=logging.INFO):
    """Send this package's log records to the console.

    The pipeline reports dropped epochs, corrupt downloads and flagged frames
    through :mod:`logging`; without this call those messages are invisible.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)-8s %(name)s: %(message)s"))
    package_logger = logging.getLogger(__name__)
    package_logger.handlers[:] = [handler]
    package_logger.setLevel(level)
    package_logger.propagate = False
    return package_logger
