"""
Layer 6: production rates from a continuum-subtracted spectrum.

Ported from ``emission-fitter/fitting.py`` with one addition: an explicit
statement of which *flux space* the solve runs in.

At fixed geometry and in the optically thin regime the model is strictly
linear in each Q, ``F_nu(lambda_i) = sum_X Q_X A_iX``, so retrieval is a
weighted linear least-squares problem with an analytic covariance.  Column X of
the design matrix is the mJy one molecule per second of species X produces in
channel i *at that channel's own heliocentric and observer distance*.

Distance correction and Q
-------------------------
The continuum may have been subtracted in distance-corrected space,
``F_dc = F x r_h^2 Delta^2``.  Q is a property of the comet, not of the
observer, so it must come out the same whichever space the continuum was fitted
in.  Two equivalent routes exist and both are implemented:

``space="physical"`` (default)
    divide the emission back by each channel's own ``distcorr_factor`` and fit
    with the standard design matrix;
``space="distcorr"``
    fit the distance-corrected emission with the design matrix multiplied
    row-wise by the same factor.

Because the factor is applied per channel and the model is linear, the two give
identical Q to machine precision (``tests/test_comspec.py`` asserts this).  What
the choice of space *does* change is upstream -- the continuum polynomial, and
therefore ``emis`` itself -- and that is what the distance-correction study in
:mod:`analysis` measures.  There is no separate "re-correction" step: dividing
by the per-channel factor before the design matrix *is* the re-correction, and
it is exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from .config import H2O_HOT_RANGE, KEY_RANGES, SPECIES, FitConfig, ModelParams
from .gasmodel import amplitude_per_unit_Q, hires_grid, species_shape_mjy, spectrum_mjy, swings_factor
from .instrument import bandpass_matrix

__all__ = ["build_design_matrix", "fit_production_rates", "FitResult", "model_curves",
           "key_range_counts", "channel_masks", "data_covariance", "STATUS_DETECTED", "STATUS_MARGINAL",
           "STATUS_UPPER_LIMIT", "STATUS_NEGATIVE", "STATUS_NOT_COVERED", "STATUS_REJECTED"]

STATUS_DETECTED = "detected"        # >= FitConfig.detection_sigma
STATUS_MARGINAL = "marginal"        # between marginal_sigma and detection_sigma: value reported, limit quoted
STATUS_UPPER_LIMIT = "upper_limit"
STATUS_NEGATIVE = "negative_fit"
STATUS_NOT_COVERED = "not_covered"
STATUS_REJECTED = "rejected"        # a covered, fitted species whose detection the review rejected (case revision)

_COLS = {"physical": ("emis_raw_mjy", "emis_raw_err_mjy"),
         "distcorr": ("emis_mjy", "emis_err_mjy")}


def channel_masks(points: pd.DataFrame, cfg: "FitConfig | None" = None,
                  space: str = "physical") -> dict:
    """
    Which channels enter the fit, and which are excluded and why.

    Single source of truth for the input filtering, so a figure marks exactly
    the channels the solve dropped.
    """
    cfg = cfg or FitConfig()
    ycol, ecol = _COLS[space]
    y = points[ycol].to_numpy(dtype=float)
    sig = points[ecol].to_numpy(dtype=float)
    finite = np.isfinite(y) & np.isfinite(sig) & (sig > 0)
    clipped = np.zeros(len(y), dtype=bool)
    if cfg.clip_sigma is not None:
        clipped = finite & (np.abs(y) > cfg.clip_sigma * sig)
    negative = np.zeros(len(y), dtype=bool)
    if cfg.drop_negative_sigma is not None:
        negative = finite & ~clipped & (y < -cfg.drop_negative_sigma * sig)
    return dict(used=finite & ~clipped & ~negative, negative=negative,
                clipped=clipped, nonfinite=~finite)


def _velocities(points: pd.DataFrame) -> np.ndarray:
    """Heliocentric radial velocity per channel [km/s]; NaN (-> Swings factor 1) when absent."""
    if "v_hel_kms" in points:
        return points["v_hel_kms"].to_numpy(dtype=float)
    return np.full(len(points), np.nan)


def data_covariance(pts: pd.DataFrame, cfg: "FitConfig", space: str = "physical"):
    """
    Full covariance of the emission channels for generalised least squares.

    Diagonal: the photometric variance of every channel (``err_raw``, physical; times the
    distance-correction factor in the corrected space).  Off-diagonal, within each band: the
    continuum-model covariance ``V C V^T`` -- ``C`` the saved covariance of the polynomial
    coefficients (``cont_cov_ij``), ``V`` the Vandermonde rows ``(lambda - lambda_ref)^k`` --
    converted from the space the continuum was fitted in (``flux_space``) to the space of the
    solve.  Its diagonal equals the continuum error the diagonal path already used, so GLS
    changes only how the channels of a band are correlated.  Returns ``None`` when GLS is off or
    the points carry no continuum covariance (then the diagonal ``emis_*_err`` are used).
    """
    if not cfg.gls or "cont_cov_00" not in pts or "err_raw" not in pts or len(pts) == 0:
        return None
    n = len(pts)
    f = pts["distcorr_factor"].to_numpy(dtype=float)
    phot = pts["err_raw"].to_numpy(dtype=float)
    if space == "distcorr":
        phot = phot * f
    if not np.all(np.isfinite(phot)) or np.any(phot <= 0):
        return None
    S = np.diag(phot ** 2)
    wl = pts["wl"].to_numpy(dtype=float)
    bands = pts["band"].to_numpy()
    fs = pts["flux_space"].astype(str).to_numpy() if "flux_space" in pts else np.full(n, "physical")
    K = 4
    for b in pd.unique(bands):
        idx = np.flatnonzero(bands == b)
        row = pts.iloc[idx[0]]
        lam_ref = float(row.get("cont_lam_ref_um", np.nan))
        if not np.isfinite(lam_ref):
            continue
        C = np.zeros((K, K))
        for i in range(K):
            for j in range(i, K):
                v = row.get(f"cont_cov_{i}{j}", np.nan)
                C[i, j] = C[j, i] = 0.0 if not np.isfinite(v) else float(v)
        V = np.stack([(wl[idx] - lam_ref) ** k for k in range(K)], axis=1)
        block = V @ C @ V.T
        cont_space = fs[idx[0]]
        if cont_space == "distcorr" and space == "physical":
            block = block / np.outer(f[idx], f[idx])
        elif cont_space == "physical" and space == "distcorr":
            block = block * np.outer(f[idx], f[idx])
        S[np.ix_(idx, idx)] += block
    return S


def key_range_counts(points: pd.DataFrame) -> dict:
    wl = points["wl"].to_numpy(dtype=float) if len(points) else np.empty(0)
    return {s: int(((wl >= r["lo"]) & (wl <= r["hi"])).sum()) for s, r in KEY_RANGES.items()}


def build_design_matrix(points: pd.DataFrame, params: ModelParams,
                        species: tuple = SPECIES, space: str = "physical"):
    """
    Model flux per unit production rate for every channel and species.

    ``A[i, k]`` is the mJy produced in channel *i* by ``Q_k = 1`` molecule per
    second, at the channel's own ``r_hel`` / ``r_obs`` and through its own
    bandpass.  In ``space="distcorr"`` each row is multiplied by that channel's
    ``distcorr_factor`` so the model lives in the same space as the data.
    """
    lam_hi = hires_grid(params)
    W = bandpass_matrix(lam_hi, points["wl"].to_numpy(), points["wlwidth"].to_numpy())
    r_h = points["r_hel"].to_numpy(dtype=float)
    delta = points["r_obs"].to_numpy(dtype=float)
    v_h = _velocities(points)
    A = np.zeros((len(points), len(species)))
    for k, s in enumerate(species):
        shape = species_shape_mjy(s, lam_hi, params)
        A[:, k] = amplitude_per_unit_Q(s, r_h, delta, params, v_h) * (W @ shape)
    if space == "distcorr":
        A = A * points["distcorr_factor"].to_numpy(dtype=float)[:, None]
    return A, lam_hi


@dataclass
class FitResult:
    """Outcome of one production-rate fit."""

    target: str
    r_ap_km: float
    phase: int
    species: tuple
    Q_fit: np.ndarray
    Q: np.ndarray
    Q_err: np.ndarray
    Q_err_formal: np.ndarray
    cov: np.ndarray
    chi2: float
    dof: int
    n_points: int
    bands_used: tuple
    covered: np.ndarray
    n_key: np.ndarray
    n_eff: np.ndarray
    status: np.ndarray
    upper_limit: np.ndarray
    Q_limit: np.ndarray
    err_scale: float
    n_dropped_negative: int = 0
    space: str = "physical"
    geometry: dict = field(default_factory=dict)
    params: ModelParams | None = None
    config: FitConfig | None = None

    #: "main" (2.7 um band), "hot" (4.6-4.9 um hot bands only) or "none"
    h2o_source: str = "none"
    n_hot_H2O: int = 0
    #: the case revision applied to the group (``revisions.CaseRevision.describe()``), or ""
    revision: str = ""

    @property
    def chi2_red(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else np.nan

    @property
    def h2o_anchored(self) -> bool:
        """True when Q(H2O) rests on the 2.7 um main band; a hot-band value is not anchored."""
        return self.h2o_source == "main"

    def caveats(self) -> list:
        out = []
        if self.revision:
            out.append("case revision (review memo 2026-09-15): " + self.revision)
        if self.h2o_source == "hot":
            out.append(f"Q(H2O) from the 4.6-4.9 um hot bands only ({self.n_hot_H2O} channels; 2.7 um "
                       "not covered): provisional -- the hot bands carry ~3 % of the water emission "
                       "and the feature is shared with CO v(1-0)")
        elif (self.h2o_source == "none" and self.config is not None and self.config.h2o_hot_max_rh_au is not None
              and self.n_hot_H2O >= H2O_HOT_RANGE["min_points"]
              and self.geometry.get("r_hel_mean", 0.0) > self.config.h2o_hot_max_rh_au):
            out.append(f"Q(H2O) hot-band fallback withheld beyond {self.config.h2o_hot_max_rh_au:g} au "
                       f"({self.n_hot_H2O} hot-band channels present)")
        if ("CO" in self.species and self.covered[self.species.index("CO")]
                and self.params is not None and self.params.co_swings
                and not np.isfinite(self.geometry.get("v_hel_mean_kms", np.nan))):
            out.append("Q(CO) without the Swings factor: heliocentric velocity unknown, the "
                       "v_h = 0 g-factor was used (up to 25 % too small if |v_h| > 10 km/s)")
        for k, s in enumerate(self.species):
            if not self.covered[k]:
                r = KEY_RANGES.get(s)
                if r is not None:
                    out.append(f"Q({s}) not covered: {self.n_key[k]} channel(s) in "
                               f"{r['lo']:.2f}-{r['hi']:.2f} um, needs {r['min_points']}")
                else:
                    out.append(f"Q({s}) unconstrained: bands not covered by the data")
            elif self.status[k] == STATUS_REJECTED:
                out.append(f"Q({s}) rejected as a spurious detection by the review (case revision); "
                           f"the fitted value {self.Q_fit[k]:+.2e} +/- {self.Q_err[k]:.1e} is not a "
                           "production rate")
            elif self.status[k] == STATUS_NEGATIVE:
                out.append(f"Q({s}) best fit is negative ({self.Q_fit[k]:+.2e}): no production "
                           f"rate derived, upper limit only")
            elif self.status[k] == STATUS_MARGINAL:
                out.append(f"Q({s}) marginal ({self.Q_fit[k] / self.Q_err[k]:.1f} sigma): value reported, "
                           f"not a detection; limit {self.Q_limit[k]:.2e}")
            elif self.status[k] == STATUS_DETECTED and self.n_eff[k] < 2.0:
                out.append(f"Q({s}) rests on {self.n_eff[k]:.1f} effective channels: a formally "
                           "significant value carried by essentially one measurement")
        if ("H2O" in self.species and "CO" in self.species
                and not self.covered[self.species.index("H2O")]
                and self.covered[self.species.index("CO")]):
            out.append("Q(CO) fitted without H2O in the model: the 4.63 um hot band is unmodelled "
                       "and its flux is absorbed into Q(CO), biasing it high")
        if self.n_dropped_negative:
            out.append(f"{self.n_dropped_negative} channel(s) dropped for being more than "
                       f"1 sigma below zero (over-subtracted continuum)")
        if np.isfinite(self.chi2_red) and self.chi2_red > 10:
            out.append(f"chi2_red = {self.chi2_red:.0f}: model and data disagree beyond the "
                       "quoted errors; treat Q errors as lower bounds")
        return out

    def Q_dict(self) -> dict:
        return {s: float(q) for s, q in zip(self.species, self.Q)}

    def ratio(self, num: str, den: str = "H2O"):
        i, j = self.species.index(num), self.species.index(den)
        if self.status[i] != STATUS_DETECTED or self.status[j] != STATUS_DETECTED:
            return np.nan, np.nan
        a, b = self.Q[i], self.Q[j]
        if not np.isfinite(a) or not np.isfinite(b) or b == 0:
            return np.nan, np.nan
        r = a / b
        va, vb, cab = self.cov[i, i], self.cov[j, j], self.cov[i, j]
        var = r ** 2 * (va / a ** 2 + vb / b ** 2 - 2 * cab / (a * b)) if a != 0 else np.nan
        return float(r), float(np.sqrt(var)) if np.isfinite(var) and var > 0 else np.nan

    def to_row(self) -> dict:
        row = dict(target=self.target, r_ap_km=self.r_ap_km, phase=self.phase,
                   n_points=self.n_points, bands_used="+".join(self.bands_used),
                   h2o_anchored=self.h2o_anchored, h2o_source=self.h2o_source,
                   n_hot_H2O=self.n_hot_H2O, fit_space=self.space,
                   chi2=self.chi2, dof=self.dof, chi2_red=self.chi2_red,
                   err_scale=self.err_scale, **self.geometry)
        row["n_dropped_negative"] = self.n_dropped_negative
        for k, s in enumerate(self.species):
            row[f"Q_{s}_status"] = str(self.status[k])
            row[f"Q_{s}"] = self.Q[k]
            row[f"Q_{s}_fit"] = self.Q_fit[k]
            row[f"Q_{s}_err"] = self.Q_err[k]
            row[f"Q_{s}_err_formal"] = self.Q_err_formal[k]
            row[f"Q_{s}_nsig"] = (self.Q_fit[k] / self.Q_err[k]
                                  if np.isfinite(self.Q_fit[k]) and np.isfinite(self.Q_err[k]) and self.Q_err[k] > 0
                                  else np.nan)
            row[f"Q_{s}_covered"] = bool(self.covered[k])
            row[f"Q_{s}_n_key"] = int(self.n_key[k])
            row[f"Q_{s}_n_eff"] = self.n_eff[k]
            row[f"Q_{s}_is_upper_limit"] = bool(self.upper_limit[k])
            row[f"Q_{s}_upper_limit"] = self.Q_limit[k]
        for num in ("CO2", "CO"):
            if num in self.species and "H2O" in self.species:
                r, e = self.ratio(num, "H2O")
                row[f"{num}_H2O"] = r
                row[f"{num}_H2O_err"] = e
        row["caveats"] = "; ".join(self.caveats())
        row["revision"] = self.revision
        return row


def fit_production_rates(points: pd.DataFrame, params: ModelParams,
                         cfg: FitConfig | None = None, space: str = "physical",
                         rev=None) -> FitResult:
    """
    Weighted linear least squares for Q(H2O), Q(CO2), Q(CO).

    Parameters
    ----------
    points : DataFrame
        Continuum-subtracted channels carrying ``wl``, ``wlwidth``, ``r_hel``,
        ``r_obs``, ``distcorr_factor``, ``band``, and the emission columns of
        both spaces.
    params : ModelParams
    cfg : FitConfig, optional
    space : {"physical", "distcorr"}
        Flux space of the solve (see the module docstring).  The retrieved Q is
        the same in both; the default is physical.
    rev : revisions.CaseRevision, optional
        The group's case revision: a species can be declared covered on a minimum number
        of its band's emission channels instead of the ``KEY_RANGES`` rule, and a
        detection can be rejected (status ``rejected``, no value reported).

    Notes
    -----
    The solve is unconstrained -- a species consistent with zero can come back
    negative and is reported as such -- and the covariance is rescaled by
    ``chi2_red`` when it exceeds 1, the standard remedy for formal photometric
    errors that understate the true scatter.
    """
    cfg = cfg or FitConfig()
    species = tuple(cfg.species)
    ycol, ecol = _COLS[space]
    masks = channel_masks(points, cfg, space)
    good = masks["used"]
    n_dropped_neg = int(masks["negative"].sum())
    pts = points.loc[good].reset_index(drop=True)
    if pts.empty:
        # Every channel was masked -- in practice an accepted band whose continuum sits above
        # all of its emission channels.  That is a continuum failure, not a measurement of Q,
        # so the group is reported as not fitted rather than as a row of "not covered" species.
        who = (f"{points['target'].iloc[0]} phase {int(points['phase'].iloc[0])}"
               if len(points) else "empty input")
        detail = (f"{n_dropped_neg} of {len(points)} channels more than "
                  f"{cfg.drop_negative_sigma:g} sigma below zero"
                  if cfg.drop_negative_sigma is not None else
                  f"all {len(points)} channels non-finite or clipped")
        raise ValueError(f"{who}: no channel survived the channel masks ({detail})")
    y = pts[ycol].to_numpy(dtype=float)
    sig = pts[ecol].to_numpy(dtype=float)

    A, _ = build_design_matrix(pts, params, species, space)

    counts = key_range_counts(pts)
    n_key = np.array([counts.get(s, 0) for s in species], dtype=int)
    if cfg.require_key_coverage:
        covered = np.array([counts.get(s, 0) >= KEY_RANGES[s]["min_points"]
                            if s in KEY_RANGES else True for s in species])
    else:
        col_max = np.abs(A).max(axis=0) if len(A) else np.zeros(len(species))
        covered = col_max > (col_max.max() * 1e-6 if col_max.max() > 0 else np.inf)
    # case revision: coverage on the band's emission channels instead of the key range (the
    # review judged the channels present to carry the band)
    if rev is not None:
        from .revisions import BAND_OF
        for s, n_min in rev.coverage().items():
            if s in species and s in BAND_OF:
                n_band = int((pts["band"].astype(str) == BAND_OF[s]).sum())
                covered[species.index(s)] = n_band >= n_min
    # H2O: the 2.7 um main band anchors the fit; only when it is not covered do the 4.6-4.9 um
    # hot bands carry Q(H2O), and the result is labelled so downstream (bright, close comets
    # such as 10P and 24P lose the 2.7 um channels to saturation and flags).
    h2o_source, n_hot = "none", 0
    if "H2O" in species:
        k = species.index("H2O")
        wl = pts["wl"].to_numpy(dtype=float)
        hr = H2O_HOT_RANGE
        n_hot = int(((wl >= hr["lo"]) & (wl <= hr["hi"])).sum())
        n_red = int(((wl >= hr["red_lo"]) & (wl <= hr["hi"])).sum())
        rh_mean = float(np.mean(pts["r_hel"])) if len(pts) else np.nan
        inside_cap = cfg.h2o_hot_max_rh_au is None or rh_mean <= cfg.h2o_hot_max_rh_au
        if covered[k]:
            h2o_source = "main"
        elif (cfg.h2o_hot_fallback and cfg.require_key_coverage and inside_cap
              and n_hot >= hr["min_points"] and n_red >= hr["min_points_red"]):
            covered[k] = True
            h2o_source = "hot"
    if len(A):
        covered &= np.abs(A).max(axis=0) > 0
    if "H2O" in species and not covered[species.index("H2O")]:
        h2o_source = "none"

    n_free = int(covered.sum())
    Q_fit = np.full(len(species), np.nan)
    Q = np.full(len(species), np.nan)
    Q_err = np.full(len(species), np.nan)
    Q_err_formal = np.full(len(species), np.nan)
    cov_full = np.full((len(species), len(species)), np.nan)
    chi2, dof, scale = np.nan, 0, 1.0

    n_eff = np.zeros(len(species))
    if len(y):
        w = (A / sig[:, None]) ** 2
        tot = w.sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            n_eff = np.where(tot > 0, tot ** 2 / (w ** 2).sum(axis=0), 0.0)

    S = data_covariance(pts, cfg, space)
    gls_used = S is not None
    if n_free and len(y) > n_free:
        if gls_used:
            L = np.linalg.cholesky(S)
            Aw = solve_triangular(L, A[:, covered], lower=True)
            yw = solve_triangular(L, y, lower=True)
        else:
            Aw = A[:, covered] / sig[:, None]
            yw = y / sig
        scl = np.abs(Aw).max(axis=0)
        scl[scl == 0] = 1.0
        As = Aw / scl
        sol, *_ = np.linalg.lstsq(As, yw, rcond=None)
        XtX_inv = np.linalg.pinv(As.T @ As)
        q = sol / scl
        cov_s = XtX_inv / np.outer(scl, scl)
        resid = As @ sol - yw
        chi2 = float(resid @ resid)
        dof = int(len(y) - n_free)
        chi2_red = chi2 / dof if dof > 0 else np.nan
        scale = float(np.sqrt(chi2_red)) if (cfg.scale_errors_by_chi2
                                             and np.isfinite(chi2_red) and chi2_red > 1) else 1.0
        idx = np.where(covered)[0]
        Q_fit[idx] = q
        Q_err_formal[idx] = np.sqrt(np.diag(cov_s))
        Q_err[idx] = Q_err_formal[idx] * scale
        cov_scaled = cov_s * scale ** 2
        for a, ia in enumerate(idx):
            for b, ib in enumerate(idx):
                cov_full[ia, ib] = cov_scaled[a, b]

    status = np.full(len(species), STATUS_NOT_COVERED, dtype=object)
    Q_limit = np.full(len(species), np.nan)
    k = cfg.upper_limit_sigma
    for i in range(len(species)):
        if not covered[i] or not np.isfinite(Q_fit[i]) or not np.isfinite(Q_err[i]) or Q_err[i] <= 0:
            continue
        nsig = Q_fit[i] / Q_err[i]
        if Q_fit[i] < 0:
            status[i] = STATUS_NEGATIVE
            Q_limit[i] = k * Q_err[i]
        elif nsig >= cfg.detection_sigma:
            status[i] = STATUS_DETECTED
            Q[i] = Q_fit[i]
        elif nsig >= cfg.marginal_sigma:
            status[i] = STATUS_MARGINAL
            Q[i] = Q_fit[i]
            Q_limit[i] = Q_fit[i] + k * Q_err[i]
        else:
            status[i] = STATUS_UPPER_LIMIT
            Q[i] = Q_fit[i]
            Q_limit[i] = Q_fit[i] + k * Q_err[i]
    if rev is not None:
        for s in rev.reject:
            if s in species:
                i = species.index(s)
                if status[i] != STATUS_NOT_COVERED:
                    status[i] = STATUS_REJECTED
                    Q[i] = np.nan
                    Q_limit[i] = np.nan
    is_ul = np.isin(status, (STATUS_UPPER_LIMIT, STATUS_NEGATIVE))

    geom = dict(
        r_hel_mean=float(np.mean(pts["r_hel"])) if len(pts) else np.nan,
        r_hel_min=float(np.min(pts["r_hel"])) if len(pts) else np.nan,
        r_hel_max=float(np.max(pts["r_hel"])) if len(pts) else np.nan,
        r_obs_mean=float(np.mean(pts["r_obs"])) if len(pts) else np.nan,
        distcorr_factor_mean=float(np.mean(pts["distcorr_factor"])) if len(pts) else np.nan,
        v_g_kms=float(np.atleast_1d(params.v_g(np.mean(pts["r_hel"])))[0]) if len(pts) else np.nan,
        jd_utc_mean=float(np.mean(pts["jd_utc"])) if "jd_utc" in pts and len(pts) else np.nan,
    )
    geom["gls"] = bool(gls_used)
    v_h = _velocities(pts)
    geom["v_hel_mean_kms"] = float(np.nanmean(v_h)) if np.isfinite(v_h).any() else np.nan
    geom["swings_CO"] = float(np.mean(np.atleast_1d(swings_factor("CO", v_h, params)))) if len(v_h) else np.nan
    return FitResult(
        target=str(points["target"].iloc[0]) if len(points) else "",
        r_ap_km=float(points["r_ap_km"].iloc[0]) if len(points) else np.nan,
        phase=int(points["phase"].iloc[0]) if len(points) else -1,
        species=species, Q_fit=Q_fit, Q=Q, Q_err=Q_err, Q_err_formal=Q_err_formal, cov=cov_full,
        chi2=chi2, dof=dof, n_points=int(len(y)),
        bands_used=tuple(sorted(pts["band"].unique())) if len(pts) else (),
        covered=covered, n_key=n_key, n_eff=n_eff, status=status, upper_limit=is_ul,
        h2o_source=h2o_source, n_hot_H2O=n_hot,
        Q_limit=Q_limit, err_scale=scale, n_dropped_negative=n_dropped_neg, space=space,
        geometry=geom, params=params, config=cfg,
        revision=rev.describe() if rev is not None else "")


def model_curves(fit: FitResult, points: pd.DataFrame, params: ModelParams,
                 n_lam: int = 3000) -> dict:
    """Dense intrinsic and instrument-convolved model curves, plus the model at each channel."""
    Q = {s: (0.0 if not np.isfinite(v) else v) for s, v in fit.Q_dict().items()}
    r_h = float(np.median(points["r_hel"])) if len(points) else np.nan
    delta = float(np.median(points["r_obs"])) if len(points) else np.nan
    v_all = _velocities(points)
    v_h = float(np.nanmedian(v_all)) if np.isfinite(v_all).any() else np.nan
    lam = np.linspace(params.lam_min_um, params.lam_max_um, n_lam)
    comp = spectrum_mjy(lam, Q, r_h, delta, params, v_h)
    o = np.argsort(points["wl"].to_numpy())
    wl_o = points["wl"].to_numpy()[o]
    ww_o = points["wlwidth"].to_numpy()[o]
    width = np.interp(lam, wl_o, ww_o, left=ww_o[0], right=ww_o[-1]) if len(wl_o) else None
    convolved = np.full_like(lam, np.nan)
    if width is not None:
        lam_hi = hires_grid(params)
        hi = spectrum_mjy(lam_hi, Q, r_h, delta, params, v_h)["total"]
        convolved = bandpass_matrix(lam_hi, lam, width) @ hi
    A, _ = build_design_matrix(points, params, fit.species, "physical")
    at_points = A @ np.nan_to_num(fit.Q)
    return dict(lam=lam, intrinsic=comp["total"], convolved=convolved,
                per_species={s: comp[s] for s in fit.species},
                at_points=at_points, r_hel_repr=r_h, r_obs_repr=delta, v_hel_repr=v_h)
