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

from .config import KEY_RANGES, SPECIES, FitConfig, ModelParams
from .gasmodel import amplitude_per_unit_Q, hires_grid, species_shape_mjy, spectrum_mjy
from .instrument import bandpass_matrix

__all__ = ["build_design_matrix", "fit_production_rates", "FitResult", "model_curves",
           "key_range_counts", "channel_masks", "STATUS_DETECTED", "STATUS_UPPER_LIMIT",
           "STATUS_NEGATIVE", "STATUS_NOT_COVERED"]

STATUS_DETECTED = "detected"
STATUS_UPPER_LIMIT = "upper_limit"
STATUS_NEGATIVE = "negative_fit"
STATUS_NOT_COVERED = "not_covered"

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
    A = np.zeros((len(points), len(species)))
    for k, s in enumerate(species):
        shape = species_shape_mjy(s, lam_hi, params)
        A[:, k] = amplitude_per_unit_Q(s, r_h, delta, params) * (W @ shape)
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

    @property
    def chi2_red(self) -> float:
        return self.chi2 / self.dof if self.dof > 0 else np.nan

    @property
    def h2o_anchored(self) -> bool:
        """True when the 2.7 um band is in the fit; without it Q(H2O) and Q(CO) are degenerate."""
        return "2.7um" in self.bands_used

    def caveats(self) -> list:
        out = []
        if not self.h2o_anchored:
            out.append("Q(H2O) not anchored at 2.7 um: constrained only by the 4.6-4.9 um hot "
                       "bands, degenerate with Q(CO)")
        for k, s in enumerate(self.species):
            if not self.covered[k]:
                r = KEY_RANGES.get(s)
                if r is not None:
                    out.append(f"Q({s}) not covered: {self.n_key[k]} channel(s) in "
                               f"{r['lo']:.2f}-{r['hi']:.2f} um, needs {r['min_points']}")
                else:
                    out.append(f"Q({s}) unconstrained: bands not covered by the data")
            elif self.status[k] == STATUS_NEGATIVE:
                out.append(f"Q({s}) best fit is negative ({self.Q_fit[k]:+.2e}): no production "
                           f"rate derived, upper limit only")
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
                   h2o_anchored=self.h2o_anchored, fit_space=self.space,
                   chi2=self.chi2, dof=self.dof, chi2_red=self.chi2_red,
                   err_scale=self.err_scale, **self.geometry)
        row["n_dropped_negative"] = self.n_dropped_negative
        for k, s in enumerate(self.species):
            row[f"Q_{s}_status"] = str(self.status[k])
            row[f"Q_{s}"] = self.Q[k]
            row[f"Q_{s}_fit"] = self.Q_fit[k]
            row[f"Q_{s}_err"] = self.Q_err[k]
            row[f"Q_{s}_err_formal"] = self.Q_err_formal[k]
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
        return row


def fit_production_rates(points: pd.DataFrame, params: ModelParams,
                         cfg: FitConfig | None = None, space: str = "physical") -> FitResult:
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
    if len(A):
        covered &= np.abs(A).max(axis=0) > 0

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

    if n_free and len(y) > n_free:
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
        if not covered[i] or not np.isfinite(Q_fit[i]) or not np.isfinite(Q_err[i]):
            continue
        if Q_fit[i] < 0:
            status[i] = STATUS_NEGATIVE
            Q_limit[i] = k * Q_err[i]
        elif Q_fit[i] < k * Q_err[i]:
            status[i] = STATUS_UPPER_LIMIT
            Q[i] = Q_fit[i]
            Q_limit[i] = Q_fit[i] + k * Q_err[i]
        else:
            status[i] = STATUS_DETECTED
            Q[i] = Q_fit[i]
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
    return FitResult(
        target=str(points["target"].iloc[0]) if len(points) else "",
        r_ap_km=float(points["r_ap_km"].iloc[0]) if len(points) else np.nan,
        phase=int(points["phase"].iloc[0]) if len(points) else -1,
        species=species, Q_fit=Q_fit, Q=Q, Q_err=Q_err, Q_err_formal=Q_err_formal, cov=cov_full,
        chi2=chi2, dof=dof, n_points=int(len(y)),
        bands_used=tuple(sorted(pts["band"].unique())) if len(pts) else (),
        covered=covered, n_key=n_key, n_eff=n_eff, status=status, upper_limit=is_ul,
        Q_limit=Q_limit, err_scale=scale, n_dropped_negative=n_dropped_neg, space=space,
        geometry=geom, params=params, config=cfg)


def model_curves(fit: FitResult, points: pd.DataFrame, params: ModelParams,
                 n_lam: int = 3000) -> dict:
    """Dense intrinsic and instrument-convolved model curves, plus the model at each channel."""
    Q = {s: (0.0 if not np.isfinite(v) else v) for s, v in fit.Q_dict().items()}
    r_h = float(np.median(points["r_hel"])) if len(points) else np.nan
    delta = float(np.median(points["r_obs"])) if len(points) else np.nan
    lam = np.linspace(params.lam_min_um, params.lam_max_um, n_lam)
    comp = spectrum_mjy(lam, Q, r_h, delta, params)
    o = np.argsort(points["wl"].to_numpy())
    wl_o = points["wl"].to_numpy()[o]
    ww_o = points["wlwidth"].to_numpy()[o]
    width = np.interp(lam, wl_o, ww_o, left=ww_o[0], right=ww_o[-1]) if len(wl_o) else None
    convolved = np.full_like(lam, np.nan)
    if width is not None:
        lam_hi = hires_grid(params)
        hi = spectrum_mjy(lam_hi, Q, r_h, delta, params)["total"]
        convolved = bandpass_matrix(lam_hi, lam, width) @ hi
    A, _ = build_design_matrix(points, params, fit.species, "physical")
    at_points = A @ np.nan_to_num(fit.Q)
    return dict(lam=lam, intrinsic=comp["total"], convolved=convolved,
                per_species={s: comp[s] for s in fit.species},
                at_points=at_points, r_hel_repr=r_h, r_obs_repr=delta)
