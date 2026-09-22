"""Beam-perturbation accounting: separating chemistry from irradiation (gap G4).

The observed rate of any atomic process is
``k_obs(phi) = k_chem + beta * phi``.  A dose-rate series at fixed (T, p,
flow) therefore identifies ``k_chem`` as the zero-dose intercept and ``beta``
as the slope, and every observed event can be assigned a

    BPI = beta * phi / (k_chem + beta * phi)

the *Beam Perturbation Index*: the probability that the event was caused by
the probe rather than by the reaction.  Events above the pre-registered BPI
threshold are excluded from causal claims.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .config import Config, DEFAULT


@dataclass
class DoseSeriesFit:
    k_chem: float
    k_chem_se: float
    beta: float
    beta_se: float
    r2: float
    n_points: int
    dose_rates: np.ndarray
    k_obs: np.ndarray
    k_obs_se: np.ndarray

    def bpi(self, dose_rate: float) -> float:
        num = self.beta * dose_rate
        den = self.k_chem + num
        return float(num / den) if den > 0 else np.nan

    def bpi_ci(self, dose_rate: float, n_draw: int = 4000,
               seed: int = 0) -> tuple:
        """CI for the BPI, propagating both fitted parameters."""
        rng = np.random.default_rng(seed)
        b = rng.normal(self.beta, max(self.beta_se, 1e-12), n_draw)
        c = rng.normal(self.k_chem, max(self.k_chem_se, 1e-12), n_draw)
        num = b * dose_rate
        den = c + num
        v = np.where(den > 0, num / den, np.nan)
        v = v[np.isfinite(v)]
        return (float(np.percentile(v, 2.5)),
                float(np.percentile(v, 97.5))) if v.size else (np.nan, np.nan)

    def as_row(self) -> Dict[str, float]:
        return {"k_chem": self.k_chem, "k_chem_se": self.k_chem_se,
                "beta": self.beta, "beta_se": self.beta_se, "r2": self.r2,
                "n_points": self.n_points}


def fit_dose_series(dose_rates: Sequence[float], k_obs: Sequence[float],
                    k_obs_se: Optional[Sequence[float]] = None
                    ) -> DoseSeriesFit:
    """Weighted linear fit of an observed process rate against dose rate."""
    phi = np.asarray(dose_rates, dtype=float)
    k = np.asarray(k_obs, dtype=float)
    ok = np.isfinite(phi) & np.isfinite(k)
    phi, k = phi[ok], k[ok]
    if k_obs_se is not None:
        se = np.asarray(k_obs_se, dtype=float)[ok]
        w = 1.0 / np.maximum(se, 1e-12) ** 2
    else:
        se = np.full_like(k, np.nan)
        w = np.ones_like(k)
    sw = w.sum()
    mx = float(np.dot(w, phi) / sw)
    my = float(np.dot(w, k) / sw)
    sxx = float(np.dot(w, (phi - mx) ** 2))
    sxy = float(np.dot(w, (phi - mx) * (k - my)))
    beta = sxy / sxx if sxx > 0 else np.nan
    inter = my - beta * mx
    resid = k - (inter + beta * phi)
    dof = max(len(k) - 2, 1)
    s2 = float(np.dot(w, resid ** 2) / dof)
    beta_se = float(np.sqrt(s2 / sxx)) if sxx > 0 else np.nan
    inter_se = float(np.sqrt(s2 * (1.0 / sw + mx ** 2 / sxx)))
    ss_tot = float(np.dot(w, (k - my) ** 2))
    r2 = 1.0 - float(np.dot(w, resid ** 2)) / ss_tot if ss_tot > 0 else np.nan
    return DoseSeriesFit(k_chem=float(inter), k_chem_se=inter_se,
                         beta=float(beta), beta_se=beta_se, r2=float(r2),
                         n_points=int(len(k)), dose_rates=phi, k_obs=k,
                         k_obs_se=se)


def annotate_events_with_bpi(events: pd.DataFrame, fit: DoseSeriesFit,
                             dose_rate: float, cfg: Config = DEFAULT
                             ) -> pd.DataFrame:
    """Attach the BPI and the pre-registered inclusion flag to every event."""
    out = events.copy()
    b = fit.bpi(dose_rate)
    out["bpi"] = b
    out["bpi_included"] = b <= cfg.analysis.bpi_exclusion_threshold
    out["dose_rate"] = dose_rate
    return out


def beam_free_extrapolation(dose_rates: Sequence[float],
                            values: Sequence[float],
                            value_se: Optional[Sequence[float]] = None
                            ) -> Dict[str, float]:
    """Zero-dose extrapolation of any dose-dependent observable.

    Used for quantities that are not rates (e.g. the dispersed fraction), for
    which the zero-dose intercept is the beam-free estimate the experiment is
    really after.
    """
    fit = fit_dose_series(dose_rates, values, value_se)
    return {"value_at_zero_dose": fit.k_chem, "se": fit.k_chem_se,
            "slope": fit.beta, "slope_se": fit.beta_se, "r2": fit.r2}


# --------------------------------------------------------------------------- #
# site-heterogeneity aware accounting
# --------------------------------------------------------------------------- #
def effective_hop_rate(cfg: Config, temp_k: float, dose_rate: float) -> float:
    """Site-occupancy-weighted effective hop rate -- what an MSD estimator sees.

    A diffusing monomer re-draws its site type after every hop, so it occupies
    site type ``s`` in proportion to ``p_s / k_s`` and the ensemble hop rate is
    the harmonic mean ``1 / sum_s (p_s / k_s)``, NOT the terrace rate.  Because
    a harmonic mean of ``k_s + beta*phi`` is not affine in ``phi``, the
    observable effective rate is **sublinear in dose**: the linear dose model
    that the Beam Perturbation Index assumes is exact only for a single
    barrier, and is biased whenever the support presents a distribution of
    them.  Quantifying that bias is one of the results of study S6.
    """
    from .truth import site_occupancy
    return site_occupancy(cfg, temp_k, dose_rate)[1]


def beam_hop_fraction(cfg: Config, temp_k: float, dose_rate: float) -> float:
    """Exact fraction of hop events caused by the beam, over all site types.

    The number of hops per unit time leaving site type ``s`` is
    ``occupancy_s * k_s``, and since occupancy is proportional to ``p_s / k_s``
    that product is proportional to the site abundance ``p_s``.  Hence

        BPI_true = sum_s p_s * (beta*phi) / (k_th,s + beta*phi)

    which is the quantity the simulator's own beam-attribution counter
    estimates, and the reference the fitted index is scored against.
    """
    from .truth import SITE_TYPES, k_hop
    s_cfg = cfg.support
    p = np.array([1.0 - s_cfg.frac_vacancy - s_cfg.frac_step,
                  s_cfg.frac_vacancy, s_cfg.frac_step])
    p = p / p.sum()
    k_th = np.array([k_hop(cfg, st, temp_k, 0.0)[0] for st in SITE_TYPES])
    kb = cfg.kinetics.beta_beam * dose_rate
    return float(np.dot(p, kb / (k_th + kb)))


def fit_beta_site_aware(cfg: Config, temp_k: float,
                        dose_rates: Sequence[float],
                        k_obs: Sequence[float]) -> Dict[str, float]:
    """Fit the beam coupling with the site distribution treated as known.

    In a real experiment the barrier distribution comes from a beam-off
    temperature series, so only the beam coupling ``beta`` is unknown.  Fitting
    the correct harmonic-mean model rather than a straight line removes the
    sublinearity bias identified in :func:`effective_hop_rate`.
    """
    from scipy import optimize
    from .truth import SITE_TYPES, k_hop
    s_cfg = cfg.support
    p = np.array([1.0 - s_cfg.frac_vacancy - s_cfg.frac_step,
                  s_cfg.frac_vacancy, s_cfg.frac_step])
    p = p / p.sum()
    k_th = np.array([k_hop(cfg, st, temp_k, 0.0)[0] for st in SITE_TYPES])
    phi = np.asarray(dose_rates, dtype=float)
    y = np.asarray(k_obs, dtype=float)

    def model(pars):
        beta, gain = pars
        keff = np.array([1.0 / np.sum(p / (k_th + beta * f)) for f in phi])
        return gain * keff

    fit = optimize.least_squares(lambda q: model(q) - y, [1e-4, 1.0],
                                 bounds=([0.0, 0.1], [1e-1, 10.0]))
    beta, gain = fit.x
    keff0 = 1.0 / float(np.sum(p / k_th))
    resid = model(fit.x) - y
    dof = max(len(y) - 2, 1)
    s2 = float(np.sum(resid ** 2) / dof)
    try:
        jac = fit.jac
        cov = s2 * np.linalg.inv(jac.T @ jac)
        beta_se = float(np.sqrt(max(cov[0, 0], 0.0)))
    except Exception:
        beta_se = np.nan
    return {"beta": float(beta), "beta_se": beta_se,
            "detection_gain": float(gain),
            "k_eff_zero_dose": keff0,
            "rmse": float(np.sqrt(np.mean(resid ** 2)))}
