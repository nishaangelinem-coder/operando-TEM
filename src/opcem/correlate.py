"""Estimators that turn synchronised streams into structure-activity statements.

Implements the seven estimators of docs/03_proposal.md section 3.7.  Two rules
are enforced in code rather than left to the analyst:

*   a lead-lag peak is never labelled significant when its magnitude is below
    ``lag_significance_sigma_mult * sigma_disp`` (gap G1);
*   an event-aligned effect is always reported next to a surrogate null built
    from the same data (gap G2).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.preprocessing import StandardScaler

from .config import Config, DEFAULT
from .transport import shift_correct, tikhonov_deconvolve


# --------------------------------------------------------------------------- #
# stream alignment
# --------------------------------------------------------------------------- #
def rate_from_ms(cfg: Config, t_ms: np.ndarray, ion_current: np.ndarray,
                 k_cal: float, background: float, kernel: np.ndarray,
                 tau_dead: float, method: str = "deconvolve",
                 bin_s: Optional[float] = None) -> Dict[str, object]:
    """Convert a raw MS channel into an estimated formation-rate trace.

    ``method`` is one of ``"deconvolve"`` (Tikhonov, GCV-selected),
    ``"shift"`` (delay-only) or ``"raw"`` (no correction at all -- included so
    that the cost of skipping the correction can be measured, study S4).
    """
    from .transport import bin_trace, gamma_kernel, moments
    conc = (np.asarray(ion_current, dtype=float) - background) / k_cal
    bs = cfg.analysis.bin_s if bin_s is None else bin_s
    if bs is not None and bs > 0:
        # Bin to the analysis grid before inverting.  The kernel dispersion
        # (~0.8 s) is the resolution limit, so nothing is lost, the noise
        # drops by sqrt(10), and the dense (exactly cross-validated) solver
        # becomes affordable on the full record.
        dt_raw = float(np.mean(np.diff(t_ms)))
        t_ms, conc = bin_trace(t_ms, conc, bs)
        kt = np.arange(len(kernel)) * dt_raw
        m, sd = moments(kt, kernel)
        # re-express the same kernel on the coarse grid from its moments
        k_shape = (m / sd) ** 2 if sd > 0 else 4.0
        theta = sd ** 2 / m if m > 0 else bs
        tau0 = max(m - k_shape * theta, 0.0)
        kernel = gamma_kernel(np.arange(0.0, max(40.0, 6 * m), bs),
                              tau0, k_shape, theta)
    dt = float(np.mean(np.diff(t_ms)))
    if method == "raw":
        return {"t": t_ms, "rate": conc, "method": method, "lam": None}
    if method == "shift":
        return {"t": t_ms, "rate": shift_correct(conc, dt, tau_dead),
                "method": method, "lam": None}
    if method == "deconvolve":
        out = tikhonov_deconvolve(conc, kernel)
        return {"t": t_ms, "rate": out["x"], "method": method,
                "lam": out["lam"], "gcv": out.get("gcv"),
                "lam_grid": out.get("lam_grid")}
    raise ValueError(f"unknown method {method!r}")


def resample(t_src: np.ndarray, y: np.ndarray, t_dst: np.ndarray
             ) -> np.ndarray:
    return np.interp(t_dst, t_src, y)


# --------------------------------------------------------------------------- #
# event-aligned analysis
# --------------------------------------------------------------------------- #
@dataclass
class EventEffect:
    delta: float
    ci_lo: float
    ci_hi: float
    n_events: int
    p_surrogate: float
    null_mean: float
    null_sd: float
    z: float
    per_event: np.ndarray
    detrended: bool
    baseline: float

    def as_row(self) -> Dict[str, float]:
        d = {k: v for k, v in asdict(self).items() if k != "per_event"}
        d["delta_rel"] = (self.delta / self.baseline
                          if self.baseline else np.nan)
        return d


def _window_means(t: np.ndarray, y: np.ndarray, t_ev: float,
                  pre_s: float, post_s: float, detrend: bool,
                  gap_s: float = 0.0
                  ) -> Optional[Tuple[float, float]]:
    """Pre- and post-window means around an event, excluding a central gap.

    ``gap_s`` blanks an interval either side of the event.  Its purpose is to
    separate a perturbation's *direct* effect from its *structure-mediated*
    one: with the gap set to the actuation duration plus a few times the
    transport dispersion, the post-window is measured after the gas
    composition has returned to baseline, so any residual rate change has to
    be carried by the catalyst's state rather than by its coverages.
    """
    pre = (t >= t_ev - gap_s - pre_s) & (t < t_ev - gap_s)
    post = (t > t_ev + gap_s) & (t <= t_ev + gap_s + post_s)
    if pre.sum() < 2 or post.sum() < 2:
        return None
    if detrend:
        # remove a local linear trend fitted on the PRE window only, so that
        # a slow global drift (e.g. progressive sintering) cannot be read as
        # an event effect.  Fitting on both windows would absorb the effect.
        tp, yp = t[pre], y[pre]
        if len(tp) < 2:
            return None
        sl, ic = np.polyfit(tp - t_ev, yp, 1)
        y_pre = float(np.mean(yp - (sl * (tp - t_ev) + ic)))
        tq = t[post]
        y_post = float(np.mean(y[post] - (sl * (tq - t_ev) + ic)))
        return y_pre, y_post
    return float(np.mean(y[pre])), float(np.mean(y[post]))


def event_aligned_effect(cfg: Config, t: np.ndarray, y: np.ndarray,
                         event_times: Sequence[float],
                         pre_s: Optional[float] = None,
                         post_s: Optional[float] = None,
                         detrend: bool = True,
                         gap_s: float = 0.0,
                         n_boot: Optional[int] = None,
                         n_surrogate: Optional[int] = None,
                         seed: int = 0) -> EventEffect:
    """Event-related average of the rate trace, with a circular-shift null.

    The estimator is the mean over events of (post-window mean minus
    pre-window mean).  Uncertainty is a percentile bootstrap over events; the
    null is built by circularly shifting the rate trace relative to the event
    times, which preserves the trace's autocorrelation and the event times'
    clustering, and therefore does not manufacture significance the way an
    i.i.d. permutation null would.
    """
    a = cfg.analysis
    pre = a.pre_window_s if pre_s is None else pre_s
    post = a.post_window_s if post_s is None else post_s
    nb = a.n_bootstrap if n_boot is None else n_boot
    ns = a.n_surrogate if n_surrogate is None else n_surrogate
    rng = np.random.default_rng(seed)

    # enforce a guard interval so that overlapping events are not double counted
    ev = np.sort(np.asarray(event_times, dtype=float))
    kept: List[float] = []
    for e in ev:
        if kept and e - kept[-1] < a.event_guard_s:
            continue
        kept.append(float(e))
    ev = np.asarray(kept)

    def collect(shift: int = 0) -> np.ndarray:
        ys = np.roll(y, shift)
        vals = []
        for e in ev:
            w = _window_means(t, ys, e, pre, post, detrend, gap_s)
            if w is not None:
                vals.append(w[1] - w[0])
        return np.asarray(vals)

    per = collect(0)
    if per.size == 0:
        return EventEffect(np.nan, np.nan, np.nan, 0, np.nan, np.nan, np.nan,
                           np.nan, per, detrend, float(np.nanmean(y)))
    delta = float(np.mean(per))
    boot = np.array([np.mean(rng.choice(per, per.size, replace=True))
                     for _ in range(nb)])
    lo, hi = np.percentile(boot, [2.5, 97.5])

    n = len(y)
    min_shift = int((max(pre, post) + gap_s)
                    / max(float(np.mean(np.diff(t))), 1e-9)) + 1
    nulls = []
    for _ in range(ns):
        sh = int(rng.integers(min_shift, max(n - min_shift, min_shift + 1)))
        v = collect(sh)
        if v.size:
            nulls.append(float(np.mean(v)))
    nulls = np.asarray(nulls)
    if nulls.size:
        p = float((np.sum(np.abs(nulls) >= abs(delta)) + 1) / (nulls.size + 1))
        nm, nsd = float(nulls.mean()), float(nulls.std(ddof=1))
        z = (delta - nm) / nsd if nsd > 0 else np.nan
    else:
        p, nm, nsd, z = np.nan, np.nan, np.nan, np.nan
    return EventEffect(delta, float(lo), float(hi), int(per.size), p, nm, nsd,
                       float(z), per, detrend, float(np.nanmean(y)))


# --------------------------------------------------------------------------- #
# lead-lag
# --------------------------------------------------------------------------- #
@dataclass
class LeadLag:
    lags: np.ndarray
    rho: np.ndarray
    peak_lag: float
    peak_rho: float
    significant: bool
    reason: str


def lead_lag(cfg: Config, x: np.ndarray, y: np.ndarray, dt: float,
             max_lag_s: float = 20.0,
             sigma_disp: Optional[float] = None) -> LeadLag:
    """Lagged cross-correlation with the transport-resolution gate of gap G1.

    A peak whose |lag| is smaller than ``lag_significance_sigma_mult`` times
    the measured dispersion is NOT reported as significant, because the
    transport blur cannot resolve it.
    """
    sd = cfg.transport.sigma_disp() if sigma_disp is None else sigma_disp
    n_lag = int(round(max_lag_s / dt))
    xs = (x - np.nanmean(x)) / (np.nanstd(x) or 1.0)
    ys = (y - np.nanmean(y)) / (np.nanstd(y) or 1.0)
    xs = np.nan_to_num(xs); ys = np.nan_to_num(ys)
    lags = np.arange(-n_lag, n_lag + 1) * dt
    rho = np.empty(len(lags))
    n = len(xs)
    for i, k in enumerate(range(-n_lag, n_lag + 1)):
        if k >= 0:
            a, b = xs[:n - k], ys[k:]
        else:
            a, b = xs[-k:], ys[:n + k]
        rho[i] = float(np.dot(a, b) / len(a)) if len(a) > 2 else np.nan
    j = int(np.nanargmax(np.abs(rho)))
    pl, pr = float(lags[j]), float(rho[j])
    gate = cfg.analysis.lag_significance_sigma_mult * sd
    if abs(pl) < gate:
        return LeadLag(lags, rho, pl, pr, False,
                       f"|lag|={abs(pl):.2f}s < {gate:.2f}s "
                       f"({cfg.analysis.lag_significance_sigma_mult}x "
                       f"sigma_disp): unresolvable by transport")
    return LeadLag(lags, rho, pl, pr, True, "resolved")


# --------------------------------------------------------------------------- #
# multivariable models
# --------------------------------------------------------------------------- #
def fit_elastic_net(df: pd.DataFrame, y: np.ndarray, cols: Sequence[str],
                    force_cols: Sequence[str] = (), cfg: Config = DEFAULT,
                    seed: int = 0) -> Dict[str, object]:
    """Elastic net on standardised descriptors, with reactor terms forced in.

    Reactor variables (T, p_CO, p_O2) are always in the design matrix and are
    left unpenalised by being fitted first and the structural terms fitted to
    the residual.  Structural coefficients are therefore partial effects, not
    effects confounded with the operating point (gap G8).
    """
    use = [c for c in cols if c in df.columns]
    forced = [c for c in force_cols if c in df.columns]
    x = df[use].to_numpy(dtype=float)
    ok = np.isfinite(x).all(axis=1) & np.isfinite(y)
    x, yy = x[ok], np.asarray(y, dtype=float)[ok]
    if forced:
        z = df[forced].to_numpy(dtype=float)[ok]
        if np.ptp(z, axis=0).max() > 0:
            base = LinearRegression().fit(z, yy)
            yy = yy - base.predict(z)
    if len(yy) < 10:
        return {"r2": np.nan, "coef": {}, "n": len(yy), "alpha": np.nan}
    sc = StandardScaler().fit(x)
    xs = sc.transform(x)
    model = ElasticNetCV(l1_ratio=cfg.analysis.elastic_net_l1_ratio, cv=5,
                         random_state=seed, max_iter=20000)
    model.fit(xs, yy)
    pred = model.predict(xs)
    ss_res = float(np.sum((yy - pred) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    return {"r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan,
            "coef": {c: float(v) for c, v in zip(use, model.coef_)},
            "alpha": float(model.alpha_), "n": int(len(yy)),
            "forced": forced}


def two_stage_least_squares(z: np.ndarray, x: np.ndarray, y: np.ndarray
                            ) -> Dict[str, float]:
    """2SLS estimate of the effect of x on y, instrumented by z.

    Valid when the instrument z (here: a *randomised* actuator setting) affects
    y only through x.  Reported next to the observational coefficient so the
    two can be compared under a known confounder (study S4b).
    """
    z = np.asarray(z, dtype=float).reshape(-1, 1)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(z).all(axis=1) & np.isfinite(x) & np.isfinite(y)
    z, x, y = z[ok], x[ok], y[ok]
    if len(y) < 10:
        return {"beta_2sls": np.nan, "beta_ols": np.nan, "first_stage_r2":
                np.nan, "n": int(len(y))}
    s1 = LinearRegression().fit(z, x)
    x_hat = s1.predict(z)
    r2_1 = float(s1.score(z, x))
    s2 = LinearRegression().fit(x_hat.reshape(-1, 1), y)
    ols = LinearRegression().fit(x.reshape(-1, 1), y)
    # analytic standard error of the 2SLS slope
    resid = y - s2.predict(x_hat.reshape(-1, 1))
    sxx = float(np.sum((x_hat - x_hat.mean()) ** 2))
    se = float(np.sqrt(np.sum(resid ** 2) / max(len(y) - 2, 1) / max(sxx,
                                                                    1e-12)))
    return {"beta_2sls": float(s2.coef_[0]), "se_2sls": se,
            "beta_ols": float(ols.coef_[0]), "first_stage_r2": r2_1,
            "n": int(len(y))}


# --------------------------------------------------------------------------- #
# state flux (hidden Markov model)
# --------------------------------------------------------------------------- #
def hmm_state_flux(df: pd.DataFrame, cols: Sequence[str], n_states: int = 4,
                   n_iter: int = 60, seed: int = 0
                   ) -> Dict[str, object]:
    """Gaussian HMM on the descriptor trajectory, fitted by EM.

    Implemented directly (diagonal covariances, Baum-Welch) to avoid a
    dependency; the point is not the inference machinery but that the
    *transition matrix* becomes an estimand, so that catalyst dynamics are
    described by fluxes between motif regimes rather than by occupancies
    alone (gap G6).
    """
    use = [c for c in cols if c in df.columns]
    x = df[use].to_numpy(dtype=float)
    ok = np.isfinite(x).all(axis=1)
    x = x[ok]
    if len(x) < 4 * n_states:
        return {"ok": False}
    x = (x - x.mean(0)) / (x.std(0) + 1e-12)
    n, d = x.shape
    rng = np.random.default_rng(seed)
    mu = x[rng.choice(n, n_states, replace=False)].copy()
    var = np.ones((n_states, d))
    pi = np.full(n_states, 1.0 / n_states)
    a = np.full((n_states, n_states), 1.0 / n_states)

    def log_emis() -> np.ndarray:
        d2 = ((x[:, None, :] - mu[None, :, :]) ** 2) / var[None, :, :]
        return -0.5 * (d2 + np.log(2 * np.pi * var[None, :, :])).sum(axis=2)

    ll_old = -np.inf
    for _ in range(n_iter):
        le = log_emis()
        # forward-backward in log space
        la = np.empty((n, n_states)); lb = np.zeros((n, n_states))
        la[0] = np.log(pi + 1e-300) + le[0]
        log_a = np.log(a + 1e-300)
        for k in range(1, n):
            la[k] = le[k] + _logsumexp(la[k - 1][:, None] + log_a, axis=0)
        for k in range(n - 2, -1, -1):
            lb[k] = _logsumexp(log_a + (le[k + 1] + lb[k + 1])[None, :],
                               axis=1)
        ll = float(_logsumexp(la[-1], axis=0))
        lg = la + lb
        lg -= _logsumexp(lg, axis=1)[:, None]
        g = np.exp(lg)
        xi = np.zeros((n_states, n_states))
        for k in range(n - 1):
            m = (la[k][:, None] + log_a
                 + (le[k + 1] + lb[k + 1])[None, :])
            m -= _logsumexp(m.reshape(-1), axis=0)
            xi += np.exp(m)
        pi = g[0] / g[0].sum()
        a = xi / np.maximum(xi.sum(axis=1, keepdims=True), 1e-300)
        w = g.sum(axis=0)
        mu = (g.T @ x) / np.maximum(w[:, None], 1e-300)
        var = np.maximum((g.T @ (x ** 2)) / np.maximum(w[:, None], 1e-300)
                         - mu ** 2, 1e-4)
        if abs(ll - ll_old) < 1e-4 * abs(ll_old if ll_old else 1.0):
            break
        ll_old = ll
    states = np.argmax(g, axis=1)
    return {"ok": True, "states": states, "trans": a, "occupancy":
            g.mean(axis=0), "loglik": ll_old, "cols": use, "mask": ok,
            "mu": mu, "dwell": 1.0 / np.maximum(1.0 - np.diag(a), 1e-9)}


def _logsumexp(x: np.ndarray, axis: int) -> np.ndarray:
    m = np.max(x, axis=axis, keepdims=True)
    out = m + np.log(np.sum(np.exp(x - m), axis=axis, keepdims=True))
    return np.squeeze(out, axis=axis)


def state_rate_contribution(states: np.ndarray, rate: np.ndarray
                            ) -> pd.DataFrame:
    """Mean rate while each HMM state is occupied, with a standard error."""
    rows = []
    for s in np.unique(states):
        m = states == s
        r = rate[m]
        rows.append({"state": int(s), "n": int(m.sum()),
                     "occupancy": float(m.mean()),
                     "mean_rate": float(np.nanmean(r)),
                     "sem_rate": float(np.nanstd(r, ddof=1)
                                       / np.sqrt(max(m.sum(), 1)))})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# power
# --------------------------------------------------------------------------- #
def detection_power(cfg: Config, t: np.ndarray, y: np.ndarray,
                    event_times: Sequence[float],
                    n_grid: Sequence[int], n_rep: int = 40,
                    seed: int = 0) -> pd.DataFrame:
    """Empirical power of the event-aligned test vs number of events.

    For each event count, events are subsampled at random ``n_rep`` times and
    the fraction of subsamples whose bootstrap CI excludes zero is reported.
    """
    rng = np.random.default_rng(seed)
    ev = np.asarray(event_times, dtype=float)
    rows = []
    for n in n_grid:
        if n > len(ev):
            continue
        hits = 0
        deltas = []
        for r in range(n_rep):
            sel = rng.choice(ev, size=n, replace=False)
            e = event_aligned_effect(cfg, t, y, sel, n_boot=400,
                                     n_surrogate=0, seed=int(rng.integers(1e6)))
            deltas.append(e.delta)
            if np.isfinite(e.ci_lo) and np.isfinite(e.ci_hi) and \
                    (e.ci_lo > 0 or e.ci_hi < 0):
                hits += 1
        rows.append({"n_events": int(n), "power": hits / n_rep,
                     "mean_delta": float(np.nanmean(deltas)),
                     "sd_delta": float(np.nanstd(deltas, ddof=1))})
    return pd.DataFrame(rows)


def ols_natural_units(df: pd.DataFrame, y: np.ndarray, cols: Sequence[str],
                      control_cols: Sequence[str] = ()) -> Dict[str, object]:
    """Unpenalised OLS in the descriptors' own units, with controls included.

    The elastic net standardises its design matrix, so its coefficients are in
    units of per-standard-deviation and cannot be compared with a known partial
    derivative.  Where a target exists, the coefficient has to be estimated in
    natural units and without shrinkage, which is what this does.  Controls are
    included in the same design matrix rather than partialled out first, so the
    reported coefficients are partial effects in the usual regression sense.
    """
    use = [c for c in cols if c in df.columns]
    ctl = [c for c in control_cols if c in df.columns]
    x = df[use + ctl].to_numpy(dtype=float)
    yy = np.asarray(y, dtype=float)
    ok = np.isfinite(x).all(axis=1) & np.isfinite(yy)
    x, yy = x[ok], yy[ok]
    # drop control columns that do not vary, which would make X singular
    keep = [i for i in range(x.shape[1]) if np.ptp(x[:, i]) > 0]
    names = [(use + ctl)[i] for i in keep]
    x = x[:, keep]
    if len(yy) < len(names) + 5:
        return {"coef": {}, "r2": np.nan, "n": int(len(yy))}
    model = LinearRegression().fit(x, yy)
    return {"coef": {n: float(v) for n, v in zip(names, model.coef_)},
            "intercept": float(model.intercept_),
            "r2": float(model.score(x, yy)), "n": int(len(yy)),
            "controls": ctl}
