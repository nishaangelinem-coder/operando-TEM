"""Gas-transport transfer function: modelling, calibration, and inversion.

This module implements the correction that gap G1 requires.  The mass
spectrometer does not see the outlet composition, it sees

.. math:: S_j(t) = (R_j * h)(t) + \\epsilon(t)

where ``h`` is the impulse response of reactor volume + heated capillary +
differentially pumped inlet.  Everything downstream of here assumes ``h`` has
been *measured*, not guessed, which is why :func:`calibrate_from_tracer` is
step 4 of the protocol in docs/04_experimental_setup.md and runs before the
primary dataset is acquired.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from scipy import optimize
from scipy.special import gammaln

from .config import Config, TransportConfig, DEFAULT


# --------------------------------------------------------------------------- #
# kernel
# --------------------------------------------------------------------------- #
def gamma_kernel(t: np.ndarray, tau0: float, k: float, theta: float
                 ) -> np.ndarray:
    """Gamma impulse response with a pure dead time, normalised to unit area.

    ``tau_dead = tau0 + k*theta`` (mean) and ``sigma_disp = theta*sqrt(k)``.
    The tanks-in-series / axial-dispersion form reproduces the
    dead-time-plus-tailing behaviour of a nanoreactor and transfer line.
    """
    u = np.asarray(t, dtype=float) - tau0
    out = np.zeros_like(u)
    m = u > 0
    out[m] = np.exp((k - 1.0) * np.log(u[m]) - u[m] / theta
                    - gammaln(k) - k * np.log(theta))
    s = out.sum()
    return out / s if s > 0 else out


def kernel_from_config(cfg: Config, dt: float, length_s: float = 30.0,
                       transport: Optional[TransportConfig] = None
                       ) -> np.ndarray:
    tr = transport or cfg.transport
    t = np.arange(0.0, length_s, dt)
    return gamma_kernel(t, tr.tau0_s, tr.gamma_k, tr.gamma_theta_s)


def moments(t: np.ndarray, h: np.ndarray) -> Tuple[float, float]:
    """Empirical (mean, standard deviation) of a discrete kernel."""
    w = h / h.sum()
    m = float(np.dot(w, t))
    v = float(np.dot(w, (t - m) ** 2))
    return m, float(np.sqrt(max(v, 0.0)))


# --------------------------------------------------------------------------- #
# tracer experiment
# --------------------------------------------------------------------------- #
@dataclass
class TracerResult:
    tau0_s: float
    gamma_k: float
    gamma_theta_s: float
    tau_dead_s: float
    sigma_disp_s: float
    tau_dead_ci: Tuple[float, float]
    sigma_disp_ci: Tuple[float, float]
    rmse: float
    n_pulses: int
    kernel: np.ndarray
    kernel_t: np.ndarray

    def as_row(self) -> Dict[str, float]:
        return {"tau_dead_s": self.tau_dead_s,
                "tau_dead_lo": self.tau_dead_ci[0],
                "tau_dead_hi": self.tau_dead_ci[1],
                "sigma_disp_s": self.sigma_disp_s,
                "sigma_disp_lo": self.sigma_disp_ci[0],
                "sigma_disp_hi": self.sigma_disp_ci[1],
                "rmse": self.rmse, "n_pulses": self.n_pulses}


def make_pulse_train(t: np.ndarray, starts: Sequence[float],
                     width_s: float = 1.0, amplitude: float = 1.0
                     ) -> np.ndarray:
    """Timestamped rectangular tracer pulses (the valve is on the master clock)."""
    u = np.zeros_like(t)
    for s in starts:
        u[(t >= s) & (t < s + width_s)] = amplitude
    return u


def calibrate_from_tracer(t: np.ndarray, u: np.ndarray, y: np.ndarray,
                          n_boot: int = 300, seed: int = 0) -> TracerResult:
    """Fit ``h`` from a measured tracer input/response pair.

    Non-linear least squares on (tau0, k, theta) plus a gain and an offset.
    Confidence intervals for the reported moments come from a residual
    bootstrap, which propagates the actual MS noise rather than assuming it.
    """
    dt = float(np.mean(np.diff(t)))
    n = len(t)
    kt = np.arange(0.0, min(40.0, n * dt), dt)

    def model(p: np.ndarray) -> np.ndarray:
        tau0, k, theta, gain, off = p
        h = gamma_kernel(kt, max(tau0, 0.0), max(k, 1.01), max(theta, 1e-3))
        return gain * np.convolve(u, h, mode="full")[:n] + off

    def resid(p: np.ndarray) -> np.ndarray:
        return model(p) - y

    y0 = float(np.median(y[:max(3, int(0.5 / dt))]))
    p0 = np.array([1.0, 2.0, 0.5, float(np.ptp(y)) / max(np.ptp(u), 1e-9), y0])
    bounds = ([0.0, 1.01, 1e-3, 0.0, -np.inf],
              [30.0, 50.0, 20.0, np.inf, np.inf])
    fit = optimize.least_squares(resid, p0, bounds=bounds, x_scale="jac")
    p = fit.x
    r = resid(p)
    rmse = float(np.sqrt(np.mean(r ** 2)))

    rng = np.random.default_rng(seed)
    taus, sigs = [], []
    for _ in range(n_boot):
        yb = model(p) + rng.choice(r, size=n, replace=True)
        try:
            fb = optimize.least_squares(lambda q: model(q) - yb, p,
                                        bounds=bounds, x_scale="jac",
                                        max_nfev=200)
            tau0, k, theta = fb.x[:3]
            taus.append(tau0 + k * theta)
            sigs.append(theta * np.sqrt(k))
        except Exception:
            continue
    taus, sigs = np.asarray(taus), np.asarray(sigs)
    tau_ci = (float(np.percentile(taus, 2.5)), float(np.percentile(taus, 97.5))) \
        if taus.size else (np.nan, np.nan)
    sig_ci = (float(np.percentile(sigs, 2.5)), float(np.percentile(sigs, 97.5))) \
        if sigs.size else (np.nan, np.nan)

    h = gamma_kernel(kt, p[0], p[1], p[2])
    return TracerResult(
        tau0_s=float(p[0]), gamma_k=float(p[1]), gamma_theta_s=float(p[2]),
        tau_dead_s=float(p[0] + p[1] * p[2]),
        sigma_disp_s=float(p[2] * np.sqrt(p[1])),
        tau_dead_ci=tau_ci, sigma_disp_ci=sig_ci, rmse=rmse,
        n_pulses=int(np.sum(np.diff((u > 0).astype(int)) > 0)),
        kernel=h, kernel_t=kt)


# --------------------------------------------------------------------------- #
# inversion
# --------------------------------------------------------------------------- #
def _conv_matrix(h: np.ndarray, n: int) -> np.ndarray:
    """Lower-triangular Toeplitz convolution operator, truncated to n x n."""
    m = min(len(h), n)
    a = np.zeros((n, n))
    for i in range(m):
        np.fill_diagonal(a[i:, :n - i], h[i])
    return a


def tikhonov_deconvolve(y: np.ndarray, h: np.ndarray,
                        lam: Optional[float] = None,
                        lam_grid: Optional[Sequence[float]] = None,
                        method: str = "auto") -> Dict[str, object]:
    """Second-difference-regularised deconvolution with a GCV-selected weight.

    Minimises ``||h*x - y||^2 + lam ||D2 x||^2``.  The default solver works in
    the frequency domain on a zero-padded signal, which is O(n log n) and so
    usable on the full 10 Hz MS record; ``method="dense"`` builds the explicit
    Toeplitz operator and is kept as an independently coded reference that the
    FFT path is checked against in the test suite.

    Generalised cross-validation selects ``lam``; the whole GCV curve is
    returned so that the regularisation choice is auditable rather than
    hidden.
    """
    n = len(y)
    if method == "auto":
        method = "dense" if n <= 2000 else "fft"
    if method == "dense":
        return _tikhonov_dense(y, h, lam, lam_grid)

    # Mirror-pad, do NOT zero-pad.  Zero padding asserts that the record
    # drops to zero just past its end; that artificial step dominates the
    # inversion and leaks back across the whole trace.  Mirroring keeps the
    # padded signal continuous at both wrap-around boundaries.
    m = int(2 ** np.ceil(np.log2(2 * n + len(h))))
    ypad = np.concatenate([y, y[::-1]])
    reps = int(np.ceil(m / len(ypad)))
    yp = np.tile(ypad, reps)[:m]
    hp = np.zeros(m); hp[:len(h)] = h
    d2 = np.zeros(m); d2[0], d2[1], d2[-1] = -2.0, 1.0, 1.0

    hf = np.fft.rfft(hp)
    yf = np.fft.rfft(yp)
    df = np.fft.rfft(d2)
    h2 = np.abs(hf) ** 2
    dd = np.abs(df) ** 2

    grid = list(lam_grid) if lam_grid is not None else \
        list(DEFAULT.analysis.tikhonov_lambda_grid)
    if lam is not None:
        grid = [lam]

    best, gcvs = None, []
    for l in grid:
        den = h2 + l * dd
        xf = np.conj(hf) * yf / np.maximum(den, 1e-30)
        x = np.fft.irfft(xf, n=m)[:n]
        resid = np.fft.irfft(hf * xf - yf, n=m)[:n]
        rss = float(np.sum(resid ** 2))
        # effective degrees of freedom of the smoother, restricted to the
        # real record length
        tr = float(np.sum(h2 / np.maximum(den, 1e-30))) * (n / m)
        g = rss / max((n - tr) ** 2, 1e-9)
        gcvs.append(g)
        if best is None or g < best[1]:
            best = (l, g, x)
    return {"x": best[2], "lam": best[0], "gcv": np.asarray(gcvs),
            "lam_grid": np.asarray(grid)}


def _tikhonov_dense(y: np.ndarray, h: np.ndarray, lam: Optional[float],
                    lam_grid: Optional[Sequence[float]]) -> Dict[str, object]:
    """Explicit-operator reference implementation (small n only)."""
    n = len(y)
    if n > 2000:
        raise ValueError("dense deconvolution is O(n^2) memory; use fft")
    a = _conv_matrix(h, n)
    d = np.zeros((n - 2, n))
    for i in range(n - 2):
        d[i, i:i + 3] = (1.0, -2.0, 1.0)
    ata, dtd = a.T @ a, d.T @ d
    aty = a.T @ y
    if lam is not None:
        return {"x": np.linalg.solve(ata + lam * dtd, aty), "lam": lam,
                "gcv": None, "lam_grid": None}
    grid = list(lam_grid) if lam_grid is not None else \
        list(DEFAULT.analysis.tikhonov_lambda_grid)
    best, gcvs = None, []
    for l in grid:
        mm = np.linalg.solve(ata + l * dtd, a.T)
        x = mm @ y
        trace = float(np.trace(a @ mm))
        rss = float(np.sum((a @ x - y) ** 2))
        g = rss / max((n - trace) ** 2, 1e-9)
        gcvs.append(g)
        if best is None or g < best[1]:
            best = (l, g, x)
    return {"x": best[2], "lam": best[0], "gcv": np.asarray(gcvs),
            "lam_grid": np.asarray(grid)}


def shift_correct(y: np.ndarray, dt: float, tau_dead: float) -> np.ndarray:
    """Delay-only correction: shift the MS trace earlier by tau_dead.

    This is the conservative alternative to deconvolution.  It removes the
    delay but not the dispersion, so it blurs event-aligned estimates without
    amplifying noise.  Both are reported side by side in the results.
    """
    k = int(round(tau_dead / dt))
    if k <= 0:
        return y.copy()
    return np.concatenate([y[k:], np.full(k, y[-1])])


def forward(rate: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Apply the transfer function (used by the virtual instruments)."""
    return np.convolve(rate, h, mode="full")[:len(rate)]


def bin_trace(t: np.ndarray, y: np.ndarray, bin_s: float
              ) -> Tuple[np.ndarray, np.ndarray]:
    """Average a trace onto a coarser uniform grid.

    The transport kernel's dispersion sets the finest timescale that can be
    recovered at all, so the MS record is binned to the analysis grid before
    inversion.  This both respects that resolution limit and averages the
    measurement noise down by the square root of the binning factor.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    edges = np.arange(t[0], t[-1] + bin_s, bin_s)
    idx = np.clip(np.digitize(t, edges) - 1, 0, len(edges) - 2)
    nb = len(edges) - 1
    sums = np.bincount(idx, weights=y, minlength=nb)
    cnts = np.bincount(idx, minlength=nb)
    keep = cnts > 0
    centres = (edges[:-1] + bin_s / 2.0)[keep]
    return centres, sums[keep] / cnts[keep]
