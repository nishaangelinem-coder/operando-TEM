"""S4 -- the cost of not correcting for transport, and S4b -- association
versus intervention.

S4 asks the question an experimentalist actually faces: if the transfer
function is not measured, what does the analysis conclude?  Three correction
levels (none, delay-only, full deconvolution) are compared on identical data,
for both the event-aligned effect and the lead-lag direction.

S4b constructs the confounding that gap G8 warns about.  Temperature is
modulated, which moves both the structure (through the hop and dissociation
rates) and the rate (through the Arrhenius factor) -- a textbook common
driver.  The observational regression of rate on structure is then compared
with a two-stage least-squares estimate that uses a *randomised* actuator as
the instrument.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import correlate, report, transport
from opcem.config import DEFAULT
from opcem.truth import Actuation, Schedule
from s3_events import DILUTE, PRE, POST, expected_delta_windowed

SEEDS = tuple(range(301, 321))   # 20 seeds: S4c was underpowered at 4
DURATION = 1200.0


def run_s4(cfg=DEFAULT, seed: int = 61):
    # The dilute loading of study S3 is used, because that is the regime in
    # which single events are temporally resolvable at all; comparing
    # correction levels where no event can be isolated would compare noise.
    cfg = K.sparse_cfg(DURATION, n_atoms=DILUTE, pre_s=PRE, post_s=POST)
    _d0, dtrue, _tau = expected_delta_windowed(cfg)
    rows, ll_rows = [], []
    for s in SEEDS:
        exp = K.build(cfg, seed=s, duration=DURATION, with_chip=False)
        ev = K.truth_event_times(exp.patch, "nucleate")
        for method in ("raw", "shift", "deconvolve"):
            rec = K.recover_rate(exp, "counterfactual", method)
            t_r, r_hat = rec["t"], rec["rate"]
            e = correlate.event_aligned_effect(cfg, t_r, r_hat, ev,
                                               pre_s=PRE, post_s=POST, seed=s)
            truth_r = np.interp(t_r, exp.patch.t, exp.patch.true_rate)
            rows.append({"seed": s, "correction": method,
                         "delta": e.delta, "ci_lo": e.ci_lo, "ci_hi": e.ci_hi,
                         "expected_delta": dtrue,
                         "recovery_ratio": e.delta / dtrue,
                         "rate_rmse": float(np.sqrt(np.nanmean(
                             (r_hat - truth_r) ** 2))),
                         "n_events": e.n_events,
                         "p_surrogate": e.p_surrogate})
            # lead-lag between the dimer count and the recovered rate
            x = np.interp(t_r, exp.desc.t.to_numpy(),
                          exp.desc.N2.to_numpy())
            dt = float(np.mean(np.diff(t_r)))
            ll = correlate.lead_lag(cfg, x, r_hat, dt, max_lag_s=20.0)
            ll_rows.append({"seed": s, "correction": method,
                            "peak_lag_s": ll.peak_lag,
                            "peak_rho": ll.peak_rho,
                            "significant": ll.significant,
                            "gate_reason": ll.reason})
    df = pd.DataFrame(rows)
    summ = df.groupby("correction").agg(
        rate_rmse=("rate_rmse", "mean"),
        delta=("delta", "mean"), delta_sd=("delta", "std"),
        recovery_ratio=("recovery_ratio", "mean"),
        recovery_ratio_sd=("recovery_ratio", "std"),
        expected_delta=("expected_delta", "first"),
        p_surrogate=("p_surrogate", "mean"),
        n_events=("n_events", "mean")).reset_index()
    ll = pd.DataFrame(ll_rows)
    ll_s = ll.groupby("correction").agg(
        peak_lag_s=("peak_lag_s", "mean"), peak_lag_sd=("peak_lag_s", "std"),
        peak_rho=("peak_rho", "mean"),
        frac_called_significant=("significant", "mean")).reset_index()

    report.save_table(summ, "table_S4a_transport_correction",
                      "Effect of the transport correction level on the "
                      "event-aligned estimate. COUNTERFACTUAL regime. "
                      "'recovery_ratio' is the estimated effect divided by "
                      "the window-averaged ground-truth effect.",
                      cfg=cfg, seed=seed)
    report.save_table(ll_s, "table_S4b_lead_lag",
                      "Lead-lag peak between the observed dimer count and the "
                      "recovered rate. Without correction the peak sits at a "
                      "lag close to the transport dead time and would be read "
                      "as the structure leading the chemistry by that amount.",
                      cfg=cfg, seed=seed)
    report.save_json({"per_seed": df, "lead_lag": ll}, "s4_transport_effect",
                     cfg=cfg, seed=seed)
    return {"correction": summ, "lead_lag": ll_s}


def _instrument_arm(cfg, s, strength: str):
    """Build one run of the interventional arm at a given actuator strength.

    A weak instrument gives an unbiased but useless estimate, so the design
    requirement is that the randomised actuator must move the structural
    variable strongly.  Both arms are run so that the requirement is
    demonstrated rather than asserted.
    """
    rng = np.random.default_rng(s if strength == "weak" else s + 5000)
    acts = []
    n_blocks = int(DURATION // 50)
    for i in range(n_blocks):
        acts.append(Actuation(t_start=i * 50.0, duration=50.0, kind="t_step",
                              delta_t_c=float(25.0 * np.sin(i / 2.2))))
    if strength == "weak":
        dur, x_o2 = 8.0, 0.20
    else:
        # A long, strong oxidising pulse: long compared with the dimer
        # lifetime and with the cluster dissociation time it has to beat.
        dur, x_o2 = 30.0, 0.60
    z_times = []
    for i in range(n_blocks):
        if rng.random() < 0.5:
            t0 = i * 50.0 + 10.0
            acts.append(Actuation(t_start=t0, duration=dur, kind="o2_pulse",
                                  x_o2=x_o2))
            z_times.append(t0)
    sch = Schedule(cfg, acts)
    exp = K.build(cfg, seed=s, duration=DURATION, schedule=sch,
                  with_chip=False)
    rec = K.recover_rate(exp, "counterfactual", "deconvolve")
    t_r, r_hat = rec["t"], rec["rate"]
    d = exp.desc
    x = np.interp(t_r, d.t.to_numpy(), d.N2.to_numpy())
    z = np.zeros_like(t_r)
    for t0 in z_times:
        z[(t_r >= t0) & (t_r < t0 + dur + 20.0)] = 1.0
    return exp, t_r, r_hat, d, x, z, len(z_times)


def run_s4b(cfg=DEFAULT, seed: int = 62):
    """Observational vs interventional estimate under a known confounder."""
    # Technical loading: this estimator does not need single events to be
    # temporally resolvable, only the population to respond to the instrument.
    cfg = K.short_cfg(DURATION)
    rows = []
    for s in SEEDS:
        for strength in ("weak", "strong"):
            _exp, t_r, r_hat, d, x, z, n_on = _instrument_arm(cfg, s, strength)
            iv = correlate.two_stage_least_squares(z, x, r_hat)
            fit = correlate.fit_elastic_net(
                d, np.interp(d.t.to_numpy(), t_r, r_hat), ["N2"],
                force_cols=["T_c", "p_co", "p_o2"], cfg=cfg, seed=s)
            rows.append({"seed": s, "instrument": strength,
                         "beta_ols_naive": iv["beta_ols"],
                         "beta_2sls": iv["beta_2sls"],
                         "se_2sls": iv["se_2sls"],
                         "first_stage_r2": iv["first_stage_r2"],
                         "beta_partial_elasticnet":
                             fit["coef"].get("N2", np.nan),
                         "n": iv["n"], "n_instrument_on": int(n_on)})
    df = pd.DataFrame(rows)
    summ = df.groupby("instrument", sort=False).agg(
        beta_ols_naive=("beta_ols_naive", "mean"),
        beta_ols_naive_sd=("beta_ols_naive", "std"),
        beta_partial_elasticnet=("beta_partial_elasticnet", "mean"),
        beta_2sls=("beta_2sls", "mean"), beta_2sls_sd=("beta_2sls", "std"),
        first_stage_r2=("first_stage_r2", "mean"),
        n_seeds=("seed", "count")).reset_index()
    report.save_table(summ, "table_S4c_observational_vs_interventional",
                      "Effect of the dimer count on the recovered rate, "
                      "estimated three ways under a deliberately confounded "
                      "temperature modulation: naive OLS, elastic net with "
                      "reactor variables partialled out, and 2SLS "
                      "instrumented by a randomised O2 pulse, at two actuator "
                      "strengths. A weak instrument gives an unbiased but "
                      "imprecise estimate, so the design requirement is a "
                      "strong one.",
                      cfg=cfg, seed=seed)
    report.save_json({"per_seed": df, "summary": summ},
                     "s4b_interventional", cfg=cfg, seed=seed)
    return summ
