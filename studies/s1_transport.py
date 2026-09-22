"""S1 -- can the transport transfer function be recovered from tracer pulses?

This is step 4 of the measurement protocol and the precondition for every
event-level statement (gap G1).  The test is run at three flow rates, because
the delay scales with the reciprocal of the flow and the calibration has to be
repeated at every set point rather than measured once.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import report, transport
from opcem.config import DEFAULT

FLOWS = (1.0, 2.0, 4.0)
N_REPEAT = 5


def run(cfg=DEFAULT, seed: int = 31) -> pd.DataFrame:
    dt = 1.0 / cfg.ms.rate_hz
    t = np.arange(0.0, 150.0, dt)
    starts = [12.0, 37.0, 62.0, 87.0, 112.0]
    u = transport.make_pulse_train(t, starts, width_s=1.0)
    rows = []
    for flow in FLOWS:
        tr_cfg = cfg.transport.scaled(flow)
        h_true = transport.gamma_kernel(np.arange(0.0, 40.0, dt),
                                        tr_cfg.tau0_s, tr_cfg.gamma_k,
                                        tr_cfg.gamma_theta_s)
        for rep in range(N_REPEAT):
            rng = np.random.default_rng(seed * 101 + int(flow * 10) + rep)
            y = 1.7 * transport.forward(u, h_true) + 0.05
            y = y * (1.0 + cfg.ms.drift_rel * 0.5
                     * np.sin(2 * np.pi * t / 90.0))
            y = y + rng.normal(0.0, cfg.ms.noise_rel * y.max(), len(t))
            fit = transport.calibrate_from_tracer(t, u, y, n_boot=200,
                                                  seed=int(rng.integers(1e6)))
            rows.append({
                "flow_sccm": flow, "replicate": rep,
                "tau_dead_true_s": tr_cfg.tau_dead(),
                "tau_dead_fit_s": fit.tau_dead_s,
                "tau_dead_ci_lo": fit.tau_dead_ci[0],
                "tau_dead_ci_hi": fit.tau_dead_ci[1],
                "sigma_disp_true_s": tr_cfg.sigma_disp(),
                "sigma_disp_fit_s": fit.sigma_disp_s,
                "sigma_disp_ci_lo": fit.sigma_disp_ci[0],
                "sigma_disp_ci_hi": fit.sigma_disp_ci[1],
                "fit_rmse": fit.rmse, "n_pulses": fit.n_pulses,
                "tau_err_s": fit.tau_dead_s - tr_cfg.tau_dead(),
                "sigma_err_s": fit.sigma_disp_s - tr_cfg.sigma_disp(),
                "tau_covered": bool(fit.tau_dead_ci[0] <= tr_cfg.tau_dead()
                                    <= fit.tau_dead_ci[1]),
                "sigma_covered": bool(fit.sigma_disp_ci[0]
                                      <= tr_cfg.sigma_disp()
                                      <= fit.sigma_disp_ci[1]),
            })
    df = pd.DataFrame(rows)
    summ = df.groupby("flow_sccm").agg(
        tau_true=("tau_dead_true_s", "first"),
        tau_fit_mean=("tau_dead_fit_s", "mean"),
        tau_fit_sd=("tau_dead_fit_s", "std"),
        tau_bias=("tau_err_s", "mean"),
        tau_coverage=("tau_covered", "mean"),
        sigma_true=("sigma_disp_true_s", "first"),
        sigma_fit_mean=("sigma_disp_fit_s", "mean"),
        sigma_fit_sd=("sigma_disp_fit_s", "std"),
        sigma_bias=("sigma_err_s", "mean"),
        sigma_coverage=("sigma_covered", "mean"),
        rmse=("fit_rmse", "mean"), n=("replicate", "count")).reset_index()
    report.save_table(summ, "table_S1_transport_calibration",
                      "Recovery of the transport transfer function from "
                      "timestamped tracer pulses at three flow rates "
                      "(5 replicates each). 'coverage' is the fraction of "
                      "bootstrap 95% intervals containing the true value.",
                      cfg=cfg, seed=seed)
    report.save_json({"per_replicate": df, "summary": summ}, "s1_transport",
                     cfg=cfg, seed=seed)
    return summ


if __name__ == "__main__":
    print(run().to_string(index=False))
