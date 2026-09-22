"""S3 -- the event-aligned estimator: when does it work, with what attenuation,
and with what power?

This study runs in the *counterfactual* regime in which the whole active area
of the reactor is the imaged field (representativeness fraction = 1) and the
product detector has matching sensitivity.  That regime is physically
unavailable -- study S7 quantifies the shortfall by about nine orders of
magnitude -- but it is the regime in which the estimator's own properties can
be measured against known truth, which is the point of a digital twin.  Every
table produced here says so.

Two things are measured.

**Where the estimator works.** The primary independent variable is Pt loading,
because the nucleation rate scales roughly as the square of the monomer count
and therefore sets the dimensionless event density ``lambda*(t_pre+t_post)``.
When that approaches one, the windows of neighbouring events overlap and no
single event's contribution can be isolated.  The boundary converts directly
into a maximum Pt loading for the experiment.

**Where the attenuation comes from.** The recovered effect is compared against
three references, which separates the estimator's unavoidable window averaging
from any loss caused by the instrument:

* ``delta_0``  -- the exact instantaneous single-event effect, TOF(2)-2*TOF(1);
* ``delta_W``  -- the same effect averaged over a post-window of length W given
  the dimer's finite lifetime tau, i.e. ``delta_0 * (tau/W) * (1-exp(-W/tau))``;
* ``delta_true`` -- the effect the same estimator recovers from the *noise-free
  true rate*, which is the achievable ceiling for any instrument.

The ratio ``delta_dec / delta_true`` is then the cost attributable to the
mass spectrometer and the transport inversion alone.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import correlate, report
from opcem.config import DEFAULT
from opcem.truth import k_diss, site_occupancy, tof_entity

LOADINGS = (6, 8, 16, 32, 90)
SEEDS = (201, 202, 203, 204)
DURATION = 1200.0
PRE = POST = 4.0
DILUTE = 8          # loading used for the controls and the power curve


def expected_delta(cfg) -> float:
    """Exact instantaneous effect of one nucleation event on the field rate."""
    t_k = cfg.reactor.t_setpoint_c + 273.15
    p_co = cfg.reactor.x_co * cfg.reactor.p_total_mbar / 1000.0
    p_o2 = cfg.reactor.x_o2 * cfg.reactor.p_total_mbar / 1000.0
    _occ, _k, mult = site_occupancy(cfg, t_k, cfg.imaging.dose_rate)
    tof1 = tof_entity(cfg, 1, "terrace", t_k, p_co, p_o2) * mult
    tof2 = tof_entity(cfg, 2, "terrace", t_k, p_co, p_o2) * mult
    return float(tof2 - 2.0 * tof1)


def expected_delta_windowed(cfg, post_s: float = POST) -> tuple:
    """Window-averaged expected effect, given the dimer's finite lifetime."""
    d0 = expected_delta(cfg)
    t_k = cfg.reactor.t_setpoint_c + 273.15
    p_o2 = cfg.reactor.x_o2 * cfg.reactor.p_total_mbar / 1000.0
    tau = 1.0 / k_diss(cfg, 2, t_k, p_o2)
    factor = (tau / post_s) * (1.0 - np.exp(-post_s / tau))
    return d0, float(d0 * factor), float(tau)


def _true_rate_on_grid(exp, bin_s: float = 1.0):
    t = np.arange(float(exp.patch.t[0]), float(exp.patch.t[-1]), bin_s)
    return t, np.interp(t, exp.patch.t, exp.patch.true_rate)


def run(cfg=DEFAULT, seed: int = 51):
    sweep_rows, ctrl_rows, power_rows = [], [], []
    for n_at in LOADINGS:
        c = K.sparse_cfg(DURATION, n_atoms=n_at, pre_s=PRE, post_s=POST)
        d0, dW, tau = expected_delta_windowed(c)
        for s in SEEDS:
            exp = K.build(c, seed=s, duration=DURATION, with_chip=False)
            lam = K.event_rate(exp.patch, "nucleate")
            ev = K.truth_event_times(exp.patch, "nucleate")
            tt, rt = _true_rate_on_grid(exp)
            e_true = correlate.event_aligned_effect(
                c, tt, rt, ev, pre_s=PRE, post_s=POST, seed=s)
            rec = K.recover_rate(exp, "counterfactual", "deconvolve")
            e_dec = correlate.event_aligned_effect(
                c, rec["t"], rec["rate"], ev, pre_s=PRE, post_s=POST, seed=s)
            base = float(np.nanmean(rec["rate"]))
            sweep_rows.append({
                "n_pt_atoms": n_at,
                "pt_density_per_nm2": n_at / c.support.field_nm ** 2,
                "seed": s, "lambda_event_per_s": lam,
                "event_density_lambda_W": lam * (PRE + POST),
                "mean_rate": base, "delta_0": d0, "delta_W": dW,
                "tau_dimer_s": tau,
                "delta_0_rel": d0 / base,
                "delta_true": e_true.delta,
                "delta_dec": e_dec.delta,
                "ci_lo": e_dec.ci_lo, "ci_hi": e_dec.ci_hi,
                "n_events": e_dec.n_events, "n_events_total": len(ev),
                "p_surrogate": e_dec.p_surrogate, "z_vs_null": e_dec.z,
                "ratio_true_over_W": e_true.delta / dW,
                "ratio_dec_over_true": (e_dec.delta / e_true.delta
                                        if e_true.delta else np.nan),
                "ratio_dec_over_W": e_dec.delta / dW,
                "detected": bool(e_dec.ci_lo > 0
                                 and e_dec.p_surrogate < c.analysis.alpha)})
            if n_at == DILUTE:
                for label, kind, obs, detr in (
                        ("truth_events_detrended", "nucleate", False, True),
                        ("truth_events_no_detrend", "nucleate", False, False),
                        ("observed_events_detrended", "nucleate", True, True),
                        ("dissociation_control", "dissoc", False, True)):
                    evx = (K.observed_event_times(exp.desc, kind) if obs
                           else K.truth_event_times(exp.patch, kind))
                    e = correlate.event_aligned_effect(
                        c, rec["t"], rec["rate"], evx, pre_s=PRE, post_s=POST,
                        detrend=detr, seed=s)
                    row = {"seed": s, "event_set": label, "delta_W": dW,
                           "n_events_total": len(evx)}
                    row.update(e.as_row())
                    row["ratio_over_W"] = e.delta / dW
                    ctrl_rows.append(row)
                pw = correlate.detection_power(
                    c, rec["t"], rec["rate"], ev,
                    n_grid=(5, 10, 20, 40, 80), n_rep=30, seed=s)
                pw["seed"] = s
                power_rows.append(pw)

    df = pd.DataFrame(sweep_rows)
    sweep = df.groupby("n_pt_atoms").agg(
        pt_density_per_nm2=("pt_density_per_nm2", "first"),
        lambda_event_per_s=("lambda_event_per_s", "mean"),
        event_density_lambda_W=("event_density_lambda_W", "mean"),
        n_events_used=("n_events", "mean"),
        mean_rate=("mean_rate", "mean"),
        delta_0_rel=("delta_0_rel", "mean"),
        delta_W=("delta_W", "first"),
        delta_true=("delta_true", "mean"),
        delta_dec=("delta_dec", "mean"), delta_dec_sd=("delta_dec", "std"),
        ratio_true_over_W=("ratio_true_over_W", "mean"),
        ratio_dec_over_true=("ratio_dec_over_true", "mean"),
        ratio_dec_over_W=("ratio_dec_over_W", "mean"),
        p_surrogate=("p_surrogate", "mean"),
        frac_detected=("detected", "mean")).reset_index()

    ctrl = pd.DataFrame(ctrl_rows)
    ctrl_s = ctrl.groupby("event_set").agg(
        n_events=("n_events", "mean"), delta=("delta", "mean"),
        delta_sd=("delta", "std"), delta_rel=("delta_rel", "mean"),
        ci_lo=("ci_lo", "mean"), ci_hi=("ci_hi", "mean"),
        p_surrogate=("p_surrogate", "mean"), z_vs_null=("z", "mean"),
        delta_W=("delta_W", "first"),
        ratio_over_W=("ratio_over_W", "mean")).reset_index()

    power = pd.concat(power_rows, ignore_index=True)
    power_s = power.groupby("n_events").agg(
        power=("power", "mean"), power_sd=("power", "std"),
        mean_delta=("mean_delta", "mean"),
        sd_delta=("sd_delta", "mean")).reset_index()

    report.save_table(sweep, "table_S3a_event_density_sweep",
                      "Recovery of the single-event effect against Pt loading, "
                      "which sets the dimensionless event density "
                      "lambda*(t_pre+t_post). COUNTERFACTUAL regime. "
                      "'ratio_dec_over_true' isolates the cost of the mass "
                      "spectrometer and the transport inversion from the "
                      "estimator's own window averaging.", cfg=cfg, seed=seed)
    report.save_table(ctrl_s, "table_S3b_event_aligned_controls",
                      f"Event-aligned effect at {DILUTE} Pt atoms per field "
                      "with its controls. The dissociation row is a sign "
                      "control; 'no_detrend' shows what progressive sintering "
                      "does to an undetrended estimate; 'observed_events' uses "
                      "the events the pipeline actually detects.",
                      cfg=cfg, seed=seed)
    report.save_table(power_s, "table_S3c_detection_power",
                      f"Empirical power of the event-aligned test against the "
                      f"number of independent events at {DILUTE} Pt atoms per "
                      "field (bootstrap 95% CI excluding zero; 4 runs x 30 "
                      "subsamples). The apparent power at very small n is a "
                      "bootstrap artefact and is excluded from interpretation.",
                      cfg=cfg, seed=seed)
    report.save_json({"sweep_per_seed": df, "sweep": sweep, "controls": ctrl,
                      "power": power}, "s3_events", cfg=cfg, seed=seed)
    return {"sweep": sweep, "controls": ctrl_s, "power": power_s}
