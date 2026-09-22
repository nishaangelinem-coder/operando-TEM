"""Unit tests for the load-bearing pieces of the pipeline.

These test properties that would silently corrupt a result if broken:
mass conservation in the chip model, kernel moments, recovery of a known
transfer function, the equivalence of two independently coded deconvolution
solvers, the provenance guard, and the lead-lag resolution gate.

Run with:  python -m pytest tests -q     (or  python tests/test_opcem.py)
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from opcem import (beam, correlate, descriptors, instruments, mobility,
                   transport, truth, vision)
from opcem.config import DEFAULT, KB_EV


def _cfg(duration=120.0):
    import copy
    c = copy.deepcopy(DEFAULT)
    c.duration_s = duration
    return c


# --------------------------------------------------------------------------- #
def test_kernel_moments_match_closed_form():
    dt = 0.05
    h = transport.kernel_from_config(DEFAULT, dt, length_s=60.0)
    t = np.arange(len(h)) * dt
    m, sd = transport.moments(t, h)
    assert abs(m - DEFAULT.transport.tau_dead()) < 0.01, (m,)
    assert abs(sd - DEFAULT.transport.sigma_disp()) < 0.01, (sd,)
    assert abs(h.sum() - 1.0) < 1e-12


def test_tracer_calibration_recovers_kernel():
    dt = 0.1
    t = np.arange(0.0, 150.0, dt)
    u = transport.make_pulse_train(t, [12.0, 50.0, 90.0, 125.0], width_s=1.0)
    h = transport.kernel_from_config(DEFAULT, dt)
    rng = np.random.default_rng(0)
    y = 1.5 * transport.forward(u, h) + 0.04
    y = y + rng.normal(0.0, 0.02 * y.max(), len(t))
    fit = transport.calibrate_from_tracer(t, u, y, n_boot=60, seed=1)
    assert abs(fit.tau_dead_s - DEFAULT.transport.tau_dead()) < 0.05
    assert abs(fit.sigma_disp_s - DEFAULT.transport.sigma_disp()) < 0.05


def test_deconvolution_solvers_agree():
    """The FFT and dense solvers must agree at a fixed regularisation weight."""
    dt = 1.0
    t = np.arange(0.0, 400.0, dt)
    h = transport.kernel_from_config(DEFAULT, dt)
    x = 1.5 + 0.4 * np.sin(2 * np.pi * t / 60.0)
    y = transport.forward(x, h)
    a = transport.tikhonov_deconvolve(y, h, lam=1.0, method="dense")["x"]
    b = transport.tikhonov_deconvolve(y, h, lam=1.0, method="fft")["x"]
    inner = slice(30, -30)          # ignore edge transients
    rel = np.max(np.abs(a[inner] - b[inner])) / np.ptp(x)
    assert rel < 0.05, rel


def test_deconvolution_beats_no_correction():
    dt = 1.0
    t = np.arange(0.0, 600.0, dt)
    h = transport.kernel_from_config(DEFAULT, dt)
    rng = np.random.default_rng(2)
    x = 1.5 + 0.5 * np.sin(2 * np.pi * t / 50.0)
    y = transport.forward(x, h) + rng.normal(0.0, 0.02, len(t))
    dec = transport.tikhonov_deconvolve(y, h)["x"]
    sh = transport.shift_correct(y, dt, DEFAULT.transport.tau_dead())
    inner = slice(20, -20)
    e_raw = np.sqrt(np.mean((y[inner] - x[inner]) ** 2))
    e_sh = np.sqrt(np.mean((sh[inner] - x[inner]) ** 2))
    e_de = np.sqrt(np.mean((dec[inner] - x[inner]) ** 2))
    assert e_de < e_sh < e_raw, (e_de, e_sh, e_raw)


# --------------------------------------------------------------------------- #
def test_chip_model_conserves_pt_mass():
    c = _cfg()
    t = np.arange(0.0, 400.0, 1.0)
    mf = truth.chip_mean_field(c, t, truth.Schedule(c), k_enc=1.3e-4,
                               n_bins=1)
    drift = abs(mf["mass"][-1] / mf["mass"][0] - 1.0)
    assert drift < 1e-6, drift


def test_patch_conserves_pt_atoms():
    c = _cfg()
    res = truth.PatchSimulator(c, seed=4, dose_rate=1e3).run(duration=200.0)
    total = np.array([np.sum(s) for s in res.sizes])
    assert np.all(total == c.n_pt_atoms), (total.min(), total.max())


def test_tof_hierarchy_dimer_is_best_per_atom():
    c = _cfg()
    t_k = c.reactor.t_setpoint_c + 273.15
    p_co = c.reactor.x_co * c.reactor.p_total_mbar / 1000.0
    p_o2 = c.reactor.x_o2 * c.reactor.p_total_mbar / 1000.0
    per_atom = {n: truth.tof_entity(c, n, "terrace", t_k, p_co, p_o2) / n
                for n in (1, 2, 3, 4, 6)}
    assert max(per_atom, key=per_atom.get) == 2, per_atom


def test_beam_contribution_is_linear_in_dose():
    c = _cfg()
    t_k = c.reactor.t_setpoint_c + 273.15
    k0 = sum(truth.k_hop(c, "terrace", t_k, 0.0))
    k1 = sum(truth.k_hop(c, "terrace", t_k, 1e4))
    assert abs((k1 - k0) - c.kinetics.beta_beam * 1e4) < 1e-9


# --------------------------------------------------------------------------- #
def test_dose_series_fit_recovers_intercept_and_slope():
    phi = np.array([1e2, 3e2, 1e3, 3e3, 1e4])
    k_chem, bet = 2.47, 2.0e-4
    rng = np.random.default_rng(1)
    k = k_chem + bet * phi + rng.normal(0.0, 0.02, len(phi))
    fit = beam.fit_dose_series(phi, k, np.full(len(phi), 0.02))
    assert abs(fit.k_chem - k_chem) < 0.1, fit.k_chem
    assert abs(fit.beta - bet) / bet < 0.15, fit.beta
    assert abs(fit.bpi(1e3) - bet * 1e3 / (k_chem + bet * 1e3)) < 0.03


def test_lead_lag_gate_rejects_unresolvable_lags():
    c = _cfg()
    dt = 1.0
    t = np.arange(0.0, 400.0, dt)
    x = np.sin(2 * np.pi * t / 40.0)
    y = np.roll(x, 1)                      # 1 s lag, below 2*sigma_disp
    ll = correlate.lead_lag(c, x, y, dt)
    assert not ll.significant
    assert "unresolvable" in ll.reason
    y2 = np.roll(x, 6)                     # 6 s lag, resolvable
    ll2 = correlate.lead_lag(c, x, y2, dt)
    assert ll2.significant


def test_event_aligned_detects_injected_step():
    """Irregularly spaced events, so the circular-shift null cannot realign."""
    c = _cfg()
    dt = 1.0
    t = np.arange(0.0, 1800.0, dt)
    rng = np.random.default_rng(3)
    y = np.full_like(t, 1.0) + rng.normal(0.0, 0.02, len(t))
    events = np.sort(rng.uniform(30.0, 1760.0, 45))
    for e in events:                        # inject a +0.1 step for 6 s
        y[(t > e) & (t <= e + 6)] += 0.1
    eff = correlate.event_aligned_effect(c, t, y, events, n_boot=300,
                                         n_surrogate=300, seed=1)
    assert 0.07 < eff.delta < 0.13, eff.delta
    assert eff.ci_lo > 0
    assert eff.p_surrogate < 0.05, eff.p_surrogate


def test_circular_shift_null_is_conservative_for_periodic_events():
    """A periodic event train against a periodic signal is NOT significant.

    Shifting the trace by one event period realigns it exactly, so the
    surrogate null contains values as large as the observed effect.  That is
    the correct answer: with a strictly periodic driver, temporal coincidence
    carries almost no information, and an i.i.d. permutation null would
    wrongly call it significant.
    """
    c = _cfg()
    dt = 1.0
    t = np.arange(0.0, 900.0, dt)
    rng = np.random.default_rng(3)
    y = np.full_like(t, 1.0) + rng.normal(0.0, 0.02, len(t))
    events = np.arange(40.0, 860.0, 30.0)
    for e in events:
        y[(t > e) & (t <= e + 6)] += 0.1
    eff = correlate.event_aligned_effect(c, t, y, events, n_boot=200,
                                         n_surrogate=300, seed=1)
    assert eff.ci_lo > 0                      # the effect itself is real
    assert eff.p_surrogate > 0.05, eff.p_surrogate   # but not distinguishable


def test_event_aligned_null_on_random_events():
    c = _cfg()
    dt = 1.0
    t = np.arange(0.0, 900.0, dt)
    rng = np.random.default_rng(5)
    y = 1.0 + rng.normal(0.0, 0.05, len(t))
    ev = np.sort(rng.uniform(20.0, 880.0, 40))
    eff = correlate.event_aligned_effect(c, t, y, ev, n_boot=300,
                                         n_surrogate=200, seed=2)
    assert eff.ci_lo < 0 < eff.ci_hi, (eff.ci_lo, eff.ci_hi)


def test_detrending_removes_a_linear_drift():
    c = _cfg()
    dt = 1.0
    t = np.arange(0.0, 900.0, dt)
    y = 2.0 - 0.001 * t                     # pure downward drift, no events
    ev = np.arange(40.0, 860.0, 30.0)
    off = correlate.event_aligned_effect(c, t, y, ev, detrend=False,
                                         n_boot=100, n_surrogate=0, seed=1)
    on = correlate.event_aligned_effect(c, t, y, ev, detrend=True,
                                        n_boot=100, n_surrogate=0, seed=1)
    assert abs(on.delta) < abs(off.delta) / 10.0, (on.delta, off.delta)


# --------------------------------------------------------------------------- #
def test_provenance_guard_blocks_snapshot3d_columns():
    import pandas as pd
    df = pd.DataFrame({"N1": [1.0], "cn_pt_support": [2.0]})
    df.attrs["provenance"] = {"N1": "census2d", "cn_pt_support": "snapshot3d"}
    descriptors.check_provenance(df, ["N1"], ["census2d"])
    try:
        descriptors.check_provenance(df, ["N1", "cn_pt_support"],
                                     ["census2d"])
    except ValueError as exc:
        assert "cn_pt_support" in str(exc)
    else:
        raise AssertionError("provenance guard did not fire")


def test_msd_noise_correction_removes_static_noise_floor():
    c = _cfg()
    sig = instruments.localisation_sigma_nm(c, c.imaging.dose_rate)
    # a perfectly immobile atom observed with localisation noise
    msd_pure_noise = 2.0 * sig ** 2
    k = descriptors.k_hop_from_msd(c, msd_pure_noise, 0.2,
                                   c.imaging.dose_rate)
    assert k == 0.0, k


def test_eta_analytic_matches_numerical_surface_average():
    c = _cfg()
    rng = np.random.default_rng(0)
    n = 40000
    x = rng.uniform(0, 40, n); y = rng.uniform(0, 40, n)
    q = 2 * np.pi / c.support.corrug_wavelength_nm
    a = c.support.corrug_amp_nm
    gx = a * q * np.cos(q * x) * np.sin(q * y)
    gy = a * q * np.sin(q * x) * np.cos(q * y)
    eta_num = 0.5 * np.mean(gx ** 2 + gy ** 2)
    assert abs(eta_num - mobility.eta_analytic(c)) / eta_num < 0.03


def test_projection_loses_out_of_plane_displacement():
    c = _cfg()
    res = truth.PatchSimulator(c, seed=9, dose_rate=1e3).run(duration=200.0)
    gt = mobility.msd_3d_truth(res.pos3d, res.ids, res.sizes, lag_frames=1,
                               box_nm=c.support.field_nm)
    assert gt["n"] > 100
    assert gt["msd_2d"] < gt["msd_3d"]
    assert abs(gt["eta"] - mobility.eta_analytic(c)) < 0.12, gt["eta"]


def test_detector_recall_rises_with_dose():
    c = _cfg()
    res = truth.PatchSimulator(c, seed=8, dose_rate=1e3).run(duration=60.0)
    recalls = []
    for dose in (3e2, 1e3, 1e4):
        fr = instruments.observe_positions(c, res,
                                           np.random.default_rng(1),
                                           dose_rate=dose)
        n_true = sum(len(s) for s in res.sizes)
        n_real = sum(int((f.true_id >= 0).sum()) for f in fr)
        recalls.append(n_real / n_true)
    assert recalls[0] < recalls[1] < recalls[2] <= 1.0, recalls


def test_two_stage_least_squares_beats_ols_under_confounding():
    rng = np.random.default_rng(4)
    n = 800
    z = rng.integers(0, 2, n).astype(float)     # randomised instrument
    u = rng.normal(size=n)                      # the confounder
    x = 0.8 * z + 1.2 * u + rng.normal(0, 0.3, n)
    beta = 0.5
    y = beta * x - 1.5 * u + rng.normal(0, 0.3, n)
    out = correlate.two_stage_least_squares(z, x, y)
    assert abs(out["beta_2sls"] - beta) < abs(out["beta_ols"] - beta)
    assert abs(out["beta_2sls"] - beta) < 0.15, out


def test_control_policies_share_the_actuation_budget():
    from opcem import control
    c = _cfg(600.0)
    for pol in (control.FixedSchedule(c, duration_s=600.0),
                control.EventTriggered(c, threshold=1.0)):
        sch = truth.Schedule(c)
        truth.PatchSimulator(c, seed=6, dose_rate=1e3, schedule=sch).run(
            duration=600.0, controller=pol)
        assert pol.log.n_actuations <= c.control.actuation_budget, pol.name


def test_hmm_recovers_two_injected_regimes():
    import pandas as pd
    rng = np.random.default_rng(1)
    n = 600
    state = (np.arange(n) // 100) % 2
    x = np.where(state == 0, 1.0, 4.0) + rng.normal(0, 0.3, n)
    df = pd.DataFrame({"a": x, "b": x * 0.5 + rng.normal(0, 0.3, n)})
    out = correlate.hmm_state_flux(df, ["a", "b"], n_states=2, seed=0)
    assert out["ok"]
    agree = max(np.mean(out["states"] == state),
                np.mean(out["states"] != state))
    assert agree > 0.9, agree


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    bad = 0
    for name, fn in fns:
        try:
            fn()
            print(f"PASS  {name}")
        except Exception as exc:
            bad += 1
            print(f"FAIL  {name}: {type(exc).__name__}: {exc}")
    print(f"\n{len(fns) - bad}/{len(fns)} passed")
    raise SystemExit(1 if bad else 0)
