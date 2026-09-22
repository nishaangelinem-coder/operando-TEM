"""Manuscript figures, rebuilt from the cached study outputs in results/raw.

Figures are print artefacts, so only the light-mode steps of the validated
categorical palette are used.  Every figure carries a legend, and every figure
has a companion table under results/tables -- which is also what the palette's
contrast warning on the aqua slot obliges (the relief rule).
"""
from __future__ import annotations

import json
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

import common as K
from opcem import (correlate, instruments, mobility, report, transport, truth,
                   vision)
from opcem.config import DEFAULT

P = report.PALETTE


def load(name: str):
    path = os.path.join(report.RAW, f"{name}.json")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return json.load(fh)["data"]


def frame(obj) -> pd.DataFrame:
    return pd.DataFrame(obj)


def _plt():
    report.style()
    import matplotlib.pyplot as plt
    return plt


# --------------------------------------------------------------------------- #
def fig1_platform() -> None:
    """Example micrograph, detection, population recovery, and the MS delay."""
    plt = _plt()
    cfg = K.short_cfg(300.0)
    exp = K.build(cfg, seed=1, duration=300.0, with_chip=False)
    rng = np.random.default_rng(7)
    k = 10                      # early, while the catalyst is still dispersed
    fov = 20.0                  # a 20 nm crop holds ~20 atoms at this density
    crop = (10.0, 10.0)
    img = instruments.render_adf_frame(cfg, exp.patch.pos3d[k],
                                       exp.patch.sizes[k], rng,
                                       crop_origin_nm=crop, field_nm=fov)
    txy, tsz = instruments.crop_truth(exp.patch.pos3d[k], exp.patch.sizes[k],
                                      crop, fov)
    det = vision.detect_atoms(cfg, img)

    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))

    a = ax[0, 0]
    ext = [0, fov, 0, fov]
    # At 200 e-/A^2 per frame the raw counts are shot-noise dominated (mean
    # ~0.4 counts per pixel), so the raw frame is visually uninformative even
    # though the atoms are detectable.  What is shown is the matched-filtered
    # image, which is what the detector operates on and what an analyst looks
    # at; the detections are the ones found on it.
    sig_pix = cfg.imaging.probe_sigma_nm / cfg.imaging.pixel_nm
    from scipy import ndimage as _nd
    disp = _nd.gaussian_filter(img - _nd.gaussian_filter(
        img, sigma=max(8.0 * sig_pix, 8.0)), sigma=sig_pix)
    a.imshow(disp, origin="lower", extent=ext, cmap="gray_r",
             vmin=np.percentile(disp, 2), vmax=np.percentile(disp, 99.9))
    if len(txy):
        a.scatter(txy[:, 0], txy[:, 1], s=70, facecolors="none",
                  edgecolors=P["truth"], linewidths=0.9,
                  label=f"ground truth (n={len(txy)})")
    if len(det.xy):
        a.scatter(det.xy[:, 0], det.xy[:, 1], s=16, marker="x",
                  color=P["observed"], linewidths=0.9,
                  label=f"detected (n={len(det.xy)})")
    a.set(xlabel="x (nm)", ylabel="y (nm)",
          title="(a) matched-filtered ADF-STEM crop, "
                f"{cfg.imaging.dose_rate:.0f} "
                r"e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$")
    a.legend(loc="upper right", fontsize=7, frameon=True, framealpha=0.88,
             facecolor=P["surface"], edgecolor="none")
    a.grid(False)
    # A rendered micrograph is the one panel a reader could mistake for data
    # if it were separated from its caption, so it carries the word on its face.
    for col, al, z in ((("#ffffff"), 0.40, 6), (P["observed"], 0.32, 7)):
        a.text(0.5, 0.5, "SIMULATED", transform=a.transAxes, fontsize=24,
               color=col, alpha=al, ha="center", va="center", rotation=30,
               zorder=z, fontweight="bold")

    b = ax[0, 1]
    tp = exp.patch.t
    b.plot(tp, exp.patch.pop[:, 0], color=P["truth"], label=r"true $N_1$")
    b.plot(exp.desc.t, exp.desc.N1, color=P["observed"], lw=1.0,
           label=r"observed $N_1$")
    b.plot(tp, exp.patch.pop[:, 1], color=P["truth"], ls="--",
           label=r"true $N_2$")
    b.plot(exp.desc.t, exp.desc.N2, color=P["observed"], ls="--", lw=1.0,
           label=r"observed $N_2$")
    b.set(xlabel="time (s)", ylabel="entities in field",
          title="(b) population recovery")
    b.legend(ncol=2, loc="upper right")

    c = ax[1, 0]
    t_ms = exp.t_ms
    ion = exp.ms_patch_cf.ion_current["CO2"]
    tr = np.interp(t_ms, exp.patch.t, exp.patch.true_rate)
    m = (t_ms > 40) & (t_ms < 130)
    c.plot(t_ms[m], tr[m], color=P["truth"], label="true rate")
    scale = np.nanmean(tr[m]) / np.nanmean(ion[m])
    c.plot(t_ms[m], ion[m] * scale, color=P["uncorrected"],
           label=r"MS ion current (scaled)")
    rec = K.recover_rate(exp, "counterfactual", "deconvolve")
    mm = (rec["t"] > 40) & (rec["t"] < 130)
    c.plot(rec["t"][mm], rec["rate"][mm], color=P["corrected"],
           label="deconvolved")
    td = cfg.transport.tau_dead()
    lo, hi = c.get_ylim()
    c.set_ylim(lo, hi + 0.22 * (hi - lo))       # headroom for the annotation
    lo, hi = c.get_ylim()
    y_ann = hi - 0.07 * (hi - lo)
    c.annotate("", xy=(104 + td, y_ann), xytext=(104, y_ann),
               arrowprops=dict(arrowstyle="<->", color=P["accent"], lw=1.1))
    c.text(104 + td / 2, y_ann - 0.075 * (hi - lo),
           rf"$\tau_\mathrm{{dead}}={td:.2f}$ s", color=P["accent"],
           ha="center", fontsize=7.5)
    c.set(xlabel="time (s)", ylabel=r"rate (molecules s$^{-1}$)",
          title="(c) transport delay and its removal")
    c.legend(loc="upper left", ncol=3, fontsize=6.8)

    d = ax[1, 1]
    doses = np.geomspace(1e2, 1e4, 40)
    d.plot(doses, [instruments.single_atom_cnr(cfg, x) for x in doses],
           color=P["truth"], label="single-atom CNR")
    d.axhline(cfg.imaging.detect_threshold_cnr, color=P["observed"], ls=":",
              label="detection threshold")
    d2 = d.twiny()
    d2.set_visible(False)
    d.set(xscale="log", xlabel=r"dose rate (e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$)",
          ylabel="contrast-to-noise ratio",
          title="(d) detection limit vs dose")
    d.legend(loc="upper left")

    report.save_figure(fig, "fig1_platform",
                       "Simulated platform output. (a) rendered ADF-STEM crop "
                       "with ground-truth and detected atom positions; (b) "
                       "recovered monomer and dimer counts against truth, "
                       "showing the size-dependent detection bias; (c) the "
                       "transport delay between the true rate and the MS "
                       "channel, and its removal by deconvolution; (d) "
                       "single-atom contrast-to-noise against dose rate.")
    plt.close(fig)


def fig2_consistency() -> None:
    """Patch kMC ensemble against the coarse-grained chip model."""
    plt = _plt()
    cfg = K.short_cfg(600.0)
    ens = truth.patch_ensemble(cfg, (1, 2, 3, 4, 5, 6, 7, 8), 600.0,
                               dose_rate=0.0, record_fps=1.0)
    mf = truth.chip_mean_field(cfg, ens["t"], truth.Schedule(cfg),
                               k_enc=K.k_enc(cfg), n_bins=1)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    a = ax[0]
    for i, (name, ls) in enumerate((("N_1", "-"), ("N_2", "--"),
                                    ("N_{4+}", ":"))):
        col = (0, 1, 3)[i]
        a.plot(ens["t"], ens["pop"][:, col], color=P["truth"], ls=ls,
               label=rf"kMC ${name}$")
        a.plot(mf["t"], mf["pop"][:, col], color=P["observed"], ls=ls, lw=1.0,
               label=rf"chip model ${name}$")
    a.set(xlabel="time (s)", ylabel="entities per field-equivalent area",
          title="(a) population trajectories")
    a.legend(ncol=2, fontsize=6.8, loc="upper right")

    b = ax[1]
    b.plot(ens["t"], ens["rate"], color=P["truth"], label="kMC ensemble mean")
    b.fill_between(ens["t"], ens["rate"] - ens["rate_sd"],
                   ens["rate"] + ens["rate_sd"], color=P["truth"], alpha=0.18,
                   lw=0, label="kMC $\\pm$1 s.d. across replicas")
    b.plot(mf["t"], mf["rate"], color=P["observed"], label="chip mean field")
    mass_err = abs(mf["mass"][-1] / mf["mass"][0] - 1.0)
    b.text(0.03, 0.06, f"Pt mass drift {mass_err:.1e}", transform=b.transAxes,
           fontsize=7, color=P["neutral"])
    b.set(xlabel="time (s)", ylabel=r"rate (molecules s$^{-1}$)",
          title="(b) rate: microscopic vs coarse-grained")
    b.legend(loc="upper right")
    report.save_figure(fig, "fig2_consistency",
                       "Consistency of the two levels of the digital twin. "
                       "The chip model has one fitted coefficient (the "
                       "encounter rate) and otherwise uses the kMC's own rate "
                       "constants; Pt mass is conserved to machine precision.")
    plt.close(fig)


def fig3_transport() -> None:
    plt = _plt()
    cfg = DEFAULT
    dt = 1.0 / cfg.ms.rate_hz
    t = np.arange(0.0, 150.0, dt)
    starts = [12.0, 37.0, 62.0, 87.0, 112.0]
    u = transport.make_pulse_train(t, starts, width_s=1.0)
    rng = np.random.default_rng(5)
    h_true = transport.kernel_from_config(cfg, dt)
    y = 1.7 * transport.forward(u, h_true) + 0.05
    y = y + rng.normal(0.0, cfg.ms.noise_rel * y.max(), len(t))
    fit = transport.calibrate_from_tracer(t, u, y, n_boot=200, seed=2)

    s1 = load("s1_transport")
    summ = frame(s1["summary"]) if s1 else None

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
    a = ax[0]
    m = t < 40
    a.plot(t[m], u[m] * y.max(), color=P["accent"], lw=1.0,
           label="tracer valve (scaled)")
    a.plot(t[m], y[m], color=P["observed"], label="MS response")
    yp = 1.7 * transport.forward(u, fit.kernel) + 0.05
    a.plot(t[m], yp[m], color=P["corrected"], ls="--", label="fitted model")
    a.set(xlabel="time (s)", ylabel="ion current (a.u.)",
          title="(a) tracer pulse and fit")
    a.legend(loc="upper right", fontsize=7)

    b = ax[1]
    kt = np.arange(len(h_true)) * dt
    b.plot(kt, h_true / h_true.max(), color=P["truth"], label="true kernel")
    b.plot(fit.kernel_t, fit.kernel / fit.kernel.max(), color=P["corrected"],
           ls="--", label="fitted kernel")
    b.set(xlim=(0, 12), xlabel="time (s)", ylabel="h(t), normalised",
          title="(b) impulse response")
    b.legend(loc="upper right")

    c = ax[2]
    if summ is not None:
        c.errorbar(summ.flow_sccm, summ.tau_fit_mean,
                   yerr=summ.tau_fit_sd.fillna(0), fmt="o", ms=5,
                   color=P["corrected"], capsize=3, label=r"fitted $\tau_{dead}$")
        c.plot(summ.flow_sccm, summ.tau_true, "s--", ms=5, color=P["truth"],
               mfc="none", label=r"true $\tau_{dead}$")
        c.errorbar(summ.flow_sccm, summ.sigma_fit_mean,
                   yerr=summ.sigma_fit_sd.fillna(0), fmt="^",
                   ms=5, color=P["observed"], capsize=3,
                   label=r"fitted $\sigma_{disp}$")
        c.plot(summ.flow_sccm, summ.sigma_true, "v--", ms=5,
               color=P["accent"], mfc="none", label=r"true $\sigma_{disp}$")
    c.set(xlabel="flow (sccm)", ylabel="time (s)",
          title="(c) recovery at three flows")
    c.legend(fontsize=6.8, loc="upper right")
    report.save_figure(fig, "fig3_transport",
                       "Calibration of the gas-transport transfer function "
                       "from timestamped tracer pulses, and recovery of its "
                       "dead time and dispersion at three flow rates.")
    plt.close(fig)


def fig4_events() -> None:
    plt = _plt()
    from s3_events import DILUTE, DURATION, PRE, POST, expected_delta_windowed
    cfg = K.sparse_cfg(DURATION, n_atoms=DILUTE, pre_s=PRE, post_s=POST)
    exp = K.build(cfg, seed=201, duration=DURATION, with_chip=False)
    rec = K.recover_rate(exp, "counterfactual", "deconvolve")
    t_r, r_hat = rec["t"], rec["rate"]
    ev = K.truth_event_times(exp.patch, "nucleate")
    pre, post = PRE, POST
    lags = np.arange(-pre, post + cfg.analysis.bin_s, cfg.analysis.bin_s)
    stack, null_stack = [], []
    rng = np.random.default_rng(3)
    for e in ev:
        seg = np.interp(e + lags, t_r, r_hat)
        base = np.mean(seg[lags < 0])
        stack.append(seg - base)
    for _ in range(300):
        e = rng.uniform(t_r[0] + pre, t_r[-1] - post)
        seg = np.interp(e + lags, t_r, r_hat)
        null_stack.append(seg - np.mean(seg[lags < 0]))
    stack = np.asarray(stack); null_stack = np.asarray(null_stack)
    m = stack.mean(0)
    se = stack.std(0, ddof=1) / np.sqrt(len(stack))
    nsd = null_stack.std(0, ddof=1)

    s3 = load("s3_events")
    power = frame(s3["power"]) if s3 else None
    sweep = frame(s3["sweep"]) if s3 else None
    _d0, dtrue, _tau = expected_delta_windowed(cfg)

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    a = ax[0]
    a.fill_between(lags, -2 * nsd, 2 * nsd, color=P["uncorrected"], alpha=0.25,
                   lw=0, label=r"random-time null, $\pm2$ s.d.")
    a.plot(lags, m, color=P["corrected"], label=f"event average (n={len(ev)})")
    a.fill_between(lags, m - 1.96 * se, m + 1.96 * se, color=P["corrected"],
                   alpha=0.3, lw=0, label="95% CI")
    a.axvline(0, color=P["neutral"], lw=0.8, ls=":")
    a.axhline(dtrue, color=P["truth"], lw=1.0, ls="--",
              label=rf"ground truth $\Delta r$ = {dtrue:.3f}")
    a.set(xlabel="time relative to nucleation event (s)",
          ylabel=r"$\Delta$ rate (molecules s$^{-1}$)",
          title="(a) event-related average")
    a.legend(fontsize=6.8, loc="upper left")

    b = ax[2]
    if power is not None:
        g = power.groupby("n_events").power.agg(["mean", "std"]).reset_index()
        b.errorbar(g.n_events, g["mean"], yerr=g["std"].fillna(0), fmt="o-",
                   ms=5, color=P["corrected"], capsize=3, label="empirical power")
        b.axhline(0.8, color=P["observed"], ls=":", label="80% power")
        n80 = g.loc[g["mean"] >= 0.8, "n_events"]
        if len(n80):
            b.axvline(n80.iloc[0], color=P["accent"], ls="--", lw=1.0,
                      label=f"n = {int(n80.iloc[0])} events")
    b.set(xscale="log", xlabel="number of independent events",
          ylabel="power", ylim=(0, 1.05), title="(c) detection power")
    b.legend(loc="lower right")

    c = ax[1]
    if sweep is not None:
        c.errorbar(sweep.event_density_lambda_W, sweep.ratio_dec_over_W,
                   yerr=(sweep.delta_dec_sd / sweep.delta_W).fillna(0),
                   fmt="o-", ms=5, color=P["corrected"], capsize=3,
                   label="deconvolved MS")
        c.plot(sweep.event_density_lambda_W, sweep.ratio_true_over_W, "s--",
               ms=5, mfc="none", color=P["truth"],
               label="noise-free true rate")
        c.axhline(1.0, color=P["neutral"], ls=":", lw=1.0,
                  label="window-averaged truth")
        for _, r in sweep.iterrows():
            c.annotate(f"{int(r.n_pt_atoms)}",
                       (r.event_density_lambda_W, r.ratio_dec_over_W),
                       textcoords="offset points", xytext=(4, 5),
                       fontsize=6.5, color=P["neutral"])
        c.set(xscale="log", xlabel=r"event density $\lambda(t_{pre}+t_{post})$",
              ylabel=r"recovered / $\Delta_W$",
              title="(b) event-density criterion")
        c.legend(fontsize=6.6, loc="lower left")
    report.save_figure(fig, "fig4_events",
                       "Event-aligned estimator in the counterfactual regime. "
                       "(a) event-related average of the deconvolved rate "
                       "around Pt1+Pt1 nucleation, against the window-averaged "
                       "ground-truth effect and a random-time null; (b) "
                       "recovery against the dimensionless event density, "
                       "labelled by Pt atoms per field; (c) empirical power "
                       "against the number of events.")
    plt.close(fig)


def fig5_correction() -> None:
    plt = _plt()
    s4 = load("s4_transport_effect")
    if s4 is None:
        return
    df = frame(s4["per_seed"])
    ll = frame(s4["lead_lag"])
    order = ["raw", "shift", "deconvolve"]
    cols = {"raw": P["uncorrected"], "shift": P["observed"],
            "deconvolve": P["corrected"]}
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))

    a = ax[0]
    for i, meth in enumerate(order):
        v = df.loc[df.correction == meth, "recovery_ratio"]
        a.bar(i, v.mean(), yerr=v.std(), color=cols[meth], width=0.62,
              capsize=4, label=meth)
        a.text(i, v.mean() + (v.std() if np.isfinite(v.std()) else 0) + 0.04,
               f"{v.mean():.2f}", ha="center", fontsize=7.5,
               color=P["neutral"])
    a.axhline(1.0, color=P["truth"], ls="--", lw=1.0, label="exact effect")
    a.set_xticks(range(len(order))); a.set_xticklabels(order)
    a.set(ylabel=r"recovered / true $\Delta r$",
          title="(a) event effect recovery")
    a.legend(fontsize=7, loc="upper left")

    b = ax[1]
    for i, meth in enumerate(order):
        v = df.loc[df.correction == meth, "rate_rmse"]
        sd = v.std() if np.isfinite(v.std()) else 0.0
        b.bar(i, v.mean(), yerr=sd, color=cols[meth], width=0.62, capsize=4)
        b.text(i, (v.mean() + sd) * 1.18, f"{v.mean():.4f}", ha="center",
               fontsize=7.5, color=P["neutral"])
    b.set_xticks(range(len(order))); b.set_xticklabels(order)
    b.set(yscale="log", ylabel=r"rate RMSE (molecules s$^{-1}$)",
          title="(b) rate reconstruction error")
    lo, hi = b.get_ylim()
    b.set_ylim(lo, hi * 2.2)

    c = ax[2]
    for i, meth in enumerate(order):
        v = ll.loc[ll.correction == meth, "peak_lag_s"]
        c.bar(i, v.mean(), yerr=v.std(), color=cols[meth], width=0.62,
              capsize=4)
    td = DEFAULT.transport.tau_dead()
    c.axhline(td, color=P["accent"], ls="--", lw=1.0,
              label=rf"$\tau_\mathrm{{dead}}$ = {td:.2f} s")
    c.axhline(0.0, color=P["truth"], ls=":", lw=1.0, label="true lag = 0")
    frac = ll.loc[ll.correction == "raw", "significant"].mean()
    c.text(0.02, 0.62, f"called significant\nin {frac:.0%} of runs",
           transform=c.transAxes, fontsize=7, color=P["neutral"])
    c.set_xticks(range(len(order))); c.set_xticklabels(order)
    c.set(ylabel="lead-lag peak (s)", title="(c) apparent lead-lag")
    lo, hi = c.get_ylim()
    c.set_ylim(lo, hi + 0.30 * (hi - lo))
    c.legend(fontsize=7, loc="upper right")
    report.save_figure(fig, "fig5_correction",
                       "What skipping the transport correction costs. Without "
                       "it the event-aligned effect is attenuated and the "
                       "lead-lag peak sits near the transport dead time, which "
                       "would be read as the structure leading the chemistry.")
    plt.close(fig)


def fig6_beam() -> None:
    plt = _plt()
    s6 = load("s6_beam")
    if s6 is None:
        return
    g = frame(s6["grouped"]); fits = frame(s6["fits"]); bpi = frame(s6["bpi"])
    gate = s6.get("usable_gate", 0.8)
    k_eff_0 = s6.get("k_eff_zero_dose", np.nan)
    lin = fits.iloc[0]; site = fits.iloc[1]; fd = fits.iloc[2]
    # JSON round-trips booleans as 0/1, and a numeric Series passed to
    # DataFrame.__getitem__ is read as a column selector, not a mask.
    usable = g.usable.astype(bool)
    use, bad = g[usable], g[~usable]
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))

    a = ax[0]
    a.errorbar(use.dose_rate, use.k_hop_obs, yerr=use.k_hop_sem.fillna(0),
               fmt="o", ms=5, color=P["observed"], capsize=3,
               label="observed (usable)")
    a.scatter(bad.dose_rate, bad.k_hop_obs, marker="x", s=28,
              color=P["uncorrected"], label="excluded by the gate")
    a.plot(g.dose_rate, g.k_eff_true, "s--", ms=4, mfc="none",
           color=P["truth"], label="true effective hop rate")
    xx = np.linspace(0, g.dose_rate.max() * 1.05, 60)
    a.plot(xx, lin.zero_dose_value + lin.beta_fit * xx, color=P["uncorrected"],
           ls="-", lw=1.1, label="linear fit")
    a.scatter([0], [lin.zero_dose_value], marker="*", s=70,
              color=P["uncorrected"], zorder=5)
    a.scatter([0], [site.zero_dose_value], marker="*", s=70,
              color=P["corrected"], zorder=5, label="site-aware intercept")
    a.axhline(k_eff_0, color=P["accent"], ls=":", lw=1.0,
              label=rf"truth at zero dose = {k_eff_0:.2f}")
    a.set(xscale="log", xlabel=r"dose rate (e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$)",
          ylabel=r"hop rate (s$^{-1}$)", title="(a) dose series")
    a.set_ylim(0, max(g.k_hop_obs.max(), g.k_eff_true.max()) * 1.55)
    a.legend(fontsize=5.8, loc="upper left", ncol=2, columnspacing=0.8,
             handlelength=1.4)

    b = ax[1]
    b.plot(bpi.dose_rate, bpi.bpi_linear_fit, "o-", ms=5,
           color=P["uncorrected"], label="linear model")
    b.fill_between(bpi.dose_rate, bpi.bpi_linear_ci_lo, bpi.bpi_linear_ci_hi,
                   color=P["uncorrected"], alpha=0.22, lw=0)
    b.plot(bpi.dose_rate, bpi.bpi_site_aware_fit, "^-", ms=5,
           color=P["corrected"], label="site-aware model")
    b.plot(bpi.dose_rate, bpi.bpi_analytic_truth, "s--", ms=5, mfc="none",
           color=P["truth"], label="analytic truth")
    b.axhline(DEFAULT.analysis.bpi_exclusion_threshold, color=P["observed"],
              ls=":", label="pre-registered threshold")
    b.set(xscale="log",
          xlabel=r"dose rate (e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$)",
          ylabel="beam perturbation index", title="(b) BPI vs truth")
    # a log axis spanning one decade puts minor ticks on top of each other
    b.set_xticks([2e3, 5e3, 1e4, 2e4])
    b.set_xticklabels(["2k", "5k", "10k", "20k"])
    b.set_xticks([], minor=True)
    b.legend(fontsize=6.2, loc="upper left")

    c = ax[2]
    c.errorbar(use.dose_rate, use.frac_dispersed,
               yerr=use.frac_dispersed_sem.fillna(0), fmt="o", ms=5,
               color=P["observed"], capsize=3, label="observed (usable)")
    c.scatter(bad.dose_rate, bad.frac_dispersed, marker="x", s=28,
              color=P["uncorrected"], label="excluded by the gate")
    c.plot(xx, fd.zero_dose_value + fd.beta_fit * xx, color=P["corrected"],
           label="zero-dose extrapolation")
    c.scatter([0], [fd.zero_dose_value], marker="*", s=90, color=P["accent"],
              zorder=5, label=f"beam-free estimate {fd.zero_dose_value:.2f}")
    c.set(xlabel=r"dose rate (e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$)",
          ylabel="dispersed fraction", title="(c) the reported observable")
    c.set_xticks([0, 5000, 10000, 15000, 20000])
    c.set_xticklabels(["0", "5k", "10k", "15k", "20k"])
    c.legend(fontsize=6.2, loc="lower left")
    report.save_figure(fig, "fig6_beam",
                       "Separating chemistry from irradiation. (a) the dose "
                       f"series, with doses below the {gate:.2f} "
                       "persistent-track gate excluded because false positives "
                       "dominate the tracking, and the two zero-dose "
                       "intercepts; (b) the Beam Perturbation Index under both "
                       "dose models against the analytic truth; (c) the "
                       "reported dispersed fraction is itself dose-dependent "
                       "and needs extrapolating.")
    plt.close(fig)


def fig7_representativeness() -> None:
    plt = _plt()
    s7a = load("s7a_budget")
    pvc = load("s7b_patch_vs_chip")
    psync = load("s7e_perturbation_sync")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7))

    a = ax[0]
    if s7a:
        tr = frame(s7a["tradeoff"])
        # Detectability and time resolution BOTH scale as the reciprocal of
        # the flow, so no flow improves one without costing the other by the
        # same factor.  Plotting one against the other shows the trade-off
        # directly and avoids a dual-axis chart, in which the two curves would
        # simply lie on top of each other.
        a.plot(tr.chip_rate_over_limit, tr.resolvable_lag_s, "o-", ms=5,
               color=P["truth"], label="operating curve")
        for _, r in tr.iterrows():
            a.annotate(f"{r.flow_sccm:g}", (r.chip_rate_over_limit,
                                            r.resolvable_lag_s),
                       textcoords="offset points", xytext=(5, -9),
                       fontsize=6.4, color=P["neutral"])
        a.axvline(3.0, color=P["observed"], ls=":",
                  label=r"detectability floor (3$\times$)")
        ok = tr[tr.chip_detectable.astype(bool)]
        if len(ok):
            a.axhspan(ok.resolvable_lag_s.min(), a.get_ylim()[1] * 4,
                      color=P["corrected"], alpha=0.10, lw=0)
            a.text(0.05, 0.92,
                   f"feasible only above\n{ok.resolvable_lag_s.min():.0f} s "
                   "time resolution",
                   transform=a.transAxes, fontsize=6.8, color=P["neutral"],
                   va="top")
        a.set(xscale="log", yscale="log",
              xlabel="chip signal / detection limit",
              ylabel="smallest resolvable lag (s)",
              title="(a) the flow trade-off (labels: sccm)")
        a.legend(fontsize=6.6, loc="lower right")

    b = ax[1]
    if s7a:
        bud = s7a["budget"]
        names = ["imaged\nfield", "whole\nchip", "MS limit\n(2 sccm)"]
        vals = [bud["field_rate_molec_per_s"], bud["chip_rate_molec_per_s"],
                bud["qms_flux_limit_molec_per_s"]]
        colours = [P["observed"], P["truth"], P["uncorrected"]]
        b.bar(names, vals, color=colours, width=0.6)
        b.set_yscale("log")
        b.set_ylim(min(vals) / 30, max(vals) * 1e5)
        for i, v in enumerate(vals):
            b.text(i, v * 3.0, f"{v:.1e}", ha="center", fontsize=7,
                   color=P["neutral"])
        b.set(ylabel=r"CO$_2$ flux (molecules s$^{-1}$)",
              title="(b) sensitivity budget")
        b.text(0.03, 0.95,
               f"field short by {bud['field_detector_shortfall']:.0e}"
               r"$\times$" + f"\n$f_{{rep}}$ = {bud['f_rep']:.1e}",
               transform=b.transAxes, fontsize=7, color=P["accent"],
               va="top")

    c = ax[2]
    if pvc:
        d = frame(pvc["summary"])
        xx = d.dose_rate.replace(0, 30)
        c.plot(xx, d.patch_dispersed_frac, "o-", ms=5, color=P["observed"],
               label="illuminated field")
        c.fill_between(xx,
                       d.patch_dispersed_frac - d.patch_sd.fillna(0),
                       d.patch_dispersed_frac + d.patch_sd.fillna(0),
                       color=P["observed"], alpha=0.22, lw=0)
        c.axhline(d.chip_dispersed_frac.iloc[0], color=P["truth"], ls="--",
                  label="un-illuminated chip")
        c.set(xscale="log",
              xlabel=r"dose rate (e$^-$ $\mathrm{\AA}^{-2}$ s$^{-1}$)",
              ylabel="dispersed fraction",
              title="(c) field vs chip state")
        c.legend(fontsize=7, loc="lower left")
    report.save_figure(fig, "fig7_representativeness",
                       "The representativeness problem. (a) the flow that "
                       "makes the product detectable is the flow that destroys "
                       "the time resolution; (b) the imaged field's product "
                       "flux against the whole chip's and the detection "
                       "limit; (c) the illuminated field's structural state "
                       "diverges systematically from the un-illuminated chip.")
    plt.close(fig)


def fig8_control() -> None:
    plt = _plt()
    s8 = load("s8_control")
    cfg = K.short_cfg(900.0)
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.7))

    a = ax[0]
    if s8:
        d = frame(s8["summary"])
        et = d[d.policy == "event_triggered"].sort_values("threshold")
        ol = d[d.policy == "open_loop"].iloc[0]
        fs = d[d.policy == "fixed_schedule"].iloc[0]
        a.errorbar(et.threshold, et.integrated_yield,
                   yerr=et.integrated_yield_sd.fillna(0), fmt="o-", ms=5,
                   color=P["corrected"], capsize=3, label="event-triggered")
        a.axhline(fs.integrated_yield, color=P["observed"], ls="--",
                  label="fixed schedule (same budget)")
        a.axhline(ol.integrated_yield, color=P["uncorrected"], ls=":",
                  label="open loop (no actuation)")
        a.set(xlabel=r"trigger threshold on $N_{\geq 3}$",
              ylabel="integrated product (molecules)",
              title="(a) yield vs trigger threshold")
        a.legend(fontsize=7, loc="lower right")

    b = ax[1]
    from opcem import control
    from opcem.truth import PatchSimulator, Schedule
    best_th = 2.0
    if s8:
        d = frame(s8["comparison"])
        best_th = float(d.best_threshold.iloc[0])
    for name, pol_f, col, ls in (
            ("open loop", lambda: control.OpenLoop(cfg), P["uncorrected"], ":"),
            ("fixed schedule",
             lambda: control.FixedSchedule(cfg, duration_s=900.0),
             P["observed"], "--"),
            (f"event-triggered (th={best_th:g})",
             lambda: control.EventTriggered(cfg, threshold=best_th),
             P["corrected"], "-")):
        pol = pol_f()
        res = PatchSimulator(cfg, seed=701, dose_rate=cfg.imaging.dose_rate,
                             schedule=Schedule(cfg)).run(duration=900.0,
                                                         controller=pol)
        w = max(1, int(round(10.0 / (res.t[1] - res.t[0]))))
        sm = np.convolve(res.true_rate, np.ones(w) / w, mode="same")
        b.plot(res.t, sm, color=col, ls=ls, label=name)
    b.set(xlabel="time (s)", ylabel=r"rate (molecules s$^{-1}$, 10 s mean)",
          title="(b) example trajectories")
    b.legend(fontsize=7, loc="upper right")
    report.save_figure(fig, "fig8_control",
                       "Closed-loop control against matched baselines. All "
                       "policies spend the same actuation budget, so the "
                       "comparison isolates when it is spent; the advantage "
                       "depends on the trigger threshold, which has to be "
                       "pre-registered.")
    plt.close(fig)


ALL = {
    "fig1": fig1_platform, "fig2": fig2_consistency, "fig3": fig3_transport,
    "fig4": fig4_events, "fig5": fig5_correction, "fig6": fig6_beam,
    "fig7": fig7_representativeness, "fig8": fig8_control,
}


def make_all(only: Optional[List[str]] = None) -> None:
    for name, fn in ALL.items():
        if only and name not in only:
            continue
        print(f"  {name} ...", end="", flush=True)
        fn()
        print(" ok")


if __name__ == "__main__":
    import sys
    make_all(sys.argv[1:] or None)
