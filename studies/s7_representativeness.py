"""S7 -- is the imaged field representative, and can spontaneous single-atom
events be seen in the product signal at all? (gap G5)

This is the study that decides what the platform can honestly claim.  Three
questions:

1.  What fraction of the reactor's active area does one field of view cover,
    and what product flux do the imaged atoms actually generate?  Compared
    with the detection limit of a quadrupole MS, this sets whether a
    spontaneous event in the field can ever appear in the signal.
2.  Does the imaged field's structural state track the chip's?  The chip is
    not illuminated, so the beam alone makes the two differ.
3.  If spontaneous-event analysis is unavailable, does
    *perturbation-synchronised* event analysis work?  A global actuator step
    drives the same structural transition in every patch of the chip at once,
    so the chip-averaged rate does respond -- and the imaged field becomes a
    sampler of that synchronised response rather than the source of the signal.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import correlate, report, truth
from opcem.config import DEFAULT
from opcem.truth import Actuation, Schedule

SEEDS = (601, 602, 603, 604)
DURATION = 1800.0
NA = 6.02214076e23
POSITIONS = ("inlet", "centre", "outlet")
FLOWS_SCCM = (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0)


def run_sensitivity(cfg=DEFAULT, seed: int = 91):
    """The sensitivity and representativeness budget, and the flow trade-off.

    A quadrupole mass spectrometer detects a *mole fraction*, so its flux
    detection limit is ``x_min * F`` and therefore falls with the total flow.
    The transport dead time, in contrast, *rises* as the reciprocal of the
    flow.  The two constraints pull in opposite directions, and the feasible
    window between them is the real design space of the experiment.
    """
    c = K.short_cfg(DURATION)
    exp = K.build(c, seed=SEEDS[0], duration=600.0, with_chip=False)
    rate_field = float(np.mean(exp.patch.true_rate))       # molecules/s
    f_rep = c.reactor.f_rep
    rate_chip = rate_field / f_rep
    from s3_events import expected_delta_windowed
    d_field = expected_delta_windowed(c)[0]

    def flux_limit(flow_sccm: float) -> float:
        return c.ms.detection_mole_fraction * flow_sccm * 7.4362e-7 * NA

    lim_ref = flux_limit(c.reactor.flow_sccm)
    budget = {
        "imaged_field_nm": c.support.field_nm,
        "reactive_zone_um2": (c.reactor.zone_length_um
                              * c.reactor.zone_width_um),
        "projected_area_fraction": c.reactor.projected_area_fraction,
        "catalyst_area_multiplier": c.reactor.catalyst_area_multiplier,
        "f_rep": f_rep,
        "pt_atoms_in_field": c.n_pt_atoms,
        "pt_atoms_on_chip": c.n_pt_atoms / f_rep,
        "field_rate_molec_per_s": rate_field,
        "chip_rate_molec_per_s": rate_chip,
        "flow_sccm": c.reactor.flow_sccm,
        "qms_flux_limit_molec_per_s": lim_ref,
        "chip_rate_over_limit": rate_chip / lim_ref,
        "field_rate_over_limit": rate_field / lim_ref,
        "one_event_rel_effect_on_field": d_field / rate_field,
        "one_event_rel_effect_on_chip": d_field / rate_chip,
        "ms_relative_noise_per_bin": c.ms.noise_rel,
        "events_needed_for_snr1_on_chip":
            (c.ms.noise_rel / (d_field / rate_chip)) ** 2,
        "field_detector_shortfall": lim_ref / rate_field,
    }
    df = pd.DataFrame([{"quantity": k, "value": v}
                       for k, v in budget.items()])

    # Sensitivity band on the single most load-bearing assumption: the ratio
    # of the catalyst film's Pt-accessible area to the projected zone area.
    # The conclusion has to survive a wide range of it, so the range is
    # reported rather than a point value.
    band = []
    for mult in (3.5e2, 3.5e3, 3.5e4, 3.5e5, 3.5e6):
        f = c.reactor.projected_area_fraction / mult
        band.append({"catalyst_area_multiplier": mult,
                     "f_rep": f,
                     "pt_atoms_on_chip": c.n_pt_atoms / f,
                     "chip_rate_molec_per_s": rate_field / f,
                     "chip_rate_over_limit": (rate_field / f) / lim_ref,
                     "field_detector_shortfall": lim_ref / rate_field})
    band = pd.DataFrame(band)

    trade = []
    for fl in FLOWS_SCCM:
        tr = c.transport.scaled(fl)
        lim = flux_limit(fl)
        trade.append({
            "flow_sccm": fl,
            "qms_flux_limit_molec_per_s": lim,
            "chip_rate_over_limit": rate_chip / lim,
            "chip_detectable": bool(rate_chip / lim >= 3.0),
            "tau_dead_s": tr.tau_dead(),
            "sigma_disp_s": tr.sigma_disp(),
            "resolvable_lag_s": (cfg.analysis.lag_significance_sigma_mult
                                 * tr.sigma_disp()),
            "field_rate_over_limit": rate_field / lim})
    trade = pd.DataFrame(trade)

    report.save_table(df, "table_S7a_representativeness_budget",
                      "Representativeness and sensitivity budget. The imaged "
                      "field's product flux is ~12 orders of magnitude below "
                      "any product detector's limit, so spontaneous "
                      "single-atom event analysis against a mass spectrometer "
                      "is excluded, not merely difficult. The whole chip is "
                      "measurable at low flow.", cfg=c, seed=seed)
    report.save_table(band, "table_S7b_frep_sensitivity",
                      "Sensitivity of the representativeness fraction and the "
                      "chip signal to the catalyst-area multiplier, which is "
                      "an estimate rather than a measurement. The field's "
                      "detector shortfall does not depend on it at all, so the "
                      "central conclusion is unchanged over four orders of "
                      "magnitude of the assumption.", cfg=c, seed=seed)
    report.save_table(trade, "table_S7c_flow_tradeoff",
                      "The flow trade-off. A quadrupole's flux detection "
                      "limit falls with flow, so low flow is needed to see the "
                      "product at all; the transport dead time and dispersion "
                      "rise as the reciprocal of the flow, so low flow "
                      "destroys the time resolution. 'resolvable_lag_s' is the "
                      "smallest lag the pipeline will call significant.",
                      cfg=c, seed=seed)
    report.save_json({"budget": budget, "tradeoff": trade, "band": band},
                     "s7a_budget", cfg=c, seed=seed)
    return {"budget": df, "frep_band": band, "tradeoff": trade}


def run_patch_vs_chip(cfg=DEFAULT, seed: int = 92) -> pd.DataFrame:
    """Does the illuminated field's structural state track the dark chip's?"""
    c = K.short_cfg(DURATION)
    rows = []
    t_grid = np.arange(0.0, DURATION + 1.0, 1.0)
    sch = Schedule(c)
    chip = truth.chip_mean_field(c, t_grid, sch, k_enc=K.k_enc(c),
                                 dose_rate=0.0)
    chip_frac = chip["pop"][:, 0] / np.maximum(chip["pop"].sum(axis=1), 1e-9)
    for dose in (0.0, 1e3, 1e4):
        for s in SEEDS:
            p = truth.PatchSimulator(c, seed=s, dose_rate=dose,
                                     schedule=Schedule(c)).run(
                duration=DURATION, record_fps=1.0)
            pf = p.pop[:, 0] / np.maximum(p.pop.sum(axis=1), 1e-9)
            pf_i = np.interp(t_grid, p.t, pf)
            rows.append({
                "dose_rate": dose, "seed": s,
                "patch_dispersed_frac": float(np.mean(pf_i)),
                "chip_dispersed_frac": float(np.mean(chip_frac)),
                "difference": float(np.mean(pf_i) - np.mean(chip_frac)),
                "rel_difference": float((np.mean(pf_i) - np.mean(chip_frac))
                                        / np.mean(chip_frac)),
                "pearson_r": float(np.corrcoef(pf_i, chip_frac)[0, 1])})
    df = pd.DataFrame(rows)
    summ = df.groupby("dose_rate").agg(
        patch_dispersed_frac=("patch_dispersed_frac", "mean"),
        patch_sd=("patch_dispersed_frac", "std"),
        chip_dispersed_frac=("chip_dispersed_frac", "first"),
        rel_difference=("rel_difference", "mean"),
        rel_difference_sd=("rel_difference", "std"),
        pearson_r=("pearson_r", "mean")).reset_index()
    report.save_table(summ, "table_S7d_patch_vs_chip",
                      "Dispersed fraction of the illuminated imaged field "
                      "versus the un-illuminated chip mean field. The chip "
                      "never sees the beam, so any non-zero dose makes the two "
                      "differ systematically, independently of detection bias.",
                      cfg=c, seed=seed)
    report.save_json({"per_seed": df, "summary": summ}, "s7b_patch_vs_chip",
                     cfg=c, seed=seed)
    return summ


def run_position_stratified(cfg=DEFAULT, seed: int = 93) -> pd.DataFrame:
    """Recovered effect at three positions along the flow axis."""
    from s3_events import DILUTE, PRE, POST, expected_delta_windowed
    c = K.sparse_cfg(DURATION, n_atoms=DILUTE, pre_s=PRE, post_s=POST)
    rows = []
    for pos_i, pos in enumerate(POSITIONS):
        cc = K.sparse_cfg(DURATION, n_atoms=DILUTE, pre_s=PRE, post_s=POST)
        frac = (pos_i + 0.5) / len(POSITIONS)
        # local temperature and product enrichment of that window
        cc.reactor.t_setpoint_c = (c.reactor.t_setpoint_c
                                   - c.reactor.t_gradient_c / 2.0
                                   + c.reactor.t_gradient_c * frac)
        enrich = 1.0 + c.reactor.product_enrichment * (frac - 0.5) * 2.0
        d_exp = expected_delta_windowed(cc)[1] * enrich
        for s in SEEDS:
            exp = K.build(cc, seed=s, duration=DURATION, with_chip=False)
            rec = K.recover_rate(exp, "counterfactual", "deconvolve")
            ev = K.truth_event_times(exp.patch, "nucleate")
            e = correlate.event_aligned_effect(cc, rec["t"], rec["rate"], ev,
                                               pre_s=PRE, post_s=POST, seed=s)
            rows.append({"position": pos, "local_T_c": cc.reactor.t_setpoint_c,
                         "enrichment": enrich, "seed": s,
                         "delta": e.delta, "ci_lo": e.ci_lo, "ci_hi": e.ci_hi,
                         "expected_delta": d_exp,
                         "recovery_ratio": e.delta / d_exp,
                         "n_events": e.n_events,
                         "mean_rate": e.baseline})
    df = pd.DataFrame(rows)
    summ = df.groupby("position", sort=False).agg(
        local_T_c=("local_T_c", "first"), enrichment=("enrichment", "first"),
        mean_rate=("mean_rate", "mean"),
        delta=("delta", "mean"), delta_sd=("delta", "std"),
        expected_delta=("expected_delta", "first"),
        recovery_ratio=("recovery_ratio", "mean"),
        n_events=("n_events", "mean")).reset_index()
    # heterogeneity test across positions (one-way ANOVA on the per-seed deltas)
    from scipy import stats
    groups = [df.loc[df.position == p, "delta"].to_numpy() for p in POSITIONS]
    f, pval = stats.f_oneway(*groups)
    summ.attrs["anova_F"] = float(f)
    summ.attrs["anova_p"] = float(pval)
    het = pd.DataFrame([{"anova_F": float(f), "anova_p": float(pval),
                         "n_per_group": len(groups[0]),
                         "delta_spread_frac":
                             float((summ.delta.max() - summ.delta.min())
                                   / summ.delta.mean())}])
    report.save_table(summ, "table_S7e_position_stratified",
                      "Event-aligned effect recovered at three pre-registered "
                      "positions along the flow axis, each with its own local "
                      "temperature and product enrichment.", cfg=c, seed=seed)
    report.save_table(het, "table_S7f_position_heterogeneity",
                      "Heterogeneity of the recovered effect across the three "
                      "positions (one-way ANOVA on per-run estimates).",
                      cfg=c, seed=seed)
    report.save_json({"per_seed": df, "summary": summ, "heterogeneity": het},
                     "s7c_positions", cfg=c, seed=seed)
    return summ


def run_perturbation_synchronised(cfg=DEFAULT, seed: int = 94) -> pd.DataFrame:
    """The route that survives a realistic representativeness fraction.

    A global O2 step drives oxidative redispersion everywhere on the chip at
    once, so the chip-averaged rate does move even though the imaged field
    contributes ~1e-12 of it.  Two things then have to be separated, and two
    design choices decide whether the analysis has any inferential power.

    **Direct versus structure-mediated.**  Raising the oxygen partial pressure
    changes the adsorbate coverages, which changes the rate immediately and has
    nothing to do with the catalyst's structure.  Measuring *during* the pulse
    therefore measures mostly that direct effect.  The structure-mediated
    component is isolated by blanking a gap of the pulse duration plus several
    times the transport dispersion and measuring the residual afterwards, when
    the gas composition is back at baseline: any rate difference that survives
    has to be carried by the catalyst's state.

    **Periodic versus randomised actuation.**  Under a circular-shift null, a
    strictly periodic pulse train realigns with itself after one period, so the
    null contains values as large as the observed effect and nothing can be
    concluded.  This is the correct answer -- with a purely periodic driver,
    temporal coincidence carries almost no information -- and it is why the
    randomised arm is the one that identifies anything.  Both arms are run.
    """
    c = K.short_cfg(DURATION)
    gap = 8.0 + 3.0 * c.transport.sigma_disp()
    rows = []
    for s in SEEDS:
        rng = np.random.default_rng(s * 3 + 1)
        periodic = list(np.arange(60.0, DURATION - 60.0, 40.0))
        # randomised arm: same number of pulses, exponential-ish spacing
        gaps = rng.uniform(20.0, 60.0, len(periodic))
        randomised, tnow = [], 60.0
        for g in gaps:
            if tnow > DURATION - 60.0:
                break
            randomised.append(float(tnow))
            tnow += g
        for arm, times in (("periodic", periodic),
                           ("randomised", randomised)):
            acts = [Actuation(t_start=t, duration=8.0, kind="o2_pulse",
                              x_o2=0.20) for t in times]
            sch = Schedule(c, acts)
            exp = K.build(c, seed=s, duration=DURATION, schedule=sch,
                          with_chip=True)
            for which, label in (("chip", "chip_rate"),
                                 ("counterfactual", "field_rate_cf")):
                rec = K.recover_rate(exp, which, "deconvolve")
                for mode, gp in (("during_pulse_direct", 0.0),
                                 ("after_pulse_structural", gap)):
                    e = correlate.event_aligned_effect(
                        c, rec["t"], rec["rate"], times, pre_s=20.0,
                        post_s=20.0, gap_s=gp, seed=s)
                    row = {"seed": s, "arm": arm, "signal": label,
                           "mode": mode, "n_pulses": len(times),
                           "gap_s": gp}
                    row.update(e.as_row())
                    rows.append(row)
            # structural response measured in the imaged field
            d = exp.desc
            for mode, gp in (("during_pulse_direct", 0.0),
                             ("after_pulse_structural", gap)):
                pre_v, post_v = [], []
                for t0 in times:
                    m0 = (d.t >= t0 - gp - 20.0) & (d.t < t0 - gp)
                    m1 = (d.t > t0 + gp) & (d.t <= t0 + gp + 20.0)
                    if m0.sum() > 2 and m1.sum() > 2:
                        pre_v.append(d.loc[m0, "frac_dispersed"].mean())
                        post_v.append(d.loc[m1, "frac_dispersed"].mean())
                if pre_v:
                    dv = np.asarray(post_v) - np.asarray(pre_v)
                    rows.append({
                        "seed": s, "arm": arm,
                        "signal": "field_dispersed_frac", "mode": mode,
                        "n_pulses": len(times), "gap_s": gp,
                        "n_events": len(dv), "delta": float(dv.mean()),
                        "baseline": float(np.mean(pre_v)),
                        "delta_rel": float(dv.mean() / np.mean(pre_v))})
    df = pd.DataFrame(rows)
    summ = df.groupby(["arm", "signal", "mode"], sort=False).agg(
        n_events=("n_events", "mean"), delta=("delta", "mean"),
        delta_sd=("delta", "std"), delta_rel=("delta_rel", "mean"),
        p_surrogate=("p_surrogate", "mean"),
        z_vs_null=("z", "mean")).reset_index()
    report.save_table(summ, "table_S7g_perturbation_synchronised",
                      "Perturbation-synchronised event analysis. A global O2 "
                      "step moves the chip-averaged rate even at a "
                      "representativeness fraction of ~5e-12. "
                      "'during_pulse_direct' measures the coverage-mediated "
                      "response; 'after_pulse_structural' blanks a gap of the "
                      "pulse duration plus three dispersions and measures what "
                      "survives, which can only be carried by the catalyst's "
                      "state. The periodic arm is not distinguishable from its "
                      "circular-shift null by construction; the randomised arm "
                      "is the one that identifies anything.", cfg=c, seed=seed)
    report.save_json({"per_seed": df, "summary": summ, "gap_s": gap},
                     "s7e_perturbation_sync", cfg=c, seed=seed)
    return summ


def run(cfg=DEFAULT, seed: int = 91):
    return {"sensitivity": run_sensitivity(cfg, seed),
            "patch_vs_chip": run_patch_vs_chip(cfg, seed + 1),
            "positions": run_position_stratified(cfg, seed + 2),
            "perturbation": run_perturbation_synchronised(cfg, seed + 3)}
