"""S2 -- descriptor fidelity, and whether dynamic descriptors earn their place.

Two questions:

1.  How accurately does the pipeline recover each element of D(t) from the
    observed streams?  Detection is size-dependent, so the *bias* matters more
    than the noise: a monomer recall below one inflates the apparent dimer
    fraction and therefore biases exactly the quantity the platform is meant
    to measure.
2.  Does adding rate-type descriptors (hop rate, nucleation rate, dimer
    lifetime) explain variance in the true rate that the census descriptors
    (counts, mean size, coordination) do not?  That is the falsifiable content
    of the central hypothesis (gap G6).

Both the median and the Viterbi size decoders are reported, because the choice
changes the dynamic descriptors by a factor of a few.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import correlate, descriptors, report
from opcem.config import DEFAULT, KB_EV
from opcem.truth import k_diss, k_hop

SEEDS = (101, 102, 103, 104, 105)
DURATION = 600.0


def truth_frame(exp) -> pd.DataFrame:
    """Ground-truth population and rate on the descriptor grid."""
    p = exp.patch
    t = exp.desc.t.to_numpy()
    out = {"t": t}
    for i, name in enumerate(("N1", "N2", "N3", "N4p")):
        out[f"{name}_true"] = np.interp(t, p.t, p.pop[:, i])
    out["rate_true"] = np.interp(t, p.t, p.true_rate)
    return pd.DataFrame(out)


def run(cfg=DEFAULT, seed: int = 41):
    cfg = K.short_cfg(DURATION)
    t_k = cfg.reactor.t_setpoint_c + 273.15
    p_o2 = cfg.reactor.x_o2 * cfg.reactor.p_total_mbar / 1000.0
    khop_true = sum(k_hop(cfg, "terrace", t_k, cfg.imaging.dose_rate))
    tau_dimer_true = 1.0 / k_diss(cfg, 2, t_k, p_o2)

    acc_rows, r2_rows, dec_rows, fc_rows = [], [], [], []
    for s in SEEDS:
        exp = K.build(cfg, seed=s, duration=DURATION, with_chip=False)
        tf = truth_frame(exp)
        d = exp.desc

        for name in ("N1", "N2", "N3", "N4p"):
            obs, tru = d[name].to_numpy(), tf[f"{name}_true"].to_numpy()
            m = np.isfinite(obs) & np.isfinite(tru)
            acc_rows.append({
                "seed": s, "descriptor": name,
                "true_mean": float(tru[m].mean()),
                "obs_mean": float(obs[m].mean()),
                "bias": float(obs[m].mean() - tru[m].mean()),
                "rel_bias": float((obs[m].mean() - tru[m].mean())
                                  / max(tru[m].mean(), 1e-9)),
                "pearson_r": float(np.corrcoef(obs[m], tru[m])[0, 1])
                if m.sum() > 3 else np.nan})
        # dispersed fraction, the quantity most exposed to detection bias
        tot_t = tf[["N1_true", "N2_true", "N3_true", "N4p_true"]].sum(axis=1)
        fd_t = (tf.N1_true / tot_t.clip(lower=1e-9)).to_numpy()
        fd_o = d.frac_dispersed.to_numpy()
        m = np.isfinite(fd_o) & np.isfinite(fd_t)
        acc_rows.append({"seed": s, "descriptor": "frac_dispersed",
                         "true_mean": float(fd_t[m].mean()),
                         "obs_mean": float(fd_o[m].mean()),
                         "bias": float(fd_o[m].mean() - fd_t[m].mean()),
                         "rel_bias": float((fd_o[m].mean() - fd_t[m].mean())
                                           / max(fd_t[m].mean(), 1e-9)),
                         "pearson_r": float(np.corrcoef(fd_o[m], fd_t[m])[0, 1])})
        for name, tv in (("k_hop", khop_true), ("tau_dimer", tau_dimer_true)):
            v = d[name].to_numpy()
            v = v[np.isfinite(v)]
            acc_rows.append({"seed": s, "descriptor": name,
                             "true_mean": tv, "obs_mean": float(v.mean()),
                             "bias": float(v.mean() - tv),
                             "rel_bias": float((v.mean() - tv) / tv),
                             "pearson_r": np.nan})

        # --- variance explained: census only vs census + dynamic ---
        y = tf.rate_true.to_numpy()
        sets = {
            "census_only": descriptors.CENSUS_COLS,
            "dynamic_only": descriptors.DYNAMIC_COLS,
            "census_plus_dynamic": (list(descriptors.CENSUS_COLS)
                                    + list(descriptors.DYNAMIC_COLS)),
            "census_matched_dim": list(descriptors.CENSUS_COLS)[:len(
                descriptors.DYNAMIC_COLS)],
        }
        for label, cols in sets.items():
            fit = correlate.fit_elastic_net(
                d, y, cols, force_cols=descriptors.REACTOR_COLS, cfg=cfg,
                seed=s)
            r2_rows.append({"seed": s, "descriptor_set": label,
                            "n_features": len(cols), "r2": fit["r2"],
                            "n": fit["n"]})

        # --- forecasting: where rate-type descriptors should earn their keep
        # The instantaneous rate is, by construction, a linear functional of
        # the instantaneous motif census, so census variables must explain it
        # almost perfectly and dynamic variables cannot add to that.  The
        # falsifiable claim of the central hypothesis is about PREDICTION:
        # transition rates and lifetimes should forecast where the catalyst is
        # going, which a census cannot.  Both tasks are therefore reported.
        for horizon in (30.0, 60.0, 120.0):
            nb = int(round(horizon / cfg.analysis.bin_s))
            if nb >= len(d) - 20:
                continue
            y_fut = np.full(len(d), np.nan)
            y_fut[:-nb] = y[nb:]
            # forecast the CHANGE, not the level: forecasting the level is
            # dominated by persistence and would hide the difference
            y_delta = y_fut - y
            for label, cols in (("census_only", descriptors.CENSUS_COLS),
                                ("census_plus_dynamic",
                                 list(descriptors.CENSUS_COLS)
                                 + list(descriptors.DYNAMIC_COLS))):
                fit = correlate.fit_elastic_net(
                    d, y_delta, cols, force_cols=descriptors.REACTOR_COLS,
                    cfg=cfg, seed=s)
                fc_rows.append({"seed": s, "horizon_s": horizon,
                                "descriptor_set": label, "r2": fit["r2"],
                                "n": fit["n"]})

        # --- decoder sensitivity ---
        for method, mf in (("median", 3), ("median", 5), ("viterbi", 3)):
            dd = descriptors.build_descriptors(
                cfg, exp.tracks, exp.frames, exp.q_pt, exp.schedule,
                dose_rate=exp.dose_rate, decode_method=method, min_frames=mf)
            tr = dd.attrs["transitions"]
            n_nuc_true = int(np.sum(exp.patch.events["kind"] == "nucleate"))
            td = dd.tau_dimer.to_numpy(); td = td[np.isfinite(td)]
            dec_rows.append({
                "seed": s, "decoder": f"{method}(w={mf})",
                "n_transitions": len(tr),
                "n_nucleations_obs": int(np.sum(tr.kind == "nucleate"))
                if len(tr) else 0,
                "n_nucleations_true": n_nuc_true,
                "tau_dimer_obs_s": float(td.mean()) if td.size else np.nan,
                "tau_dimer_true_s": tau_dimer_true})

    acc = pd.DataFrame(acc_rows)
    acc_s = acc.groupby("descriptor").agg(
        true_mean=("true_mean", "mean"), obs_mean=("obs_mean", "mean"),
        bias=("bias", "mean"), rel_bias=("rel_bias", "mean"),
        rel_bias_sd=("rel_bias", "std"),
        pearson_r=("pearson_r", "mean")).reset_index()
    r2 = pd.DataFrame(r2_rows)
    r2_s = r2.groupby("descriptor_set").agg(
        n_features=("n_features", "first"), r2_mean=("r2", "mean"),
        r2_sd=("r2", "std"), n_seeds=("seed", "count")).reset_index()
    fc = pd.DataFrame(fc_rows)
    fc_s = fc.groupby(["horizon_s", "descriptor_set"]).agg(
        r2_mean=("r2", "mean"), r2_sd=("r2", "std"),
        n_seeds=("seed", "count")).reset_index()
    piv = fc_s.pivot(index="horizon_s", columns="descriptor_set",
                     values="r2_mean").reset_index()
    piv["r2_gain_from_dynamic"] = (piv["census_plus_dynamic"]
                                   - piv["census_only"])

    dec = pd.DataFrame(dec_rows)
    dec_s = dec.groupby("decoder").agg(
        n_transitions=("n_transitions", "mean"),
        nucleations_obs=("n_nucleations_obs", "mean"),
        nucleations_true=("n_nucleations_true", "mean"),
        tau_dimer_obs_s=("tau_dimer_obs_s", "mean"),
        tau_dimer_true_s=("tau_dimer_true_s", "first")).reset_index()
    dec_s["nucleation_recall"] = (dec_s.nucleations_obs
                                  / dec_s.nucleations_true)
    dec_s["tau_dimer_rel_err"] = (dec_s.tau_dimer_obs_s
                                  / dec_s.tau_dimer_true_s - 1.0)

    report.save_table(acc_s, "table_S2a_descriptor_accuracy",
                      "Recovered versus true descriptors, averaged over "
                      f"{len(SEEDS)} independent {DURATION:.0f} s runs at the "
                      "reference dose. Bias, not noise, is the dominant error: "
                      "size-dependent detection systematically under-counts "
                      "monomers.", cfg=cfg, seed=seed)
    report.save_table(r2_s, "table_S2b_variance_explained",
                      "Variance in the true turnover rate explained by each "
                      "descriptor set, with reactor variables (T, p_CO, p_O2) "
                      "forced into every model so structural coefficients are "
                      "partial effects.", cfg=cfg, seed=seed)
    report.save_table(piv, "table_S2c_forecast_skill",
                      "Variance explained in the FUTURE CHANGE of the true "
                      "rate, r(t+h) - r(t), by census descriptors alone "
                      "versus census plus dynamic (rate- and lifetime-type) "
                      "descriptors. This is the prediction task in which "
                      "dynamic descriptors are expected to add information "
                      "that a census cannot carry.", cfg=cfg, seed=seed)
    report.save_table(dec_s, "table_S2d_decoder_sensitivity",
                      "Sensitivity of the dynamic descriptors to the size "
                      "decoder. The choice changes the observed transition "
                      "count by a factor of a few and must be pre-registered.",
                      cfg=cfg, seed=seed)
    report.save_json({"accuracy": acc, "r2": r2, "forecast": fc,
                      "decoder": dec}, "s2_descriptors", cfg=cfg, seed=seed)
    return {"accuracy": acc_s, "r2": r2_s, "forecast": piv,
            "decoder": dec_s}
