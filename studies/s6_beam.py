"""S6 -- separating chemistry from irradiation (gap G4).

A dose-rate series at fixed temperature, pressure and flow.  Three results,
two of them cautionary.

**The usable dose window is bounded from below by detection, not by damage.**
Below roughly 10^3 e-/A^2/s the matched-filter detector runs close to its noise
floor, false positives outnumber real atoms across a 40 nm field, and a few
chains of coincidentally linked false positives dominate the mean squared
displacement.  The persistence filter removes most of them but not all.  A
diagnostic that needs no ground truth -- the fraction of detections belonging to
a track of at least five frames -- identifies the unusable doses, and the
extrapolation is restricted accordingly.  Because beam damage bounds the window
from above, the accessible range is about one decade, which is what limits the
precision of any zero-dose extrapolation.

**The reference is the effective hop rate, not the terrace rate.** A monomer
occupies each support site type in proportion to ``p_s / k_s``, so the
ensemble hop rate an MSD estimator measures is a harmonic mean over site types.

**The linear dose model is therefore misspecified.** A harmonic mean of
``k_s + beta*phi`` is sublinear in dose, so a straight-line fit biases both the
zero-dose intercept and the Beam Perturbation Index.  Fitting the correct
site-aware model -- with the barrier distribution taken as known from a
beam-off temperature series, which is how it would be obtained experimentally
-- removes that bias.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import beam, report
from opcem.config import DEFAULT
from opcem.truth import k_hop

DOSES = (1e2, 3e2, 1e3, 2e3, 3e3, 6e3, 1e4, 2e4)
PERSISTENCE_GATE = 0.80
SEEDS = (501, 502, 503)
DURATION = 400.0


def run(cfg=DEFAULT, seed: int = 81):
    c = K.short_cfg(DURATION)
    t_k = c.reactor.t_setpoint_c + 273.15

    per = []
    for dose in DOSES:
        for s in SEEDS:
            exp = K.build(c, seed=s, duration=DURATION, dose_rate=dose,
                          with_chip=False)
            kh = exp.desc.k_hop.to_numpy(); kh = kh[np.isfinite(kh)]
            kn = exp.desc.k_nuc.to_numpy(); kn = kn[np.isfinite(kn)]
            fd = exp.desc.frac_dispersed.to_numpy()
            fd = fd[np.isfinite(fd)]
            ev = exp.patch.events
            nuc = ev[ev["kind"] == "nucleate"]
            hc = exp.patch.hop_counts
            per.append({
                "dose_rate": dose, "seed": s,
                "persistent_frac": float(exp.desc.persistent_frac.iloc[0]),
                "k_hop_obs": float(kh.mean()) if kh.size else np.nan,
                "k_nuc_obs": float(kn.mean()) if kn.size else np.nan,
                "frac_dispersed_obs": float(fd.mean()) if fd.size else np.nan,
                "k_eff_true": beam.effective_hop_rate(c, t_k, dose),
                "k_terrace_true": float(sum(k_hop(c, "terrace", t_k, dose))),
                "bpi_analytic": beam.beam_hop_fraction(c, t_k, dose),
                "bpi_simulated": float(hc["beam"]
                                       / max(sum(hc.values()), 1)),
                "bpi_simulated_nucleations":
                    float(nuc["beam"].mean()) if len(nuc) else np.nan})
    df = pd.DataFrame(per)

    def sem(v):
        v = np.asarray(v, dtype=float)
        return float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 \
            else np.nan

    g = df.groupby("dose_rate").agg(
        persistent_frac=("persistent_frac", "mean"),
        k_hop_obs=("k_hop_obs", "mean"), k_hop_sem=("k_hop_obs", sem),
        k_nuc_obs=("k_nuc_obs", "mean"), k_nuc_sem=("k_nuc_obs", sem),
        frac_dispersed=("frac_dispersed_obs", "mean"),
        frac_dispersed_sem=("frac_dispersed_obs", sem),
        k_eff_true=("k_eff_true", "first"),
        k_terrace_true=("k_terrace_true", "first"),
        bpi_analytic=("bpi_analytic", "first"),
        bpi_simulated=("bpi_simulated", "mean")).reset_index()
    g["usable"] = g.persistent_frac >= PERSISTENCE_GATE
    g["k_hop_rel_err_vs_keff"] = g.k_hop_obs / g.k_eff_true - 1.0

    use = g[g.usable].copy()
    fit_lin = beam.fit_dose_series(use.dose_rate, use.k_hop_obs,
                                   use.k_hop_sem)
    fit_site = beam.fit_beta_site_aware(c, t_k, use.dose_rate.to_numpy(),
                                        use.k_hop_obs.to_numpy())
    fd_ex = beam.beam_free_extrapolation(use.dose_rate, use.frac_dispersed,
                                         use.frac_dispersed_sem)
    k_eff_0 = beam.effective_hop_rate(c, t_k, 0.0)

    fits = pd.DataFrame([
        {"model": "linear in dose (standard BPI assumption)",
         "zero_dose_value": fit_lin.k_chem, "se": fit_lin.k_chem_se,
         "reference": k_eff_0,
         "rel_error": fit_lin.k_chem / k_eff_0 - 1.0,
         "beta_fit": fit_lin.beta, "beta_se": fit_lin.beta_se,
         "beta_true": c.kinetics.beta_beam, "r2_or_rmse": fit_lin.r2},
        {"model": "site-aware harmonic mean (barriers known)",
         "zero_dose_value": fit_site["k_eff_zero_dose"]
         * fit_site["detection_gain"], "se": np.nan,
         "reference": k_eff_0,
         "rel_error": (fit_site["k_eff_zero_dose"]
                       * fit_site["detection_gain"]) / k_eff_0 - 1.0,
         "beta_fit": fit_site["beta"], "beta_se": fit_site["beta_se"],
         "beta_true": c.kinetics.beta_beam,
         "r2_or_rmse": fit_site["rmse"]},
        {"model": "dispersed fraction, linear extrapolation",
         "zero_dose_value": fd_ex["value_at_zero_dose"], "se": fd_ex["se"],
         "reference": np.nan, "rel_error": np.nan,
         "beta_fit": fd_ex["slope"], "beta_se": fd_ex["slope_se"],
         "beta_true": np.nan, "r2_or_rmse": fd_ex["r2"]}])

    bpi_rows = []
    for dose in use.dose_rate:
        lo, hi = fit_lin.bpi_ci(dose, seed=seed)
        row = g.loc[g.dose_rate == dose].iloc[0]
        b_site = beam.beam_hop_fraction(c, t_k, dose)   # with fitted beta
        import copy as _copy
        c2 = _copy.deepcopy(c)
        c2.kinetics.beta_beam = fit_site["beta"]
        bpi_rows.append({
            "dose_rate": dose,
            "bpi_linear_fit": fit_lin.bpi(dose),
            "bpi_linear_ci_lo": lo, "bpi_linear_ci_hi": hi,
            "bpi_site_aware_fit": beam.beam_hop_fraction(c2, t_k, dose),
            "bpi_analytic_truth": row.bpi_analytic,
            "bpi_simulated_truth": row.bpi_simulated,
            "included_at_threshold":
                bool(beam.beam_hop_fraction(c2, t_k, dose)
                     <= cfg.analysis.bpi_exclusion_threshold)})
    bpi = pd.DataFrame(bpi_rows)

    report.save_table(g, "table_S6a_dose_series",
                      "Observed atomic-process rates against electron dose "
                      f"rate at fixed T, p and flow ({len(SEEDS)} runs per "
                      "point). 'persistent_frac' is the ground-truth-free "
                      "contamination diagnostic; doses below the "
                      f"{PERSISTENCE_GATE:.2f} gate are excluded because false "
                      "positives dominate the tracking. 'k_eff_true' is the "
                      "site-occupancy-weighted effective hop rate, which is "
                      "what an MSD estimator measures.", cfg=c, seed=seed)
    report.save_table(fits, "table_S6b_dose_models",
                      "Zero-dose extrapolation under the standard "
                      "linear-in-dose assumption and under the correct "
                      "site-aware harmonic-mean model, fitted only to the "
                      "usable dose window.", cfg=c, seed=seed)
    report.save_table(bpi, "table_S6c_beam_perturbation_index",
                      "Beam Perturbation Index from both dose models against "
                      "the analytic truth and the simulator's own "
                      "beam-attribution counter. The last column applies the "
                      f"pre-registered threshold of "
                      f"{cfg.analysis.bpi_exclusion_threshold}.",
                      cfg=c, seed=seed)
    report.save_json({"per_seed": df, "grouped": g, "fits": fits, "bpi": bpi,
                      "usable_gate": PERSISTENCE_GATE,
                      "k_eff_zero_dose": k_eff_0},
                     "s6_beam", cfg=c, seed=seed)
    return {"dose_series": g, "fits": fits, "bpi": bpi}
