"""S8 -- does closing the loop help, and is the gain from intelligence or from
doing more? (gap G7)

Three policies share an identical actuation budget, so the comparison isolates
*when* the actuation is spent rather than *how much*.  The event-trigger
threshold is swept, because a badly chosen threshold spends the budget early
and gives away the advantage -- which is a result about the method, not a
failure to be hidden.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import control, report
from opcem.config import DEFAULT
from opcem.truth import PatchSimulator, Schedule

SELECT_SEEDS = tuple(range(701, 711))   # 10 threshold-selection seeds
EVAL_SEEDS = tuple(range(711, 725))     # 14 held-out evaluation seeds
SEEDS = SELECT_SEEDS + EVAL_SEEDS
DURATION = 900.0
THRESHOLDS = (1.0, 2.0, 3.0, 5.0, 8.0, 12.0)


def _run_one(c, seed, policy) -> dict:
    sch = Schedule(c)
    res = PatchSimulator(c, seed=seed, dose_rate=c.imaging.dose_rate,
                         schedule=sch).run(duration=DURATION,
                                           controller=policy)
    pop = res.pop
    disp = pop[:, 0] / np.maximum(pop.sum(axis=1), 1e-9)
    rew = control.reward(c, res.t, res.true_rate, policy.log.n_actuations,
                         c.imaging.dose_rate, float(disp[-20:].mean()))
    rew["dimer_atom_fraction"] = float(np.mean(
        2 * pop[:, 1] / np.maximum(
            pop[:, 0] + 2 * pop[:, 1] + 3 * pop[:, 2] + 5 * pop[:, 3], 1e-9)))
    rew["first_actuation_s"] = (float(policy.log.fired_t[0])
                                if policy.log.fired_t else np.nan)
    rew["median_actuation_s"] = (float(np.median(policy.log.fired_t))
                                 if policy.log.fired_t else np.nan)
    return rew


def run(cfg=DEFAULT, seed: int = 111):
    c = K.short_cfg(DURATION)
    rows = []
    for s in SEEDS:
        rows.append({"policy": "open_loop", "threshold": np.nan, "seed": s,
                     **_run_one(c, s, control.OpenLoop(c))})
        rows.append({"policy": "fixed_schedule", "threshold": np.nan,
                     "seed": s,
                     **_run_one(c, s, control.FixedSchedule(
                         c, duration_s=DURATION))})
        for th in THRESHOLDS:
            rows.append({"policy": "event_triggered", "threshold": th,
                         "seed": s,
                         **_run_one(c, s, control.EventTriggered(
                             c, threshold=th))})
    df = pd.DataFrame(rows)
    df["split"] = np.where(df.seed.isin(SELECT_SEEDS), "select", "eval")
    summ = df.groupby(["policy", "threshold"], dropna=False).agg(
        integrated_yield=("integrated_yield", "mean"),
        integrated_yield_sd=("integrated_yield", "std"),
        mean_rate=("mean_rate", "mean"),
        final_rate=("final_rate", "mean"),
        final_rate_sd=("final_rate", "std"),
        final_dispersed_frac=("final_dispersed_frac", "mean"),
        dimer_atom_fraction=("dimer_atom_fraction", "mean"),
        n_actuations=("n_actuations", "mean"),
        first_actuation_s=("first_actuation_s", "mean"),
        median_actuation_s=("median_actuation_s", "mean"),
        n_seeds=("seed", "count")).reset_index()
    base = summ.loc[summ.policy == "open_loop", "integrated_yield"].iloc[0]
    fixed = summ.loc[summ.policy == "fixed_schedule",
                     "integrated_yield"].iloc[0]
    summ["yield_vs_open_loop"] = summ.integrated_yield / base
    summ["yield_vs_fixed_schedule"] = summ.integrated_yield / fixed

    # Threshold selection and evaluation are separated.  The trigger
    # threshold is chosen on one half of the simulator seeds and the paired
    # comparison against the fixed schedule is then made on the other half.
    # Reporting the best-on-the-same-data threshold would be selection on the
    # test set, and the gain it shows is not an estimate of anything.
    from scipy import stats
    sel = df[df.split == "select"]
    et_sel = sel[sel.policy == "event_triggered"].groupby(
        "threshold").integrated_yield.mean()
    best_th = float(et_sel.idxmax())

    ev = df[df.split == "eval"]
    a = ev[(ev.policy == "event_triggered")
           & (ev.threshold == best_th)].sort_values("seed")
    b = ev[ev.policy == "fixed_schedule"].sort_values("seed")
    o = ev[ev.policy == "open_loop"].sort_values("seed")
    t_fs, p_fs = stats.ttest_rel(a.integrated_yield.to_numpy(),
                                 b.integrated_yield.to_numpy())
    t_ol, p_ol = stats.ttest_rel(a.integrated_yield.to_numpy(),
                                 o.integrated_yield.to_numpy())
    t_fr, p_fr = stats.ttest_rel(a.final_rate.to_numpy(),
                                 b.final_rate.to_numpy())
    d_fs = a.integrated_yield.to_numpy() - b.integrated_yield.to_numpy()
    # in-sample value, reported only so the optimism can be seen
    a_in = sel[(sel.policy == "event_triggered")
               & (sel.threshold == best_th)]
    b_in = sel[sel.policy == "fixed_schedule"]
    cmp = pd.DataFrame([{
        "threshold_selected_on": "select split (6 seeds)",
        "best_threshold": best_th,
        "in_sample_yield_ratio": float(a_in.integrated_yield.mean()
                                       / b_in.integrated_yield.mean()),
        "held_out_event_yield": float(a.integrated_yield.mean()),
        "held_out_fixed_yield": float(b.integrated_yield.mean()),
        "held_out_open_yield": float(o.integrated_yield.mean()),
        "held_out_yield_ratio_vs_fixed":
            float(a.integrated_yield.mean() / b.integrated_yield.mean()),
        "paired_mean_difference": float(d_fs.mean()),
        "paired_sd": float(d_fs.std(ddof=1)),
        "paired_t_vs_fixed": float(t_fs), "paired_p_vs_fixed": float(p_fs),
        "paired_t_vs_open": float(t_ol), "paired_p_vs_open": float(p_ol),
        "held_out_final_rate_event": float(a.final_rate.mean()),
        "held_out_final_rate_fixed": float(b.final_rate.mean()),
        "paired_t_final_rate": float(t_fr),
        "paired_p_final_rate": float(p_fr),
        "n_pairs": len(d_fs),
        # Minimum detectable effect at 80% power for this paired design, so a
        # null can be read as an exclusion of a given size rather than only as
        # a failure to reject.
        "mde_80pct_absolute": float(2.8 * d_fs.std(ddof=1)
                                    / np.sqrt(len(d_fs))),
        "mde_80pct_relative": float(2.8 * d_fs.std(ddof=1)
                                    / np.sqrt(len(d_fs))
                                    / b.integrated_yield.mean()),
        "actuations_event": float(a.n_actuations.mean()),
        "actuations_fixed": float(b.n_actuations.mean())}])

    report.save_table(summ, "table_S8a_control_policies",
                      "Closed-loop control against matched baselines. All "
                      "policies share the same actuation budget "
                      f"({c.control.actuation_budget} O2 pulses), so the "
                      "comparison isolates when the budget is spent. "
                      f"{len(SEEDS)} paired runs per row, split into "
                      f"{len(SELECT_SEEDS)} threshold-selection seeds and "
                      f"{len(EVAL_SEEDS)} held-out evaluation seeds.",
                      cfg=c, seed=seed)
    report.save_table(cmp, "table_S8b_control_paired_test",
                      "Held-out evaluation of the closed-loop policy. The "
                      "trigger threshold was chosen on the selection seeds "
                      "and the paired comparison made on disjoint evaluation "
                      "seeds; the in-sample ratio is shown alongside so the "
                      "selection optimism is visible.", cfg=c, seed=seed)
    report.save_json({"per_run": df, "summary": summ, "comparison": cmp},
                     "s8_control", cfg=c, seed=seed)
    return {"summary": summ, "comparison": cmp}
