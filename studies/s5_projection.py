"""S5 -- how biased is mobility measured from projected images? (gap G3)

The continuous image stream has no depth information, so the mean-squared
displacement it reports is the projected one.  The bias is quantified against
the simulator's own 3D coordinates, and the projection-aware correction -- whose
only 3D input is the sparse snapshot stream, never an interpolation of it -- is
validated.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import mobility, report
from opcem.config import DEFAULT

SEEDS = (401, 402, 403, 404)
DURATION = 400.0
AMPLITUDES = (0.0, 0.35, 0.70, 1.05)


def run(cfg=DEFAULT, seed: int = 71):
    rows = []
    for amp in AMPLITUDES:
        c = K.short_cfg(DURATION)
        c.support.corrug_amp_nm = amp
        eta_an = mobility.eta_analytic(c)
        for s in SEEDS:
            exp = K.build(c, seed=s, duration=DURATION, with_chip=False)
            gt = mobility.msd_3d_truth(exp.patch.pos3d, exp.patch.ids,
                                       exp.patch.sizes, lag_frames=1,
                                       box_nm=c.support.field_nm)
            est = mobility.eta_from_snapshots(c, exp.snapshots)
            obs_msd = exp.desc.msd_1s.to_numpy()
            obs_msd = float(np.nanmean(obs_msd))
            corr = mobility.correct_msd(obs_msd, est["eta"],
                                        est.get("eta_se", 0.0))
            rows.append({
                "corrug_amp_nm": amp, "seed": s,
                "eta_analytic": eta_an,
                "eta_truth_msd": gt["eta"],
                "eta_from_snapshots": est["eta"],
                "eta_snapshot_se": est.get("eta_se", np.nan),
                "msd_3d_truth_nm2": gt["msd_3d"],
                "msd_2d_truth_nm2": gt["msd_2d"],
                "projection_bias_frac": (gt["msd_2d"] - gt["msd_3d"])
                / gt["msd_3d"] if gt["msd_3d"] else np.nan,
                "msd_2d_observed_nm2": obs_msd,
                "msd_3d_corrected_nm2": corr["msd_3d_est"],
                "corrected_error_frac": (corr["msd_3d_est"] - gt["msd_3d"])
                / gt["msd_3d"] if gt["msd_3d"] else np.nan,
                "n_pairs": gt["n"], "n_snapshots": len(exp.snapshots),
                "sigma_z_nm": exp.snapshots[0]["sigma_z_nm"]
                if exp.snapshots else np.nan})
    df = pd.DataFrame(rows)
    summ = df.groupby("corrug_amp_nm").agg(
        eta_analytic=("eta_analytic", "first"),
        eta_truth_msd=("eta_truth_msd", "mean"),
        eta_from_snapshots=("eta_from_snapshots", "mean"),
        eta_snapshot_sd=("eta_from_snapshots", "std"),
        projection_bias_frac=("projection_bias_frac", "mean"),
        msd_2d_observed_nm2=("msd_2d_observed_nm2", "mean"),
        msd_3d_truth_nm2=("msd_3d_truth_nm2", "mean"),
        corrected_error_frac=("corrected_error_frac", "mean"),
        corrected_error_sd=("corrected_error_frac", "std"),
        sigma_z_nm=("sigma_z_nm", "first")).reset_index()
    report.save_table(summ, "table_S5_projection_bias",
                      "Projected versus true 3D single-atom mean-squared "
                      "displacement, against support corrugation amplitude. "
                      "'eta_from_snapshots' is estimated only from the "
                      "intermittent 3D stream, whose depth precision is "
                      "sigma_z; the corrected error column is what remains "
                      "after applying it.", cfg=cfg, seed=seed)
    report.save_json({"per_seed": df, "summary": summ}, "s5_projection",
                     cfg=cfg, seed=seed)
    return summ
