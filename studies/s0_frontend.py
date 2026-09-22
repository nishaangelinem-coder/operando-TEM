"""S0 -- validation of the perception front end.

Checks that the fast analytic observation model used by every other study
reproduces what the pixel-level detector actually achieves, over a dose
series.  Any study that relied on the surrogate without this check would be
reporting the surrogate's assumptions rather than a measurement pipeline's
performance.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import common as K
from opcem import instruments, report, vision
from opcem.config import DEFAULT

DOSES = (3e2, 1e3, 3e3, 1e4)


def pixel_level(cfg, patch, dose, rng, crops, frames_idx, tol_nm=0.15):
    fov = cfg.imaging.render_field_nm
    n1_t = n1_d = n2_t = n2_d = 0
    n_det = n_match = 0
    loc, amp_err = [], []
    for k in frames_idx:
        for crop in crops:
            img = instruments.render_adf_frame(
                cfg, patch.pos3d[k], patch.sizes[k], rng, dose_rate=dose,
                crop_origin_nm=crop)
            txy, tsz = instruments.crop_truth(patch.pos3d[k], patch.sizes[k],
                                              crop, fov)
            det = vision.detect_atoms(cfg, img, dose_rate=dose)
            n_det += len(det.xy)
            if len(tsz) == 0:
                continue
            m = vision.match_detections(txy, det.xy, tol_nm=tol_nm)
            pairs = np.asarray(m["pairs"]).reshape(-1, 2)
            n_match += len(pairs)
            matched = set(pairs[:, 0].tolist())
            for i, sz in enumerate(tsz):
                if sz == 1:
                    n1_t += 1; n1_d += int(i in matched)
                else:
                    n2_t += 1; n2_d += int(i in matched)
            if len(pairs):
                loc.append(m["rmse_nm"])
                amp_err.append(det.amp[pairs[:, 1]] - tsz[pairs[:, 0]])
    ae = np.concatenate(amp_err) if amp_err else np.array([np.nan])
    return {"recall_n1": n1_d / max(n1_t, 1), "recall_n2p": n2_d / max(n2_t, 1),
            "precision": n_match / max(n_det, 1),
            "loc_rmse_pm": float(np.nanmean(loc) * 1000) if loc else np.nan,
            "amp_bias_atoms": float(np.nanmean(ae)),
            "amp_sd_atoms": float(np.nanstd(ae))}


def surrogate_level(cfg, patch, dose, seed):
    frames = instruments.observe_positions(
        cfg, patch, np.random.default_rng(seed), dose_rate=dose)
    n1_t = n1_d = n2_t = n2_d = 0
    tp = fp = 0
    bias = []
    for k, fr in enumerate(frames):
        tsz = np.asarray(patch.sizes[k])
        lut = {int(a): int(b) for a, b in zip(patch.ids[k], tsz)}
        real = fr.true_id >= 0
        tp += int(real.sum()); fp += int((~real).sum())
        seen = set(fr.true_id[real].tolist())
        for i, tid in enumerate(patch.ids[k]):
            if int(tid) in seen:
                if tsz[i] == 1:
                    n1_d += 1
                else:
                    n2_d += 1
        n1_t += int((tsz == 1).sum()); n2_t += int((tsz >= 2).sum())
        if real.any():
            bias.append(fr.amp[real]
                        - np.array([lut[int(x)] for x in fr.true_id[real]]))
    b = np.concatenate(bias) if bias else np.array([np.nan])
    return {"recall_n1": n1_d / max(n1_t, 1), "recall_n2p": n2_d / max(n2_t, 1),
            "precision": tp / max(tp + fp, 1),
            "loc_rmse_pm": instruments.localisation_sigma_nm(cfg, dose) * 1000,
            "amp_bias_atoms": float(np.nanmean(b)),
            "amp_sd_atoms": float(np.nanstd(b))}


def run(cfg=DEFAULT, seed: int = 21) -> pd.DataFrame:
    from opcem.truth import PatchSimulator
    patch = PatchSimulator(cfg, seed=seed, dose_rate=1e3).run(duration=120.0)
    crops = [(x, y) for x in (4.0, 16.0, 28.0) for y in (4.0, 16.0, 28.0)]
    frames_idx = [60, 180, 300, 420, 540]
    rows = []
    for dose in DOSES:
        rng = np.random.default_rng(seed * 31 + int(np.log10(dose) * 10))
        px = pixel_level(cfg, patch, dose, rng, crops, frames_idx)
        su = surrogate_level(cfg, patch, dose, seed * 17)
        for kk in px:
            rows.append({"dose_rate_e_A2_s": dose, "quantity": kk,
                         "pixel_level": px[kk], "surrogate": su[kk],
                         "abs_diff": abs(px[kk] - su[kk])})
        rows.append({"dose_rate_e_A2_s": dose,
                     "quantity": "single_atom_cnr",
                     "pixel_level": instruments.single_atom_cnr(cfg, dose),
                     "surrogate": instruments.single_atom_cnr(cfg, dose),
                     "abs_diff": 0.0})
    df = pd.DataFrame(rows)
    wide = df.pivot_table(index="dose_rate_e_A2_s", columns="quantity",
                          values=["pixel_level", "surrogate"])
    wide.columns = [f"{b}__{a}" for a, b in wide.columns]
    wide = wide.reset_index()
    report.save_table(df, "table_S0_frontend_validation",
                      "Pixel-level detector versus the fast analytic "
                      "observation model, over a dose series. The surrogate "
                      "was calibrated at the reference dose (1e3 e-/A2/s) "
                      "only; agreement away from it is a test, not a fit.",
                      cfg=cfg, seed=seed)
    report.save_json({"long": df, "wide": wide}, "s0_frontend", cfg=cfg,
                     seed=seed)
    return df


if __name__ == "__main__":
    print(run().to_string(index=False))
