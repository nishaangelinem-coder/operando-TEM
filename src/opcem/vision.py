"""Perception front end: atom detection from ADF frames, and linking.

Two entry points:

*   :func:`detect_atoms` works on rendered pixels and is what a real pipeline
    would run.  Its measured precision, recall and localisation error are
    reported in study S0 and are what justify the parameters of the fast
    analytic observation model in :mod:`opcem.instruments`.
*   :func:`link_frames` turns per-frame detections into tracks by gated
    nearest-neighbour assignment, and reports the fraction of links that are
    correct against the ground-truth identities, so that the descriptor
    uncertainties downstream are honest.

Cluster sizes are read from *integrated* ADF intensity rather than from
resolving the individual atoms in a sub-nanometre cluster, which is how atom
counting is done in practice.  The counting precision is therefore
shot-noise-limited and absolute in atom units.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage

from .config import Config, DEFAULT
from .instruments import (ObservedFrame, dose_per_frame, single_atom_cnr)


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #
@dataclass
class Detection:
    xy: np.ndarray          # (n, 2) nm
    amp: np.ndarray         # (n,) atoms (from integrated intensity)
    peak: np.ndarray        # (n,) peak significance in sigma


def detect_atoms(cfg: Config, img: np.ndarray,
                 dose_rate: Optional[float] = None,
                 fps: Optional[float] = None,
                 threshold_sigma: float = 5.0) -> Detection:
    """Matched-filter blob detection with sub-pixel centroid refinement.

    The matched filter is a Gaussian of the probe width, which is the optimal
    linear detector for a point source in additive noise; the detection
    threshold is set on the robust noise level of the filtered image, so the
    achieved false-positive rate follows from ``threshold_sigma``.
    Amplitudes are integrated intensities inside a 3-sigma disk with a local
    annulus background, converted to atom counts using the known dose and the
    ADF cross-section, i.e. by ordinary ADF atom counting.
    """
    im = cfg.imaging
    sig_pix = im.probe_sigma_nm / im.pixel_nm
    npix = img.shape[0]

    bg = ndimage.gaussian_filter(img, sigma=max(8.0 * sig_pix, 8.0))
    resid = img - bg
    filt = ndimage.gaussian_filter(resid, sigma=sig_pix)
    # At single-atom dose the frame holds <1 count per pixel, so a robust
    # (MAD-based) spread of the filtered image collapses to zero and grossly
    # under-estimates the noise.  The noise is therefore computed analytically
    # from the shot-noise variance and the matched filter's l2 norm,
    # sum(g^2) = 1/(4*pi*sigma^2); mean(img) is a conservative (upper-bound)
    # estimate of the per-pixel variance because it includes the atoms.
    var_per_pixel = max(float(np.mean(img)), 1e-12)
    noise = float(np.sqrt(var_per_pixel / (4.0 * np.pi * sig_pix ** 2)))

    fp_size = int(max(3, round(4 * sig_pix)))
    mx = ndimage.maximum_filter(filt, size=fp_size)
    peaks = (filt == mx) & (filt > threshold_sigma * noise)
    ys, xs = np.nonzero(peaks)
    if xs.size == 0:
        return Detection(np.zeros((0, 2)), np.zeros(0), np.zeros(0))

    d = dose_per_frame(cfg, dose_rate, fps)
    counts_per_atom = im.sigma_adf_pt_A2 * d
    r_in = max(2, int(np.ceil(3.0 * sig_pix)))
    r_out = r_in + 3
    yy0, xx0 = np.mgrid[-r_out:r_out + 1, -r_out:r_out + 1]
    rr = np.hypot(xx0, yy0)
    disk, ring = rr <= r_in, (rr > r_in) & (rr <= r_out)

    out_xy, out_amp, out_pk = [], [], []
    for y0, x0 in zip(ys, xs):
        if not (r_out <= y0 < npix - r_out and r_out <= x0 < npix - r_out):
            continue
        patch = img[y0 - r_out:y0 + r_out + 1, x0 - r_out:x0 + r_out + 1]
        # The annulus background must be estimated with the MEAN, not the
        # median: for sparse Poisson data the median is 0 and the aperture sum
        # is then biased upward by the whole background it contains.
        local_bg = float(np.mean(patch[ring]))
        amp = float(patch[disk].sum() - disk.sum() * local_bg) \
            / max(counts_per_atom, 1e-12)
        # centroid on the matched-filtered, background-free image
        fpatch = filt[y0 - r_in:y0 + r_in + 1, x0 - r_in:x0 + r_in + 1]
        wgt = np.clip(fpatch, 0.0, None)
        if wgt.sum() <= 0:
            continue
        gy, gx = np.mgrid[y0 - r_in:y0 + r_in + 1, x0 - r_in:x0 + r_in + 1]
        cy = float((wgt * gy).sum() / wgt.sum())
        cx = float((wgt * gx).sum() / wgt.sum())
        out_xy.append((cx * im.pixel_nm, cy * im.pixel_nm))
        out_amp.append(amp)
        out_pk.append(float(filt[y0, x0] / noise))
    return Detection(np.asarray(out_xy).reshape(-1, 2),
                     np.asarray(out_amp), np.asarray(out_pk))


def match_detections(truth_xy: np.ndarray, det_xy: np.ndarray,
                     tol_nm: float = 0.25) -> Dict[str, object]:
    """Greedy one-to-one matching of detections to ground truth."""
    if det_xy.size == 0 or truth_xy.size == 0:
        return {"precision": 0.0, "recall": 0.0, "rmse_nm": np.nan,
                "n_match": 0, "pairs": []}
    d = np.linalg.norm(truth_xy[:, None, :] - det_xy[None, :, :], axis=2)
    pairs: List[Tuple[int, int]] = []
    used_t, used_d = set(), set()
    order = np.dstack(np.unravel_index(np.argsort(d, axis=None), d.shape))[0]
    for i, j in order:
        if d[i, j] > tol_nm:
            break
        if i in used_t or j in used_d:
            continue
        used_t.add(int(i)); used_d.add(int(j)); pairs.append((int(i), int(j)))
    if not pairs:
        return {"precision": 0.0, "recall": 0.0, "rmse_nm": np.nan,
                "n_match": 0, "pairs": []}
    err = np.array([d[i, j] for i, j in pairs])
    return {"precision": len(pairs) / len(det_xy),
            "recall": len(pairs) / len(truth_xy),
            "rmse_nm": float(np.sqrt(np.mean(err ** 2))),
            "n_match": len(pairs), "pairs": pairs}


# --------------------------------------------------------------------------- #
# linking / tracking
# --------------------------------------------------------------------------- #
@dataclass
class Tracks:
    """Per-frame track assignment.

    ``label[k]`` gives the track label of each detection in frame ``k``
    (-1 = unassigned).  ``xy[k]``, ``amp[k]`` mirror the observed frames.
    """
    t: np.ndarray
    label: List[np.ndarray]
    xy: List[np.ndarray]
    amp: List[np.ndarray]
    link_accuracy: float
    n_tracks: int


def link_frames(frames: Sequence[ObservedFrame], max_disp_nm: float = 0.8,
                score_truth: bool = True) -> Tracks:
    """Gated nearest-neighbour linking between consecutive frames.

    A detection is linked to the closest detection in the previous frame within
    ``max_disp_nm``; otherwise a new track is started.  ``link_accuracy`` is
    the fraction of links that connect detections with the same ground-truth
    identity, and is reported so that descriptor error budgets are explicit.
    """
    labels: List[np.ndarray] = []
    next_label = 0
    good = bad = 0
    prev_xy = prev_lab = prev_id = None
    for fr in frames:
        n = len(fr.xy)
        lab = -np.ones(n, dtype=np.int64)
        if prev_xy is not None and n and len(prev_xy):
            d = np.linalg.norm(fr.xy[:, None, :] - prev_xy[None, :, :], axis=2)
            taken = set()
            order = np.dstack(np.unravel_index(np.argsort(d, axis=None),
                                               d.shape))[0]
            for i, j in order:
                if d[i, j] > max_disp_nm:
                    break
                if lab[i] != -1 or j in taken:
                    continue
                lab[i] = prev_lab[j]
                taken.add(int(j))
                if score_truth:
                    if fr.true_id[i] == prev_id[j] and fr.true_id[i] >= 0:
                        good += 1
                    else:
                        bad += 1
        for i in range(n):
            if lab[i] == -1:
                lab[i] = next_label
                next_label += 1
        labels.append(lab)
        prev_xy, prev_lab, prev_id = fr.xy, lab, fr.true_id
    acc = good / (good + bad) if (good + bad) else np.nan
    return Tracks(t=np.array([f.t for f in frames]), label=labels,
                  xy=[f.xy for f in frames], amp=[f.amp for f in frames],
                  link_accuracy=float(acc), n_tracks=next_label)
