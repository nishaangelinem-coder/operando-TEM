"""Shared construction of one synchronised virtual experiment.

Produces, on a single time base:

*   the patch kMC ground truth (the only level with atom trajectories),
*   the observed detections, tracks and descriptor table,
*   intermittent 3D snapshots,
*   the chip-level mean-field rate that the mass spectrometer integrates,
*   two MS traces: the realistic one (chip rate) and an explicitly
    counterfactual one in which the entire active area IS the imaged field.

The counterfactual trace is the object that most of the estimator validation
needs, and it is labelled as a counterfactual everywhere it is used, because
the representativeness fraction of a real chip makes it physically
unavailable (study S7 quantifies by how much).
"""
from __future__ import annotations

import copy
import json
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

from opcem import (beam, control, correlate, descriptors, instruments,
                   mobility, report, transport, truth, vision)   # noqa: E402
from opcem.config import Config, DEFAULT                          # noqa: E402

K_ENC_CACHE = os.path.join(report.RAW, "k_enc_calibration.json")


def k_enc(cfg: Config = DEFAULT, recompute: bool = False) -> float:
    """Load (or compute once and cache) the chip model's encounter coefficient."""
    if not recompute and os.path.exists(K_ENC_CACHE):
        with open(K_ENC_CACHE) as fh:
            d = json.load(fh)
        return float(d.get("k_enc", d.get("data", {}).get("k_enc")))
    cal = truth.calibrate_k_enc(cfg)
    os.makedirs(report.RAW, exist_ok=True)
    payload = {"_provenance": {"simulated": True,
                               "banner": report.SIMULATED_BANNER,
                               "git": report.git_sha(), "seed": cfg.seed},
               **cal}
    with open(K_ENC_CACHE, "w") as fh:
        json.dump(payload, fh, indent=2)
    return float(cal["k_enc"])


@dataclass
class VirtualExperiment:
    cfg: Config
    seed: int
    dose_rate: float
    patch: truth.PatchResult
    frames: List[instruments.ObservedFrame]
    tracks: vision.Tracks
    q_pt: np.ndarray
    desc: pd.DataFrame
    snapshots: List[Dict[str, object]]
    schedule: truth.Schedule
    t_ms: np.ndarray
    kernel: np.ndarray
    ms_chip: instruments.MSTrace
    ms_patch_cf: instruments.MSTrace
    chip: Dict[str, np.ndarray]

    @property
    def dt_ms(self) -> float:
        return float(np.mean(np.diff(self.t_ms)))


def build(cfg: Config = DEFAULT, seed: int = 1,
          duration: Optional[float] = None,
          dose_rate: Optional[float] = None,
          schedule: Optional[truth.Schedule] = None,
          controller=None,
          with_chip: bool = True,
          decode_method: str = "median",
          min_frames: int = 3) -> VirtualExperiment:
    dur = cfg.duration_s if duration is None else duration
    dr = cfg.imaging.dose_rate if dose_rate is None else dose_rate
    sch = schedule or truth.Schedule(cfg)
    rng = np.random.default_rng(seed * 7919 + 13)

    patch = truth.PatchSimulator(cfg, seed=seed, dose_rate=dr,
                                 schedule=sch).run(duration=dur,
                                                   controller=controller)
    frames = instruments.observe_positions(cfg, patch, rng, dose_rate=dr)
    q = instruments.eels_oxidation_proxy(cfg, patch, rng, dose_rate=dr)
    tracks = vision.link_frames(frames)
    desc = descriptors.build_descriptors(cfg, tracks, frames, q, sch,
                                         dose_rate=dr,
                                         decode_method=decode_method,
                                         min_frames=min_frames)
    snaps = instruments.snapshot3d(cfg, patch, rng)

    dt_ms = 1.0 / cfg.ms.rate_hz
    t_ms = np.arange(0.0, dur, dt_ms)
    kern = transport.kernel_from_config(cfg, dt_ms)

    patch_rate_ms = np.interp(t_ms, patch.t, patch.true_rate)
    if with_chip:
        # The chip mean field is smooth on the MS sampling interval, so it is
        # integrated on a 1 Hz grid and interpolated up; this is ~10x cheaper
        # and the interpolation error is far below the MS noise floor.
        t_chip = np.arange(0.0, dur + 1.0, 1.0)
        chip = truth.chip_mean_field(cfg, t_chip, sch, k_enc=k_enc(cfg),
                                     dose_rate=0.0)
        chip_rate = np.interp(t_ms, t_chip, chip["rate"])
        chip["rate_ms"] = chip_rate
    else:
        chip = {"t": t_ms, "rate": patch_rate_ms}
        chip_rate = patch_rate_ms

    ms_chip = instruments.virtual_qms(cfg, t_ms, chip_rate, kern,
                                      np.random.default_rng(seed * 104729 + 7))
    ms_cf = instruments.virtual_qms(cfg, t_ms, patch_rate_ms, kern,
                                    np.random.default_rng(seed * 104729 + 11))
    return VirtualExperiment(cfg=cfg, seed=seed, dose_rate=dr, patch=patch,
                             frames=frames, tracks=tracks, q_pt=q, desc=desc,
                             snapshots=snaps, schedule=sch, t_ms=t_ms,
                             kernel=kern, ms_chip=ms_chip,
                             ms_patch_cf=ms_cf, chip=chip)


def recover_rate(exp: VirtualExperiment, which: str = "counterfactual",
                 method: str = "deconvolve",
                 tau_dead: Optional[float] = None,
                 kernel: Optional[np.ndarray] = None) -> Dict[str, object]:
    """Invert an MS channel back to an estimated formation rate."""
    cfg = exp.cfg
    trace = exp.ms_patch_cf if which == "counterfactual" else exp.ms_chip
    td = cfg.transport.tau_dead() if tau_dead is None else tau_dead
    kern = exp.kernel if kernel is None else kernel
    return correlate.rate_from_ms(cfg, trace.t, trace.ion_current["CO2"],
                                  cfg.ms.k_cal["CO2"],
                                  cfg.ms.background["CO2"], kern, td,
                                  method=method)


def truth_event_times(patch: truth.PatchResult, kind: str = "nucleate"
                      ) -> np.ndarray:
    ev = patch.events
    return np.sort(ev["t"][ev["kind"] == kind].astype(float))


def observed_event_times(desc: pd.DataFrame, kind: str = "nucleate"
                         ) -> np.ndarray:
    tr = desc.attrs.get("transitions")
    if tr is None or not len(tr):
        return np.zeros(0)
    return np.sort(tr.loc[tr.kind == kind, "t"].to_numpy(dtype=float))


def short_cfg(duration: float, n_atoms: Optional[int] = None) -> Config:
    c = copy.deepcopy(DEFAULT)
    c.duration_s = duration
    if n_atoms is not None:
        c.n_pt_atoms = n_atoms
    return c


def sparse_cfg(duration: float, n_atoms: int = 10,
               pre_s: float = 4.0, post_s: float = 4.0) -> Config:
    """Configuration in which single atomic events are temporally resolvable.

    Event-aligned analysis requires the event rate to be small compared with
    the reciprocal of the analysis window: with lambda*(t_pre+t_post) of order
    one, the pre-window of one event overlaps the post-windows of many others
    and no single event's contribution survives.  Because the nucleation rate
    scales roughly as the square of the monomer count, the criterion is a
    constraint on Pt loading -- which is a design prescription for the
    experiment, not a property of the estimator.  Study S3 measures where the
    boundary lies.
    """
    c = short_cfg(duration, n_atoms)
    c.analysis.pre_window_s = pre_s
    c.analysis.post_window_s = post_s
    return c


def event_rate(patch, kind: str = "nucleate") -> float:
    """Observed ground-truth event rate (events per second)."""
    ev = patch.events
    n = int(np.sum(ev["kind"] == kind))
    span = float(patch.t[-1] - patch.t[0]) if len(patch.t) > 1 else 1.0
    return n / max(span, 1e-9)
