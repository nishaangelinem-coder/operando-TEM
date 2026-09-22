"""Assembly of the descriptor vector D(t) from observed streams.

Design rules enforced here (docs/04_experimental_setup.md section 4):

1.  every column carries a provenance tag from :data:`opcem.config.PROVENANCE`;
2.  the intermittent 3D stream is never interpolated onto the continuous grid;
3.  rate-type descriptors are estimated with an explicit noise correction, so
    that localisation error does not masquerade as atom mobility.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import Config, DEFAULT, PROVENANCE
from .instruments import ObservedFrame, localisation_sigma_nm
from .truth import Schedule
from .vision import Tracks


def sizes_from_amp(amp: np.ndarray) -> np.ndarray:
    """Assign an integer atom count to each detection from its ADF amplitude."""
    return np.clip(np.rint(amp), 1, None).astype(int)


# --------------------------------------------------------------------------- #
# mobility
# --------------------------------------------------------------------------- #
def msd_from_tracks(tracks: Tracks, lag_frames: int = 1,
                    monomer_only: bool = True) -> Tuple[float, int]:
    """Mean-squared *projected* displacement at a given frame lag.

    Returns (MSD in nm^2, number of displacement pairs).  This is a projected
    quantity: the image stream cannot report the out-of-plane component
    (gap G3), so this systematically under-estimates the 3D MSD.
    """
    by_label: Dict[int, List[Tuple[int, np.ndarray, float]]] = {}
    for k in range(len(tracks.label)):
        s = sizes_from_amp(tracks.amp[k])
        for i, lab in enumerate(tracks.label[k]):
            by_label.setdefault(int(lab), []).append(
                (k, tracks.xy[k][i], float(s[i])))
    vals: List[float] = []
    for _lab, seq in by_label.items():
        seq.sort(key=lambda z: z[0])
        for a in range(len(seq) - lag_frames):
            k0, p0, s0 = seq[a]
            k1, p1, s1 = seq[a + lag_frames]
            if k1 - k0 != lag_frames:
                continue
            if monomer_only and not (s0 == 1 and s1 == 1):
                continue
            vals.append(float(np.sum((p1 - p0) ** 2)))
    if not vals:
        return float("nan"), 0
    return float(np.mean(vals)), len(vals)


def k_hop_from_msd(cfg: Config, msd_nm2: float, dt_s: float,
                   dose_rate: Optional[float] = None) -> float:
    """Invert a 2-D random walk for the hop rate, correcting localisation noise.

    For an unbiased square-lattice walk of step ``a``, MSD(dt) = k_hop a^2 dt.
    Independent localisation error adds ``2 sigma_loc^2`` to every measured
    squared displacement; omitting that term is the standard way to
    manufacture spurious mobility at low dose.
    """
    a = cfg.support.lattice_nm
    sig = localisation_sigma_nm(cfg, dose_rate)
    corrected = msd_nm2 - 2.0 * sig ** 2
    return float(max(corrected, 0.0) / (a ** 2 * dt_s))


# --------------------------------------------------------------------------- #
# transitions
# --------------------------------------------------------------------------- #
def viterbi_sizes(amp: np.ndarray, cnr1: float, n_max: int = 6,
                  switch_log_prior: float = -6.0) -> np.ndarray:
    """Decode the most likely integer size sequence of one track.

    Emission log-likelihood of observing amplitude ``a`` from an ``n``-atom
    entity is ``-0.5 * ((a - n) * cnr1)**2`` because the atom-counting error is
    Gaussian with standard deviation ``1/cnr1`` atoms.  A constant
    ``switch_log_prior`` penalises any size change, which is what suppresses
    the noise-driven 1<->2 flicker that a per-frame threshold produces.

    A temporal median filter is the cruder alternative; the two are compared in
    study S2 because the choice changes the dynamic descriptors by a factor of
    a few, and that sensitivity has to be reported rather than buried.
    """
    a = np.asarray(amp, dtype=float)
    states = np.arange(1, n_max + 1)
    n_t, n_s = len(a), len(states)
    if n_t == 0:
        return np.zeros(0, dtype=int)
    emis = -0.5 * ((a[:, None] - states[None, :]) * cnr1) ** 2
    delta = np.empty((n_t, n_s)); psi = np.zeros((n_t, n_s), dtype=int)
    delta[0] = emis[0]
    for k in range(1, n_t):
        trans = np.full((n_s, n_s), switch_log_prior)
        np.fill_diagonal(trans, 0.0)
        tot = delta[k - 1][:, None] + trans
        psi[k] = np.argmax(tot, axis=0)
        delta[k] = emis[k] + tot[psi[k], np.arange(n_s)]
    path = np.empty(n_t, dtype=int)
    path[-1] = int(np.argmax(delta[-1]))
    for k in range(n_t - 2, -1, -1):
        path[k] = psi[k + 1, path[k + 1]]
    return states[path]


def decode_tracks(tracks: Tracks, cnr1: float, method: str = "viterbi",
                  min_frames: int = 3,
                  switch_log_prior: float = -6.0
                  ) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """Return {track label: (times, decoded integer sizes)}."""
    by_label: Dict[int, List[Tuple[float, float]]] = {}
    for k in range(len(tracks.label)):
        for i, lab in enumerate(tracks.label[k]):
            by_label.setdefault(int(lab), []).append(
                (float(tracks.t[k]), float(tracks.amp[k][i])))
    out: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for lab, seq in by_label.items():
        seq.sort(key=lambda z: z[0])
        ts = np.array([z[0] for z in seq])
        amp = np.array([z[1] for z in seq])
        if method == "viterbi":
            ss = viterbi_sizes(amp, cnr1, switch_log_prior=switch_log_prior)
        elif method == "median":
            raw = sizes_from_amp(amp)
            w = min_frames
            ss = np.array([int(np.median(raw[max(0, i - w // 2):i + w // 2 + 1]))
                           for i in range(len(raw))])
        else:
            raise ValueError(f"unknown decode method {method!r}")
        out[lab] = (ts, ss)
    return out


def transitions_from_decoded(
        decoded: Dict[int, Tuple[np.ndarray, np.ndarray]]) -> pd.DataFrame:
    """Read size transitions off decoded per-track size sequences."""
    rows = []
    for lab, (ts, ss) in decoded.items():
        ch = np.flatnonzero(np.diff(ss) != 0)
        for i in ch:
            before, after = int(ss[i]), int(ss[i + 1])
            rows.append({"t": float(ts[i + 1]), "label": lab,
                         "size_from": before, "size_to": after,
                         "kind": "nucleate" if (before == 1 and after == 2)
                         else ("grow" if after > before else "shrink")})
    return pd.DataFrame(rows, columns=["t", "label", "size_from", "size_to",
                                       "kind"]).sort_values("t").reset_index(
        drop=True)


def detect_transitions(tracks: Tracks, min_frames: int = 2
                       ) -> pd.DataFrame:
    """Observed size transitions of each track (the estimator's 'events').

    A transition is recorded when a track's median size over ``min_frames``
    consecutive frames differs from its median over the preceding
    ``min_frames``.  Median filtering is essential: single-frame amplitude
    noise of ~0.3 atoms would otherwise generate a transition every few
    frames.
    """
    by_label: Dict[int, List[Tuple[float, int]]] = {}
    for k in range(len(tracks.label)):
        s = sizes_from_amp(tracks.amp[k])
        for i, lab in enumerate(tracks.label[k]):
            by_label.setdefault(int(lab), []).append((float(tracks.t[k]),
                                                      int(s[i])))
    rows = []
    for lab, seq in by_label.items():
        if len(seq) < 2 * min_frames:
            continue
        seq.sort(key=lambda z: z[0])
        ts = np.array([z[0] for z in seq])
        ss = np.array([z[1] for z in seq])
        for a in range(min_frames, len(seq) - min_frames + 1):
            before = int(np.median(ss[a - min_frames:a]))
            after = int(np.median(ss[a:a + min_frames]))
            if after != before:
                rows.append({"t": ts[a], "label": lab, "size_from": before,
                             "size_to": after,
                             "kind": "nucleate" if (before == 1 and after == 2)
                             else ("grow" if after > before else "shrink")})
    df = pd.DataFrame(rows, columns=["t", "label", "size_from", "size_to",
                                     "kind"])
    if len(df):
        # keep only the first transition of each contiguous run per track
        df = df.sort_values(["label", "t"]).reset_index(drop=True)
        keep = np.ones(len(df), dtype=bool)
        for i in range(1, len(df)):
            if (df.label[i] == df.label[i - 1]
                    and df.t[i] - df.t[i - 1] < 1.0
                    and df.size_to[i] == df.size_to[i - 1]):
                keep[i] = False
        df = df[keep].sort_values("t").reset_index(drop=True)
    return df


def dimer_lifetimes_decoded(
        decoded: Dict[int, Tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    """Residence times in the size-2 state, from decoded size sequences.

    Only runs that both start and end inside the observation window are used,
    so the estimate is not truncated by censoring at the edges.
    """
    out = []
    for _lab, (ts, ss) in decoded.items():
        i = 0
        while i < len(ss):
            if ss[i] == 2:
                j = i
                while j + 1 < len(ss) and ss[j + 1] == 2:
                    j += 1
                if i > 0 and j < len(ss) - 1:
                    out.append(float(ts[j + 1] - ts[i]))
                i = j + 1
            else:
                i += 1
    return np.asarray(out)


def dimer_lifetimes(tracks: Tracks) -> np.ndarray:
    """Observed residence times of tracks in the size-2 state (seconds)."""
    by_label: Dict[int, List[Tuple[float, int]]] = {}
    for k in range(len(tracks.label)):
        s = sizes_from_amp(tracks.amp[k])
        for i, lab in enumerate(tracks.label[k]):
            by_label.setdefault(int(lab), []).append((float(tracks.t[k]),
                                                      int(s[i])))
    out = []
    for _lab, seq in by_label.items():
        seq.sort(key=lambda z: z[0])
        run_start = None
        for t, s in seq:
            if s == 2 and run_start is None:
                run_start = t
            elif s != 2 and run_start is not None:
                out.append(t - run_start)
                run_start = None
    return np.asarray(out)


# --------------------------------------------------------------------------- #
# assembly
# --------------------------------------------------------------------------- #
def _precompute(tracks: Tracks, decoded: Dict[int, Tuple[np.ndarray,
                                                          np.ndarray]],
                min_track_frames: int = 5) -> Dict[str, object]:
    """Precompute displacement pairs and dimer runs once, not per bin.

    Recomputing these inside every analysis bin is what makes a naive
    implementation quadratic in the record length; the descriptors are
    identical either way.
    """
    # map (label, frame) -> decoded size, and per-frame arrays
    frame_of_t = {float(tt): k for k, tt in enumerate(tracks.t)}
    size_lut: Dict[Tuple[int, int], int] = {}
    for lab, (ts, ss) in decoded.items():
        for tt, sv in zip(ts, ss):
            k = frame_of_t.get(float(tt))
            if k is not None:
                size_lut[(int(lab), k)] = int(sv)

    # Persistence filter.  A false positive is an independent noise excursion,
    # so it appears in one frame and is never linked forward; a real atom
    # persists.  Requiring a minimum track length therefore removes almost all
    # false positives, which is essential at low dose where they outnumber real
    # detections (see study S6).  The surviving/total detection ratio is kept
    # as an observable contamination diagnostic that needs no ground truth.
    track_len = {int(lab): len(ts) for lab, (ts, _ss) in decoded.items()}
    persistent = {lab for lab, n in track_len.items()
                  if n >= min_track_frames}
    n_det = sum(len(l) for l in tracks.label)
    n_persist = sum(int(np.sum([int(l) in persistent for l in lab]))
                    for lab in tracks.label) if n_det else 0
    persistent_frac = n_persist / n_det if n_det else np.nan

    # consecutive-frame displacement pairs, flagged monomer/monomer
    disp_t: List[float] = []
    disp_d2: List[float] = []
    disp_mono: List[bool] = []
    for k in range(len(tracks.label) - 1):
        lab0 = tracks.label[k]
        lab1 = tracks.label[k + 1]
        pos0 = {int(l): tracks.xy[k][i] for i, l in enumerate(lab0)}
        for i, l in enumerate(lab1):
            if int(l) not in persistent:
                continue
            p0 = pos0.get(int(l))
            if p0 is None:
                continue
            p1 = tracks.xy[k + 1][i]
            s0 = size_lut.get((int(l), k), 1)
            s1 = size_lut.get((int(l), k + 1), 1)
            disp_t.append(float(tracks.t[k + 1]))
            disp_d2.append(float(np.sum((p1 - p0) ** 2)))
            disp_mono.append(bool(s0 == 1 and s1 == 1))

    # dimer residence runs, keeping only runs bounded on both sides
    run_end: List[float] = []
    run_dur: List[float] = []
    for _lab, (ts, ss) in decoded.items():
        if int(_lab) not in persistent:
            continue
        i = 0
        while i < len(ss):
            if ss[i] == 2:
                j = i
                while j + 1 < len(ss) and ss[j + 1] == 2:
                    j += 1
                if i > 0 and j < len(ss) - 1:
                    run_end.append(float(ts[j + 1]))
                    run_dur.append(float(ts[j + 1] - ts[i]))
                i = j + 1
            else:
                i += 1

    order = np.argsort(disp_t) if disp_t else np.zeros(0, dtype=int)
    oe = np.argsort(run_end) if run_end else np.zeros(0, dtype=int)
    return {
        "size_lut": size_lut,
        "disp_t": np.asarray(disp_t)[order] if disp_t else np.zeros(0),
        "disp_d2": np.asarray(disp_d2)[order] if disp_t else np.zeros(0),
        "disp_mono": np.asarray(disp_mono)[order] if disp_t
        else np.zeros(0, bool),
        "run_end": np.asarray(run_end)[oe] if run_end else np.zeros(0),
        "run_dur": np.asarray(run_dur)[oe] if run_end else np.zeros(0),
        "persistent": persistent,
        "persistent_frac": persistent_frac,
    }


def build_descriptors(cfg: Config, tracks: Tracks,
                      frames: Sequence[ObservedFrame],
                      q_pt: Optional[np.ndarray],
                      schedule: Schedule,
                      dose_rate: Optional[float] = None,
                      bin_s: Optional[float] = None,
                      window_s: float = 10.0,
                      decode_method: str = "median",
                      min_frames: int = 3,
                      min_track_frames: int = 5,
                      switch_log_prior: float = -6.0) -> pd.DataFrame:
    """Assemble D(t) on the analysis grid.

    Census columns are per-bin averages; dynamic (rate-type) columns are
    estimated on a trailing window of ``window_s`` because a 1 s bin at 5 fps
    contains too few displacement pairs to estimate a rate.
    """
    from .instruments import single_atom_cnr
    bs = cfg.analysis.bin_s if bin_s is None else bin_s
    dr = cfg.imaging.dose_rate if dose_rate is None else dose_rate
    cnr1 = single_atom_cnr(cfg, dr)
    t_all = tracks.t
    dt_frame = float(np.mean(np.diff(t_all))) if len(t_all) > 1 else 1.0
    t0, t1 = float(t_all[0]), float(t_all[-1])
    edges = np.arange(t0, t1 + bs, bs)
    centres = edges[:-1] + bs / 2.0

    decoded = decode_tracks(tracks, cnr1, method=decode_method,
                            min_frames=min_frames,
                            switch_log_prior=switch_log_prior)
    pre = _precompute(tracks, decoded, min_track_frames=min_track_frames)
    persistent = pre["persistent"]
    trans = transitions_from_decoded(
        {lab: v for lab, v in decoded.items() if int(lab) in persistent})
    size_lut = pre["size_lut"]
    tr_t = trans.t.to_numpy(dtype=float) if len(trans) else np.zeros(0)
    tr_kind = trans.kind.to_numpy() if len(trans) else np.zeros(0, dtype=object)

    sig_loc = localisation_sigma_nm(cfg, dr)
    rows = []
    for c in centres:
        idx = np.flatnonzero((t_all >= c - bs / 2) & (t_all < c + bs / 2))
        if idx.size == 0:
            continue
        n1 = n2 = n3 = n4 = 0
        sz_all: List[int] = []
        for k in idx:
            s = np.array([size_lut.get((int(lab), int(k)), 1)
                          for lab in tracks.label[k]], dtype=int)
            n1 += int(np.sum(s == 1)); n2 += int(np.sum(s == 2))
            n3 += int(np.sum(s == 3)); n4 += int(np.sum(s >= 4))
            sz_all.extend(s.tolist())
        nf = idx.size
        n1, n2, n3, n4 = n1 / nf, n2 / nf, n3 / nf, n4 / nf
        sz = np.asarray(sz_all, dtype=float)
        n_ent = max(n1 + n2 + n3 + n4, 1e-9)
        mean_size = float(sz.mean()) if sz.size else np.nan
        cn = float(np.sum(sz * (sz - 1)) / max(np.sum(sz), 1e-9)) \
            if sz.size else np.nan

        # --- dynamic descriptors on a trailing window, from precomputes ---
        lo, hi = c - window_s, c
        a = int(np.searchsorted(pre["disp_t"], lo, "left"))
        b = int(np.searchsorted(pre["disp_t"], hi, "right"))
        if b > a:
            mono = pre["disp_mono"][a:b]
            d2 = pre["disp_d2"][a:b][mono]
            if d2.size:
                msd = float(np.mean(d2))
                khop = float(max(msd - 2.0 * sig_loc ** 2, 0.0)
                             / (cfg.support.lattice_nm ** 2 * dt_frame))
            else:
                msd, khop = np.nan, np.nan
        else:
            msd, khop = np.nan, np.nan

        ra = int(np.searchsorted(pre["run_end"], lo, "left"))
        rb = int(np.searchsorted(pre["run_end"], hi, "right"))
        tau_d = float(np.mean(pre["run_dur"][ra:rb])) if rb > ra else np.nan

        if tr_t.size:
            ta = int(np.searchsorted(tr_t, lo, "left"))
            tb = int(np.searchsorted(tr_t, hi, "right"))
            kinds = tr_kind[ta:tb]
            knuc = float(np.sum(kinds == "nucleate")) / window_s
            kgrow = float(np.sum(kinds == "grow")) / window_s
            kdis = float(np.sum(kinds == "shrink")) / window_s
        else:
            knuc = kgrow = kdis = 0.0

        cond = schedule.conditions(float(c))
        qi = float(np.nanmean(q_pt[idx])) if q_pt is not None else np.nan
        rows.append({
            "t": float(c), "N1": n1, "N2": n2, "N3": n3, "N4p": n4,
            "mean_size": mean_size, "cn_pt_pt": cn,
            "frac_dispersed": n1 / n_ent,
            "msd_1s": msd, "k_hop": khop, "k_nuc": knuc, "k_grow": kgrow,
            "k_diss": kdis, "tau_dimer": tau_d, "q_pt": qi,
            "T_c": cond["t_c"], "p_co": cond["p_co"], "p_o2": cond["p_o2"],
            "dose_rate": dr, "f_rep": cfg.reactor.f_rep,
            "persistent_frac": pre["persistent_frac"],
        })
    df = pd.DataFrame(rows)
    df.attrs["provenance"] = {c: PROVENANCE.get(c, "derived")
                              for c in df.columns}
    df.attrs["link_accuracy"] = tracks.link_accuracy
    df.attrs["persistent_frac"] = pre["persistent_frac"]
    df.attrs["transitions"] = trans
    df.attrs["decoded"] = decoded
    return df


CENSUS_COLS = ["N1", "N2", "N3", "N4p", "mean_size", "cn_pt_pt",
               "frac_dispersed"]
DYNAMIC_COLS = ["k_hop", "k_nuc", "k_grow", "k_diss", "tau_dimer", "msd_1s"]
REACTOR_COLS = ["T_c", "p_co", "p_o2"]
SPECTRO_COLS = ["q_pt"]


def check_provenance(df: pd.DataFrame, cols: Sequence[str],
                     allow: Sequence[str]) -> None:
    """Raise if a model is about to consume a stream it must not (rule 1/2)."""
    prov = df.attrs.get("provenance", {})
    bad = [c for c in cols if prov.get(c, "derived") not in allow]
    if bad:
        raise ValueError(
            f"columns {bad} have provenance "
            f"{[prov.get(c) for c in bad]}, which this estimator does not "
            f"accept (allowed: {list(allow)})")
