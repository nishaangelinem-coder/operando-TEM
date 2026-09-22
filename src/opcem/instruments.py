"""Virtual instruments: forward models of what each detector would record.

Every function here maps ground truth (from :mod:`opcem.truth`) to a realistic
detector output.  Two fidelity levels are provided and are never mixed
silently:

*   :func:`render_adf_frame` renders actual pixels, including Poisson shot
    noise, probe blur, scan jitter, drift and window background.  It is used to
    validate the detection front end (:mod:`opcem.vision`) and to produce the
    example micrograph figure.
*   :func:`observe_positions` is a fast analytic surrogate that emits noisy
    projected coordinates directly, with the dose-dependent localisation
    precision, missed detections and false positives that the pixel-level
    pipeline actually achieves.  The bulk studies use this level, and the
    equivalence of the two levels is checked in study S0.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import ndimage

from .config import Config, DEFAULT
from .truth import PatchResult, SITE_TYPES


# --------------------------------------------------------------------------- #
# ADF-STEM
# --------------------------------------------------------------------------- #
def dose_per_frame(cfg: Config, dose_rate: Optional[float] = None,
                   fps: Optional[float] = None) -> float:
    """Electron dose accumulated in one frame (e-/A^2)."""
    dr = cfg.imaging.dose_rate if dose_rate is None else dose_rate
    f = cfg.imaging.fps if fps is None else fps
    return dr / f


def localisation_sigma_nm(cfg: Config, dose_rate: Optional[float] = None,
                          fps: Optional[float] = None) -> float:
    """Single-atom localisation precision, scaling as dose^(-1/2)."""
    im = cfg.imaging
    d = dose_per_frame(cfg, dose_rate, fps)
    d_ref = im.ref_dose_rate / im.fps
    return im.loc_sigma_nm_at_ref * np.sqrt(d_ref / max(d, 1e-9))


def single_atom_cnr(cfg: Config, dose_rate: Optional[float] = None,
                    fps: Optional[float] = None) -> float:
    """Contrast-to-noise ratio of one Pt atom in one frame."""
    im = cfg.imaging
    d = dose_per_frame(cfg, dose_rate, fps)
    probe_area = np.pi * (2.0 * im.probe_sigma_nm * 10.0) ** 2   # A^2
    n_bg = im.bg_scatter_frac * d * probe_area
    s_atom = im.sigma_adf_pt_A2 * d
    return s_atom / np.sqrt(max(n_bg + s_atom, 1e-12))


def render_adf_frame(cfg: Config, pos3d: np.ndarray, sizes: np.ndarray,
                     rng: np.random.Generator,
                     dose_rate: Optional[float] = None,
                     fps: Optional[float] = None,
                     drift_nm: Tuple[float, float] = (0.0, 0.0),
                     crop_origin_nm: Tuple[float, float] = (0.0, 0.0),
                     field_nm: Optional[float] = None) -> np.ndarray:
    """Render one ADF-STEM frame in detected electrons per pixel.

    Rendering is done on a CROP of the imaged field (default 8 nm) because
    resolving single atoms requires ~0.2 A sampling, and sampling the whole
    40 nm field that finely would be a 2000 x 2000 array per frame.  The
    analytic observation model is used for the full field instead.

    The z coordinate is *discarded* by the projection, which is the physical
    origin of gap G3: a continuous image stream cannot report depth.
    """
    im = cfg.imaging
    L = im.render_field_nm if field_nm is None else field_nm
    npix = int(round(L / im.pixel_nm))
    d = dose_per_frame(cfg, dose_rate, fps)
    pix_area_a2 = (im.pixel_nm * 10.0) ** 2

    # support / window background, with slowly varying thickness
    bg = im.bg_scatter_frac * d * pix_area_a2
    field = np.full((npix, npix), bg, dtype=float)
    w = rng.normal(size=(npix, npix))
    w = ndimage.gaussian_filter(w, sigma=npix / 12.0)
    w /= max(np.std(w), 1e-12)
    field *= (1.0 + im.window_bg_level * w)

    sig_total = im.sigma_adf_pt_A2 * d          # electrons per Pt atom
    sig_pix = im.probe_sigma_nm / im.pixel_nm
    jit = rng.normal(0.0, im.scan_jitter_nm, size=npix)   # one per scan line
    atoms = np.zeros((npix, npix), dtype=float)
    x0, y0 = crop_origin_nm
    for (x, y, _z), n in zip(pos3d, sizes):
        xi = (x + drift_nm[0] - x0) / im.pixel_nm
        yi = (y + drift_nm[1] - y0) / im.pixel_nm
        if not (0.0 <= xi < npix - 1 and 0.0 <= yi < npix - 1):
            continue
        xi += jit[int(yi)] / im.pixel_nm
        ix, iy = int(np.floor(xi)), int(np.floor(yi))
        if not (0 <= ix < npix - 1 and 0 <= iy < npix - 1):
            continue
        fx, fy = xi - ix, yi - iy
        for dx in (0, 1):
            for dy in (0, 1):
                wgt = (fx if dx else 1 - fx) * (fy if dy else 1 - fy)
                atoms[iy + dy, ix + dx] += wgt * sig_total * n
    # gaussian_filter conserves the sum, so each atom keeps its intended
    # integrated signal of sig_total * n electrons
    atoms = ndimage.gaussian_filter(atoms, sigma=sig_pix, mode="nearest")
    return rng.poisson(np.clip(field + atoms, 0.0, None)).astype(float)


def crop_truth(pos3d: np.ndarray, sizes: np.ndarray,
               crop_origin_nm: Tuple[float, float], field_nm: float,
               margin_nm: float = 0.2):
    """Ground-truth atoms fully inside a render crop (for scoring detection)."""
    x0, y0 = crop_origin_nm
    m = ((pos3d[:, 0] > x0 + margin_nm)
         & (pos3d[:, 0] < x0 + field_nm - margin_nm)
         & (pos3d[:, 1] > y0 + margin_nm)
         & (pos3d[:, 1] < y0 + field_nm - margin_nm))
    return pos3d[m][:, :2] - np.array([x0, y0]), np.asarray(sizes)[m]


# --------------------------------------------------------------------------- #
# fast analytic observation model
# --------------------------------------------------------------------------- #
@dataclass
class ObservedFrame:
    t: float
    xy: np.ndarray            # (n_det, 2) observed projected positions (nm)
    amp: np.ndarray           # (n_det,) apparent atom count (from ADF intensity)
    true_id: np.ndarray       # ground-truth entity id, -1 for a false positive
    q_pt: float               # EELS oxidation-state proxy for the field
    dose_rate: float


def observe_positions(cfg: Config, res: PatchResult, rng: np.random.Generator,
                      dose_rate: Optional[float] = None,
                      include_drift: bool = True) -> List[ObservedFrame]:
    """Fast surrogate for the pixel pipeline: noisy projected coordinates.

    Applies, in order: projection (z discarded), stage drift, scan jitter,
    dose-dependent localisation noise, size-dependent detection probability,
    and Poisson false positives over the field.
    """
    im = cfg.imaging
    dr = im.dose_rate if dose_rate is None else dose_rate
    sig = localisation_sigma_nm(cfg, dr)
    cnr1 = single_atom_cnr(cfg, dr)
    L = cfg.support.field_nm
    out: List[ObservedFrame] = []
    for k, t in enumerate(res.t):
        pos, sizes, ids = res.pos3d[k], res.sizes[k], res.ids[k]
        xy = pos[:, :2].copy()
        if include_drift:
            xy[:, 0] += im.drift_nm_per_s * t
            xy[:, 1] += 0.4 * im.drift_nm_per_s * t
        xy += rng.normal(0.0, im.scan_jitter_nm, size=xy.shape)
        # localisation precision improves as sqrt(signal)
        xy += rng.normal(0.0, 1.0, size=xy.shape) * (sig / np.sqrt(
            np.maximum(sizes, 1))[:, None])
        # Atom counting from integrated ADF intensity.  The shot-noise limit is
        # ABSOLUTE in atom units and set by the single-atom CNR, so a large
        # cluster is counted proportionally more precisely than a monomer.
        amp_noisy = np.asarray(sizes, dtype=float) + rng.normal(
            0.0, 1.0 / max(cnr1, 1e-6), size=len(sizes))
        # Detection thresholds the NOISY amplitude, not the true one.  This
        # reproduces two things measured on rendered pixels in study S0: the
        # ~49% monomer recall at the reference dose, and the +0.23 atom
        # Eddington bias of the surviving amplitudes (only upward noise
        # excursions of a marginal atom clear the threshold).
        keep = amp_noisy * cnr1 > im.detect_threshold_cnr
        xy_k, amp_k = xy[keep], amp_noisy[keep]
        ids_k = np.asarray(ids)[keep]
        fp_density = im.false_positive_per_nm2 * (
            im.ref_dose_rate / max(dr, 1e-9)) ** im.fp_dose_exponent
        n_fp = rng.poisson(fp_density * L * L)
        if n_fp:
            xy_k = np.vstack([xy_k, rng.uniform(0, L, size=(n_fp, 2))])
            amp_k = np.concatenate([amp_k, rng.uniform(0.6, 1.4, size=n_fp)])
            ids_k = np.concatenate([ids_k, -np.ones(n_fp, dtype=np.int64)])
        out.append(ObservedFrame(t=float(t), xy=xy_k, amp=amp_k,
                                 true_id=ids_k, q_pt=np.nan, dose_rate=dr))
    return out


# --------------------------------------------------------------------------- #
# EELS
# --------------------------------------------------------------------------- #
def eels_oxidation_proxy(cfg: Config, res: PatchResult,
                         rng: np.random.Generator,
                         dose_rate: Optional[float] = None) -> np.ndarray:
    """Field-averaged Pt oxidation-state proxy q_Pt from the Pt M4,5 white lines.

    Ground truth: a Pt atom is more oxidised the smaller its entity and the
    deeper its support site (vacancy-anchored Pt is the most oxidised).  Only
    the field average is returned, because at single-atom concentration the
    per-atom core-loss signal is below the noise floor at tolerable dose --
    which is itself one of the findings this platform has to respect.
    """
    dr = cfg.imaging.dose_rate if dose_rate is None else dose_rate
    snr = cfg.imaging.eels_snr_at_ref * np.sqrt(
        dr / cfg.imaging.ref_dose_rate)
    q_site = {"terrace": 1.6, "vacancy": 2.4, "step": 1.9}
    out = np.empty(len(res.t))
    for k in range(len(res.t)):
        sizes = np.asarray(res.sizes[k], dtype=float)
        sites = np.asarray(res.sites[k])
        if sizes.size == 0:
            out[k] = np.nan
            continue
        base = np.array([q_site[SITE_TYPES[s]] for s in sites])
        q_atom = base * np.power(sizes, -0.35)       # larger => more metallic
        w = sizes / sizes.sum()
        out[k] = float(np.dot(w, q_atom))
    return out + rng.normal(0.0, np.mean(out) / max(snr, 1e-6), size=out.shape)


# --------------------------------------------------------------------------- #
# intermittent 3D snapshots
# --------------------------------------------------------------------------- #
def snapshot3d(cfg: Config, res: PatchResult, rng: np.random.Generator,
               interval_s: Optional[float] = None) -> List[Dict[str, object]]:
    """Multi-tilt ptychographic 3D snapshots at a coarse interval.

    Lateral precision is the imaging precision; depth precision is set by the
    measured depth resolution (FWHM) of the reconstruction, converted to a
    standard deviation.  This is the ONLY stream that carries z, and it is
    tagged ``snapshot3d`` so that no continuous-time estimator can consume it
    by accident (design rule 2).
    """
    im = cfg.imaging
    iv = im.snapshot3d_interval_s if interval_s is None else interval_s
    sig_xy = localisation_sigma_nm(cfg, im.dose_rate)
    sig_z = im.depth_resolution_nm / 2.3548
    snaps = []
    next_t = 0.0
    for k, t in enumerate(res.t):
        if t + 1e-9 < next_t:
            continue
        next_t = t + iv
        pos = res.pos3d[k].copy()
        pos[:, :2] += rng.normal(0.0, sig_xy, size=pos[:, :2].shape)
        pos[:, 2] += rng.normal(0.0, sig_z, size=pos[:, 2].shape)
        snaps.append({"t": float(t), "pos3d": pos,
                      "sizes": np.asarray(res.sizes[k]).copy(),
                      "sigma_xy_nm": sig_xy, "sigma_z_nm": sig_z,
                      "provenance": "snapshot3d"})
    return snaps


# --------------------------------------------------------------------------- #
# quadrupole mass spectrometer
# --------------------------------------------------------------------------- #
@dataclass
class MSTrace:
    t: np.ndarray
    ion_current: Dict[str, np.ndarray]
    rate_hz: float
    provenance: str = "ms"


def _pink_noise(n: int, rng: np.random.Generator) -> np.ndarray:
    """Unit-variance 1/f noise via spectral shaping."""
    f = np.fft.rfftfreq(n, d=1.0)
    amp = np.ones_like(f)
    amp[1:] = 1.0 / np.sqrt(f[1:])
    amp[0] = 0.0
    ph = rng.normal(size=amp.shape) + 1j * rng.normal(size=amp.shape)
    x = np.fft.irfft(amp * ph, n=n)
    sd = np.std(x)
    return x / sd if sd > 0 else x


def virtual_qms(cfg: Config, t_grid: np.ndarray, rate_molec_s: np.ndarray,
                kernel: np.ndarray, rng: np.random.Generator,
                inlet_flux: Optional[Dict[str, float]] = None,
                ar_tracer: Optional[np.ndarray] = None) -> MSTrace:
    """Forward model of the online QMS.

    ``rate_molec_s`` is the *true* CO2 formation rate integrated over the whole
    reactive zone, in the same arbitrary flux units as ``inlet_flux``.  It is
    convolved with the measured transport kernel, scaled by the per-species
    calibration factor, and corrupted by 1/f baseline drift, Gaussian noise,
    a constant background and N2 interference on the m/z 28 (CO) channel.
    """
    ms = cfg.ms
    n = len(t_grid)
    if inlet_flux is None:
        inlet_flux = {"CO": 1.0, "O2": 2.0, "Ar": 0.0}

    def conv(x: np.ndarray) -> np.ndarray:
        y = np.convolve(x, kernel, mode="full")[:n]
        return y

    co2 = conv(rate_molec_s)
    co = conv(np.full(n, inlet_flux["CO"]) - rate_molec_s)
    o2 = conv(np.full(n, inlet_flux["O2"]) - 0.5 * rate_molec_s)
    ar = conv(ar_tracer if ar_tracer is not None else
              np.full(n, inlet_flux.get("Ar", 0.0)))

    out: Dict[str, np.ndarray] = {}
    for name, sig in (("CO2", co2), ("CO", co), ("O2", o2), ("Ar", ar)):
        val = ms.k_cal[name] * sig + ms.background[name]
        val = val * (1.0 + ms.drift_rel * _pink_noise(n, rng))
        scale = np.mean(np.abs(val)) if np.mean(np.abs(val)) > 0 else 1.0
        val = val + rng.normal(0.0, ms.noise_rel * scale, size=n)
        if name == "CO":
            val = val + ms.n2_interference
        out[name] = val
    return MSTrace(t=t_grid, ion_current=out, rate_hz=1.0 / float(
        np.mean(np.diff(t_grid))) if n > 1 else ms.rate_hz)
