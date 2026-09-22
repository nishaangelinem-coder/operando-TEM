"""Projection-aware mobility estimation (gap G3).

A continuous image stream reports only the projected coordinates, so the
mean-squared displacement it yields is

.. math:: MSD_{2D} = MSD_{3D} - \\langle \\Delta z^2 \\rangle

and is therefore biased low.  The bias is not a constant: it depends on the
support's out-of-plane corrugation, which is exactly what the intermittent 3D
snapshots can measure even though they are far too sparse to track an atom.

The correction implemented here is deliberately minimal: estimate the
anisotropy ratio ``eta = <dz^2> / <dx^2 + dy^2>`` from the snapshot stream's
height statistics, then report ``MSD_3D = (1 + eta) * MSD_2D`` with the
uncertainty of ``eta`` propagated.  This keeps the 3D-derived quantity tagged
as such and never interpolates the snapshot stream onto the continuous grid.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from .config import Config, DEFAULT


def msd_3d_truth(pos_by_frame: Sequence[np.ndarray],
                 ids_by_frame: Sequence[np.ndarray],
                 sizes_by_frame: Sequence[np.ndarray],
                 lag_frames: int = 1, monomer_only: bool = True,
                 box_nm: Optional[float] = None) -> Dict[str, float]:
    """Ground-truth 3D and projected MSD from the simulator's own coordinates.

    Available only inside the digital twin; it is what the projection-aware
    estimator is validated against.

    ``box_nm`` applies the minimum-image convention to the in-plane
    displacement.  The patch simulator uses periodic boundaries, so an atom
    leaving one edge reappears at the other; without minimum-imaging those
    rare wrap events contribute displacements of order the box size and
    swamp the mean squared displacement entirely.  The observed-track
    estimator is immune because its linking gate rejects such jumps.
    """
    idx = [{int(i): k for k, i in enumerate(ids)} for ids in ids_by_frame]
    d3, d2, dz = [], [], []
    for k in range(len(pos_by_frame) - lag_frames):
        m0, m1 = idx[k], idx[k + lag_frames]
        for pid, i0 in m0.items():
            i1 = m1.get(pid)
            if i1 is None:
                continue
            if monomer_only and not (sizes_by_frame[k][i0] == 1
                                     and sizes_by_frame[k + lag_frames][i1] == 1):
                continue
            a = pos_by_frame[k][i0]
            b = pos_by_frame[k + lag_frames][i1]
            d = b - a
            if box_nm:
                d[:2] -= box_nm * np.round(d[:2] / box_nm)
            d3.append(float(np.sum(d ** 2)))
            d2.append(float(d[0] ** 2 + d[1] ** 2))
            dz.append(float(d[2] ** 2))
    if not d3:
        return {"msd_3d": np.nan, "msd_2d": np.nan, "eta": np.nan, "n": 0}
    m3, m2, mz = float(np.mean(d3)), float(np.mean(d2)), float(np.mean(dz))
    return {"msd_3d": m3, "msd_2d": m2, "msd_z": mz,
            "eta": mz / m2 if m2 > 0 else np.nan,
            "bias_frac": (m2 - m3) / m3 if m3 > 0 else np.nan,
            "n": len(d3)}


def eta_from_snapshots(cfg: Config, snapshots: Sequence[Dict[str, object]],
                       lattice_nm: Optional[float] = None) -> Dict[str, float]:
    """Estimate the out-of-plane anisotropy from intermittent 3D snapshots.

    The snapshots cannot follow an atom, but they do sample the support height
    field.  Fitting the local height gradient from the (x, y, z) point cloud
    gives <|grad z|^2>, and for a lattice hop of length ``a`` in a random
    in-plane direction,

        <dz^2> = 0.5 * <|grad z|^2> * a^2 ,   <dx^2 + dy^2> = a^2

    so ``eta = 0.5 * <|grad z|^2>``.  The depth noise of the reconstruction
    inflates the apparent gradient and is subtracted using the reported
    ``sigma_z``.
    """
    a = cfg.support.lattice_nm if lattice_nm is None else lattice_nm
    grads: List[float] = []
    for snap in snapshots:
        pos = np.asarray(snap["pos3d"], dtype=float)
        sig_z = float(snap.get("sigma_z_nm", 0.0))
        if len(pos) < 6:
            continue
        # local plane fit around each point using its k nearest neighbours
        for i in range(len(pos)):
            d = np.linalg.norm(pos[:, :2] - pos[i, :2], axis=1)
            order = np.argsort(d)[1:6]
            if len(order) < 4:
                continue
            p = pos[order]
            amat = np.column_stack([p[:, 0] - pos[i, 0],
                                    p[:, 1] - pos[i, 1],
                                    np.ones(len(p))])
            try:
                sol, *_ = np.linalg.lstsq(amat, p[:, 2] - pos[i, 2], rcond=None)
            except np.linalg.LinAlgError:
                continue
            g2 = float(sol[0] ** 2 + sol[1] ** 2)
            # subtract the variance the depth noise alone would produce
            spread = float(np.mean(d[order] ** 2))
            g2_noise = 2.0 * sig_z ** 2 / max(spread, 1e-9)
            grads.append(max(g2 - g2_noise, 0.0))
    if not grads:
        return {"eta": np.nan, "eta_se": np.nan, "n": 0}
    g = np.asarray(grads)
    eta = 0.5 * float(np.mean(g))
    se = 0.5 * float(np.std(g, ddof=1) / np.sqrt(len(g)))
    return {"eta": eta, "eta_se": se, "n": int(len(g)),
            "mean_grad2": float(np.mean(g))}


def eta_analytic(cfg: Config) -> float:
    """Analytic anisotropy for the configured sinusoidal corrugation.

    For ``z = A sin(qx) sin(qy)`` with ``q = 2 pi / lambda``, averaging
    ``|grad z|^2 = A^2 q^2 [cos^2(qx) sin^2(qy) + sin^2(qx) cos^2(qy)]`` over
    the surface gives ``A^2 q^2 / 2``, hence ``eta = A^2 q^2 / 4``.
    """
    s = cfg.support
    q = 2.0 * np.pi / s.corrug_wavelength_nm
    return float(s.corrug_amp_nm ** 2 * q ** 2 / 4.0)


def correct_msd(msd_2d: float, eta: float, eta_se: float = 0.0
                ) -> Dict[str, float]:
    """Projection-corrected 3D MSD with propagated uncertainty."""
    val = (1.0 + eta) * msd_2d
    return {"msd_3d_est": float(val),
            "msd_3d_se": float(msd_2d * eta_se),
            "eta_used": float(eta)}
