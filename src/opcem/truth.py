"""Ground-truth generative model (the 'digital twin').

Two coupled levels, as required by docs/04_experimental_setup.md section 5.1:

*   :class:`PatchSimulator` -- a spatially resolved, discrete-time kinetic Monte
    Carlo model of Pt entities inside ONE imaged field of view.  This is the
    only level that produces atom trajectories, so it is the only level the
    microscope can see.  It is illuminated by the electron beam.

*   :func:`chip_mean_field` -- a deterministic population-balance model of the
    whole reactive zone, position-binned along the flow axis so that the
    temperature and partial-pressure gradients of the real chip are present.
    This is what the mass spectrometer actually integrates over.  It is NOT
    illuminated (dose rate zero outside the imaged field), which is the source
    of the systematic patch/chip discrepancy analysed in study S7.

Nothing in this module is a measurement.  It is a forward model whose outputs
are consumed by :mod:`opcem.instruments`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .config import Config, KB_EV, DEFAULT

SITE_TYPES = ("terrace", "vacancy", "step")


# --------------------------------------------------------------------------- #
# thermodynamic / kinetic helpers
# --------------------------------------------------------------------------- #
def langmuir_coverages(cfg: Config, temp_k: float, p_co: float,
                       p_o2: float) -> Tuple[float, float]:
    """Competitive Langmuir coverages for CO (molecular) and O (dissociative).

    Adsorption constants are corrected from the reference temperature with a
    van't Hoff factor using the configured adsorption enthalpies.
    """
    a = cfg.activity
    t_ref = a.t_ref_k
    def corr(k_ref: float, dh: float) -> float:
        return k_ref * np.exp(-(dh / KB_EV) * (1.0 / temp_k - 1.0 / t_ref))
    k_co = corr(a.k_ads_co, a.dh_ads_co)
    k_o2 = corr(a.k_ads_o2, a.dh_ads_o2)
    s_o = np.sqrt(max(k_o2 * p_o2, 0.0))
    den = 1.0 + k_co * p_co + s_o
    return (k_co * p_co) / den, s_o / den


def tof_entity(cfg: Config, size: int, site: str, temp_k: float,
               p_co: float, p_o2: float) -> float:
    """Ground-truth turnover frequency (molecules/s) of one Pt_n entity."""
    a = cfg.activity
    th_co, th_o = langmuir_coverages(cfg, temp_k, p_co, p_o2)
    if size <= 4:
        pref, e_act = a.a_n[size], a.e_act[size]
    else:
        pref, e_act = a.a_n_large_per_atom * size, a.e_act_large
    return (pref * th_co * th_o * np.exp(-e_act / (KB_EV * temp_k))
            * a.site_mult[site])


def k_hop(cfg: Config, site: str, temp_k: float, dose_rate: float
          ) -> Tuple[float, float]:
    """Return (thermal, beam) contributions to the monomer hop rate (1/s)."""
    k = cfg.kinetics
    barrier = k.e_hop_terrace + cfg.support.e_trap[site]
    k_th = k.nu0_hop * np.exp(-barrier / (KB_EV * temp_k))
    k_bm = k.beta_beam * dose_rate
    return k_th, k_bm


def k_diss(cfg: Config, size: int, temp_k: float, p_o2: float) -> float:
    """Dissociation rate Pt_n -> Pt_(n-1) + Pt_1 (1/s), oxygen-accelerated."""
    if size < 2:
        return 0.0
    k = cfg.kinetics
    e_b = k.e_bind.get(size, k.e_bind_large)
    ox = 1.0 + k.gamma_ox * np.sqrt(max(p_o2, 0.0))
    return k.nu0_diss * np.exp(-e_b / (KB_EV * temp_k)) * ox


def corrugation_z(cfg: Config, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Support surface height, which makes lateral diffusion 3-dimensional."""
    s = cfg.support
    w = 2.0 * np.pi / s.corrug_wavelength_nm
    return s.corrug_amp_nm * np.sin(w * x) * np.sin(w * y)


# --------------------------------------------------------------------------- #
# actuator schedule
# --------------------------------------------------------------------------- #
@dataclass
class Actuation:
    """One actuator action on the reactor."""
    t_start: float
    duration: float
    kind: str = "o2_pulse"
    x_o2: Optional[float] = None
    delta_t_c: float = 0.0


class Schedule:
    """Time-resolved reactor conditions produced by a list of actuations."""

    def __init__(self, cfg: Config, actuations: Sequence[Actuation] = ()):
        self.cfg = cfg
        self.actuations = list(actuations)

    def conditions(self, t: float) -> Dict[str, float]:
        r = self.cfg.reactor
        x_o2, x_co, dt_c = r.x_o2, r.x_co, 0.0
        for a in self.actuations:
            if a.t_start <= t < a.t_start + a.duration:
                if a.kind == "o2_pulse":
                    x_o2 = a.x_o2 if a.x_o2 is not None else \
                        self.cfg.control.o2_pulse_x_o2
                elif a.kind == "co_lean":
                    x_co = r.x_co * 0.25
                elif a.kind == "t_step":
                    dt_c = a.delta_t_c
        p_tot_bar = r.p_total_mbar / 1000.0
        return {"x_co": x_co, "x_o2": x_o2,
                "p_co": x_co * p_tot_bar, "p_o2": x_o2 * p_tot_bar,
                "t_c": r.t_setpoint_c + dt_c}

    def is_active(self, t: float) -> bool:
        return any(a.t_start <= t < a.t_start + a.duration
                   for a in self.actuations)


# --------------------------------------------------------------------------- #
# patch-level kMC
# --------------------------------------------------------------------------- #
@dataclass
class PatchResult:
    """Ground truth for one imaged field of view."""
    t: np.ndarray                       # kMC macro-step times (s)
    true_rate: np.ndarray               # molecules/s from the imaged atoms
    pop: np.ndarray                     # (T,4) counts of N1,N2,N3,N4+
    pos3d: List[np.ndarray]             # per frame (n_entities,3) nm
    sizes: List[np.ndarray]             # per frame entity sizes
    sites: List[np.ndarray]             # per frame site-type indices
    ids: List[np.ndarray]               # per frame entity ids (for tracking truth)
    events: "np.ndarray"                # structured array of ground-truth events
    hop_counts: Dict[str, int]          # thermal vs beam-attributed hops
    monomer_tracks: Dict[int, List[Tuple[float, float, float, float]]]
    conditions: np.ndarray              # (T,3) t_c, p_co, p_o2
    cfg: Config = field(repr=False, default_factory=lambda: DEFAULT)

    @property
    def n_frames(self) -> int:
        return len(self.t)


EVENT_DTYPE = np.dtype([("t", "f8"), ("kind", "U12"), ("size_from", "i4"),
                        ("size_to", "i4"), ("beam", "?"), ("ids", "i8")])


class PatchSimulator:
    """Discrete-time kMC of Pt entities in one imaged field of view.

    The time step is chosen so that the largest per-entity transition
    probability stays small (checked at construction); transitions are then
    drawn as independent Bernoulli trials, which is the standard tau-leaping
    approximation and is far cheaper than exact Gillespie for the ~10^5 hops
    a single run contains.
    """

    def __init__(self, cfg: Config = DEFAULT, seed: int = 0,
                 dose_rate: Optional[float] = None,
                 schedule: Optional[Schedule] = None,
                 dt: float = 0.02, n_atoms: Optional[int] = None):
        self.cfg = cfg
        self.rng = np.random.default_rng(seed)
        self.dose_rate = cfg.imaging.dose_rate if dose_rate is None else dose_rate
        self.schedule = schedule or Schedule(cfg)
        self.dt = dt
        self.n_atoms = cfg.n_pt_atoms if n_atoms is None else n_atoms
        self._check_dt()

    # -- setup -------------------------------------------------------------- #
    def _check_dt(self) -> None:
        t_k = self.cfg.reactor.t_setpoint_c + 273.15
        k_th, k_bm = k_hop(self.cfg, "terrace", t_k, self.dose_rate)
        p = (k_th + k_bm) * self.dt
        if p > 0.2:
            raise ValueError(
                f"dt={self.dt}s gives hop probability {p:.3f} per step; "
                "reduce dt (tau-leaping validity)")

    def _draw_sites(self, n: int) -> np.ndarray:
        s = self.cfg.support
        p = np.array([1.0 - s.frac_vacancy - s.frac_step,
                      s.frac_vacancy, s.frac_step])
        return self.rng.choice(3, size=n, p=p / p.sum())

    # -- main loop ---------------------------------------------------------- #
    def run(self, duration: Optional[float] = None,
            record_fps: Optional[float] = None,
            controller=None,
            control_every_s: float = 1.0) -> PatchResult:
        cfg = self.cfg
        dur = cfg.duration_s if duration is None else duration
        fps = cfg.imaging.fps if record_fps is None else record_fps
        L = cfg.support.field_nm
        a_latt = cfg.support.lattice_nm

        # entity state arrays (index = entity slot)
        n_ent = self.n_atoms
        size = np.ones(n_ent, dtype=np.int32)
        xy = self.rng.uniform(0.0, L, size=(n_ent, 2))
        site = self._draw_sites(n_ent)
        alive = np.ones(n_ent, dtype=bool)
        eid = np.arange(n_ent, dtype=np.int64)
        next_id = n_ent

        rec_every = max(1, int(round(1.0 / (fps * self.dt))))
        # The controller is polled at the analysis-bin rate, not every kMC
        # step: Layer 3 has one analysis bin to decide (docs/04 section 4), so
        # polling faster than that would not be implementable on real hardware.
        ctrl_every = max(1, int(round(control_every_s / self.dt)))
        n_steps = int(round(dur / self.dt))

        t_rec, rate_rec, pop_rec, cond_rec = [], [], [], []
        pos_rec, size_rec, site_rec, id_rec = [], [], [], []
        events: List[tuple] = []
        hop_counts = {"thermal": 0, "beam": 0}
        tracks: Dict[int, List[Tuple[float, float, float, float]]] = {}

        cap_r2 = cfg.kinetics.capture_radius_nm ** 2
        for step in range(n_steps + 1):
            t = step * self.dt
            if controller is not None and step % ctrl_every == 0:
                controller.observe(t, size[alive], site[alive], self.schedule)
            cond = self.schedule.conditions(t)
            t_k = cond["t_c"] + 273.15
            p_co, p_o2 = cond["p_co"], cond["p_o2"]

            idx = np.flatnonzero(alive)

            # ---- turnover (instantaneous, no state change) ----
            rate = 0.0
            for i in idx:
                rate += tof_entity(cfg, int(size[i]), SITE_TYPES[site[i]],
                                   t_k, p_co, p_o2)

            if step % rec_every == 0:
                z = corrugation_z(cfg, xy[idx, 0], xy[idx, 1])
                t_rec.append(t)
                rate_rec.append(rate)
                s_i = size[idx]
                pop_rec.append([int(np.sum(s_i == 1)), int(np.sum(s_i == 2)),
                                int(np.sum(s_i == 3)), int(np.sum(s_i >= 4))])
                cond_rec.append([cond["t_c"], p_co, p_o2])
                pos_rec.append(np.column_stack([xy[idx], z]))
                size_rec.append(s_i.copy())
                site_rec.append(site[idx].copy())
                id_rec.append(eid[idx].copy())
                for k, i in enumerate(idx):
                    if size[i] == 1:
                        tracks.setdefault(int(eid[i]), []).append(
                            (t, xy[i, 0], xy[i, 1], float(z[k])))

            if step == n_steps:
                break

            # ---- dissociation of clusters ----
            clus = idx[size[idx] >= 2]
            if clus.size:
                kd = np.array([k_diss(cfg, int(size[i]), t_k, p_o2)
                               for i in clus])
                fire = self.rng.random(clus.size) < kd * self.dt
                for i in clus[fire]:
                    n_from = int(size[i])
                    size[i] = n_from - 1
                    # emit a free monomer next to the parent
                    ang = self.rng.uniform(0, 2 * np.pi)
                    nx = (xy[i, 0] + 1.3 * a_latt * np.cos(ang)) % L
                    ny = (xy[i, 1] + 1.3 * a_latt * np.sin(ang)) % L
                    slot = self._alloc(alive)
                    if slot is None:
                        continue
                    alive[slot] = True
                    size[slot] = 1
                    xy[slot] = (nx, ny)
                    site[slot] = self._draw_sites(1)[0]
                    eid[slot] = next_id
                    next_id += 1
                    events.append((t, "dissoc", n_from, n_from - 1, False,
                                   int(eid[i])))
                    if size[i] == 1:
                        site[i] = self._draw_sites(1)[0]

            # ---- monomer hops ----
            idx = np.flatnonzero(alive)
            mon = idx[size[idx] == 1]
            if mon.size:
                kth = np.empty(mon.size)
                for j, i in enumerate(mon):
                    kth[j] = k_hop(cfg, SITE_TYPES[site[i]], t_k, 0.0)[0]
                kbm = cfg.kinetics.beta_beam * self.dose_rate
                p_tot = (kth + kbm) * self.dt
                fire = self.rng.random(mon.size) < p_tot
                hop = mon[fire]
                if hop.size:
                    # attribute each hop to thermal or beam origin
                    frac_beam = kbm / (kth[fire] + kbm)
                    is_beam = self.rng.random(hop.size) < frac_beam
                    hop_counts["beam"] += int(is_beam.sum())
                    hop_counts["thermal"] += int((~is_beam).sum())
                    d = self.rng.integers(0, 4, size=hop.size)
                    dx = np.where(d == 0, a_latt, np.where(d == 1, -a_latt, 0.0))
                    dy = np.where(d == 2, a_latt, np.where(d == 3, -a_latt, 0.0))
                    xy[hop, 0] = (xy[hop, 0] + dx) % L
                    xy[hop, 1] = (xy[hop, 1] + dy) % L
                    site[hop] = self._draw_sites(hop.size)
                    # ---- encounter / merge, checked only for atoms that moved
                    for k, i in enumerate(hop):
                        if not alive[i] or size[i] != 1:
                            continue
                        others = np.flatnonzero(alive)
                        others = others[others != i]
                        if others.size == 0:
                            continue
                        d2 = ((xy[others, 0] - xy[i, 0]) ** 2
                              + (xy[others, 1] - xy[i, 1]) ** 2)
                        cap = cap_r2 * np.maximum(size[others], 1) ** (2.0 / 3.0)
                        hit = others[d2 < cap]
                        if hit.size == 0:
                            continue
                        j = hit[np.argmin(d2[np.isin(others, hit)])]
                        n_from = int(size[j])
                        size[j] = n_from + 1
                        kind = "nucleate" if n_from == 1 else "grow"
                        events.append((t, kind, n_from, n_from + 1,
                                       bool(is_beam[k]), int(eid[j])))
                        alive[i] = False

        ev = np.array(events, dtype=EVENT_DTYPE) if events else \
            np.empty(0, dtype=EVENT_DTYPE)
        return PatchResult(
            t=np.asarray(t_rec), true_rate=np.asarray(rate_rec),
            pop=np.asarray(pop_rec, dtype=float), pos3d=pos_rec,
            sizes=size_rec, sites=site_rec, ids=id_rec, events=ev,
            hop_counts=hop_counts, monomer_tracks=tracks,
            conditions=np.asarray(cond_rec), cfg=cfg)

    @staticmethod
    def _alloc(alive: np.ndarray) -> Optional[int]:
        free = np.flatnonzero(~alive)
        return int(free[0]) if free.size else None


# --------------------------------------------------------------------------- #
# chip-level mean field
# --------------------------------------------------------------------------- #
def site_occupancy(cfg: Config, temp_k: float, dose_rate: float
                   ) -> Tuple[np.ndarray, float, float]:
    """Time-averaged site-type occupancy of a diffusing monomer.

    In the patch kMC a monomer re-draws its site type after every hop, so its
    steady-state occupancy of site type ``s`` is proportional to
    ``p_s / k_s``: a deep site is entered as often as its abundance dictates
    but left more slowly.  Returns (occupancy, effective hop rate, effective
    activity multiplier), which is what the coarse-grained chip model needs.
    """
    s_cfg = cfg.support
    p = np.array([1.0 - s_cfg.frac_vacancy - s_cfg.frac_step,
                  s_cfg.frac_vacancy, s_cfg.frac_step])
    p = p / p.sum()
    k = np.array([sum(k_hop(cfg, st, temp_k, dose_rate)) for st in SITE_TYPES])
    resid = p / k
    occ = resid / resid.sum()
    k_eff = 1.0 / resid.sum()
    mult = float(np.dot(occ, [cfg.activity.site_mult[st] for st in SITE_TYPES]))
    return occ, float(k_eff), mult


def chip_mean_field(cfg: Config, t_grid: np.ndarray, schedule: Schedule,
                    k_enc: float, n_max: int = 14, n_bins: int = 12,
                    dose_rate: float = 0.0,
                    n_atoms: Optional[int] = None) -> Dict[str, np.ndarray]:
    """Size-resolved population balance for the whole reactive zone.

    A Smoluchowski-type balance restricted to monomer addition and single-atom
    loss, using exactly the dissociation rate constants of the patch kMC:

    .. math::
        \\dot n_1 = -2 a_{11} n_1^2 - \\sum_{k\\ge2} a_{1k} n_1 n_k
                    + \\sum_{k\\ge2} k_d(k)\\, n_k
        \\dot n_k = a_{1,k-1} n_1 n_{k-1} - a_{1k} n_1 n_k
                    - k_d(k) n_k + k_d(k{+}1) n_{k+1}

    with :math:`a_{1k} = k_{enc}\\,(k_{hop}/k_{hop}^{ref})\\,k^{2/3}`, matching
    the capture-radius scaling of the kMC.  The scheme conserves
    :math:`\\sum_k k\\,n_k` exactly in the continuum limit; the residual drift
    is reported by the caller.

    ``k_enc`` is the single free coefficient and is fitted to an ensemble of
    patch kMC runs by :func:`calibrate_k_enc`.

    The zone is split into ``n_bins`` slices along the flow axis, each with its
    own temperature and product-enrichment factor, so the chip carries the
    gradients that make one imaged field unrepresentative (gap G5).  The chip
    is not illuminated: ``dose_rate`` defaults to zero.

    Because the physical chip holds O(10^13) Pt atoms, the population is
    treated as deterministic; the chip rate therefore responds only to global
    drivers, never to single-atom stochasticity.  That is the central
    quantitative obstacle analysed in study S7.
    """
    r = cfg.reactor
    n0 = cfg.n_pt_atoms if n_atoms is None else n_atoms
    nb = n_bins
    frac = (np.arange(nb) + 0.5) / nb
    t_c_bin = r.t_setpoint_c - r.t_gradient_c / 2.0 + r.t_gradient_c * frac
    enrich = 1.0 + r.product_enrichment * (frac - 0.5) * 2.0

    sizes = np.arange(1, n_max + 1)
    y = np.zeros((nb, n_max))
    y[:, 0] = n0

    _, k_ref, _ = site_occupancy(cfg, r.t_setpoint_c + 273.15, 0.0)

    nt = len(t_grid)
    out_pop = np.zeros((nt, 4))
    out_rate = np.zeros(nt)
    out_rate_bin = np.zeros((nt, nb))
    out_mass = np.zeros(nt)

    dt = float(np.mean(np.diff(t_grid))) if nt > 1 else 1.0
    sub = max(1, int(np.ceil(dt / 0.02)))
    h = dt / sub

    cap = sizes ** (2.0 / 3.0)
    for it, t in enumerate(t_grid):
        cond = schedule.conditions(t)
        p_co, p_o2 = cond["p_co"], cond["p_o2"]
        dt_c = cond["t_c"] - r.t_setpoint_c
        for b in range(nb):
            t_k = t_c_bin[b] + dt_c + 273.15
            _, k_eff, mult = site_occupancy(cfg, t_k, dose_rate)
            a1k = k_enc * (k_eff / k_ref) * cap
            kd = np.array([k_diss(cfg, int(n), t_k, p_o2) for n in sizes])
            kd[0] = 0.0
            yb = y[b]
            for _ in range(sub):
                n1 = yb[0]
                flux_up = a1k * n1 * yb            # k -> k+1 (index k-1 -> k)
                flux_up[-1] = 0.0                  # top bin does not grow out
                flux_down = kd * yb                # k -> k-1
                dy = np.zeros(n_max)
                dy -= flux_up
                dy[1:] += flux_up[:-1]
                dy -= flux_down
                dy[:-1] += flux_down[1:]
                # monomer bookkeeping: nucleation removes two monomers,
                # growth removes one, each dissociation releases one
                # nucleation removes a SECOND monomer; every growth event
                # k -> k+1 (k >= 2) also consumes one monomer; every
                # dissociation releases one.  Omitting the growth term is the
                # classic way to break mass conservation in this balance.
                dy[0] += (-flux_up[0] - float(np.sum(flux_up[1:]))
                          + float(np.sum(flux_down[1:])))
                yb = np.maximum(yb + h * dy, 0.0)
            y[b] = yb
        rt = np.zeros(nb)
        for b in range(nb):
            t_k = t_c_bin[b] + dt_c + 273.15
            _, _, mult = site_occupancy(cfg, t_k, dose_rate)
            tofs = np.array([tof_entity(cfg, int(n), "terrace", t_k, p_co, p_o2)
                             for n in sizes]) * mult
            rt[b] = float(np.dot(y[b], tofs)) * enrich[b]
        out_rate[it] = rt.mean()
        out_rate_bin[it] = rt
        ym = y.mean(axis=0)
        out_pop[it] = [ym[0], ym[1], ym[2], ym[3:].sum()]
        out_mass[it] = float(np.dot(ym, sizes))

    return {"t": t_grid, "rate": out_rate, "rate_bin": out_rate_bin,
            "pop": out_pop, "t_c_bin": t_c_bin, "enrich_bin": enrich,
            "mass": out_mass, "dist": y}


def patch_ensemble(cfg: Config, seeds: Sequence[int], duration: float,
                   dose_rate: float = 0.0,
                   schedule: Optional[Schedule] = None,
                   record_fps: float = 1.0) -> Dict[str, np.ndarray]:
    """Average several patch kMC replicas onto a common grid."""
    runs = [PatchSimulator(cfg, seed=s, dose_rate=dose_rate,
                           schedule=schedule).run(duration=duration,
                                                  record_fps=record_fps)
            for s in seeds]
    t = runs[0].t
    rate = np.mean([r.true_rate for r in runs], axis=0)
    pop = np.mean([r.pop for r in runs], axis=0)
    rate_sd = np.std([r.true_rate for r in runs], axis=0, ddof=1)
    return {"t": t, "rate": rate, "rate_sd": rate_sd, "pop": pop,
            "n_replicas": len(runs), "runs": runs}


def calibrate_k_enc(cfg: Config = DEFAULT, seeds: Sequence[int] = (1, 2, 3, 4,
                                                                   5, 6, 7, 8),
                    duration: float = 600.0,
                    grid: Sequence[float] = tuple(np.geomspace(2e-5, 4e-4, 25))
                    ) -> Dict[str, float]:
    """Fit the single encounter coefficient of the chip model to patch kMC.

    The objective is the squared error on the ensemble-mean monomer count and
    rate trajectories of beam-off patch runs.  A one-parameter fit of a
    coarse-grained model to the microscopic model is the standard
    consistency step; the residual is reported so that the quality of the
    coarse-graining is visible (Table 3).
    """
    ens = patch_ensemble(cfg, seeds, duration, dose_rate=0.0, record_fps=1.0)
    t = ens["t"]
    sched = Schedule(cfg)
    best = None
    for k_enc in grid:
        mf = chip_mean_field(cfg, t, sched, k_enc=k_enc, dose_rate=0.0,
                             n_bins=1)
        # normalise each channel by its own scale so neither dominates
        e_pop = np.mean(((mf["pop"][:, 0] - ens["pop"][:, 0])
                         / max(ens["pop"][:, 0].max(), 1e-9)) ** 2)
        e_rate = np.mean(((mf["rate"] - ens["rate"])
                          / max(ens["rate"].max(), 1e-9)) ** 2)
        obj = e_pop + e_rate
        if best is None or obj < best["objective"]:
            best = {"k_enc": float(k_enc), "objective": float(obj),
                    "rmse_pop_frac": float(np.sqrt(e_pop)),
                    "rmse_rate_frac": float(np.sqrt(e_rate))}
    best.update({"n_replicas": len(seeds), "duration_s": duration})
    return best
