"""Layer 3: closing the loop onto the chemistry (gap G7).

Three policies with a *shared actuation budget*, which is what separates
"the controller is intelligent" from "the controller did more":

*   :class:`OpenLoop` -- no actuation at all (the reference trajectory).
*   :class:`FixedSchedule` -- the budget spent at evenly spaced times.
*   :class:`EventTriggered` -- the same budget spent when a pre-registered
    precursor signature fires, subject to a refractory period.

The precursor is deliberately simple and pre-registered rather than learned
in-sample: the count of entities of three or more atoms, averaged over a
trigger window, crossing a threshold.  Larger clusters are the irreversible
sink in this system, so their appearance is the earliest observable that
deactivation is under way; an O2 pulse raises the oxidative dissociation rate
and drives them back towards monomers and dimers.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from .config import Config, DEFAULT
from .truth import Actuation, Schedule


@dataclass
class PolicyLog:
    fired_t: List[float] = field(default_factory=list)
    observed_n3: List[float] = field(default_factory=list)
    observed_t: List[float] = field(default_factory=list)

    @property
    def n_actuations(self) -> int:
        return len(self.fired_t)


class BasePolicy:
    """Common interface: ``observe`` is called by the simulator every step."""

    name = "base"

    def __init__(self, cfg: Config = DEFAULT):
        self.cfg = cfg
        self.log = PolicyLog()

    def observe(self, t: float, sizes: np.ndarray, sites: np.ndarray,
                schedule: Schedule) -> None:  # pragma: no cover - interface
        raise NotImplementedError

    def budget_used(self) -> int:
        return self.log.n_actuations


class OpenLoop(BasePolicy):
    name = "open_loop"

    def observe(self, t, sizes, sites, schedule) -> None:
        return None


class FixedSchedule(BasePolicy):
    """Spend the whole budget on an even grid, decided before the run."""

    name = "fixed_schedule"

    def __init__(self, cfg: Config = DEFAULT, duration_s: Optional[float] = None,
                 budget: Optional[int] = None):
        super().__init__(cfg)
        dur = cfg.duration_s if duration_s is None else duration_s
        b = cfg.control.actuation_budget if budget is None else budget
        self.times = list(np.linspace(dur / (b + 1), dur * b / (b + 1), b)) \
            if b > 0 else []
        self._installed = False

    def observe(self, t, sizes, sites, schedule) -> None:
        if self._installed:
            return
        c = self.cfg.control
        for tt in self.times:
            schedule.actuations.append(
                Actuation(t_start=float(tt), duration=c.o2_pulse_duration_s,
                          kind="o2_pulse", x_o2=c.o2_pulse_x_o2))
            self.log.fired_t.append(float(tt))
        self._installed = True


class EventTriggered(BasePolicy):
    """Fire an O2 pulse when the deactivation precursor crosses threshold.

    The observation is the *observed* cluster count, i.e. what the detection
    front end reports, not the ground truth; the controller therefore inherits
    the detection bias, which is the realistic situation.
    """

    name = "event_triggered"

    def __init__(self, cfg: Config = DEFAULT, budget: Optional[int] = None,
                 threshold: Optional[float] = None,
                 detect_recall_n3: float = 1.0,
                 rng: Optional[np.random.Generator] = None):
        super().__init__(cfg)
        c = cfg.control
        self.budget = c.actuation_budget if budget is None else budget
        self.threshold = c.trigger_n3_threshold if threshold is None \
            else threshold
        self.window = c.trigger_window_s
        self.refractory = c.refractory_s
        self.recall = detect_recall_n3
        self.rng = rng or np.random.default_rng(0)
        self._hist: List[tuple] = []
        self._last_fire = -np.inf

    def observe(self, t, sizes, sites, schedule) -> None:
        n3 = float(np.sum(np.asarray(sizes) >= 3))
        if self.recall < 1.0:
            n3 = float(self.rng.binomial(int(n3), self.recall))
        self._hist.append((t, n3))
        self._hist = [h for h in self._hist if h[0] >= t - self.window]
        self.log.observed_t.append(t)
        self.log.observed_n3.append(n3)
        if self.log.n_actuations >= self.budget:
            return
        if t - self._last_fire < self.refractory:
            return
        if len(self._hist) < 3:
            return
        if float(np.mean([h[1] for h in self._hist])) >= self.threshold:
            c = self.cfg.control
            schedule.actuations.append(
                Actuation(t_start=float(t), duration=c.o2_pulse_duration_s,
                          kind="o2_pulse", x_o2=c.o2_pulse_x_o2))
            self.log.fired_t.append(float(t))
            self._last_fire = t


def reward(cfg: Config, t: np.ndarray, rate: np.ndarray,
           n_actuations: int, dose_rate: float,
           final_dispersed_frac: float) -> Dict[str, float]:
    """Integrated product, and the components of the control objective.

    Reported as separate columns rather than collapsed into one number, so
    that a policy cannot look good by trading away a term the reader cares
    about.
    """
    dt = float(np.mean(np.diff(t))) if len(t) > 1 else 1.0
    yield_molec = float(np.nansum(rate) * dt)
    return {"integrated_yield": yield_molec,
            "mean_rate": float(np.nanmean(rate)),
            "final_rate": float(np.nanmean(rate[-max(1, len(rate) // 20):])),
            "n_actuations": int(n_actuations),
            "dose_rate": float(dose_rate),
            "final_dispersed_frac": float(final_dispersed_frac)}
