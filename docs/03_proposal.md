# 03 — Proposal: Catalytic Event Microscopy (CEM)

## 3.1 One-sentence statement

**Catalytic Event Microscopy** treats an operando electron-microscopy experiment not as
a movie to be interpreted but as an *event-driven causal inference problem*, in which
atomic transitions are the treatments, the transport-corrected product rate is the
outcome, the electron beam and the transport field are named confounders, and the
reactor's actuators are available both as instruments for identification and as
controls for stabilisation.

## 3.2 Central hypothesis

> Catalytic activity of a dynamic single-atom Pt catalyst is governed not by the
> instantaneous population of structural motifs but by the **transition rates and
> lifetimes** among a small number of recurring atomic motifs; therefore a descriptor
> set containing rates and lifetimes will explain measured turnover better than a
> census descriptor set of equal dimension, and an intervention that changes those
> rates will change turnover in the predicted direction.

This is falsifiable in three separate ways: (H1) the dynamic descriptor set fails to
beat the census set; (H2) the event-aligned effect is indistinguishable from the
surrogate null; (H3) an intervention that provably changes motif lifetimes fails to
move the rate in the predicted direction.

## 3.3 Reframing: why "event" is the right unit

| Framing | Unit of analysis | Failure mode |
|---|---|---|
| Movie framing | the whole time series | confounding by T, p; visual overfitting |
| Census framing | a population snapshot | dilution by motif duty cycle (G6) |
| **Event framing** | one atomic transition | requires many events, but supports replication, a null model, and an effect size |

The event framing imports the machinery of a field that already solved the
"one-noisy-trace" problem — event-related averaging in electrophysiology and
event-study methods in econometrics — and it makes the transport delay a *nuisance
parameter to be calibrated* rather than an unstated assumption.

## 3.4 The three layers of the platform

```
                 ┌──────────────────────────────────────────────┐
  Layer 3        │  CONTROL: event-triggered actuation policy    │
  (act)          │  reward = ∫ product − λ_act·actuation − λ_D·dose │
                 └──────────────▲────────────────┬───────────────┘
                                │ D(t), r̂(t)     │ setpoints
                 ┌──────────────┴────────────────▼───────────────┐
  Layer 2        │  INFERENCE: transport deconvolution, event    │
  (infer)        │  alignment, HMM state flux, BPI, regression   │
                 └──────────────▲────────────────▲───────────────┘
                                │ atoms, spectra │ MS channels
                 ┌──────────────┴────────────────┴───────────────┐
  Layer 1        │  SENSE: STEM/4D-STEM · EELS · QMS · electrical │
  (measure)      │  all on ONE hardware time base                 │
                 └───────────────────────────────────────────────┘
```

The non-negotiable design rule is the **single time base**: every stream is timestamped
by the same clock (see `04_experimental_setup.md` §3), because every inference in
Layer 2 is a statement about relative timing.

## 3.5 Descriptor vector

For each analysis bin *t* (default 1 s) the platform emits

D(t) = [ N₁, N₂, N₃, N₄₊, d̄, CN_Pt–Pt, CN_Pt–sup, q_Pt, k_hop, k_nuc, k_diss, τ_dimer, θ_ads, f_rep, φ_dose, T, p_CO, p_O₂, V, I ]

split by provenance into three blocks that are **never mixed silently**:

* **Census (2D-observable):** N₁, N₂, N₃, N₄₊, d̄, CN_Pt–Pt (projected)
* **Dynamic (2D-observable, requires tracking):** k_hop, k_nuc, k_diss, τ_dimer
* **Spectroscopic / reactor:** q_Pt, θ_ads, T, p, V, I, φ_dose
* **3D-derived (intermittent only):** CN_Pt–sup, true-depth-corrected mobility

and each element carries an uncertainty and a provenance tag.

## 3.6 Outcome vector

Y(t) = [ r_prod(t), X_reactant(t), S_prod(t) ] from the deconvolved, calibrated MS,
normalised by the Pt inventory **in the reactive zone**, not the Pt dispensed onto the
chip.

## 3.7 Estimators (what is actually computed)

1. **Transfer function.** Tracer pulse ⇒ h(t) fitted as a gamma / axial-dispersion
   kernel; report τ_dead, σ_disp with CIs per set point.
2. **Deconvolution.** Tikhonov-regularised, with the regularisation weight chosen by
   generalised cross-validation; delay-shift-only fallback reported alongside.
3. **Event-aligned effect.** Δr = mean(r̂ | post) − mean(r̂ | pre), averaged over events,
   with a circular-shift surrogate null and a BCa bootstrap CI over events (and a
   chip-level cluster bootstrap when chips > 1).
4. **Lead–lag.** ρ(τ) = corr[X(t), r̂(t+τ)], interpreted only for |τ| > 2σ_disp.
5. **Multivariable model.** Elastic-net and a state-space variant, with T and p always
   forced into the model so that structural coefficients are read as partial effects.
6. **State flux.** A Gaussian HMM on D(t) gives motif states and a transition matrix;
   the per-state rate contribution is estimated by regressing r̂ on state occupancy.
7. **Beam Perturbation Index.** From the dose-rate series: k_obs(φ) = k_chem + βφ,
   BPI = βφ/(k_chem + βφ).
8. **Interventional estimate.** In the randomised arm, the actuator is the instrument;
   the effect of the structural variable on the rate is estimated by two-stage least
   squares, reported next to the observational elastic-net coefficient.

## 3.8 Control policy (Layer 3)

Action space (discrete, safety-clamped): `hold`, `O₂-pulse(Δt)`, `CO-lean(Δt)`,
`T-step(±ΔT)`, `bias-pulse(V, Δt)`.
Observation: D(t) with its uncertainty, plus r̂(t).
Reward: ∫ r_prod dt − λ_act·(actuation cost) − λ_dose·φ − λ_risk·P(irreversible sintering).
Trigger: a pre-registered precursor signature (rising k_nuc, falling CN_Pt–sup,
growing N₃) crossing a threshold learned from the HMM's pre-failure state.

The claim being tested at Layer 3 is narrow and checkable: **an event-triggered policy
achieves higher integrated product and better end-of-run dispersion than both open-loop
and a fixed periodic schedule delivering the same total actuation budget.** The
"same budget" control is what distinguishes intelligence from merely doing more.

## 3.9 What this project delivers, and what it does not

**Delivers.** The full Layer-1 specification, a complete and executable Layer-2/Layer-3
software stack, a ground-truth digital twin of the experiment, and a quantitative
validation of every estimator against known truth, including its failure modes and its
power as a function of event count.

**Does not deliver.** Any experimental measurement. No physical microscope, nanoreactor
or mass spectrometer was operated. All numbers reported in the manuscript are outputs of
the simulator specified in `04_experimental_setup.md` §5 and implemented in
`src/opcem/`. This is stated in the manuscript abstract, title, and integrity
declaration, and is the reason the manuscript is framed as a *methodology and
in-silico validation* paper rather than an experimental one.

**Why this is still worth publishing.** The estimators in §3.7 have failure modes that
are invisible without ground truth. An experimental paper cannot show that the
event-aligned estimator inverts its sign when the transport delay is uncorrected,
because in an experiment the true answer is unknown. A digital twin can, and that
result is a prerequisite for trusting the same pipeline on real data. The deliverable
is therefore a *pre-registrable analysis protocol with characterised operating limits*.

## 3.10 Model systems and reactions (for the experimental phase that follows)

| System | Why | First reaction |
|---|---|---|
| Pt₁/CeO₂ (100) nanocubes | oxygen-vacancy chemistry; redispersion documented and facet-dependent [R12, R14] | CO oxidation |
| Pt₁/TiO₂ | reducible support, charge transfer | CO oxidation, then C₂H₄ hydrogenation |
| Pt₁/γ-Al₂O₃ | irreversible single-atom→cluster conversion documented with MS [R11]; industrial relevance [R7] | CO oxidation |
| Pt₁/N-doped carbon | fuel-cell relevant | ORR (liquid cell) |
| Pt–M dual atom (Fe, Co, Sn) | ensemble-controlled selectivity | CO₂ reduction |

CO oxidation on Pt/CeO₂ is the designated first demonstration because it has the
largest existing operando-TEM + MS baseline [R1, R2, R3] against which the platform's
event-level output can be cross-checked.

## 3.11 Titles

Primary: **"Catalytic Event Microscopy: transport-corrected, beam-aware causal inference
from operando electron microscopy of single-atom platinum catalysts."**

Alternative (control-forward): "Closing the loop on single atoms: event-triggered
control of dynamic Pt catalysts under operando electron microscopy."
