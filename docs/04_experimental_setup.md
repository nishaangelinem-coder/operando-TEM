# 04 — Experimental Setup: Hardware, Software, and the Digital Twin

This document has two halves. **§1–§4** specify the *physical* platform, to the level of
detail needed for procurement and for a facility proposal. **§5–§7** specify the
*digital twin* that stands in for that platform in this study, including the exact
parameter values used, so that every simulated number in the manuscript is traceable.

---

## 1. Physical architecture

```
 MFC bank ──┐
 (CO, O2,   ├─► 4-port fast switching valve ─► MEMS nanoreactor ─► heated capillary ─► QMS
  He, Ar,   │        (t_sw < 50 ms)              (in TEM column)      (150 °C)        │
  13CO)  ───┘                │                        │                              │
                             │                        ├─ ADF-STEM / 4D-STEM detector  │
                        TTL trigger out               ├─ EELS spectrometer            │
                             │                        └─ on-chip T sensor, heater I/V │
                             ▼                                    │                   │
                    ┌────────────────────────────────────────────────────────────┐    │
                    │  FPGA / DAQ master clock (10 MHz, all streams timestamped) │◄───┘
                    └────────────────────────────────────────────────────────────┘
                             │
                             ▼
                    Edge compute node (GPU) ── Layer 2 inference ── Layer 3 policy
                             │                                          │
                             └────────── setpoints back to MFC / heater / potentiostat
```

### 1.1 Microscope and detectors

| Item | Specification | Purpose |
|---|---|---|
| Probe-corrected STEM | 200/300 kV, cold FEG, probe ≤ 0.8 Å, monochromated option | single-atom ADF imaging; 60–80 kV low-voltage option to reduce knock-on |
| Direct-electron pixelated detector | ≥ 2 kHz full-frame, 4D-STEM capable | ptychography, low-dose phase contrast |
| ADF detector with fast DAC readout | per-pixel timestamping | continuous 2D atom tracking |
| EELS spectrometer | energy resolution ≤ 30 meV (monochromated) / ≤ 0.8 eV standard; dual-EELS | Pt M₄,₅ & O K core loss (oxidation state); low-loss gas fingerprint; vibrational adsorbate bands |
| Beam blanker | electrostatic, ≤ 1 µs rise, TTL-driven | stroboscopic duty cycling, beam-blanked control periods |
| Stage | piezo, drift ≤ 0.2 nm/min after settling | multi-position stratified sampling |

### 1.2 Nanoreactor (the critical component)

| Parameter | Target | Rationale |
|---|---|---|
| Geometry | unidirectional flow channel, separate inlet/outlet | defines a conversion (G1, G5) |
| Windows | SiN_x, 8–15 nm, ≥ 3 viewing windows along the flow axis | position-stratified design (G5) |
| Channel | ~4–6 µm deep × ~30–50 µm wide, reactive zone ~200 µm long | keeps internal volume ≤ ~100 nL to limit σ_disp |
| Pressure | up to 1 bar (1000 mbar) | technical relevance; matches [R1] |
| Temperature | RT–800 °C, on-chip Pt-resistor thermometry, 4-wire | local T is a required descriptor |
| T uniformity | mapped, not assumed; target ≤ 5 °C across reactive zone | microheaters have real in-plane gradients [R6] |
| Dead volume, chip→MS | minimised; ≤ 20 nL internal + capillary | τ_dead is the dominant nuisance parameter |
| Electrical feedthroughs | ≥ 6, for electrocatalysis variant (WE/CE/RE + sense) | ORR/CO₂RR extension |

### 1.3 Gas handling and product analytics

| Item | Specification |
|---|---|
| Mass-flow controllers | 0.1–10 sccm, 1 % FS, ≤ 100 ms settling; one per gas |
| Switching valve | 4-port, pneumatic, < 50 ms, TTL-triggered, **timestamp published to the master clock** |
| Tracer supply | Ar, He, and isotopically labelled ¹³CO / ¹⁸O₂ for channel-unambiguous pulses |
| Transfer line | fused-silica capillary, 150 °C jacket, minimal length | prevents condensable loss (H₂O, oxygenates) |
| QMS | 1–200 amu, ≥ 10 Hz multi-channel, secondary-electron multiplier, differentially pumped inlet |
| Optional PTR-MS | for oxygenate selectivity in CO₂RR work |
| Calibration standards | certified CO/O₂/CO₂/Ar mixtures spanning the expected composition range |

### 1.4 Electrochemistry variant

Potentiostat with ≥ 10 kHz sampling, floating ground, EIS 0.1 Hz–100 kHz, current
resolution ≤ 1 pA; liquid-cell chip with the same three-window geometry.

---

## 2. Measurement protocol (pre-registered order)

| # | Step | Output | Gap addressed |
|---|---|---|---|
| 0 | Chip QC: window thickness map, heater calibration, leak check | chip acceptance record | — |
| 1 | Catalyst deposition; Pt inventory in reactive zone by ICP-MS on sacrificial chips + ADF atom counting | n_Pt, f_rep | G5 |
| 2 | MS blank, fragmentation pattern, background, T-dependent baseline | I_bg, cracking matrix | G1 |
| 3 | MS calibration against certified mixtures at each (T, p) | k_j | G1 |
| 4 | **Tracer-pulse transfer-function measurement at every (T, p, flow) set point** | h(t), τ_dead, σ_disp | **G1** |
| 5 | Beam-off reaction map: conversion vs T and composition | intrinsic kinetics baseline | G4 |
| 6 | Low-dose operando run: synchronised STEM + EELS + MS + reactor logs | primary dataset | — |
| 7 | Dose-rate series at fixed (T, p, flow), 5 points ≥ 1 decade | BPI, k_chem | **G4** |
| 8 | Position-stratified repeats: inlet / centre / outlet windows | heterogeneity test | **G5** |
| 9 | Programmed perturbations (T ramps, CO/O₂ steps, isotope pulses, randomised arm) | interventional estimate | **G8** |
| 10 | Intermittent 3D: multi-tilt ptychographic snapshots at defined times | depth prior | **G3** |
| 11 | Ensemble cross-validation outside the microscope: operando XAS / DRIFTS / Raman on the same batch | external validity | G5 |

Step 4 before step 6 is not negotiable: without h(t) the primary dataset cannot be
analysed at event resolution.

---

## 3. Synchronisation: the single time base

* One FPGA/DAQ device (e.g. 10 MHz counter) is the **master clock**.
* It emits a periodic TTL "frame sync" consumed by the STEM scan generator, the EELS
  spectrometer trigger, the QMS acquisition, the potentiostat, and the valve controller.
* Every stream writes `(master_tick, local_index, payload)`. No stream is trusted to
  report its own wall-clock time.
* Valve actuations, beam-blank transitions and setpoint changes are logged as
  *events on the same tick axis*, to sub-10 ms.
* Acceptance test: a deliberate valve step must appear in the reactor log and in the
  QMS channel with a difference equal to the independently measured τ_dead to within
  σ_disp.

Nominal sampling targets: ADF/STEM 1–30 fps (dose-limited), EELS 1 spectrum per frame
or per N frames, QMS 10 Hz, reactor scalars 10 Hz, trigger log event-driven at 100 Hz
resolution.

---

## 4. Software architecture (Layer 2 + Layer 3)

```
src/opcem/
  truth.py        ground-truth generative model (digital twin only)
  instruments.py  virtual ADF-STEM, EELS, QMS forward models + noise
  transport.py    h(t) model, tracer calibration, delay estimation, Tikhonov deconvolution
  vision.py       atom detection (LoG + refinement), nearest-neighbour linking, tracking
  descriptors.py  D(t) assembly with provenance tags and uncertainties
  correlate.py    event alignment, surrogate nulls, bootstrap, lead-lag, elastic net, HMM, 2SLS
  beam.py         dose-series fit, Beam Perturbation Index
  control.py      event-triggered policy, fixed-schedule and open-loop baselines
  report.py       tables and figures
```

Design rules enforced in code:

1. **Provenance tags.** Every descriptor column is tagged `census2d`, `dynamic2d`,
   `spectro`, `reactor` or `snapshot3d`. Any estimator that mixes `snapshot3d` into a
   continuous-time model must declare it.
2. **No silent interpolation** of the intermittent 3D stream onto the 1 s grid.
3. **Delay uncertainty propagates.** Lags are reported with the transfer-function CI
   folded in; `correlate.py` refuses to label a lag significant if |τ| < 2σ_disp.
4. **Seeded RNG everywhere**, seed recorded in every output file.
5. **Pre-registered thresholds** live in `config.py`, not in analysis scripts.

Runtime budget for closed-loop operation: detection + descriptor extraction must
complete within one analysis bin (1 s) on a single GPU; the policy is a lookup on the
HMM posterior, so its cost is negligible.

---

## 5. The digital twin (what is actually run in this study)

### 5.1 Ground-truth generative model (`truth.py`)

A continuous-time Markov (kinetic Monte Carlo, Gillespie) model of Pt species on a
support lattice, coupled to a Langmuir–Hinshelwood CO-oxidation microkinetic model.

**State.** A list of Pt entities, each with size n ∈ {1,2,3,…}, a 3D position on/in the
support slab, and a support-site type s ∈ {terrace, O-vacancy, step}.

**Reaction channels.**

| Channel | Rate law | Notes |
|---|---|---|
| Hop of a monomer | k_hop = ν₀ exp(−E_a(s)/k_BT) + β·φ | φ = electron dose rate; **this is the beam confound** |
| Dimerisation Pt₁+Pt₁→Pt₂ | diffusion-limited encounter within r_cap | |
| Growth Pt_n+Pt₁→Pt_{n+1} | capture radius ∝ n^{1/3} | |
| Dissociation Pt_n→Pt_{n−1}+Pt₁ | k_diss(n) = ν₀ exp(−E_b(n)/k_BT)·(1+γ·p_O₂^{1/2}) | oxidative redispersion |
| Trap/untrap at O-vacancy | E_a raised by E_trap on a vacancy site | anchoring |

**Turnover.** Each entity contributes a site-specific turnover frequency

  TOF_n(T, p) = A_n · θ_CO · θ_O · exp(−E_act,n / k_B T)

with Langmuir coverages θ_CO, θ_O from p_CO, p_O₂ and adsorption energies that depend
on n. The **ground-truth intrinsic activity ordering is set by A_n and E_act,n**, and in
the baseline configuration the *dimer* is the most productive motif per Pt atom —
this is the known answer that the pipeline must recover without being told.

Total instantaneous rate: R(t) = Σ_entities TOF_{n,s}(T(x), p(x)) — i.e. the rate is
position-weighted, which is what creates the representativeness problem (G5).

### 5.2 Virtual instruments (`instruments.py`)

| Stream | Forward model |
|---|---|
| ADF-STEM frame | Gaussian probe blur (σ = 0.45 Å), per-atom intensity ∝ Z^1.7·n_z, Poisson shot noise at the requested dose, scan-line jitter, slow stage drift, SiN window background with 1/f spatial noise |
| Projected positions | true 3D positions projected along the optic axis ⇒ depth information destroyed (implements G3) |
| EELS | Pt M₄,₅ white-line ratio mapped from ground-truth oxidation state + Gaussian noise scaled by dose; O K pre-peak from vacancy concentration |
| QMS | R(t) → species fluxes → convolution with h(t) → ×k_j calibration → + 1/f baseline drift + Gaussian noise + channel cross-talk (m/z 28 N₂ interference) → 10 Hz sampling |
| Reactor scalars | T(x) field with an imposed inlet→outlet gradient; p, flow, and the actuator log |

### 5.3 Transport transfer function (`transport.py`)

h(t) is modelled as a gamma kernel with a pure dead time:

  h(t) = 0 for t < τ_0 ;  h(t) = (t−τ_0)^{k−1} e^{−(t−τ_0)/θ} / [Γ(k) θ^k] for t ≥ τ_0

so that τ_dead = τ_0 + kθ and σ_disp = θ√k. This is the standard tanks-in-series /
axial-dispersion form and is sufficient to reproduce the dead-time-plus-tailing
behaviour of a nanoreactor + heated capillary + differentially pumped inlet.

### 5.4 Parameter values used

All values are in `src/opcem/config.py`. The physically-anchored ones and their
justification:

| Parameter | Value | Anchor |
|---|---|---|
| Reactor pressure | 700 mbar | within the demonstrated operando-TEM envelope [R1] |
| Temperature | 250–400 °C, baseline 300 °C | CO-oxidation light-off range for Pt |
| τ_dead | 3.202 s at 2 sccm, scaling as 1/flow | order of magnitude of nanoreactor + capillary + QMS inlet |
| σ_disp | 0.838 s at 2 sccm | consistent with a ~100 nL cell plus capillary at a few sccm |
| QMS rate | 10 Hz | typical multi-channel QMS |
| Analysis bin | 1.0 s | matches QMS and Layer-3 budget |
| Dose-rate series | 1×10² to 1×10⁴ e⁻ Å⁻² s⁻¹, 5 points | brackets low-dose practice [R10] far below manipulation regimes [R15] |
| Lateral localisation σ | 0.25 Å at reference dose, ∝ dose^(−1/2) | achievable single-atom precision |
| Depth "resolution" | 6.6 nm for single-projection snapshots | measured value [R8] |
| Product enrichment, inlet→outlet | 18 % | within the 12–21 % range reported for cell non-uniformity [R4] |
| Projected area fraction | 1.6×10⁻⁷ (one 40 × 40 nm field of a 200 × 50 µm zone) | computed from the chip geometry above |
| Catalyst area multiplier | 3.5×10⁴ | Pt-accessible area of a ~100 nm film of 5 nm nanocrystals (equivalently a few µg of a ~50 m² g⁻¹ powder) relative to the projected zone area. **Estimated, not measured**; study S7 reports the conclusion's sensitivity to it over four decades |
| f_rep | 4.6×10⁻¹² | projected area fraction divided by the catalyst area multiplier |
| QMS detection limit | 1 ppm mole fraction | a quadrupole measures a *concentration*, so its flux limit is x_min × F and scales with the flow |
| Hop barrier, terrace | 1.32 eV | chosen so that the monomer hop rate is a few s⁻¹ at 300 °C, i.e. resolvable at 5 frames s⁻¹ |
| Beam coupling β | 2.0×10⁻⁴ s⁻¹ per e⁻ Å⁻² s⁻¹ | chosen so that the beam contributes ~11 % of hops at the reference dose and ~49 % at 10⁴ e⁻ Å⁻² s⁻¹ — a deliberately strong confound |
| Support corrugation | A = 0.70 nm, λ = 5.0 nm | gives η = A²q²/4 = 0.19, i.e. a ~16 % projection bias in mean-squared displacement |

**These are plausible, cited-where-possible values for a *simulation*, not
measurements.** No claim is made that the twin reproduces any specific published
dataset.

---

## 6. Study list

| ID | Question | Key output |
|---|---|---|
| S0 | Does the fast observation model reproduce the pixel-level detector? | Table S0 |
| S1 | Can h(t) be recovered from tracer pulses, at every flow? | Table S1 |
| S2 | How biased is each descriptor, and do dynamic descriptors earn their place? | Tables S2a–S2d |
| S3 | When does the event-aligned estimator work, and where does its attenuation come from? | Tables S3a–S3c |
| S4 | What does skipping the transport correction cost? | Tables S4a–S4b |
| S4b | Does a multivariable observational estimate recover a known coefficient, and is a gas pulse a usable instrument? | Table S4c |
| S5 | How biased is projected mobility, and can 3D snapshots correct it? | Table S5 |
| S6 | Can the chemical rate be separated from the beam's, and over what dose window? | Tables S6a–S6c |
| S7 | Is the imaged field representative, and can a single event ever be seen in the product? | Tables S7a–S7g |
| S8 | Does closed-loop control beat a budget-matched schedule out of sample? | Tables S8a–S8b |

---

## 7. Reproducibility

`python studies/run_all.py` regenerates every table and figure from seed `20260922`.
Each output file records the code version (git SHA), the seed, and the full config dict.
