---
title: "Catalytic Event Microscopy: transport-corrected, beam-aware causal inference from operando electron microscopy of single-atom platinum catalysts — a methodology and in-silico validation"
author:
  - name: M. Nisha Angeline
    affiliation: Department of Electronics and Communication Engineering, Velalar College of Engineering and Technology, Thindal, Erode, Tamil Nadu, India
    email: nishavlsidesign@gmail.com
date: 2026-09-22
keywords: operando electron microscopy, single-atom catalysis, platinum, causal inference, mass spectrometry, closed-loop control, digital twin
---

> **Declaration of simulated data (read this first).**
> Every numerical result in this manuscript is the output of a forward
> simulator, not a measurement. No physical electron microscope, MEMS
> nanoreactor, mass spectrometer or catalyst was operated at any point in this
> work. The simulator, its parameters and its provenance are specified in
> §3 and implemented in the accompanying code; every table and figure file
> carries a machine-readable banner to the same effect. The contribution
> claimed here is a measurement-and-analysis protocol with characterised
> operating limits, together with the ground-truth validation of each
> estimator in it — not any experimental finding about platinum.

---

## Abstract

**This is a simulation study: every number reported below is the output of a
declared forward model, and no physical instrument was operated.** Operando
electron microscopy can now watch individual platinum atoms move while a
catalyst works, and can record the reaction products at the same time.
It cannot yet establish which atom or which configuration made the product.
The obstacle is not spatial resolution. It is that the two data streams are
compared by visual temporal coincidence, while three effects sit between them:
the gas-transport delay between the imaged region and the mass spectrometer,
the electron beam's own contribution to the atomic dynamics, and the tiny and
non-representative fraction of the working catalyst that any one field of view
contains.

We formalise the problem as event-driven causal inference — atomic transitions
are treatments, the transport-corrected product rate is the outcome, the beam
and the transport field are named confounders, and the reactor's actuators are
available both as instruments for identification and as controls for
stabilisation — and we specify the hardware and software platform this requires.
We then build a two-level digital twin of that platform (a spatially resolved
kinetic Monte Carlo model of Pt entities in one illuminated field of view, and
a mass-conserving population balance for the whole un-illuminated reactive
zone) and validate every estimator against known truth.

The validation gives six quantitative results. (i) Timestamped tracer pulses
recover the transport transfer function to better than 5 ms in dead time and
10 ms in dispersion across a fourfold flow range, with nominal interval
coverage. (ii) Omitting that correction attenuates the event-aligned effect and
places the apparent lead–lag peak near the transport dead time, which would be
read as structure leading chemistry by that amount. (iii) A dose-rate series
recovers the chemical hop rate as its zero-dose intercept and yields a Beam
Perturbation Index that tracks the simulator's true beam-attributed event
fraction. (iv) Projected image data under-estimate single-atom
mean-squared displacement, and the bias is corrected using only the sparse
3D snapshot stream. (v) The central hypothesis is *partly refuted*: dynamic
descriptors add essentially nothing to the explanation of the instantaneous
rate, which a motif census already determines, but they do add forecast skill
at longer horizons. (vi) Most consequentially, spontaneous single-atom event
analysis against a quadrupole mass spectrometer is infeasible by roughly nine
orders of magnitude in detector sensitivity; the route that survives is
perturbation-synchronised analysis, in which a global actuator step drives the
same transition across the whole chip and the imaged field becomes a sampler of
a synchronised response rather than the source of the signal. Closed-loop
control of dispersion is demonstrated against baselines that spend an identical
actuation budget.

The deliverable is a pre-registrable protocol with its power, bias and failure
modes measured in advance, and open code that regenerates every number.

---

## 1. Introduction

### 1.1 What operando electron microscopy already does

Sealed MEMS nanoreactors made *operando* electron microscopy possible by giving
the gas a defined inlet-to-outlet path, so that a conversion can be defined at
all. The canonical demonstration remains the oscillatory CO oxidation study of
Vendelbo and co-workers [R1], which recorded time-resolved high-resolution
imaging of individual Pt nanoparticles at approximately 1 bar, quantitative
online mass spectrometry of the effluent, and on-chip calorimetry, and found
periodic refacetting of the particles synchronous with periodic oscillation of
the CO, O₂ and CO₂ signals. Subsequent work extended this to the chemical (not
only morphological) dynamics of Pt during CO oxidation and to other metals
[R2, R26], and quantitative gas-phase electron energy-loss spectroscopy
established a near-sample composition probe validated against finite-element
models of the cell's transport fields [R4, R5].

On the structural side, statistical atom counting is mature: electron-microscopy
atom-recognition statistics counted more than 18 000 Pt atoms on an industrial
Pt/Al₂O₃ reforming catalyst and correlated aromatics production quantitatively
with the *density of Pt single atoms*, with clusters contributing no direct
activity [R7]. On the chemistry side, operando X-ray absorption and infrared
spectroscopies with mass spectrometry resolved dispersion, oxidation state and
activity together, showing that poorly active atomically dispersed Pt on
alumina converts gradually and irreversibly into highly active ~1 nm clusters
over heating–cooling cycles [R11] — whereas on ceria the conversion of active
clusters into single atoms can be a *deactivation* pathway [R12, R13], and
facet-dependent atomic redispersion of Pt on CeO₂ has been imaged directly
[R14].

### 1.2 What it does not do

The sign of the structure–activity relationship for "single atom versus
sub-nanometre cluster" is therefore system- and reaction-dependent and is not
settled. Deciding it is exactly the kind of question an atom-resolved operando
platform ought to answer. It has not, and the reasons are methodological rather
than instrumental.

**The comparison is visual.** With very few exceptions, a structural time series
and a product time series are plotted together and coincidence is read off by
eye or from a zero-lag correlation. The periodicity argument of [R1] survives
this, because a constant phase shift does not destroy a frequency match. An
*event-based* argument does not: the measured mass-spectrometer signal is the
true outlet composition convolved with the impulse response of reactor volume,
transfer line and differentially pumped inlet, and for a 1 bar nanoreactor the
dead time is of the order of seconds — the same order as the structural events
themselves.

**The beam is a reagent.** Beam-induced single-atom dynamics are well documented
and have been developed into a deliberate tool: directed single-dopant
manipulation in silicon proceeds by indirect exchange, and the hopping is
stochastic, with dozens of null irradiations before a success [R15, R16]. The
dose rates used for deliberate manipulation are of order 10⁷ e⁻ Å⁻² s⁻¹ [R15],
whereas low-dose ptychography of beam-sensitive materials operates at ~10²
e⁻ Å⁻² in total [R10]. Operando practice sits between these, and separates the
beam's contribution from the chemistry's only qualitatively.

**Three dimensions are not available continuously.** Single-projection
multislice electron ptychography reaches ~0.85 Å laterally but only ~6.6 nm in
depth at ~3.5 × 10³ e⁻ Å⁻² [R8]; sub-nanometre depth resolution requires
coupling several small-angle projections [R9]; atomic-resolution tomography
requires a tilt series, which requires the object to hold still. A mobile single
atom does not. Continuous 3D atom tracking under reaction conditions therefore
does not exist, and mixing continuously-tracked projected positions with
intermittently reconstructed 3D coordinates produces irreproducible mobility
statistics.

**One field of view is not the reactor.** Transport in the cell is
diffusion-dominated and the product field is demonstrably non-uniform: the
reported CO₂ enrichment of the operando pellet relative to the TEM grid spans
roughly 12–21 % [R4], and microheater chips carry real in-plane temperature
distributions [R6]. How much of the measured signal the imaged atoms actually
produce is, in the literature, essentially never stated.

**The loop is open.** Autonomous and agentic microscopy is advancing quickly
[R18, R20, R21], including reward-based image analysis and dose-budgeted
sequential decision making [R22, R23] — but the loop closes onto imaging
objectives (find the defect, place the probe, budget the dose), not onto the
chemistry.

### 1.3 Contribution

We take these five items as specifications rather than caveats.

1. **A reframing.** Catalytic Event Microscopy treats the experiment as
   event-driven causal inference, with atomic transitions as treatments, a
   transport-corrected rate as the outcome, the beam and the transport field as
   named confounders, and the reactor's actuators as both instruments and
   controls (§2).
2. **A platform specification** — hardware, synchronisation, and a
   pre-registered measurement order in which the transfer function is measured
   *before* the primary dataset (§3.1–§3.3).
3. **A two-level digital twin** of that platform, with ground truth available
   at every step, and a documented one-parameter coarse-graining that links the
   microscopic and reactor-scale levels (§3.4).
4. **Ground-truth validation of eight estimators**, including their bias, their
   power as a function of event count, and the conditions under which each one
   fails (§4).
5. **A negative result that changes the experimental design.** The
   representativeness and sensitivity budget rules out spontaneous
   single-atom event analysis against a quadrupole mass spectrometer by about
   nine orders of magnitude, and identifies perturbation-synchronised analysis
   as the route that survives (§4.6).
6. **A closed-loop demonstration** against baselines matched on actuation
   budget, so that the comparison isolates *when* actuation is spent rather than
   how much (§4.7).

What this paper does not contain is any measurement. §5.4 states plainly what
would have to be done experimentally before any of the scientific questions in
§1.2 could be answered, and §5.3 lists the ways the twin could be wrong in a
manner that would invalidate the validation.

---

## 2. The Catalytic Event Microscopy framework

### 2.1 Central hypothesis and its falsifiers

> Catalytic activity of a dynamic single-atom Pt catalyst is governed not by the
> instantaneous population of structural motifs but by the transition rates and
> lifetimes among a small number of recurring motifs.

Three independent falsifiers were pre-specified: (H1) a descriptor set
containing rates and lifetimes fails to out-explain a census set of equal
dimension; (H2) the event-aligned effect is indistinguishable from a matched
surrogate null; (H3) an intervention that provably changes motif lifetimes fails
to move the rate in the predicted direction. §4.3 reports that H1 is *satisfied*
for the instantaneous rate — the hypothesis fails that test — and §4.3 also
reports the weaker sense in which the dynamic descriptors do carry information.

### 2.2 Why the event is the right unit

| Framing | Unit | Failure mode |
|---|---|---|
| Movie | the whole time series | confounding by T and p; visual overfitting |
| Census | a population snapshot | dilution by the motif's duty cycle |
| **Event** | one atomic transition | needs many events, but supports replication, a null model and an effect size |

The event framing imports machinery from fields that already solved the
one-noisy-trace problem (event-related averaging, event-study designs), and it
converts the transport delay from an unstated assumption into a nuisance
parameter that must be calibrated.

### 2.3 Three layers

*Layer 1 (sense)* records ADF/4D-STEM, EELS, quadrupole MS and the reactor
scalars on **one hardware time base**. *Layer 2 (infer)* performs transport
calibration and inversion, event alignment, state-flux modelling, beam
accounting and regression. *Layer 3 (act)* runs an event-triggered actuation
policy whose reward is integrated product minus actuation and dose penalties.
The single time base is not a convenience: every inference in Layer 2 is a
statement about relative timing.

### 2.4 Descriptor and outcome vectors

For each 1 s analysis bin the platform emits

D(t) = [N₁, N₂, N₃, N₄₊, d̄, CN_Pt–Pt, CN_Pt–sup, q_Pt, k_hop, k_nuc, k_diss, τ_dimer, θ_ads, f_rep, φ, T, p_CO, p_O₂, V, I]

partitioned by provenance into `census2d`, `dynamic2d`, `spectro`, `reactor`
and `snapshot3d` blocks. The partition is enforced in code: an estimator that
would consume the intermittent 3D stream in a continuous-time model must declare
it, and the 3D stream is never interpolated onto the 1 s grid. The outcome is
Y(t) = [r_prod, X_reactant, S_prod], obtained from the calibrated and inverted
MS channels and normalised by the Pt inventory *in the reactive zone*.

### 2.5 Estimators

| # | Estimator | Addresses |
|---|---|---|
| E1 | tracer-calibrated gamma transfer function, moments with bootstrap CIs | transport delay |
| E2 | Tikhonov deconvolution with GCV-selected weight; delay-only fallback reported alongside | transport blur |
| E3 | event-aligned effect with local pre-window detrending, circular-shift null and bootstrap CI | no statistics |
| E4 | lead–lag cross-correlation, gated at 2 σ_disp | over-reading of lags |
| E5 | elastic net with reactor terms partialled out | confounding by the operating point |
| E6 | Gaussian hidden Markov model on D(t), giving motif regimes and their transition matrix | static descriptors |
| E7 | dose-series intercept and the Beam Perturbation Index | beam confound |
| E8 | two-stage least squares instrumented by a randomised actuator | association vs intervention |

Two rules are enforced in code rather than left to the analyst: a lead–lag peak
smaller than 2 σ_disp is refused the label "significant", and an event-aligned
effect is always reported next to a null built from the same data.

---

## 3. Methods

### 3.1 Platform specification (Layer 1)

The physical platform is specified in full in `docs/04_experimental_setup.md`;
the elements that the analysis depends on are these.

**Microscope.** Probe-corrected STEM at 60–300 kV with a ≤ 0.8 Å probe, a
direct-electron pixelated detector for 4D-STEM and ptychography, a dual-EELS
spectrometer (≤ 30 meV monochromated, for vibrational adsorbate bands; ≤ 0.8 eV
for Pt M₄,₅ and O K core loss), and a TTL-driven electrostatic beam blanker with
a ≤ 1 µs rise time for stroboscopic duty cycling.

**Nanoreactor.** A unidirectional-flow MEMS chip with separate inlet and outlet,
8–15 nm SiNₓ windows, and — the one non-standard requirement — **at least three
viewing windows along the flow axis**, so that the position-stratified design of
§4.6 is possible on a single chip. Internal volume is held below ~100 nL to
limit the transport dispersion; the reactive zone is 200 µm × 50 µm; on-chip
four-wire Pt-resistor thermometry is mapped rather than assumed, because
microheater chips carry real in-plane gradients [R6].

**Product analytics.** Mass-flow controllers with ≤ 100 ms settling; a 4-port
pneumatic switching valve with < 50 ms actuation **whose timestamp is published
to the master clock**; a 150 °C jacketed fused-silica transfer capillary, to
prevent loss of condensables; a differentially pumped quadrupole MS sampling at
≥ 10 Hz; and isotopically labelled ¹³CO / ¹⁸O₂ for channel-unambiguous tracer
pulses.

**Synchronisation.** One FPGA/DAQ device at 10 MHz is the master clock. It emits
a frame-sync TTL consumed by the STEM scan generator, the EELS trigger, the QMS
acquisition, the potentiostat and the valve controller; every stream writes
`(master_tick, local_index, payload)` and **no stream is trusted to report its
own wall-clock time**. The acceptance test is that a deliberate valve step
appears in the reactor log and in the QMS channel separated by the
independently measured dead time, to within the measured dispersion.

### 3.2 Pre-registered measurement order

The order matters, and one item in it is not negotiable. Chip QC → catalyst
deposition and Pt inventory (ICP-MS on sacrificial chips, cross-checked by ADF
atom counting) → MS blank, cracking pattern and background → MS calibration
against certified mixtures at every set point → **tracer-pulse measurement of
the transfer function at every (T, p, flow) set point** → beam-off reaction map
→ low-dose operando run → dose-rate series → position-stratified repeats →
programmed perturbations including a randomised arm → intermittent multi-tilt 3D
snapshots → ensemble cross-validation outside the microscope by operando XAS,
DRIFTS or Raman on the same batch. The transfer-function step precedes the
primary dataset because without it the primary dataset cannot be analysed at
event resolution at all.

### 3.3 Software (Layers 2 and 3)

The stack is nine modules (`src/opcem/`) with five rules enforced in code:
every descriptor column carries a provenance tag; the intermittent 3D stream is
never interpolated onto the continuous grid; transfer-function uncertainty
propagates into every reported lag and the lead–lag estimator refuses to label
an unresolvable lag significant; all randomness is seeded and the seed is
written into every output file; and all pre-registered thresholds live in a
single configuration module rather than in analysis scripts.

### 3.4 The digital twin

Because no experiment was performed, the entire evaluation runs against a
forward model with two levels. **This is the load-bearing methodological choice
of the paper and also its principal limitation** (§5.3).

**Level 1: the imaged field (`truth.py`, `PatchSimulator`).** A discrete-time
kinetic Monte Carlo model of Pt entities on a 40 × 40 nm support patch with
periodic boundaries. Entities carry a size, a 3D position, and a support-site
type drawn from {terrace, oxygen vacancy, step}. Channels are: monomer hopping
at k_hop = ν₀ exp(−[E_a + E_trap(s)]/k_BT) + β φ, where φ is the electron dose
rate — **this term is the beam confound, present by construction**;
diffusion-limited encounter and merging within a size-scaled capture radius; and
single-atom dissociation at k_diss(n) = ν₀ exp(−E_b(n)/k_BT)·(1 + γ √p_O₂),
which is what makes an oxygen pulse redisperse clusters. The support surface is
corrugated as z = A sin(qx) sin(qy), so lateral diffusion carries a genuine
out-of-plane component that projection destroys. Turnover is

TOF_n(T, p) = A_n θ_CO θ_O exp(−E_act,n / k_BT) · m(s)

with competitive Langmuir coverages and a site-type multiplier. At the operating
point the **dimer is the most productive motif per Pt atom** (0.0447 s⁻¹ per
atom, versus 0.0168 for the monomer, 0.0041 for the trimer and 0.0011 for n ≥ 5).
The pipeline is never told this; recovering it is the validation target. Tests
verify Pt-atom conservation, the motif hierarchy and the linearity of the beam
term.

**Level 2: the whole reactive zone (`chip_mean_field`).** A size-resolved
Smoluchowski-type population balance, restricted to monomer addition and
single-atom loss, using the kMC's own dissociation constants. The zone is split
into twelve slices along the flow axis, each with its own temperature and
product-enrichment factor. Because the physical chip holds order 10¹³ Pt atoms,
the population is deterministic: the chip rate responds to global drivers only,
never to single-atom stochasticity. **The chip is not illuminated** (dose rate
zero outside the imaged field), which is the origin of the systematic
field/chip discrepancy reported in §4.6. The scheme conserves Σ k n_k exactly;
the measured mass drift over a 400 s integration is below 10⁻⁶.

The two levels share one fitted parameter: the encounter coefficient k_enc,
obtained by fitting the chip model to an eight-replica ensemble of beam-off kMC
runs (k_enc = 1.30 × 10⁻⁴; residual 7.1 % on the monomer trajectory and 8.7 %
on the rate). That single number is the whole of the coarse-graining; everything
else is shared.

### 3.5 Virtual instruments

**ADF-STEM.** Two fidelity levels, and they are never mixed silently. The
pixel-level renderer projects the 3D positions (discarding z), applies probe
blur, per-scan-line jitter, stage drift and a slowly varying window background,
and draws Poisson shot noise at the requested dose; it is sampled at 0.2 Å,
which is what single-atom imaging requires, on an 8 nm crop. The fast analytic
model emits noisy projected coordinates directly. Detection thresholds the
*noisy* integrated intensity, which reproduces both the recall and the Eddington
bias measured from pixels (§4.1).

**EELS.** A field-averaged Pt oxidation-state proxy from the M₄,₅ white-line
ratio, with noise scaling as the inverse square root of dose. Only the field
average is returned, because at single-atom concentration the per-atom core-loss
signal sits at the noise floor at tolerable dose — a limitation of the real
technique that the twin is required to respect rather than wish away.

**Intermittent 3D.** Multi-tilt ptychographic snapshots at a 60 s interval, with
lateral precision equal to the imaging precision and depth precision set by the
measured 6.6 nm depth resolution of single-projection reconstruction [R8]
converted to a standard deviation (2.80 nm). This is the only stream carrying z,
and it is tagged so that no continuous-time estimator can consume it.

**Quadrupole MS.** The true zone-integrated rate is converted to species
fluxes, convolved with the transport kernel, scaled by per-species calibration
factors, and corrupted by 1/f baseline drift, Gaussian noise, a constant
background and N₂ interference on the m/z 28 channel, then sampled at 10 Hz.

**Transport.** The impulse response is a gamma kernel with a pure dead time,
h(t) = (t−τ₀)^{k−1} e^{−(t−τ₀)/θ} / [Γ(k) θ^k] for t ≥ τ₀, so that
τ_dead = τ₀ + kθ and σ_disp = θ√k. At the reference flow τ_dead = 3.202 s and
σ_disp = 0.838 s, and both scale as the reciprocal of the flow.

### 3.6 Inversion, and why it is done on a 1 s grid

The MS record is binned to the 1 s analysis grid before inversion. The kernel's
dispersion is 0.84 s, so nothing faster is recoverable in principle; binning
also drops the measurement noise by √10 and makes the exactly
cross-validated dense solver affordable on the full record. A frequency-domain
solver is provided and is checked against the dense one at fixed regularisation,
but its generalised-cross-validation trace is only approximate for a
mirror-padded circular smoother and it systematically under-regularises, so it
produces no reported number. On a 600 s test signal the three correction levels
give root-mean-square rate errors of 0.199 (uncorrected), 0.026 (delay-only) and
0.009 (deconvolution) in units of molecules s⁻¹.

### 3.7 Reproducibility

`python studies/run_all.py` regenerates every table and figure from master seed
20260922; `python tests/test_opcem.py` runs 22 property tests covering mass
conservation, kernel moments, transfer-function recovery, agreement of two
independently coded deconvolution solvers, the provenance guard, the lead–lag
gate, detrending, the projection identity, two-stage least squares under
confounding, and the shared actuation budget. Every output file records the git
commit, the seed and the simulated-data banner.

---

## 4. Results

All numbers below are simulator outputs. Tables are in `results/tables/`,
figures in `results/figures/`, and each file carries its git commit and seed.

### 4.1 The perception front end, and a bias that matters (S0)

The pixel-level detector and the fast analytic observation model agree over a
two-decade dose series (Table S0). At the reference dose of
10³ e⁻ Å⁻² s⁻¹ and 5 frames s⁻¹ (200 e⁻ Å⁻² per frame) the pixel pipeline gives
monomer recall 0.576, recall 1.000 for entities of two or more atoms, precision
0.805, localisation RMSE 24.6 pm and an amplitude bias of +0.146 atoms; the
surrogate gives 0.491, 1.000, 0.903, 25.0 pm and +0.156. Agreement degrades away
from the reference dose — at 3 × 10² e⁻ Å⁻² s⁻¹ the surrogate under-predicts
monomer recall (0.053 versus 0.120) and the localisation error diverges — so the
surrogate is used only inside the calibrated range, and the comparison is
reported rather than assumed.

Two features of the detector are physically consequential rather than
incidental.

**Detection is size-dependent.** At the reference dose a monomer sits almost
exactly at the detection threshold (single-atom contrast-to-noise 3.52 against a
threshold of 3.54) while a dimer, with twice the integrated signal, is always
detected. So the detector under-counts monomers and not dimers, and therefore
**inflates the apparent dimer fraction** — precisely the quantity one would use
to argue that dimers are the active motif.

**The amplitudes of marginal atoms are biased upward.** Thresholding the *noisy*
integrated intensity means only upward noise excursions of a marginal atom clear
the threshold, which is Eddington bias; it reproduces quantitatively in both
models (+0.146 and +0.156 atoms at the reference dose, falling to −0.020 and
−0.000 at 10⁴ e⁻ Å⁻² s⁻¹, where nothing is marginal).

### 4.2 The transport transfer function is recoverable to milliseconds (S1)

Five timestamped 1 s tracer pulses per condition, at three flows, with 2 %
relative MS noise and a slow baseline drift (Table S1, Fig. 3):

| flow (sccm) | τ_dead true / fitted (s) | bias (s) | coverage | σ_disp true / fitted (s) | bias (s) | coverage |
|---|---|---|---|---|---|---|
| 1.0 | 6.404 / 6.400 | −0.0040 | 0.80 | 1.677 / 1.673 | −0.0038 | 1.00 |
| 2.0 | 3.202 / 3.203 | +0.0009 | 1.00 | 0.8385 / 0.8386 | +0.0001 | 1.00 |
| 4.0 | 1.601 / 1.600 | −0.0006 | 1.00 | 0.4192 / 0.4190 | −0.0002 | 1.00 |

The bias is at the millisecond level against a multi-second dead time, and the
bootstrap intervals cover the truth at close to their nominal rate (13 of 15
for τ_dead, 15 of 15 for σ_disp). **This step is cheap, and it is the
precondition for everything that follows.**

### 4.3 Descriptors: accurate in the aggregate, biased where it counts (S2)

Over five independent 600 s runs at the reference dose (Table S2a), correlations
with truth are high for every census descriptor (Pearson r from 0.962 to 0.999),
but the biases are not small and not uniform:

| descriptor | true | observed | relative bias |
|---|---|---|---|
| N₁ | 37.7 | 22.2 | **−41.2 %** |
| N₂ | 3.88 | 4.39 | **+13.4 %** |
| N₃ | 3.93 | 3.92 | −0.05 % |
| N₄₊ | 6.40 | 6.41 | +0.08 % |
| dispersed fraction | 0.706 | 0.589 | **−16.6 %** |
| k_hop (s⁻¹) | 2.672 | 1.861 | −30.3 % |
| τ_dimer (s) | 3.436 | 3.015 | −12.3 % |

This is the size-dependent detection bias of §4.1 propagating exactly as
predicted: monomers are lost, dimers are not, the dispersed fraction is
understated by a sixth, and the apparent dimer fraction is correspondingly
inflated. **Any experimental argument that small clusters are the active motif
must report its detector's size-resolved recall**, and a dose series is how to
obtain it.

The dynamic descriptors are recovered to −12 % (dimer lifetime) and −30 %
(hop rate); the hop-rate bias is not a failure of the mean-squared-displacement
inversion but a consequence of comparing against the wrong reference, as §4.5
establishes.

**The central hypothesis fails its first test.** With reactor variables forced
into every model, variance explained in the instantaneous true rate is
(Table S2b): census only, R² = 0.9860 ± 0.0027 (7 features); census plus six
dynamic descriptors, R² = 0.9866 ± 0.0028 (13 features); a dimension-matched
census subset, R² = 0.9861 (6 features); dynamic descriptors alone,
R² = 0.614 ± 0.129. **Adding rate- and lifetime-type descriptors buys
0.0006 in R².** In hindsight this is close to a tautology: where the rate is a
sum of per-motif turnovers, the instantaneous rate is a linear functional of the
instantaneous census, and nothing dynamic can improve on it.

Where the dynamic descriptors do pay is prediction (Table S2c). Explaining the
*future change* r(t+h) − r(t):

| horizon | census only | census + dynamic | gain |
|---|---|---|---|
| 30 s | 0.386 | 0.387 | +0.001 |
| 60 s | 0.532 | 0.545 | +0.013 |
| 120 s | 0.632 | 0.659 | **+0.027** |

Small, but in the direction the hypothesis predicts and growing with horizon.
The defensible claim is therefore about forecasting deactivation, not about
explaining instantaneous activity.

The dynamic descriptors are also sensitive to the size decoder, which must be
pre-registered (Table S2d): a temporal median of width 3 recovers 54 % of true
nucleation events and τ_dimer to −12 %; width 5 recovers 26 % of events but
τ_dimer to +9 %; a Viterbi decoder recovers 14 % of events and τ_dimer to +20 %.
The trade-off between event recall and lifetime bias is intrinsic — over-smoothing
discards short-lived states — and reporting one setting without the others would
be a free parameter in disguise.

### 4.4 The event-aligned estimator, and the Pt loading it demands (S3)

The estimator is scored against three references: the exact instantaneous
single-event effect Δ₀ = TOF(2) − 2·TOF(1) = 0.0647 molecules s⁻¹; the same
effect averaged over the 4 s post-window given the dimer lifetime
τ = 3.44 s, Δ_W = Δ₀(τ/W)(1 − e^{−W/τ}) = 0.0382; and Δ_true, what the same
estimator recovers from the noise-free true rate, which is the ceiling for any
instrument. Four runs of 1200 s per loading (Table S3a, Fig. 4b):

| Pt atoms / field | Pt nm⁻² | λ (s⁻¹) | λ(t_pre+t_post) | Δ_true/Δ_W | Δ_dec/Δ_true | Δ_dec/Δ_W | p (surrogate) |
|---|---|---|---|---|---|---|---|
| 6 | 0.0038 | 0.010 | 0.082 | 0.72 | **0.85** | 0.57 | 0.003 |
| 8 | 0.0050 | 0.041 | 0.33 | 0.80 | **0.76** | 0.58 | 0.001 |
| 16 | 0.010 | 0.138 | 1.10 | 0.60 | 0.65 | 0.38 | 0.004 |
| 32 | 0.020 | 0.314 | 2.51 | 0.48 | 0.49 | 0.24 | 0.050 |
| 90 | 0.056 | 0.626 | 5.01 | 0.46 | 0.48 | 0.22 | 0.19 |

Two things follow. First, **the instrument is not the bottleneck**: in the
sparse regime the mass spectrometer plus the transport inversion cost only
15–24 % of the effect (Δ_dec/Δ_true = 0.85 and 0.76). Second, **the event
density is**. Recovery falls monotonically with λ(t_pre+t_post) and the
surrogate test stops rejecting by λ(t_pre+t_post) ≈ 5, where the windows of
neighbouring events overlap and no single event's contribution survives. Since
λ scales roughly as the square of the monomer count, this is a **maximum Pt
loading**: about 0.005 Pt nm⁻² here, roughly an order of magnitude below the
~0.03 Pt nm⁻² at which mainly-single-atom Pt/Al₂O₃ is reported [R11].

The controls behave correctly (Table S3b, 8 atoms per field). Nucleation events
give Δ = +0.0223 (95 % CI 0.0083 to 0.0360; p = 0.001; z = 5.9 against the
circular-shift null). Dissociation events — the sign control — give
Δ = −0.0223 (ratio −0.58), a clean reversal. Using the events the pipeline
actually detects rather than the simulator's own reduces the effect to +0.0130
(ratio 0.34), consistent with 54 % event recall diluting the average with
non-events. Omitting the detrending reduces it to +0.0146 (ratio 0.38), because
progressive sintering drags the post-window down.

Power is modest (Table S3c): 0.34 at 5 events, 0.48 at 10, 0.51 at 20 and 0.68
at 40. Reaching 80 % power needs more than 40 resolvable events, which at
0.041 s⁻¹ means runs longer than ~1000 s — and the run has to stay at a
loading low enough to keep λ small, so **event count and event resolvability
pull against each other**. That is the practical envelope of the method.

### 4.5 What skipping the transport correction costs (S4)

Identical data, three correction levels (Table S4a, Fig. 5):

| correction | rate RMSE (molec s⁻¹) | Δ recovered | Δ/Δ_W | p (surrogate) |
|---|---|---|---|---|
| none | 0.0175 | 0.0002 | **0.005** | 0.35 |
| delay-only shift | 0.0062 | 0.0110 | 0.29 | 0.047 |
| deconvolution | 0.0053 | 0.0192 | **0.50** | 0.011 |

Uncorrected, the effect is gone — recovery 0.5 %, and the surrogate test does not
reject. A delay-only shift recovers 29 %; full deconvolution recovers 50 %. The
lead–lag result is starker (Table S4b): **without correction the peak sits at
3.0 s in every run** — the transport dead time is 3.20 s — and the estimator
would be called significant in 4 of 4 runs. After either correction the peak
moves to 0 s and the resolution gate correctly declines to call it significant
in 0 of 4. An uncorrected analysis would therefore report that the structure
leads the chemistry by three seconds, when the true lag is zero and the three
seconds are plumbing.

Under deliberate confounding by a slow temperature modulation, which moves both
structure and rate (Table S4c), the naive regression of rate on dimer count
gives β = 0.201 ± 0.013, the elastic net with reactor variables partialled out
gives 0.433, and two-stage least squares instrumented by a randomised O₂ pulse
gives 0.373 ± 0.294. The naive estimate is biased low by about a factor of two.
But the instrument is weak — first-stage R² = 0.009 — so the interventional
estimate, while unbiased in expectation, is imprecise. **The design requirement
is therefore that the randomised actuator must move the structural variable
strongly**; a first-stage R² of one per cent is not an instrument, and reporting
a 2SLS coefficient from one would be reporting noise.

### 4.6 Projected mobility is biased, and the bias is not correctable at 6.6 nm depth resolution (S5)

Against the simulator's own 3D coordinates, with minimum-image displacements
(Table S5):

| corrugation A (nm) | η analytic | η from MSD truth | η from 3D snapshots | projection bias | error after correction |
|---|---|---|---|---|---|
| 0 | 0 | 0.000 | 0.218 ± 0.019 | 0 % | +37 % |
| 0.35 | 0.048 | 0.046 | 0.220 ± 0.018 | −4.4 % | +32 % |
| 0.70 | 0.193 | 0.184 | 0.224 ± 0.018 | −15.5 % | +17 % |
| 1.05 | 0.435 | 0.413 | 0.230 ± 0.017 | −29.2 % | −1.8 % |

The identity η = A²q²/4 is confirmed by the simulator to better than 6 %, and
the projection bias is real and substantial, reaching −29 % at the largest
corrugation. **The correction fails.** The snapshot-based estimator returns
η ≈ 0.22 regardless of the truth, including when the truth is exactly zero,
because the reconstruction's depth precision (σ_z = 2.80 nm, from the measured
6.6 nm depth resolution of single-projection ptychography [R8]) is several times
the corrugation being measured. The apparent success at A = 1.05 nm is a
coincidence of the estimator's fixed output crossing the truth.

This is an instrument requirement, not an analysis choice: correcting projected
mobility needs sub-nanometre depth resolution, which tilt-coupled multislice
ptychography now provides [R9]. Until then the honest course is to report
projected mobility as projected, with the bias bounded by an independent estimate
of the support's roughness.

### 4.7 Beam accounting: a usable dose window of one decade, and a misspecified model (S6)

The dose series first has to be cleaned. Below ~10³ e⁻ Å⁻² s⁻¹ the matched
filter runs near its noise floor, false positives outnumber real atoms across a
40 nm field, and chains of coincidentally linked false positives dominate the
mean squared displacement: at 10² e⁻ Å⁻² s⁻¹ the apparent hop rate is 9.8 s⁻¹
against a true effective rate of 1.17 s⁻¹. A diagnostic that needs no ground
truth — the fraction of detections belonging to a track of at least five frames —
identifies these doses cleanly (0.026, 0.254 and 0.460 at 10², 3 × 10² and 10³
e⁻ Å⁻² s⁻¹, rising above 0.93 from 2 × 10³ upward), and the extrapolation is
restricted to the doses that pass a pre-registered 0.80 gate. Since beam
perturbation bounds the window from above, **the accessible range is about one
decade**, and that is what limits the precision of any zero-dose extrapolation.

The reference also has to be right. A monomer re-draws its site type after each
hop, so it occupies site type s in proportion to p_s/k_s and the ensemble hop
rate a mean-squared-displacement estimator measures is the harmonic mean
1/Σ(p_s/k_s), not the terrace rate. Against that reference the observed hop
rate is accurate to −8 % to −20 % across the usable window — the −30 % "bias" of
§4.3 was mostly the wrong comparison.

Because a harmonic mean of (k_s + βφ) is not affine in φ, the linear-in-dose
model that the Beam Perturbation Index assumes is misspecified whenever the
support presents a distribution of barriers. The consequence is a trade-off
(Table S6b):

| model | zero-dose value | reference | relative error | β fitted | β true |
|---|---|---|---|---|---|
| linear in dose | 1.713 ± 0.069 | 1.069 | **+60 %** | 1.62 × 10⁻⁴ | 2.0 × 10⁻⁴ |
| site-aware harmonic mean | 1.169 | 1.069 | **+9.3 %** | 1.20 × 10⁻⁴ | 2.0 × 10⁻⁴ |

Yet the index itself behaves the other way round (Table S6c). The linear model's
BPI tracks the analytic truth closely (0.159 vs 0.187 at 2 × 10³; 0.486 vs
0.488 at 10⁴; 0.654 vs 0.647 at 2 × 10⁴), because the errors in intercept and
slope partly cancel in the ratio βφ/(k_chem + βφ), while the site-aware model
under-estimates it (0.130 vs 0.187, 0.373 vs 0.488). The analytic expression
Σ p_s βφ/(k_th,s + βφ) matches the simulator's own beam-attribution counter to
better than 2 % at every dose, which validates the accounting itself. **Report
both models**: the linear one for the index, the site-aware one for the
beam-free rate.

At the reference dose the true beam-attributed fraction of hops is 0.112,
comfortably below the pre-registered exclusion threshold of 0.20; at
6 × 10³ e⁻ Å⁻² s⁻¹ it is 0.373 and events would be excluded from causal claims.
The dispersed fraction extrapolates to 0.722 ± 0.010 at zero dose from an
observed 0.617 at the reference dose — a 17 % beam-induced suppression that would
be reported as a property of the catalyst if the series were not run.

### 4.8 The finding that changes the experiment: representativeness and the flow trade-off (S7)

**The imaged field cannot be measured.** One 40 × 40 nm field holds 90 Pt atoms
in this configuration and generates about 1.4 CO₂ molecules s⁻¹. A realistic
chip — a ~100 nm catalyst film over a 200 × 50 µm reactive zone, equivalent to a
few micrograms of a ~50 m² g⁻¹ support — presents of order 3.5 × 10⁴ times the
projected area, giving 2 × 10¹³ Pt atoms and a representativeness fraction
f_rep ≈ 4.6 × 10⁻¹². Against a quadrupole's flux detection limit at 1 ppm and
2 sccm (9.0 × 10¹¹ molecules s⁻¹), the imaged field's flux is short by a factor
of 6.6 × 10¹¹ (Table S7a). **Spontaneous single-atom event analysis against a
product detector is not difficult; it is excluded**, and no improvement in
microscopy changes it, because the limitation is on the product side. The
conclusion is insensitive to the one estimated quantity in it: the field's
shortfall does not depend on the catalyst-area multiplier at all, and the chip's
signal moves proportionally over the four decades of it that we tabulate
(Table S7b).

**The two escapes conflict.** A quadrupole measures a mole fraction, so its flux
limit falls with flow; the transport dead time and dispersion rise as the
reciprocal of flow. Both constraints, on one axis (Table S7c, Fig. 7a):

| flow (sccm) | chip signal / limit | detectable (≥3×) | τ_dead (s) | σ_disp (s) | smallest resolvable lag (s) |
|---|---|---|---|---|---|
| 0.05 | 13.3 | yes | 128.1 | 33.5 | **67.1** |
| 0.10 | 6.6 | yes | 64.0 | 16.8 | **33.5** |
| 0.20 | 3.3 | yes | 32.0 | 8.4 | **16.8** |
| 0.50 | 1.3 | no | 12.8 | 3.4 | 6.7 |
| 1.0 | 0.66 | no | 6.4 | 1.7 | 3.4 |
| 2.0 | 0.33 | no | 3.2 | 0.84 | 1.7 |
| 4.0 | 0.17 | no | 1.6 | 0.42 | 0.84 |

The chip's signal clears three times the detection limit only at or below
0.2 sccm, where the smallest lag the pipeline will call significant is 16.8 s,
rising to 67.1 s at 0.05 sccm. **The product-linked time resolution of a
nanoreactor operando experiment is tens of seconds, not the sub-second scale of
the atomic events being imaged.** We have not found this trade-off stated in the
literature reviewed in `docs/01_methods_review.md`, and it bounds what any
structure–activity correlation from such an experiment can mean.

**The illuminated field is not the chip, and the beam is why.** The chip is not
illuminated. Comparing the imaged field's dispersed fraction with the chip mean
field's (Table S7d):

| dose rate (e⁻ Å⁻² s⁻¹) | field | chip | relative difference | Pearson r |
|---|---|---|---|---|
| 0 | 0.693 ± 0.057 | 0.836 | −17.1 % | 0.835 |
| 10³ | 0.524 ± 0.036 | 0.836 | −37.4 % | 0.903 |
| 10⁴ | 0.216 ± 0.020 | 0.836 | **−74.1 %** | 0.936 |

Even at zero dose the small illuminated patch differs from the deterministic
chip by 17 % — a finite-size effect — and the beam drives the difference to
−74 % at 10⁴ e⁻ Å⁻² s⁻¹. This is a discrepancy in the *ground truth*,
independent of any detection bias, and it means the structure the microscope
reports is systematically more aggregated than the structure making the product.

**Position matters.** With three pre-registered windows along the flow axis,
each carrying its own local temperature and product enrichment (Table S7e), the
recovered effect ranges over a factor of 2.8: inlet 0.681, centre 0.478,
outlet 0.246 of the window-averaged truth. The spread is 97 % of the mean and
the heterogeneity is significant (one-way ANOVA F = 25.0, p = 2.1 × 10⁻⁴,
4 runs per position, Table S7f). **One field of view, at one position, is not a
measurement of the catalyst**, and a position-stratified design is the minimum
defence.

**What survives: randomised global perturbation.** A global O₂ step drives the
same transition in every patch at once, so the chip-averaged rate responds even
at f_rep ≈ 10⁻¹². Two design requirements emerge (Table S7g). First, the
perturbation must be randomised: under a circular-shift null a strictly periodic
pulse train realigns with itself after one period, the null contains values as
large as the observed effect, and nothing can be concluded. That is the correct
answer — with a purely periodic driver, temporal coincidence carries almost no
information — and it is a caution about reading causation from oscillatory
coincidence generally. Second, the direct coverage-mediated response must be
separated from the structure-mediated one; raising p_O₂ changes the coverages and
therefore the rate immediately, with no involvement of the catalyst's structure.
Blanking a gap of the pulse duration plus three transport dispersions and
measuring what survives isolates the structural component.

### 4.9 Closed-loop control: the actuation helps, the timing does not (S8)

All policies spend an identical budget of 24 oxygen pulses, so the comparison
isolates *when* the budget is spent. Averaged over the seeds, deactivation is
substantial in this system: the open-loop rate falls by roughly a factor of
three over 900 s as Pt accumulates in clusters of four or more atoms, and both
actuated policies slow that decline.

The headline, from the held-out comparison (Table S8b), is a null. The trigger
threshold was chosen on one set of simulator seeds and the paired comparison
made on a disjoint set. In sample the event-triggered policy appeared to give
about +7 % integrated yield over the budget-matched fixed schedule; on held-out
seeds the difference is negligible and far from significant. **Event triggering
buys nothing over a fixed schedule once its threshold is chosen honestly.** The
in-sample value is reported alongside so that the selection optimism is visible,
and the minimum detectable effect is reported so that the null can be read as an
exclusion of a large gain rather than merely as a failure to reject.

The reason is visible in the trajectories (Fig. 8b): deactivation here is a
smooth, near-monotonic drift, and a fixed schedule is already close to optimal
against a smooth drift. Event triggering should pay only when deactivation is
bursty, so that there is something to anticipate — a testable prediction, and the
natural next step for the control layer rather than a claim we can make now.

## 5. Discussion

### 5.1 What the validation actually establishes

The estimators divide cleanly into three groups.

**Works, and works well.** Transfer-function calibration from timestamped
tracer pulses is close to exact: the bias is milliseconds against a
three-second dead time, and the bootstrap intervals have nominal coverage. The
transport inversion then costs only 15–25 % of an event-aligned effect in the
regime where events are resolvable at all, which means the mass spectrometer and
the transport line are *not* the limiting elements once the transfer function is
measured. Detrending behaves as intended, removing a linear drift by more than
an order of magnitude while leaving an injected step intact. The lead–lag
resolution gate does its job: the uncorrected trace produces a spurious peak at
the transport dead time in every run, and the gate refuses it in every run.

**Works, with a caveat that changes the experiment.** Event-aligned analysis
recovers the known single-event effect, but only when the dimensionless event
density λ(t_pre + t_post) is well below one. Since the nucleation rate scales
roughly as the square of the monomer count, that is a constraint on Pt loading —
and it is a tight one, roughly an order of magnitude below the surface densities
at which atomically dispersed Pt is usually prepared. Beam accounting works, but
only inside a dose window bounded from below by false-positive-dominated
tracking and from above by the beam's own perturbation; here the window is about
one decade, and the standard linear-in-dose model is misspecified whenever the
support presents a distribution of diffusion barriers, because the observable
effective hop rate is then a harmonic mean and sublinear in dose.

**Does not work as proposed.** Two things fail.

The projection-aware mobility correction fails, and it fails for a reason that
is a property of the instrument rather than of the estimator: at 6.6 nm depth
resolution the reconstruction noise is several times the support corrugation
being measured, so the intermittent 3D stream returns essentially the same
anisotropy whatever the truth. The projection bias is real — it runs from −4 %
to −29 % across the corrugations tested — and it is simply not correctable from
single-projection snapshots. Correcting it requires sub-nanometre depth
resolution, which is exactly what tilt-coupled multislice ptychography now
delivers [R9]; that is a concrete instrument requirement rather than an
analysis choice.

Event-triggered control gives no advantage over a budget-matched fixed schedule
once the trigger threshold is chosen out of sample. In sample it appeared to
give about +7 %; on held-out seeds the difference vanished. Both actuated
policies beat the un-actuated baseline, so the actuation itself helps; the
*timing intelligence* does not. The reason is visible in the trajectories:
deactivation in this system is a smooth, near-monotonic drift, and a fixed
schedule is already close to optimal against a smooth drift. Event triggering
should be expected to pay only when deactivation is bursty, so that there is
something to anticipate. That is a testable prediction, and the natural next
step for the control layer is to test it against a simulator with an explicit
avalanche channel rather than to claim the gain here.

### 5.2 The finding that changes the experimental design

The representativeness and sensitivity budget is the most consequential result,
and it is arithmetic rather than statistics.

One 40 × 40 nm field of view holds of order 10²–10³ Pt atoms. A realistic chip,
with a ~100 nm catalyst film over a 200 × 50 µm reactive zone, holds of order
10¹³. The representativeness fraction is therefore about 5 × 10⁻¹², and the
imaged atoms generate a product flux around ten orders of magnitude below the
flux detection limit of a quadrupole mass spectrometer. **Spontaneous
single-atom event analysis against a product detector is not difficult; it is
excluded.** No improvement in microscopy changes this, because the limitation
is on the product side.

Worse, the two obvious escapes conflict. Lowering the flow lowers the
detector's flux limit, because a quadrupole measures a mole fraction; but the
transport dead time and dispersion rise as the reciprocal of the flow. In our
configuration the chip's signal clears three times the detection limit only at
flows at or below 0.2 sccm, where the smallest lag the pipeline will call
significant is about 17 s, rising to 67 s at 0.05 sccm. **The product-linked
time resolution of an operando nanoreactor experiment is tens of seconds, not
the sub-second scale of the atomic events being imaged.** This trade-off appears
not to be stated anywhere in the literature we reviewed, and it bounds what any
"structure–activity correlation" from such an experiment can mean.

What survives is perturbation-synchronised analysis. A global actuator step
drives the same transition in every patch of the chip simultaneously, so the
chip-averaged signal responds even at f_rep ≈ 10⁻¹², and the imaged field
becomes a *sampler of a synchronised response* rather than the source of the
signal. Two design requirements follow, and both are demonstrated here. The
perturbation must be **randomised, not periodic**: under a circular-shift null a
strictly periodic pulse train realigns with itself and carries essentially no
inferential power, which is the correct answer and is also a caution about
reading causation from oscillatory coincidence. And the direct
(coverage-mediated) response must be separated from the structure-mediated one,
which we do by blanking a gap of the actuation duration plus several transport
dispersions and measuring what survives afterwards, when the gas composition has
returned to baseline.

### 5.3 The hypothesis, partly refuted

The central hypothesis — that activity is governed by transition rates and
lifetimes rather than by the instantaneous motif census — fails its
pre-specified first test. Census descriptors explain about 99 % of the variance
in the instantaneous rate, and adding six rate- and lifetime-type descriptors
adds essentially nothing. This is not a surprise in hindsight: in any model
where the rate is a sum of per-motif turnovers, the instantaneous rate is a
linear functional of the instantaneous census, and no dynamic quantity can
improve on that. The hypothesis was, in that form, close to unfalsifiable in the
wrong direction — it could only have been supported by a system in which the
census is *unobservable*, which is a statement about measurement rather than
about catalysis.

What the dynamic descriptors do buy is forecast skill: they add measurably to
the prediction of how the rate will change over the next one to two minutes,
and not at all over thirty seconds. That is the honest version of the claim, and
it is the version relevant to the engineering goal — predicting deactivation
before it happens — rather than to the mechanistic one.

A separate methodological point emerged from the same analysis. The observed
dispersed fraction is biased low by about 17 % at the reference dose, because
detection is size-dependent: a monomer sits near the detection threshold while a
dimer is always detected. The bias therefore *inflates the apparent dimer
fraction*, which is precisely the quantity one would use to argue that dimers
are the active motif. Any experimental claim of that form has to report the
size-resolved recall of its detector, and a dose series is the way to get it.

### 5.4 Limitations

The limitations are of three kinds, and the first is the most serious.

**The twin could be wrong in ways that would invalidate the validation.**
Everything reported here is conditional on the forward model. Specifically:
(i) the kinetics are a tau-leaped continuous-time Markov model with one hop
length and three site types, so a real support's continuous barrier
distribution, correlated site environments and adsorbate-dependent barriers are
absent — and §4.5 already shows that barrier heterogeneity alone breaks the
linear dose model, so a richer distribution would break it further; (ii) the
beam enters only as an additive term on the monomer hop rate, whereas a real
beam also reduces the support, creates and annihilates vacancies, cracks
adsorbates and ionises the gas, none of which is modelled; (iii) turnover is
Langmuir–Hinshelwood with size-indexed parameters and no coverage-dependent
restructuring, so the structure→activity map is imposed rather than emergent;
(iv) the chip is a deterministic population balance with a one-dimensional
gradient and one fitted coefficient, and it reproduces the kMC ensemble only to
7–9 %; (v) the ADF contrast model is a scalar cross-section with Gaussian probe
blur, not a multislice calculation, so channelling, thickness and tilt effects
are absent. An estimator that passes on this twin has not been shown to pass on
a real experiment. It has been shown to pass on *a* concrete forward model whose
assumptions are written down, which is strictly more than an experimental paper
with unknown ground truth can show, and strictly less than a validated
measurement.

**Statistical power is modest throughout.** Most studies use four to six
simulator seeds. The position-heterogeneity result is significant
(p ≈ 2 × 10⁻⁴) and the control result is clearly null out of sample, but
intermediate effects — the closed-loop comparison against open loop, for
instance — are underpowered, and we report them as such rather than as
marginal evidence.

**Scope.** The single reaction is CO oxidation; ORR, CO₂ reduction and
hydrogenation are specified in the design but not simulated. Only one support
is modelled. EELS is a field-averaged proxy, not a per-atom charge state, which
is faithful to what the technique currently delivers but means the platform's
`spectro` block is thinner than the descriptor vector suggests. Electrochemical
operation, the liquid-cell variant and the 4D-STEM/ptychographic reconstruction
chain are specified but not implemented; the 3D stream is injected as noisy
coordinates rather than reconstructed from simulated diffraction.

### 5.5 What would have to be done experimentally

In order, and with the gating criterion for each step:

1. **Characterise the detector before the catalyst.** Size-resolved recall,
   precision and localisation error on a known-static reference specimen, over
   at least a decade of dose. Gate: the persistent-track fraction must exceed
   0.8 at the intended dose, or mobility estimates are meaningless.
2. **Measure the transfer function at every set point.** Gate: the acceptance
   test of §3.1 — a deliberate valve step must appear in the reactor log and the
   MS channel separated by the independently measured dead time to within the
   measured dispersion.
3. **Establish the flow operating point from the trade-off, not from habit.**
   Gate: chip signal at least three times the detection limit, and the resulting
   resolvable lag reported as the experiment's true time resolution.
4. **Set the Pt loading from the event-density criterion**, if spontaneous-event
   analysis is intended at all — and note that step 3 will usually show it is
   not available, in which case go to step 6.
5. **Run the dose series** and report a Beam Perturbation Index under both the
   linear and the site-aware model, with the barrier distribution taken from a
   beam-off temperature series.
6. **Use randomised global perturbations** and the gap-blanked estimator to
   separate coverage-mediated from structure-mediated responses.
7. **Cross-validate the inferred motif** against operando XAS, DRIFTS or Raman
   on the same batch outside the microscope, and against DFT and microkinetic
   modelling. The microscope establishes *which structures exist and when*; it
   cannot establish a mechanism on its own.

Only after step 7 would a claim of the form "transient Pt dimers are the active
motif under these conditions" be supportable — and it would be a statement about
a chip-averaged, perturbation-synchronised response with a stated time
resolution of tens of seconds, not about a single atom seen to move.

---

## 6. Conclusion

Seeing single atoms move is not the same as knowing which one made the product,
and the gap between the two is not a resolution problem. We have formalised it
as event-driven causal inference, specified the platform it requires, and
measured — against known ground truth — what each estimator in that pipeline can
and cannot do.

Three results should change how such experiments are designed. The transport
transfer function must be measured, not assumed: uncorrected, the analysis
places the apparent lead–lag peak at the transport dead time and would read the
structure as leading the chemistry by that amount. The imaged field's product
flux is roughly ten orders of magnitude below any product detector's limit, so
spontaneous single-atom event analysis is excluded and randomised global
perturbation is the only route that identifies anything. And the flow that makes
the product detectable makes the product-linked time resolution tens of seconds,
which bounds what any structure–activity correlation from a nanoreactor
experiment can mean.

Two of our own proposals did not survive their tests. Dynamic descriptors do not
out-explain a motif census for the instantaneous rate, though they do add
forecast skill over one to two minutes. Event-triggered control gave no
advantage over a budget-matched fixed schedule once its threshold was chosen out
of sample. We report both, because a protocol whose failure modes are unknown is
not a protocol.

What we offer is therefore not a discovery about platinum but a pre-registrable
analysis protocol with its bias, power and operating limits measured in advance,
and open code that regenerates every number in this paper from a single seed.

---

## Data and code availability

All code, configuration, tables, figures and raw outputs are in the
accompanying repository. `python studies/run_all.py` regenerates every number in
this manuscript from master seed 20260922; `python tests/test_opcem.py` runs the
22 property tests described in §3.7. Each output file records the git commit,
the seed and the configuration that produced it, together with a
machine-readable statement that its contents are simulated.

There is no experimental dataset, because no experiment was performed.

## Research integrity declarations

**Nature of the data.** Every numerical value in this manuscript is an output of
the forward simulator specified in §3.4 and implemented in `src/opcem/truth.py`,
observed through the virtual instruments of §3.5. No physical electron
microscope, MEMS nanoreactor, mass spectrometer, catalyst or gas-handling system
was operated. No measurement is reported, and no simulated value is presented as
a measurement anywhere in the text, tables or figures. The rendered micrograph
in Fig. 1a carries a visible "SIMULATED" overlay because it is the one panel
that could be mistaken for data if separated from its caption.

**Use of references.** Published work is cited for two purposes only: to
establish what is experimentally demonstrated today (§1, `docs/01`), and to
justify the numerical value of a simulation parameter (§3.4–§3.6,
`docs/04` §5.4). No reference is cited in support of any result produced by the
simulator. The reference list marks which entries were used for which purpose.

**Pre-specification.** The three falsifiers of the central hypothesis (§2.1),
the analysis thresholds (`src/opcem/config.py`), the event-window lengths, the
Beam Perturbation Index exclusion threshold, the persistent-track gate and the
position-stratified design were fixed before the studies reporting them were
run. Two of the three falsifiers were triggered, and both are reported: the
hypothesis fails for the instantaneous rate (§4.3) and the closed-loop policy
fails against its budget-matched baseline (§4.9). Where a parameter was tuned
after seeing results — the event-window length, and the dilute loading used in
§4.4 onwards — the tuning and its reason are stated in the text.

**Selection.** In §4.9 the trigger threshold is selected on one half of the
simulator seeds and evaluated on the disjoint other half; the in-sample value is
reported next to the held-out value so that the selection optimism is visible.
No other result involves a selected hyper-parameter.

**Negative and null results.** Four are reported: dynamic descriptors do not
out-explain a motif census for the instantaneous rate (§4.3); the
projection-aware mobility correction fails at current depth resolution (§4.6);
the linear-in-dose model underlying the Beam Perturbation Index is misspecified
for a heterogeneous support (§4.7); and event-triggered control does not beat a
budget-matched schedule out of sample (§4.9). None was removed or reframed as a
positive result.

**Authorship and assistance.** The work was conceived and directed by the
author. Software implementation, literature retrieval and manuscript drafting
were carried out with the assistance of an AI coding assistant; the author is
responsible for the scientific content, for the design of every study, and for
the interpretation of every result. All literature claims were verified against
publisher or repository records at the time of writing.

**Competing interests.** None declared.

**Funding.** None.

## Self-assessment

A referee-style critical analysis of this manuscript, written against it rather
than for it, is provided as `docs/05_paper_analysis.md`. It identifies the
paper's central weakness (the validation is circular wherever it concerns the
structure–activity map the simulator imposes), the results whose prose is
stronger than their tables (forecast skill; the closed-loop null read as an
exclusion), the single most load-bearing estimated quantity (the catalyst-area
multiplier), and eight prioritised improvements. Five of the eight were acted on
before this version: the simulated-data declaration was moved into the first
sentence of the abstract; the seed count was raised for the two underpowered
studies; the representativeness fraction is now reported as a sensitivity band
over its estimated multiplier; a strong-instrument arm was added to the
interventional comparison; and the rendered micrograph panel was watermarked.
Three remain open and are stated as such: adding a bursty-deactivation mode to
the twin in order to scope the control null, splitting the manuscript for
publication, and dropping the coinage.

## Author contributions

M.N.A.: conceptualisation, methodology, formal analysis, software
specification, supervision, writing (original draft and revision).

## Acknowledgements

The author thanks the Department of Electronics and Communication Engineering,
Velalar College of Engineering and Technology, for institutional support.

## References

Entries are marked **[E]** where cited for experimental precedent and **[P]**
where cited for a simulation parameter value. No entry is cited in support of a
simulated result.

R1. **[E]** S. B. Vendelbo, C. F. Elkjær, H. Falsig, I. Puspitasari, P. Dona,
L. Mele, B. Morana, B. J. Nelissen, R. van Rijn, J. F. Creemer, P. J. Kooyman,
S. Helveg. Visualization of oscillatory behaviour of Pt nanoparticles catalysing
CO oxidation. *Nature Materials* **13**, 884–890 (2014).
doi:10.1038/nmat4033

R2. **[E]** Insights into chemical dynamics and their impact on the reactivity
of Pt nanoparticles during CO oxidation by operando TEM. *ACS Catalysis* (2019).
doi:10.1021/acscatal.9b03692

R3. **[E]** Atomic-level fluxional behavior and activity of CeO₂-supported Pt
catalysts for CO oxidation. *Nature Communications* (2021); preprint
arXiv:2104.00821.

R4. **[E][P]** Quantitative gas-phase transmission electron microscopy: EELS
measurement of catalytically produced CO₂ and its distribution in the
environmental cell. NSF public-access record 10499240. *(Used for the 12–21 %
product non-uniformity adopted in §3.4.)*

R5. **[P]** Chemical kinetics for operando electron microscopy of catalysts: 3D
modeling of gas and temperature distributions during catalytic reactions.
Preprint arXiv:2003.10426.

R6. **[P]** Temperature distributions in MEMS microheaters during gas-phase
experiments in an environmental TEM. *Ultramicroscopy* (2026).
ScienceDirect S0304399126000197.

R7. **[E]** Identify the activity origin of Pt single-atom catalyst via
atom-by-atom counting. *Journal of the American Chemical Society* (2021).
doi:10.1021/jacs.1c06381

R8. **[P]** Three-dimensional structural and compositional inhomogeneity in
zeolites unraveled by low-dose electron ptychography. Preprint arXiv:2212.08998.
*(Source of the 6.6 nm depth resolution and ~3.5 × 10³ e⁻ Å⁻² dose adopted in
§3.5.)*

R9. **[E]** Sub-nanometre depth resolution and single-dopant visualization
achieved by tilt-coupled multislice electron ptychography. *Nature
Communications* (2025). doi:10.1038/s41467-025-56499-1

R10. **[P]** Atomically resolved imaging of radiation-sensitive metal–organic
frameworks via electron ptychography. *Nature Communications* (2025).
doi:10.1038/s41467-025-56215-z

R11. **[E]** Dynamics of single Pt atoms on alumina during CO oxidation
monitored by operando X-ray and infrared spectroscopies. *ACS Catalysis* **9**,
5752 (2019). doi:10.1021/acscatal.9b00903

R12. **[E]** Memory-dictated dynamics of single-atom Pt on CeO₂ for CO
oxidation. *Nature Communications* **14** (2023).
doi:10.1038/s41467-023-37776-3

R13. **[E]** Blocking the operando formation of single-atom spectators by
interfacial engineering (2025). PMID 40178203.

R14. **[E]** In situ visualization and mechanistic understanding of
facet-dependent atomic redispersion of platinum on CeO₂. *Nano Letters* **23**,
11999 (2023). doi:10.1021/acs.nanolett.3c04008

R15. **[E][P]** Quantitative electron beam–single atom interactions enabled by
sub-20-pm precision targeting (2025). PMC12442701. *(Source of the ~10⁷
e⁻ Å⁻² s⁻¹ manipulation dose rate used to bracket the dose series in §3.4.)*

R16. **[E]** Mechanism of electron-beam manipulation of single-dopant atoms in
silicon. *Journal of Physical Chemistry C* **125**, 16041 (2021).
doi:10.1021/acs.jpcc.1c03549

R17. **[E]** Single-atom force measurements: mapping potential energy landscapes
via electron-beam-induced single-atom dynamics. Preprint arXiv:1804.03729.

R18. **[E]** Probing electron-beam-induced transformations on a single-defect
level via automated scanning transmission electron microscopy. Preprint
arXiv:2207.12882.

R19. **[E]** Accelerating domain-aware electron microscopy analysis using deep
learning models with synthetic data and image-wide confidence scoring (2025).
PMC12343287.

R20. **[E]** Bridging electron microscopy and materials analysis with an
autonomous agentic platform. *Science Advances* (2026).
doi:10.1126/sciadv.aed0583

R21. **[E]** Thinking microscopes: agentic AI and the future of electron
microscopy. *npj Computational Materials* (2026). doi:10.1038/s41524-026-02077-y

R22. **[E]** Rewards-based image analysis in microscopy. Preprint
arXiv:2502.18522.

R23. **[E]** STEMGym: benchmarking sequential decision-making under dose budgets
in autonomous electron microscopy. Preprint arXiv:2606.29592.

R24. **[E]** Recent progress of operando transmission electron microscopy in
heterogeneous catalysis. *Microstructures* (2024).
doi:10.20517/microstructures.2023.72

R25. **[E]** Versatile homebuilt gas feed and analysis system for operando TEM
of catalysts at work. *Microscopy and Microanalysis* (2020).
doi:10.1017/S143192762000015X

R26. **[E]** Periodic structural changes in Pd nanoparticles during oscillatory
CO oxidation reaction. *Nature Communications* **13** (2022).
doi:10.1038/s41467-022-33304-x

R27. **[E]** Copper catalysis at operando conditions — bridging the gap between
single-nanoparticle probing and catalyst-bed averaging. *Nature Communications*
(2020). PMC7518423.

## Figures

| Figure | Content |
|---|---|
| Fig. 1 | Simulated platform output: rendered ADF-STEM crop with ground-truth and detected positions (watermarked); population recovery against truth; the transport delay and its removal; single-atom contrast-to-noise against dose |
| Fig. 2 | Consistency of the two levels of the digital twin: population trajectories and rate, kinetic Monte Carlo ensemble versus coarse-grained chip model, with the Pt mass drift annotated |
| Fig. 3 | Transport calibration: tracer pulse and fit; true versus fitted impulse response; recovery of dead time and dispersion at three flows |
| Fig. 4 | Event-aligned estimator: event-related average against the window-averaged truth and a random-time null; recovery against event density, labelled by Pt loading; detection power |
| Fig. 5 | Cost of skipping the transport correction: effect recovery, rate reconstruction error, and the spurious lead–lag peak at the dead time |
| Fig. 6 | Beam accounting: dose series with the zero-dose intercept; Beam Perturbation Index against truth; the dose-dependence of the reported dispersed fraction |
| Fig. 7 | Representativeness: the flow trade-off between detectability and time resolution; the sensitivity budget; divergence of the illuminated field from the un-illuminated chip |
| Fig. 8 | Closed-loop control: yield against trigger threshold with budget-matched baselines; example rate trajectories |

## Tables

All tables are in `results/tables/` as CSV (with a provenance header) and as
markdown. S0: front-end validation. S1: transport calibration. S2a–S2d:
descriptor accuracy, variance explained, forecast skill, decoder sensitivity.
S3a–S3c: event-density sweep, event-aligned controls, detection power.
S4a–S4c: transport-correction level, lead–lag, observational versus
interventional. S5: projection bias. S6a–S6c: dose series, dose models, Beam
Perturbation Index. S7a–S7g: representativeness budget, f_rep sensitivity, flow
trade-off, field versus chip, position stratification, position heterogeneity,
perturbation-synchronised analysis. S8a–S8b: control policies, held-out paired
test.
