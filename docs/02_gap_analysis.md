# 02 — Gap Analysis

Each gap below is stated as (i) what is currently done, (ii) why it is insufficient,
(iii) what would close it, and (iv) how this project tests that closure. Gaps are
ordered by how much they limit the *inferential* strength of an operando-TEM claim,
not by how hard they are.

---

## G1 — Correlation is argued from visual coincidence, not from a transport-corrected estimate

**Current practice.** A structural time series from the image stream and a product time
series from the mass spectrometer are plotted on the same axis, and coincidence is read
off by eye or from a zero-lag correlation [R1, R2].

**Why insufficient.** The measured MS signal is
S_j(t) = [R_j * h](t) + ε(t), where h is the impulse response of reactor + transfer line
+ inlet. For a 1-bar MEMS nanoreactor at typical flows, the mean dead time τ_dead is of
order seconds and the dispersion σ_disp of order a second. Structural events in a
single-atom catalyst occur on 0.1–10 s scales. **The delay is therefore of the same
order as the phenomenon**, so an uncorrected comparison can invert the apparent
direction of causation. Periodic-oscillation arguments [R1] survive this because a
constant phase shift does not destroy a frequency match; *event-based* arguments do not.

**Closure.** Measure h directly with timestamped non-reactive tracer pulses at every
(T, P, flow) set point; report τ_dead ± and σ_disp ±; align or deconvolve before
correlating; and propagate the residual delay uncertainty into every reported lag, so
that no lag smaller than that uncertainty is interpreted.

**Test in this project.** Study S1 calibrates h from tracer pulses and checks the
recovered τ_dead and σ_disp against the simulator's ground truth. Study S4 shows that
without the correction the apparent lead–lag peak moves to the transport dead time
and is called significant in every run.

---

## G2 — No statistical event-level estimator exists for structure→activity coupling

**Current practice.** One movie, one trace, one qualitative statement.

**Why insufficient.** A single structural event followed by a rate change is one
sample. Single-atom systems are stochastic: beam-induced hopping is explicitly reported
as stochastic, with dozens of null events before a success [R15]. Without replication
and a null model, any single coincidence has an unquantified false-positive rate.

**Closure.** Treat each structural transition as an *event* in the survival-analysis
sense: detect it, set t = 0, extract the transport-corrected product trace in
pre/post windows, average over many events, and test against a matched null built from
time-shifted or event-time-permuted surrogates. Report an effect size with a
confidence interval and the number of independent events and chips.

**Test.** Study S3 implements event-aligned averaging with a circular-shift null and
bootstrap CIs, and reports detection power as a function of event count.

---

## G3 — "3D atom tracking during reaction" is asserted beyond what the physics allows

**Current practice.** Proposals routinely promise time-resolved 3D atomic coordinates
under reaction conditions.

**Why insufficient.** Depth resolution in a single projection is ~6.6 nm at ~3.5×10³
e⁻/Å² [R8]; sub-nm depth requires multi-tilt coupling [R9]; sub-Å depth requires a full
tilt series, which requires a static object. A mobile single atom is not static. A
platform that mixes continuously-tracked 2D projected positions with intermittently
reconstructed 3D coordinates and calls both "3D trajectories" will produce
irreproducible mobility statistics, because projected displacement systematically
*underestimates* true displacement and does so anisotropically.

**Closure.** Declare a strict data typology — *continuous 2D projected trajectories*
versus *intermittent 3D snapshots* — and provide a projection-aware estimator that
infers 3D mobility parameters from 2D observations with the depth information entering
only as a prior from the intermittent snapshots. Report the projection bias explicitly.

**Test.** Study S5 quantifies the bias of 2D-projected MSD relative to ground-truth 3D
MSD and validates the projection-aware correction.

---

## G4 — Beam perturbation is controlled qualitatively, not quantified per event

**Current practice.** "Low dose was used." Occasionally a dose series is shown.

**Why insufficient.** Beam-induced and thermal/chemical hopping are the same
observable. If the total observed hopping rate is k_obs = k_chem(T, p) + k_beam(φ),
then a claim that atom mobility drives activity is confounded unless k_beam is
separated. Dose rates used for deliberate single-atom manipulation (~10⁷ e⁻ Å⁻² s⁻¹
[R15]) are ~5 orders of magnitude above low-dose ptychography practice (~10² e⁻ Å⁻²
total [R10]), so the confound is real and its magnitude is dose-dependent.

**Closure.** A per-event **Beam Perturbation Index**,
BPI = k_beam(φ) / [k_beam(φ) + k_chem], estimated from a dose-rate series at fixed
(T, p, flow) by linear regression of observed event rate on dose rate, with the
zero-dose intercept giving k_chem. Every reported atomic event carries its BPI; events
above a pre-registered BPI threshold are excluded from causal claims.

**Test.** Study S6 runs a five-point dose series, fits the intercept, and compares the
recovered k_chem and BPI to ground truth.

---

## G5 — Spatial representativeness of the imaged region is unverified

**Current practice.** One field of view is imaged; the MS reports the whole chip.

**Why insufficient.** Transport in the cell is diffusion-dominated and the product
field is demonstrably non-uniform — reported CO₂ enrichment of the operando pellet
relative to the grid spans ~12–21 % [R4] — and microheater chips carry in-plane
temperature distributions [R6]. The imaged Pt atoms may contribute a vanishing and
unrepresentative fraction of the measured signal.

**Closure.** Report the **representativeness fraction** f_rep = (imaged active
area)/(total active area) and a **position-stratified design**: pre-registered
inlet/centre/outlet fields of view, with the event-aligned effect estimated separately
per position and a heterogeneity test across positions. Keep conversion below ~10–15 %
so that the outlet composition reflects intrinsic kinetics rather than transport
limitation.

**Test.** Study S7 runs a three-position stratified design with a transport gradient in
the simulator and reports the between-position heterogeneity of the recovered effect.

---

## G6 — Descriptors are static; catalysis is dynamic

**Current practice.** Size, mean coordination number, dispersion, single-atom density —
all *census* variables, regressed across samples [R7].

**Why insufficient.** If the active entity is a transient motif, a census taken at any
instant mixes productive and unproductive configurations, and the regression
coefficient is diluted by the fraction of time the motif exists. The relevant
quantities are *rates and lifetimes*: hopping rate k_hop, dimer lifetime τ_dimer,
nucleation rate k_nuc, inter-state flux.

**Closure.** A dynamic descriptor vector D(t) that contains residence times,
transition rates and lifetimes alongside the census variables, fitted by a
hidden-Markov / state-flux model so that the *transition* structure, not only the
occupancy, is estimated.

**Test.** Studies S2 and S4 compare the variance in rate explained by static census
descriptors alone against the static + dynamic set, on the same simulated data.

---

## G7 — The loop is open: microscopy observes, it does not steer

**Current practice.** Autonomous microscopy exists and is advancing quickly, but the
loop closes onto *image quality* or *structural targets* — find defects, position the
probe, budget dose [R18, R20, R21, R22, R23]. The reward is an imaging objective.

**Why insufficient.** The scientifically valuable loop closes onto *chemistry*: detect
a precursor of deactivation and act on the gas/potential/temperature to prevent it.
Nothing in the current literature demonstrates a microscope changing reaction
conditions in response to a detected atomic event in order to protect a measured
product rate.

**Closure.** An event-triggered controller whose observation is the descriptor vector
D(t) plus the transport-corrected rate, whose action space is the reactor's actuators
(gas composition step, pulsed oxidation, temperature, potential), and whose reward is
integrated product formation minus an actuation and dose penalty.

**Test.** Study S8 compares open-loop, fixed-schedule, and event-triggered closed-loop
policies on integrated product yield and end-of-run dispersion, on matched simulator
seeds.

---

## G8 — Causal language outruns the evidence (cross-cutting)

Granger tests, lagged correlations and regression coefficients are routinely described
as establishing that a structure "is the active site". None of these identifies a causal
effect in the presence of a common driver — and in this system there is an obvious
common driver: temperature and gas composition move both the structure and the rate.
Any honest platform must (i) distinguish association from intervention, (ii) exploit
the fact that it *can* intervene, and (iii) state the assumptions under which an
intervention-based estimate is causal.

**Closure.** Use the actuators as instruments: randomise a small perturbation
(interventional design), and estimate the effect of structure on rate from the
randomised arm. Report the association-only estimate and the interventional estimate
side by side.

**Test.** Study S4b estimates a coefficient whose true value is known from the
configuration, three ways — univariate, multivariable with controls, and two-stage
least squares instrumented by a randomised actuator — under a deliberately
confounded simulator setting, and reports the instrument's first-stage strength.

---

## Gap-to-contribution map

| Gap | Contribution | Study |
|---|---|---|
| G1 transport blur | Tracer-calibrated transfer function + delay propagation | S1, S4 |
| G2 no statistics | Event-aligned estimator + surrogate null + bootstrap | S3 |
| G3 3D overclaim | Strict 2D/3D typology + projection-aware mobility | S5 |
| G4 beam confound | Beam Perturbation Index from dose-series intercept | S6 |
| G5 representativeness | f_rep budget + flow trade-off + position-stratified design + perturbation-synchronised route | S7 |
| G6 static descriptors | Dynamic descriptor set + HMM state-flux model | S2, S4 |
| G7 open loop | Event-triggered closed-loop controller | S8 |
| G8 causal overclaim | Estimation against known coefficient targets; instrument-strength diagnostic | S4b |
