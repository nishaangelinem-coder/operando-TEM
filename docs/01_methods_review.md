# 01 — Review of Existing Methods for Atom-Resolved Operando Catalysis

**Scope.** This review covers the measurement modalities that a "structure–state–activity"
platform for single-atom Pt catalysts must integrate, and states for each one what it
delivers today and where it stops. It is deliberately organised by *capability axis*
rather than by paper, because the research gap addressed in this project is an
*integration* gap, not a resolution gap.

Reference keys below (e.g. `[R1]`) resolve to `docs/references.md`.

---

## 1. Environmental / operando electron microscopy

### 1.1 Two hardware lineages

| Lineage | Principle | Typical pressure | Strengths | Limits |
|---|---|---|---|---|
| Differentially pumped ETEM | Gas admitted into the pole-piece gap, pumped in stages | ≤ 20–50 mbar | Open geometry, EELS of the gas volume is easy, large fields of view, wide temperature range | Pressure far below technical catalysis; gas path is not a defined plug-flow reactor; conversion is hard to define |
| Closed MEMS gas cell / "nanoreactor" | Catalyst sealed between electron-transparent SiN or SiO_x windows with an integrated microheater | up to ~1 bar (and above in liquid variants) | True 1-bar operando, *unidirectional* inlet→outlet flow, direct coupling to downstream analytics | Window scattering degrades resolution/contrast; small reactive volume ⇒ low absolute conversion; thermal and concentration gradients across the channel |

The closed nanoreactor is the lineage that made *operando* (as opposed to merely
*in situ*) microscopy possible, because it supplies a defined inlet–outlet transport
path against which a conversion can be defined [R1, R2].

### 1.2 The reference operando-TEM experiment

The canonical demonstration remains the oscillatory CO oxidation study of Vendelbo
and co-workers [R1], which combined, in one instrument and on one time base:

* time-resolved HRTEM of individual Pt nanoparticles inside a MEMS nanoreactor at
  ~1 bar,
* quantitative online mass spectrometry of the gas leaving the nanoreactor,
* on-chip reaction calorimetry.

The result — periodic refacetting of the Pt particles *synchronous* with periodic
oscillation of the CO/O₂/CO₂ signals — is still the clearest existing example of a
structure–activity correlation established inside a microscope. Its importance for
the present project is methodological: it established that (i) simultaneous
imaging + product analysis is physically achievable, and (ii) the argument for
structure–activity coupling rests on *temporal* coincidence of two independently
recorded traces.

Subsequent work extended this to chemical (not only morphological) dynamics of Pt
during CO oxidation, and to Pd and other metals, with conversion measured online
[R2, R3].

### 1.3 Quantitative gas-phase EELS as a near-sample chemical probe

Crozier and co-workers established EELS of the gas volume as a *quantitative*
composition probe inside the environmental cell, and validated it against finite-element
models of the cell that couple fluid flow, conduction/radiation heat transfer,
multicomponent diffusion and reaction kinetics [R4, R5]. Two results from that work
are directly load-bearing for the present design:

* Under typical operando-TEM conditions **mass transport is diffusion-dominated and
  heat transport conduction-dominated**, with radiation becoming significant above
  ~400 °C.
* The product (CO₂) concentration is **not uniform**: the reported enrichment of the
  operando pellet relative to the TEM grid ranges from ~21 % down to ~12 % as the CO₂
  mole fraction rises from 5 % to 73 % [R4].

The second point is the quantitative statement of the *spatial representativeness*
problem: the molecules the microscope is looking at and the molecules the mass
spectrometer counts are drawn from partially different populations.

### 1.4 Local temperature

MEMS microheaters do not deliver an isothermal reaction zone. Finite-element and
calibration studies of microheater chips show measurable in-plane temperature
distributions, and dedicated local-thermometry methods for nanoreactors have been
developed for exactly this reason [R5, R6]. In the open-cell geometry studied in
[R4] the grid temperature tracked the set point to better than ~2 °C, but this is a
property of that geometry and is not transferable to a closed chip without
re-calibration.

---

## 2. Atomic-scale structural metrology

### 2.1 Two-dimensional single-atom imaging

Aberration-corrected ADF-STEM routinely resolves isolated heavy atoms on light
supports, and *statistical* atom counting is now a mature quantitative method:
electron-microscopy-based atom-recognition statistics counted >18 000 Pt atoms on an
industrial Pt/Al₂O₃ reforming catalyst and correlated aromatics production
quantitatively with the **density of Pt single atoms**, with clusters contributing no
direct activity [R7]. This is the strongest existing example of an atom-counted
descriptor being regressed against a rate.

Its limitation is equally instructive: the counting is *ex situ* and *static*. It
yields a population census, not a trajectory, and the correlation is across samples,
not across time within one sample.

### 2.2 Three-dimensional information

True 3D atom coordinates are the weakest link. The state of the art is:

| Method | Lateral res. | Depth res. | Dose | Time cost |
|---|---|---|---|---|
| Multislice electron ptychography (single projection) | ~0.85 Å | ~6.6 nm | ~3.5×10³ e⁻/Å² | seconds–minutes |
| Focused-probe ptychography, beam-sensitive materials | ~2 Å | not depth-resolved | ~10² e⁻/Å² | seconds |
| Tilt-coupled multislice ptychography | sub-Å | **sub-nm** | higher | minutes (multi-tilt) |
| AET / tilt-series tomography | sub-Å | sub-Å | very high | minutes–hours |

Sources: [R8, R9, R10]. The essential point for this project is that **lateral
resolution and depth resolution are not on the same footing**. Single-projection
ptychography gives ~0.85 Å laterally but only ~6.6 nm in depth [R8]; sub-nm depth
resolution requires coupling several small-angle projections [R9], which costs
acquisition time and dose. Atomic-resolution tomography requires a tilt series, which
requires the structure to be *static* for the duration of the series — precisely the
assumption that a dynamic single-atom catalyst violates.

Consequence: **"continuous 3D atom tracking during reaction" does not currently exist.**
What exists is (a) continuous 2D projected tracking, and (b) intermittent
3D snapshots. Any honest platform must keep these two data types
typologically distinct.

### 2.3 Electronic state and adsorbates

Core-loss EELS (Pt M₄,₅; O K; Ce M₄,₅; Ti L₂,₃) and low-loss/vibrational EELS give
oxidation-state and bonding proxies, and monochromated vibrational EELS can in
principle report adsorbate signatures. In practice, at single-atom concentrations the
core-loss signal from one Pt atom is close to the noise floor at tolerable dose, so
per-atom charge state is currently a *statistical* rather than a per-atom observable.

---

## 3. Product-side analytics

### 3.1 Online mass spectrometry

Quadrupole MS is the default. Its practical envelope:

* sampling 1–10 Hz for multi-channel monitoring,
* detection limit ~ppm with a differentially pumped inlet,
* **response is convolved**, not instantaneous: the measured signal is the true
  outlet composition convolved with the impulse response of the reactor volume plus
  transfer line plus inlet, and offset by a dead time.

The convolution is the usual unstated assumption. In [R1] the correlation argument is
carried by the *periodicity* of both traces, which is robust to a constant delay. For
*event-based* correlation — "did this specific atomic rearrangement change the rate?" —
the delay and the dispersion must be measured, because the delay is typically of the
same order as, or larger than, the structural event rate.

### 3.2 Complementary ensemble spectroscopies

Operando XAS + DRIFTS + MS on Pt/Al₂O₃ during CO oxidation has resolved dispersion,
oxidation state and activity simultaneously, showing that poorly active atomically
dispersed Pt converts gradually and irreversibly into highly active ~1 nm clusters
over heating–cooling cycles [R11]. On CeO₂, the direction can reverse: clusters
disperse to single atoms under oxidising conditions, and for propane oxidation the
operando conversion of active Pt nanoclusters into single atoms is a *deactivation*
pathway [R12, R13]. Facet-dependent atomic redispersion of Pt on CeO₂ has been imaged
directly [R14], and fluxional behaviour of CeO₂-supported Pt has been correlated with
turnover frequency and with surface oxygen-vacancy concentration [R3].

These results matter because they establish the *scientific* content of the problem:
the sign of the structure–activity relationship for "single atom versus small cluster"
is system- and reaction-dependent and is not settled. That is what an event-resolved
method should be able to decide.

---

## 4. Electron-beam perturbation

The beam is a reagent. Beam-induced single-atom dynamics are well documented and have
even been developed into a deliberate tool: directed single-dopant manipulation in
silicon proceeds by an indirect-exchange mechanism, and the hopping is *stochastic* —
in one reported case a directed P hop occurred only after 68 ineffective 10-s spot
irradiations [R15, R16]. Beam-induced dynamics have been used to map potential-energy
landscapes [R17], and automated STEM has been used to probe beam-induced
transformations defect-by-defect [R18].

Quantitative context: single-atom targeting experiments operate at dose rates of order
10⁷ e⁻ Å⁻² s⁻¹ [R15], whereas low-dose ptychography of beam-sensitive materials
operates at ~10² e⁻ Å⁻² *total* [R10]. The two regimes are separated by many orders of
magnitude, which is the operating room a dose-aware protocol can exploit.

For a catalysis experiment the beam can additionally reduce the support, generate or
annihilate oxygen vacancies, crack adsorbates, and ionise the gas. Consequently, an
observed "atomic event" in an operando movie has two candidate causes — chemistry and
irradiation — and existing practice separates them only qualitatively ("we used low
dose").

---

## 5. Data-analysis and autonomous-microscopy practice

Deep-learning atom finding and segmentation of (S)TEM images is standard; real-time
multi-object tracking of in-situ TEM video has been demonstrated; synthetic-data
training pipelines with confidence scoring are established; and autonomous/agentic
microscope orchestration is an active and fast-moving area [R18, R19, R20, R21].
Reward-based ("rewards-driven") image analysis and dose-budgeted sequential
decision-making benchmarks for autonomous EM have appeared very recently [R22, R23].

What is mature: perception (find atoms), automation (drive the scope), and
retrospective analysis.

What is not established: using the perception output to *close a loop onto the
chemistry* — i.e. changing gas, potential or temperature in response to a detected
atomic event in order to steer a measured product rate.

---

## 6. Summary of the state of the art

1. Simultaneous atomic-resolution imaging and online product analysis at ~1 bar is
   demonstrated and reproducible [R1, R2].
2. Quantitative gas-phase EELS provides a near-sample cross-check of composition, and
   validated FEM models exist for the cell's transport fields [R4, R5].
3. Atom counting can be regressed against activity, but so far only statically and
   across samples [R7].
4. Continuous 3D atom tracking under reaction does not exist; 3D is intermittent and
   depth-limited [R8, R9, R10].
5. The sign and identity of the active Pt motif (single atom vs. sub-nm cluster) is
   system-dependent and contested [R3, R7, R11, R12, R13].
6. Beam-induced atom dynamics are real, stochastic, and only qualitatively controlled
   in operando practice [R15, R16, R17].
7. Correlation between structure and activity is, with very few exceptions, argued
   from *visual temporal coincidence* of two traces rather than from a
   transport-corrected, event-aligned statistical estimate.

Item 7, together with items 4 and 6, defines the gap analysed in `02_gap_analysis.md`.
