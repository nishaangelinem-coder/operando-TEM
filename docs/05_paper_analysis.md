# 05 — Critical Analysis of the Manuscript

A referee-style assessment of the paper this repository produces, written
against it rather than for it. Section numbers refer to `manuscript/manuscript.md`.

---

## 1. Summary of what is claimed

A methodology paper. It reframes operando electron microscopy of single-atom
catalysts as event-driven causal inference, specifies the platform and the
pre-registered measurement order that framing requires, builds a two-level
forward model of that platform, and measures the bias, power and failure modes
of eight estimators against known ground truth. It reports three design-changing
constraints (transport correction, the sensitivity/time-resolution trade-off,
the event-density limit on Pt loading) and two refutations of its own proposals
(dynamic descriptors do not out-explain a census; event-triggered control does
not beat a budget-matched schedule out of sample).

No experimental measurement is claimed, and this is stated in the title, the
abstract, a declaration before the abstract, every table and figure file, and
the repository README.

---

## 2. The strongest contributions

**2.1 The sensitivity/time-resolution trade-off (§4.8, Table S7b) is the most
valuable thing in the paper**, and it is arithmetic that anyone could have done
and nobody did. A quadrupole detects a mole fraction, so its flux limit falls
with flow; the transport dead time rises as the reciprocal of flow. The window
in which the chip's product is detectable is the window in which the smallest
resolvable lag is tens of seconds. This bounds what every published
"structure–activity correlation" from a nanoreactor experiment can mean, and it
is not a statement about any particular instrument's quality.

**2.2 The lead–lag artefact (Table S4b) is a clean, falsifiable, reproducible
demonstration.** Without the transport correction the peak lands at 3.0 s
against a dead time of 3.20 s, in 4 of 4 runs, and would be reported as
significant in 4 of 4. That is a specific, checkable prediction about what an
uncorrected published analysis would have found, and it is the kind of result
that changes practice because it is cheap to act on.

**2.3 The event-density criterion (Table S3a) converts a statistical
requirement into a synthesis specification.** λ(t_pre + t_post) ≲ 0.3 becomes a
maximum Pt surface density, roughly an order of magnitude below where
atomically dispersed Pt is usually prepared. That is actionable in a way that
"more statistics are needed" is not.

**2.4 The negative results are reported properly.** The in-sample versus
held-out closed-loop comparison (+6.9 % then −0.7 %, p = 0.90) is exactly the
comparison most papers omit, and reporting it is the paper's strongest signal of
good faith. The same is true of the hypothesis failing its own pre-specified
first test.

**2.5 The size-dependent detection bias (§4.1, §4.3) is a real and general
hazard.** A monomer at the detection threshold and a dimer above it means the
apparent dimer fraction is inflated — in exactly the direction that would support
a "dimers are the active site" conclusion. The paper's demand that detectors
report size-resolved recall is well founded and easy to comply with.

---

## 3. The central weakness

**Everything is validated on the author's own forward model, and the forward
model contains the answer.**

This is the objection a referee will lead with, and it is correct as far as it
goes. The simulator *imposes* that the dimer is the most productive motif; the
pipeline then recovers something whose sign and magnitude the simulator fixed.
There is no sense in which the paper discovers anything about platinum, and the
paper says so — but the danger is subtler than fabrication. A twin built by the
same person who built the estimators will tend to contain the failure modes
those estimators handle and to omit the ones they do not. Two specific
instances:

* The beam enters only as an additive term on the monomer hop rate. A real beam
  also reduces the support, creates and annihilates oxygen vacancies, cracks
  adsorbates and ionises the gas. The Beam Perturbation Index is *defined*
  against the one mechanism the twin implements, so its validation is close to
  circular. What survives the objection is narrower and should be stated as
  such: the dose-series intercept method is shown to be unbiased for an
  additive single-channel perturbation and to be misspecified as soon as the
  barrier distribution is non-degenerate — and the second half of that is the
  informative half, because it is a *failure* the author did not design in.
* Turnover is Langmuir–Hinshelwood with size-indexed parameters and no
  coverage-dependent restructuring. The structure→activity map is therefore
  imposed rather than emergent, which means the paper cannot speak to the case
  that most interests the field: a motif whose activity depends on the
  adsorbate configuration it is currently carrying.

**The honest reading** is that the paper establishes *necessary* conditions, not
sufficient ones. An estimator that fails on the twin will certainly fail on real
data. An estimator that passes has passed one concrete, fully specified test.
The paper's own §5.4 says this; a referee will want it said more prominently
than in a limitations subsection, and arguably in the abstract.

---

## 4. Specific methodological criticisms

**4.1 The counterfactual regime carries most of the estimator validation, and it
is physically impossible.** Studies S3, S4 and S7d all run at f_rep = 1 with a
detector of matching sensitivity — a regime the paper itself shows is excluded by
about ten orders of magnitude. The labelling is scrupulous, but a reader could
still come away thinking the event-aligned estimator has been validated for use.
It has been validated as an *estimator*; it has not been shown to have an
application. The paper should probably lead §4.4 with that tension rather than
resolve it two sections later.

**4.2 Power is thin.** Four to six simulator seeds throughout. The
position-heterogeneity ANOVA (p ≈ 2 × 10⁻⁴) and the null control result are
robust; the closed-loop comparison against open loop (p = 0.34) and the
observational/interventional contrast (2SLS s.d. 0.294 on a point estimate of
0.373) are not. Since seeds are free, there is no defence for this other than
compute budget, and the paper should either raise the seed count or drop the
underpowered comparisons rather than report them as suggestive.

**4.3 The weak instrument undermines the causal-inference claim.** First-stage
R² = 0.009 is a textbook weak instrument. The paper notes this and turns it into
a design requirement, which is the right move, but it means the paper's headline
framing — "the actuators are available as instruments for identification" — is
demonstrated only in principle. A referee is entitled to ask for a simulated
actuator strong enough to make the 2SLS estimate precise, which would test
whether the estimator works when its assumptions hold.

**4.4 One free parameter links the two levels of the twin, and the fit is only
7–9 %.** The chip mean field reproduces the kMC ensemble to 7.1 % on the monomer
trajectory and 8.7 % on the rate. Every chip-level number — the whole
representativeness argument — inherits that. The argument is robust because it
turns on ten orders of magnitude rather than ten per cent, but the paper should
say so explicitly instead of leaving the reader to notice.

**4.5 The catalyst-area multiplier (3.5 × 10⁴) is the single most load-bearing
number in the paper and it is an estimate.** It converts the projected area
ratio into f_rep and therefore sets the ten-orders-of-magnitude claim. It is
justified from a plausible film thickness and specific surface area, but it is
not measured, and a factor of 30 either way is defensible. The claim survives
easily — the shortfall is ~10¹¹, so even three orders of magnitude of error
changes nothing qualitatively — but the paper should present a sensitivity band
rather than a point value.

**4.6 The projection-correction failure is reported as an instrument
requirement, which is generous to the estimator.** An alternative reading is
that the estimator was poorly designed: a plane fit through five neighbours with
2.8 nm depth noise was never going to resolve a 0.7 nm corrugation, and this
could have been predicted analytically before it was run. The paper would be
stronger for saying that the analytic noise floor was computed and the estimator
run anyway to confirm it, rather than presenting the failure as a discovery.

**4.7 "Catalytic Event Microscopy" is a new name for an assembly of existing
methods.** Event-related averaging, circular-shift surrogates, Tikhonov
deconvolution with GCV, elastic nets, Gaussian HMMs and 2SLS are all standard.
The contribution is the assembly, the enforcement of two rules in code, and the
measurement of the operating envelope — which is a real contribution, but the
coinage oversells it. A referee will say so.

---

## 5. Internal consistency check

| Claim | Supported by | Verdict |
|---|---|---|
| Transfer function recoverable to ms | Table S1, 15 fits, coverage 13/15 and 15/15 | Supported |
| Uncorrected analysis invents a 3 s lead | Table S4b, 4/4 runs, dead time 3.20 s | Supported |
| Instrument costs only 15–25 % of the effect | Table S3a, Δ_dec/Δ_true = 0.85, 0.76 | Supported in the sparse regime only; stated as such |
| Event density limits Pt loading | Table S3a, monotonic over 60× in λ | Supported |
| Imaged field is ~10¹¹× below detection | Table S7a | Supported, conditional on the area multiplier (§4.5) |
| Time resolution is tens of seconds | Table S7b | Supported |
| Field diverges from chip under the beam | Table S7c, −17 % to −74 % | Supported |
| Dynamic descriptors add ~0 to instantaneous R² | Table S2b, ΔR² = 0.0006 | Supported |
| Dynamic descriptors add forecast skill | Table S2c, +0.027 at 120 s | Weakly supported; 5 seeds, small effect |
| Closed loop does not beat a matched schedule | Table S8b, p = 0.90 held out | Supported as a null; underpowered to exclude a small gain |
| BPI tracks truth under the linear model | Table S6c, within 15 % and 1 % | Supported inside the usable dose window |
| Projection bias not correctable at 6.6 nm | Table S5, η flat at 0.22 vs truth 0–0.41 | Supported |

No claim in the manuscript is unsupported by its own tables. Two are weaker than
the prose suggests (forecast skill, and the closed-loop null as an exclusion
rather than a failure to detect).

---

## 6. Research-integrity assessment

**Good.** Simulated provenance is declared in the title, in a pre-abstract
box, in §1.3, in §3.4, in §5.4, in every CSV header, every markdown table, every
figure caption sidecar and every JSON payload. The reference list distinguishes
what is cited for experimental precedent from what is cited for parameter
values, and states that no reference supports any simulated result. Code, seeds
and a one-command regeneration path are provided. Two of the author's own
hypotheses are reported as refuted. Threshold selection in the control study is
separated into selection and evaluation splits, with the in-sample optimism
shown alongside.

**Residual risks.** (i) Figures reproduce realistic micrographs; anyone
extracting `fig1_platform.png` from the repository without its caption sidecar
could mistake panel (a) for data. A visible watermark on rendered-micrograph
panels would close this. (ii) The parameter table in
`docs/04_experimental_setup.md` §5.4 mixes values anchored to publications with
plausible simulation choices in one table; the anchoring column distinguishes
them, but a stricter separation would be safer. (iii) The paper's framing invites
citation as though the protocol were validated in an experiment; the abstract's
first sentence should carry the word "simulated".

---

## 7. Journal fit

| Venue | Fit | Assessment |
|---|---|---|
| *Ultramicroscopy* / *Microscopy and Microanalysis* | **Best fit** | Methodology, instrument requirements, quantitative error budgets; a simulation-only methods paper is normal there |
| *ACS Catalysis* / *Journal of Catalysis* | Poor as submitted | Will be rejected for containing no catalysis measurement, however careful the framing |
| *npj Computational Materials* / *Digital Discovery* | Good | Digital twin plus open pipeline plus negative results is squarely in scope |
| *Nature Catalysis* / *Nature Materials* | No | Requires an experimental finding |
| *Review of Scientific Instruments* | Good for the platform half | The hardware specification and synchronisation acceptance test would stand alone |
| *Physical Review Materials* | Plausible | If reframed around the kinetic model and the transport/resolution bound |

**Recommendation:** split. The sensitivity/time-resolution bound, the lead–lag
artefact and the event-density criterion are a short, high-impact methods paper
for a microscopy or instrumentation venue. The full estimator suite, the digital
twin and the control study are a longer computational paper. Submitting the
combined manuscript to a catalysis journal would waste the strongest result.

---

## 8. Anticipated referee objections, and the defensible replies

**"This is a simulation dressed as a platform paper."**
Reply: the platform specification stands independently of the simulation — the
synchronisation acceptance test, the measurement order and the flow trade-off are
derivable from instrument parameters alone. The simulation is what allows the
estimators' bias and power to be quantified, which an experiment cannot do
because it has no ground truth. Concede that no estimator is thereby validated
for real data, and say so in the abstract.

**"The dimer result is circular."**
Reply: agreed, and it is not offered as a result. The validation target is the
estimator's recovery ratio against a known value, not the identity of the active
motif. The claims that are *not* circular are the ones about failure: the
event-density limit, the projection-correction failure, the dose-model
misspecification and the closed-loop null were not designed in.

**"Why should I believe the area multiplier?"**
Reply: present the shortfall as a function of it. The conclusion is unchanged
over three orders of magnitude, which is far wider than the plausible range.

**"Four seeds is not a study."**
Reply: valid. Raise the seed count for the underpowered comparisons, or remove
them.

**"The closed-loop null may just be low power."**
Reply: valid, and the paper should state the minimum detectable effect rather
than only the p-value. The point estimate is −0.7 % with a paired s.d. of about
115 on yields of ~860, so effects above roughly 10 % would have been detected;
the paper can exclude a large gain but not a small one.

**"Nothing here could not have been said in a comment on the existing
literature."**
Reply: the trade-off bound could have been. The quantified error budgets,
the operating envelope, the power curves and the two refutations could not.

---

## 9. What would most improve the paper

In priority order:

1. **Move the simulated-data declaration into the first sentence of the
   abstract**, not just above it.
2. **Raise the seed count** to at least 20 for S4c and S8, and report a minimum
   detectable effect for every null.
3. **Present f_rep as a sensitivity band** over the catalyst-area multiplier.
4. **Add a strong-instrument arm** to S4b so the 2SLS estimator is tested where
   its assumptions hold.
5. **Add a bursty-deactivation mode to the twin** and re-run S8; that is the
   condition under which event-triggered control should win, and testing it
   converts a null into a scoped claim.
6. **Watermark the rendered micrograph panels.**
7. **Split the manuscript** as in §7.
8. **Drop or downgrade "Catalytic Event Microscopy"** as a coinage; the assembly
   is the contribution and does not need a name.

---

## 10. Overall assessment

A careful, unusually self-critical methods paper whose best result — that the
product-linked time resolution of a nanoreactor operando experiment is tens of
seconds, and that a single field of view contributes a part in 10¹¹ of the
signal — is an arithmetic bound that constrains a whole literature. Its
weaknesses are real and mostly fixable: the validation is circular where it
concerns the imposed structure–activity map, several comparisons are
underpowered, and the framing occasionally implies more than simulation can
deliver. Its integrity practices are better than the norm, and its two
self-refutations are the reason to trust the rest.

**Verdict as submitted:** major revision for a microscopy or computational
venue; reject for a catalysis venue. **Verdict after the eight items in §9:**
a strong paper, with a genuinely load-bearing negative result.
