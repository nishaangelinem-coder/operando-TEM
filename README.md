# Catalytic Event Microscopy (`opcem`)

**Transport-corrected, beam-aware causal inference from operando electron
microscopy of single-atom platinum catalysts — a methodology and in-silico
validation.**

---

## ⚠️ Every number in this repository is simulated

No physical electron microscope, MEMS nanoreactor, mass spectrometer or catalyst
was operated in this work. All results are outputs of the forward simulator in
`src/opcem/truth.py` and the virtual instruments in `src/opcem/instruments.py`.
Every table, figure and raw output file carries a machine-readable banner saying
so, together with the git commit and the random seed that produced it.

The contribution is a **measurement-and-analysis protocol with characterised
operating limits**, validated against known ground truth — not an experimental
finding about platinum.

---

## What problem this addresses

Operando electron microscopy can watch individual Pt atoms move while a catalyst
works, and can record the products at the same time. It cannot yet establish
*which* atom or configuration made the product. The obstacle is not spatial
resolution; it is that the two data streams are compared by visual temporal
coincidence while three effects sit between them:

1. the **gas-transport delay** between the imaged region and the mass
   spectrometer (seconds — the same order as the atomic events);
2. the **electron beam's own contribution** to the atomic dynamics;
3. the **tiny, non-representative fraction** of the working catalyst that one
   field of view contains.

`opcem` formalises the experiment as event-driven causal inference, specifies the
hardware and software it requires, and measures what each estimator in it can
and cannot do.

## Headline results from the in-silico validation

| Result | Where |
|---|---|
| Tracer pulses recover the transport transfer function to < 5 ms in dead time and < 10 ms in dispersion over a 4× flow range | `docs/results.md`, Table S1 |
| Skipping the transport correction drives the recovered event effect to ~0 and puts the apparent lead–lag peak at the dead time | Table S4a/S4b |
| Event-aligned analysis needs event density λ(t_pre+t_post) ≲ 0.3, i.e. a **maximum Pt loading** | Table S3a |
| The imaged field's product flux is ~7×10¹¹ times below a quadrupole's detection limit — spontaneous single-atom event analysis is *excluded*, not merely hard | Table S7a |
| The flow that makes the chip's product detectable makes the smallest resolvable lag ~17–67 s | Table S7b |
| The illuminated field's dispersed fraction differs from the un-illuminated chip's by −17 % at zero dose and −74 % at 10⁴ e⁻ Å⁻² s⁻¹ | Table S7c |
| Dynamic descriptors add ~0 to explaining the *instantaneous* rate (census already gives R² ≈ 0.99) but do add forecast skill at 60–120 s | Tables S2b/S2c |
| Event-triggered control shows +6.9 % in sample and **−0.7 %, p = 0.90, out of sample** against a budget-matched fixed schedule | Table S8b |

Several of these are negative results. They are the point.

## Layout

```
docs/       01 methods review · 02 gap analysis · 03 proposal ·
            04 experimental setup (hardware + digital twin) · 05 paper analysis
src/opcem/  config · truth (kMC + chip mean field) · instruments ·
            transport · vision · descriptors · correlate · beam ·
            mobility · control · report
studies/    s0…s8 + run_all.py + figures.py
tests/      22 property tests
results/    tables/ (csv + md) · figures/ (png + pdf + captions) · raw/ (json)
manuscript/ manuscript.md
```

## Reproduce

```bash
pip install -r requirements.txt
python tests/test_opcem.py            # 22 property tests
python studies/run_all.py             # every table and figure, ~1 h
python studies/run_all.py s3 s6       # selected studies
python studies/run_all.py figures     # figures only, from cached json
```

Master seed `20260922`. Each output records the git commit and its seed.

## The two-level digital twin

* **`PatchSimulator`** — discrete-time kinetic Monte Carlo of Pt entities in one
  40 × 40 nm illuminated field: site-dependent hopping with an explicit
  beam term `β·φ`, diffusion-limited aggregation, oxygen-accelerated
  dissociation, a corrugated support so that projection genuinely loses
  out-of-plane motion, and a Langmuir–Hinshelwood turnover in which the **dimer
  is the most productive motif per Pt atom**. The pipeline is never told this.
* **`chip_mean_field`** — mass-conserving size-resolved population balance for
  the whole, **un-illuminated** reactive zone, position-binned along the flow
  axis. The two levels share one fitted coefficient (the encounter rate).

## Author

Dr. M. Nisha Angeline — Professor and Head, Department of Electronics and
Communication Engineering, Velalar College of Engineering and Technology,
Thindal, Erode.

## Licence

MIT (see `LICENSE`).
