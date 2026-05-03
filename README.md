# NEXUS — Nonlinear EXact Unified Simulation of the Pupillary Light Reflex

A high-fidelity computational neuroscience framework for simulating the human **Pupillary Light Reflex (PLR)** as a nonlinear delayed negative-feedback control system. The simulation is grounded in the Longtin-Milton (1989) Delay-Differential Equation and is designed to validate a novel method for non-invasively estimating **Locus Coeruleus (LC) tonic arousal state** from dynamical pupillometry.

---

## Scientific Background

The standard clinical view of pupillometry treats the pupil as a static tape measure: *how wide is it?* This project argues that the **dynamics** of the pupil's recovery after a light pulse , how fast or slow it returns to baseline,  contain far richer neurological information than the static diameter alone.

The Locus Coeruleus (LC), the brain's primary noradrenergic nucleus, directly modulates the gain `G` of the pupil's reflexive control loop via its sympathetic projections to the iris dilator muscle. This means:

- **High LC tonic output** (hyperarousal) → low gain `G` → fast, overdamped recovery
- **Low LC tonic output** (fatigue, drowsiness) → high gain `G` → slow, critically damped recovery

As gain `G` increases toward a critical threshold, the system approaches a **Hopf bifurcation**, which physically manifests as two observables:

1. **Critical Slowing Down (CSD):** The exponential recovery time constant `τ_return` diverges as `G → G_c ≈ 2.33`
2. **Hippus:** Spontaneous ~1 Hz pupillary oscillations that emerge as a stable limit cycle beyond the bifurcation

By measuring `τ_return` from a brief, standardized light pulse, this framework proposes that the underlying LC gain state can be recovered by numerically inverting the exact analytical relationship (via the Lambert W function) — providing a **continuous, non-invasive biomarker** of brainstem arousal.

---

## Key Results

### Figure 1 — Critical Slowing Down in the Time Domain
The pupil's response to an identical 200 ms light pulse changes dramatically as `G` increases. At low gain the pupil snaps back immediately. As `G` approaches the Hopf boundary, recovery slows by orders of magnitude and the trace develops underdamped ringing before settling.

![Figure 1](output/Figure_1_Traces.png)

---

### Figure 2 — Recovery Time Constant vs. Loop Gain
The exact analytical theory (Lambert W) perfectly predicts the simulated nonlinear DDE recovery times. The Valid Measurement Regime spans the Node-Focus Transition (`G ≈ 0.36`, where the system becomes underdamped) to the Hopf Bifurcation (`G ≈ 2.33`, where the system becomes oscillatory). This defines the precise operating envelope of a wearable LC pupillometer.

![Figure 2](output/Figure_2_Tau_vs_G.png)

---

### Figure 3 — Spectral Emergence of Hippus
Power spectral density computed via Welch's method shows that at `G = 0.97`, a clear spectral peak emerges within the clinical hippus band (1–3 Hz), centred at `f_c ≈ 1.07 Hz` — matching the exact Hopf frequency predicted by the DDE characteristic equation for the Longtin-Milton parameters. The Hill function nonlinearity traps the instability into a finite-amplitude stable limit cycle, replicating true biological hippus.

![Figure 3](output/Figure_3_Spectrum.png)

---

### Figure 4 — Noise Robustness
The exponential envelope extraction algorithm (using the Hilbert transform) remains accurate under realistic levels of retinal flux noise (σ up to 20%). Measurement variance grows modestly with noise but tracking of the underlying `τ_return` trend is preserved, supporting the feasibility of ambulatory use.

![Figure 4](output/Figure_4_Noise.png)

---

## Mathematical Foundations

The governing equation is the **nonlinear Longtin-Milton DDE**:

```
τ_iris · dA/dt = -A(t) + A* + γ · [ f(Φ(t−δ)) − f(A*) ]
```

Where:
- `A(t)` — pupil area (mm²)
- `Φ(t) = A(t) · (1 + stimulus)` — retinal light flux
- `f(Φ) = c·θⁿ / (θⁿ + Φⁿ)` — Hill feedback function
- `γ = G / |f′(A*)|` — gain scaling ensuring linearized loop gain is exactly `G`
- `δ = 0.300 s` — neural conduction delay
- `τ_iris = 0.311 s` — iris time constant
- `A* = 15.0 mm²` — resting pupil area

The recovery time constant is given by the exact **Lambert W inversion** of the DDE characteristic equation:

```
τ_return(G) = −1 / Re( W₀(−G·δ/τ · exp(δ/τ)) / δ − 1/τ )
```

This replaces earlier simplified approximations (e.g., `τ ≈ τ_iris/(1−G)`) which assumed zero delay and positive feedback, and incorrectly placed the bifurcation at `G = 1`.

---

## Project Structure

| File | Description |
|---|---|
| `prl_model.py` | Core nonlinear DDE model with Hill function saturation |
| `stimulus.py` | Adaptive stimulus protocol generator with IPI scaling |
| `analysis.py` | Hilbert envelope extraction, Lambert W τ-estimation, hippus detection |
| `main.py` | Full parameter sweep across G and noise levels |
| `figures.py` | figure generation |

---

## Installation & Usage

```bash
pip install -r requirements.txt
```

```bash
# Step 1: Run parameter sweep and generate simulation data
python main.py

# Step 2: Render figures from the simulation data
python figures.py
```

All outputs (figures + CSVs) are saved to the `output/` directory.

---

## Dependencies

```
numpy
scipy
pandas
matplotlib
seaborn
```

---

## References

- Longtin, A. & Milton, J.G. (1989). *Modelling autonomous oscillations in the human pupil light reflex using nonlinear delay-differential equations.* Bulletin of Mathematical Biology, 51(5), 605–624.
- Longtin, A. & Milton, J.G. (1989). *Insight into the transfer function, gain, and oscillation onset for the pupil light reflex using nonlinear delay-differential equations.* Biological Cybernetics, 61(1), 51–58.
- Aston-Jones, G. & Cohen, J.D. (2005). *An integrative theory of locus coeruleus-norepinephrine function.* Annual Review of Neuroscience, 28, 403–450.
