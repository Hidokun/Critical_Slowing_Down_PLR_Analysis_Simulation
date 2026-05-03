# Critical Slowing Down in the Pupillary Light Reflex
### A Computational Framework for Non-Invasive Locus Coeruleus Arousal Estimation

A high-fidelity simulation of the human Pupillary Light Reflex (PLR) as a **nonlinear delayed negative-feedback control system**, grounded in the Longtin-Milton (1989) Delay-Differential Equation. The project computationally validates a novel method for estimating **Locus Coeruleus (LC) tonic arousal state** from dynamical pupillometry — without any invasive procedure.

---

## The Problem

Standard clinical pupillometry treats the pupil as a static tape measure: *how wide is it?*

This completely discards the dynamical structure of the reflex. The **Locus Coeruleus (LC)** — the brain's primary noradrenergic nucleus — directly modulates the gain `G` of the pupil's reflexive control loop via its sympathetic projections to the iris dilator. This creates a direct mapping:

| Brain State | LC Tonic Output | Loop Gain `G` | Pupil Recovery |
|---|---|---|---|
| Hyperarousal / stress | High | Low | Fast, overdamped |
| Focused / alert | Moderate | Moderate | Normal |
| Fatigue / drowsiness | Low | High | Slow, ringing |
| Near bifurcation | Very low | `→ 2.33` | Spontaneous oscillations (Hippus) |

As `G` increases toward the critical threshold `G_c ≈ 2.33`, the system approaches a **Hopf bifurcation**, producing two measurable phenomena:

1. **Critical Slowing Down (CSD):** The exponential recovery time constant `τ_return` diverges — recovery becomes progressively slower
2. **Hippus:** Spontaneous ~1 Hz pupillary micro-oscillations emerge as a stable limit cycle

By measuring `τ_return` from a standardized 200 ms light pulse and inverting the exact DDE analytical relationship (via the Lambert W function), the underlying LC gain state `G` — and thus arousal level — can be recovered continuously and non-invasively.

---

## Key Results

### Figure 1 — Critical Slowing Down in the Time Domain
An identical 200 ms light pulse produces dramatically different recovery trajectories depending on `G`. At low gain the pupil snaps back immediately. As `G` approaches the Hopf boundary, recovery slows by orders of magnitude and the trace develops underdamped oscillatory ringing.

![Figure 1 — PLR Recovery Traces](figures/Figure_1_Traces.png)

---

### Figure 2 — Recovery Time Constant vs. Loop Gain
The exact analytical Lambert W theory (dashed line) precisely predicts the nonlinear DDE simulation results (orange dots). The **Valid Measurement Regime** spans from the Node-Focus Transition (`G ≈ 0.36`, where the system becomes underdamped) to the Hopf Bifurcation (`G ≈ 2.33`, where stability is lost). This defines the precise operating envelope of a non-invasive LC monitoring device.

![Figure 2 — τ_return vs G](figures/Figure_2_Tau_vs_G.png)

---

### Figure 3 — Spectral Emergence of Hippus
Welch power spectral density computed from simulated long-duration traces. Near the bifurcation (`G = 0.97`), a clear spectral peak emerges within the clinical hippus band (1–3 Hz), centred at `f_c ≈ 1.07 Hz` — matching the exact Hopf frequency predicted by the DDE characteristic equation. The Hill function nonlinearity bounds the instability into a finite-amplitude stable limit cycle, replicating true biological hippus.

![Figure 3 — Power Spectrum and Hippus](figures/Figure_3_Spectrum.png)

---

### Figure 4 — Noise Robustness
The Hilbert envelope extraction algorithm accurately tracks `τ_return` under realistic levels of retinal flux noise (σ up to 20%). Variance grows modestly with noise but the underlying trend is preserved, supporting feasibility of ambulatory wearable use.

![Figure 4 — Noise Robustness](figures/Figure_4_Noise.png)

---

## Mathematical Model

The governing equation is the **nonlinear Longtin-Milton DDE**:

```
τ_iris · dA/dt = −A(t) + A* + γ · [ f(Φ(t−δ)) − f(A*) ]
```

| Symbol | Meaning | Value |
|---|---|---|
| `A(t)` | Pupil area | mm² |
| `Φ(t) = A(t)·(1+stim)` | Retinal light flux | mm² |
| `f(Φ) = c·θⁿ/(θⁿ+Φⁿ)` | Hill negative feedback | — |
| `γ = G / \|f′(A*)\|` | Gain scaling factor | — |
| `δ` | Neural conduction delay | 0.300 s |
| `τ_iris` | Iris time constant | 0.311 s |
| `A*` | Resting pupil area | 15.0 mm² |

The recovery time constant is given by the exact **Lambert W inversion** of the DDE characteristic equation:

```
τ_return(G) = −1 / Re( W₀(−G·δ/τ · exp(δ/τ)) / δ  −  1/τ )
```

This replaces the classical simplified approximation `τ ≈ τ_iris/(1−G)`, which assumed zero delay and positive feedback, and incorrectly placed the bifurcation at `G = 1`.

---

## Project Structure

| File | Description |
|---|---|
| `prl_model.py` | Nonlinear DDE model with Hill function saturation |
| `stimulus.py` | Adaptive pulse protocol generator with IPI scaling |
| `analysis.py` | Hilbert envelope extraction, Lambert W τ-estimation, hippus detection |
| `main.py` | Full parameter sweep across G and noise levels |
| `figures.py` | Figure generation |

---

## Installation & Usage

```bash
pip install -r requirements.txt
```

```bash
# Step 1: Run the parameter sweep to generate simulation data
python main.py

# Step 2: Render figures from simulation data
python figures.py
```

All outputs (figures + CSVs) are saved to the `output/` directory.

---

## References

- Longtin, A. & Milton, J.G. (1989). *Modelling autonomous oscillations in the human pupil light reflex using nonlinear delay-differential equations.* Bulletin of Mathematical Biology, 51(5), 605–624.
- Longtin, A. & Milton, J.G. (1989). *Insight into the transfer function, gain, and oscillation onset for the pupil light reflex using nonlinear delay-differential equations.* Biological Cybernetics, 61(1), 51–58.
- Aston-Jones, G. & Cohen, J.D. (2005). *An integrative theory of locus coeruleus-norepinephrine function.* Annual Review of Neuroscience, 28, 403–450.
