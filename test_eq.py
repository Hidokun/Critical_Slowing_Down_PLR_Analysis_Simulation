from prl_model import PLRModel
from stimulus import generate_protocol
import numpy as np

# Equilibrium check
model = PLRModel(G=0.5, noise_sigma=0.0)
print(f"Nonlinear DDE Equilibrium: A_current = {model.A_current:.4f}, A* = {model.A_star:.4f}")
print()

# First-passage diagnostic
for G in [0.30, 0.60, 0.85, 0.95, 0.98]:
    t, stimulus, pulse_onsets = generate_protocol(G, num_pulses=1)
    model = PLRModel(G, noise_sigma=0.0)
    A = model.simulate(stimulus)
    
    idx_start = int((pulse_onsets[0] + 0.200) / 0.001)
    t_seg = t[idx_start:] - t[idx_start]
    A_seg = A[idx_start:]
    
    within_2pct = np.where(np.abs(A_seg - 15.0) < 0.30)[0]  # 2% of A*=15
    if len(within_2pct) > 0:
        print(f"G={G:.2f}  T_return(2%) = {t_seg[within_2pct[0]]:.3f}s")
    else:
        print(f"G={G:.2f}  never returned within 2%")