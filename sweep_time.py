import numpy as np
from scipy.special import lambertw

tau_iris = 0.311
delta = 0.300
G_values = np.round(np.arange(0.50, 2.10, 0.05), 2)
pulse_duration = 0.200
num_pulses = 10

total_time = 0
breakdown = []

for G in G_values:
    arg = -(G * delta / tau_iris) * np.exp(delta / tau_iris)
    W = lambertw(arg, k=0)
    lam = (W / delta) - (1.0 / tau_iris)
    alpha = np.real(lam)
    tau_pred = -1.0 / alpha
    
    T_IPI = max(5.0 * tau_pred, 0.5)
    block_time = num_pulses * (pulse_duration + T_IPI)
    total_time += block_time
    breakdown.append((G, round(tau_pred, 3), round(T_IPI, 3), round(block_time, 1)))

print(f"{'G':>6}  {'tau_pred (s)':>12}  {'T_IPI (s)':>10}  {'Block (s)':>10}")
print("-" * 46)
for row in breakdown:
    print(f"{row[0]:>6.2f}  {row[1]:>12.3f}  {row[2]:>10.3f}  {row[3]:>10.1f}")

print(f"\nTotal sweep time: {total_time:.1f} s  ({total_time/60:.1f} min)")