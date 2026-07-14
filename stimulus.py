import numpy as np
from scipy.special import lambertw

def generate_protocol(G, num_pulses=10, dt=0.001, pulse_duration=0.200, intensity=0.3):
    """
    Generates the stimulus protocol for a given loop gain G, using an adaptive
    inter-pulse interval (IPI) set to 5*tau_pred to guarantee recovery to within
    1% of A* before the next pulse arrives.
    """
    tau_iris = 0.311
    delta = 0.300
    
    # Calculate theoretical tau via exact DDE Lambert W solution
    z = (G * delta / tau_iris) * np.exp(delta / tau_iris)
    W = lambertw(-z, k=0)
    lam = np.real(W) / delta - 1.0 / tau_iris
    
    if lam >= 0:
        # Post-bifurcation (G >= Gc): System is in a limit cycle, no exponential decay
        T_IPI = 10.0  # Fixed observation interval
    else:
        tau_pred = -1.0 / lam
        # Cap max IPI at 60 seconds to prevent RAM overflow near the Hopf boundary
        T_IPI = min(max(5.0 * tau_pred, 2.0), 60.0)

    t_baseline = 2.0
    t_total = t_baseline + (num_pulses * (pulse_duration + T_IPI))
    num_steps = int(np.round(t_total / dt))
    
    t = np.arange(num_steps) * dt
    stimulus = np.zeros(num_steps)
    pulse_onsets = []
    
    current_time = t_baseline
    for _ in range(num_pulses):
        pulse_onsets.append(current_time)
        
        start_idx = int(np.round(current_time / dt))
        end_idx = int(np.round((current_time + pulse_duration) / dt))
        
        if start_idx < num_steps and end_idx <= num_steps:
            stimulus[start_idx:end_idx] = intensity
            
        current_time += (pulse_duration + T_IPI)
        
    return t, stimulus, np.array(pulse_onsets)