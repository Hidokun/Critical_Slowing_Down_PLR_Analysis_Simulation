import numpy as np

def generate_protocol(G, num_pulses=10, dt=0.001, pulse_duration=0.200, intensity=5.0):
    """
    Generates the stimulus protocol for a given loop gain G, using an adaptive
    inter-pulse interval (IPI) to guarantee full recovery to baseline.
    
    Parameters:
        G (float): Loop gain, used to predict recovery time and scale IPI.
        num_pulses (int): Number of consecutive light pulses (default: 10).
        dt (float): Integration step size in seconds (default: 0.001).
        pulse_duration (float): Duration of the square light pulse in seconds (default: 0.200).
        intensity (float): Multiplicative increase in effective retinal flux during the pulse.
        
    Returns:
        t (np.ndarray): Time vector for the simulation.
        stimulus (np.ndarray): Stimulus amplitude vector.
        pulse_onsets (np.ndarray): Array of times (in seconds) when pulses begin.
    """
    tau_iris = 0.311  # (s) from Longtin & Milton (1989)
    
    # Compute adaptive Inter-Pulse Interval (T_IPI)
    if G < 1.0:
        tau_pred = tau_iris / (1.0 - G)
        # Floor the IPI at 2.0 seconds for low G, otherwise 10x predicted tau
        T_IPI = max(2.0, 10.0 * tau_pred)
    else:
        # Near or past the bifurcation (oscillatory regime), fix IPI to a long window
        # since exponential recovery prediction breaks down.
        T_IPI = 10.0 
        
    t_baseline = 2.0  # 2 seconds of pure baseline before the first pulse
    
    # Calculate total duration
    t_total = t_baseline + num_pulses * (pulse_duration + T_IPI)
    
    # Generate time array
    num_steps = int(np.round(t_total / dt))
    t = np.arange(num_steps) * dt
    stimulus = np.zeros(num_steps)
    
    pulse_onsets = []
    
    # Populate the stimulus array with square pulses
    current_time = t_baseline
    for _ in range(num_pulses):
        pulse_onsets.append(current_time)
        
        # Convert time to array indices
        start_idx = int(np.round(current_time / dt))
        end_idx = int(np.round((current_time + pulse_duration) / dt))
        
        # Inject stimulus intensity
        stimulus[start_idx:end_idx] = intensity
        
        # Advance to next pulse onset
        current_time += (pulse_duration + T_IPI)
        
    return t, stimulus, np.array(pulse_onsets)