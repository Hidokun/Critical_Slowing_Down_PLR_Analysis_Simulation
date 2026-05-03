import numpy as np
from scipy.special import lambertw


def generate_protocol(G, num_pulses=10, dt=0.001,
                      pulse_duration=0.200, intensity=0.3):
    """
    Generates the stimulus protocol for a given loop gain G, using an adaptive
    inter-pulse interval (IPI) set to 5 * tau_pred to guarantee recovery to
    within 1% of A* before the next pulse arrives.

    tau_pred is computed from the exact Lambert W eigenvalue formula (Eq. 11),
    valid across the full operative range G in [0.15, 2.1].

    Parameters:
        G (float): Loop gain.
        num_pulses (int): Number of consecutive light pulses (default: 10).
        dt (float): Integration step size in seconds (default: 0.001).
        pulse_duration (float): Duration of the square light pulse in seconds
                                (default: 0.200).
        intensity (float): Normalised stimulus intensity (default: 0.3).

    Returns:
        t (np.ndarray): Time vector for the simulation.
        stimulus (np.ndarray): Stimulus amplitude vector.
        pulse_onsets (np.ndarray): Array of pulse onset times in seconds.
    """
    tau_iris = 0.311
    delta    = 0.300

    z       = G * delta / tau_iris * np.exp(delta / tau_iris)
    W       = lambertw(-z, k=0)
    lam     = np.real(W) / delta - 1.0 / tau_iris
    tau_pred = -1.0 / lam
    T_IPI   = max(5.0 * tau_pred, 2.0)

    t_baseline = 2.0
    t_total    = t_baseline + num_pulses * (pulse_duration + T_IPI)

    num_steps = int(np.round(t_total / dt))
    t         = np.arange(num_steps) * dt
    stimulus  = np.zeros(num_steps)

    pulse_onsets  = []
    current_time  = t_baseline

    for _ in range(num_pulses):
        pulse_onsets.append(current_time)
        start_idx = int(np.round(current_time / dt))
        end_idx   = int(np.round((current_time + pulse_duration) / dt))
        stimulus[start_idx:end_idx] = intensity
        current_time += (pulse_duration + T_IPI)

    return t, stimulus, np.array(pulse_onsets)