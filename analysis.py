import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import welch, hilbert

def exponential_recovery(t, epsilon, tau_return, A_star=15.0):
    """
    Exponential decay function for fitting the re-dilation curve.
    A(t) = A* + epsilon * exp(-t / tau_return)
    """
    return A_star + epsilon * np.exp(-t / tau_return)

def extract_tau_return(t, A, pulse_onsets, G, pulse_duration=0.200, dt=0.001, A_star=15.0):
    """
    Extracts the recovery time constant (tau_return) from simulated PLR traces.
    
    Parameters:
        t (np.ndarray): Time vector.
        A (np.ndarray): Pupil area vector.
        pulse_onsets (np.ndarray): Array of stimulus onset times.
        G (float): Loop gain (used for predictions and boundary logic).
        pulse_duration (float): Duration of the stimulus pulse.
        dt (float): Integration step size.
        A_star (float): Resting pupil area.
        
    Returns:
        mean_tau (float): Averaged valid tau_return across pulses.
        std_tau (float): Standard deviation of valid tau_return across pulses.
        valid_fits (int): Number of pulses that produced a valid fit.
    """
    tau_iris = 0.311
    
    # Exclude near/post-bifurcation regime from exponential fitting
    if G >= 1.0:
        return np.nan, np.nan, 0
        
    tau_pred = tau_iris / (1.0 - G)
    T_IPI = max(2.0, 10.0 * tau_pred)
    
    taus = []
    residuals = []
    
    for onset in pulse_onsets:
        t_offset = onset + pulse_duration
        idx_offset = int(np.round(t_offset / dt))
        
        # Define maximum window end: T_IPI - 50 ms
        t_max_end = t_offset + T_IPI - 0.050
        idx_max_end = int(np.round(t_max_end / dt))
        
        # Ensure we don't index out of bounds
        idx_max_end = min(idx_max_end, len(t) - 1)
        
        # Use the full window, no premature truncation
        idx_end = idx_max_end
            
        # If the window is too short
        if (idx_end - idx_offset) * dt < 0.050:
            continue
            
        t_fit = t[idx_offset:idx_end] - t[idx_offset]  # Shift time to 0
        A_fit = A[idx_offset:idx_end]
        
        # Hilbert Envelope Extraction
        A_zero_mean = A_fit - A_star
        envelope = np.abs(hilbert(A_zero_mean))
        
        # Initial guesses
        epsilon_guess = envelope[0]
        epsilon_guess = max(epsilon_guess, 1e-5)
        tau_guess = min(max(tau_pred, 1e-4), 59.9)
        
        p0 = [epsilon_guess, tau_guess]
        
        # Bounds: epsilon strictly positive (amplitude), tau in (0, 60] seconds
        bounds = ([1e-5, 1e-5], [np.inf, 60.0])
        
        try:
            # Fit exponential envelope: eps * exp(-t / tau)
            fit_func = lambda t_val, eps, tau: eps * np.exp(-t_val / tau)
            popt, pcov = curve_fit(fit_func, t_fit, envelope, p0=p0, bounds=bounds)
            
            eps_fit, tau_fit = popt
            
            # Calculate residual (Mean Squared Error) against the envelope
            env_pred = fit_func(t_fit, eps_fit, tau_fit)
            mse = np.mean((envelope - env_pred)**2)
            
            # Exclusion threshold: since envelope is smooth, 0.5 is mathematically robust
            if mse < 0.5:  
                taus.append(tau_fit)
                residuals.append(mse)
                
        except RuntimeError:
            # Curve fit failed to converge
            continue
            
    if len(taus) == 0:
        return np.nan, np.nan, 0
        
    return np.mean(taus), np.std(taus), len(taus)


def detect_hippus(t, A, fs=1000.0, threshold_power=0.1):
    """
    Detects spontaneous pupillary oscillations (hippus) using power spectral density.
    Identifies if sustained spectral power in the 1-3 Hz band exceeds a threshold.
    
    Parameters:
        t (np.ndarray): Time vector.
        A (np.ndarray): Pupil area vector.
        fs (float): Sampling frequency in Hz (default 1000.0 for dt=1ms).
        threshold_power (float): Minimum power to classify as hippus state.
        
    Returns:
        is_hippus (bool): True if hippus detected.
        band_power (float): Integrated power in the 1-3 Hz band.
    """
    # Remove mean to avoid massive DC component
    A_zero_mean = A - np.mean(A)
    
    # Compute PSD using Welch's method
    # nperseg defines the frequency resolution. 4 seconds gives 0.25 Hz resolution.
    nperseg = int(4.0 * fs)
    if len(A_zero_mean) < nperseg:
        nperseg = len(A_zero_mean)
        
    freqs, psd = welch(A_zero_mean, fs, nperseg=nperseg)
    
    # Integrate power in the 1 to 3 Hz band
    band_mask = (freqs >= 1.0) & (freqs <= 3.0)
    band_power = np.trapezoid(psd[band_mask], freqs[band_mask])
    
    is_hippus = bool(band_power > threshold_power)
    
    return is_hippus, band_power


def compute_pipr(A_baseline_mean, A_post_mean):
    """
    Computes the Post-Illumination Pupil Response (PIPR) differential metric.
    R = ((A_PIPR - A0) / A0) * 100%
    
    Parameters:
        A_baseline_mean (float): Mean area over 2s pre-stimulus window (A0).
        A_post_mean (float): Mean area over 5-7s post-stimulus window (A_PIPR).
        
    Returns:
        R (float): Normalized percentage change.
    """
    return ((A_post_mean - A_baseline_mean) / A_baseline_mean) * 100.0