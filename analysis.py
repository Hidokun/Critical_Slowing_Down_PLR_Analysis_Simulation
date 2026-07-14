import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import welch, butter, filtfilt

def constrained_damped_oscillator(t, y0, e, tau, omega):
    """
    4-parameter morphological fit. 
    Phase is algebraically locked (maximum displacement forced at t=0).
    """
    return y0 + e * np.exp(-t / tau) * np.cos(omega * t)

def extract_tau_return(t, A, pulse_onsets, G, pulse_duration=0.200, dt=0.001, A_star=12.0):
    if G >= 2.10:
        return np.nan, np.nan, 0

    # 4.0 Hz Butterworth preserves the 1.07 Hz biological envelope
    fs = 1.0 / dt
    b, a = butter(2, 4.0, btype='low', fs=fs)
    padlen = min(3 * max(len(a), len(b)), len(A) - 1)
    A_smooth = filtfilt(b, a, A, padlen=padlen)

    taus = []

    for i, onset in enumerate(pulse_onsets):
        if i + 1 < len(pulse_onsets):
            end_t = pulse_onsets[i+1] - 0.05
        else:
            end_t = t[-1]
            
        pre_stim_mask = (t >= onset - 1.0) & (t < onset)
        if np.any(pre_stim_mask):
            empirical_baseline = np.mean(A[pre_stim_mask])
        else:
            empirical_baseline = A_star

        window_mask = (t >= onset) & (t < end_t)
        t_window = t[window_mask]
        A_raw_window = A[window_mask]
        A_smooth_window = A_smooth[window_mask]

        if len(t_window) < 100:
            continue

        search_mask = (t_window - onset) <= 1.5
        if not np.any(search_mask):
            continue
            
        peak_idx = np.argmin(A_smooth_window[search_mask])
        t_peak = t_window[peak_idx]

        redil_mask = (t_window >= t_peak)
        t_fit = t_window[redil_mask] - t_peak  
        A_raw_fit = A_raw_window[redil_mask]
        A_smooth_fit = A_smooth_window[redil_mask] 

        if len(t_fit) < 50:
            continue

        try:
            displacement = np.min(A_raw_fit) - empirical_baseline
            
            p0_guess = [empirical_baseline, displacement, 1.0, 6.7]
            
            # Restore the broad, honest physical bounds
            bounds_strict = (
                [empirical_baseline - 0.2, -np.inf, 0.05, 3.0], 
                [empirical_baseline + 0.2, 0.0, 15.0, 10.0]
            )
            
            # STEP 1: Guided Initializer with tuned armor
            popt_guess, _ = curve_fit(
                constrained_damped_oscillator, 
                t_fit, 
                A_smooth_fit, 
                p0=p0_guess, 
                bounds=bounds_strict, 
                method='trf',
                loss='soft_l1',
                f_scale=0.1,  # Forces robust rejection of noise-induced wobbles
                maxfev=5000
            )
            
            tau_guide = popt_guess[2]
            omega_guide = popt_guess[3]

            # STEP 2: The Exact Extraction with tuned armor
            p0_exact = [empirical_baseline, displacement, tau_guide, omega_guide]
            
            popt_exact, _ = curve_fit(
                constrained_damped_oscillator, 
                t_fit, 
                A_raw_fit, 
                p0=p0_exact, 
                bounds=bounds_strict, 
                method='trf',
                loss='soft_l1',
                f_scale=0.1,  # The master key. Rejects CCD static unconditionally.
                maxfev=5000
            )
            
            taus.append(popt_exact[2])
            
        except RuntimeError:
            pass

    if len(taus) > 0:
        return np.mean(taus), np.std(taus), len(taus)
    else:
        return np.nan, np.nan, 0

def detect_hippus(t, A, fs=1000.0, nperseg=4000):
    A_zero_mean = A - np.mean(A)
    freqs, psd = welch(A_zero_mean, fs, nperseg=nperseg)
    band_mask = (freqs >= 0.1) & (freqs <= 5.0)
    valid_psd = psd[band_mask]
    if len(valid_psd) == 0:
        return False, 0.0
    peak_power = np.max(valid_psd)
    return peak_power > 0.5, peak_power