import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import welch, hilbert
from scipy.special import lambertw

def exponential_recovery(t, epsilon, tau_return, A_star=15.0):
    """
    Exponential decay function for fitting the re-dilation curve.
    A(t) = A* + epsilon * exp(-t / tau_return)
    """
    return A_star + epsilon * np.exp(-t / tau_return)
def extract_tau_return(t, A, pulse_onsets, G, pulse_duration=0.200, dt=0.001, A_star=15.0):
    """
    Extracts the recovery time constant (tau_return) from simulated PLR traces.
    Uses Hilbert envelope for oscillatory recoveries (G > 0.5) and direct
    absolute deviation for monotonic recoveries (G <= 0.5).
    """
    tau_iris = 0.311

    if G >= 2.318:  # The true Hopf bifurcation is around G=2.318 for tau=0.311
        return np.nan, np.nan, 0

    # Exact theoretical envelope decay using Lambert W function
    delta = 0.300
    arg = - (G * delta / tau_iris) * np.exp(delta / tau_iris)
    W = lambertw(arg, k=0)
    lam = (W / delta) - (1.0 / tau_iris)
    alpha = np.real(lam)
    
    tau_pred = -1.0 / alpha
    T_IPI = max(5.0 * tau_pred, 0.5)

    taus = []
    residuals = []

    for onset in pulse_onsets:
        t_offset = onset + pulse_duration
        idx_offset = int(np.round(t_offset / dt))

        t_max_end = t_offset + T_IPI - 0.050
        idx_max_end = int(np.round(t_max_end / dt))
        idx_max_end = min(idx_max_end, len(t) - 1)
        idx_end = idx_max_end

        if (idx_end - idx_offset) * dt < 0.050:
            continue

        t_fit = t[idx_offset:idx_end] - t[idx_offset]
        A_fit = A[idx_offset:idx_end]
        A_zero_mean = A_fit - A_star

        zero_crossings = np.where(np.diff(np.sign(A_zero_mean)))[0]

        if len(zero_crossings) >= 2:
            envelope = np.abs(hilbert(A_zero_mean))
            # Hilbert transform has edge artifacts, so clip the last 10%
            clip_idx = int(len(envelope) * 0.9)
            t_fit = t_fit[:clip_idx]
            A_fit = A_fit[:clip_idx]
            envelope = envelope[:clip_idx]
            mse_threshold = 0.1
        else:
            envelope = np.abs(A_zero_mean)
            mse_threshold = 0.1

        # Find the minimum of A_fit (constriction nadir) and start fitting from there
        nadir_idx = np.argmin(A_fit)
        if nadir_idx > 5:
            t_fit = t_fit[nadir_idx:] - t_fit[nadir_idx]
            envelope = envelope[nadir_idx:]

        if envelope[0] < 0.05:
            continue

        epsilon_guess = float(np.clip(envelope[0], 1e-5, A_star))
        tau_guess = float(np.clip(tau_pred, 1e-4, 59.9))
        p0 = [epsilon_guess, tau_guess]
        bounds = ([1e-5, 1e-5], [A_star * 2.0, 60.0])

        try:
            fit_func = lambda t_val, eps, tau: eps * np.exp(-t_val / tau)
            popt, _ = curve_fit(fit_func, t_fit, envelope,
                                p0=p0, bounds=bounds, maxfev=8000)
            eps_fit, tau_fit = popt

            env_pred = fit_func(t_fit, eps_fit, tau_fit)
            mse = np.mean((envelope - env_pred) ** 2)

            if mse < mse_threshold:
                taus.append(tau_fit)
                residuals.append(mse)

        except RuntimeError:
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