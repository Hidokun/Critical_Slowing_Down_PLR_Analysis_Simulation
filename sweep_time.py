import numpy as np
import pandas as pd
from prl_model import PLRModel
from stimulus import generate_protocol
from analysis import extract_tau_return

def run_monte_carlo_sweep(G_values, noise_levels, num_trials=50, dt=0.001):
    results = []
    
    # 1. Optical Spatial Integration Factor
    # A standard clinical camera tracks the pupil boundary using ~400 perimeter pixels.
    # Therefore, 2D CCD static is attenuated when integrated into a 1D area measurement.
    N_pixels = 400.0
    integration_attenuation = 1.0 / np.sqrt(N_pixels)

    for sigma_obs in noise_levels:
        print(f"--- Running Monte Carlo for Observational Noise: {sigma_obs*100}% ---")
        
        # Apply the optical attenuation to the raw CCD noise
        effective_1d_noise = sigma_obs * integration_attenuation

        for G in G_values:
            t, stimulus, pulse_onsets = generate_protocol(
                G=G, num_pulses=1, dt=dt, pulse_duration=0.200, intensity=0.3
            )

            mask = t <= 15.0
            t = t[mask]
            stimulus = stimulus[mask]

            tau_trials = []

            for trial in range(num_trials):
                # Inject the physically correct effective 1D noise, NOT the raw 2D noise
                model = PLRModel(G=G, proc_noise_sigma=0.005, obs_noise_sigma=effective_1d_noise, dt=dt)
                A = model.simulate(stimulus)

                tau_mean, _, valid_fits = extract_tau_return(t, A, pulse_onsets, G, dt=dt)

                if valid_fits > 0 and np.isfinite(tau_mean):
                    tau_trials.append(tau_mean)

            if len(tau_trials) > 0:
                ensemble_mean = np.mean(tau_trials)
                ensemble_std = np.std(tau_trials)
                ci_lower = np.percentile(tau_trials, 2.5)
                ci_upper = np.percentile(tau_trials, 97.5)
            else:
                ensemble_mean = np.nan
                ensemble_std = np.nan
                ci_lower = np.nan
                ci_upper = np.nan

            results.append({
                'G': G,
                'sigma': sigma_obs,  # Record the nominal CCD noise for the legend
                'tau_return_mean': ensemble_mean,
                'tau_return_std': ensemble_std,
                'tau_return_ci_lower': ci_lower,
                'tau_return_ci_upper': ci_upper,
                'valid_trials': len(tau_trials),
                'total_trials': num_trials
            })

            print(f"G={G:.2f} | Valid: {len(tau_trials)}/{num_trials} | Mean Tau: {ensemble_mean:.3f}")

    return pd.DataFrame(results)

if __name__ == "__main__":
    import os
    np.random.seed(42)
    G_vals = np.linspace(0.15, 2.25, 40)
    noise_vals = [0.0, 0.01, 0.05, 0.10, 0.20]
    
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    
    print("Initiating Monte Carlo Noise Sweep...")
    df_sweep = run_monte_carlo_sweep(G_vals, noise_vals, num_trials=50)
    
    sweep_path = os.path.join(output_dir, 'sweep_results_monte_carlo.csv')
    df_sweep.to_csv(sweep_path, index=False)
    print(f"Monte Carlo sweep complete. Results saved to {sweep_path}")