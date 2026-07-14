import os
import numpy as np
import pandas as pd

from prl_model import PLRModel
from stimulus import generate_protocol
from analysis import extract_tau_return, detect_hippus
from figures import plot_figure_1, plot_figure_2, plot_figure_4

def run_parameter_sweep(output_dir="output"):
    """
    Executes the full G sweep across all noise levels as defined in Section 3.
    """
    print("Starting parameter sweep...")
    G_values = np.round(np.arange(0.10, 2.30, 0.05), 2)  # Sweep up to 0.99*Gc = 2.295 per Section 3.2
    sigma_values = [0.0, 0.01, 0.05, 0.10, 0.20]
    
    results = []
    
    for sigma in sigma_values:
        print(f"  Simulating noise level sigma = {sigma}...")
        for G in G_values:
            # 1. Generate protocol
            t, stimulus, pulse_onsets = generate_protocol(G, num_pulses=10)
            
            # 2. Simulate model with strict OBSERVATIONAL noise (camera error)
            model = PLRModel(G, obs_noise_sigma=sigma)
            A = model.simulate(stimulus)
            
            # 3. Analyze trace (Explicitly pass pulse_duration to sync the dynamic window)
            tau_mean, tau_std, valid_fits = extract_tau_return(t, A, pulse_onsets, G, pulse_duration=0.200)
            is_hippus, band_power = detect_hippus(t, A)
            
            # 4. Store results
            results.append({
                "G": G,
                "sigma": sigma,
                "tau_return_mean": tau_mean,
                "tau_return_std": tau_std,
                "valid_fits": valid_fits,
                "is_hippus": is_hippus,
                "hippus_power": band_power
            })
            
    df_results = pd.DataFrame(results)
    output_path = os.path.join(output_dir, "sweep_results.csv")
    df_results.to_csv(output_path, index=False)
    print(f"Parameter sweep complete. Results saved to {output_path}")

def generate_trace_examples(output_dir="output"):
    """
    Generates specific, isolated traces for Figure 1.
    """
    print("Generating trace examples for Figure 1...")
    
    # --- Figure 1: Single pulse recoveries across G ---
    G_examples = [0.15,0.30,0.85,1.20,1.70,2.10, 2.35]
    trace_data = []

    dt_fig1   = 0.001
    t_pre     = 1.0    # 1 s baseline before pulse
    pulse_dur = 0.200
    t_post    = 12.0   # 12 s post-pulse — captures G=0.98 recovery
    stim_amp  = 0.3   # matches intensity in generate_protocol

    n_pre   = int(t_pre / dt_fig1)
    n_pulse = int(pulse_dur / dt_fig1)
    n_post  = int(t_post / dt_fig1)
    n_total = n_pre + n_pulse + n_post

    stimulus_fig1 = np.zeros(n_total)
    stimulus_fig1[n_pre : n_pre + n_pulse] = stim_amp

    t_raw     = np.arange(n_total) * dt_fig1
    t_shifted = t_raw - t_pre  # zero at stimulus onset

    for G in G_examples:
        model = PLRModel(G=G, obs_noise_sigma=0.0)
        A = model.simulate(stimulus_fig1)
        valid_idx = (t_shifted >= -0.5) & (t_shifted <= 12.0)
        for ti, ai in zip(t_shifted[valid_idx], A[valid_idx]):
            trace_data.append({"G": G, "time": round(ti, 4), "Area": ai})

    df_traces = pd.DataFrame(trace_data)
    trace_path = os.path.join(output_dir, "fig1_traces.csv")
    df_traces.to_csv(trace_path, index=False)

    print(f"Trace examples complete. Results saved to {trace_path}")

if __name__ == "__main__":
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    run_parameter_sweep(output_dir)
    generate_trace_examples(output_dir)

    print("Generating figures...")
    plot_figure_1(output_dir)
    plot_figure_2(output_dir)
    plot_figure_4(output_dir)
    print("Figure generation complete.")