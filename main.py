import os
import numpy as np
import pandas as pd

from prl_model import PLRModel
from stimulus import generate_protocol
from analysis import extract_tau_return, detect_hippus

def run_parameter_sweep(output_dir="output"):
    """
    Executes the full G sweep across all noise levels as defined in Section 3.
    """
    print("Starting parameter sweep...")
    G_values = np.round(np.arange(0.10, 1.00, 0.01), 2)
    sigma_values = [0.0, 0.01, 0.05, 0.10, 0.20]
    
    results = []
    
    for sigma in sigma_values:
        print(f"  Simulating noise level sigma = {sigma}...")
        for G in G_values:
            # 1. Generate protocol
            t, stimulus, pulse_onsets = generate_protocol(G, num_pulses=10)
            
            # 2. Simulate model
            model = PLRModel(G, noise_sigma=sigma)
            A = model.simulate(stimulus)
            
            # 3. Analyze trace
            tau_mean, tau_std, valid_fits = extract_tau_return(t, A, pulse_onsets, G)
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
    Generates specific, isolated traces for Figure 1 (time domain transitions)
    and Figure 3 (spectral/hippus emergence).
    """
    print("Generating trace examples for Figures 1 and 3...")
    
    # --- Figure 1: Single pulse recoveries across G ---
    G_examples = [0.30, 0.60, 0.85, 0.95, 0.98]
    trace_data = []
    
    for G in G_examples:
        # Generate exactly 1 pulse with a long enough baseline to settle
        t, stimulus, pulse_onsets = generate_protocol(G, num_pulses=1)
        model = PLRModel(G, noise_sigma=0.0) # Zero noise for clean qualitative overlay
        A = model.simulate(stimulus)
        
        # Shift time so stimulus onset is at t=0
        t_shifted = t - pulse_onsets[0]
        
        # Only keep data from -1.0s to 10.0s (or max available)
        valid_idx = (t_shifted >= -1.0) & (t_shifted <= 10.0)
        
        for ti, ai in zip(t_shifted[valid_idx], A[valid_idx]):
            trace_data.append({"G": G, "time": ti, "Area": ai})
            
    df_traces = pd.DataFrame(trace_data)
    trace_path = os.path.join(output_dir, "fig1_traces.csv")
    df_traces.to_csv(trace_path, index=False)
    
    # --- Figure 3: Spontaneous Hippus (Oscillatory regime) ---
    # Simulate a 60-second window with ambient noise (sigma=0.05), NO stimulus
    hippus_data = []
    dt = 0.001
    t_hippus = np.arange(0, 60.0, dt)
    zero_stimulus = np.zeros_like(t_hippus)
    
    for G in [0.50, 0.97]:
        model = PLRModel(G, noise_sigma=0.05)
        A = model.simulate(zero_stimulus)
        
        for ti, ai in zip(t_hippus, A):
            hippus_data.append({"G": G, "time": ti, "Area": ai})
            
    df_hippus = pd.DataFrame(hippus_data)
    hippus_path = os.path.join(output_dir, "fig3_hippus_traces.csv")
    df_hippus.to_csv(hippus_path, index=False)
    
    print(f"Trace examples complete. Results saved to {trace_path} and {hippus_path}")

if __name__ == "__main__":
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    run_parameter_sweep(output_dir)
    generate_trace_examples(output_dir)