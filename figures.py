import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import welch


plt.style.use('default')
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'lines.linewidth': 2,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

def plot_figure_1(output_dir="output"):
    """
    Figure 1: PLR traces at varying G overlaid.
    Shows the qualitative transition from fast damped return to critical slowing down.
    """
    trace_path = os.path.join(output_dir, "fig1_traces.csv")
    if not os.path.exists(trace_path):
        print(f"Skipping Fig 1: {trace_path} not found.")
        return

    df = pd.read_csv(trace_path)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = sns.color_palette("rocket", n_colors=len(df['G'].unique()))
    
    # Plot each G trace
    for idx, G in enumerate(sorted(df['G'].unique())):
        data = df[df['G'] == G]
        ax.plot(data['time'], data['Area'], label=f'$G = {G:.2f}$', color=colors[idx])
        
    # Formatting and annotations
    ax.axvspan(0, 0.200, color='gray', alpha=0.2, label='Stimulus Pulse (200 ms)')
    ax.axhline(15.0, color='k', linestyle='--', linewidth=1, alpha=0.5, label='Resting Area $A^*$')
    
    # Text annotation for the degenerate trace
    ax.text(0.3, 11.5, "Degenerate return\n($\tau_{return} < 200$ ms)", 
            color=colors[0], fontsize=10, va='center')
    
    ax.set_xlim(-0.5, 6.0)
    ax.set_ylim(10.0, 15.5)
    ax.set_xlabel("Time from stimulus onset (s)")
    ax.set_ylabel(r"Pupil Area ($mm^2$)")
    ax.set_title("PLR Recovery Dynamics Approaching Bifurcation")
    ax.legend(loc='lower right', framealpha=0.9)
    
    fig.savefig(os.path.join(output_dir, "Figure_1_Traces.png"))
    plt.close(fig)

def plot_figure_2(output_dir="output"):
    """
    Figure 2: tau_return vs G.
    Empirical curve from simulation overlaid with theoretical tau/(1-G).
    """
    sweep_path = os.path.join(output_dir, "sweep_results.csv")
    if not os.path.exists(sweep_path):
        print(f"Skipping Fig 2: {sweep_path} not found.")
        return

    df = pd.read_csv(sweep_path)
    # Use zero-noise data for the baseline validation
    df_clean = df[df['sigma'] == 0.0].dropna(subset=['tau_return_mean'])
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Theoretical curve (linearised approximation)
    G_theory = np.linspace(0.1, 0.99, 100)
    tau_iris = 0.311
    tau_theory = tau_iris / (1 - G_theory)
    
    ax.plot(G_theory, tau_theory, 'k--', linewidth=1.5, 
            label=r'Theoretical: $\tau_{return} \approx \tau_{iris}/(1-G)$')
    
    # Simulated points (full nonlinear DDE)
    ax.scatter(df_clean['G'], df_clean['tau_return_mean'], 
               color='#d95f02', s=40, zorder=5, label='Simulated DDE (Zero Noise)')
    
    # Validity shading
    degen_bound = 0.36
    ax.axvspan(degen_bound, 1.0, color='#4daf4a', alpha=0.08, label='Valid Measurement Regime')
    ax.axvline(degen_bound, color='#e41a1c', linestyle=':', alpha=0.8, 
               label=r'Degeneration Bound ($\tau_{return} = 200$ ms)')
    
    ax.set_yscale('log')
    ax.set_xlim(0.1, 1.0)
    ax.set_ylim(0.1, 50.0)
    ax.set_xlabel("Loop Gain ($G$)")
    ax.set_ylabel(r"Recovery Time Constant $\tau_{return}$ (s)")
    ax.set_title("Critical Slowing Down: Simulation vs Theory")
    ax.legend(loc='upper left')
    ax.grid(True, which="major", ls="-", alpha=0.2)
    ax.grid(True, which="minor", ls=":", alpha=0.1)
    
    fig.savefig(os.path.join(output_dir, "Figure_2_Tau_vs_G.png"))
    plt.close(fig)

def plot_figure_3(output_dir="output"):
    """
    Figure 3: Power spectrum at G=0.50 vs G=0.97.
    Shows the emergence of the 1/(2*delta) spectral peak near bifurcation.
    """
    hippus_path = os.path.join(output_dir, "fig3_hippus_traces.csv")
    if not os.path.exists(hippus_path):
        print(f"Skipping Fig 3: {hippus_path} not found.")
        return

    df = pd.read_csv(hippus_path)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = {'0.5': '#377eb8', '0.97': '#e41a1c'}
    labels = {'0.5': 'G = 0.50 (Stable state)', '0.97': 'G = 0.97 (Near bifurcation)'}
    
    for G_val in [0.50, 0.97]:
        # Exact floating point matching can be tricky, cast to string key
        data = df[np.isclose(df['G'], G_val)]['Area'].values
        if len(data) == 0:
            continue
            
        data_zero_mean = data - np.mean(data)
        
        # Welch's method PSD (dt=0.001 -> fs=1000 Hz)
        freqs, psd = welch(data_zero_mean, fs=1000.0, nperseg=4000)
        
        mask = (freqs >= 0.1) & (freqs <= 5.0)
        ax.plot(freqs[mask], psd[mask], label=labels[str(G_val)], 
                color=colors[str(G_val)], linewidth=2, alpha=0.8)
        
    # Annotate clinical and theoretical bounds
    ax.axvspan(1.0, 3.0, color='gray', alpha=0.15, label='Clinical Hippus Band (1-3 Hz)')
    f_c = 1.0 / (2.0 * 0.300) # 1 / (2*delta)
    ax.axvline(f_c, color='k', linestyle='--', linewidth=1.5, 
               label=r'Theoretical Hopf Freq ($f_c = 1.67$ Hz)')
    
    ax.set_xlim(0.1, 5.0)
    ax.set_yscale('log')
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel(r"Power Spectral Density ($mm^4$/Hz)")
    ax.set_title("Emergence of Spontaneous Pupillary Oscillations (Hippus)")
    ax.legend(loc='upper right')
    ax.grid(True, which="both", ls=":", alpha=0.2)
    
    fig.savefig(os.path.join(output_dir, "Figure_3_Spectrum.png"))
    plt.close(fig)

def plot_figure_4(output_dir="output"):
    """
    Figure 4: Noise robustness.
    Fitted tau_return across four sigma levels with error shading.
    """
    sweep_path = os.path.join(output_dir, "sweep_results.csv")
    if not os.path.exists(sweep_path):
        print(f"Skipping Fig 4: {sweep_path} not found.")
        return

    df = pd.read_csv(sweep_path)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    sigmas = [0.01, 0.05, 0.10, 0.20]
    colors = sns.color_palette("viridis", n_colors=len(sigmas))
    
    # Plot true zero-noise curve as baseline
    df_clean = df[df['sigma'] == 0.0].dropna(subset=['tau_return_mean'])
    ax.plot(df_clean['G'], df_clean['tau_return_mean'], 'k-', 
            linewidth=2, alpha=0.7, label='Zero Noise Baseline')
    
    for idx, sigma in enumerate(sigmas):
        df_sig = df[(df['sigma'] == sigma) & (df['valid_fits'] > 0)].dropna(subset=['tau_return_mean'])
        if len(df_sig) == 0:
            continue
            
        label_str = rf'$\sigma = {int(sigma*100)}\%$'
        
        # Plot mean
        ax.plot(df_sig['G'], df_sig['tau_return_mean'], 
                label=label_str, color=colors[idx], marker='o', markersize=4, linestyle='none')
        
        # Add 1-standard-deviation shading
        ax.fill_between(df_sig['G'], 
                        df_sig['tau_return_mean'] - df_sig['tau_return_std'],
                        df_sig['tau_return_mean'] + df_sig['tau_return_std'],
                        color=colors[idx], alpha=0.15)
                        
    ax.set_yscale('log')
    # Focus the x-axis on the valid measurement regime
    ax.set_xlim(0.36, 0.95)
    ax.set_ylim(0.2, 30.0)
    ax.set_xlabel("Loop Gain ($G$)")
    ax.set_ylabel(r"Estimated $\tau_{return}$ (s) $\pm 1$ SD")
    ax.set_title("Measurement Robustness Under Retinal Flux Noise")
    ax.legend(loc='upper left', ncol=2)
    ax.grid(True, which="major", ls="-", alpha=0.2)
    ax.grid(True, which="minor", ls=":", alpha=0.1)
    
    fig.savefig(os.path.join(output_dir, "Figure_4_Noise.png"))
    plt.close(fig)

if __name__ == "__main__":
    output_dir = "output"
    if not os.path.exists(output_dir):
        print(f"Error: {output_dir} directory not found. Run main.py first.")
    else:
        plot_figure_1(output_dir)
        plot_figure_2(output_dir)
        plot_figure_3(output_dir)
        plot_figure_4(output_dir)