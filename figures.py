import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.special import lambertw

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

def _find_existing_file(output_dir, candidates):
    for name in candidates:
        path = os.path.join(output_dir, name)
        if os.path.exists(path):
            return path
    return None

def _load_sweep_dataframe(output_dir):
    sweep_path = _find_existing_file(output_dir, ['sweep_results_monte_carlo.csv', 'sweep_results.csv'])
    if sweep_path is None:
        return None, None
    return pd.read_csv(sweep_path), sweep_path

def plot_figure_1(output_dir='output'):
    """ Figure 1: PLR traces at varying G overlaid. """
    trace_path = _find_existing_file(output_dir, ['fig1_traces.csv', 'fig1_traces_new.csv'])
    if trace_path is None:
        print(f"Skipping Fig 1: no trace file found in {output_dir}.")
        return
        
    df = pd.read_csv(trace_path)
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = sns.color_palette("rocket", n_colors=len(df['G'].unique()))
    
    for idx, G in enumerate(sorted(df['G'].unique())):
        data = df[df['G'] == G]
        if G >= 2.318:
            lbl = f"G = {G:.2f} (Post-bifurcation limit cycle)"
        else:
            lbl = f"G = {G:.2f}"
        ax.plot(data['time'], data['Area'], label=lbl, color=colors[idx])
        
    ax.axvspan(0, 0.200, color='gray', alpha=0.2, label='Dynamical Block Stimulus (200 ms)')
    ax.axhline(12.0, color='k', linestyle='--', linewidth=1, alpha=0.5, label=r'Mesopic Baseline $A^* = 12.0$ mm$^2$')
    
    ax.set_xlim(-0.5, 12.0)
    ax.set_ylim(8.0, 13.5)
    ax.set_xlabel('Time from stimulus onset (s)')
    ax.set_ylabel(r'Pupil Area (mm$^2$)')
    ax.set_title('PLR Recovery Dynamics Approaching Bifurcation')
    
    # De-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='lower right', framealpha=0.9, fontsize=9)
    
    fig.savefig(os.path.join(output_dir, 'Figure_1_Traces.png'))
    plt.close(fig)

def plot_figure_2(output_dir='output'):
    """ Figure 2: tau_return vs G. Empirical curve vs theoretical Lambert W. """
    df, sweep_path = _load_sweep_dataframe(output_dir)
    if df is None:
        print(f"Skipping Fig 2: no sweep file found in {output_dir}.")
        return
        
    df_clean = df[df['sigma'] == 0.0].dropna(subset=['tau_return_mean'])
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Theoretical curve: Exact DDE via Lambert W
    G_theory = np.linspace(0.1, 2.38, 200)
    tau_iris = 0.311
    delta = 0.300
    tau_theory = np.zeros_like(G_theory)
    
    for i, g_val in enumerate(G_theory):
        arg = - (g_val * delta / tau_iris) * np.exp(delta / tau_iris)
        W = lambertw(arg, k=0)
        lam = (W / delta) - (1.0 / tau_iris)
        tau_theory[i] = -1.0 / np.real(lam)
        
    ax.plot(G_theory, tau_theory, 'k--', linewidth=1.5, label=r'Theoretical Lambert $W_0$ Exact DDE')
    
    # Simulated points
    ax.scatter(df_clean['G'], df_clean['tau_return_mean'], color='#d95f02', s=40, zorder=5, label='Simulated DDE (Zero Observational Noise)')
    
    degen_bound = 0.15
    kinematic_bound = 2.10
    
    ax.axvspan(degen_bound, kinematic_bound, color='#4daf4a', alpha=0.08, label=r'Operative Window $G \in [0.15, 2.1]$')
    ax.axvline(degen_bound, color='#e41a1c', linestyle=':', alpha=0.8, label=r'Injectivity Limit $G^* = 0.145$')
    ax.axvline(kinematic_bound, color='#ff7f00', linestyle=':', alpha=0.8, label=r'Kinematic Ceiling $G = 2.1$')
    ax.axvline(2.318, color='#984ea3', linestyle='-.', alpha=0.8, label=r'Hopf Bifurcation $G_c = 2.318$')
    
    ax.set_yscale('log')
    ax.set_xlim(0.1, 2.5)
    ax.set_ylim(0.1, 100.0)
    ax.set_xlabel('Loop Gain $G$')
    ax.set_ylabel(r'Envelope Recovery Time Constant $\tau_{return}$ (s)')
    ax.set_title('Critical Slowing Down: Exact DDE Theory vs Simulation')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, which='major', ls='-', alpha=0.2)
    ax.grid(True, which='minor', ls=':', alpha=0.1)
    
    fig.savefig(os.path.join(output_dir, 'Figure_2_Tau_vs_G.png'))
    plt.close(fig)

def plot_figure_4(output_dir='output'):
    """ Figure 4: Noise robustness. Uses either Monte Carlo or standard sweep outputs. """
    df, _ = _load_sweep_dataframe(output_dir)
    if df is None:
        print(f"Skipping Fig 4: no sweep file found in {output_dir}.")
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    sigmas = [0.01, 0.05, 0.10, 0.20]
    
    # User-defined absolute contrast sequence
    colors = ['red', 'blue', 'black', 'yellow']

    df_clean = df[df['sigma'] == 0.0].dropna(subset=['tau_return_mean'])
    if len(df_clean) > 0:
        # Baseline is Black, slightly thicker
        ax.plot(df_clean['G'], df_clean['tau_return_mean'], 'k-', linewidth=2.5, alpha=0.9, label='Zero Noise Baseline')

    for idx, sigma in enumerate(sigmas):
        df_sig = df[df['sigma'] == sigma].copy()

        if {'tau_return_ci_lower', 'tau_return_ci_upper'}.issubset(df_sig.columns):
            df_valid = df_sig.dropna(subset=['tau_return_mean', 'tau_return_ci_lower', 'tau_return_ci_upper'])
            if len(df_valid) == 0:
                continue
            label_str = rf'Obs. Noise $\sigma_{{obs}}$ = {int(sigma*100)}%'
            
            # Distinct colors for each scatter plot and variance band
            ax.plot(df_valid['G'], df_valid['tau_return_mean'], label=label_str, color=colors[idx], marker='o', markersize=4, linestyle='none')
            ax.fill_between(df_valid['G'], df_valid['tau_return_ci_lower'], df_valid['tau_return_ci_upper'], color=colors[idx], alpha=0.15)
            
        elif 'valid_fits' in df_sig.columns:
            df_valid = df_sig[(df_sig['valid_fits'] > 0)].dropna(subset=['tau_return_mean'])
            if len(df_valid) == 0:
                continue
            label_str = rf'Obs. Noise $\sigma_{{obs}}$ = {int(sigma*100)}%'
            
            ax.plot(df_valid['G'], df_valid['tau_return_mean'], label=label_str, color=colors[idx], marker='o', markersize=4, linestyle='none')
            ax.fill_between(df_valid['G'], df_valid['tau_return_mean'] - df_valid['tau_return_std'], df_valid['tau_return_mean'] + df_valid['tau_return_std'], color=colors[idx], alpha=0.15)

    ax.set_yscale('log')
    ax.set_xlim(0.15, 2.35)
    ax.set_ylim(0.1, 100.0)
    ax.set_xlabel('Loop Gain $G$')
    ax.set_ylabel(r'Estimated $\tau_{return}$ (s)')
    ax.set_title('Measurement Robustness Under CCD Observational Noise')

    ax.legend(loc='upper left', ncol=2)
    ax.grid(True, which='major', ls='-', alpha=0.2)
    ax.grid(True, which='minor', ls=':', alpha=0.1)

    fig.savefig(os.path.join(output_dir, 'Figure_4_Noise.png'))
    plt.close(fig)

if __name__ == '__main__':
    output_dir = 'output'
    if not os.path.exists(output_dir):
        print(f"Error: {output_dir} directory not found. Run main.py first.")
    else:
        plot_figure_1(output_dir)
        plot_figure_2(output_dir)
        plot_figure_4(output_dir)