import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import zeta
from scipy.stats import poisson, geom
import os
from scipy.stats import poisson



def plot_degree_evolution_separate(storage, stats, arrival_times, n0, m0, tmax,
                                   model_name='BA_pref', 
                                   plot_theory=True,
                                   plot_rescaled=True,
                                   save_dir='./figures'):
    """
    Create TWO separate figures for ANY network model with theoretical predictions:
    1. Raw degree plots (grid 2x2, one per tau)
    2. Combined rescaled/data collapse plot
    
    Theoretical predictions:
    - BA_pref: k(t) ≈ m₀√(t/tau)
    - BA_random: k(t) ≈ m₀[log(n₀+t-1) - log(n₀+tau-1) + 1]
    - fixed_pref: k(t) ≈ (2m₀/n₀)·t  (for large n₀, t ≥ n₀)
    
    Parameters
    ----------
    storage : dict
        Raw trajectories
    stats : dict
        Statistics per arrival time
    arrival_times : list
        Arrival times to plot
    n0 : int
        Initial nodes (network size for fixed_pref)
    m0 : int
        Connections per step
    tmax : int
        Maximum time
    model_name : str
        'BA_pref', 'BA_random', or 'fixed_pref'
    plot_theory : bool
        Whether to plot theoretical prediction
    plot_rescaled : bool
        Whether to create rescaled/collapse plot
    save_dir : str
        Output directory
        
    Returns
    -------
    fig_raw : matplotlib.figure.Figure
    fig_rescaled : matplotlib.figure.Figure or None
    """
    os.makedirs(save_dir, exist_ok=True)
    
    n_times = len(arrival_times)
    colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray']
    
    model_titles = {
        'BA_pref': 'BA with Preferential Attachment',
        'BA_random': 'BA with Random Attachment',
        'fixed_pref': 'Fixed-Size with Preferential Attachment'
    }
    title_main = model_titles.get(model_name, model_name)
    
    fig_raw, axes_raw = plt.subplots(2, 2, figsize=(14, 10))
    fig_raw.suptitle(f'{title_main} — Raw Degree Evolution k(t)', 
                     fontsize=16, y=0.995)
    axes_raw = axes_raw.flatten()
    
    for i, tau in enumerate(arrival_times):
        mean = stats[tau]["mean"]
        std = stats[tau]["std"]
        t_vals = np.arange(tau, tmax + 1)
        
        ax = axes_raw[i]
        
        ax.plot(t_vals, mean, label="Empirical mean k(t)", 
               color=colors[i % len(colors)], linewidth=2.5)
        ax.fill_between(t_vals, mean - std, mean + std,
                       alpha=0.3, label="±1 std", 
                       color=colors[i % len(colors)])
        
        if plot_theory:
            if model_name == 'BA_pref':
                theory_raw = m0 * np.sqrt(t_vals / tau)
                label_theory = r"Theory: $m_0\sqrt{t/\tau}$"
                
            elif model_name == 'BA_random':
                theory_raw = m0 * (np.log(n0 + t_vals - 1) - np.log(n0 + tau - 1) + 1)
                label_theory = r"Theory: $m_0[\log(n_0+t-1) - \log(n_0+\tau-1) + 1]$"
                
            else:
                theory_raw = None
                label_theory = None
            
            if theory_raw is not None:
                ax.plot(t_vals, theory_raw, '--', 
                       label=label_theory, 
                       color="black", linewidth=2)
        
        ax.set_xlabel("t", fontsize=11)
        ax.set_ylabel("k(t)", fontsize=11)
        ax.set_title(f"tau = {tau}", fontsize=12, fontweight='bold')
        ax.legend(loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for i in range(n_times, len(axes_raw)):
        axes_raw[i].set_visible(False)
    
    plt.tight_layout()
    
    # Save figure 1
    path_raw = os.path.join(save_dir, f'{model_name}_raw_degree_evolution.png')
    fig_raw.savefig(path_raw, dpi=300, bbox_inches='tight')
    print(f"✓ Raw degree plots saved: {path_raw}")
    
    fig_rescaled = None
    
    if plot_rescaled:
        fig_rescaled, ax_rescaled = plt.subplots(1, 1, figsize=(12, 8))
        
        if model_name == 'BA_pref':
            fig_rescaled.suptitle(f'{title_main} — Data Collapse', 
                                 fontsize=16, y=0.995)
            
            # Universal curve
            t_theory = np.arange(1, tmax + 1)
            theory_universal = m0 * np.sqrt(t_theory)
            ax_rescaled.plot(t_theory, theory_universal, '--', 
                            label=r"Theory: $m_0\sqrt{t}$ (universal)", 
                            color="black", linewidth=3.5, zorder=10)
            
            ylabel = r"$k'(t) = \sqrt{\tau} \cdot k(t)$"
            subtitle = "All arrival times collapse onto the universal curve"
            
            # Rescaling function
            rescale_func = lambda mean, tau: np.sqrt(tau) * mean
            
        elif model_name == 'BA_random':
            fig_rescaled.suptitle(f'{title_main} — Data Collapse', 
                                 fontsize=16, y=0.995)
            
            # Universal curve
            t_theory = np.arange(1, tmax + 1)
            theory_universal = m0 * np.log(n0 + t_theory - 1)
            ax_rescaled.plot(t_theory, theory_universal, '--', 
                            label=r"Theory: $m_0\log(n_0+t-1)$ (universal)", 
                            color="black", linewidth=3.5, zorder=10)
            
            ylabel = r"$k''(t) = k(t) + m_0\log(n_0+\tau-1) - m_0$"
            subtitle = "All arrival times collapse onto the universal curve"
            
            # Rescaling function
            rescale_func = lambda mean, tau: mean + m0 * np.log(n0 + tau - 1) - m0
            
            
        else:
            # Fallback
            fig_rescaled.suptitle(f'{title_main} — Rescaled Degree', 
                                 fontsize=16, y=0.995)
            ylabel = r"$k'(t)$"
            subtitle = "Rescaled degree evolution"
            rescale_func = lambda mean, tau: mean
        
        # Plot curves for all tau
        for i, tau in enumerate(arrival_times):
            mean = stats[tau]["mean"]
            t_vals = np.arange(tau, tmax + 1)
            
            # Apply rescaling
            k_rescaled = rescale_func(mean, tau)
            
            ax_rescaled.plot(t_vals, k_rescaled, 
                            label=f"tau = {tau}", 
                            color=colors[i % len(colors)], 
                            linewidth=2.5, 
                            alpha=0.85)
        
        ax_rescaled.set_xlabel("t", fontsize=13)
        ax_rescaled.set_ylabel(ylabel, fontsize=13)
        ax_rescaled.set_title(subtitle, fontsize=13, fontweight='bold', pad=20)
        ax_rescaled.legend(loc='best', fontsize=11, framealpha=0.95, ncol=2)
        ax_rescaled.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure 2
        path_rescaled = os.path.join(save_dir, f'{model_name}_rescaled_data_collapse.png')
        fig_rescaled.savefig(path_rescaled, dpi=300, bbox_inches='tight')
        print(f"✓ Rescaled plot saved: {path_rescaled}")
        
        plt.show()
    else:
        plt.show()
    
    return fig_raw, fig_rescaled


def plot_best_models_grid(stats, arrival_times, tmax, results_df, fit_params,
                          model_name='model', metric_name='RMSE',
                          save_path=None, figsize=(18, 12)):
    """
    Plot best fitted models in 2x2 grid.
    """
    n_times = min(len(arrival_times), 4)
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle(f'Best Model Fits: {model_name}', fontsize=18, fontweight='bold', y=0.995)
    axes = axes.flatten()
    
    for idx, tau in enumerate(arrival_times[:4]):
        ax = axes[idx]
        
        mean = stats[tau]["mean"]
        std = stats[tau]["std"]
        t_vals = np.arange(tau, tmax + 1)
        
        row = results_df[results_df['tau'] == tau]
        if len(row) == 0:
            ax.text(0.5, 0.5, f'No data for tau={tau}',
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        best_model = row['best_model'].values[0]
        aic = row['AIC'].values[0]
        metric_value = row[metric_name].values[0]
        
        best_fit = fit_params[tau][best_model]
        predictions = best_fit['predictions']
        
        # Plot with error bars (subsample)
        step = max(1, len(t_vals) // 100)
        ax.errorbar(t_vals[::step], mean[::step], yerr=std[::step],
                   fmt='o', label='Data', color='steelblue',
                   alpha=0.6, markersize=4, capsize=3, elinewidth=1.5, zorder=3)
        
        ax.plot(t_vals, predictions, '-', label=f'{best_model}',
               color='darkred', linewidth=3, zorder=5)
        
        ax.set_xlabel('Time t', fontsize=12, fontweight='bold')
        ax.set_ylabel('Degree k(t)', fontsize=12, fontweight='bold')
        ax.set_title(f'tau = {tau}', fontsize=14, fontweight='bold')
        
        textstr = f'{best_model}\nAIC = {aic:.1f}\n{metric_name} = {metric_value:.4f}'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))
        
        ax.legend(loc='lower right', fontsize=10, framealpha=0.95)
        ax.grid(True, alpha=0.35, linestyle='--')
    
    for idx in range(n_times, 4):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()
    
    return fig, axes


def plot_no_growth(times, k_all, n0, m0, n_display=10, 
                           t_min=None, figsize=(18, 8), save_fig=True):
    """
    Create two plots:
    1. Individual degree trajectories (sample of nodes)
    2. Mean degree with ±1 std band
    
    Parameters:
    -----------
    times : array
        Time points
    k_all : array
        Degrees of all vertices (shape: n0 x T+1)
    n0 : int
        Total number of vertices
    m0 : int
        Number of edges added per time step
    n_display : int
        Number of individual trajectories to show
    t_min : int, optional
        Minimum time for x-axis (default: 0)
    figsize : tuple
        Figure size
    save_fig : bool
        If True, save figure to ./figures directory
    """
    if t_min is None:
        t_min = 0
    
    # theoretical curve
    k_theory = 2 * m0 * times / n0
    
    # compute mean and std across all nodes
    k_mean = k_all.mean(axis=0)
    k_std = k_all.std(axis=0)
    
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    ax1 = axes[0]
    
    # select n_display random nodes to display
    rng = np.random.default_rng(42)
    display_indices = rng.choice(n0, size=n_display, replace=False)
    
    for idx in display_indices:
        ax1.plot(times, k_all[idx], lw=1, alpha=0.5)
    
    # theoretical line
    ax1.plot(times, k_theory, "k--", lw=2.5, label=r"Theory: $2m_0 t / n_0$")
    
    ax1.set_xlim(t_min, times[-1])
    ax1.set_xlabel("time $t$", fontsize=12)
    ax1.set_ylabel(r"degree $k_i(t)$", fontsize=12)
    ax1.set_title(f"Individual trajectories ({n_display} nodes shown)", fontsize=13)
    ax1.legend(loc="best", framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[1]
    
    ax2.fill_between(times, k_mean - k_std, k_mean + k_std,
                     color='salmon', alpha=0.4, label=r"$\pm 1$ std")
    
    # plot mean
    ax2.plot(times, k_mean, color='red', lw=2, label=r"Empirical mean $\bar{k}(t)$")
    
    # theoretical line
    ax2.plot(times, k_theory, "k--", lw=2.5, label=r"Theory: $2m_0 t / n_0$")
    
    ax2.set_xlim(t_min, times[-1])
    ax2.set_xlabel("time $t$", fontsize=12)
    ax2.set_ylabel(r"degree $k(t)$", fontsize=12)
    ax2.set_title(f"Mean degree across all $n_0={n0}$ nodes", fontsize=13)
    ax2.legend(loc="best", framealpha=0.9)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    if save_fig:
        os.makedirs("./figures", exist_ok=True)
        
        T = times[-1]
        filename = f"./figures/no_growth_n0_{n0}_m0_{m0}_T_{T}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {filename}")
    
    plt.show()


def plot_best_model_no_growth(times, k_all, results_dict, all_fit_params, 
                              show_std=True, save_fig=True, figsize=(10, 6)):
    """
    Plot empirical data, best fitted model, and theoretical curve for no-growth model.
    
    Parameters
    ----------
    times : array
        Time points
    k_all : array
        Degrees of all vertices (shape: n0 x T+1)
    results_dict : dict
        Results from model_selection_no_growth
    all_fit_params : dict
        All fitted models from model_selection_no_growth
    show_std : bool
        If True, show ±1 std band around empirical mean
    save_fig : bool
        If True, save figure to ./figures directory
    figsize : tuple
        Figure size
    """
    
    n0 = results_dict['n0']
    m0 = results_dict['m0']
    t_min = results_dict['t_min']
    best_model = results_dict['best_model']
    best_info = results_dict['best_info']
    
    if best_info is None:
        print("Error: No valid fit found!")
        return
    
    k_mean = k_all.mean(axis=0)
    k_std = k_all.std(axis=0)
    
    # Theoretical curve (full range)
    k_theory = 2 * m0 * times / n0
    
    # Best model predictions
    t_fit = best_info['t_vals']
    y_fit = best_info['predictions']
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if show_std:
        ax.fill_between(times, k_mean - k_std, k_mean + k_std,
                        color='lightblue', alpha=0.3, label=r'$\pm 1$ std')
    
    ax.plot(times, k_mean, 'o', color='steelblue', markersize=2, 
            alpha=0.6, label='Empirical mean $\\bar{k}(t)$')
    
    ax.plot(t_fit, y_fit, '-', color='red', lw=2.5, 
            label=f'Best fit: {best_model}')
    
    ax.plot(times, k_theory, 'k--', lw=2, 
            label=r'Theory: $k(t) = 2m_0 t / n_0$')
    
    # Vertical line at t_min
    ax.axvline(t_min, color='gray', linestyle=':', lw=1.5, alpha=0.7,
               label=f'Fit range: $t \\geq {t_min}$')
    
    ax.set_xlabel('time $t$', fontsize=12)
    ax.set_ylabel(r'degree $k(t)$', fontsize=12)
    
    param_str = ", ".join([f"{name}={val:.4f}" 
                          for name, val in zip(best_info['param_names'], 
                                              best_info['params'])])
    title = f"No-growth PA: $n_0={n0}$, $m_0={m0}$ | Best: {best_model}\n"
    title += f"AIC={best_info['AIC']:.1f}, RMSE={best_info['RMSE']:.4f} | {param_str}"
    ax.set_title(title, fontsize=11)
    
    ax.legend(loc='best', framealpha=0.9, fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_fig:
        os.makedirs("./figures", exist_ok=True)
        T = times[-1]
        filename = f"./figures/no_growth_fit_n0_{n0}_m0_{m0}_T_{T}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {filename}")
    
    plt.show()

def plot_degree_distributions(
    all_deg_pref, 
    all_deg_rand, 
    degrees_no_growth,
    n0_ba=None, 
    m0_ba=None, 
    tmax_ba=None,
    n0_ng=None,
    m0_ng=None,
    tmax_ng=None,
    log_scale=True,
    adaptive_binning=True,
    save_fig=True,
    figsize=(14, 5)
):
    """
    Plot degree distributions for BA preferential, BA random, and no-growth models
    in separate subplots.
    """
    
    all_deg_pref = np.array(all_deg_pref)
    all_deg_rand = np.array(all_deg_rand)
    degrees_no_growth = np.array(degrees_no_growth)
    
    # Compute distributions for BA models 
    deg_pref, counts_pref = np.unique(all_deg_pref, return_counts=True)
    prob_pref = counts_pref / counts_pref.sum()
    
    deg_rand, counts_rand = np.unique(all_deg_rand, return_counts=True)
    prob_rand = counts_rand / counts_rand.sum()
    
    # Compute distribution for no-growth model
    if adaptive_binning:
        # Freedman-Diaconis rule: bin_width = 2 * IQR / n^(1/3)
        n_samples = len(degrees_no_growth)
        q75, q25 = np.percentile(degrees_no_growth, [75, 25])
        iqr = q75 - q25
        
        bin_width = 2 * iqr / np.cbrt(n_samples)
        bin_width = max(1, int(np.round(bin_width)))  # At least 1
        
        # Create bins
        k_min = degrees_no_growth.min()
        k_max = degrees_no_growth.max()
        bins = np.arange(k_min, k_max + bin_width, bin_width)
        
        # Compute histogram
        counts_ng, bin_edges = np.histogram(degrees_no_growth, bins=bins, density=False)
        prob_ng = counts_ng / counts_ng.sum()
        deg_ng = (bin_edges[:-1] + bin_edges[1:]) / 2  # Bin centers
        
        binning_label = f"Freedman-Diaconis (width={bin_width})"
    else:
        # Standard binning (width = 1)
        deg_ng, counts_ng = np.unique(degrees_no_growth, return_counts=True)
        prob_ng = counts_ng / counts_ng.sum()
        binning_label = "standard (width=1)"
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # === BA Preferential ===
    ax1 = axes[0]
    if log_scale:
        ax1.loglog(deg_pref, prob_pref, 'o', markersize=2, alpha=0.7, color='blue')
    else:
        ax1.plot(deg_pref, prob_pref, 'o', markersize=2, alpha=0.7, color='blue')
    
    title1 = "BA Preferential"
    if n0_ba and m0_ba and tmax_ba:
        title1 += f"\n$n_0={n0_ba}$, $m_0={m0_ba}$, $t_{{max}}={tmax_ba}$"
        title1 += f"\n$N_{{samples}}={len(all_deg_pref)}$"
    ax1.set_title(title1, fontsize=11)
    ax1.set_xlabel('degree $k$', fontsize=12)
    ax1.set_ylabel('$P(k)$', fontsize=12)
    ax1.grid(True, alpha=0.3, which='both' if log_scale else 'major')
    
    # === BA Random (expected geometric: x linear, y log) ===
    ax2 = axes[1]
    ax2.plot(deg_rand, prob_rand, 's', markersize=2, alpha=0.7, color='green')
    
    if log_scale:
        ax2.set_xscale('linear')  
        ax2.set_yscale('log')     
    else:
        ax2.set_xscale('linear')
        ax2.set_yscale('linear')
    
    title2 = "BA Random (Geometric-like)"
    if n0_ba and m0_ba and tmax_ba:
        title2 += f"\n$n_0={n0_ba}$, $m_0={m0_ba}$, $t_{{max}}={tmax_ba}$"
        title2 += f"\n$N_{{samples}}={len(all_deg_rand)}$"
    ax2.set_title(title2, fontsize=11)
    ax2.set_xlabel('degree $k$', fontsize=12)
    ax2.set_ylabel('$P(k)$', fontsize=12)
    ax2.grid(True, alpha=0.3, which='both' if log_scale else 'major')
    
    # === No-growth ===
    ax3 = axes[2]
    ax3.plot(deg_ng, prob_ng, '^', markersize=2, alpha=0.7, color='red')
    ax3.set_xscale('linear')
    ax3.set_yscale('linear')

    title3 = f"No-growth PA ({binning_label})"
    if n0_ng and m0_ng and tmax_ng:
        title3 += f"\n$n_0={n0_ng}$, $m_0={m0_ng}$, $t_{{max}}={tmax_ng}$"
        title3 += f"\n$N_{{samples}}={len(degrees_no_growth)}$"
    ax3.set_title(title3, fontsize=11)
    ax3.set_xlabel('degree $k$', fontsize=12)
    ax3.set_ylabel('$P(k)$', fontsize=12)
    ax3.grid(True, alpha=0.3, which='major')
    
    plt.tight_layout()
    
    if save_fig:
        os.makedirs('./figures', exist_ok=True)
        scale_str = 'loglog' if log_scale else 'linear'
        bin_str = 'FD' if adaptive_binning else 'standard'
        filename = f'./figures/degree_dist_separate_{scale_str}_{bin_str}.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {filename}")
    
    plt.show()

def plot_degree_distribution_with_fit(degrees, fit_results, figsize=(8, 6), point_size=8):
    """
    Plot empirical degree distribution vs fitted best model (only PMF).
    Rescales the model to match empirical mass in the fitting tail.
    """
    degrees = np.array(degrees, dtype=np.int64)
    best_model = fit_results['best_model']

    # Empirical distribution
    unique_k, counts = np.unique(degrees, return_counts=True)
    N = len(degrees)
    pmf_empirical = counts / N

    # Tail range used for fitting
    k_min_tail = fit_results.get('k_min_tail', unique_k.min())
    k_max_tail = fit_results.get('k_max_tail', unique_k.max())

    mask_tail = (unique_k >= k_min_tail) & (unique_k <= k_max_tail)
    tail_mass = pmf_empirical[mask_tail].sum()

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter(unique_k, pmf_empirical, alpha=0.7, s=point_size,
               label='Empirical', color='black',
               zorder=3, edgecolors='white', linewidths=0.3)

    k_model = None
    pmf_model = None
    label = ""

    # =====================================================================
    #                           POISSON
    # =====================================================================
    if best_model == 'poisson':
        lam = fit_results['poisson']['lambda']
        se = fit_results['poisson'].get('se_lambda')

        k_model = np.arange(k_min_tail, k_max_tail + 1)
        pmf_model = poisson.pmf(k_model, lam)

        # rescale to tail
        pmf_model /= pmf_model.sum()

        label = f"Poisson (λ={lam:.3f}"
        if se is not None:
            label += f"±{se:.3f}"
        label += ")"

    # =====================================================================
    #                           GAUSSIAN
    # =====================================================================
    if best_model == 'gaussian':
        mu = fit_results['gaussian']['mu']
        se_mu = fit_results['gaussian'].get('se_mu')
        sigma = fit_results['gaussian']['sigma']
        se_sigma = fit_results['gaussian'].get('se_sigma')

        k_model = np.arange(k_min_tail, k_max_tail + 1, dtype=float)
        log_unnorm = -0.5 * ((k_model - mu) / sigma) ** 2
        a = np.max(log_unnorm)
        pmf_model = np.exp(log_unnorm - a)
        pmf_model /= pmf_model.sum()   # discrete normalization

        pmf_model_scaled = pmf_model * tail_mass


        label = f"Gaussian (μ={mu:.3f}"
        if se_mu is not None:
            label += f"±{se_mu:.3f}"
        label += f", σ={sigma:.3f}"
        if se_sigma is not None:
            label += f"±{se_sigma:.3f}"
        label += ")"


    # =====================================================================
    #                     INFINITE DISPLACED GEOMETRIC
    # =====================================================================
    elif best_model == 'displaced_geometric':
        q = fit_results['displaced_geometric']['q']
        se = fit_results['displaced_geometric'].get('se_q')

        k_model = np.arange(k_min_tail, k_max_tail + 1)
        pmf_model = (1 - q) * (q ** (k_model - k_min_tail))

        label = f"Displaced Geometric (q={q:.3f}"
        if se is not None:
            label += f"±{se:.3f}"
        label += f", k≥{k_min_tail})"

    # =====================================================================
    #                    TRUNCATED GEOMETRIC 
    # =====================================================================
    elif best_model == 'geometric_trunc':
        q = fit_results['geometric_trunc']['q']
        se = fit_results['geometric_trunc'].get('se_q')

        k_model = np.arange(k_min_tail, k_max_tail + 1)
        L = k_max_tail - k_min_tail + 1
        pmf_model = (1 - q) * (q ** (k_model - k_min_tail)) / (1 - q**L)

        label = f"Geometric Trunc (q={q:.4f}"
        if se is not None:
            label += f"±{se:.4f}"
        label += f", k∈[{k_min_tail},{k_max_tail}])"

    # =====================================================================
    #                     ZETA (HURWITZ)
    # =====================================================================
    elif best_model == 'zeta':
        gamma = fit_results['zeta']['gamma']
        se = fit_results['zeta'].get('se_gamma')

        k_model = np.arange(k_min_tail, k_max_tail + 1)
        Z = np.sum((k_model) ** (-gamma))
        pmf_model = (k_model ** (-gamma)) / Z

        label = f"Zeta (γ={gamma:.3f}"
        if se is not None:
            label += f"±{se:.3f}"
        label += f", k≥{k_min_tail})"

    # =====================================================================
    #                    RIGHT TRUNCATED ZETA (FINITE)
    # =====================================================================
    elif best_model == 'zeta_trunc':
        gamma = fit_results['zeta_trunc']['gamma']
        se = fit_results['zeta_trunc'].get('se_gamma')

        k_min_fit = fit_results['zeta_trunc'].get('k_min_fit', k_min_tail)
        k_max_fit = fit_results['zeta_trunc'].get('k_max_fit', k_max_tail)

        k_model = np.arange(k_min_fit, k_max_fit + 1)
        Z = np.sum(k_model ** (-gamma))
        pmf_model = (k_model ** (-gamma)) / Z

        label = f"Zeta Trunc (γ={gamma:.3f}"
        if se is not None:
            label += f"±{se:.3f}"
        label += f", k∈[{k_min_fit},{k_max_fit}])"

    
    if k_model is not None and pmf_model is not None:
        pmf_model_scaled = pmf_model * tail_mass
        ax.plot(k_model, pmf_model_scaled, 'r-', linewidth=2.0,
                label=label, alpha=0.9, zorder=4)

    
    geometric_models = {'displaced_geometric', 'geometric_trunc'}
    if best_model in geometric_models:
        # Random attachment / geometric: x linear, y log
        ax.set_xscale('linear')
        ax.set_yscale('log')
    elif best_model in {'gaussian'}:
        # Gaussian: x linear, y linear
        ax.set_xscale('linear')
        ax.set_yscale('linear')
    else:
        # Power-law-like (zeta, zeta_trunc, ecc.): log-log
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax.set_xlabel('Degree k', fontsize=13, fontweight='bold')
    ax.set_ylabel('P(k)', fontsize=13, fontweight='bold')
    ax.set_title('Degree Distribution (PMF)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, framealpha=0.9, loc='best')
    ax.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    return fig
