from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import zeta, gammaln
from scipy.stats import poisson, geom

try:
    import numdifftools as ndt
    HAS_NUMDIFFTOOLS = True
except ImportError:
    HAS_NUMDIFFTOOLS = False
    print("WARNING: numdifftools not installed. Standard errors will not be computed.")
    print("Install with: pip install numdifftools")


def estimate_initial_parameters(t_vals, y_vals, model_name):
    """
    Smart initial parameter estimation using transformed regressions.
    """
    t_vals = np.maximum(t_vals, 1e-10)
    y_vals = np.maximum(y_vals, 1e-10)
    
    y_mean = np.mean(y_vals)
    y_max = np.max(y_vals)
    y_min = np.min(y_vals)
    t_min = t_vals[0]
    t_max = t_vals[-1]
    
    if model_name == 'Model 0':
        lr = LinearRegression(fit_intercept=False)
        lr.fit(t_vals.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        p0 = [a_init]
        
    elif model_name == 'Model 1':
        sqrt_t = np.sqrt(t_vals)
        lr = LinearRegression(fit_intercept=False)
        lr.fit(sqrt_t.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        p0 = [a_init]
        
    elif model_name == 'Model 2':
        log_t = np.log(t_vals)
        log_y = np.log(y_vals)
        lr = LinearRegression()
        lr.fit(log_t.reshape(-1, 1), log_y)
        b_init = np.clip(lr.coef_[0], 0.1, 1.5)
        a_init = max(np.exp(lr.intercept_), 1e-6)
        p0 = [a_init, b_init]
        
    elif model_name == 'Model 3':
        log_y = np.log(y_vals)
        lr = LinearRegression()
        lr.fit(t_vals.reshape(-1, 1), log_y)
        c_init = np.clip(lr.coef_[0], -1e-4, 1e-4)
        a_init = max(np.exp(lr.intercept_), 1e-6)
        p0 = [a_init, c_init]
        

    elif model_name == 'Model 4':
        d1_init = -t_min + 1.0
        x = np.log(t_vals - t_min + 1.0)
        lr = LinearRegression(fit_intercept=False)
        lr.fit(x.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        p0 = [a_init, d1_init]


    elif model_name == 'Model 0+':
        lr = LinearRegression()
        lr.fit(t_vals.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        d_init = lr.intercept_
        p0 = [a_init, d_init]
        
    elif model_name == 'Model 1+':
        sqrt_t = np.sqrt(t_vals)
        lr = LinearRegression()
        lr.fit(sqrt_t.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        d_init = lr.intercept_
        p0 = [a_init, d_init]
        
    elif model_name == 'Model 2+':
        log_t = np.log(t_vals)
        log_y = np.log(np.maximum(y_vals, 1e-10))
        lr_log = LinearRegression()
        lr_log.fit(log_t.reshape(-1, 1), log_y)
        b_init = np.clip(lr_log.coef_[0], 0.1, 1.5)
        a_init = max(np.exp(lr_log.intercept_), 1e-6)
        y_pred_no_d = a_init * np.power(t_vals, b_init)
        d_init = np.median(y_vals - y_pred_no_d)
        p0 = [a_init, b_init, d_init]
        
    elif model_name == 'Model 3+':
        mid_idx = len(t_vals) // 2
        window = len(t_vals) // 4
        t_mid = t_vals[mid_idx-window:mid_idx+window]
        y_mid = y_vals[mid_idx-window:mid_idx+window]
        log_y_mid = np.log(np.maximum(y_mid, 1e-10))
        lr = LinearRegression()
        lr.fit(t_mid.reshape(-1, 1), log_y_mid)
        c_init = np.clip(lr.coef_[0], -1e-4, 1e-4)
        a_init = max(np.exp(lr.intercept_), 1e-6)
        d_init = y_min * 0.1
        p0 = [a_init, c_init, d_init]
        
    elif model_name == 'Model 4+':
        d1_init = -t_min + 1.0
        x = np.log(t_vals - t_min + 1.0)
        lr = LinearRegression()
        lr.fit(x.reshape(-1, 1), y_vals)
        a_init = max(lr.coef_[0], 1e-6)
        d2_init = lr.intercept_
        p0 = [a_init, d1_init, d2_init]

        
    else:
        raise ValueError(f"Unknown model: {model_name}")
    
    return p0


def model_selection_unweighted(storage, stats, arrival_times, tmax, model_name='model'):
    """
    Unweighted model selection with AIC.
    
    Parameters
    ----------
    storage : dict
        Raw trajectories
    stats : dict
        Statistics
    arrival_times : list
        Arrival times
    tmax : int
        Maximum time
    model_name : str
        Network model name for display
    
    Returns
    -------
    results_df : pd.DataFrame
        Best models per tau
    all_fit_params : dict
        All fit details
    """
    
    # Model ensemble
    models = {
        'Model 0':  lambda t, a: a * t,
        'Model 1':  lambda t, a: a * np.sqrt(t),
        'Model 2':  lambda t, a, b: a * np.power(t, b),
        'Model 3':  lambda t, a, c: a * np.exp(c * t),
        'Model 4':  lambda t, a, d1: a * np.log(np.maximum(t + d1, 1e-10)),
        'Model 0+': lambda t, a, d: a * t + d,
        'Model 1+': lambda t, a, d: a * np.sqrt(t) + d,
        'Model 2+': lambda t, a, b, d: a * np.power(t, b) + d,
        'Model 3+': lambda t, a, c, d: a * np.exp(c * t) + d,
        'Model 4+': lambda t, a, d1, d2: a * np.log(np.maximum(t + d1, 1e-10)) + d2
    }
    
    param_names = {
        'Model 0': ['a'], 'Model 1': ['a'], 'Model 2': ['a', 'b'],
        'Model 3': ['a', 'c'], 'Model 4': ['a', 'd1'],
        'Model 0+': ['a', 'd'], 'Model 1+': ['a', 'd'],
        'Model 2+': ['a', 'b', 'd'], 'Model 3+': ['a', 'c', 'd'],
        'Model 4+': ['a', 'd1', 'd2']
    }
    
    results = []
    all_fit_params = {}
    
    print(f"\n{'='*90}")
    print(f"MODEL SELECTION (Unweighted): {model_name}")
    print(f"{'='*90}")
    
    for tau in arrival_times:
        mean = stats[tau]["mean"]
        t_vals = np.arange(tau, tmax + 1)
        n = len(t_vals)
        
        print(f"\ntau = {tau}")
        print("-" * 90)
        
        tau_fits = {}
        best_aic = np.inf
        best_model = None
        
        for model_key, model_func in models.items():
            try:
                p0 = estimate_initial_parameters(t_vals, mean, model_key)
                if model_key == "Model 4":
                    d1_init = -tau + 1.0
                    x = np.log(t_vals - tau + 1.0)
                    a_init = np.polyfit(x, mean, 1)[0]
                    a_init = max(a_init, 1e-8)
                    p0 = [a_init, d1_init]

                elif model_key == "Model 4+":
                    d1_init = -tau + 1.0
                    x = np.log(t_vals - tau + 1.0)
                    coeffs = np.polyfit(x, mean, 1)
                    a_init = max(coeffs[0], 1e-8)
                    d2_init = coeffs[1]
                    p0 = [a_init, d1_init, d2_init]
                                
                # ---- add bounds only for log models to make them converge better ----
                if model_key == "Model 4":
                    bounds = (
                        [0.0, -tau + 1e-6],   # a >= 0, d1 > -tau
                        [np.inf, 0.0]         # upper bounds
                    )

                elif model_key == "Model 4+":
                    bounds = (
                        [0.0, -tau + 1e-6, -np.inf],  # a>=0, d1>-tau
                        [np.inf, 0.0, np.inf]
                    )

                else:
                    bounds = (-np.inf, np.inf)

                popt, pcov = curve_fit(
                    model_func, t_vals, mean,
                    p0=p0,
                    bounds=bounds,
                    maxfev=50000,
                    method='trf'
                )

                
                if np.any(np.isnan(popt)) or np.any(np.isinf(popt)):
                    raise ValueError("NaN/Inf in parameters")
                
                y_pred = model_func(t_vals, *popt)
                
                if np.any(np.isnan(y_pred)) or np.any(np.isinf(y_pred)):
                    raise ValueError("NaN/Inf in predictions")
                
                residuals = mean - y_pred
                RSS = np.sum(residuals**2)
                RMSE = np.sqrt(RSS / n)
                
                p = len(param_names[model_key])
                AIC = n * np.log(2 * np.pi) + n * np.log(RSS / n) + n + 2 * (p + 1)
                
                perr = np.sqrt(np.diag(pcov))
                
                tau_fits[model_key] = {
                    'params': popt,
                    'param_names': param_names[model_key],
                    'uncertainties': perr,
                    'AIC': AIC,
                    'RSS': RSS,
                    'RMSE': RMSE,
                    'predictions': y_pred,
                    'func': model_func
                }
                
                if AIC < best_aic:
                    best_aic = AIC
                    best_model = model_key
                
                print(f"{model_key:12s}: AIC={AIC:10.2f}, RMSE={RMSE:.4f} ✓")
                
            except Exception as e:
                print(f"{model_key:12s}: FAILED → {str(e)[:40]}")
                tau_fits[model_key] = None
        
        all_fit_params[tau] = tau_fits
        
        if best_model is not None:
            best_info = tau_fits[best_model]
            
            print(f"\n{'='*90}")
            print(f"Best: {best_model}, AIC={best_info['AIC']:.2f}, RMSE={best_info['RMSE']:.4f}")
            
            result_row = {
                'tau': tau,
                'best_model': best_model,
                'AIC': best_info['AIC'],
                'RMSE': best_info['RMSE']
            }
            
            for name, val, err in zip(best_info['param_names'],
                                     best_info['params'],
                                     best_info['uncertainties']):
                result_row[f'{name}'] = val
                result_row[f'{name}_err'] = err
            
            results.append(result_row)
    
    results_df = pd.DataFrame(results)
    
    print(f"\n{'='*90}")
    print(f"SUMMARY: {model_name}")
    print(f"{'='*90}")
    print(results_df.to_string(index=False))
    print(f"{'='*90}\n")
    
    return results_df, all_fit_params


def model_selection_no_growth(times, k_mean, n0, m0, t_min=None):
    """
    Model selection for no-growth model using mean degree k_mean(t).
    
    Parameters
    ----------
    times : array
        Time points
    k_mean : array
        Mean degree over all nodes at each time
    n0 : int
        Number of vertices (fixed)
    m0 : int
        Edges added per time step
    t_min : int, optional
        Minimum time to include in fit (default: n0, as recommended)
    
    Returns
    -------
    results_dict : dict
        Best model info
    all_fit_params : dict
        All model fit details
    """
    
    if t_min is None:
        t_min = n0
    
    # Filter data for t >= t_min
    mask = times >= t_min
    t_vals = times[mask]
    y_vals = k_mean[mask]
    n = len(t_vals)
    
    # Theoretical value from model
    k_theory = 2 * m0 * t_vals / n0
    
    # Model ensemble
    models = {
        'Model 0':  lambda t, a: a * t,
        'Model 1':  lambda t, a: a * np.sqrt(t),
        'Model 2':  lambda t, a, b: a * np.power(t, b),
        'Model 3':  lambda t, a, c: a * np.exp(c * t),
        'Model 4':  lambda t, a, d1: a * np.log(np.maximum(t + d1, 1e-10)),
        'Model 0+': lambda t, a, d: a * t + d,
        'Model 1+': lambda t, a, d: a * np.sqrt(t) + d,
        'Model 2+': lambda t, a, b, d: a * np.power(t, b) + d,
        'Model 3+': lambda t, a, c, d: a * np.exp(c * t) + d,
        'Model 4+': lambda t, a, d1, d2: a * np.log(np.maximum(t + d1, 1e-10)) + d2
    }
    
    param_names = {
        'Model 0': ['a'], 'Model 1': ['a'], 'Model 2': ['a', 'b'],
        'Model 3': ['a', 'c'], 'Model 4': ['a', 'd1'],
        'Model 0+': ['a', 'd'], 'Model 1+': ['a', 'd'],
        'Model 2+': ['a', 'b', 'd'], 'Model 3+': ['a', 'c', 'd'],
        'Model 4+': ['a', 'd1', 'd2']
    }
        
    
    all_fit_params = {}
    best_aic = np.inf
    best_model = None
    
    print(f"\n{'='*90}")
    print(f"MODEL SELECTION: No-growth PA model (n0={n0}, m0={m0})")
    print(f"Fitting range: t ∈ [{t_min}, {times[-1]}], n_points = {n}")
    print(f"Theoretical prediction: k(t) = 2*m0*t/n0 = {2*m0/n0:.6f}*t")
    print(f"{'='*90}\n")
    
    for model_key, model_func in models.items():
        try:
            p0 = estimate_initial_parameters(t_vals, y_vals, model_key)
            
            # Use higher maxfev for complex models
            max_iterations = 50000
            
            popt, pcov = curve_fit(
                model_func, t_vals, y_vals,
                p0=p0,
                maxfev=max_iterations,
                method='lm'
            )
            
            if np.any(np.isnan(popt)) or np.any(np.isinf(popt)):
                raise ValueError("NaN/Inf in parameters")
            
            y_pred = model_func(t_vals, *popt)
            
            if np.any(np.isnan(y_pred)) or np.any(np.isinf(y_pred)):
                raise ValueError("NaN/Inf in predictions")
            
            residuals = y_vals - y_pred
            RSS = np.sum(residuals**2)
            RMSE = np.sqrt(RSS / n)
            
            p = len(param_names[model_key])
            AIC = n * np.log(2 * np.pi) + n * np.log(RSS / n) + n + 2 * (p + 1)
            
            perr = np.sqrt(np.diag(pcov))
            
            all_fit_params[model_key] = {
                'params': popt,
                'param_names': param_names[model_key],
                'uncertainties': perr,
                'AIC': AIC,
                'RSS': RSS,
                'RMSE': RMSE,
                'predictions': y_pred,
                'func': model_func,
                't_vals': t_vals
            }
            
            if AIC < best_aic:
                best_aic = AIC
                best_model = model_key
            
            # Print parameter values
            param_str = ", ".join([f"{name}={val:.6f}±{err:.6f}" 
                                  for name, val, err in zip(param_names[model_key], popt, perr)])
            print(f"{model_key:12s}: AIC={AIC:10.2f}, RMSE={RMSE:.6f}")
            print(f"              {param_str}")
            
        except Exception as e:
            error_msg = str(e)
            if len(error_msg) > 60:
                error_msg = error_msg[:60] + "..."
            print(f"{model_key:12s}: FAILED → {error_msg}")
            all_fit_params[model_key] = None
    
    print(f"\n{'='*90}")
    if best_model is not None:
        best_info = all_fit_params[best_model]
        print(f"BEST MODEL: {best_model}")
        print(f"AIC = {best_info['AIC']:.2f}, RMSE = {best_info['RMSE']:.6f}")
        
        param_str = ", ".join([f"{name}={val:.6f}±{err:.6f}" 
                              for name, val, err in zip(best_info['param_names'], 
                                                       best_info['params'], 
                                                       best_info['uncertainties'])])
        print(f"Parameters: {param_str}")
        
        # Compare with theoretical value for Model 0
        if best_model == 'Model 0':
            a_fitted = best_info['params'][0]
            a_theory = 2 * m0 / n0
            print(f"\nTheoretical slope: a = 2*m0/n0 = {a_theory:.6f}")
            print(f"Fitted slope:      a = {a_fitted:.6f}")
            print(f"Relative difference: {abs(a_fitted - a_theory)/a_theory * 100:.2f}%")
    else:
        print("WARNING: No model could be fitted successfully!")
    
    print(f"{'='*90}\n")
    
    results_dict = {
        'best_model': best_model,
        'best_info': all_fit_params[best_model] if best_model else None,
        'n0': n0,
        'm0': m0,
        't_min': t_min,
        't_max': times[-1],
        'theoretical_slope': 2 * m0 / n0
    }
    
    return results_dict, all_fit_params



def fit_degree_distributions(
    degrees,
    model_names=None,
    verbose=True,
    compute_se=True,
    k_min_tail=None,
    k_max_tail=None
):
    """
    Fit discrete distributions to degree sequence using MLE.
    Models are fitted on the range [k_min_tail, k_max_tail].
    If not provided, the full range of empirical degrees is used.
    """
    degrees = np.array(degrees, dtype=np.int64)
    N = len(degrees)

    # Full distribution
    unique_k_full, counts_full = np.unique(degrees, return_counts=True)
    k_min_full = unique_k_full.min()
    k_max_full = unique_k_full.max()

    # Range for the fit
    if k_min_tail is None:
        k_min_tail = k_min_full
    if k_max_tail is None:
        k_max_tail = k_max_full

    mask_tail = (unique_k_full >= k_min_tail) & (unique_k_full <= k_max_tail)
    unique_k_tail = unique_k_full[mask_tail]
    counts_tail = counts_full[mask_tail]
    mask_tail_z = (unique_k_full >= k_min_tail)
    unique_k_z = unique_k_full[mask_tail_z]
    counts_tail_z = counts_full[mask_tail_z]

    if len(unique_k_tail) < 5:
        raise ValueError("Too few distinct degrees in tail range.")

    # avg value for initialization 
    mean_deg_tail = np.sum(unique_k_tail * counts_tail) / np.sum(counts_tail)

    if model_names is None:
        model_names = ['poisson', 'displaced_geometric', 'geometric_trunc', 'zeta', 'zeta_trunc', 'gaussian']

    compute_se = compute_se and HAS_NUMDIFFTOOLS

    # ============ NLL FUNCTIONS (TUTTI SULLA TAIL) ============

    def nll_poisson(params):
        """Poisson on tail [k_min_tail, k_max_tail]"""
        lambda_param = params[0]
        if lambda_param <= 0:
            return np.inf
        log_pmf = (unique_k_tail * np.log(lambda_param)
                   - lambda_param
                   - gammaln(unique_k_tail + 1))
        return -np.sum(counts_tail * log_pmf)
    
    
        
    def nll_geometric_trunc(params):
        """
        Truncated geometric on k ∈ [k_min_tail, k_max_tail].
        P(k) = (1-q)/(1-q^L) * q^(k - k_min_tail)
        """
        q = params[0]
        if q <= 0 or q >= 1:
            return np.inf

        L = k_max_full - k_min_full + 1

        try:
            log_norm = np.log(1 - q) - np.log(1 - q**L)
            log_pmf = log_norm + (unique_k_tail - k_min_tail) * np.log(q)
            return -np.sum(counts_tail * log_pmf)
        except FloatingPointError:
            return np.inf

        
        
    def nll_displaced_geometric(params):
        """
        Infinite displaced geometric with support starting at k_min_tail:
        P(k) = (1-q) * q^(k - k_min_tail)
        """
        q = params[0]
        if q <= 0 or q >= 1:
            return np.inf

        try:
            log_pmf = np.log(1 - q) + (unique_k_tail - k_min_full) * np.log(q)
            return -np.sum(counts_tail * log_pmf)
        except FloatingPointError:
            return np.inf
        
    def nll_gaussian(params):
        """
        Discrete Gaussian on k ∈ [k_min_tail, k_max_tail].
        P(k) = exp(-(k-mu)^2/(2*sigma^2)) / Z
        """
        mu, sigma = params
        if sigma <= 0:
            return np.inf

        k = unique_k_tail.astype(float)
        log_unnorm = -0.5 * ((k - mu) / sigma) ** 2
        k_range = np.arange(k_min_tail, k_max_tail + 1, dtype=float)
        log_unnorm_range = -0.5 * ((k_range - mu) / sigma) ** 2
        
        a = np.max(log_unnorm_range)
        Z = np.log(np.sum(np.exp(log_unnorm_range - a))) + a

        log_pmf = log_unnorm - Z

        return -np.sum(counts_tail * log_pmf)

    

    def nll_zeta(params):
        """
        Left-truncated Zeta (Hurwitz) distribution.
        Support: k = k_min, k_min+1, ..., infinity.

        p(k) = k^{-gamma} / zeta(gamma, k_min)
        """
        gamma = params[0]
        if gamma <= 1:
            return np.inf  # divergence

        try:
            # Normalization constant = Hurwitz zeta(gamma, k_min)
            Z = zeta(gamma, k_min_tail)
            if not np.isfinite(Z) or Z <= 0:
                return np.inf

            # log p(k) = -gamma*log(k) - log(Z)
            log_pmf = -gamma * np.log(unique_k_z) - np.log(Z)

            return -np.sum(counts_tail_z * log_pmf)

        except Exception:
            return np.inf

    def nll_zeta_trunc(params):
        """Zeta truncated at k_max_tail."""
        gamma = params[0]
        if gamma <= 0:
            return np.inf
        try:
            k_range = np.arange(k_min_tail, k_max_tail + 1)
            norm_const = np.sum(k_range ** (-gamma))
            if not np.isfinite(norm_const) or norm_const <= 0:
                return np.inf
            log_pmf = -gamma * np.log(unique_k_tail) - np.log(norm_const)
            return -np.sum(counts_tail * log_pmf)
        except Exception:
            return np.inf

    # ==================== STANDARD ERRORS =======================

    def compute_standard_errors(nll_func, params_mle):
        if not compute_se:
            return None
        try:
            hessian_ndt = ndt.Hessian(nll_func, full_output=False)
            H = hessian_ndt(params_mle)
            cov_matrix = np.linalg.inv(H)
            se = np.sqrt(np.diag(cov_matrix))
            return se
        except Exception:
            return None

    # ========================== FIT MODELS ================================

    results = {}
    aic_values = {}

    if 'poisson' in model_names:
        res = minimize(nll_poisson, x0=[mean_deg_tail],
                       method='L-BFGS-B', bounds=[(1e-9, None)])
        if res.success:
            se_lambda = compute_standard_errors(nll_poisson, res.x)
            results['poisson'] = {
                'lambda': res.x[0],
                'se_lambda': se_lambda[0] if se_lambda is not None else None,
                'nll': res.fun,
                'AIC': 2 * 1 + 2 * res.fun,
                'n_params': 1
            }
            aic_values['poisson'] = results['poisson']['AIC']

    if 'geometric_trunc' in model_names:
        q_init = (mean_deg_tail - k_min_tail) / (mean_deg_tail - k_min_tail + 1)
        q_init = np.clip(q_init, 1e-9, 1 - 1e-9)
        
        res = minimize(nll_geometric_trunc, x0=[q_init],
                       method='L-BFGS-B', bounds=[(1e-9, 1 - 1e-9)])
        if res.success:
            se_q = compute_standard_errors(nll_geometric_trunc, res.x)
            results['geometric_trunc'] = {
                'q': res.x[0],
                'se_q': se_q[0] if se_q is not None else None,
                'nll': res.fun,
                'AIC': 2 * 1 + 2 * res.fun, 
                'n_params': 1,
                'k_max_fit': k_max_tail
            }
            aic_values['geometric_trunc'] = results['geometric_trunc']['AIC']
    
    if 'displaced_geometric' in model_names:
        q_init = (mean_deg_tail - k_min_tail) / (mean_deg_tail - k_min_tail + 1)
        q_init = np.clip(q_init, 1e-9, 1 - 1e-9)
        
        res = minimize(nll_displaced_geometric, x0=[q_init],
                       method='L-BFGS-B', bounds=[(1e-9, 1 - 1e-9)])
        if res.success:
            se_q = compute_standard_errors(nll_displaced_geometric, res.x)
            results['displaced_geometric'] = {
                'q': res.x[0],
                'se_q': se_q[0] if se_q is not None else None,
                'nll': res.fun,
                'AIC': 2 * 1 + 2 * res.fun, 
                'n_params': 1
            }
            aic_values['displaced_geometric'] = results['displaced_geometric']['AIC']


    if 'gaussian' in model_names:
        mu0 = mean_deg_tail
        sigma0 = np.std(np.repeat(unique_k_tail, counts_tail))
        sigma0 = max(sigma0, 1e-6)

        res = minimize(
            nll_gaussian,
            x0=[mu0, sigma0],
            method='L-BFGS-B',
            bounds=[(k_min_tail, k_max_tail), (1e-6, None)]
        )

        if res.success:
            se = compute_standard_errors(nll_gaussian, res.x)
            se_mu = se[0] if se is not None else None
            se_sigma = se[1] if se is not None else None

            results['gaussian'] = {
                'mu': res.x[0],
                'sigma': res.x[1],
                'se_mu': se_mu,
                'se_sigma': se_sigma,
                'nll': res.fun,
                'AIC': 2 * 2 + 2 * res.fun,  # 4 parameters: mu, sigma, k_min_fit, k_max_fit
                'n_params': 2
            }
            aic_values['gaussian'] = results['gaussian']['AIC']


    if 'zeta' in model_names:
        res = minimize(nll_zeta, x0=[3.0],
                       method='L-BFGS-B', bounds=[(1.01, 5.0)])
        if res.success:
            se_gamma = compute_standard_errors(nll_zeta, res.x)
            results['zeta'] = {
                'gamma': res.x[0],
                'se_gamma': se_gamma[0] if se_gamma is not None else None,
                'nll': res.fun,
                'AIC': 2 * 1 + 2 * res.fun, # 2 parameters: gamma, k_min_fit
                'n_params': 1,
                'k_min_fit': k_min_tail
            }
            aic_values['zeta'] = results['zeta']['AIC']

    if 'zeta_trunc' in model_names:
        res = minimize(nll_zeta_trunc, x0=[3.0],
                       method='L-BFGS-B', bounds=[(0.1, 5.0)])
        if res.success:
            se_gamma = compute_standard_errors(nll_zeta_trunc, res.x)
            results['zeta_trunc'] = {
                'gamma': res.x[0],
                'se_gamma': se_gamma[0] if se_gamma is not None else None,
                'k_min_fit': k_min_tail,
                'k_max_fit': k_max_tail,
                'nll': res.fun,
                'AIC': 2 * 1 + 2 * res.fun, # 3 parameters: gamma, k_min_fit, k_max_fit
                'n_params': 1
            }
            aic_values['zeta_trunc'] = results['zeta_trunc']['AIC']

    if len(aic_values) == 0:
        print("ERROR: No models fitted successfully!")
        return None

    best_model = min(aic_values, key=aic_values.get)
    best_aic = aic_values[best_model]
    for model in results:
        results[model]['delta_AIC'] = results[model]['AIC'] - best_aic

    N_used = int(np.sum(counts_tail))
    
    results['best_model'] = best_model
    results['N'] = N
    results['N_used'] = N_used
    results['k_min_full'] = k_min_full
    results['k_max_full'] = k_max_full
    results['k_min_tail'] = k_min_tail
    results['k_max_tail'] = k_max_tail
    results['mean_degree'] = mean_deg_tail

    if verbose:
        print("\n" + "="*95)
        print("MAXIMUM LIKELIHOOD ESTIMATION - FAIR MODEL SELECTION")
        print(f"Total N = {N} | Used N = {N_used} in range k ∈ [{k_min_tail}, {k_max_tail}]")
        print("="*95)
        print(f"{'Model':<18} {'Parameters':<40} {'AIC':<12} {'ΔAIC':<12}")
        print("-"*95)

        for model in model_names:
            if model not in results:
                continue
            r = results[model]

            # ---------- PARAMETER FORMATTING ----------
            if model == 'poisson':
                ps = f"λ = {r['lambda']:.4f}"
                if r.get('se_lambda') is not None:
                    ps += f" (±{r['se_lambda']:.4f})"

            elif model == 'displaced_geometric':  
                # Infinite displaced geometric: p(k) = (1-q)/q * q^k
                ps = f"q = {r['q']:.4f}"
                if r.get('se_q') is not None:
                    ps += f" (±{r['se_q']:.4f})"

            elif model == 'geometric_trunc':
                # Truncated geometric: P(k) = (1-q)/(1-q^L) * q^(k-k_min)
                ps = f"q = {r['q']:.4f}"
                if r.get('se_q') is not None:
                    ps += f" (±{r['se_q']:.4f})"
                ps += f", k <= {r['k_max_fit']}]"

            elif model == 'zeta':
                # Left-truncated zeta: P(k) = k^{-γ} / ζ(γ, k_min)
                ps = f"γ = {r['gamma']:.4f}"
                if r.get('se_gamma') is not None:
                    ps += f" (±{r['se_gamma']:.4f})"
                ps += f", k_min={r['k_min_fit']}"

            elif model == 'gaussian':
                ps  = f"μ = {r['mu']:.2f}"
                if r.get('se_mu') is not None:
                    ps += f" (±{r['se_mu']:.2f})"
                ps += f", σ = {r['sigma']:.2f}"
                if r.get('se_sigma') is not None:
                    ps += f" (±{r['se_sigma']:.2f})"

            elif model == 'zeta_trunc':
                # Fully truncated zeta: k_min <= k <= k_max
                ps = f"γ = {r['gamma']:.4f}"
                if r.get('se_gamma') is not None:
                    ps += f" (±{r['se_gamma']:.4f})"
                ps += f", k∈[{r['k_min_fit']},{r['k_max_fit']}]"

            else:
                ps = " "

            marker = " ★" if model == best_model else ""
            print(f"{model:<18} {ps:<40} {r['AIC']:<12.2f} {r['delta_AIC']:<12.2f}{marker}")

        print("="*95 + "\n")


    return results
