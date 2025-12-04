import numpy as np
from typing import List, Dict, Literal, Tuple
import multiprocessing as mp
import os
import random

def get_optimal_workers():
    """
    Auto-detect optimal number of workers based on available resources.
    Conservative: leaves 1-2 cores free for system.
    """
    n_cpus = mp.cpu_count()
    # Leave some cores free
    if n_cpus <= 4:
        return max(1, n_cpus - 1)
    elif n_cpus <= 8:
        return n_cpus - 1
    else:
        return n_cpus - 2


def _simulate_ba_preferential(n0: int, m0: int, tmax: int, arrival_times: List[int]):
    """BA with preferential attachment - optimized with stubs"""
    N_final = n0 + tmax
    degrees = np.zeros(N_final, dtype=int)

    # Initial complete graph
    stubs = []
    for u in range(n0):
        for v in range(u + 1, n0):
            stubs.extend([u, v])
            degrees[u] += 1
            degrees[v] += 1

    vertex_time_series = {i: {} for i in arrival_times}

    for t in range(1, tmax + 1):
        new_vertex = n0 + t - 1
        degrees_new = 0
        L = len(stubs)

        targets = set()
        while len(targets) < m0:
            v = stubs[np.random.randint(0, L)]
            if v != new_vertex:
                targets.add(v)

        for v in targets:
            degrees[v] += 1
            degrees_new += 1
            stubs.extend([v, new_vertex])

        degrees[new_vertex] = degrees_new

        if t in arrival_times:
            vertex_time_series[t][new_vertex] = []

        for i in arrival_times:
            for vid in vertex_time_series[i]:
                vertex_time_series[i][vid].append(degrees[vid])

    return vertex_time_series, degrees

def _simulate_ba_random(n0: int, m0: int, tmax: int, arrival_times: List[int]):
    """BA random attachment - optimized with arrays"""
    
    N_final = n0 + tmax
    degrees = np.zeros(N_final, dtype=int)
    degrees[:n0] = n0 - 1
    vertex_time_series = {i: {} for i in arrival_times}

    for t in range(1, tmax + 1):
        new_vertex = n0 + t - 1
        targets = np.array(random.sample(range(n0 + t - 1), min(m0, n0 + t - 1)))

        degrees[targets] += 1
        degrees[new_vertex] = len(targets)

        if t in arrival_times:
            vertex_time_series[t][new_vertex] = []

        for i in arrival_times:
            for vid in vertex_time_series[i]:
                vertex_time_series[i][vid].append(degrees[vid])

    return vertex_time_series, degrees


def simulate_network_single_run(n0: int, m0: int, tmax: int,
                               arrival_times: List[int],
                               model: Literal['BA_pref', 'BA_random', 'fixed_pref'] = 'BA_pref'):
    if model == 'BA_random':
        return _simulate_ba_random(n0, m0, tmax, arrival_times)
    else:
        return _simulate_ba_preferential(n0, m0, tmax, arrival_times)

def _run_single_simulation(args):
    """
    Wrapper for parallel execution. Must be picklable (top-level function).
    Returns: (run_index, time_series, degrees)
    """
    run_idx, n0, m0, tmax, arrival_times, model, seed = args

    np.random.seed(seed)

    ts, degrees = simulate_network_single_run(n0, m0, tmax, arrival_times, model)
    cleaned = {i: list(ts[i].values())[0] for i in arrival_times if i in ts}

    return run_idx, cleaned, degrees

def init_storage(arrival_times: List[int], tmax: int, n_runs: int, model: str):
    """Create storage for multi-run statistics"""
    storage = {}
    for i in arrival_times:
                length = tmax - i + 1
                storage[i] = np.zeros((n_runs, length), dtype=float)

    return storage


def compute_statistics(storage: Dict[int, np.ndarray]):
    """Compute mean, variance, std per time step"""
    stats = {}
    for i, arr in storage.items():
        stats[i] = {
            "mean": np.mean(arr, axis=0),
            "var": np.var(arr, axis=0),
            "std": np.std(arr, axis=0),
            "n": arr.shape[0]
        }
    return stats



def run_network_multi(
    n_runs: int,
    n0: int,
    m0: int,
    tmax: int,
    arrival_times: List[int],
    model: Literal['BA_pref', 'BA_random'] = 'BA_pref',
    n_jobs: int = -1,
    save_data: bool = False,
    verbose: bool = True):
    """
    Parallelized network simulation using multiprocessing.
    Parameters
    ----------
    n_runs : int
        Number of simulation runs
    n0 : int
        Initial nodes
    m0 : int
        Connections per step
    tmax : int
        Time steps
    arrival_times : List[int]
        Times to track
    model : str
        'BA_pref', 'BA_random'
    n_jobs : int
        Number of parallel workers:
        - -1: use all cores
        - -2: use all cores but one
        - positive int: use exactly that many
        Default: -1
    save_data : bool
        Save results to ./data/
    verbose : bool
        Print progress

    Returns
    -------
    storage : dict
        Raw trajectories
    stats : dict
        Statistics
    """

    # Determine number of workers
    if n_jobs == -1:
        n_workers = get_optimal_workers()
    elif n_jobs < -1:
        n_workers = max(1, mp.cpu_count() + n_jobs + 1)
    else:
        n_workers = min(n_jobs, mp.cpu_count())

    if verbose:
        print(f"\n{'='*70}")
        print(f"Running {n_runs} simulations of {model} on {n_workers} cores")
        print(f"{'='*70}")

    # Prepare arguments for parallel execution
    base_seed = np.random.randint(0, 2**31)
    args_list = [
        (i, n0, m0, tmax, arrival_times, model, base_seed + i)
        for i in range(n_runs)
    ]

    # Initialize storage
    storage = init_storage(arrival_times, tmax, n_runs, model)
    all_degrees = []
    degrees_final = None

    # Parallel execution with memory control
    with mp.Pool(processes=n_workers, maxtasksperchild=10) as pool:
        results = pool.imap(_run_single_simulation, args_list, chunksize=1)

        for run_idx, cleaned, degrees in results:
            # Store trajectories
            for i in arrival_times:
                if i in cleaned:
                    traj = cleaned[i]
                    storage[i][run_idx, :len(traj)] = traj

            # Accumulate degrees
            all_degrees.extend(degrees)
            if run_idx == n_runs - 1:
                degrees_final = degrees

            # Progress
            if verbose and (run_idx + 1) % max(1, n_runs // 10) == 0:
                print(f"  Progress: {run_idx + 1}/{n_runs} ({100*(run_idx+1)//n_runs}%)")

    if verbose:
        print(f"{'='*70}\n")

    # Compute statistics
    stats = compute_statistics(storage)

    # Save data
    if save_data:
        os.makedirs('./data', exist_ok=True)
        deg_seq = np.sort(degrees_final)[::-1]
        np.savetxt(f'./data/degree_sequence_{model}.txt', deg_seq, fmt='%d')

        all_deg = np.array(all_degrees)
        degrees_unique, counts = np.unique(all_deg, return_counts=True)
        prob = counts / counts.sum()
        np.savetxt(f'./data/degree_dist_{model}.txt',
                  np.column_stack([degrees_unique, prob]),
                  fmt=['%d', '%.8f'])

        for i in arrival_times[:4]:
            if i in stats:
                t = np.arange(i, tmax + 1)
                mean = stats[i]["mean"]
                np.savetxt(f'./data/ts_{model}_tau{i}.txt',
                          np.column_stack([t[:len(mean)], mean]),
                          fmt=['%d', '%.6f'])

    return storage, stats, all_degrees

def save_no_growth_outputs(times, k_all, vertex_indices=None, prefix="no_growth"):
    """
    Save:
    1) degree sequence at final time,
    2) time series for selected vertices,
    matching Lab 1.1 requirements.
    """

    os.makedirs("./data", exist_ok=True)

    final_degrees = k_all[:, -1]   # degree of every vertex at time T
    np.savetxt(
        f"./data/degree_sequence_{prefix}.txt",
        final_degrees,
        fmt="%d"
    )

    if vertex_indices is None:
        # Pick four vertices to track
        vertex_indices = [0, 1, 2, 3]

    for v in vertex_indices:
        data = np.column_stack([times, k_all[v]])
        np.savetxt(
            f"./data/ts_{prefix}_v{v}.txt",
            data,
            fmt=["%d", "%d"]
        )

def no_growth_pa(
    n0: int,
    m0: int,
    T: int,
    seed: int = None
):
    """
    Preferential attachment without growth (Model B-like).
    simple graph: no multiedges, no self loops.

    """

    rng = np.random.default_rng(seed)

    degrees = np.zeros(n0, dtype=int)
    adj = [set() for _ in range(n0)]
    k_all = np.zeros((n0, T+1), dtype=int)
    k_all[:, 0] = degrees

    for t in range(1, T+1):

        i = rng.integers(n0)

        for _ in range(m0):

            candidates = [j for j in range(n0) if j != i and j not in adj[i]]

            if not candidates:
                continue

            total_deg = degrees.sum()

            if total_deg == 0:
                j = rng.choice(candidates)
            else:
                cand_deg = degrees[candidates].astype(float)
                s = cand_deg.sum()

                if s > 0:
                    probs = cand_deg / s
                    j = rng.choice(candidates, p=probs)
                else:
                    j = rng.choice(candidates)

            adj[i].add(j)
            adj[j].add(i)
            degrees[i] += 1
            degrees[j] += 1

        k_all[:, t] = degrees

    return np.arange(T+1), k_all
