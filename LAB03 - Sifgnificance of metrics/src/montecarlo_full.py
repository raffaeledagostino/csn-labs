import os
import numpy as np
import pandas as pd
import networkx as nx
import multiprocessing as mp
from tqdm import tqdm

from src.utils import LANG_DICT, load_stats
from src.graph_ops import er_constant, sw_model, edges_extraction

def process_language(lang, T):
    """
    Perform Monte Carlo simulations to compute p-values for a single language.
    
    Runs T simulations comparing the clustering coefficient of real dependency networks
    against two null models: Erdős-Rényi (ER) random graphs and Switching Model (SW) graphs.
    The p-value represents the fraction of simulations where the null model's clustering
    exceeds the real network's clustering coefficient.
    
    Args:
        lang (str): Two-letter ISO 639-1 language code (e.g., 'en', 'de', 'fr').
        T (int): Number of Monte Carlo simulation iterations to perform.
    
    Returns:
        tuple: A tuple containing:
            - str: Full language name from LANG_DICT.
            - float: P-value for ER model (fraction of times C_ER > C_real).
            - float: P-value for SW model (fraction of times C_SW > C_real).
            - float: Mean fraction of effective switches across all SW simulations.
            - float: Standard error of the mean for effective switches.
            - list: List of clustering coefficients from all SW simulations.
            - list: List of clustering coefficients from all ER simulations.
    """
    clustering_dataset = pd.read_csv("./data/clustering_dataset.csv")
    clustering_real = float(clustering_dataset[clustering_dataset["Language"] == LANG_DICT[lang]]["Clustering"].iloc[0])
    
    counter_er = 0
    counter_sw = 0

    store_C_sw = []
    store_C_er = []

    effective_switches_sw_iter = []
    
    for _ in range(T):  
        
        degrees = load_stats("./data/dependency_network", lang)
        nodes_er, edges_er = er_constant(degrees[0], degrees[1])
        G_er = nx.Graph()
        G_er.add_nodes_from(nodes_er)
        G_er.add_edges_from(edges_er)
        
        nodes_sw, edges_sw = edges_extraction("./data/dependency_network", lang)
        print(f"Processing {lang}: Simulation {_+1}/{T}")
        edges_sw, effective_switches_sw = sw_model(edges_sw)
        effective_switches_sw_iter.append(effective_switches_sw)
        G_sw = nx.Graph()
        G_sw.add_nodes_from(nodes_sw)
        G_sw.add_edges_from(edges_sw)
        
        C_er = nx.average_clustering(G_er)
        C_sw = nx.average_clustering(G_sw)

        store_C_sw.append(C_sw)
        store_C_er.append(C_er)

        
        if C_er > clustering_real:
            counter_er += 1
        if C_sw > clustering_real:
            counter_sw += 1

    p_value_er = counter_er / T
    p_value_sw = counter_sw / T


    return (LANG_DICT[lang], p_value_er, p_value_sw, np.mean(effective_switches_sw_iter), np.std(effective_switches_sw_iter)/np.sqrt(T), store_C_sw, store_C_er)

def main():
    """
    Main execution function for parallel Monte Carlo simulations across all languages.
    
    Orchestrates the parallel processing of multiple languages using multiprocessing.
    Each language undergoes T Monte Carlo simulations to compute p-values comparing
    real dependency networks against null models (ER and SW). Results are saved to CSV.
    
    Process:
        1. Initialize simulation parameters (T iterations per language)
        2. Create multiprocessing pool using all available CPU cores
        3. Process all languages in parallel with progress tracking
        4. Aggregate results into DataFrame
        5. Save results to CSV file
    
    Output:
        Creates './data/p_values_dataset.csv' containing:
            - Language: Full language name
            - P-value_ER: P-value for Erdős-Rényi model
            - P-value_SW: P-value for Switching Model
            - Mean_Effective_Switches_SW: Mean fraction of successful switches
            - Std_Mean_Effective_Switches_SW: Standard error of effective switches
            - C_SW_Simulations: List of clustering coefficients from SW simulations
            - C_ER_Simulations: List of clustering coefficients from ER simulations
    """
    T = 10000
    languages = list(LANG_DICT.keys())

    num_cores = mp.cpu_count()
    print(f"Running with {num_cores} cores")

    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(
            pool.starmap(
                process_language, 
                [(lang, T) for lang in languages]
            ),
            total=len(languages),
            desc="Processing languages"
        ))

    lang_p_values = pd.DataFrame(results, columns=["Language", "P-value_ER", "P-value_SW", "Mean_Effective_Switches_SW", "Std_Mean_Effective_Switches_SW", "C_SW_Simulations", "C_ER_Simulations"])

    output_path = "./data/p_values_dataset.csv"
    lang_p_values.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    main()