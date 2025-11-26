import os
import numpy as np
from src.utils import LANG_DICT
import random

def er_constant(N, E):
    """
    Generate an Erdős-Rényi random graph with a fixed number of edges.
    
    Creates a random graph by repeatedly sampling pairs of nodes until
    the desired number of edges is reached. Self-loops and double edges are not allowed.
    
    Args:
        N (int): Number of nodes in the graph.
        E (int): Number of edges to generate.
    
    Returns:
        tuple: A tuple containing:
            - list: List of node IDs from 1 to N.
            - set: Set of edge tuples (i, j) where i and j are node IDs.
    """
    nodes = list(range(1, N + 1))
    edges = set()
    while len(edges) < E:
        i, j = np.random.randint(1, N+1, size=2)
        if i != j:
            edge = tuple(sorted([i, j]))
            edges.add(edge)
    return (nodes, edges)

def sw_model(set_edges):
    """
    Apply the switching model to randomize a graph while preserving degree distribution.
    
    Performs Q * |E| edge rewiring attempts, where Q = log(|E|). Each attempt selects
    two random edges (u,v) and (x,y) and tries to replace them with (u,x) and (v,y)
    if the new edges don't already exist in the graph.
    
    Args:
        set_edges (set): Set of edge tuples representing the graph edges.
    
    Returns:
        tuple: A tuple containing:
            - set: Updated set of edges after switching.
            - float: Fraction of successful switches performed.
    """

    effective_switches = 0
    edges_list = list(set_edges)
    num_edges = len(edges_list)
    Q = int(round(np.log(num_edges)))
    total_attempts = Q * num_edges

    for _ in range(total_attempts):

        i1, i2 = random.sample(range(num_edges), 2)
        edge1 = edges_list[i1]
        edge2 = edges_list[i2]

        u, v = edge1
        x, y = edge2

        if u in (x, y) or v in (x, y):
            continue

        random_choice = random.random()

        if random_choice < 0.5:
            new_edge1 = tuple(sorted((u, x)))
            new_edge2 = tuple(sorted((v, y)))
        else:
            new_edge1 = tuple(sorted((u, y)))
            new_edge2 = tuple(sorted((v, x)))

        if new_edge1 in set_edges or new_edge2 in set_edges:
            continue
        
        set_edges.add(new_edge1)
        set_edges.add(new_edge2)
        set_edges.discard(edge1)
        set_edges.discard(edge2)

        edges_list[i1] = new_edge1
        edges_list[i2] = new_edge2

        effective_switches += 1
           
    return (set_edges, effective_switches/total_attempts)



def parse_ud_conllu(data_folder: str = "./data/PUD", lang: str = "en"):
    """
    Parse Universal Dependencies CoNLL-U format files to extract dependency graph structure.
    
    Reads a CoNLL-U formatted file and extracts syntactic dependency relationships
    between words, creating a directed graph where edges represent head-dependent relations.
    
    Args:
        data_folder (str): Path to the directory containing CoNLL-U files. 
                          Defaults to "./data/PUD".
        lang (str): Language code (e.g., 'en', 'de', 'fr'). Defaults to "en".
    
    Returns:
        tuple: A tuple containing:
            - set: Set of unique words (nodes) in the dependency graph.
            - set: Set of edge tuples (parent_word, child_word) representing dependencies.
    """
    with open(data_folder + f"/{lang}_pud-ud-test.conllu", "r", encoding="utf-8") as f:

        lines = f.readlines()

        id_to_word_map = {}
        sentence_id = 0

        for line in lines:
            if line.startswith("#"):
                continue
            if len(line) < 7:
                sentence_id += 1
                continue

            words = line.split("\t")
            if "-" in words[0] or "." in words[0]:
                continue
            
            if words[2] == "_":
                id_to_word_map[(sentence_id, int(words[0]))] = words[1]
            else:
                id_to_word_map[(sentence_id, int(words[0]))] = words[2]

        edges = set()
        sentence_id = 0
        for line in lines:
            if line.startswith("#"):
                continue
            if len(line) < 7:
                sentence_id += 1
                continue

            words = line.split("\t")
            if "-" in words[0] or "." in words[0]:
                continue

            parent_id = int(words[6])
            if parent_id != 0:
                child_id = int(words[0])
                parent_word = id_to_word_map[(sentence_id, parent_id)]
                child_word = id_to_word_map[(sentence_id, child_id)]
                if parent_word != child_word:
                    edges.add((parent_word, child_word))

        nodes = set()
        for parent, child in edges:
            nodes.add(parent)
            nodes.add(child)

        return (nodes, edges)

def edges_extraction(base_path: str = None, lang: str = None, fixed_path: str = None) -> set[tuple[str, str]]:
    """
    Load dependency network edges from a text file.
    
    Reads a file containing edge tuples (one per line) and returns them as a set.
    Edges are normalized by sorting to ensure consistency. The first line (header)
    is skipped.
    
    Args:
        base_path (str, optional): Base directory path for language-specific files.
        lang (str, optional): Language code used to construct the file path.
        fixed_path (str, optional): Direct file path. If provided, overrides base_path and lang.
    
    Returns:
        tuple: A tuple containing:
            - set: Set of unique nodes (words) in the network.
            - set: Set of normalized edge tuples (word1, word2) where word1 <= word2.
    """
    import ast

    if fixed_path:
        path = fixed_path
    else:
        path = f"{base_path}/{LANG_DICT[lang]}_dependency_network_edges.txt"

    edges = set()
    nodes = set()
    with open(path, "r", encoding="utf-8") as f:
        _ = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                tup = ast.literal_eval(line)
            except Exception:
                continue
            if isinstance(tup, tuple) and len(tup) == 2:
                normalized_edge = tuple(sorted(tup))
                edges.add(normalized_edge)
                nodes.add(tup[0])
                nodes.add(tup[1])
    return nodes, set(edges)

def adj_extraction(base_path: str = None, lang: str = None, fixed_path: str = None) -> dict:
    """
    Load adjacency list representation of a graph from a text file.
    
    Parses a file containing adjacency list entries in the format ([node], [neighbors])
    and constructs a dictionary mapping each node to its list of neighbors. The first
    line (header) is skipped.
    
    Args:
        base_path (str, optional): Base directory path for language-specific files.
        lang (str, optional): Language code used to construct the file path.
        fixed_path (str, optional): Direct file path. If provided, overrides base_path and lang.
    
    Returns:
        dict: Dictionary mapping each node to a list of its neighboring nodes.
              Format: {node_name: [neighbor1, neighbor2, ...]}.
    """
    import ast

    if fixed_path:
        path = fixed_path
    else:
        path = f"{base_path}/{LANG_DICT[lang]}_dependency_network_adj.txt"

    adj_dict = dict()
    with open(path, "r", encoding="utf-8") as f:
        _ = f.readline()
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                parsed = ast.literal_eval(line)
                
                if isinstance(parsed, tuple) and len(parsed) == 2:
                    node_list = parsed[0]
                    neighbors_list = parsed[1]
                    
                    if isinstance(node_list, list) and len(node_list) == 1:
                        node_name = node_list[0]
                        adj_dict[node_name] = neighbors_list
                        
            except Exception as e:
                print(f"Error parsing line: {line[:80]}..., error: {e}")
                continue
    
    return adj_dict
