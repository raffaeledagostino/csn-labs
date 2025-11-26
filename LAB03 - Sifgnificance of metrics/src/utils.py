import numpy as np
from matplotlib.ticker import MaxNLocator, FuncFormatter

LANG_DICT = {
    "en": "english",
    "ar": "arabic",
    "cs": "czech",
    "de": "german",
    "es": "spanish",
    "fi": "finnish",
    "fr": "french",
    "gl": "galician",
    "hi": "hindi",
    "id": "indonesian",
    "is": "icelandic",
    "it": "italian",
    "ja": "japanese",
    "ko": "korean",
    "pl": "polish",
    "pt": "portuguese",
    "ru": "russian",
    "sv": "swedish",
    "th": "thai",
    "tr": "turkish",
    "zh": "chinese",
}

def load_stats(
    base_path: str = None, lang: str = None, fixed_path: str = None
) -> tuple[int, int]:
    """
    Load basic network statistics from a dependency network adjacency file.
    
    Reads the first line of the adjacency list file which contains metadata
    about the network structure: the total number of nodes and edges.
    
    Args:
        base_path (str, optional): Base directory path for language-specific files.
        lang (str, optional): Language code (e.g., 'en', 'de') used to construct 
                             the file path via LANG_DICT.
        fixed_path (str, optional): Direct file path. If provided, overrides 
                                   base_path and lang parameters.
    
    Returns:
        tuple[int, int]: A tuple containing:
            - int: Total number of nodes in the network.
            - int: Total number of edges in the network.
    """
    if fixed_path:
        path = fixed_path
    else:
        path = f"{base_path}/{LANG_DICT[lang]}_dependency_network_adj.txt"

    with open(path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
        nodes, edges = map(int, first_line.split())
        return (nodes, edges)

def format_pvalue(val, T=10000):
    """
    Format p-values with uncertainty estimates for display.
    
    Applies special formatting rules for extreme p-values and calculates
    standard error for intermediate values based on binomial proportion variance.
    
    Args:
        val (float): P-value to format, expected to be in the range [0, 1].
        T (int, optional): Number of Monte Carlo simulations used to compute 
                          the p-value. Defaults to 1000.
    
    Returns:
        str: Formatted p-value string:
            - ">1-10^-4" for val = 1.0 (all simulations extreme)
            - "<10^-4" for val = 0.0 (no simulations extreme)
            - "value ± stderr" for intermediate values, where stderr = sqrt(val*(1-val)/T)
    """
    if val == 1.0:
        return ">1-10^-4"
    elif val == 0.0:
        return "<10^-4"
    else:
        std_error = np.sqrt(val * (1 - val) / T)
        return f"{val:.4f} ± {std_error:.4f}"

def format_xaxis_clean(ax, sample_values):
    """
    Format x-axis with scientific notation and clean tick labels.
    
    Automatically scales axis values by appropriate power of 10 and displays
    the scale factor separately to reduce label clutter. Useful for histograms
    with very small or very large values.
    
    Args:
        ax (matplotlib.axes.Axes): The axis object to format.
        sample_values (array-like): Sample data values used to determine 
                                   appropriate scaling factor.
    
    Side Effects:
        - Sets major tick locator to MaxNLocator with up to 5 bins
        - Formats tick labels as scaled values (e.g., "0.50" instead of "5.0e-3")
        - Adds scale factor annotation (e.g., "×10^-3") near the axis if exp != 0
        - Hides default offset text
    """
    vals = np.asarray(sample_values, dtype=float)
    vmax = np.nanmax(np.abs(vals))
    
    if not np.isfinite(vmax) or vmax == 0:
        exp = 0
    else:
        exp = int(np.floor(np.log10(vmax)))
    
    scale = 10.0 ** exp
    
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))
    
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x/scale:.2f}"))
    
    ax.xaxis.get_offset_text().set_visible(False)
    
    if exp != 0:  
        ax.text(
            1.02, -0.00,
            rf"$\times\!10^{{{exp}}}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
            color="#333333"
        )
