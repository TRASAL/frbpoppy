"""Define covenience functions."""
import matplotlib.pyplot as plt
import numpy as np
import os


def hist(parameter, bin_type='lin', n_bins=25, norm='max', edges=True,
         bins=None):
    """Bin up a parameter either in a lin or log space.

    Why is this not a standard option in numpy or matplotlib?

    Args:
        parameter (array): To be binned
        bin_type (str): Either 'lin' or 'log'
        n_bins (int): Number of bins. Can be overriden internally
        norm (bool): Whether to normalise to 'max' or 'prob'

    Returns:
        tuple: bin centers, values per bin

    """
    # Drop NaN-values
    parameter = parameter[~np.isnan(parameter)]

    # Determine number of bins
    if n_bins != 25:
        pass
    elif len(parameter) < 50:
        n_bins = 15
    elif len(parameter) > 500:
        n_bins = 50

    # Determine type of binning
    if bin_type == 'lin':
        _bins = n_bins
    elif bin_type == 'log':
        min_f = np.log10(np.min(parameter[parameter != 0]))
        max_f = np.log10(max(parameter))
        _bins = np.logspace(min_f, max_f, n_bins)

    # Allow for custom bins
    if bins is not None:
        _bins = bins

    # Allow for probability weighting
    weights = None
    if norm == 'prob':
        weights = np.ones(len(parameter)) / len(parameter)

    # Bin
    n, bin_edges = np.histogram(parameter, bins=_bins, weights=weights)

    if norm == 'max':
        n = n/max(n)  # Normalise

    # Centre bins
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Ensure there are edges on the outer bins of the histograms
    if edges:
        if bin_type == 'lin':
            bin_dif = np.diff(bins)[-1]
            bins = np.insert(bins, 0, bins[0] - bin_dif)
            bins = np.insert(bins, len(bins), bins[-1] + bin_dif)
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)
        else:
            bin_dif = np.diff(np.log10(bins))[-1]
            bins = np.insert(bins, 0, 10**(np.log10(bins[0])-bin_dif))
            bins = np.insert(bins, len(bins), 10**(np.log10(bins[-1])+bin_dif))
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)

    return bins, n


def plot_aa_style(cols=1):
    """Use a plotting style specifically for Astronomy & Astrophyiscs."""
    # Use a style file
    dir = os.path.abspath(os.path.dirname(__file__))
    style_file = os.path.join(dir, './aa.mplstyle')
    plt.style.use(style_file)

    # Check whether interactive backend is present
    if os.name == 'posix' and "DISPLAY" not in os.environ:
        # On a cluster
        plt.switch_backend('Agg')
        # Use just the basic latex package allowing \text to be used
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    if cols == 2:
        plt.rcParams["figure.figsize"] = (5.75373, 3.556)


def rel_path(path):
    """Return the relative path name for this directory."""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), path)
