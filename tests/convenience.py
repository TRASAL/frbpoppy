"""Define covenience functions."""
import matplotlib.pyplot as plt
import os


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
