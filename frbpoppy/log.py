# -*- coding: future_fstrings -*-
"""Quick and dirty logging functions."""
import inspect
import sys
import numpy as np


def pprint(*s, output=True):
    """Hack to make for more informative print statements."""
    f = inspect.stack()[1][1].split('/')[-1]
    m = '{:13.13} |'.format(f)

    if output:
        print(m, *s)
    else:
        lines = []
        for e in s:
            lines.append('\n'.join([f'{m} {f}' for f in e.split('\n')]))
        return '\n'.join(lines)


def progressbar(it, prefix="", size=69, file=sys.stdout):
    """Progressbar from adapted from Stack Overflow.

    Args:
        it (generator): range of values
        prefix (str): Words displayed before the progress bar
        size (int): Display width
        file: Where to direct output

    Returns:
        type: Description of returned object.

    """
    count = len(it)
    size -= len(prefix)

    def show(j):
        x = int((size)*j/count)
        print(f'{prefix} [{"#"*x}{"."*(size-x)}] {j}/{count}')

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)

    file.flush()


def hist(parameter, bins='lin', n_bins=25, norm=True):
    """Bin up a parameter either in a lin or log space.

    Why is this not a standard option in numpy or matplotlib?

    Args:
        parameter (array): To be binned
        bins (str): Either 'lin' or 'log'
        n_bins (int): Number of bins. Can be overriden internally
        norm (bool): Whether to normalise

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
    if bins == 'lin':
        bins = n_bins
    elif bins == 'log':
        min_f = np.log10(min(parameter))
        max_f = np.log10(max(parameter))
        bins = np.logspace(min_f, max_f, n_bins)

    # Bin
    n, bins = np.histogram(parameter, bins=bins)

    if norm:
        n = n/max(n)  # Normalise
    bin_centres = (bins[:-1] + bins[1:]) / 2

    return bin_centres, n
