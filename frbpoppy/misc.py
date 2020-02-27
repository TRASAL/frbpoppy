"""Convenience functions."""
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


def hist(parameter, bin_type='lin', n_bins=25, norm='max', edges=True,
         bins=None):
    """Bin up a parameter either in a lin or log space.

    Why is this not a standard option in numpy or matplotlib?

    Args:
        parameter (array): To be binned
        bin_type (str): Either 'lin' or 'log' or 'ln'
        n_bins (int): Number of bins. Can be overriden internally
        norm (bool): Whether to normalise to 'max' or 'prob' or none

    Returns:
        tuple: bin centers, values per bin

    """
    if isinstance(parameter, list):
        parameter = np.array(parameter)

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
    elif bin_type == 'ln':
        min_f = np.log(np.min(parameter[parameter != 0]))
        max_f = np.log(max(parameter))
        _bins = np.logspace(min_f, max_f, n_bins, base=np.e)

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


def calc_lognormal_input(lognormal_mean, lognormal_std):
    """Calculate the mean and std of a lognormal distribution.

    See
    https://stackoverflow.com/questions/48014712/get-lognormal-random-number
    -given-log10-mean-and-log10-standard-deviation/48016650#48016650
    for more info
    """
    normal_std = np.sqrt(np.log(1 + (lognormal_std/lognormal_mean)**2))
    normal_mean = np.log(lognormal_mean) - normal_std**2 / 2
    return normal_mean, normal_std
