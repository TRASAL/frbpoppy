"""Define an easy to use histogram function."""


def hist(parameter, bin_type='lin', n_bins=25, norm='max', edges=True):
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
        bins = n_bins
    elif bin_type == 'log':
        min_f = np.log10(np.min(parameter[parameter != 0]))
        max_f = np.log10(max(parameter))
        bins = np.logspace(min_f, max_f, n_bins)

    weights = None
    if norm == 'prob':
        weights = np.ones(len(parameter)) / len(parameter)

    # Bin
    n, bin_edges = np.histogram(parameter, bins=bins, weights=weights)

    if norm == 'max':
        n = n/max(n)  # Normalise

    # Centre bins
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Ensure there are edges on the outer bins of the histograms
    if edges:
        if bin_type == 'lin':
            bin_diff = np.diff(bins)[-1]
            bins = np.insert(bins, 0, bins[0] - bin_diff)
            bins = np.insert(bins, len(bins), bins[-1] + bin_diff)
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)
        else:
            bin_diff = np.diff(np.log10(bins))[-1]
            bins = np.insert(bins, 0, 10**(np.log10(bins[0]) - bin_diff))
            bins = np.insert(bins, len(bins), 10**(np.log10(bins[-1]) + bin_diff))
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)

    return bins, n
