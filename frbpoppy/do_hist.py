"""Bin up Pandas DataFrames into histograms ready for Bokeh plotting."""
import numpy as np
import pandas as pd


def histogram(dfs, n_bins=50, log=False, mc=False, cum=False):
    """
    Quick function to 'histogram' each column of each dataframe.

    Args:
        dfs (list): List of dataframes
        n_bins (int): Number of bins
        log (bool): Bin in log space or not
        mc (bool): Whether running for a Monte Carlo plot
    Returns:
        hists (list): List of histogramed dataframes
    """
    np.seterr(divide='ignore', invalid='ignore')
    cols = ['ra', 'dec', 'dist_co', 'gb', 'gl', 'gx', 'gy', 'gz', 'z', 'dm',
            'dm_host', 'dm_igm', 'dm_mw', 'lum_bol', 'si', 'w_arr', 'w_int',
            'fluence', 's_peak', 'snr', 't_dm', 'T_sky', 'T_sys',
            'w_eff', 'time']
    # Determine bin limits
    limits = {}

    for c in cols:

        low = 1e99
        high = -1e99

        for df in dfs:
            if c not in df:
                continue

            col = pd.to_numeric(df[c], errors='coerce')
            dlow = col.min()
            dhigh = col.max()

            if isinstance(dlow, str) or isinstance(dhigh, str):
                continue

            if dlow < low:
                low = dlow

            if dhigh > high:
                high = dhigh

        limits[c] = (low, high)

    hists = []

    # Bin each dataframe
    for df in dfs:

        hist = pd.DataFrame()

        for c in cols:

            low = limits[c][0]
            high = limits[c][1]

            if low == 1e99 or high == -1e99:
                continue

            if c not in df:
                continue

            if df[c].nunique() == 1 and df[c].iloc[0] == 'None':
                continue

            if log:
                if low == 0:
                    low = 1e-3
                if high == 0:
                    high = 1e-3
                bins = np.geomspace(low, high, num=n_bins)
            else:
                bins = np.linspace(low, high, n_bins)

            col = df[c].apply(pd.to_numeric, errors='coerce')
            col = col.dropna()
            h, _ = np.histogram(col, bins=bins)

            # Normalise
            h = h/sum(h)

            # Cumulative
            if cum:
                h = [sum(h[i:]) for i in range(len(h))]

            hist[f'{c}'] = pd.Series(h)
            hist[f'{c}_x'] = pd.Series(bins[:-1])

        hist['population'] = df['population'].iloc[0]

        hists.append(hist)

    # Ugly, but will have to do
    # Make the bottom of the bins line up
    for c in cols:
        min_c = 1e99
        for hist in hists:
            if c in hist:
                h = hist[c]
                m = np.min(h[h.to_numpy().nonzero()[0]])
                if m < min_c:
                    min_c = m

        for hist in hists:
            if c in hist:
                hist.loc[hist[c] == 0., c] = 10**(int(np.log10(min_c)) - 1)

    return hists
