"""Bin up Pandas DataFrames into histograms ready for Bokeh plotting."""
import numpy as np
import pandas as pd


def histogram(dfs, n_bins=50, log=False, mc=False):
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
    # Determine all column names
    columns = [df.columns.values for df in dfs]
    cols = list(min(columns, key=len))
    for c in ['id', 'in_par', 'frbcat']:
        if c in cols:
            cols.remove(c)

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

        hist = pd.DataFrame(np.nan,
                            index=np.arange(n_bins-1),
                            columns=['empty'])

        if not mc:
            hist['color'] = df['color'][0]
            hist['population'] = df['population'][0]

        if log or not mc:
            hist['bottom'] = 10**(round(np.log10(1/len(df))) - 1)
        else:
            hist['bottom'] = 0

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
            h = [e/h.sum() for e in h]

            if mc:
                hist['id'] = df['id'].iloc[0]
                hist['in_par'] = df['in_par'].iloc[0]
                hist['out_par'] = c
                hist['survey'] = df['survey'].iloc[0]
                hist['frbs'] = len(col)  # number of frbs in population
                hist['frbcat'] = df['frbcat'].iloc[0]
                hist['top'] = pd.Series(h)
                hist['left'] = pd.Series(bins[:-1])
                hist['right'] = pd.Series(bins[1:])
            else:
                hist[c] = pd.Series(h)
                hist[f'{c}_left'] = pd.Series(bins[:-1])
                hist[f'{c}_right'] = pd.Series(bins[1:])

        del hist['empty']
        hists.append(hist)

    # Make the bottom of the bins line up
    bottom = 1e99
    for hist in hists:
        if hist['bottom'][0] < bottom:
            bottom = hist['bottom'][0]

    for hist in hists:
        hist['bottom'] = bottom

    if mc:
        return pd.concat(hists)
    else:
        return hists
