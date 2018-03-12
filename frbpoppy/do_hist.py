"""Bin up Pandas DataFrames into histograms ready for Bokeh plotting."""
import numpy as np
import pandas as pd


def histogram(dfs):
    """
    Quick function to 'histogram' each column of each dataframe

    Args:
        dfs (list): List of dataframes
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

        for c in cols:

            hist = pd.DataFrame(np.nan,
                                index=np.arange(14),
                                columns=['empty'])

            low = limits[c][0]
            high = limits[c][1]

            if low == 1e99 or high == -1e99:
                continue

            if c not in df:
                continue

            if df[c].nunique() == 1 and df[c].iloc[0] == 'None':
                continue

            if high - low > 1500:
                if low == 0:
                    low = 1e-3
                bins = np.geomspace(low, high, num=15)
            else:
                bins = np.linspace(low, high, 15)

            col = df[c].convert_objects(convert_numeric=True)
            col = col.dropna()
            h, _ = np.histogram(col, bins=bins)

            # Normalise
            h = [e/h.sum() for e in h]

            hist['id'] = df['id'].iloc[0]
            hist['in_par'] = df['in_par'].iloc[0]
            hist['out_par'] = c
            hist['top'] = pd.Series(h)
            hist['bottom'] = 0.
            hist['left'] = pd.Series(bins[:-1])
            hist['right'] = pd.Series(bins[1:])
            hist['survey'] = df['survey'].iloc[0]
            hist['frbs'] = len(col)  # number of frbs in population
            hist['frbcat'] = df['frbcat'].iloc[0]

            del hist['empty']

            hists.append(hist)

    return pd.concat(hists)
