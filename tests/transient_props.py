"""Create a 4D graph of radio transient objects."""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from frbpoppy import Frbcat
import frbpoppy.galacticops as go
from convenience import plot_aa_style, rel_path


def import_frbcat():
    cat = Frbcat()
    cat.filter(one_per_frb=True, repeat_bursts=True, repeaters=True)
    df = cat.df

    # Remove Pushichino events
    df = df[~df.telescope.str.startswith('pushchino')]

    db = pd.DataFrame()

    dm_igm = df['dm'] - df['dm_mw']
    db['z'] = dm_igm / 950
    db['dist_co'] = go.Redshift(db['z']).dist_co() * 1e3  # Gpc -> kpc
    db['pseudo_lum'] = (df.fluence / df.w_eff) * db.dist_co**2
    # Observing bandwidth
    db['bw'] = df['bandwidth']*1e-3  # MHz - > GHz
    # Pulse width
    db['w_eff'] = df['w_eff'] * 1e-3  # ms -> s
    # Object type
    db['obj'] = df['obj']

    return db


def plot_data(df):
    plot_aa_style(cols=2)
    plt.rcParams["figure.figsize"] = (5.75373, 5.75373)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for obj, mdf in df.groupby(['obj']):
        ax.scatter(np.log10(mdf.w_eff),
                   np.log10(mdf.bw),
                   np.log10(mdf.pseudo_lum),
                   label=obj)

    def power_of_10(n):
        c = 10 ** round(np.log10(n))
        if type(n) == float:
            c = float(c)
        return np.isclose(c, n) or 10 * c == n

    # Set axis properties
    def log_tick_formatter(val, pos=None):
        new_val = 10**val
        if power_of_10(new_val):
            return "{:.3g}".format(new_val)
        else:
            return ''

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))

    # Set labels
    ax.set_xlabel(r'Pulse Width (s)')
    ax.set_ylabel(r'Bandwidth (GHz)')
    ax.set_zlabel(r'S$_{\text{peak}}$D$^2$ (Jy kpc$^2$)')
    plt.legend()

    # Save figure
    plt.show()
    # plt.savefig(rel_path('./plots/transients.pdf'))


if __name__ == '__main__':
    test_df = pd.DataFrame({'obj': ['a', 'a', 'b'],
                            'pseudo_lum': [1e40, 1e41, 1e40],
                            'w_eff': [1, 100, 0.1],
                            'bw': [100, 1, 100]})
    df = import_frbcat()
    plot_data(df)
