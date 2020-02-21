"""TODO."""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

from frbpoppy import Frbcat
import frbpoppy.galacticops as go
from convenience import plot_aa_style, rel_path, hist

SCALE = 'lin'

def import_frbcat():
    cat = Frbcat()
    cat.filter(one_per_frb=True, one_offs=False, repeaters=True,
               repeat_bursts=True)
    df = cat.df

    # Remove Pushichino events
    df = df[~df.telescope.str.startswith('pushchino')]

    # Keep CHIME
    df = df[df.telescope.str.startswith('chime')]

    db = pd.DataFrame()
    db['utc'] = df.utc
    dm_igm = df['dm'] - df['dm_mw']
    db['z'] = dm_igm / 950
    db['dist_co'] = go.Redshift(db['z']).dist_co() * 1e3  # Gpc -> kpc
    db['pseudo_lum'] = (df.fluence / df.w_eff) * db.dist_co**2
    # Observing bandwidth
    db['bw'] = df['bandwidth']*1e-3  # MHz - > GHz
    # Pulse width
    db['w_eff'] = df['w_eff'] * 1e-3  # ms -> s
    # Object type
    db['type'] = df['type']

    # Add exposure time information
    exposure = {'FRB180725.J0613+67': np.nan,
                'FRB180814.J0422+73': 0.09,
                'FRB180908.J1232+74': 4/(.5*(53+36)),
                'FRB180916.J0158+65': 10/23,
                'FRB181017.J1705+68': 2/20,
                'FRB181017.J18+81': 3/(.5*(55+159)),
                'FRB181030.J1054+73': 2/(.5*(27+19)),
                'FRB181119.J12+65': 3/19,
                'FRB181128.J0456+63': 2/16,
                'FRB190110.J1353+48': 3/23,
                'FRB190116.J1249+27': 2/8,
                'FRB190117.J2207+17': 5/19,
                'FRB190208.J1855+46': 2/20,
                'FRB190209.J0937+77': 2/(.5*(34+28)),
                'FRB190212.J02+20': 2/17,
                'FRB190222.J2052+69': 2/20,
                'FRB190417.J1939+59': 3/29,
                'FRB190604.J1435+53': 2/30,
                'FRB190907.J08+46': 3/23}

    db['rate'] = np.nan
    db['frb_name'] = df['frb_name']
    for name, value in exposure.items():
        db.loc[db.frb_name == name, 'rate'] = value

    return db


def plot_w_eff(df):
    plot_aa_style(cols=1)

    fig = plt.figure()

    for obj, mdf in df.groupby(['type']):
        plt.plot(*hist(mdf.w_eff*1e3, bin_type=SCALE, norm=False, bins=5),
                 label=obj)

    plt.xlabel(r'Pulse Width (ms)')
    if SCALE == 'log':
      plt.xscale('log')
    plt.legend()
    plt.savefig(rel_path('./plots/w_eff_frbcat.pdf'))


def plot_w_eff_rate(df):
    plot_aa_style(cols=1)

    fig = plt.figure()

    df = df.sort_values('utc')
    df['n_burst'] = np.where(df.duplicated('frb_name'), 'other', 'first')
    import IPython; IPython.embed()
    for obj, mdf in df.groupby(['n_burst']):
        plt.scatter(mdf.w_eff*1e3, mdf.rate, label=obj, marker='x')

    plt.xlabel(r'Pulse Width (ms)')
    plt.ylabel(r'Rate (/hour)')

    if SCALE == 'log':
      plt.xscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(rel_path('./plots/w_eff_rate_frbcat.pdf'))


if __name__ == '__main__':
    df = import_frbcat()
    plot_w_eff(df)
    plot_w_eff_rate(df)
