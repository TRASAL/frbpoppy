"""Create a 4D graph of radio transient objects."""
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import urllib.request

from frbpoppy import Frbcat, pprint
import frbpoppy.galacticops as go

from convenience import plot_aa_style, rel_path, set_axes_equal


def import_frbcat():
    """Import frbcat."""
    cat = Frbcat()
    cat.filter(one_offs=True, repeaters=True, repeat_bursts=True)
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
    db['type'] = df['type']

    return db


def bw_from_epn():
    """Scrape bandwidth from EPN website."""
    pprint('Getting data from EPN')
    bw_db = {'name': [], 'bw': []}
    url = 'http://www.epta.eu.org/epndb/list.php'

    def clean(l):
        return float(l.split('>')[1].split(' ')[0])

    with urllib.request.urlopen(url) as resp:
        data = resp.read().decode().replace('&nbsp', ' ')
        for line in data.split('\n'):
            if line.startswith('<li>'):
                name = line.split()[0].split('>')[-1]
                freqs = [clean(l) for l in line.split('<a ')[1:]]
                bw = (max(freqs) - min(freqs)) / 1e3
                bw_db['name'].append(name)
                bw_db['bw'].append(bw)

    pprint('Finished getting data from EPN')
    bw_db = pd.DataFrame.from_dict(bw_db)
    return bw_db


def import_pulsars(use_epn=True):
    """Import pulsar info from ATNF."""
    from psrqpy import QueryATNF

    query = QueryATNF()
    df = query.pandas

    db = pd.DataFrame()

    def calc_bandwidth(row):
        bw = 500

        s400 = ~np.isnan(row.S400)
        s1400 = ~np.isnan(row.S1400)
        s2000 = ~np.isnan(row.S2000)

        if s400 & s1400:
            bw += (1400 - 400)
        if (s400 & s2000) or (s400 & s1400 & s2000):
            bw += (2000 - 400)
        if s1400 & s2000:
            bw += (2000 - 1400)
        if sum((s400, s1400, s2000)) == 0:
            bw = np.nan

        return bw

    # Observing bandwidth
    db['bw'] = df.apply(calc_bandwidth, axis=1)*1e-3  # MHz - > GHz

    def calc_flux(row):
        return max([row.S400, row.S1400, row.S2000])

    # Pseudo luminosity in Jy*kpc**2
    db['pseudo_lum'] = df.apply(calc_flux, axis=1) * 1e-3 * df.DIST_DM**2
    # Pulse width
    db['w_eff'] = df['W10'] * 1e-3  # ms -> s
    # Object type
    db['type'] = 'pulsar'

    # Find bw from EPN if available
    if use_epn:
        db['name'] = df.NAME
        epn_db = bw_from_epn()
        epn_db['bw'] += 0.5  # Set a default of 500 MHz
        merged_db = pd.merge(db, epn_db, on='name')
        merged_db['bw'] = merged_db[['bw_x', 'bw_y']].apply(np.max, axis=1)
        merged_db.drop(['bw_x', 'bw_y'], inplace=True, axis=1)
        return merged_db
    else:
        return db


def import_magnetars():
    """Import magnetar info."""
    m_db = {'bw': [41, 290],
            'pseudo_lum': [1.1*3.3**2, 0.6*8.3**2],
            'w_eff': [0.6, 1.8]}
    m_db['type'] = ['magnetar' for i in range(len(m_db['bw']))]
    return pd.DataFrame.from_dict(m_db)


def import_crab_giants():
    """Import crab giant pulse info."""
    crab_db = {'bw': [8.1],
               'pseudo_lum': [4e6*2**2],
               'w_eff': [2e-3]}
    crab_db['type'] = ['giant pulse' for i in range(len(crab_db['bw']))]
    return pd.DataFrame.from_dict(crab_db)


def plot_data(df):
    """Plot transient properities in a DataFrame."""
    plot_aa_style(cols=2)
    plt.rcParams["figure.figsize"] = (5.75373, 5.75373)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for obj, mdf in df.groupby(['type']):
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

    # Equal figure proportions
    # set_axes_equal(ax)

    # Save figure
    plt.tight_layout()
    plt.savefig(rel_path('./plots/transients.pdf'))
    plt.show()


if __name__ == '__main__':
    test = pd.DataFrame({'obj': ['a', 'a', 'b'],
                         'pseudo_lum': [1e40, 1e41, 1e40],
                         'w_eff': [1, 100, 0.1],
                         'bw': [100, 1, 100]})

    pulsars = import_pulsars(use_epn=True)
    frbs = import_frbcat()
    magnetars = import_magnetars()
    giants = import_crab_giants()

    tot_df = pd.concat([frbs, pulsars, magnetars, giants], sort=False)

    plot_data(tot_df)
