"""Create a 4D graph of radio transient objects."""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import urllib.request
from matplotlib.markers import MarkerStyle

from frbpoppy import Frbcat, pprint
import frbpoppy.galacticops as go

from tests.convenience import rel_path


def import_frbcat():
    """Import frbcat."""
    cat = Frbcat(frbpoppy=False)
    cat.frbpoppify()
    cat.filter(one_offs=True, repeaters=True, repeat_bursts=False)
    df = cat.df

    # Remove Pushichino events
    df = df[~df.telescope.str.startswith('pushchino')]

    db = pd.DataFrame()

    dm_igm = df['dm'] - df['dm_mw']
    db['z'] = dm_igm / 950
    db['dist_co'] = go.Redshift(db['z']).dist_co() * 1e6  # Gpc -> kpc
    db['pseudo_lum'] = (df.fluence / df.w_eff) * db.dist_co**2
    # Observing bandwidth
    db['bw'] = df['bandwidth']  # MHz
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
                bw = (max(freqs) - min(freqs))
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
    db['bw'] = df.apply(calc_bandwidth, axis=1)  # MHz

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
        epn_db['bw'] += 500  # Set a default of 500 MHz
        merged_db = pd.merge(db, epn_db, on='name')
        merged_db['bw'] = merged_db[['bw_x', 'bw_y']].apply(np.max, axis=1)
        merged_db.drop(['bw_x', 'bw_y'], inplace=True, axis=1)
        return merged_db
    else:
        return db


def import_magnetars():
    """Import magnetar info."""
    m_db = {'bw': [41*1e3, 290*1e3],
            'pseudo_lum': [1.1*3.3**2, 0.6*8.3**2],
            'w_eff': [0.6, 1.8]}
    m_db['type'] = ['magnetar' for i in range(len(m_db['bw']))]
    return pd.DataFrame.from_dict(m_db)


def import_crab_giants():
    """Import crab giant pulse info."""
    crab_db = {'bw': [8.1*1e3],
               'pseudo_lum': [4e6*2**2],
               'w_eff': [2e-3]}
    crab_db['type'] = ['giant pulse' for i in range(len(crab_db['bw']))]
    return pd.DataFrame.from_dict(crab_db)


def plot_data(df):
    """Plot transient properities in a DataFrame."""
    # plot_aa_style(cols=2)
    plt.rcParams["figure.figsize"] = (5.75373, 5.75373*3)
    plt.rcParams["font.family"] = "Times"
    plt.rcParams['pdf.fonttype'] = 42

    fig = plt.figure()
    ax = fig.add_subplot(111, projection=Axes3D.name)

    # Clean data
    df = df[df.bw < 1e10]

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    markers = list(MarkerStyle.markers.keys())
    markers = ['*', '^', 's', '.', 's']

    c = 0
    for obj, mdf in df.groupby(['type']):
        print(f'Plotting {obj}s')
        color = colors[c]
        marker = markers[c]

        c += 1
        alpha = 1
        linewidth = 1
        markersize = 30
        if len(mdf.w_eff) > 5:
            alpha = 0.3
            linewidth = 0.7
            markersize = 10
        if len(mdf.w_eff) > 100:
            alpha = 0.1
            linewidth = 0.2
            markersize = 5

        ax.scatter(np.log10(mdf.w_eff*1e3),
                   np.log10(mdf.bw),
                   np.log10(mdf.pseudo_lum),
                   label=obj,
                   color=color,
                   s=markersize,
                   marker=marker)

        for i in range(len(mdf)):
            ax.plot([np.log10(mdf.w_eff.iloc[i]*1e3) for n in range(2)],
                    [np.log10(mdf.bw.iloc[i]) for n in range(2)],
                    [np.log10(mdf.pseudo_lum.iloc[i]), -4],
                    color=color,
                    alpha=alpha,
                    linewidth=linewidth)

    def power_of_10(n):
        c = 10 ** round(np.log10(n))
        if type(n) == float:
            c = float(c)
        return np.isclose(c, n) or 10 * c == n

    # Set axis properties
    def log_tick_formatter(val, pos=None):
        new_val = 10**val
        return r"$10^{%02d}$" % (val)
        if power_of_10(new_val):
            return "{:.3g}".format(new_val)
        else:
            return ''

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))

    # Set labels
    ax.set_xlabel(r'Pulse Width (ms)')
    ax.set_ylabel(r'Bandwidth (MHz)')
    ax.set_zlabel(r'S$_{peak}$D$^2$ (Jy kpc$^2$)')
    plt.legend()

    # Equal figure proportions
    ax.set_xlim(-1, 4)
    ax.set_ylim(1, 6)
    ax.set_zlim(-4, 14)

    # Get rid of colored axes planes
    # First remove fill
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    # Now set color to white (or whatever is "invisible")
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    # Save figure
    ax.view_init(azim=-52, elev=10)
    plt.tight_layout()

    plt.savefig(rel_path('./plots/transients.pdf'), transparent=True)
    # plt.show()


if __name__ == '__main__':

    RELOAD = False
    csv_path = rel_path('./plots/transients.csv')

    if RELOAD:
        test = pd.DataFrame({'obj': ['a', 'a', 'b'],
                             'pseudo_lum': [1e40, 1e41, 1e40],
                             'w_eff': [1, 100, 0.1],
                             'bw': [800, 1400, 1000]})

        pulsars = import_pulsars(use_epn=True)
        frbs = import_frbcat()
        magnetars = import_magnetars()
        giants = import_crab_giants()

        tot_df = pd.concat([frbs, pulsars, magnetars, giants], sort=False)

        tot_df.to_csv(csv_path)
    else:
        tot_df = pd.read_csv(csv_path)

    plot_data(tot_df)
