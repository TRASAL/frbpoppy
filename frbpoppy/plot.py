"""
Plot FRB populations with a Bokeh server. Run script with:

    $ bokeh serve --show code/plot.py --args <pop_example.csv>

in which all csv-files with populations can be given after ``--args``, and as
well as the optional arguments of ``-noshow`` and ``-nofrbcat``, to
respectively not show the resulting plot, and to not overplot frbcat
"""

import numpy as np
import os
import pandas as pd
import sys
sys.path.append("..")

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div, Panel, Tabs
from bokeh.models.widgets import Select
from bokeh.palettes import Category10, viridis
from bokeh.plotting import figure

from frbpoppy.log import pprint
from frbpoppy.frbcat import Frbcat
from frbpoppy.paths import paths

# Number of dataframes/populations
num_df = 0


def histogram(dfs):
    """
    Quick function to 'histogram' each column of each dataframes.

    Args:
        dfs (list): List of dataframes
    Returns:
        hists (list): List of histogramed dataframes
    """
    # Determine all column names
    columns = [df.columns.values for df in dfs]
    cols = min(columns, key=len)

    # Determine bin limits
    limits = {}

    for c in cols:

        low = 1e99
        high = -1e99

        for df in dfs:
            if c not in df:
                continue

            dlow = df[c].min()
            dhigh = df[c].max()

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

        hist = pd.DataFrame(np.nan, index=np.arange(49), columns=['empty'])
        hist['color'] = df['color'][0]
        hist['population'] = df['population'][0]
        hist['bottom'] = 10**(round(np.log10(1/len(df))) - 1)

        for c in cols:

            low = limits[c][0]
            high = limits[c][1]

            if low == 1e99 or high == -1e99:
                continue

            if c not in df:
                continue

            if df[c].nunique() == 1 and df[c][0] == 'None':
                continue

            bins = np.linspace(low, high, 50)

            if high - low > 1000:
                if low == 0:
                    low = 1e-3
                bins = np.geomspace(low, high, num=50)

            col = df[c].apply(pd.to_numeric, errors='coerce')
            col = col.dropna()
            h, _ = np.histogram(col, bins=bins)

            # Normalise
            h = [e/h.sum() for e in h]

            hist[c] = pd.Series(h)
            hist[c+'_left'] = pd.Series(bins[:-1])
            hist[c+'_right'] = pd.Series(bins[1:])

        del hist['empty']
        hists.append(hist)

    return hists


def plot_pop(files=[], frbcat=True, hist_axis='log'):
    """
    Plot populations in browser using Bokeh.

    Args:
        files (list): List of population files to plot (currently only works
                      with csv files - file an issue if you would like more
                      options)
        frbcat (bool): Whether to plot frbcat parameters. Defaults to True
        hist_log (bool): Whether to plot using a linear or log axis

    """
    # Configure colours
    mc = len(files)
    if frbcat:
        mc += 1
    if mc > 10:
        colours = viridis(mc)
    else:
        colours = Category10[10][:mc]

    # Dataframes
    dfs = []

    def read(path=None):
        """
        Mini-function to read in data

        Args:
            path (str): Path to file to read
        """
        global num_df

        if os.path.isfile(path):
            try:
                df = pd.read_csv(path)
            except ValueError:
                return
            name = f.split('_')[-1].split('.')[0]
        else:
            m = 'Skipping population {} - contains no sources'.format(f)
            pprint(m)
            return

        # Reduce population size if it's too large
        if df.shape[0] > 30000:
            df = df.iloc[::10]

        df['population'] = name
        df['color'] = colours[num_df]
        # Sidestepping issue in Bokeh
        df['lum_bol'] = df['lum_bol'] / 1e30
        dfs.append(df)
        num_df += 1

    # Get dataframe from files
    for f in files:
        read(path=f)

    # Import frbcat
    if frbcat:
        df = Frbcat().df
        num = len(dfs)
        df['color'] = colours[num]
        dfs.append(df)

    # Create axis options
    axis_map = {
        'Comoving Distance (Gpc)': 'dist_co',
        'Declination (°)': 'dec',
        'Dispersion Measure - Host (pc/cm^3)': 'dm_host',
        'Dispersion Measure - IGM (pc/cm^3)': 'dm_igm',
        'Dispersion Measure - Milky Way (pc/cm^3)': 'dm_mw',
        'Dispersion Measure (pc/cm^3)': 'dm',
        'Fluence (Jy*ms)': 'fluence',
        'Galactic Latitude (degrees)': 'gb',
        'Galactic Longitude (degrees)': 'gl',
        'Galactic X (Gpc)': 'gx',
        'Galactic Y (Gpc)': 'gy',
        'Galactic Z (Gpc)': 'gz',
        'Luminosity - Bolometric (10^30 ergs/s)': 'lum_bol',
        'Peak Flux Density (Jy)': 's_peak',
        'Pulse Width - Effective (ms)': 'w_eff',
        'Pulse Width - Intrinsic (ms)': 'w_int',
        'Redshift': 'z',
        'Right Ascension (°)': 'ra',
        'Signal to Noise Ratio': 'snr',
        }

    x_axis = Select(title='',
                    options=sorted(axis_map.keys()),
                    value='Galactic Longitude (degrees)')

    y_axis = Select(title='',
                    options=sorted(axis_map.keys()),
                    value='Galactic Latitude (degrees)')

    # Set up tools
    props = [("pop", "@population"), ("x", "@x"), ("y", "@y")]
    hover = HoverTool(tooltips=props)
    sc_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

    # Create scatter plot
    sp = figure(plot_height=700,
                plot_width=700,
                active_scroll='wheel_zoom',
                toolbar_location='right',
                tools=sc_tools)

    # Stop labels falling off
    sp.min_border_left = 80

    # Create Column Data Sources for interacting with the plot
    props = dict(x=[], y=[], color=[], population=[])
    sp_sources = [ColumnDataSource(props) for df in dfs]

    # Plot scatter plot for the populations
    for source in sp_sources:
        sp.circle(x='x',
                  y='y',
                  source=source,
                  size=7,
                  alpha=0.6,
                  color='color',
                  legend='population')

    # Set up tools
    props = [("pop", "@population"), ("frac", "@top")]
    hover = HoverTool(tooltips=props)
    hp_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

    # Create histogram plot
    hp = figure(plot_height=700,
                plot_width=700,
                active_scroll='wheel_zoom',
                toolbar_location='right',
                tools=hp_tools,
                x_axis_type=hist_axis,
                y_axis_type="log")

    # Create Column Data Sources for interacting with the plot
    hists = histogram(dfs)
    props = dict(top=[], left=[], right=[], bottom=[], color=[], population=[])
    hp_sources = [ColumnDataSource(props) for hist in hists]

    # Plot histogram values
    for source in hp_sources:
        hp.quad(bottom='bottom',
                left='left',
                right='right',
                top='top',
                color='color',
                line_color='color',
                legend='population',
                alpha=0.4,
                source=source)

    # Interactive goodness
    def update():
        x_name = axis_map[x_axis.value]
        y_name = axis_map[y_axis.value]

        sp.xaxis.axis_label = x_axis.value
        sp.yaxis.axis_label = y_axis.value
        hp.xaxis.axis_label = x_axis.value
        hp.yaxis.axis_label = 'Fraction'

        for i, source in enumerate(sp_sources):

            # Create an empty data set
            cols = [x_name, y_name, 'color', 'population']
            empty = pd.DataFrame(np.nan, index=[0], columns=cols)

            # Ensure columns are present in each population
            if x_name not in dfs[i] or y_name not in dfs[i]:
                df = empty
            else:
                df = dfs[i][cols]
                df = df.replace('None', np.nan)
                df[x_name].apply(pd.to_numeric, errors='coerce')
                df[y_name].apply(pd.to_numeric, errors='coerce')
                df = df.dropna()

            # Update data
            source.data = dict(
                x=df[x_name],
                y=df[y_name],
                color=df['color'],
                population=df['population']
            )

        for i, source in enumerate(hp_sources):
            # Create an empty data set
            cols = [x_name,
                    x_name + '_left',
                    x_name + '_right',
                    'bottom',
                    'color',
                    'population']
            empty = pd.DataFrame(np.nan, index=[0], columns=cols)

            # Ensure columns are present in each population
            if x_name not in hists[i]:
                df = empty
            else:
                df = hists[i]
                df = df[cols]
                for e in cols[:4]:
                    df[e].apply(pd.to_numeric)
                df = df.dropna()

            # Update data
            source.data = dict(
                top=df[x_name],
                left=df[x_name + '_left'],
                right=df[x_name + '_right'],
                bottom=df['bottom'],
                color=df['color'],
                population=df['population']
            )

    # What to interact with
    controls = [x_axis, y_axis]

    for control in controls:
        control.on_change('value', lambda attr, old, new: update())

    # Layout options
    sizing_mode = 'fixed'

    cwd = os.path.dirname(__file__)

    def path(p):
        d = os.path.join(cwd, 'plot_config/{}.html'.format(p))
        return open(d).read()

    text_top = Div(text=path('text_top'))
    text_bottom = Div(text=path('text_bottom'))

    sidebar = [text_top, x_axis, y_axis, text_bottom]
    s = widgetbox(sidebar, width=350)
    tab1 = Panel(child=sp, title='scatter')
    tab2 = Panel(child=hp, title='histogram')
    tabs = Tabs(tabs=[tab1, tab2])
    L = layout([[s, tabs]], sizing_mode=sizing_mode)

    sp.legend.click_policy = 'hide'
    hp.legend.click_policy = 'hide'
    update()  # initial load of the data
    curdoc().add_root(L)
    curdoc().title = 'frbpoppy'


# Parse system arguments
# I know ArgumentParser is nicer, but bokeh only works with argv

args = sys.argv

frbcat = True
if '-nofrbcat' in args:
    frbcat = False

files = []
for a in args:
    if a.endswith('.csv'):
        files.append(a)

# Check whether populations have been given as input
if len(files) == 0:
    pprint('Nothing to plot: plot arguments are empty')
else:
    plot_pop(files=files, frbcat=frbcat)
