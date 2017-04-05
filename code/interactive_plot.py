"""
Plot FRB populations with a Bokeh server. Run script with:

    $ bokeh serve --show code/interactive_plot.py --args <pop_example.csv>

in which all csv-files with populations can be given after ``--args``, and as
well as the optional arguments of ``-noshow`` and ``-nofrbcat``, to
respectively not show the resulting plot, and to not overplot frbcat
"""

import numpy as np
import os
import pandas as pd
import sys

from bokeh.io import curdoc
from bokeh.layouts import row, column, layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div
from bokeh.models.widgets import Select
from bokeh.palettes import Category10
from bokeh.plotting import figure

from frbcat import get_frbcat
from log import pprint

# Number of dataframes/populations
num_df = 0


def plot_pop(files=[], frbcat=True):
    """
    Function to plot populations in browser using Bokeh

    Args:
        files (list): List of population files to plot (currently only works
                      with csv files - file an issue if you would like more
                      options)
        frbcat (bool): Whether to plot frbcat parameters. Defaults to True
    """

    # Configure colours
    mc = len(files)
    if frbcat:
        mc += 1
    colours = Category10[10][:mc]

    # Dataframes
    dfs = []

    def read(path=None):
        '''
        Mini-function to read in data

        Args:
            path (str): Path to file to read
        '''

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

        df['population'] = name
        df['color'] = colours[num_df]
        dfs.append(df)
        num_df += 1

    # Get dataframe from files
    for f in files:
        read(path=f)

    # Import frbcat
    if frbcat:
        df = get_frbcat()
        num = len(dfs)
        df['color'] = colours[num]
        dfs.append(df)

    # Create Column Data Sources for interacting with the plot
    props = dict(x=[], y=[], color=[], population=[])
    sources = [ColumnDataSource(props) for df in dfs]

    # Create axis options
    axis_map = {
        'Dispersion Measure - Host (pc/cm^3)': 'dm_host',
        'Dispersion Measure - IGM (pc/cm^3)': 'dm_igm',
        'Dispersion Measure - Milky Way (pc/cm^3)': 'dm_mw',
        'Dispersion Measure (pc/cm^3)': 'dm',
        'Distance (Gpc)': 'dist',
        'Fluence (Jy*ms)': 'fluence',
        'Galactic Latitude (degrees)': 'gb',
        'Galactic Longitude (degrees)': 'gl',
        'Galactic X (Gpc)': 'gx',
        'Galactic Y (Gpc)': 'gy',
        'Galactic Z (Gpc)': 'gz',
        'Luminosity - Bolometric': 'lum_bol',
        'Peak Flux Density (W/(m**2*Hz))': 's_peak',
        'Pulse Width - Effective (ms)': 'w_eff',
        'Pulse Width - Intrinsic (ms)': 'w_int',
        'Redshift': 'z',
        'Signal to Noise Ratio': 'snr',
        'Spectral Index': 'si',
        }

    x_axis = Select(title='',
                    options=sorted(axis_map.keys()),
                    value='Galactic Latitude (degrees)')

    # Add histogram option
    axis_map['Fraction'] = 'number'

    y_axis = Select(title='',
                    options=sorted(axis_map.keys()),
                    #value='Number')
                    value='Galactic Longitude (degrees)')

    # What to display while hovering
    props = [("Pop", "@population"), ("X", "@x"), ("Y", "@y")]
    hover = HoverTool(tooltips=props)

    # Set up tools
    tools_to_show = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

    # Set up title
    title = '|'
    for df in dfs:
        name = df['population'].iloc[0]
        num = df.shape[0]
        title += ' {}: {} |'.format(name, num)

    # Create main plot
    p = figure(plot_height=600,
               plot_width=600,
               title=title,
               active_scroll='wheel_zoom',
               toolbar_location='above',
               tools=tools_to_show,
               webgl=True)

    for source in sources:

        # Plot scatter plot for the populations
        p.circle(x='x',
                 y='y',
                 source=source,
                 size=7,
                 alpha=0.7,
                 color='color',
                 legend='population')

    # How to filter data
    def filter_data(df, x_name, y_name, bin_edges):

        # How to plot a histogram
        if y_name == 'number':

            # Initialise different histograms
            xn = []
            yn = []

            # Find x-axis values in that population
            dfx = df[x_name]

            # Check there are values
            dfx = pd.to_numeric(dfx, errors='coerce')
            dfx.dropna(inplace=True)

            print(dfx)
            # Ensure bins overlap if existing
            print(bin_edges)
            if len(bin_edges) > 0:
                be = bin_edges
            else:
                be = 15

            hist, bin_edges = np.histogram(dfx, bins=be)

            # Normalise
            hist = [h/dfx.count() for h in hist]

            # Ugly hack for plotting a histogram with points
            m = 15
            xn = []
            yn = []
            cn = []
            pn = []

            for i, v in enumerate(hist):

                if i < (len(hist) - 1):
                    upper_y = hist[i+1]
                else:
                    upper_y = 0

                dy = upper_y - hist[i]
                dx = bin_edges[i+1] - bin_edges[i]

                xh = [bin_edges[i] + a*dx/m for a in range(m)]
                xh.extend([bin_edges[i+1] for _ in range(m)])
                yh = [v for _ in range(m)]
                yh.extend([v + a*dy/m for a in range(m)])

                xn.extend(xh)
                yn.extend(yh)
                cn.extend(df['color'].iloc[0] for i in range(len(xh)))
                pn.extend(df['population'].iloc[0] for i in range(len(xh)))

            props = {x_name: xn, 'number': yn, 'color': cn, 'population': pn}

            df = pd.DataFrame(props)

        # Only return relevant values
        df = df[[x_name, y_name, 'color', 'population']]
        df = pd.to_numeric(df, errors='coerce')
        df = df.dropna()

        return df, bin_edges

    # Interactive goodness
    def update():
        x_name = axis_map[x_axis.value]
        y_name = axis_map[y_axis.value]

        p.xaxis.axis_label = x_axis.value
        p.yaxis.axis_label = y_axis.value

        bins = []

        for i, source in enumerate(sources):

            # Create an empty data set
            cols = [x_name, y_name, 'color', 'population']
            empty = pd.DataFrame(np.nan, index=[0], columns=cols)

            # Ensure columns are present in each population
            if x_name not in dfs[i]:
                df = empty
            elif y_name not in dfs[i] and y_name != 'number':
                df = empty
            else:
                # Apply filtering
                df, bins = filter_data(dfs[i], x_name, y_name, bins)

            # Don't plot if empty
            #if df.shape[0] > 0:

            # Update data
            source.data = dict(
                x=df[x_name],
                y=df[y_name],
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
    L = layout([[s, p]], sizing_mode=sizing_mode)

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
