from collections import OrderedDict
import math
import numpy as np
import os
import pandas as pd

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.models.widgets import RadioButtonGroup, Select
from bokeh.palettes import Category10
from bokeh.plotting import figure

from log import pprint

# Number of dataframes/populations
num_df = 0

def plot_pop(pops=[], files=[], show=True):
    """
    Function to plot populations in browser using Bokeh

    Args:
        pops (list): List of population objects to plot
        files (list): List of population files to plot (currently only works
                      with csv files - file an issue if you would like more
                      options)
        show (boolean): Show plot or not. Defaults to True
    """

    lp = len(pops)
    lf = len(files)

    # Configure colours
    colours = OrderedDict()

    # Dataframes
    dfs = []

    def read(path=None,pop=None):
        '''
        Mini-function to read in data

        Args:
            path (str): Path to file to read
            pop (str): Population input
        '''
        global num_df

        if path:
            if os.path.isfile(path):
                df = pd.read_csv(path)
                name = f.split('_')[-1].split('.')[0]
                df[name] = True
            else:
                pprint('Skipping population {} - contains no sources'.format(f))
                return

        if pop:
            v = pop.values()
            if v:
                df = pd.read_csv(StringIO(v))
                name = pop.name
                df[name] = True
            else:
                name = pop.name
                pprint('Skipping population {} - contains no sources'.format(name))
                return

        dfs.append(df)
        colours[name] = Category10[10][num_df]
        num_df += 1

    # Check whether populations have been given as input
    if lp == 0 and lf == 0:
        print('Nothing to plot: plot arguments are empty')
        return

    # Get dataframe from populations
    elif lp > 0:
        for p in pops:
            read(pop=p)

    # Get dataframe from files
    elif lf > 0:
        for f in files:
            read(path=f)

    # Combine populations to single dataframe, using updating the initial
    # population if observation values have been found in the second one
    df = dfs[0]

    if num_df > 1:
        for d in dfs:
            df = pd.merge(df,
                          d,
                          on=['gx', 'gy', 'gz'],
                          how='outer',
                          suffixes=('_old',''))

            for c in df:
                if c + '_old' in df:
                    df[c] = df[c].fillna(df[c + '_old'])
                    del df[c + '_old']

    # Add colours to each source
    for n in colours:
        df.loc[(df[n].notnull()), 'color'] = colours[n]
        df.loc[(df[n].notnull()), 'population'] = n

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
        'Peak Flux Density (W/m**2)': 's_peak',
        'Pulse Width - Effective (ms)': 'w_eff',
        'Pulse Width - Intrinsic (ms)': 'w_int',
        'Redshift': 'z',
        'Signal to Noise Ratio': 'snr',
        'Spectral Index': 'si',
        }

    x_axis = Select(title='X-Axis',
                    options=sorted(axis_map.keys()),
                    value='Galactic X (Gpc)')

    # Add histogram option
    axis_map['Number'] = 'number'

    y_axis = Select(title='Y-Axis',
                    options=sorted(axis_map.keys()),
                    value='Number')
                    #value='Galactic Y (Gpc)')

    # Create Column Data Source for interacting with the plot
    source = ColumnDataSource(data=dict(x=[], y=[], color=[]))

    # What to display while hovering
    hover = HoverTool(tooltips=[("Pop", "$population"), ("X", "@x"), ("Y", "@y")])
    # Set up tools
    tools_to_show = ['box_zoom', 'pan', 'save',hover,'reset','wheel_zoom']

    # Create main plot
    p = figure(plot_height=600,
               plot_width=600,
               title='',
               active_scroll='wheel_zoom',
               toolbar_location='above',
               tools=tools_to_show)

    # Plot scatter plot
    p.circle(x='x',
             y='y',
             source=source,
             size=7,
             color='color',
             alpha=0.5)

    # Interactive goodness
    def update():
        x_name = axis_map[x_axis.value]
        y_name = axis_map[y_axis.value]

        p.xaxis.axis_label = x_axis.value
        p.yaxis.axis_label = y_axis.value

        num_pop = len(df)
        num_subpop =len(df[(df.pmsurv.notnull())])

        # Set up title
        title = '|'
        for n in colours:
            num_pop = df[(df[n].notnull())][x_name].shape[0]
            title += ' {}: {} |'.format(n.capitalize(), num_pop)
        p.title.text = title

        # How to plot a histogram
        if y_name == 'number':

            # Initialise different histograms
            xn = []
            yn = []
            cn = []
            bin_edges = []

            # For each population
            for n in colours:

                # Find x-axis values in that population
                dfx = df[(df[n].notnull())][x_name]

                # Check there are values
                dfx = pd.to_numeric(dfx, errors='coerce')
                dfx.dropna(inplace=True)

                # Ensure bins overlap
                if len(bin_edges) > 0:
                    be = bin_edges
                else:
                    be = 15
                hist, bin_edges = np.histogram(dfx, bins=be)

                xh = [val for val in bin_edges for _ in (0, 1)]
                yh = [val for val in hist for _ in (0, 1)]
                xn.extend(xh)
                yn.append(0)
                yn.extend(yh)
                yn.append(0)
                cn.extend([colours[n] for i in range(len(xh))])

            db = pd.DataFrame({x_name: xn, 'number': yn, 'color': cn})

            #p.line(x='x', y='y', source=source)
        else:
            db = df

        source.data = dict(
            x=db[x_name],
            y=db[y_name],
            color=db['color']
        )

    # What to interact with
    controls = [x_axis, y_axis]

    for control in controls:
        control.on_change('value', lambda attr, old, new: update())

    # Layout options
    sizing_mode = 'scale_width'

    inputs = widgetbox(*controls, width=100)
    l = layout([
        [inputs, p],
    ], sizing_mode=sizing_mode)

    if show:
        update()  # initial load of the data
        curdoc().add_root(l)
        curdoc().title = 'frbpoppy'

plot_pop(pops=[], files=['/media/DATA/phd/frbpoppy/data/results/population_initial.csv',
                             '/media/DATA/phd/frbpoppy/data/results/population_pmsurv.csv'])
