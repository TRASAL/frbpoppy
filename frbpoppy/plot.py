"""Allow FRB populations to be explored interactively.

Can be run with:

    $ bokeh serve --show code/plot.py --args <pop_example.csv>

in which all csv-files with populations can be given after ``--args``, and as
well as the optional arguments of ``-noshow`` and ``-nofrbcat``, to
respectively not show the resulting plot, and to not overplot frbcat
"""

import numpy as np
import os
import pandas as pd
import sys

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div, Panel, Tabs
from bokeh.models.widgets import Select
from bokeh.palettes import Category10, viridis
from bokeh.plotting import figure

from frbpoppy.frbcatpoppy import Frbcat
from frbpoppy.do_hist import histogram
from frbpoppy.log import pprint
from frbpoppy import unpickle


class Tab():
    """Gather all elements needed for a plot."""

    def __init__(self):
        """Initializing."""
        self.fig = None
        self.sources = []
        self.name = ''


class Plot():
    """Gather plotting options."""

    def __init__(self, files=[], frbcat=True):
        """Initializing."""
        # From arguments
        self.files = files
        self.frbcat = frbcat

        # Predefined
        self.height = 700  # Plot height
        self.width = 700  # Plot width

        # Initializing arguments set later in the code
        self.x_axis = None
        self.y_axis = None
        self.n_df = 0
        self.dfs = []
        self.tabs = []

        # Set parameters
        self.params = {'Comoving Distance (Gpc)': 'dist_co',
                       'Declination (°)': 'dec',
                       'Dispersion Measure - Host (pc/cm^3)': 'dm_host',
                       'Dispersion Measure - IGM (pc/cm^3)': 'dm_igm',
                       'Dispersion Measure - Milky Way (pc/cm^3)': 'dm_mw',
                       'Dispersion Measure (pc/cm^3)': 'dm',
                       'Fluence (Jy*ms)': 'fluence',
                       'Galactic Latitude (°)': 'gb',
                       'Galactic Longitude (°)': 'gl',
                       'Galactic X (Gpc)': 'gx',
                       'Galactic Y (Gpc)': 'gy',
                       'Galactic Z (Gpc)': 'gz',
                       'Luminosity - Bolometric (10^30 ergs/s)': 'lum_bol',
                       'Peak Flux Density (Jy)': 's_peak',
                       'Pulse Width - Effective (ms)': 'w_eff',
                       'Pulse Width - Intrinsic (ms)': 'w_int',
                       'Redshift': 'z',
                       'Right Ascension (°)': 'ra',
                       'Spectral Index': 'si',
                       'Signal to Noise Ratio': 'snr',
                       'Time (days)': 'time'}

        # Running plotting
        self.set_colours()
        self.set_widgets()
        self.get_data()
        self.make_scatter()
        self.make_histogram(kind='lin')
        self.make_histogram(kind='log')
        self.make_histogram(kind='cum')
        self.set_layout()

    def set_colours(self):
        """Determine which colours need to be used."""
        # Ensure number of overplots is known
        n = len(self.files)
        if self.frbcat:
            n += 1
        if n > 10:
            self.colours = viridis(n)
        else:
            self.colours = Category10[10][:n]

    def get_data(self):
        """Read in populations."""
        # Read in files
        for f in self.files:

            # Check whether file exists
            if os.path.isfile(f):
                try:
                    df = unpickle(f).frbs.to_df()
                except ValueError:
                    pprint(f'Unpacking {f} seemed to have failed.')
                    continue
                if '.' in f:
                    name = f.split('/')[-1].split('.')[0]
                    if '_for_plotting' in name:
                        name = name.split('_for_plotting')[0]
                else:
                    name = f

            # If things haven't worked
            if df is None:
                m = 'Skipping population {} - contains no sources'.format(f)
                pprint(m)
                continue

            # Downsample population size if it's too large
            if df.shape[0] > 10000:
                pprint(f'Downsampling population {f} (else too big to plot)')
                df = df.sample(n=10000)

            df['population'] = name
            df['color'] = self.colours[self.n_df]
            df['lum_bol'] = df['lum_bol'] / 1e30  # Sidestepping Bokeh issue

            if df.empty:
                m = 'Skipping population {} - contains no sources'.format(f)
                pprint(m)
                continue
            else:
                self.dfs.append(df)
                self.n_df += 1

        # Add on frbcat
        if self.frbcat:
            df = Frbcat(frbpoppy=True).df
            # Filter by survey if wished
            if isinstance(self.frbcat, str):
                if df['survey'].str.match(self.frbcat).any():
                    df = df[df.survey == self.frbcat]
                elif df['telescope'].str.match(self.frbcat).any():
                    df = df[df.telescope == self.frbcat]
                else:
                    m = 'Your chosen input for frbcat is not found.'
                    raise ValueError(m)

                df['population'] = f'frbcat {self.frbcat}'

            df['color'] = self.colours[len(self.dfs)]
            self.dfs.append(df)

    def set_widgets(self):
        """Set up widget details."""
        self.x_axis = Select(title='',
                             options=sorted(self.params.keys()),
                             value='Right Ascension (°)')

        self.y_axis = Select(title='',
                             options=sorted(self.params.keys()),
                             value='Declination (°)')

    def make_scatter(self):
        """Set up a scatter plot."""
        # Initializing plot
        tab = Tab()
        tab.name = 'Scatter'

        # Set up interactive tools
        props = [("pop", "@population"), ("x", "@x"), ("y", "@y")]
        hover = HoverTool(tooltips=props)
        tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create scatter plot
        tab.fig = figure(plot_height=self.height,
                         plot_width=self.width,
                         active_scroll='wheel_zoom',
                         toolbar_location='right',
                         tools=tools)

        # Stop labels falling off
        tab.fig.min_border_left = 80

        # Create Column Data Sources for interacting with the plot
        props = dict(x=[], y=[], color=[], population=[])
        tab.sources = [ColumnDataSource(props) for df in self.dfs]

        # Plot scatter plot of FRB populations
        for source in tab.sources:
            tab.fig.circle(x='x',
                           y='y',
                           source=source,
                           size=7,
                           alpha=0.6,
                           color='color',
                           legend='population')

        self.tabs.append(tab)

    def make_histogram(self, kind='lin'):
        """Set up a histogram plot."""
        # Initializing plot
        tab = Tab()

        if kind == 'lin':
            tab.name = 'Hist (Lin)'
            axis_type = 'linear'
            log = False
            cum = False
        elif kind == 'log':
            tab.name = 'Hist (Log)'
            axis_type = 'log'
            log = True
            cum = False
        elif kind == 'cum':
            tab.name = 'Hist (Cum)'
            axis_type = 'log'
            log = True
            cum = True

        # Set up interactive tools
        tools = ['box_zoom', 'pan', 'save', 'reset', 'wheel_zoom']

        # Create histogram plot
        tab.fig = figure(plot_height=self.height,
                         plot_width=self.width,
                         active_scroll='wheel_zoom',
                         toolbar_location='right',
                         tools=tools,
                         x_axis_type=axis_type,
                         y_axis_type="log")

        # Create Column Data Sources for interacting with the plot
        hists = histogram(self.dfs, log=log, cum=cum)

        props = dict(x=[], y=[], population=[])
        tab.sources = [ColumnDataSource(props) for hist in hists]

        # Plot histogram values
        for i, source in enumerate(tab.sources):
            tab.fig.step(x='x',
                         y='y',
                         color=self.colours[i],
                         alpha=0.8,
                         legend='population',
                         line_width=2.5,
                         source=source,
                         mode='before')

        if kind == 'lin':
            self.hists_lin = hists
        elif kind == 'log':
            self.hists_log = hists
        elif kind == 'cum':
            self.hists_cum = hists

        self.tabs.append(tab)

    def update(self):
        """Update plots when interacted with."""
        for tab in self.tabs:
            x_abr = self.params[self.x_axis.value]
            y_abr = self.params[self.y_axis.value]

            tab.fig.xaxis.axis_label = self.x_axis.value
            tab.fig.yaxis.axis_label = self.y_axis.value

            if tab.name.startswith('Hist'):
                tab.fig.yaxis.axis_label = 'Fraction'

            for i, source in enumerate(tab.sources):

                cols = [f'{x_abr}_x', f'{x_abr}', 'population']

                if tab.name == 'Scatter':
                    cols = [x_abr, y_abr, 'color', 'population']
                    dfs = self.dfs
                elif tab.name == "Hist (Lin)":
                    dfs = self.hists_lin
                elif tab.name == "Hist (Log)":
                    dfs = self.hists_log
                elif tab.name == "Hist (Cum)":
                    dfs = self.hists_cum

                # Ensure columns are present in each population
                if (x_abr not in dfs[i] or
                   (tab.name == 'Scatter' and y_abr not in dfs[i])):
                    df = pd.DataFrame(np.nan, index=[0], columns=cols)
                else:
                    # Clean up data
                    df = dfs[i][cols]
                    df = df.replace('None', np.nan)
                    df = df.dropna()

                # Update data
                if tab.name == 'Scatter':
                    source.data = dict(x=df[x_abr],
                                       y=df[y_abr],
                                       color=df['color'],
                                       population=df['population'])
                else:
                    source.data = dict(x=df[f'{x_abr}_x'],
                                       y=df[f'{x_abr}'],
                                       population=df['population']
                                       )

    def set_layout(self):
        """Create the plot layout."""
        # What to interact with
        for control in [self.x_axis, self.y_axis]:
            control.on_change('value', lambda attr, old, new: self.update())

        # Set up sidebar
        cwd = os.path.dirname(__file__)

        def path(p):
            d = os.path.join(cwd, 'plot_config/{}.html'.format(p))
            return open(d).read()

        text_top = Div(text=path('text_top'))
        text_bottom = Div(text=path('text_bottom'))

        sidebar = [text_top, self.x_axis, self.y_axis, text_bottom]
        s = widgetbox(sidebar, width=380)

        # Set up tabs
        panels = []
        for tab in self.tabs:
            panels.append(Panel(child=tab.fig, title=tab.name))
            tab.fig.legend.click_policy = 'hide'
        tabs = Tabs(tabs=panels, width=self.width)

        # Add sidebar and tabs
        L = layout([[s, tabs]])

        # Initial load of data
        self.update()

        # Showtime
        curdoc().title = 'frbpoppy'
        curdoc().add_root(L)


# Parse system arguments
# (I know ArgumentParser is nicer, but bokeh only works with argv)
args = sys.argv

# Whether to plot the frbcat population
if '-frbcat' in args:
    frbcat = args[args.index('-frbcat') + 1]
    if frbcat == 'True':
        frbcat = True
    elif frbcat == 'False':
        frbcat = False
else:
    frcat = True

# Which files to plot
files = []
for a in args:
    a = a.strip('"')
    if a.endswith('.p'):
        files.append(a)

# Check whether populations have been given as input
if len(files) == 0:
    pprint('Nothing to plot: plot arguments are empty')
else:
    Plot(files=files, frbcat=frbcat)
