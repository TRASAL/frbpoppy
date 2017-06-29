"""Use Bokeh to plot a population or the results of a Monte Carlo."""

from bokeh.io import curdoc, output_file
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Panel, Tabs
from bokeh.models.widgets import RangeSlider, Select
from bokeh.palettes import Spectral11
from bokeh.plotting import figure
import os
import pandas as pd
import sqlite3


class Plot():
    """Class to plot a population or the output of a Monte Carlo."""

    def __init__(self, files=[]):
        """Initializing."""
        self.surveys = ['A', 'B']

        self.pars = {'dec': 'Declination (°)',
                     'dist': 'Distance (Gpc)',
                     'dm_host': 'Dispersion Measure - Host (pc/cm^3)',
                     'dm_igm': 'Dispersion Measure - IGM (pc/cm^3)',
                     'dm_mw': 'Dispersion Measure - Milky Way (pc/cm^3)',
                     'dm': 'Dispersion Measure (pc/cm^3)',
                     'fluence': 'Fluence (Jy*ms)',
                     'gb': 'Galactic Latitude (degrees)',
                     'gl': 'Galactic Longitude (degrees)',
                     'gx': 'Galactic X (Gpc)',
                     'gy': 'Galactic Y (Gpc)',
                     'gz': 'Galactic Z (Gpc)',
                     'lum_bol': 'Luminosity - Bolometric (10^30 ergs/s)',
                     'ra': 'Right Ascension (°)',
                     's_peak': 'Peak Flux Density (Jy)',
                     'snr': 'Signal to Noise Ratio',
                     'w_eff': 'Pulse Width - Effective (ms)',
                     'w_int': 'Pulse Width - Intrinsic (ms)',
                     'z': 'Redshift'}

        self.inv_pars = {v: k for k, v in self.pars.items()}

    def pop(self):
        """Plot a populations."""
        return 'TODO'

    def path(self, s):
        """Return the path to a file in the results folder."""
        return os.path.join(os.path.dirname(__file__), '../data/results/' + s)

    def import_df(self):
        """Import a sql database into a pandas dataframe."""
        conn = sqlite3.connect(self.path('ks.db'))
        df = pd.read_sql_query("select * from pars;", conn)
        cols = [c for c in df.columns if c.startswith('level')]
        df = df.drop(cols, axis=1)
        return df

    def mc(self):
        """Plot the results of a Monte Carlo run."""
        # Group parameters
        self.in_pars = ['w_int']
        self.out_pars = ['dec', 'dist', 'dm', 'fluence', 'ra', 'snr', 's_peak',
                         'w_eff', 'z']

        # Import goodness-of-fit data
        df = self.import_df()
        df = df[(df.survey != 'APERTIF') & (df.survey != 'WHOLESKY')]
        self.surveys = df.survey.unique().tolist()

        # Setup observed parameter choice
        out_opt = [self.pars[o] for o in self.out_pars]
        out_sel = Select(title="Observed parameter:",
                         value=out_opt[-2],
                         options=out_opt)

        # Setup input parameter choice
        in_opt = [self.pars[o] for o in self.in_pars]
        in_sel = Select(title="In parameter:",
                        value=in_opt[0],
                        options=in_opt)

        # Setup input parameters
        sliders = []
        for p in self.in_pars:
            # TODO convert min and max to values obtained from table
            t = self.pars[p]
            r = RangeSlider(start=0, end=10, range=(0, 9), step=.1, title=t)
            sliders.append(r)

        # Set up tools
        props = [("survey", "@survey")]
        hover = HoverTool(tooltips=props)
        sc_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create scatter plot
        sp = figure(plot_height=600,
                    plot_width=600,
                    active_scroll='wheel_zoom',
                    toolbar_location='right',
                    tools=sc_tools,
                    y_axis_type="log"
                    )

        # Create a Column Data Source for interacting with the plot
        props = dict(x=[], y=[], survey=[])
        sp_sources = [ColumnDataSource(props) for s in self.surveys]

        # Set colours
        colours = Spectral11[:len(self.surveys)]

        # Plot ks values various surveys
        for i, source in enumerate(sp_sources):
            sp.line(x='x',
                    y='y',
                    source=source,
                    line_width=5,
                    line_color=colours[i],
                    alpha=0.7,
                    legend='survey')

        # Interactive goodness
        def update():

            # Update axes
            x_name = self.inv_pars[in_sel.value]
            y_name = self.inv_pars[out_sel.value]
            sp.xaxis.axis_label = in_sel.value
            sp.yaxis.axis_label = 'P-value ({})'.format(out_sel.value)

            # Update data
            for i, source in enumerate(sp_sources):

                # Refilter data
                name = self.surveys[i]
                test = (df.survey == name)
                length = df[test].shape[0]
                x = df[test][x_name + '_max'].tolist()
                y = df[test]['ks_' + y_name].tolist()
                s = [name for j in range(length)]

                source.data = dict(x=x, y=y, survey=s)

        # What to interact with
        controls = [out_sel, in_sel]

        for control in controls:
            control.on_change('value', lambda attr, old, new: update())

        # Layout options
        sizing_mode = 'fixed'
        sidebar = [in_sel, out_sel, *sliders]
        s = widgetbox(sidebar, width=350)
        tab1 = Panel(child=sp, title='scatter')
        tabs = Tabs(tabs=[tab1])
        L = layout([[s, tabs]], sizing_mode=sizing_mode)

        # Make legend clickable
        sp.legend.click_policy = 'hide'

        update()  # Initial load of the data
        curdoc().add_root(L)
        curdoc().title = 'frbpoppy'


Plot().mc()
