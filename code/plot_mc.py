"""Use Bokeh to plot a population or the results of a Monte Carlo."""

from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div, Panel, Span, Tabs
from bokeh.models.widgets import CheckboxButtonGroup, RangeSlider
from bokeh.models.widgets import Select, Slider
from bokeh.palettes import Spectral11
from bokeh.plotting import figure

import os
import numpy as np
import pandas as pd
import sqlite3


class Plot:
    """Class to plot a population or the output of a Monte Carlo."""

    def __init__(self, files=[]):
        """Initializing."""
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
                     'w_int_min': 'Pulse Width - Intrinsic - Minimum (ms)',
                     'w_int_max': 'Pulse Width - Intrinsic - Maximum (ms)',
                     'z': 'Redshift'}

        self.inv_pars = {v: k for k, v in self.pars.items()}

    def pop(self, files=[]):
        """Plot a population."""
        return 'TODO'

    def path(self, s, where='results'):
        """Return the path to a file in the results folder."""
        if where == 'results':
            r = os.path.join(os.path.dirname(__file__), '../data/results/' + s)
        if where == 'html':
            cwd = os.path.dirname(__file__)
            d = os.path.join(cwd, 'plot_config/{}.html'.format(s))
            r = open(d).read()
        return r

    def import_df(self, query=None):
        """Import a sql database into a pandas dataframe."""
        # Options DataFrame
        if query:
            p = self.path('hists.db')
        else:
            query = 'select * from pars;'
            p = self.path('ks_3.db')

        # Get DataFrame
        conn = sqlite3.connect(p)
        df = pd.read_sql_query(query, conn)

        # Clean up DataFrame
        cols = [c for c in df.columns if c.startswith('level')]
        df = df.drop(cols, axis=1)
        df = df.fillna('')

        # Only return entries if present in frbcat
        df = df[(df.survey != 'APERTIF') & (df.survey != 'WHOLESKY')]

        return df

    def mc(self):
        """Plot the results of a Monte Carlo run."""
        # Import goodness-of-fit data
        df = self.import_df()
        self.surveys = df.survey.unique().tolist()

        # Group parameters
        self.in_pars = ['w_int_min', 'w_int_max', 'dm_host']
        self.out_pars = ['dec', 'dist', 'dm', 'fluence', 'ra', 'snr', 's_peak',
                         'w_eff', 'z']

        # ----- Set widget options -----

        # Optional frbcat plotting
        cat_sel = CheckboxButtonGroup(labels=['frbcat'], active=[0])
        # Optional survey plotting
        sur_sel = CheckboxButtonGroup(labels=self.surveys, active=[-1])

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
        single_opt = []
        dual_opt = []
        single_sel = []
        dual_sel = []
        exempt = ['index', 'survey', 'telescope', 'path', 'in_par', 'id']

        for c in df:
            if c.endswith('min'):
                dual_opt.append(c[:-4])
            elif c.endswith('max'):
                continue
            elif c.startswith('ks_'):
                exempt.append(c)
            elif c not in exempt:
                single_opt.append(c)

        for p in single_opt:

            # Minimum and maximum values
            mi = min(df[p])
            ma = max(df[p])

            # Stepsize
            s = np.sort(df[p].unique())
            st = s[1] - s[0]

            # Set standard value
            v = df[p].value_counts().idxmax()

            r = Slider(start=mi, end=ma, step=st, value=v, title=p)
            single_sel.append(r)

        for p in dual_opt:

            # Minimum and maximum values
            mi = min(df[p + '_min'])
            ma = max(df[p + '_max'])

            # Stepsize
            s = np.sort(df[p + '_min'].unique())
            st = s[1] - s[0]

            # Set standard values
            mi_v = df[p + '_min'].value_counts().idxmax()
            ma_v = df[p + '_max'].value_counts().idxmax()

            r = RangeSlider(start=mi,
                            end=ma,
                            range=(mi_v, ma_v),
                            step=st,
                            title=p)
            dual_sel.append(r)

        #  ----- Set plot options -----

        # Set up tools
        props = [("survey", "@survey")]
        hover = HoverTool(tooltips=props)
        sc_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create scatter plot
        lp = figure(plot_height=600,
                    plot_width=600,
                    active_scroll='wheel_zoom',
                    toolbar_location='right',
                    tools=sc_tools,
                    y_axis_type="log"
                    )

        # Create a Column Data Source for interacting with the plot
        props = dict(x=[], y=[], survey=[])
        lp_sources = [ColumnDataSource(props) for s in self.surveys]

        # Set colours
        colours = Spectral11[:len(self.surveys) + 1]

        # Plot ks values for various surveys
        for i, source in enumerate(lp_sources):
            lp.line(x='x',
                    y='y',
                    source=source,
                    line_width=5,
                    line_color=colours[i],
                    alpha=0.7,
                    legend='survey')

        # Add a line to link to the histogram plots
        span = Span(location=0,
                    dimension='height', line_color='green',
                    line_dash='dashed', line_width=3)
        lp.add_layout(span)

        # Set up tools
        props = [("survey", "@survey")]
        hover = HoverTool(tooltips=props)
        hp_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create histogram plot
        hp = figure(plot_height=600,
                    plot_width=600,
                    active_scroll='wheel_zoom',
                    toolbar_location='right',
                    tools=hp_tools)

        # Create Column Data Sources for interacting with the plot
        props = dict(top=[], left=[], right=[], bottom=[], survey=[])
        opt = [*self.surveys].append('frbcat')
        hp_sources = [ColumnDataSource(props) for s in self.surveys]
        cat_sources = [ColumnDataSource(props) for s in self.surveys]

        # Plot histogram values
        for i, source in enumerate(hp_sources):
            hp.quad(bottom='bottom',
                    left='left',
                    right='right',
                    top='top',
                    color=colours[i],
                    line_color=colours[i],
                    line_dash='solid',
                    legend='survey',
                    alpha=0.4,
                    source=source)

        for i, source in enumerate(cat_sources):
            hp.quad(bottom='bottom',
                    left='left',
                    right='right',
                    top='top',
                    color=colours[i],
                    line_color=colours[i],
                    line_dash='dashed',
                    legend='survey',
                    line_width=5,
                    line_alpha=0.6,
                    alpha=0,
                    source=source)

        # Interactive goodness
        def update():

            # Get values from sliders
            val = {}
            for s in single_sel:
                val[s.title] = s.value
            for s in dual_sel:
                val[s.title + '_min'] = s.range[0]
                val[s.title + '_max'] = s.range[1]

            # Update axes
            x_name = self.inv_pars[in_sel.value]
            y_name = self.inv_pars[out_sel.value]
            lp.xaxis.axis_label = in_sel.value
            lp.yaxis.axis_label = 'P-value ({})'.format(out_sel.value)
            span.location = val[x_name]

            # Update data
            for i, source in enumerate(lp_sources):

                # Refilter data
                survey = self.surveys[i]
                filt = val.copy()
                filt['survey'] = survey
                filt['in_par'] = x_name
                del filt[x_name]

                t = [df[f] == filt[f] for f in filt]
                df['test'] = pd.DataFrame(t).all()

                x = df[(df.test == True)][x_name].tolist()
                y = df[(df.test == True)]['ks_' + y_name].tolist()

                length = df[(df.test == True)].shape[0]
                s = [survey for j in range(length)]

                source.data = dict(x=x, y=y, survey=s)

            # Update axes
            hp.xaxis.axis_label = out_sel.value

            # Get the id of the set of input values
            test = ((df.test == True) & (df[x_name] == val[x_name]))
            iden = df[test].iloc[0]['id']

            # Update histogram data
            for i, source in enumerate(hp_sources):

                # Import data
                name = self.surveys[i]
                query = "SELECT * FROM pars WHERE id='{}' and in_par='{}'"
                query += " and out_par='{}' and survey='{}';"
                query = query.format(iden, x_name, y_name, name)
                dh = self.import_df(query=query)
                # 1: frbcat population (solid), 0: virtual population (dashed)
                dc = dh[(dh.frbcat == 1)]
                dh = dh[(dh.frbcat == 0)]

                # Prepare data for plotting
                top = dh.top.tolist()
                bottom = dh.bottom.tolist()
                left = dh.left.tolist()
                right = dh.right.tolist()
                length = dh.shape[0]
                s = [name for j in range(length)]

                # Plot data
                source.data = dict(top=top,
                                   bottom=bottom,
                                   left=left,
                                   right=right,
                                   survey=s)

                cat_source = cat_sources[i]

                # Prepare data for plotting
                top = dc.top.tolist()
                bottom = dc.bottom.tolist()
                left = dc.left.tolist()
                right = dc.right.tolist()
                length = dc.shape[0]
                s = [name + ' (🐱)' for j in range(length)]

                # Plot data
                cat_source.data = dict(top=top,
                                       bottom=bottom,
                                       left=left,
                                       right=right,
                                       survey=s)

        # What to interact with
        for control in [in_sel, out_sel, *single_sel]:
            control.on_change('value', lambda attr, old, new: update())
        for control in dual_sel:
            control.on_change('range', lambda attr, old, new: update())
        for checkbox in [cat_sel, sur_sel]:
            checkbox.on_change('active', lambda attr, old, new: update())

        # Get text
        text_top = Div(text=self.path('mc_top', where='html'))
        text_bottom = Div(text=self.path('mc_bottom', where='html'))

        # Layout options
        sizing_mode = 'fixed'
        sidebar = [cat_sel, sur_sel, in_sel, out_sel, *single_sel, *dual_sel]
        s = widgetbox(sidebar, width=350)
        tab1 = Tabs(tabs=[Panel(child=lp, title='p-value')])
        tab2 = Tabs(tabs=[Panel(child=hp, title='Histograms')])
        L = layout([[s, lp, hp]], sizing_mode=sizing_mode)

        # Make legend clickable
        lp.legend.click_policy = 'hide'
        hp.legend.click_policy = 'hide'

        update()  # Initial load of the data
        curdoc().add_root(L)
        curdoc().title = 'frbpoppy'


Plot().mc()
