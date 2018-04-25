"""Use Bokeh to plot a population or the results of a Monte Carlo."""

from bokeh.io import curdoc, show
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div, Span
from bokeh.models.widgets import CheckboxButtonGroup, RangeSlider
from bokeh.models.widgets import Select, Slider
from bokeh.palettes import Spectral11
from bokeh.plotting import figure

import os
import numpy as np
import pandas as pd
import sqlite3

from frbpoppy.paths import paths
from frbpoppy.log import pprint


class Plot:
    """Class to plot a population or the output of a Monte Carlo."""

    def __init__(self, files=[]):
        """Initializing."""
        self.pars = {'dec': 'Declination (°)',
                     'dist': 'Distance (Gpc)',
                     'dm_host': 'Dispersion Measure - Host (pc/cm^3)',
                     'dm_igm': 'Dispersion Measure - IGM (pc/cm^3)',
                     'dm_igm_slope': 'Dispersion Measure - IGM - Slope (pc/cm^3)',
                     'dm_mw': 'Dispersion Measure - Milky Way (pc/cm^3)',
                     'dm': 'Dispersion Measure (pc/cm^3)',
                     'fluence': 'Fluence (Jy*ms)',
                     'freq_min': 'Source frequency - minimum (Hz)',
                     'freq_max': 'Source frequency - maximum (Hz)',
                     'gb': 'Galactic Latitude (degrees)',
                     'gl': 'Galactic Longitude (degrees)',
                     'gx': 'Galactic X (Gpc)',
                     'gy': 'Galactic Y (Gpc)',
                     'gz': 'Galactic Z (Gpc)',
                     'lum_bol': 'Luminosity - Bolometric (10^30 ergs/s)',
                     'lum_bol_min': 'Luminosity - Bolometric - Minimum (10^30 ergs/s)',
                     'lum_bol_max': 'Luminosity - Bolometric - Maximum (10^30 ergs/s)',
                     'lum_bol_slope': 'Luminosity - Bolometric - Slope (10^30 ergs/s)',
                     'n_day': 'Number of seed sources (#/day)',
                     'ra': 'Right Ascension (°)',
                     'rep': 'Repeater Fraction',
                     's_peak': 'Peak Flux Density (Jy)',
                     'si_mean': 'Spectral Index - Mean',
                     'si_sigma': 'Spectral Index - Sigma',
                     'snr': 'Signal to Noise Ratio',
                     'w_eff': 'Pulse Width - Effective (ms)',
                     'w_int': 'Pulse Width - Intrinsic (ms)',
                     'w_int_min': 'Pulse Width - Intrinsic - Minimum (ms)',
                     'w_int_max': 'Pulse Width - Intrinsic - Maximum (ms)',
                     'z': 'Redshift'}

        self.inv_pars = {v: k for k, v in self.pars.items()}

    def path(self, s, where='results'):
        """Return the path to a file in the results folder."""
        if where == 'results':
            r = os.path.join(paths.results(), s)
        if where == 'html':
            cd = paths.code
            d = os.path.join(cd, f'plot_config/{s}.html')
            r = open(d).read()
        return r

    def import_df(self, query=None, loc=None):
        """Import a sql database into a pandas dataframe."""
        # Options DataFrame
        if not query:
            query = 'select * from pars;'
        if not loc:
            loc = 'ks_test_small.db'
        elif loc == 'ks':
            loc = 'ks_test_small.db'
        elif loc == 'hist':
            loc = 'hists_test_small.db'
        elif loc == 'rate':
            loc = 'rates_test_small.db'

        p = self.path(loc)
        conn = sqlite3.connect(p)
        df = pd.read_sql_query(query, conn)

        # Clean up DataFrame
        cols = [c for c in df.columns if c.startswith('level')]
        df = df.drop(cols, axis=1)
        df = df.fillna('')

        # Only return entries if present in frbcat
        # TODO Check whether this is still relevant
        df = df[(df.survey != 'APERTIF') & (df.survey != 'WHOLESKY')]

        return df

    def mc(self):
        """Plot the results of a Monte Carlo run."""
        # Import goodness-of-fit data
        df = self.import_df()
        # Ugly fix due to rouuding errors
        df['w_int_max'] = round(df['w_int_max'], 8)
        self.surveys = sorted(df.survey.unique().tolist())

        # Group parameters
        self.in_pars = df['in_par'].unique().tolist()

        self.out_pars = ['dec', 'dist', 'dm', 'fluence', 'ra', 'snr', 's_peak',
                         'w_eff', 'z']

        # ----- Set widget options -----

        # Optional frbcat plotting
        cat_sel = CheckboxButtonGroup(labels=['frbcat'], active=[0])

        # Optional survey plotting
        if len(self.surveys) >= 7:
            selection = [6, 7]
        else:
            selection = [i for i in range(len(self.surveys))]
        sur_sel = CheckboxButtonGroup(labels=self.surveys, active=selection)

        # Setup observed parameter choice
        out_opt = [self.pars[o] for o in self.out_pars]
        out_sel = Select(title="Observed parameter:",
                         value=out_opt[-2],
                         options=out_opt)

        # Setup input parameter choice
        in_opt = [self.pars[o] for o in self.in_pars]
        in_sel = Select(title="In parameter:",
                        value=in_opt[-1],
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
            if len(s) > 1:
                st = s[1] - s[0]
            else:
                st = 0.00001
                ma += st

            pprint(f"{p:>15} {mi:>15} {ma:>15} {st:>15}")

            # Set standard value
            v = df[p].value_counts().idxmax()

            r = Slider(start=mi,
                       end=ma,
                       step=st,
                       value=v,
                       title=p,
                       callback_policy="mouseup")

            single_sel.append(r)

        for p in dual_opt:

            title = p

            # For the logrithmic sliders
            if p == 'freq' or p == 'lum_bol':
                p_min = np.log10(df[p + '_min'])
                p_max = np.log10(df[p + '_max'])
            else:
                p_min = df[p + '_min']
                p_max = df[p + '_max']

                # title += ' (10^x)'

            # Minimum and maximum values
            mi = min(p_min)
            ma = max(p_max)

            # Stepsize
            s = np.sort(p_min.unique())
            if len(s) > 1:
                st = s[1] - s[0]
            else:
                st = 0.00001
                ma += st

            # Set standard values
            mi_v = p_min.value_counts().idxmax()
            ma_v = p_max.value_counts().idxmax()

            r = RangeSlider(start=mi,
                            end=ma,
                            value=(mi_v, ma_v),
                            step=st,
                            title=title,
                            callback_policy="mouseup")
            dual_sel.append(r)

        #  ----- Set plot options -----

        # Left most plot with ks-test
        # --------------------------
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
            lp.circle(x='x',
                      y='y',
                      source=source,
                      size=5,
                      color=colours[i],
                      alpha=0.7,
                      legend='survey')

        # Add a line to link to the histogram plots
        span = Span(location=0,
                    dimension='height', line_color='green',
                    line_dash='dashed', line_width=3)
        lp.add_layout(span)

        # Centre plot with Histograms
        # ---------------------------
        # Set up tools
        props = [("survey", "@survey"), ("frbs", "@frbs")]
        hover = HoverTool(tooltips=props)
        hp_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create histogram plot
        hp = figure(plot_height=600,
                    plot_width=600,
                    active_scroll='wheel_zoom',
                    toolbar_location='right',
                    tools=hp_tools)

        # Create Column Data Sources for interacting with the plot
        props = dict(top=[], left=[], right=[], bottom=[], survey=[], frbs=[])
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

        # Right most plot with rates
        # --------------------------
        # Set up tools
        props = [("Survey", '@surveys'), ("Value", '@frbs_day')]
        hover = HoverTool(tooltips=props)
        rp_tools = ['box_zoom', 'pan', 'save', hover, 'reset', 'wheel_zoom']

        # Create bar plot
        rp = figure(plot_height=600,
                    plot_width=600,
                    active_scroll='wheel_zoom',
                    toolbar_location='right',
                    tools=rp_tools,
                    y_range=self.surveys,
                    x_axis_type="log"
                    )

        # Create a Column Data Source for interacting with the plot
        props = dict(surveys=[], srcs_day=[], frbs_day=[], baseline=[])
        rp_source = ColumnDataSource(props)

        rp.hbar(y='surveys',
                left='baseline',
                right='frbs_day',
                height=0.7,
                source=rp_source)

        rp.xaxis.axis_label = 'FRBs per day'
        rp.y_range.range_padding = 0.1
        rp.ygrid.grid_line_color = None

        # Interactive goodness
        def update():

            # Get values from sliders
            val = {}
            for s in single_sel:
                val[s.title] = s.value
            for s in dual_sel:
                val[s.title + '_min'] = s.value[0]
                val[s.title + '_max'] = s.value[1]

            # Get active surveys
            act_sur = [self.surveys[i] for i in sur_sel.active]

            # Update axes
            x_name = self.inv_pars[in_sel.value]
            y_name = self.inv_pars[out_sel.value]
            lp.xaxis.axis_label = in_sel.value
            lp.yaxis.axis_label = 'P-value ({})'.format(out_sel.value)
            span.location = val[x_name]

            # For the logarithmic sliders
            for v in val:
                if ('freq' in v or 'lum_bol' in v) and 'slope' not in v:
                    val[v] = float(10**val[v])

            # Update data
            for i, source in enumerate(lp_sources):

                # Refilter data
                survey = self.surveys[i]
                filt = val.copy()
                filt['survey'] = survey

                del filt[x_name]

                query = 'SELECT * FROM pars WHERE '

                # An utterly and completely ridiculous way to solve SQL
                # rounding bug with decimals
                for f in filt:
                    value = str(filt[f])
                    if '.' not in value:
                        if 'e' in value:
                            n, p = value.split('e')
                            value = n + '.0e' + p
                        if f in ['w_int_max', 'w_int_min']:
                            value += '.0'
                    if f == 'si_sigma':
                        if value != '0.0' and len(value.split('.')[-1]) != 2:
                            value += '0'
                    if f == 'si_mean':
                        if len(value.split('.')[-1]) == 2:
                            value = value[:-1]
                    if f in ['freq_max']:
                        pre, post = value.split('.')
                        value = str(round(float(value[:18]), 15-len(pre)))

                    query += "{} LIKE '{}' AND ".format(f, value)

                query = query[:-4]  # Remove the last AND

                dk = self.import_df(query=query)

                x = dk[x_name].tolist()
                y = dk['ks_' + y_name].tolist()
                length = dk.shape[0]
                s = [survey for j in range(length)]

                # Prevent superfluous graphs
                if all(e == '' for e in y) or survey not in act_sur:
                    x = y = s = []

                source.data = dict(x=x, y=y, survey=s)

            # Update axes
            hp.xaxis.axis_label = out_sel.value

            # Get the id of the set of input values
            err = 0.0000001
            val_up = val[x_name] + err
            val_down = val[x_name] - err
            val_test = ((val_down <= dk[x_name]) & (dk[x_name] <= val_up))
            par_test = (dk['in_par'] == x_name)
            test = (val_test & par_test)
            try:
                iden = dk[test].iloc[0]['id']
                pprint('Parameters found')
            except IndexError:
                pprint('No parameters found')
                pprint('SQL constraints:', query)
                pprint('Looking for:', x_name)
                pprint('With value:', val[x_name])
                pprint('Finding:', *[repr(r) for r in dk[x_name].values])
                pprint('Looking in:\n', dk)
                iden = ''

            # Update histogram data
            for i, source in enumerate(hp_sources):

                # Import data
                name = self.surveys[i]
                query = "SELECT * FROM pars WHERE id='{}' and "
                query += "in_par='{}' and out_par='{}' and "
                query += "survey='{}';"
                query = query.format(iden, x_name, y_name, name)
                dh = self.import_df(query=query, loc='hist')

                if 'frbs' not in dh:
                    dh['frbs'] = '?'

                dc = dh[(dh.frbcat == 1)]  # 1: frbcat population (dashed)
                dh = dh[(dh.frbcat == 0)]  # 0: virtual population (solid)

                # Prepare data for plotting
                top = dh.top.tolist()
                bottom = dh.bottom.tolist()
                left = dh.left.tolist()
                right = dh.right.tolist()
                length = dh.shape[0]
                s = [name for j in range(length)]
                frbs = dh.frbs.tolist()

                # Prevent superfluous graphs
                if name not in act_sur:
                    top = []
                    bottom = []
                    left = []
                    right = []
                    s = []
                    frbs = []

                # Plot data
                source.data = dict(top=top,
                                   bottom=bottom,
                                   left=left,
                                   right=right,
                                   survey=s,
                                   frbs=frbs)

                cat_source = cat_sources[i]

                # Prepare frbcat data for plotting
                top = dc.top.tolist()
                bottom = dc.bottom.tolist()
                left = dc.left.tolist()
                right = dc.right.tolist()
                length = dc.shape[0]
                s = [name + ' (🐱)' for j in range(length)]
                frbs = dc.frbs.tolist()

                # Prevent superfluous graphs
                if name not in act_sur or not cat_sel.active:
                    top = bottom = left = right = s = frbs = []

                # Plot data
                cat_source.data = dict(top=top,
                                       bottom=bottom,
                                       left=left,
                                       right=right,
                                       survey=s,
                                       frbs=frbs)

            # Update rate plot
            # Construct query to import sql database
            query = "SELECT * FROM pars WHERE id='{}' and in_par='{}';"
            query = query.format(iden, x_name)
            dr = self.import_df(query=query, loc='rate')

            frb_rates = dr.frbs_per_day.tolist()
            survs = dr.survey.tolist()
            min_frbs = min(dr.frbs_per_day.tolist()) / 10.
            baseline = [min_frbs for r in frb_rates]

            # Update source
            rp_source.data = dict(surveys=survs,
                                  frbs_day=frb_rates,
                                  baseline=baseline)

        # What to interact with
        for control in [in_sel, out_sel, *single_sel, *dual_sel]:
            control.on_change('value', lambda attr, old, new: update())
        for checkbox in [cat_sel, sur_sel]:
            checkbox.on_change('active', lambda attr, old, new: update())

        # Get text
        text_top = Div(text=self.path('mc_top', where='html'))
        text_bottom = Div(text=self.path('mc_bottom', where='html'))

        # Layout options
        sizing_mode = 'fixed'
        left_sb = widgetbox([cat_sel, sur_sel, in_sel, out_sel], width=600)
        centre_sb = widgetbox(single_sel, width=600)
        right_sb = widgetbox(dual_sel, width=600)

        grid = [[lp, hp, rp], [left_sb, centre_sb, right_sb]]
        page = layout(grid, sizing_mode=sizing_mode)

        # Make legend clickable
        lp.legend.click_policy = 'hide'
        hp.legend.click_policy = 'hide'

        update()  # Initial load of the data
        curdoc().add_root(page)
        curdoc().title = 'frbpoppy'


# Plot().mc()
