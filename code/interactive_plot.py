import numpy as np
import pandas as pd

from bokeh.plotting import figure
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource
from bokeh.models.widgets import Slider, Select
from bokeh.io import curdoc

di = pd.read_csv('/media/DATA/phd/frbpoppy/data/results/population_initial.csv')
di['color'] = 'grey'

ds = pd.read_csv('/media/DATA/phd/frbpoppy/data/results/population_wholesky.csv')
ds = ds[:500]
ds['color'] = 'purple'

# Combine populations to single dataframe, using updating the initial
# population if observation values have been found in the second one
df = pd.merge(di, ds, on=['gx', 'gy', 'gz'], how='outer', suffixes=('_old',''))
for c in df:
    if c + '_old' in df:
        df[c] = df[c].fillna(df[c + '_old'])
        del df[c + '_old']

axis_map = {
    'Distance (Gpc)': 'dist',
    'Dispersion Measure (pc/cm^3)': 'dm',
    'Dispersion Measure - Host (pc/cm^3)': 'dm_host',
    'Dispersion Measure - IGM (pc/cm^3)': 'dm_igm',
    'Dispersion Measure - Milky Way(pc/cm^3)': 'dm_mw',
    'Fluence (Jy*ms)': 'fluence',
    'Galactic Latitude (degrees)': 'gb',
    'Galactic Longitude (degrees)': 'gl',
    'Galactic X (Gpc)': 'gx',
    'Galactic Y (Gpc)': 'gy',
    'Galactic Z (Gpc)': 'gz',
    'Luminosity - Bolometric': 'lum_bol',
    'Peak Flux Density (W/m**2)': 's_peak',
    'Spectral Index': 'si',
    'Signal to Noise Ratio': 'snr',
    'Pulse Width - Effective (ms)': 'w_eff',
    'Pulse Width - Intrinsic (ms)': 'w_int',
    'Redshift': 'z',
    }

x_axis = Select(title='X-Axis',
                options=sorted(axis_map.keys()),
                value='Galactic X (Gpc)')

#x_log = Slider(title='X-Axis Log Scale', start=0, end=1, value=0, step=1)

y_axis = Select(title='Y-Axis',
                options=sorted(axis_map.keys()),
                value='Galactic Y (Gpc)')


# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[], y=[], color=[]))

# Create main plot
p = figure(plot_height=600, plot_width=600, title='', active_scroll='wheel_zoom',toolbar_location='above')
p.circle(x='x', y='y', source=source, size=7, color='color', line_color=None)


def update():
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    #p.x_axis_type= x_log.value

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value

    p.title.text = '%d sources' % len(df)

    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        color=df['color']
    )


controls = [x_axis, y_axis]

for control in controls:
    control.on_change('value', lambda attr, old, new: update())

sizing_mode = 'scale_width'

inputs = widgetbox(*controls, width=300)
l = layout([
    [inputs, p],
], sizing_mode=sizing_mode)

update()  # initial load of the data

curdoc().add_root(l)
curdoc().title = 'frbpoppy'
