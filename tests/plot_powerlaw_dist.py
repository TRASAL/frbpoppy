"""Test powerlaw distribution."""
import numpy as np
from bokeh.plotting import figure, show
from frbpoppy import powerlaw

ps = [powerlaw(1e38, 1e43, 1) for i in range(int(1e5))]

min_f = np.log10(min(ps))
max_f = np.log10(max(ps))
log_bins = np.logspace(min_f, max_f, 50)
hist, edges = np.histogram(ps, bins=log_bins)

p = figure(x_axis_type='log', y_axis_type='log')

p.quad(top=hist,
       bottom=1e-12,
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5)

show(p)
