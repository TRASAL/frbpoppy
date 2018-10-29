import numpy as np
from bokeh.plotting import figure, show

from frbpoppy import unpickle

pop = unpickle('apertif')

w_effs = pop.get('w_eff')
t_dms = pop.get('t_dm')

ratios = [e[0]/e[1] for e in zip(w_effs, t_dms)]

# Bin up
min_f = np.log10(min(ratios))
max_f = np.log10(max(ratios))
log_bins = np.logspace(min_f, max_f, 50)
hist, edges = np.histogram(ratios, bins=log_bins)
cum_hist = [sum(hist[i:]) for i in range(len(hist))]

# Plot
p = figure(title='Ratio tau / tau_dm', x_axis_type='log', y_axis_type='log')
p.yaxis.axis_label = 'Fraction'

# Find optimum lower limit
m = min(e for e in cum_hist if e > 0)
bottom = 10**(np.log10(m) - 1)

p.quad(top=hist,
       bottom=bottom,
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5,
       legend='Apertif')

show(p)
