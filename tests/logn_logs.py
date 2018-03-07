"""Check the log N log S slope of a population."""
import numpy as np
from scipy import optimize
from bokeh.plotting import figure, show

from frbpoppy.population import unpickle

pop_aper = unpickle('apertif')

# Get fluences
fluences = []
for src in pop_aper.sources:
    fluences.append(src.frbs[0].fluence)

# Bin up
log_bins = np.logspace(np.log10(min(fluences)), np.log10(max(fluences)), 10)
hist, edges = np.histogram(fluences, bins=log_bins)
cum_hist = [sum(hist[i:]) for i in range(len(hist))]

# Plot
p = figure(title='Log N / Log S', x_axis_type='log', y_axis_type='log')
p.xaxis.axis_label = 'Fluence (Jy*ms)'
p.yaxis.axis_label = 'Cumulative Number'

p.quad(top=cum_hist,
       bottom=10**(np.log10(min(hist)) - 1),
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5,
       legend='Apertif')

# Fit
from lmfit import Model

xs = 10**((np.log10(edges[:-1]) + np.log10(edges[1:])) / 2)


def powerlaw(x, norm=1, alpha=-1):
    return norm*x**(alpha)


model = Model(powerlaw)
params = model.make_params()

params['alpha'].min = -5.0
params['alpha'].max = -0.5

result = model.fit(cum_hist, params, x=xs)

alpha = result.params['alpha'].value
alpha_err = result.params['alpha'].stderr
norm = result.params['norm'].value
ys = [powerlaw(x, norm=norm, alpha=alpha) for x in xs]

p.line(xs, ys, line_width=3, line_alpha=0.6, legend=f'α = {alpha:.3}')

show(p)
