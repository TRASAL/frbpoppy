"""Check the log N log S slope of a population."""
import numpy as np
from bokeh.plotting import figure, show

from frbpoppy.population import unpickle

MAKE = False
FIT = False
POP = 'lum_neg'
FLUENCE_LIMIT_LOW = False  # False or a float value

if MAKE:
    from frbpoppy.do_populate import generate
    from frbpoppy.do_survey import observe

    days = 7
    n_per_day = 5000

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_dist_pars=[1e40, 1e50, -1.5],
                          z_max=2.5,
                          pulse=[0.1, 10],
                          repeat=0.0)

    # Observe FRB population
    pop_aper = observe(population, 'APERTIF')

else:
    pop_aper = unpickle(POP)

# Get fluences
fluences = []
for src in pop_aper.sources:
    f = src.frbs[0].fluence
    if FLUENCE_LIMIT_LOW:
        if f > 10**(np.log10(FLUENCE_LIMIT_LOW) - 1):
            fluences.append(f)
    else:
        fluences.append(f)

# Bin up
log_bins = np.logspace(np.log10(min(fluences)), np.log10(max(fluences)), 50)
hist, edges = np.histogram(fluences, bins=log_bins)
cum_hist = [sum(hist[i:]) for i in range(len(hist))]

# Plot
p = figure(title='Log N / Log S', x_axis_type='log', y_axis_type='log')
p.xaxis.axis_label = 'Fluence (Jy*ms)'
p.yaxis.axis_label = 'Cumulative Number'

# Find optimum lower limit
m = min(e for e in cum_hist if e > 0)
bottom = 10**(np.log10(m) - 1)

p.quad(top=cum_hist,
       bottom=bottom,
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5,
       legend=POP)

if FIT:
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
