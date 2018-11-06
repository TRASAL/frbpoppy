"""Check the log N log S slope of a population."""
import numpy as np
from bokeh.plotting import figure, show

from frbpoppy.population import unpickle

MAKE = False

if MAKE:
    from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

    # Generate an FRB population
    days = 14
    population = CosmicPopulation(days*5000,
                                  z_max=0.01,
                                  lum_range=[1e40, 1e40],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='constant',
                                  days=days,
                                  dm_host_model='normal',
                                  dm_host_mu=0,
                                  dm_host_sigma=0,
                                  dm_igm_index=0,
                                  dm_igm_sigma=0,
                                  dm_mw_model='zero',
                                  emission_range=[10e6, 10e9],
                                  lum_index=0,
                                  pulse_model='uniform',
                                  pulse_range=[1., 1.],
                                  pulse_mu=1.,
                                  pulse_sigma=0.,
                                  repeat=0.)
    population.name = 'test'
    population.save()

    # Setup a survey
    survey = Survey('PERFECT', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'lognlogs'
    surv_pop.save()

else:
    surv_pop = unpickle('lognlogs')


parms = np.array(surv_pop.get('fluence'))
min_p = min(parms)
max_p = max(parms)

# Bin up
min_f = np.log10(min(parms))
max_f = np.log10(max(parms))
log_bins = np.logspace(min_f, max_f, 50)
hist, edges = np.histogram(parms, bins=log_bins)
cum_hist = [sum(hist[i:]) for i in range(len(hist))]

# Plot
p = figure(title='Log N / Log S', x_axis_type='log', y_axis_type='log')
p.xaxis.axis_label = 'Fluence (Jy ms)'
p.yaxis.axis_label = 'Cumulative Number'

# Find optimum lower limit
m = min(e for e in cum_hist if e > 0)
bottom = 10**(np.log10(m) - 1)

p.quad(top=cum_hist,
       bottom=bottom,
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5,
       legend='perfect')

alpha, alpha_err, norm = surv_pop.calc_logn_logs(parameter='fluence',
                                                 min_p=min_p,
                                                 max_p=max_p)
print(alpha, alpha_err, norm)

xs = 10**((np.log10(edges[:-1]) + np.log10(edges[1:])) / 2)
# import IPython; IPython.embed()
xs = xs[xs>=min_p]
xs = xs[xs<=max_p]
ys = [norm*x**(alpha) for x in xs]
p.line(xs, ys, line_width=3, line_alpha=0.6,
       legend=f'α = {alpha:.3} ± {round(abs(alpha_err), 2)}')

show(p)
