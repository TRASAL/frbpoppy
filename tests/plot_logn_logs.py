"""Check the log N log S slope of a population."""
import numpy as np
from bokeh.plotting import figure, show

from frbpoppy.population import unpickle

MAKE = True

if MAKE:
    from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

    # Generate an FRB population
    days = 30
    population = CosmicPopulation(days*5000,
                                  z_max=0.1,
                                  lum_range=[1e36, 1e36],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='constant',
                                  days=days)
    population.name = 'cosmic'
    population.save()

    # Setup a survey
    survey = Survey('PERFECT', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'perfect'
    surv_pop.save()

else:
    surv_pop = unpickle('perfect')


f_lim = surv_pop.survey.calc_fluence_limit()
fluences = [f for f in surv_pop.get('fluence') if f >= f_lim]
fluences = surv_pop.get('fluence')

# Bin up
min_f = np.log10(min(fluences))
max_f = np.log10(max(fluences))
log_bins = np.logspace(min_f, max_f, 50)
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
       legend='perfect')

alpha, alpha_err, norm = surv_pop.calc_logn_logs()
print(alpha, alpha_err, norm)
f_lim = surv_pop.survey.calc_fluence_limit()

xs = 10**((np.log10(edges[:-1]) + np.log10(edges[1:])) / 2)
xs = [x for x in xs if x > f_lim]
ys = [norm*x**(alpha) for x in xs]
p.line(xs, ys, line_width=3, line_alpha=0.6, legend=f'α = {alpha:.3}')

show(p)
