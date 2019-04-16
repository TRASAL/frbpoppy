"""Check the log N log S slope of a population."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy.population import unpickle

MAKE = False

if MAKE:

    # Generate an FRB population
    days = 14
    population = CosmicPopulation(days*5000,
                                  z_max=0.01,
                                  lum_range=[1e40, 1e40],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='vol_co',
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
                                  pulse_sigma=0.)
    population.name = 'test'

    # Setup a survey
    survey = Survey('perfect', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'lognlogslocal'
    surv_pop.save()

else:
    surv_pop = unpickle('lognlogslocal')


parms = surv_pop.frbs.fluence
min_p = min(parms)
max_p = max(parms)

# Bin up
min_f = np.log10(min(parms))
max_f = np.log10(max(parms))
log_bins = np.logspace(min_f, max_f, 50)
hist, edges = np.histogram(parms, bins=log_bins)
n_gt_s = np.cumsum(hist[::-1])[::-1]

# Calculate alpha
alpha, alpha_err, norm = surv_pop.calc_logn_logs(parameter='fluence',
                                                 min_p=min_p,
                                                 max_p=max_p)

print(alpha, alpha_err, norm)
xs = 10**((np.log10(edges[:-1]) + np.log10(edges[1:])) / 2)
xs = xs[xs >= min_p]
xs = xs[xs <= max_p]
ys = [norm*x**(alpha) for x in xs]


fig = plt.figure()
ax = fig.add_subplot(111)

plt.step(edges[:-1], n_gt_s, where='post')
plt.plot(xs, ys, linestyle='--',
         label=rf'$\alpha$ = {alpha:.3} ± {round(abs(alpha_err), 2)}')

plt.xlabel('S (Jy ms)')
plt.ylabel('N (>S)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plots/logn_logs_local.pdf')
