"""Check the log N log F slope of a population."""
import numpy as np
import matplotlib.pyplot as plt
import os

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy.population import unpickle

MAKE = True


if MAKE:
    population = CosmicPopulation.simple(1e5, generate=True)
    survey = Survey('perfect')
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'lognlogflocal'
    surv_pop.save()
else:
    surv_pop = unpickle('lognlogflocal')

# Get parameter
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

# Change working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

plt.style.use('./aa.mplstyle')
fig = plt.figure()
ax = fig.add_subplot(111)

plt.step(edges[:-1], n_gt_s, where='post')
plt.plot(xs, ys, linestyle='--',
         label=rf'$\alpha$ = {alpha:.3} $\pm$ {round(abs(alpha_err), 2)}')

plt.xlabel('Fluence (Jy ms)')
plt.ylabel(r'N(${>}\text{S}_{\text{min}}$)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plots/logn_logf_local.pdf')
