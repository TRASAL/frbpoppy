"""Replicate figure 3 of Niino (2018).

This figure compares input luminosity vs apparent luminosity for a variety of
cases.
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Adapt, Survey, SurveyPopulation
from frbpoppy import unpickle

CREATE = False
OBSERVE = False

if CREATE:
    days = 28
    n_per_day = 5000

    # Generate population with standard candles
    pop_sc = CosmicPopulation(n_per_day*days,
                              days=days,
                              name='sc',
                              dm_host_model='normal',
                              dm_host_mu=0,
                              dm_host_sigma=0,
                              dm_igm_index=0,
                              dm_mw_model='zero',
                              emission_range=[10e6, 10e9],
                              lum_range=[1e33, 1e33],
                              lum_index=0,
                              n_model='constant',
                              pulse_model='uniform',
                              pulse_range=[1., 1.],
                              pulse_mu=1.,
                              pulse_sigma=0.,
                              repeat=0.,
                              si_mu=0.,
                              si_sigma=0.,
                              z_max=0.1)
    pop_sc.save()

    # Generate population with a powerlaw
    pop_pl = copy.deepcopy(pop_sc)
    pop_pl = Adapt(pop_pl).lum_bol(1e33, 1e40, -2)
    pop_pl.name = 'pl'
    pop_pl.save()

if OBSERVE:
    if not CREATE:
        pop_sc = unpickle('sc')
        pop_pl = unpickle('pl')

    # Create Survey
    perfect = Survey('PERFECT-SMALL', gain_pattern='perfect', sidelobes=0)

    # Observe populations
    pop_sc_obs = SurveyPopulation(pop_sc, perfect)
    pop_sc_obs.name = 'sc-obs'
    pop_sc_obs.save()
    pop_pl_obs = SurveyPopulation(pop_pl, perfect)
    pop_pl_obs.name = 'pl-obs'
    pop_pl_obs.save()

else:
    pop_sc_obs = unpickle('sc-obs')
    pop_pl_obs = unpickle('pl-obs')


fig = plt.figure()
ax = fig.add_subplot(111)

for pop in (pop_sc_obs, pop_pl_obs):

    lum_org = np.array(pop.get('lum_bol'))

    # Create Survey
    perfect = Survey('PERFECT-SMALL', gain_pattern='airy', sidelobes=0)

    def multiply_with_int_pro(i):
        return i*perfect.intensity_profile(sidelobes=0)

    vec_int_pro = np.vectorize(multiply_with_int_pro)

    lum_app = vec_int_pro(lum_org)

    mini = min(lum_app)
    maxi = max(lum_org)
    log_bins = 10 ** np.linspace(np.log10(mini), np.log10(maxi), 50)

    color = next(ax._get_lines.prop_cycler)['color']
    # Original luminosity function
    n, bins = np.histogram(lum_org, bins=log_bins, density=False)
    bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
    plt.step(bincentres, n, where='mid', label=pop.name, color=color)

    # Apparent luminosity functions
    n, bins = np.histogram(lum_app, bins=log_bins, density=False)
    bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
    plt.step(bincentres, n, where='mid', label=f'{pop.name[:-3]}app',
             color=color,
             linestyle='--')

plt.xlabel('Luminosity')
plt.xscale('log')
plt.ylabel(r'Number')
plt.yscale('log')
plt.legend()

plt.tight_layout()
plt.savefig('plots/luminosity_functions.pdf')
