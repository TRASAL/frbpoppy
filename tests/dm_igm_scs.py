"""Replicate middle plot of figure 3 of Niino (2018).

This figure compares input luminosity vs apparent luminosity for a variety of
cases.
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import unpickle
import frbpoppy.distributions as dis

CREATE = False
OBSERVE = True

pop = {}

if CREATE:
    days = 56
    n_per_day = 5000

    for n in (35, 36, 37):
        p = float(f'1e{n}')

        if n == 35:
            # Generate population with standard candles
            pop[n] = CosmicPopulation(n_per_day*days,
                                      days=days,
                                      name=f'sc-{n}',
                                      dm_host_model='normal',
                                      dm_host_mu=0,
                                      dm_host_sigma=0,
                                      dm_igm_index=1200,
                                      dm_igm_sigma=0,
                                      dm_mw_model='zero',
                                      emission_range=[10e6, 10e9],
                                      lum_range=[p, p],
                                      lum_index=0,
                                      n_model='sfr',
                                      pulse_model='uniform',
                                      pulse_range=[1., 1.],
                                      pulse_mu=1.,
                                      pulse_sigma=0.,
                                      repeat=0.,
                                      si_mu=0.,
                                      si_sigma=0.,
                                      z_max=2.5)
            pop[n].save()

        else:
            pop[n] = copy.deepcopy(pop[35])
            n_gen = len(pop[35].frbs.lum_bol)
            pop[n].frbs.lum_bol = dis.powerlaw(10**n, 10**n, 0, n_gen)
            pop[n].name = f'sc-{n}'
            pop[n].save()

pop_obs = {}

if OBSERVE:

    for n in (35, 36, 37):

        if not CREATE:
            pop[n] = unpickle(f'sc-{n}')

        # Create Survey
        perfect = Survey('perfect-small', gain_pattern='gaussian',
                         n_sidelobes=8)

        # Observe populations
        pop_obs[n] = SurveyPopulation(pop[n], perfect)
        pop_obs[n].name = f'sc-{n}-obs'
        pop_obs[n].rates()
        pop_obs[n].save()

else:
    for n in (35, 36, 37):
        pop_obs[n] = unpickle(f'sc-{n}-obs')


f, (ax1) = plt.subplots(1, 1)

for n in pop_obs:

    pop = pop_obs[n]
    limit = 1e-9

    s_peak = pop.frbs.s_peak
    dm_igm = pop.frbs.dm_igm

    mini = min(s_peak)
    maxi = max(s_peak)

    bins_s_peak = np.logspace(np.log10(mini), np.log10(maxi), 50)

    dm_igm = dm_igm[(s_peak > limit)]

    print(f'{len(dm_igm)} FRBs in graph of {pop.name}')

    n, bins = np.histogram(dm_igm, bins=50)
    n = n/max(n)
    bincentres = (bins[:-1] + bins[1:]) / 2
    ax1.step(bincentres, n, where='mid', label=pop.name)

# TODO Add in plot of frbcat, parkes, dm_igm
# parkes_dm_igm = Frbcat().df[('telescope' == 'parkes')]['dm_igm']

ax1.set_xlabel('$DM_{IGM}$')
ax1.set_ylabel(r'Fraction')
ax1.set_ylim([0, 1])
ax1.legend()

plt.tight_layout()
plt.savefig(f'plots/dm_igm_sc.pdf')
