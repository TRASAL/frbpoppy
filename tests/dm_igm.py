"""Replicate figure 3 of Niino (2018).

This figure compares input luminosity vs apparent luminosity for a variety of
cases.
"""
import copy
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Adapt, Survey, SurveyPopulation
from frbpoppy import unpickle, plot, Frbcat

CREATE = False
OBSERVE = False

if CREATE:
    days = 56
    n_per_day = 5000

    # Generate population with standard candles
    pop_sc = CosmicPopulation(n_per_day*days,
                              days=days,
                              name='sc-sfr',
                              dm_host_model='normal',
                              dm_host_mu=0,
                              dm_host_sigma=0,
                              dm_igm_index=1200,
                              dm_igm_sigma=0,
                              dm_mw_model='zero',
                              emission_range=[10e6, 10e9],
                              lum_range=[1e36, 1e36],
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
    pop_sc.save()

    # Generate population with a powerlaw
    pop_pl = copy.deepcopy(pop_sc)
    pop_pl = Adapt(pop_pl).lum_bol(1e36, 1e50, -2)
    pop_pl.name = 'pl-sfr'
    pop_pl.save()

if OBSERVE:
    if not CREATE:
        pop_sc = unpickle('sc-sfr')
        pop_pl = unpickle('pl-sfr')

    # Create Survey
    perfect = Survey('PERFECT-SMALL', gain_pattern='airy', sidelobes=0)

    # Observe populations
    pop_sc_obs = SurveyPopulation(pop_sc, perfect)
    pop_sc_obs.name = 'sc-sfr-obs'
    pop_sc_obs.rates()
    pop_sc_obs.save()
    pop_pl_obs = SurveyPopulation(pop_pl, perfect)
    pop_pl_obs.name = 'pl-sfr-obs'
    pop_pl_obs.rates()
    pop_pl_obs.save()

else:
    pop_sc_obs = unpickle('sc-sfr-obs')
    pop_pl_obs = unpickle('pl-sfr-obs')


limits = np.logspace(np.log10(1e-21), np.log10(1e0), 40)
k = 0

for limit in limits:

    f, (ax1, ax2) = plt.subplots(1, 2)

    for pop in (pop_sc_obs, pop_pl_obs):

        s_peak = np.array(pop.get('s_peak'))
        dm_igm = np.array(pop.get('dm_igm'))
        lum = np.array(pop.get('lum_bol'))

        mini = min(s_peak)
        maxi = max(s_peak)

        bins_s_peak = np.logspace(np.log10(mini), np.log10(maxi), 50)

        dm_igm = dm_igm[(s_peak > limit)]

        print(f'{len(dm_igm)} FRBs in graph of {pop.name}')

        n, bins = np.histogram(s_peak[(s_peak >= limit)], bins=bins_s_peak)
        bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
        ax1.step(bincentres, n, where='mid', label=pop.name)

        n, bins = np.histogram(dm_igm, bins=50)
        n = n/max(n)
        bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
        ax2.step(bincentres, n, where='mid', label=pop.name)

    # TODO Add in plot of frbcat, parkes, dm_igm
    # parkes_dm_igm = Frbcat().df[('telescope' == 'parkes')]['dm_igm']

    ax1.axvline(x=limit, color='r')

    ax1.set_xscale('log')
    ax1.set_xlabel('$S_{peak}$')
    ax1.set_ylabel(r'Number')
    ax1.set_ylim([0, 21000])

    ax2.set_xlabel('$DM_{IGM}$')
    ax2.set_ylabel(r'Fraction')
    ax2.set_ylim([0, 1])
    ax2.legend()

    plt.tight_layout()
    plt.savefig(f'plots/dm_igm_{k}.png')
    plt.clf()
    k += 1
