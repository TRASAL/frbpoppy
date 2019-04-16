"""Plot how beampatterns can change the IGM DM distribution."""

import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import unpickle

CREATE = False
OBSERVE = True
SIDELOBES = [0, 1, 8]

if CREATE:
    days = 30
    n_per_day = 5000

    # Generate population with standard candles
    pop = CosmicPopulation(n_per_day*days,
                           days=days,
                           name='simple',
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
                           si_mu=0.,
                           si_sigma=0.,
                           z_max=2.5)
    pop.save()

pop_obs = {}

if OBSERVE:

    if not CREATE:
        pop = unpickle(f'simple')

    for n in SIDELOBES:

        # Create Survey
        survey = Survey('perfect-small', gain_pattern='airy', n_sidelobes=n)

        # Observe populations
        pop_obs[n] = SurveyPopulation(pop, survey)
        pop_obs[n].name = f'obs-airy-{n}'
        pop_obs[n].rates()
        pop_obs[n].save()

else:
    for n in SIDELOBES:
        pop_obs[n] = unpickle(f'obs-airy-{n}')

f, (ax1) = plt.subplots(1, 1)

for sidelobe in SIDELOBES:

    pop = pop_obs[sidelobe]
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
    ax1.step(bincentres, n, where='mid', label=f'airy-{sidelobe}')

ax1.set_xlabel('$DM_{IGM}$')
ax1.set_ylabel(r'Fraction')
ax1.set_ylim([0, 1])
ax1.legend()

plt.tight_layout()
plt.savefig(f'plots/dm_igm_sidelobes.pdf')
