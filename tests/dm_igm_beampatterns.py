"""Plot how beampatterns can change the IGM DM distribution."""

import numpy as np
import matplotlib.pyplot as plt
import os

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

BEAMPATTERNS = ['perfect', 'airy', 'gaussian']

# Generate population with standard candles
pop = CosmicPopulation(5e5,
                       days=1,
                       name='standard',
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
                       w_model='uniform',
                       w_range=[1., 1.],
                       w_mu=1.,
                       w_sigma=0.,
                       si_mu=0.,
                       si_sigma=0.,
                       z_max=2.5)

pop_obs = {}

survey = Survey('perfect-small', n_sidelobes=0)

for pattern in BEAMPATTERNS:

    survey.gain_pattern = pattern

    # Observe populations
    pop_obs[pattern] = SurveyPopulation(pop, survey)
    pop_obs[pattern].name = f'obs-{pattern}'
    pop_obs[pattern].rates()
    pop_obs[pattern].save()

# Change working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
plt.style.use('./aa.mplstyle')
f, (ax1) = plt.subplots(1, 1)

for p in BEAMPATTERNS:

    pop = pop_obs[p]
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
    ax1.step(bincentres, n, where='mid', label=p)

ax1.set_xlabel(r'DM$_{\text{IGM}}$')
ax1.set_ylabel(r'Fraction')
ax1.set_ylim([0, 1])
ax1.legend()
plt.tight_layout()
plt.savefig('./plots/dm_igm_beampatterns.pdf')
