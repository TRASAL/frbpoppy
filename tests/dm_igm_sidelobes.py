"""Plot how sidelobes can change the IGM DM distribution."""

import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

from convenience import plot_aa_style, rel_path

SIDELOBES = [0, 1, 8]

# Generate population with standard candles
pop = CosmicPopulation(5e5,
                       days=1,
                       name='standard',
                       dm_host_model='gaussian',
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

survey = Survey('perfect-small', gain_pattern='airy')

for sidelobe in SIDELOBES:

    survey.n_sidelobes = sidelobe

    # Observe populations
    pop_obs[sidelobe] = SurveyPopulation(pop, survey)
    pop_obs[sidelobe].name = f'obs-{sidelobe}'
    pop_obs[sidelobe].rates()
    pop_obs[sidelobe].save()

plot_aa_style()
f, (ax1) = plt.subplots(1, 1)

for p in SIDELOBES:

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
plt.savefig(rel_path(f'./plots/dm_igm_sidelobes.pdf'))
