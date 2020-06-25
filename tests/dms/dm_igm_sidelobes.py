"""Plot how sidelobes can change the IGM DM distribution."""

import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

from convenience import plot_aa_style, rel_path

SIDELOBES = [0, 1, 8]

# Generate population with standard candles
pop = CosmicPopulation(n_srcs=5e5, n_days=1, name='standard')
pop.set_dist(model='sfr', z_max=2.5, H_0=67.74, W_m=0.3089, W_v=0.6911)
pop.set_dm_igm(model='ioka', slope=1200, sigma=0)
pop.set_dm(mw=False, host=False, igm=True)
pop.set_emission_range(low=10e6, high=10e9)
pop.set_lum(model='constant', value=1e36)
pop.set_w(model='constant', value=1)
pop.set_si(model='constant', value=0)
pop.generate()

pop_obs = {}

survey = Survey('perfect-small')
survey.beam_size_fwhm = 10.
survey.set_beam(model='airy')

for sidelobe in SIDELOBES:

    survey.set_beam(model='airy', n_sidelobes=sidelobe)

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
