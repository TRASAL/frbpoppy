# -*- coding: future_fstrings -*-
"""Script to reproduce Fig. 6 of Connor et al. 2017.

Plot N(>S) over log S (S being the flux density) for various spectral indices.
"""
import copy
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

CREATE = False
OBSERVE = False
SIS = (-2, 0, 2)  # Spectral indices

pop = {}

if CREATE:
    days = 14
    n_per_day = 5000
    n_tot = n_per_day*days

    for si in SIS:

        if si == min(SIS):

            pop[si] = CosmicPopulation(n_tot,
                                       days=days,
                                       name=f'si-{si}',
                                       dm_host_model='normal',
                                       dm_host_mu=0,
                                       dm_host_sigma=0,
                                       dm_igm_index=1200,
                                       dm_igm_sigma=0,
                                       dm_mw_model='zero',
                                       emission_range=[10e6, 10e9],
                                       lum_range=[1e40, 1e40],
                                       lum_index=0,
                                       n_model='vol_co',
                                       w_model='uniform',
                                       w_range=[1., 1.],
                                       w_mu=1.,
                                       w_sigma=0.,
                                       repeat=0.,
                                       si_mu=si,
                                       si_sigma=0.,
                                       z_max=2.5)
            pop[si].save()

        else:
            pop[si] = copy.deepcopy(pop[min(SIS)])
            pop[si].frbs.si = np.random.normal(si, 0, n_tot)
            pop[si].name = f'si-{si}'
            pop[si].save()

pop_obs = {}

if OBSERVE or CREATE:

    for si in SIS:

        if not CREATE:
            pop[si] = unpickle(f'si-{si}')

        # Create Survey
        perfect = Survey('perfect', gain_pattern='perfect')

        # Observe populations
        pop_obs[si] = SurveyPopulation(pop[si], perfect)
        pop_obs[si].name = f'si-{si}-obs'
        pop_obs[si].rates()
        pop_obs[si].save()

else:
    for si in SIS:
        pop_obs[si] = unpickle(f'si-{si}-obs')


# Plot log N and alpha versus log S
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

min_s = 1e99
max_s = -1e99

for si in SIS:

    pop = pop_obs[si]

    s_peak = pop.frbs.s_peak

    # Bin up
    number, bins = np.histogram(np.log10(s_peak), bins=500)  # N(S)
    n_gt_s = np.cumsum(number[::-1])[::-1]  # N(>S) from N(S)
    x = bins[:-1]  # log(S)
    y = np.log10(n_gt_s)  # log(N(>S))

    ax1.step(x, y, where='pre', label=fr"$\gamma$ of {si}")

    # Plot alpha
    # Calculate derivative
    der = np.diff(y) / np.diff(x)
    bin_centres = (x[:-1] + x[1:]) / 2

    # Smooth function
    derhat = savgol_filter(der, 51, 3)
    ax2.step(bin_centres, derhat, where='mid')

    if min(bin_centres) <= min_s:
        min_s = min(bin_centres)
    if max(bin_centres) >= max_s:
        max_s = max(bin_centres)

# Add a -3/2 slope
x = np.linspace(min_s, max_s, 1000)
y = -1.5*x
y -= min(y)
y += min(np.log10(n_gt_s))
x = x[y <= max(np.log10(n_gt_s))]
y = y[y <= max(np.log10(n_gt_s))]
ax1.step(x, y, where='mid', color='grey', alpha=0.5)

# Plot alpha over log S

# Plot a Euclidean line
x = np.linspace(min_s, max_s, 1000)
y = np.ones_like(x) * -1.5
ax2.step(x, y, where='mid', color='grey', alpha=0.5)


ax1.set_ylabel(r'log N(>S$_{\text{peak}}$)')
ax1.legend()
ax2.set_xlabel(r'log S$_{\text{peak}}$')
ax2.set_ylabel(r'$\alpha$')
ax2.set_ylim(ax2.get_ylim()[::-1])

plt.tight_layout()
plt.savefig(f'plots/logn_logs_si.pdf')
