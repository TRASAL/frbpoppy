"""Script to reproduce Fig. 6 of Connor et al. 2017.

Plot N(>S) over log S (S being the flux density) for various spectral indices.
"""
import copy
import numpy as np
import math
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle
from frbpoppy import Adapt

CREATE = False
OBSERVE = False
SIS = (-2, 0, 2)  # Spectral indices

pop = {}

if CREATE:
    days = 14
    n_per_day = 5000

    for si in SIS:

        if si == min(SIS):

            pop[si] = CosmicPopulation(n_per_day*days,
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
                                       pulse_model='uniform',
                                       pulse_range=[1., 1.],
                                       pulse_mu=1.,
                                       pulse_sigma=0.,
                                       repeat=0.,
                                       si_mu=si,
                                       si_sigma=0.,
                                       z_max=2.5)
            pop[si].save()

        else:
            pop[si] = copy.deepcopy(pop[min(SIS)])
            pop[si] = Adapt(pop[si]).si(si, 0.)
            pop[si].name = f'si-{si}'
            pop[si].save()

pop_obs = {}

if OBSERVE or CREATE:

    for si in SIS:

        if not CREATE:
            pop[si] = unpickle(f'si-{si}')

        # Create Survey
        perfect = Survey('PERFECT', gain_pattern='perfect')

        # Observe populations
        pop_obs[si] = SurveyPopulation(pop[si], perfect)
        pop_obs[si].name = f'si-{si}-obs'
        pop_obs[si].rates()
        pop_obs[si].save()

else:
    for si in SIS:
        pop_obs[si] = unpickle(f'si-{si}-obs')


def calc_logn_logs(parms):
    """Unbiased JP method."""
    f_0 = min(parms)
    n = len(parms)
    alpha = -1/((1/n)*sum([math.log(f/f_0) for f in parms]))
    alpha *= (n-1)/n  # Removing bias in alpha
    alpha_err = n*alpha/((n-1)*(n-2)**0.5)
    norm = n / (f_0**alpha)  # Normalisation at lowest parameter
    return alpha, alpha_err, norm


# Plot log N and alpha versus log S
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

min_s = 1e99
max_s = -1e99

for si in SIS:

    pop = pop_obs[si]

    s_peak = np.array(pop.get('s_peak'))

    # Bin up
    min_f = min(np.log10(s_peak))
    max_f = max(np.log10(s_peak))
    hist, bins = np.histogram(np.log10(s_peak), bins=100)
    cum_hist = np.cumsum(hist[::-1])[::-1]

    bin_centres = (bins[1:] + bins[:-1]) / 2

    if min(bin_centres) <= min_s:
        min_s = min(bin_centres)
    if max(bin_centres) >= max_s:
        max_s = max(bin_centres)

    ax1.step(bin_centres, np.log10(cum_hist), where='mid',
             label=fr"$\gamma$ of {si}")

    # Plot alpha
    n_points = 100
    frac = 0.5
    frac_width = frac*(bins[-1] - bins[0])
    s_peak.sort()
    points = np.linspace(bins[0], bins[-1], n_points)
    bins = []
    for point in points:
        low = point - 0.5*frac_width
        high = point + 0.5*frac_width
        bins.append((low, high))

    alphas = []
    alpha_errs = []
    s_peaks = []

    for b in bins:
        low, high = 10**b[0], 10**b[1]
        section = s_peak[s_peak >= low]
        section = section[section <= high]
        log_s_peak = (b[0] + b[1]) / 2.
        if len(section) <= 2:
            continue
        alpha, alpha_err, norm = calc_logn_logs(section)
        alphas.append(alpha)
        alpha_errs.append(alpha_err)
        s_peaks.append(log_s_peak)

    ax2.step(s_peaks, alphas, where='mid')


ax1.set_ylabel(r'log N(>S)')
ax1.legend()

# Add a -3/2 slope
x = np.linspace(min_s, max_s, 1000)
y = -1.5*x
y -= min(y)
y += min(np.log10(cum_hist))
x = x[y <= max(np.log10(cum_hist))]
y = y[y <= max(np.log10(cum_hist))]
ax1.step(x, y, where='mid', color='grey', alpha=0.5)

# Plot alpha over log S

# Plot a Euclidean line
x = np.linspace(min_s, max_s, 1000)
y = np.ones_like(x) * -1.5
ax2.step(x, y, where='mid', color='grey', alpha=0.5)
ax2.set_xlabel(r'log S$_{\text{peak}}$')
ax2.set_ylabel(r'$\alpha$')
ax2.set_ylim(ax2.get_ylim()[::-1])

plt.tight_layout()
plt.savefig(f'plots/si_logn_logs.pdf')
