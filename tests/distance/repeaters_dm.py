"""Plot DM/SNR distributions of repeater populations."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot
from frbpoppy import split_pop, pprint, hist

from tests.convenience import plot_aa_style, rel_path

DAYS = 4
INTERACTIVE_PLOT = False
PLOTTING_LIMIT_N_SRCS = 0
SNR = False

r = CosmicPopulation.simple(n_srcs=int(1e5), n_days=DAYS, repeaters=True)
r.set_dist(z_max=0.01)
r.set_lum(model='powerlaw', low=1e35, high=1e45, power=-1.7,
          per_source='different')
r.set_time(model='poisson', rate=3)
r.set_dm_igm(model='ioka', slope=1000, std=0)
r.set_dm(mw=False, igm=True, host=False)
r.set_w('constant', value=1)

r.generate()

# Set up survey
survey = Survey('perfect', n_days=DAYS)
survey.set_beam(model='perfect')
survey.snr_limit = 1e6

surv_pop = SurveyPopulation(r, survey)

pprint(f'{r.n_bursts()}:{surv_pop.n_bursts()}')
pprint(f'{surv_pop.n_sources()} sources detected')

if r.n_bursts() < PLOTTING_LIMIT_N_SRCS:
    pprint(f'Not sufficient FRB sources for plotting')
    exit()

# Split population into seamingly one-off and repeater populations
mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
pop_rep, pop_one = split_pop(surv_pop, mask)
pop_rep.name += ' (> 1 burst)'
pop_one.name += ' (1 burst)'

if INTERACTIVE_PLOT:
    plot(r, pop_rep, pop_one, frbcat=False, mute=False)

# Plot dm distribution
if SNR:
    plot_aa_style(cols=2)
    f, (ax1, ax2) = plt.subplots(1, 2)
else:
    plot_aa_style(cols=1)
    f, ax1 = plt.subplots(1, 1)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

pops = (r, pop_rep, pop_one)

for i, pop in enumerate(pops):

    # Distinguish populations
    if pop.name.endswith('(1 burst)'):
        label = '1 burst'
        linestyle = 'solid'
    elif pop.name.endswith('(> 1 burst)'):
        label = '$>$1 burst'
        linestyle = 'dashed'
    else:
        label = 'cosmic'
        linestyle = 'dashdot'

    pprint(f'Number of bursts in {label}: {pop.n_bursts()}')

    # Do stuff with data
    dm = pop.frbs.dm
    x, y = hist(dm)
    x *= 200  # Normalise x-axis z=0.01, z=2

    # Plot DM distributions
    ax1.step(x, y, where='mid', linestyle=linestyle, label=label,
             color=colors[i])

    # Plot fluence distributions
    snr = pop.frbs.snr

    if snr is None:
        continue

    if not SNR:
        continue

    try:
        ax2.step(*hist(snr, bin_type='log'), where='mid', linestyle=linestyle,
                 color=colors[i])
    except ValueError:
        pprint('Zero sources available to plot')
        continue

ax1.set_xlabel(r'DM$_{\textrm{ex}}$ ($\textrm{pc}\ \textrm{cm}^{-3}$)')
ax1.set_ylabel('Fraction')

if SNR:
    ax2.set_xlabel(r'SNR')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.yaxis.tick_right()
    plt.figlegend(loc='upper center', ncol=len(pops), framealpha=1)
else:
    plt.figlegend(loc='upper center', ncol=3, framealpha=1, prop={'size': 8},
                  bbox_to_anchor=(0.5, 1.07), bbox_transform=ax1.transAxes)

# Test the difference between the distributions
rep = pops[1].frbs.dm[:20]  # First 20 repeaters
one = pops[2].frbs.dm[:min([200, len(pops[2].frbs.dm)])]  # First x one-offs
print(f'KS test: {ks_2samp(rep, one)}')

plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_dm_dist.pdf'))
plt.clf()
