"""Plot DM/SNR distributions of repeater populations."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot
from frbpoppy import split_pop, pprint

from convenience import hist, plot_aa_style, rel_path

DAYS = 1
INTERACTIVE_PLOT = True
PLOTTING_LIMIT_N_SRCS = 0
SNR = True

r = CosmicPopulation.simple(n_srcs=int(1e5), n_days=DAYS, repeaters=True)
r.set_dist(z_max=0.01)
r.set_lum(model='powerlaw', low=1e35, high=1e45, power=2,
          per_source='different')
r.set_time(model='poisson', lam=3)
r.set_dm_igm(model='ioka', slope=1000, sigma=0)
r.set_dm(mw=False, igm=True, host=False)
r.set_w('constant', value=1)
# r.set_w(model='lognormal', mu=np.log(1), sigma=np.log(2))

r.generate()

# Set up survey
survey = Survey('perfect', n_days=DAYS)
survey.set_beam(model='perfect')
# survey.t_samp = 1
survey.snr_limit = 1e5

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
    plot(r, pop_rep, pop_one, frbcat=False)

# Plot dm distribution
if SNR:
    plot_aa_style(cols=2)
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
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

    ax2.step(*hist(snr, bin_type='log'), where='mid', linestyle=linestyle,
             color=colors[i])

ax1.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
ax1.set_ylabel('Fraction')
# ax1.set_xlim([0, 10])

if SNR:
    ax2.set_xlabel(r'SNR')
    plt.xscale('log')
    plt.yscale('log')
    plt.figlegend(loc='upper center', ncol=len(pops), framealpha=1)
else:
    plt.legend()

plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_dm_dist.pdf'))
plt.clf()
