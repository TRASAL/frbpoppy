"""Plot DM/SNR distributions of repeater populations."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, plot
from frbpoppy import split_pop, pprint

from convenience import hist, plot_aa_style, rel_path

DAYS = 1
INTERACTIVE_PLOT = False
PLOTTING_LIMIT_N_FRBS = 0
SNR = False

r = RepeaterPopulation.simple(int(1e6))
r.lum_min = 1e40
r.lum_max = 1e45
r.lum_pow = 0
r.lum_rep_model = 'independent'
r.z_max = 0.01
r.times_rep_model = 'poisson'
r.n_days = DAYS

# Set DM distributions
r.dm_host_model = 'gaussian'
r.dm_host_mu = 0
r.dm_host_sigma = 0
r.dm_igm_index = 1000
r.dm_igm_sigma = 0
r.dm_mw_model = 'zero'

r.generate()

survey = Survey('perfect', strategy='regular', n_days=DAYS)
survey.beam_pattern = 'perfect'
survey.snr_limit = 1e16

surv_pop = SurveyPopulation(r, survey)

# Check whether sufficient frbs
n_frbs = len(surv_pop.frbs.index)


def n_bursts(pop):
    """Count the number of bursts."""
    return np.count_nonzero(~np.isnan(pop.frbs.time))


pprint(f'{n_bursts(r)}:{n_bursts(surv_pop)}')
pprint(f'{n_frbs} sources detected')
if n_frbs < PLOTTING_LIMIT_N_FRBS:
    pprint(f'Not sufficient FRBs for plotting')
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

    dm = pop.frbs.dm
    pprint(f'Number of bursts in {label}: {n_bursts(pop)}')
    x, y = hist(dm)
    # x *= 190
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
plt.savefig(rel_path(f'plots/rep_dist.pdf'))
plt.clf()
