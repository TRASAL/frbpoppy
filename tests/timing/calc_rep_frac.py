"""Calculate fraction of repeaters in detected frbs."""
from copy import deepcopy
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, merge_pop
from frbpoppy import pprint, log10normal, hist

from tests.convenience import plot_aa_style, rel_path

MAX_DAYS = 50
N_SRCS = int(1e5)
RATE = 0.1  # per day
PLOT_SNR = False


def setup_pop(n_srcs):
    pop = CosmicPopulation.simple(n_srcs, n_days=MAX_DAYS,
                                  repeaters=True)
    pop.set_dist(z_max=0.01)
    pop.set_lum(model='powerlaw', low=1e30, high=1e40, power=-1)
    pop.set_w(model='constant', value=1)
    pop.set_dm(mw=False, igm=True, host=False)
    return pop


def adapt_pop(pop, dist_type):
    if dist_type == 'rep':
        pop.set_time(model='poisson', rate=RATE)
        pop.generate()
        return pop

    elif dist_type == 'dist':
        rate_dist = log10normal(RATE, 2, N_SRCS)
        pop.set_time(model='poisson', rate=rate_dist)
        pop.generate()
        return pop

    elif dist_type == 'mix':
        # 50% one offs, 50% repeaters
        repeaters = setup_pop(N_SRCS/2)
        repeaters.set_time(model='poisson', rate=RATE)
        repeaters.name = 'repeaters'
        repeaters.generate()
        one_offs = setup_pop(N_SRCS/2)
        one_offs.set_time(model='single')
        one_offs.name = 'one-offs'
        one_offs.generate()
        return merge_pop(repeaters, one_offs, random=True)

    else:
        raise ValueError('Unknown distribution type')


r = setup_pop(N_SRCS)

# Set up surveys
survey = Survey('perfect', n_days=MAX_DAYS)
survey.mount_type = 'transit'
survey.t_obs = 60*60*24
survey.set_beam(model='perfect')
survey.snr_limit = 1e-2

chime_survey = deepcopy(survey)
chime_survey.set_beam(model='chime')
chime_survey.snr_limit = 1e-4

# Set up plot style
plot_aa_style(cols=2)
f, (ax1, ax2) = plt.subplots(1, 2)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dotted']
dist_types = ['rep', 'dist', 'mix']


for i, dist_type in enumerate(dist_types):

    r = adapt_pop(r, dist_type)
    surv_pops = [SurveyPopulation(r, s) for s in (survey, chime_survey)]

    color = colors[i]

    # Plot rates
    rate = np.array([RATE])

    if dist_type == 'mix':
        rate = np.array([RATE, 0])

    bins = np.linspace(0, 0.5, 20)
    rates, values = hist(rate, bins=bins, norm='prob')

    if dist_type == 'dist':
        rate_dist = log10normal(RATE, 2, N_SRCS)
        rates, values = hist(rate_dist, bin_type='lin', norm='prob')

    ax2.step(rates, values, where='mid', label=dist_type, color=color)

    for ii, surv_pop in enumerate(surv_pops):

        linestyle = linestyles[ii]

        if PLOT_SNR:
            plt.clf()
            snr = surv_pop.frbs.snr
            plt.step(*hist(snr, bin_type='lin', norm='prob'), where='mid',
                     color=color, linestyle=linestyle)
            # plt.xscale('log')
            plt.yscale('log')
            plt.show()
            exit()

        # See how fraction changes over time
        n_pointings = survey.n_pointings
        days = np.linspace(0, MAX_DAYS, (MAX_DAYS*n_pointings)+1)
        fracs = []
        for day in days:
            t = surv_pop.frbs.time.copy()
            time = np.where(t < day, t, np.nan)
            n_rep = ((~np.isnan(time)).sum(1) > 1).sum()
            n_one_offs = ((~np.isnan(time)).sum(1) == 1).sum()
            frac = n_rep / (n_rep + n_one_offs)
            fracs.append(frac)

        pprint(f'{dist_type} has:')
        pprint(f' - n_bursts: {surv_pop.n_bursts()}')
        pprint(f' - n_rep: {n_rep}')
        pprint(f' - n_one_offs: {n_one_offs}')

        ax1.plot(days, fracs, color=color, linestyle=linestyle)


# Further plot details
ax1.set_xlabel(r'Time (days)')
ax1.set_ylabel(r'$N_{\textrm{repeaters}}/N_{\textrm{detections}}$')
ax1.set_xlim(0, max(days))
ax1.set_yscale('log')
ax1.set_ylim(1e-2, 1e0)

ax2.set_xlabel(r'Poisson Rate (day$^{-1}$)')
ax2.set_ylabel(r'Fraction')
# ax2.set_xscale('log')
ax2.set_xlim(-0.01, 0.5)
ax2.set_ylim(0, 1)
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
ax2.set_ylim(0, 1.1)

# Legend details
plt.rcParams['legend.title_fontsize'] = 8

# Add distribution types
elements = []
for i, dist_type in enumerate(dist_types):
    line = Line2D([0], [0], color=colors[i])
    label = dist_type
    elements.append((line, label))
lines, labels = zip(*elements)
ax2.legend(lines, labels, title='Populations', prop={'size': 8},
           loc='upper right')

# Add line styles
elements = []
elements.append((Line2D([0], [0], color='gray'), 'perfect'))
elements.append((Line2D([0], [0], color='gray', linestyle='dashed'), 'chime'))
lines, labels = zip(*elements)
ax1.legend(lines, labels, title='Beam patterns', prop={'size': 8},
           loc='lower right')

# # Add legends to plot
# plt.gca().add_artist(leg_1)
# plt.gca().add_artist(leg_2)

# Save figure
plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_frac.pdf'))
plt.clf()
