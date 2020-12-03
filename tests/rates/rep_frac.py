"""Calculate repeater fraction over time for various rate distributions."""
from copy import deepcopy
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, merge_pop
from frbpoppy import log10normal, hist, clustered

from tests.convenience import plot_aa_style, rel_path

MAX_DAYS = 100
N_SRCS = int(1e5)
RATE = 0.1  # per day
PLOT_SNR = False


def setup_pop(n_srcs):
    pop = CosmicPopulation.simple(n_srcs, n_days=MAX_DAYS,
                                  repeaters=True)
    pop.set_dist(z_max=0.01)
    pop.set_lum(model='powerlaw', per_source='different', low=1e35, high=1e40,
                power=-1.7)
    pop.set_w(model='constant', value=1)
    pop.set_dm(mw=False, igm=True, host=False)
    pop.set_dm_igm(model='ioka', slope=1000, std=0)
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

    elif dist_type == 'weibull':
        pop.set_time(model='clustered', r=RATE, k=0.34)
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


def get_frac_over_time(surv_pop, survey):
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
    print(f'# one-offs: {n_one_offs}')
    print(f'# repeaters: {n_rep}')
    return days, fracs


r_pop = setup_pop(N_SRCS)

# Set up surveys
survey = Survey('perfect', n_days=MAX_DAYS)
survey.mount_type = 'transit'
survey.t_obs = 60*60*24
survey.set_beam(model='perfect')
survey.snr_limit = 10000

chime_survey = deepcopy(survey)
chime_survey.set_beam(model='chime-frb')
chime_survey.snr_limit = 1

# Set up plot style
plot_aa_style(cols=2)
f, (ax1, ax2) = plt.subplots(1, 2)

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
linestyles = ['solid', 'dotted']
dist_types = ['rep', 'dist', 'mix']

# Further plot details
ax1.set_xlabel(r'Poisson Rate (day$^{-1}$)')
ax1.set_ylabel(r'Fraction')
ax1.set_xlim(-0.01, 0.5)
# ax1.set_ylim(0, 1)
# ax1.set_yscale('log')

ax2.set_xlabel(r'Time (days)')
ax2.set_ylabel(r'$f_{\text{rep}}$')
ax2.set_xlim(0, MAX_DAYS)
ax2.set_yscale('log')
ax2.set_ylim(1e-1, 1e0)
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()

# Legend details
plt.rcParams['legend.title_fontsize'] = 8

ax1_elements = []
ax2_elements = [(Line2D([0], [0], color='gray'), 'perfect')]

for i, dist_type in enumerate(dist_types):
    print(f'Distribution type: {dist_type}')
    color = colors[i]

    # Plot rates
    rate = np.array([RATE])

    if dist_type == 'mix':
        rate = np.array([RATE, 0])

    bins = np.linspace(0, 0.5, 40)
    rates, values = hist(rate, bins=bins, norm='prob')

    if dist_type == 'dist':
        rate_dist = log10normal(RATE, 2, N_SRCS)
        rates, values = hist(rate_dist, bins=bins, norm='prob')

    if dist_type == 'weibull':
        rate_dist = clustered(n_srcs=N_SRCS, n_days=MAX_DAYS, z=0)
        rates, values = hist(rate_dist, bins=bins, norm='prob')

    ax1.step(rates, values, where='mid', label=dist_type, color=color)

    r = adapt_pop(r_pop, dist_type)
    surv_pop = SurveyPopulation(r, survey)
    linestyle = linestyles[0]

    # See how fraction changes over time
    days, fracs = get_frac_over_time(surv_pop, survey)

    ax2.plot(days, fracs, color=color, linestyle=linestyle)

    # Add distribution types
    line = Line2D([0], [0], color=colors[i])
    label = dist_type
    ax1_elements.append((line, label))
    lines, labels = zip(*ax1_elements)
    ax1.legend(lines, labels, title='Populations', prop={'size': 8},
               loc='upper right')

    # Add line styles
    lines, labels = zip(*ax2_elements)
    ax2.legend(lines, labels, title='Beam patterns', prop={'size': 8},
               loc='lower right')

    # Save figure
    plt.tight_layout()
    plt.savefig(rel_path(f'plots/rep_frac_{i}.pdf'))


# Plot surveys with chime beam components
ax2_elements.append((Line2D([0], [0], color='gray', linestyle='dashed'),
                     'chime-frb'))
# Add line styles
lines, labels = zip(*ax2_elements)
ax2.legend(lines, labels, title='Beam patterns', prop={'size': 8},
           loc='lower right')

for i, dist_type in enumerate(dist_types):
    r = adapt_pop(r_pop, dist_type)
    surv_pop = SurveyPopulation(r, chime_survey)
    linestyle = linestyles[1]
    days, fracs = get_frac_over_time(surv_pop, survey)
    ax2.plot(days, fracs, color=colors[i], linestyle=linestyle)

# Save figure
plt.tight_layout()
plt.savefig(rel_path('plots/rep_frac.pdf'))
plt.clf()
