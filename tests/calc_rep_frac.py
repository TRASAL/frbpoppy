"""Calculate fraction of repeaters in detected frbs."""
from bisect import bisect
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, merge_pop
from frbpoppy import pprint

from convenience import plot_aa_style, rel_path, hist

MAX_DAYS = 100
N_SRCS = 1e4
LAM = 0.4  # per day
PLOT_SNR = False

# Define standard setup
def set_up_pop(n_srcs):
    r = CosmicPopulation.simple(n_srcs, n_days=MAX_DAYS, repeaters=True)
    r.set_dist(z_max=0.01)
    r.set_lum(model='powerlaw', low=1e30, high=1e40, power=0)
    r.set_w(model='constant', value=1)
    r.set_dm(mw=False, igm=True, host=False)
    return r

r = set_up_pop(N_SRCS)

# Set up survey
survey = Survey('perfect', n_days=MAX_DAYS)
survey.mount_type = 'transit'
survey.t_obs = 60*60*24
survey.set_beam(model='chime')
survey.snr_limit = 1e-5
# import IPython; IPython.embed()

# Set up plot style
plot_aa_style(cols=2)
f, (ax1, ax2) = plt.subplots(1, 2)

for s in ['rep', 'dist', 'mix']:  #, 'mix', 'dist'):

    if s == 'rep':
        r.set_time(model='regular', lam=LAM)
        r.generate()
    elif s == 'dist':
        r.set_time(model='poisson', lam=LAM)
        r.generate()
    elif s == 'mix':
        # 50% one offs, 50% repeaters
        repeaters = set_up_pop(N_SRCS/2)
        repeaters.set_time(model='regular', lam=LAM)
        repeaters.name = 'repeaters'
        repeaters.generate()
        one_offs = set_up_pop(N_SRCS/2)
        one_offs.set_time(model='single')
        one_offs.name = 'one-offs'
        one_offs.generate()
        r = merge_pop(repeaters, one_offs, random=True)
        del one_offs
        del repeaters

    surv_pop = SurveyPopulation(r, survey)

    if PLOT_SNR:
        plt.clf()
        snr = surv_pop.frbs.snr
        plt.step(*hist(snr, bin_type='log', norm='prob'), where='mid')
        plt.xscale('log')
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

    pprint(f'{s} has:')
    pprint(f' - n_bursts: {surv_pop.n_bursts()}')
    pprint(f' - n_rep: {n_rep}')
    pprint(f' - n_one_offs: {n_one_offs}')

    ax1.plot(days, fracs)

    # Plot rates
    dt = 1/np.diff(r.frbs.time/(1+r.frbs.z)[:, np.newaxis]).flatten()
    bins = np.logspace(np.log10(LAM)-1, np.log10(LAM)+2, 20)
    bins, values = hist(dt, bins=bins, norm='prob')
    if s == 'mix':   # One-offs have no time difference
        values /= 2
    ax2.step(bins, values, where='mid', label=s)

# Further plot details
ax1.set_xlabel(r'Time (days)')
ax1.set_ylabel(r'$N_{\textrm{repeaters}}/N_{\textrm{detections}}$')
ax1.set_xlim(0, max(days))
ax1.set_ylim(0, 1.1)

ax2.set_xlabel(r'Rate (/day)')
ax2.set_ylabel(r'Fraction')
ax2.set_xscale('log')
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
ax2.set_ylim(0, 1.1)

# Legend details
plt.figlegend(loc='upper center', ncol=3, framealpha=1)

plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_frac.pdf'))
plt.clf()
