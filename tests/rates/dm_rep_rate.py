"""Plot redshift versus observed repetition rate."""
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, pprint, split_pop

from tests.convenience import plot_aa_style, rel_path
from tests.chime.obs_rep_frac import get_chime_rep

n_srcs = 1e5
n_days = 4
RATE = 3
CHIME_HIST = False
SCATTER = False

pop = CosmicPopulation.simple(n_srcs, n_days=n_days, repeaters=True)
pop.set_dist(z_max=2.0)
pop.set_w(model='constant', value=1)
pop.set_dm(mw=False, igm=True, host=False)
pop.set_dm_igm(model='ioka', slope=1000, std=0)
pop.set_time(model='poisson', rate=RATE)
pop.set_time(model='regular', rate=RATE)
pop.set_w('constant', value=1)
pop.generate()

survey = Survey('perfect', n_days=n_days)
survey.snr_limit = 1e6

plot_aa_style(cols=1)
f, ax1 = plt.subplots(1, 1)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

lum_funcs = ['neg. pl.', 'flat pl.', 'std. candle']
for i, lum_func in enumerate(lum_funcs):

    if lum_func == 'neg. pl.':
        pop.set_lum(model='powerlaw', per_source='different', low=1e40,
                    high=1e45, power=-1.7)
    elif lum_func == 'flat pl.':
        pop.set_lum(model='powerlaw', per_source='different', low=1e40,
                    high=1e45, power=0)
    elif lum_func == 'gauss':
        pop.set_lum(model='gauss', per_source='different', mean=1e42,
                    std=1e10)
    elif lum_func == 'std. candle':
        pop.set_lum(model='constant', value=1e42)

    pop.generate()
    surv_pop = SurveyPopulation(pop, survey)
    # Split population into seamingly one-off and repeater populations
    mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
    pop_ngt1, pop_nle1 = split_pop(surv_pop, mask)

    pprint(f'{pop_ngt1.n_bursts()} repeater bursts')
    if pop_ngt1.n_bursts() < 10:
        pprint(f'Insufficient FRB sources for plotting')

    # Plot
    if SCATTER:
        def sample(n):
            return np.random.randint(low=0, high=pop_ngt1.n_sources(), size=n)

        def accept(ii):
            return surv_pop.frbs.dm[ii] < 1000

        ii = sample(25)
        mask = accept(ii)
        reject, = np.where(~mask)
        while reject.size > 0:
            fill = sample(reject.size)
            mask = accept(fill)
            ii[reject[mask]] = fill[mask]
            reject = reject[~mask]

        n_bursts = np.sum(~np.isnan(surv_pop.frbs.time[ii]), axis=1)
        ax1.scatter(pop_ngt1.frbs.dm[ii], n_bursts,
                    label=lum_func, color=colors[i])
    else:
        bins = np.linspace(min(pop_ngt1.frbs.dm), max(pop_ngt1.frbs.dm), 25)
        n_bursts = []
        for bin in zip(bins[:-1], bins[1:]):
            low, high = bin
            ii = np.where((pop_ngt1.frbs.dm >= low) & (pop_ngt1.frbs.dm < high))
            n_burst_per_src = np.sum(~np.isnan(pop_ngt1.frbs.time[ii]), axis=1)
            n_bursts.append(np.mean(n_burst_per_src))

        ax1.step(bins[:-1], np.array(n_bursts), where='mid',
                 label=lum_func, color=colors[i])


ax1.set_ylabel(r'$\bar{n}_{\text{obs.~bursts, frbpoppy}}$')

plt.legend()

df = get_chime_rep()
chime_dms = []
chime_n_bursts = []

for i, group in df.groupby('name'):
    chime_dms.append(np.mean(group.dm) - np.mean(group.ne2001))
    chime_n_bursts.append(len(group))

chime_dms = np.array(chime_dms)
chime_n_bursts = np.array(chime_n_bursts)
chime_n_days = (max(group.timestamp) - min(group.timestamp)).days

ax2 = ax1.twinx()

if CHIME_HIST:
    n_bursts = []
    for bin in zip(bins[:-1], bins[1:]):
        low, high = bin
        ii = np.where((chime_dms >= low) & (chime_dms < high))
        n_burst_per_src = chime_n_bursts[ii]
        n_bursts.append(np.mean(n_burst_per_src))

    ax2.step(bins[:-1], np.array(n_bursts), where='mid', label='chime',
             color=colors[len(lum_funcs)])
else:
    ax2.scatter(chime_dms, chime_n_bursts, label='chime', marker='x',
                color=colors[len(lum_funcs)])


ax2.set_ylabel(r'$n_{\text{obs.~bursts, chime}}$', color=colors[len(lum_funcs)])
ax2.tick_params(axis='y', labelcolor=colors[len(lum_funcs)])

ax1.set_xlabel(r'DM$_{\textrm{ex}}$ ($\textrm{pc}\ \textrm{cm}^{-3}$)')

# # Add simulated chime burstiness
# surv_pop = unpickle('cosmic_chime')
# # Split population into seamingly one-off and repeater populations
# mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
# pop_ngt1, pop_nle1 = split_pop(surv_pop, mask)
#
# n_bursts = np.sum(~np.isnan(pop_ngt1.frbs.time), axis=1)
# ax2.scatter(pop_ngt1.frbs.dm, n_bursts, label='frbpoppy', color=colors[len(lum_funcs)+1])

plt.tight_layout()
plt.savefig(rel_path(f'plots/dm_rep_rate.pdf'))
plt.clf()
