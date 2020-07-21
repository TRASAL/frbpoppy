"""Plot redshift versus observed repetition rate."""
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, pprint

from tests.convenience import plot_aa_style, rel_path
from tests.chime.obs_rep_frac import get_chime_rep

n_srcs = 1e5
n_days = 4
RATE = 3
CHIME_AVG = False

pop = CosmicPopulation.simple(n_srcs, n_days=n_days, repeaters=True)
pop.set_dist(z_max=2)
pop.set_lum(model='powerlaw', per_source='different', low=1e40, high=1e45,
            power=-1.7)
pop.set_w(model='constant', value=1)
pop.set_dm(mw=False, igm=True, host=False)
pop.set_dm_igm(model='ioka', slope=1000, std=0)
pop.set_time(model='poisson', rate=RATE)
pop.generate()

survey = Survey('perfect', n_days=n_days)
survey.snr_limit = 1e4

surv_pop = SurveyPopulation(pop, survey)

pprint(f'{surv_pop.n_bursts()} bursts')
if surv_pop.n_bursts() < 10:
    pprint(f'Insufficient FRB sources for plotting')
    exit()

plot_aa_style(cols=1)
f, ax1 = plt.subplots(1, 1)

# Plot
bins = np.linspace(min(surv_pop.frbs.dm), max(surv_pop.frbs.dm), 25)
n_bursts = []
for bin in zip(bins[:-1], bins[1:]):
    low, high = bin
    ii = np.where((surv_pop.frbs.dm >= low) & (surv_pop.frbs.dm < high))
    n_burst_per_src = np.sum(~np.isnan(surv_pop.frbs.time[ii]), axis=1)
    n_bursts.append(np.mean(n_burst_per_src))


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

ax1.step(bins[:-1], np.array(n_bursts) / n_days, where='mid',
         label='sim obs', color=colors[0])
ax1.set_ylabel(r'$\bar{r}_{\text{obs.~bursts, frbpoppy}}$ (day $^{-1}$)',
               color=colors[0])
ax1.tick_params(axis='y', labelcolor=colors[0])
# ax1.set_yscale('log')


df = get_chime_rep()
chime_dms = []
chime_n_bursts = []

for i, group in df.groupby('name'):
    chime_dms.append(np.mean(group.dm))
    chime_n_bursts.append(len(group))

chime_dms = np.array(chime_dms)
chime_n_bursts = np.array(chime_n_bursts)
chime_n_days = (max(group.timestamp) - min(group.timestamp)).days

ax2 = ax1.twinx()

# if CHIME_AVG:
#     n_bursts = []
#     for bin in zip(bins[:-1], bins[1:]):
#         low, high = bin
#         ii = np.where((chime_dms >= low) & (chime_dms < high))
#         n_burst_per_src = chime_n_bursts[ii]
#         n_bursts.append(np.mean(n_burst_per_src))
#
#     ax2.step(bins[:-1], np.array(n_bursts), where='mid', label='chime',
#              color=colors[1])

ax2.scatter(chime_dms, chime_n_bursts/chime_n_days, label='chime', marker='x',
            color=colors[1])


ax2.set_ylabel(r'$r_{\text{obs.~bursts, chime}}$ (day $^{-1}$)', color=colors[1])
ax2.tick_params(axis='y', labelcolor=colors[1])
# ax2.set_yscale('log')

ax1.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
# ax1.set_xscale('log')

# plt.figlegend(loc='upper center', ncol=3, framealpha=1, prop={'size': 8},
#                   bbox_to_anchor=(0.5, 1.07), bbox_transform=ax1.transAxes)

plt.tight_layout()
plt.savefig(rel_path(f'plots/dm_rep_rate.pdf'))
plt.clf()
