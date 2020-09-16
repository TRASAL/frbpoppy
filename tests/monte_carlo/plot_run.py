"""Plot distributions from a single Monte Carlo run."""
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from matplotlib.offsetbox import AnchoredText
import numpy as np
import os
from functools import reduce

from frbpoppy import poisson_interval, unpickle, hist

from frbcat import TNS

from tests.rates.alpha_real import EXPECTED
from tests.convenience import plot_aa_style, rel_path

from run_mc import Run

r = Run(survey_name='htru',
        li=-1.5,
        si=-1.1,
        alpha=-2.0)

# Find closest population match
path = rel_path('../data/populations/mc')
options = os.listdir(path)
alphas = np.zeros(len(options))
lis = np.zeros(len(options))
sis = np.zeros(len(options))

for i, line in enumerate(options):
    if line.split('_')[-1].split('.')[0] != r.survey_name:
        continue
    values = map(line.split('_').__getitem__, [2, 4, 6])
    alphas[i], lis[i], sis[i] = [float(v) for v in values]

alpha_i = np.where(np.abs(alphas-r.alpha) == np.abs(alphas-r.alpha).min())
li_i = np.where(np.abs(lis-r.li) == np.abs(lis-r.li).min())
si_i = np.where(np.abs(sis-r.si) == np.abs(sis-r.si).min())
i = reduce(np.intersect1d, (alpha_i, li_i, si_i))[0]
r.alpha, r.li, r.si = map(options[i].split('_').__getitem__, [2, 4, 6])

f = f'mc/complex_alpha_{r.alpha}_lum_{r.li}_si_{r.si}_{r.survey_name}'
surv_pop = unpickle(f)
df = TNS(repeaters=False).df
import IPython; IPython.embed()

plot_aa_style()
plt.rcParams["figure.figsize"] = (5.75373, 5.75373)
fig, axes = plt.subplots(2, 2)

# Rates
ax = axes[0, 0]
# Frbcat
n_frbs, n_days = EXPECTED[r.survey_name]
errs = [e/n_days for e in poisson_interval(n_frbs)]
ax.errorbar(n_frbs/n_days, 0.6, xerr=np.array([errs]).T, fmt='x')
# Frbpoppy
errs = [e/surv_pop.n_days for e in poisson_interval(surv_pop.n_sources())]
ax.errorbar(surv_pop.n_sources()/surv_pop.n_days, 0.4, xerr=np.array([errs]).T,
            fmt='x')

# Other axes properties
ax.set_xscale('log')
ax.set_xlabel(r'Rate (day$^{-1}$)')
ax.set_yticks([0.6, 0.4])
ax.set_ylim([0, 1])
ax.set_yticklabels(['frbcat', 'frbpoppy'])

# Plot information
axes[0, 1].set_axis_off()
text = f'{r.survey_name} \n li={r.li} \n si={r.si} \n alpha={r.alpha}'
axes[0, 1].text(.5, .7, text,
                horizontalalignment='center',
                verticalalignment='center',
                transform=axes[0, 1].transAxes,
                linespacing=1.6)


# DM distributions
ax = axes[1, 0]
ax.step(*hist(surv_pop.frbs.dm), where='mid')
mask = df.telescope.str.lower() == r.survey_name
ax.step(*hist(df[mask].dm), where='mid')

# Compare distributions
ks = ks_2samp(surv_pop.frbs.dm, df[mask].dm)
text = fr'$p={round(ks[1], 2)}$'
if ks[1] < 0.01:
    text = fr'$p={ks[1]:.2e}$'
anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
ax.add_artist(anchored_text)
ax.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
ax.set_ylabel('Fraction')

# SNR distribution
ax = axes[1, 1]
ax.step(*hist(surv_pop.frbs.snr, bin_type='log'), where='mid')
ax.step(*hist(df[mask].snr, bin_type='log'), where='mid')

# Compare distributions
ks = ks_2samp(surv_pop.frbs.snr, df[mask].snr)
text = fr'$p={round(ks[1], 2)}$'
if ks[1] < 0.01:
    text = fr'$p={ks[1]:.2e}$'
anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
ax.add_artist(anchored_text)
ax.set_xlabel(r'SNR')
ax.set_xscale('log')
ax.set_ylabel('Fraction')

# Final settings
plt.tight_layout()
p = f'./plots/run.pdf'
plt.savefig(rel_path(p))
