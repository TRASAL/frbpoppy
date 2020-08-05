"""Compare the input versus output poisson rate.

TODO: Finish script
"""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, log10normal, Survey, SurveyPopulation
from frbpoppy import hist

from tests.convenience import plot_aa_style, rel_path

N_DAYS = 50
N_SRCS = int(1e3)
RATE = 0.1  # per day
SINGLE_INPUT_RATE = True
init_surv_time_frac = 0.1

# Set up a population
pop = CosmicPopulation.simple(N_SRCS, n_days=N_DAYS, repeaters=True)
pop.set_dist(z_max=0.01)
pop.set_lum(model='powerlaw', low=1e30, high=1e40, power=-1)
pop.set_w(model='constant', value=1)
pop.set_dm(mw=False, igm=True, host=False)

# Add a distribution of Poisson burst rates
if SINGLE_INPUT_RATE:
    rate_dist = RATE
else:
    rate_dist = log10normal(RATE, 2, N_SRCS)
pop.set_time(model='poisson', rate=rate_dist)
pop.generate()

# Survey the high fluences
survey = Survey('perfect', n_days=N_DAYS)
survey.set_beam(model='perfect')
survey.snr_limit = 1e-6
survey.t_obs = 60*60  # seconds

# Check the burst rate
surv_pop = SurveyPopulation(pop, survey)
time = surv_pop.frbs.time

# Set up plot
plot_aa_style()

if isinstance(rate_dist, np.ndarray):
    min_rate = np.log10(np.min(rate_dist[rate_dist != 0]))
    max_rate = np.log10(max(rate_dist))
else:
    min_rate = np.log10(RATE) - 1
    max_rate = np.log10(RATE) + 1
    rate_dist = np.array(rate_dist)
bins = np.logspace(min_rate, max_rate, 20)

# Plot original distribution
rates, values = hist(rate_dist, bins=bins, norm='max')
plt.step(rates, values, where='mid', label='original')

# Plot the first frac of time
mask = (time < init_surv_time_frac*N_DAYS) & ~np.isnan(time)
half_n_bursts = np.count_nonzero(mask, axis=1)
half_rate = half_n_bursts / (init_surv_time_frac*N_DAYS)
rates, values = hist(half_rate, bins=bins, norm='max')
plt.step(rates, values, where='mid', zorder=2,
         label=fr'{int(init_surv_time_frac*100)}\% of survey time')

# Plot the full rate
n_bursts = np.count_nonzero(~np.isnan(time), axis=1)
rate = n_bursts / N_DAYS
rates, values = hist(rate, bins=bins, norm='max')
plt.step(rates, values, where='mid', zorder=1,
         label=r'100\% of survey time')

# Set up plot details
plt.xlabel('Rate (per day)')
plt.ylabel('Fraction')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/rate_intrinsic_vs_observed.pdf'))
