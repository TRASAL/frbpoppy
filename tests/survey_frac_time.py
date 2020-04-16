"""Input versus output poisson rate."""
from copy import deepcopy

from frbpoppy import CosmicPopulation, log10normal, Survey, SurveyPopulation
from frbpoppy import plot, powerlaw


N_DAYS = 50
N_SRCS = int(1e3)
init_surv_time_frac = 0.05

# Set up a population
pop = CosmicPopulation.simple(N_SRCS, n_days=N_DAYS, repeaters=True)
pop.name = 'cosmic'
pop.set_dist(z_max=0.01)
pop.set_dm(mw=False, igm=True, host=False)

# Add a distribution of luminosities
# lum_dist = powerlaw(low=1e30, high=1e40, power=-1, shape=N_SRCS)
# pop.set_lum(model='gauss', mean=lum_dist, std=1e2)
pop.set_lum(model='powerlaw', low=1e30, high=1e40, power=-1)

# Add a distribution of Poisson burst rates
rate_dist = log10normal(0.1, 2, N_SRCS)  # per day
pop.set_time(model='poisson', rate=rate_dist)

# Add a distribution of width rates
# w_eff_mean_dist = log10normal(1, 2, N_SRCS)
# pop.set_w(model='gauss', mean=w_eff_mean_dist, std=2)
pop.set_w(model='constant', value=1)

# Generate population
pop.generate()

# Survey the high fluences
survey = Survey('perfect', n_days=N_DAYS)
survey.mount_type = 'transit'
survey.set_beam(model='perfect')
survey.snr_limit = 1e1
survey.t_obs = 60*60  # seconds

# Check the burst rate
full_surv_pop = SurveyPopulation(pop, survey)
full_surv_pop.name = 'full'
time = full_surv_pop.frbs.time
print(f'# in full survey population: {full_surv_pop.n_bursts()}')

# Fractional mask
frac_surv_pop = deepcopy(full_surv_pop)
frac_surv_pop.name = 'fraction'
mask = (frac_surv_pop.frbs.time < init_surv_time_frac*N_DAYS)
frac_surv_pop.frbs.apply(mask)
print(f'# in fractional survey population: {frac_surv_pop.n_bursts()}')

plot(pop, full_surv_pop, frac_surv_pop)
