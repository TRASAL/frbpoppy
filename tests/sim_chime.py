"""Simulate large CHIME population."""
from frbpoppy import LargePopulation, CosmicPopulation, log10normal, Survey
from frbpoppy import plot

N_SRCS = 100
N_DAYS = 14
RATE = 1  # per day

r = CosmicPopulation(N_SRCS, n_days=N_DAYS, repeaters=True)
r.set_dist(model='vol_co', z_max=0.01)
r.set_dm_host(model='gauss', mean=100, std=0)
r.set_dm_igm(model='ioka', slope=1000, std=0)
r.set_dm(mw=False, igm=True, host=True)
r.set_emission_range(low=100e6, high=10e9)
r.set_lum(model='powerlaw', per_source='different', low=1e35, high=1e45,
          power=-1)
r.set_si(model='gauss', mean=-1.4, std=0)
r.set_w(model='log10normal', per_source='different', mean=0.01, std=5)
# rate = log10normal(RATE, 2, int(N_SRCS))
# r.set_time(model='poisson', rate=rate)
r.set_time(model='regular', rate=RATE)

s = Survey('chime', n_days=N_DAYS)
s.set_beam(model='perfect-small')

surv_pop = LargePopulation(r, s, run=True).pops[0]
surv_pop.name = 'chime'
surv_pop.save()

# import IPython; IPython.embed()
print(surv_pop.source_rate)
print(surv_pop.burst_rate)
# plot(surv_pop)
