"""Survey brightness distribution."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

from pointings import plot_coordinates

MAX_DAYS = 10

r = CosmicPopulation.simple(1e4, n_days=MAX_DAYS, repeaters=True)
r.set_dist(z_max=0.01)
r.set_direction(min_dec=-5, max_dec=5)
r.set_time(model='regular', lam=1)
r.set_lum(model='powerlaw', low=1e45, high=1e45, power=0,
          per_source='different')
r.set_dm_host(model='gauss', mu=100, sigma=100)
# In between Cordes review & Petroff review
r.set_dm_igm(model='ioka', slope=950, sigma=0)
r.set_dm_mw(model='ne2001')
r.set_dm(mw=True, igm=True, host=True)
r.generate()

# Set up survey
survey = Survey('apertif', n_days=MAX_DAYS)
survey.set_beam(model='perfect')
survey.set_pointings(ra=[100, 200, 300], dec=[0, 0, 0])

surv_pop = SurveyPopulation(r, survey)

frbs = surv_pop.frbs

# plot_coordinates(*survey.pointings)
# import IPython; IPython.embed()
plot(surv_pop, mute=False)
