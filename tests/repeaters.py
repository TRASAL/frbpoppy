"""Play around with a simple repeater population."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot, pprint

MAX_DAYS = 1
PLOT = False

# Set up a fairly simple repeater population
r = CosmicPopulation.simple(1e5, n_days=MAX_DAYS, repeaters=True)
r.set_dist(z_max=0.01)
r.set_direction(min_dec=-10, max_dec=10)
r.set_time(model='regular', lam=1)
r.set_lum(model='powerlaw', low=1e45, high=1e45, power=0,
          per_source='different')
r.set_dm_host(model='gauss', mu=100, sigma=100)
r.set_dm_igm(model='ioka', slope=950, sigma=0)
r.set_dm_mw(model='ne2001')
r.set_dm(mw=True, igm=True, host=True)
r.generate()

# Set up survey
survey = Survey('apertif', n_days=MAX_DAYS)
survey.set_beam(model='perfect')

# Construct SurveyPopulation
surv_pop = SurveyPopulation(r, survey)

# Give some information on the detections
frac_sources = surv_pop.n_sources()/r.n_sources()
frac_bursts = surv_pop.n_bursts()/r.n_bursts()
pprint(f'{100*frac_sources:.4}% ({surv_pop.n_sources()}) sources seen')
pprint(f'{100*frac_bursts:.4}% ({surv_pop.n_bursts()}) bursts seen')

# Plot results
if PLOT:
    plot(surv_pop, mute=False)
