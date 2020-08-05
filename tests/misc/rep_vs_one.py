"""Generate a repeater population and split into repeaters and one-offs."""
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import split_pop, pprint, plot

DAYS = 1

r = CosmicPopulation.simple(int(1e4), n_days=DAYS, repeaters=True)
r.set_time(model='regular', rate=2)
r.set_lum(model='powerlaw', low=1e40, high=1e45, per_source='different')
r.generate()

survey = Survey('chime', n_days=DAYS)
survey.set_beam(model='perfect')
surv_pop = SurveyPopulation(r, survey)

# Split population into seamingly one-off and repeater populations
mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
pop_ngt1, pop_nle1 = split_pop(surv_pop, mask)
pop_ngt1.name += ' (> 1 burst)'
pop_nle1.name += ' (1 burst)'

pops = [pop_nle1, pop_ngt1]

pprint(f'{surv_pop.n_sources()} sources detected')
pprint(f'{surv_pop.n_bursts()} bursts detected')

plot(surv_pop)
