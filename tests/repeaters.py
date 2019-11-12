"""Test repeater population."""
import numpy as np

from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, plot
from frbpoppy import split_pop, pprint

DAYS = 5

r = RepeaterPopulation.simple(int(1e5))
r.times_rep_model = 'clustered'
r.days = DAYS
r.generate()

survey = Survey('apertif', strategy='regular', n_days=DAYS)
survey.gain_pattern = 'apertif'
survey.snr_limit = 1.

pops = []
surv_pop = SurveyPopulation(r, survey)

# Split population into seamingly one-off and repeater populations
mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
pop_ngt1, pop_nle1 = split_pop(surv_pop, mask)
pop_ngt1.name += ' (> 1 burst)'
pop_nle1.name += ' (1 burst)'

pops.append(pop_nle1)
pops.append(pop_ngt1)

pprint(f'{len(surv_pop.frbs.index)} FRBs detected')

if len(surv_pop.frbs.index) > 100:
    plot(*pops, frbcat=False, mute=False, show=True)
