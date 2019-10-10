"""Test survey strategies."""
import numpy as np
from frbpoppy import (CosmicPopulation, RepeaterPopulation, Survey,
                      SurveyPopulation, plot)

pop = RepeaterPopulation.simple(1e4)
pop.times_rep_model = 'clustered'
pop.generate()
s = Survey('chime', n_days=1, strategy='follow-up')
surv_pop = SurveyPopulation(pop, s)
print(surv_pop.n_bursts())
if surv_pop.n_bursts() > 1:
    plot(surv_pop, frbcat=False)
