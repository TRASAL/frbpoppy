"""Test repeaters."""
import numpy as np
from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, pprint

N_DAYS = 14

# Repeater population
n_gen = 1e4
cosmic = RepeaterPopulation.simple(n_gen)
cosmic.n_days = N_DAYS
cosmic.lum_min = 1e40
cosmic.lum_max = 1e45
cosmic.lum_rep_model = 'independent'
cosmic.times_rep_model = 'clustered'
cosmic.generate()

# Survey
survey = Survey('chime')
survey.n_days = N_DAYS
survey.beam_pattern = 'chime'

# Survey Population
surv_pop = SurveyPopulation(cosmic, survey)

# Check detections
frbs = surv_pop.frbs
frac_sources = surv_pop.n_sources()/np.count_nonzero(~np.isnan(cosmic.frbs.ra))
frac_bursts = surv_pop.n_bursts()/np.count_nonzero(~np.isnan(cosmic.frbs.time))
pprint(f'% sources det: {100*frac_sources:.4}% ({surv_pop.n_sources()} seen)')
pprint(f'% bursts det: {100*frac_bursts:.4}% ({surv_pop.n_bursts()} seen)')
