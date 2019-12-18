"""Test repeaters."""
import numpy as np
from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, pprint

# Repeater population
n_gen = 1e4
cosmic = RepeaterPopulation.simple(n_gen)
cosmic.n_days = 2
cosmic.lum_rep_model = 'independent'
cosmic.generate()

# Survey
survey = Survey('perfect')
survey.n_days = 1
survey.snr_limit = 0
survey.beam_pattern = 'gaussian'

# Survey Population
surv_pop = SurveyPopulation(cosmic, survey)

# Check detections
frbs = surv_pop.frbs
frac_sources = surv_pop.n_sources()/np.count_nonzero(~np.isnan(cosmic.frbs.ra))
frac_bursts = surv_pop.n_bursts()/np.count_nonzero(~np.isnan(cosmic.frbs.time))
pprint(f'% sources det: {100*frac_sources:.4}%')
pprint(f'% bursts det: {100*frac_bursts:.4}%')
