"""Test repeater survey detections."""
from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, plot, pprint

cosmic = RepeaterPopulation.simple(1e4)
cosmic.n_days = 2
cosmic.lum_rep_model = 'independent'
cosmic.generate()
survey = Survey('perfect-small')
survey.n_days = 2
# survey.snr_limit = 5e9
survey.strategy = 'regular'
surv_pop = SurveyPopulation(cosmic, survey)
frbs = surv_pop.frbs
pprint(f'# sources: {surv_pop.n_sources()}')
pprint(f'# bursts: {surv_pop.n_bursts()}')
# plot(surv_pop, mute=False)
