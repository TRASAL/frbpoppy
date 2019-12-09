"""Survey brightness distribution."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

cosmic_pop = CosmicPopulation(1e4, name='simple', n_days=0.23,
                             lum_range=[1e45, 1e45], z_max=0.05,
                             lum_index=0, w_range=[1.,1.], w_sigma=0,
                             dm_mw_model='zero', dm_host_mu=0, dm_host_sigma=0,
                             dm_igm_sigma=None)
survey = Survey('perfect')
print(survey.gain_pattern)
surv_pop = SurveyPopulation(cosmic_pop, survey)
plot(surv_pop, frbcat=False)
