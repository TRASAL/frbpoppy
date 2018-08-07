"""Testing new implementation of frbpoppy classes."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

pop = CosmicPopulation(1000, days=2, lum_range=[1e40, 1e45])
survey = Survey('HTRU', gain_pattern='parkes')
surv_pop = SurveyPopulation(pop, survey)
rates = surv_pop.rates(scale_area=False)
print(rates)
