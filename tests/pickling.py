"""Test the pickling of a population."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

# Generate an FRB population
days = 1
population = CosmicPopulation(days*5000, days=days, name='test-cosmic')

# Setup a survey
survey = Survey('APERTIF')

# Observe the FRB population
surv_pop = SurveyPopulation(population, survey, name='test-survey')

# Attempt to pickle populations
population.save(extention='p')
surv_pop.save(extention='p')
