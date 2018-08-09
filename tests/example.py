"""Short example of how frbpoppy works."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

# Generate an FRB population
population = CosmicPopulation(10000, days=2, name='example')

# Setup a survey
survey = Survey('APERTIF')

# Observe the FRB population
result = SurveyPopulation(population, survey)

# Plot populations
plot(population, result)
