"""Short example of how frbpoppy works."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

# Generate an FRB population
population = CosmicPopulation(5000*7, days=7, name='example_rates')

# Setup a survey
survey = Survey('PERFECT', gain_pattern='perfect')

# Observe the FRB population
result = SurveyPopulation(population, survey)

# Print rates
print(result.rates())

# Plot populations
plot(result, frbcat=False, mute=False)
