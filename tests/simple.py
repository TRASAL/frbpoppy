"""Short example of how frbpoppy works."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

# Generate an FRB population
cosmic_pop = CosmicPopulation.simple(1e4, generate=True)

# Setup a survey
survey = Survey('perfect')

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)

# Check the detection rates
print(survey_pop.rates())

# Plot populations
plot(cosmic_pop, survey_pop, frbcat=False, mute=False)
