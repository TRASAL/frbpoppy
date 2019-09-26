"""Short example of how frbpoppy works.

The first time you run frbpoppy, a series of cosmological databases will be
constructed to set up subsequent runs. This first run can take ~2h on a 4 core
machine. Subsequent runs will take mere seconds.
"""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

# Generate an FRB population
cosmic_pop = CosmicPopulation(1e5, name='example', days=0.23)

# Setup a survey
survey = Survey('htru')

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)

# Check the detection rates
print(survey_pop.rates())

# Plot populations in a browser
plot(cosmic_pop, survey_pop, frbcat='parkes')
