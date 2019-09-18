# -*- coding: future_fstrings -*-
"""Short example of how frbpoppy works."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

PLOT = False

# Generate an FRB population
cosmic_pop = CosmicPopulation(1e5, name='example', days=0.23)

# Setup a survey
survey = Survey('htru')

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)

# Check the detection rates
print(survey_pop.rates())

# Plot populations
if PLOT:
    plot(cosmic_pop, survey_pop, frbcat='parkes')
