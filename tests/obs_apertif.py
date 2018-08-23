"""Test predictions for an apertif population."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, paths, plot

old = paths.results()
paths.results(old + 'apertif_predictions/')

days = 90
n_per_day = 5000

# Generate FRB population
pop = CosmicPopulation(n_per_day*days,
                       days=days,
                       lum_range=[1e45, 1e50],
                       lum_index=-1.5)

iab_survey = Survey('APERTIF-IAB-12', gain_pattern='apertif')
iab_pop = SurveyPopulation(pop, iab_survey)

fly_survey = Survey('APERTIF-FLY', gain_pattern='apertif')
fly_pop = SurveyPopulation(pop, fly_survey)

# Plot FRB populations
plot(pop, iab_pop, fly_pop, mute=False)
