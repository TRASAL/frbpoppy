"""How to access frb population parameters."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

cosmic_pop = CosmicPopulation.simple(1e3, repeaters=True, generate=True)
dm = cosmic_pop.frbs.dm  # Get dispersion measure values

survey_pop = SurveyPopulation(cosmic_pop, Survey('perfect'))
survey_dm = survey_pop.frbs.dm  # Also works for SurveyPopulations

# Inbuilt functions help with detection rates
print(f'Number of sources: {survey_pop.n_sources()}')
print(f'Number of bursts: {survey_pop.n_bursts()}')
print(f'Number of repeating sources: {survey_pop.n_repeaters()}')
print(f'Number of one-off sources: {survey_pop.n_oneoffs()}')
