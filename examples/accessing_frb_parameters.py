"""How to access frb population parameters."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

cosmic_pop = CosmicPopulation.complex(1e5, generate=True)
dm = cosmic_pop.frbs.dm  # Get dispersion measure values

survey_pop = SurveyPopulation(cosmic_pop, Survey('apertif'))
survey_dm = survey_pop.frbs.dm  # Also works for SurveyPopulations
