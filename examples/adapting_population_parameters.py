"""Example of changing parameters."""
from frbpoppy import CosmicPopulation, Survey

# Set up a population with arguments such as z_max ...
cosmic_pop = CosmicPopulation(1e4, z_max=0.01, generate=False)
# ... or adapt the population with e.g.
cosmic_pop.z_max = 2.5
# Generate the population
cosmic_pop.generate()

# Setup a survey with arguments ...
survey = Survey('apertif', beam_pattern='airy', n_sidelobes=2)
# ... or adapt the survey later with
survey.snr_limit = 2

# For a full list of available arguments or parameters check the classes as
# defined in /frbpoppy/*.py
