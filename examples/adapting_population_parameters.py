"""Example of changing parameters."""
from frbpoppy import CosmicPopulation, Survey

# Set up a population
cosmic_pop = CosmicPopulation(1e4, generate=False)
# ... or adapt the population per parameter, e.g.
cosmic_pop.set_dist(z_max=2.5)
# Or to adapt the luminosity
cosmic_pop.set_lum(model='powerlaw', low=1e44, high=1e45)
# Generate the population
cosmic_pop.generate()

# Setup a survey with arguments ...
survey = Survey('apertif', beam_pattern='airy', n_sidelobes=2)
# ... or adapt the survey later with
survey.snr_limit = 2

# For a full list of available arguments or parameters check the classes as
# defined in /frbpoppy/*.py
