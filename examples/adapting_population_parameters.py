"""Example of changing parameters."""
from frbpoppy import CosmicPopulation, Survey
import frbpoppy.gen_dists as gd

# Set up a population
cosmic_pop = CosmicPopulation(1e4, generate=False)
# ... or adapt the population per parameter, e.g.
cosmic_pop.set_dist(z_max=2.5)
# Or to adapt the luminosity
cosmic_pop.set_lum(model='powerlaw', low=1e44, high=1e45)
# Generate the population
cosmic_pop.generate()

# Setup a survey
survey = Survey('apertif')
# and adapt the survey with
survey.set_beam('airy', n_sidelobes=2)
survey.snr_limit = 2

# For a full list of available arguments or parameters check the classes as
# defined in /frbpoppy/*.py

# Some more difficult examples

# Use a convolution of two distributions
cosmic_pop = CosmicPopulation.simple(n_srcs=1e4, repeaters=True)
# Draw the mean value per source from a normal distribution
mean_dist = gd.trunc_norm(mean=1, std=2, shape=int(1e4))
# And use those means as a _means_ to generate from a new Gaussian
# distribution per source
cosmic_pop.set_w(model='gauss', per_source='different',
                 mean=mean_dist, std=0.01*mean_dist)
# And generate afterwards
cosmic_pop.generate()
