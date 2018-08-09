"""Generate a standard population."""

from frbpoppy import CosmicPopulation, plot

days = 100

pop = CosmicPopulation(5000*days,
                       days=days,
                       lum_index=-0.5,
                       name='standard')
