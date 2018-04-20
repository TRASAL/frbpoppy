"""Code to inspect a population."""
import os.path

from frbpoppy.population import unpickle
from frbpoppy.paths import paths
from frbpoppy.log import pprint

pop_path = os.path.join(paths.results(), 'population_initial.p')

pop = unpickle(pop_path)

pprint(pop)

for src in pop.sources:
    pprint(src.dm_host, src.dm_igm, src.dm_mw)
