"""Code to inspect a population."""
import os
import sys
sys.path.append("..")
from frbpoppy.population import Population, unpickle

pop_path = './simple_pop.p'

pop = unpickle(pop_path)

print(pop)

for src in pop.sources:
    print(src.dm_host, src.dm_igm, src.dm_mw)
