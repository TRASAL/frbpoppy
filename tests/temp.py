# -*- coding: future_fstrings -*-
"""Survey brightness distribution."""
from frbpoppy import plot, unpickle

large = unpickle('alpha_large_-2.5_askap-fly')

plot(large, frbcat='askap-fly')
