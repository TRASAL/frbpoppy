"""Code to inspect a population."""
from frbpoppy import unpickle, plot

pop = unpickle('perfect')
plot(pop, mute=False, frbcat=False)
