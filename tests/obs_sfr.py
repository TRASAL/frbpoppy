"""Test FRB sources following star forming rate."""
from frbpoppy.do_populate import generate
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle

MAKE = True

if MAKE:
    days = 7
    n_per_day = 5000

    # Generate population following a constant number density / comoving volume
    pop_cst = generate(n_per_day*days,
                       days=days,
                       z_max=6.0,
                       n_model='vol_co',
                       name='vol_co')

    # Generate population following star forming rate
    pop_sfr = generate(n_per_day*days,
                       days=days,
                       z_max=8.,
                       n_model='sfr',
                       name='sfr')

else:
    pop_cst = unpickle('vol_co')
    pop_sfr = unpickle('sfr')

# Plot populations
plot(pop_cst, pop_sfr)
