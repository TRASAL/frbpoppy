"""Test FRB sources following star forming rate."""
from frbpoppy.do_populate import generate
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle

MAKE = False

if MAKE:
    days = 7
    n_per_day = 5000

    # Generate population following a constant number density / comoving volume
    pop_cst = generate(n_per_day*days,
                       days=days,
                       lum_dist_pars=[1e40, 1e50, 0.],
                       z_max=4.,
                       pulse=[0.1, 10],
                       repeat=0.0,
                       n_model='constant',
                       name='constant')

    # Generate population following star forming rate
    pop_sfr = generate(n_per_day*days,
                       days=days,
                       lum_dist_pars=[1e40, 1e50, 0.],
                       z_max=4.,
                       pulse=[0.1, 10],
                       repeat=0.0,
                       n_model='sfr',
                       name='sfr')

else:
    pop_cst = unpickle('constant')
    pop_sfr = unpickle('sfr')

# Plot populations
plot(pop_cst, pop_sfr)
