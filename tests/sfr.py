"""Test FRB sources following star forming rate."""
from frbpoppy.do_populate import generate
from frbpoppy.do_plot import plot

days = 7
n_per_day = 5000

# Generate population following a constant number density per comoving volume
pop_cst = generate(n_per_day*days,
                   days=days,
                   lum_dist_pars=[1e40, 1e50, -1.5],
                   z_max=4.,
                   pulse=[0.1, 10],
                   repeat=0.0,
                   n_model='constant',
                   name='constant')

# Generate population following star forming rate
pop_sfr = generate(n_per_day*days,
                   days=days,
                   lum_dist_pars=[1e40, 1e50, -1.5],
                   z_max=4.,
                   pulse=[0.1, 10],
                   repeat=0.0,
                   n_model='sfr',
                   name='sfr')

# Plot populations
plot(pop_cst, pop_sfr)
