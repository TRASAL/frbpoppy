"""Test whether distributions for simple populations hold true."""
from frbpoppy.do_plot import plot
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe

days = 1
n_per_day = 20000

# Generate FRB population
pop = generate(n_per_day*days,
               days=days,
               lum_dist_pars=[1e42, 1e42, -1.0],
               z_max=0.1,
               dm_pars=[0, 1200],
               electron_model='zero',
               emission_pars=[10e6, 10e9],
               pulse=[5, 5],
               si_pars=[0., 0.],
               repeat=0.0)

# Observe FRB population
surv_pop = observe(pop, 'APERTIF', pattern='tophat')

# # Plot populations
# plot(pop, surv_pop)
