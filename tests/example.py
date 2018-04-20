"""Short example of how frbpoppy works."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot

days = 7
n_per_day = 5000

# Generate FRB population
population = generate(n_per_day*days,
                      days=days,
                      lum_dist_pars=[1e42, 1e50, 0],
                      z_max=2.5,
                      pulse=[0.1, 10],
                      repeat=0.0)

# Observe FRB population
result = observe(population, 'APERTIF')

# Plot populations
plot(population, result)
