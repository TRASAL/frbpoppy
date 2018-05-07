"""Testing a perfect survey detecting everything."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot

days = 3
n_per_day = 5000

# Generate FRB population
population = generate(n_per_day*days,
                      days=days,
                      lum_dist_pars=[1e45, 1e45, 0.],
                      z_max=5.0,
                      pulse=[1, 1],
                      si_pars=[0., 0.],
                      repeat=0.0)

# Observe FRB population
result = observe(population, 'PERFECT')

# Plot populations
plot(population, result)
