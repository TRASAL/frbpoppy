"""Testing a perfect survey detecting everything."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot

days = 3
n_per_day = 5000

# Generate FRB population
population = generate(n_per_day*days, days=days)

# Observe FRB population
result = observe(population, 'PERFECT', gain_pattern='perfect')

# Plot populations
plot(population, result)
