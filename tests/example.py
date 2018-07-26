"""Short example of how frbpoppy works."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot

days = 1
n_per_day = 5000

# Generate FRB population
population = generate(n_per_day*days,
                      days=days,
                      name='example')

# Observe FRB population
result = observe(population, 'APERTIF')

# Plot populations
plot(population, result)
