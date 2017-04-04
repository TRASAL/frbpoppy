from populate import generate
from dosurvey import observe
from bokeh_server import plot

# Generate FRB population
population = generate(100)

# Observe FRB population
survey_population = observe(population, 'PMSURV')

# Plot populations
plot(population, survey_population)
