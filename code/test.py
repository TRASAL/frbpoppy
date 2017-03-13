from populate import generate
from dosurvey import observe
from plot import plot_pop
from log import pprint

# Generate FRB population
population = generate(10)

# Observe FRB population
survey_population = observe(population, 'WHOLESKY')

plot_pop(pops=[population, survey_population], show=False)
