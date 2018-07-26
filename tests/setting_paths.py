"""How to change paths."""
from frbpoppy.do_plot import plot
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.log import pprint
from frbpoppy.paths import paths
import os.path

days = 1
n_per_day = 10000

pprint('Old path:', paths.populations())

# Change path to which results are saved
paths.results(os.path.expanduser("~/Downloads/frbpoppy/"))

pprint('New path:', paths.populations())

# Generate FRB population
pop = generate(n_per_day*days,
               days=days,
               lum_range=[1e40, 1e46])

# Observe FRB population
surv_pop = observe(pop, 'APERTIF')

# Plot populations
plot(pop, surv_pop)
