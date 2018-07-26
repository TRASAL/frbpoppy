"""Test predictions for an apertif population."""
from frbpoppy.do_plot import plot
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.paths import paths

old = paths.results()
paths.results(old + 'apertif_predictions/')

days = 90
n_per_day = 5000

# Generate FRB population
pop = generate(n_per_day*days,
               days=days,
               lum_range=[1e45, 1e50],
               lum_index=-1.5)

# Observe FRB population
surv_pop = observe(pop, 'APERTIF-FLY', gain_pattern='airy')

# Observe FRB population
iab_pop = observe(pop, 'APERTIF-IAB-12', gain_pattern='airy')

# Plot FRB populations
plot(pop, surv_pop, iab_pop, mute=False)
