"""Test predictions for an apertif population."""
from frbpoppy.do_plot import plot
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.paths import paths

old = paths.results()
paths.results(old + 'apertif_predictions/')

days = 1
n_per_day = 10000

# Generate FRB population
pop = generate(n_per_day*days,
               days=days,
               lum_dist_pars=[1e39, 1e45, -1.0],
               z_max=2.0,
               repeat=0.0)

# Observe FRB population
surv_pop = observe(pop, 'APERTIF-IAB-10', pattern='tophat')

# Save populations
pop.save()
surv_pop.save()

# Plot results
files = [paths.results()+'population_initial.csv',
         paths.results()+'population_apertif.csv']

plot(files=files, mute=False)
