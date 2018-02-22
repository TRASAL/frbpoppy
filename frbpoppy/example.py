from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from do_plot import plot

days = 30
n_per_day = 5000

# Generate FRB population
population = generate(n_per_day*days,
                      days=days,
                      lum_dist_pars=[1e40, 1e50, -1.5],
                      z_max=2.5,
                      pulse=[0.1, 10],
                      repeat=0.0)

# Observe FRB populations
surveys = ['APERTIF',
           'HTRU',
           'UTMOST-1D']

results = []

for s in surveys:
    results.append(observe(population, s))

# Plot populations
plot(population, *results, mute=False)
